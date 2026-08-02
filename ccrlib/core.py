"""
ccrlib.core
===========
A compact counterparty-credit-exposure Monte-Carlo engine used to study how well
current exposure (CE) and mark-to-market (MTM) proxy the change in Effective EPE
(EEPE) under market shocks.

Design
------
* Risk factors: geometric Brownian motion (equity / FX / commodity) and a
  Vasicek-style Ornstein-Uhlenbeck short rate (interest rate). Optionally, each
  option-bearing asset carries a SEPARATE stochastic *implied volatility* factor
  (mean-reverting log-OU) used only to reprice options -- distinct from the
  *realized* volatility that diffuses the underlying.
* Common random numbers: draws are built once and reused across every scenario so
  Monte-Carlo noise cancels in the deltas (dCE, dMTM, dEEPE).
* Scenarios ("shocks") perturb the initial risk-factor levels, optionally scale
  the realized simulation volatility, and optionally shift implied-vol levels.

Everything here is deliberately transparent (single netting set, flat discount
curve, non-discounted EE) -- it is a teaching / research engine, not a booking
system.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field
from math import erf, sqrt

# --------------------------------------------------------------------------- #
#  Small analytics
# --------------------------------------------------------------------------- #
_SQRT2 = sqrt(2.0)
_erf = np.vectorize(erf)

def norm_cdf(x):
    return 0.5 * (1.0 + _erf(np.asarray(x, float) / _SQRT2))

def norm_pdf(x):
    return np.exp(-0.5 * np.asarray(x, float) ** 2) / np.sqrt(2 * np.pi)

def bs_price(S, K, r, q, sig, tau, call=True):
    """Black-Scholes; vectorised over S and/or sig. Handles tau -> 0 (intrinsic)."""
    S = np.asarray(S, float); sig = np.asarray(sig, float)
    if tau <= 1e-9:
        return np.maximum(S - K, 0.0) if call else np.maximum(K - S, 0.0)
    vs = sig * np.sqrt(tau)
    d1 = (np.log(S / K) + (r - q + 0.5 * sig * sig) * tau) / vs
    d2 = d1 - vs
    if call:
        return S * np.exp(-q * tau) * norm_cdf(d1) - K * np.exp(-r * tau) * norm_cdf(d2)
    return K * np.exp(-r * tau) * norm_cdf(-d2) - S * np.exp(-q * tau) * norm_cdf(-d1)

def bachelier_price(F, K, sig_n, tau, DF, call=True):
    """Normal-model option (for interest-rate optionality). sig_n is absolute vol."""
    F = np.asarray(F, float); sig_n = np.asarray(sig_n, float)
    if tau <= 1e-9:
        payoff = np.maximum(F - K, 0.0) if call else np.maximum(K - F, 0.0)
        return DF * payoff
    s = sig_n * np.sqrt(tau)
    d = (F - K) / s
    if call:
        return DF * ((F - K) * norm_cdf(d) + s * norm_pdf(d))
    return DF * ((K - F) * norm_cdf(-d) + s * norm_pdf(d))


# --------------------------------------------------------------------------- #
#  Risk-factor and instrument definitions
# --------------------------------------------------------------------------- #
@dataclass
class Asset:
    name: str
    process: str            # 'gbm' or 'ou'
    spot: float
    vol: float              # realized / diffusion vol (lognormal for gbm; absolute for ou)
    q: float = 0.0          # carry/yield (gbm): drift = r_disc - q
    kappa: float = 0.10     # ou mean-reversion speed
    theta: float = 0.03     # ou long-run level
    # optional stochastic implied vol attached to this asset (for options):
    imp0: float | None = None
    iv_kappa: float = 3.0
    iv_nu: float = 0.8      # vol-of-vol (log-OU)
    iv_rho: float = -0.5    # leverage correlation to own spot Brownian


@dataclass
class Trade:
    asset: str
    kind: str               # 'fwd','swap','call','put','rate_call','rate_put'
    strike: float
    notional: float = 1.0
    maturity: float = 1.0
    duration: float = 1.0   # annuity/DV01 scaling for rate instruments


@dataclass
class Scenario:
    """An instantaneous market shock applied to the initial state."""
    spot_mult: dict = field(default_factory=dict)   # gbm: multiplicative level shock (1+ds)
    rate_add: dict = field(default_factory=dict)     # ou: additive level shock
    vol_mult: float = 1.0                             # scales REALIZED diffusion vol of spots/rate
    iv_add: dict = field(default_factory=dict)        # additive implied-vol LEVEL shock (base implied vol)
    iv_vov_mult: float = 1.0                           # scales the VOL-OF-VOL (diffusion vol of implied vol)
    iv_persist: bool = True                           # shift iv mean-reversion target too


R_DISC = 0.03   # flat discount rate


# --------------------------------------------------------------------------- #
#  Simulator
# --------------------------------------------------------------------------- #
class Simulator:
    def __init__(self, assets, corr, times, n_paths=20000, seed=0):
        self.assets = assets
        self.names = [a.name for a in assets]
        self.idx = {a.name: i for i, a in enumerate(assets)}
        self.corr = np.asarray(corr, float)
        self.times = np.asarray(times, float)
        self.dt = float(times[1] - times[0])
        self.nsteps = len(times) - 1
        self.opt_assets = [a for a in assets if a.imp0 is not None]
        self.n = n_paths
        self._build_draws(seed)

    def _build_draws(self, seed):
        rng = np.random.default_rng(seed)
        h = self.n // 2
        N = len(self.assets)
        Zs = rng.standard_normal((h, self.nsteps, N))
        Zs = np.concatenate([Zs, -Zs], axis=0)               # antithetic
        L = np.linalg.cholesky(self.corr)
        self.Zs = Zs @ L.T                                    # correlated spot/rate shocks
        no = len(self.opt_assets)
        Eta = rng.standard_normal((h, self.nsteps, max(no, 1)))
        self.Eta = np.concatenate([Eta, -Eta], axis=0)        # independent iv shocks

    def simulate(self, sc: Scenario):
        n = self.Zs.shape[0]
        sdt = np.sqrt(self.dt)
        factors = {}
        # initialise levels with shocks
        lvl = {}
        for a in self.assets:
            if a.process == 'gbm':
                lvl[a.name] = np.full(n, a.spot * sc.spot_mult.get(a.name, 1.0))
            else:  # ou
                lvl[a.name] = np.full(n, a.spot + sc.rate_add.get(a.name, 0.0))
        paths = {name: np.empty((self.nsteps + 1, n)) for name in self.names}
        for name in self.names:
            paths[name][0] = lvl[name]
        # implied-vol paths
        ivp = {}
        iv_theta = {}
        for j, a in enumerate(self.opt_assets):
            i0 = a.imp0 + sc.iv_add.get(a.name, 0.0)
            iv_theta[a.name] = i0 if sc.iv_persist else a.imp0
            ivp[a.name] = np.empty((self.nsteps + 1, n))
            ivp[a.name][0] = np.full(n, i0)
        lnsig = {a.name: np.log(ivp[a.name][0]) for a in self.opt_assets}

        for k in range(1, self.nsteps + 1):
            for a in self.assets:
                i = self.idx[a.name]
                Z = self.Zs[:, k - 1, i]
                if a.process == 'gbm':
                    sig = a.vol * sc.vol_mult                       # diffusion (width) scales with the vol shock
                    # Convexity uses the BASE vol so that shocking the simulation volatility changes only the
                    # dispersion (envelope width), NOT the central/median path. This keeps the LEVEL shock
                    # identical between methodology A (vol shocked) and B (vol fixed) -- they differ only in width.
                    drift = (R_DISC - a.q - 0.5 * a.vol * a.vol) * self.dt
                    lvl[a.name] = lvl[a.name] * np.exp(drift + sig * sdt * Z)
                else:  # ou (exact)
                    sig = a.vol * sc.vol_mult
                    e_a = np.exp(-a.kappa * self.dt)
                    sd = sig * np.sqrt((1 - np.exp(-2 * a.kappa * self.dt)) / (2 * a.kappa))
                    lvl[a.name] = lvl[a.name] * e_a + a.theta * (1 - e_a) + sd * Z
                paths[a.name][k] = lvl[a.name]
            for j, a in enumerate(self.opt_assets):
                i = self.idx[a.name]
                Wv = a.iv_rho * self.Zs[:, k - 1, i] + np.sqrt(1 - a.iv_rho ** 2) * self.Eta[:, k - 1, j]
                lnsig[a.name] = lnsig[a.name] + a.iv_kappa * (np.log(iv_theta[a.name]) - lnsig[a.name]) * self.dt \
                    + (a.iv_nu * sc.iv_vov_mult) * sdt * Wv       # vol-of-vol scaled by the diffusion-vol shock
                ivp[a.name][k] = np.exp(lnsig[a.name])
        return paths, ivp


# --------------------------------------------------------------------------- #
#  Pricing
# --------------------------------------------------------------------------- #
def price_portfolio(trades, assets, paths, ivp, times):
    """Return netted portfolio value V of shape (nsteps+1, n_paths)."""
    amap = {a.name: a for a in assets}
    nsteps = len(times) - 1
    V = np.zeros_like(next(iter(paths.values())))
    for tr in trades:
        a = amap[tr.asset]
        for k in range(nsteps + 1):
            tau = tr.maturity - times[k]
            if tau <= 1e-9:
                continue
            S = paths[tr.asset][k]
            if tr.kind == 'fwd':
                val = tr.notional * (S * np.exp(-a.q * tau) - tr.strike * np.exp(-R_DISC * tau))
            elif tr.kind == 'swap':
                val = tr.notional * tr.duration * (S - tr.strike)          # linear in rate
            elif tr.kind in ('call', 'put'):
                iv = ivp[tr.asset][k] if tr.asset in ivp else a.vol
                val = tr.notional * bs_price(S, tr.strike, R_DISC, a.q, iv, tau, call=(tr.kind == 'call'))
            elif tr.kind in ('rate_call', 'rate_put'):
                iv = ivp[tr.asset][k] if tr.asset in ivp else a.vol
                DF = np.exp(-R_DISC * tau)
                val = tr.notional * tr.duration * bachelier_price(S, tr.strike, iv, tau, DF,
                                                                  call=(tr.kind == 'rate_call'))
            else:
                raise ValueError(tr.kind)
            V[k] = V[k] + val
    return V


# --------------------------------------------------------------------------- #
#  Exposure metrics
# --------------------------------------------------------------------------- #
def exposure_metrics(V, times, eepe_window=1.0):
    """Return dict with MTM, CE, EE(t), EffEE(t), EPE, EEPE, PFE95(t)."""
    E = np.maximum(V, 0.0)
    EE = E.mean(axis=1)
    EffEE = np.maximum.accumulate(EE)
    win = times <= eepe_window + 1e-9
    w = np.diff(times[win], prepend=0.0)
    Tw = times[win][-1]
    EPE = float(np.sum(EE[win] * w) / Tw)
    EEPE = float(np.sum(EffEE[win] * w) / Tw)
    MTM = float(V[0, 0])
    return dict(MTM=MTM, CE=max(MTM, 0.0), EE=EE, EffEE=EffEE,
                EPE=EPE, EEPE=EEPE, PFE95=np.percentile(E, 95, axis=1), times=times)
