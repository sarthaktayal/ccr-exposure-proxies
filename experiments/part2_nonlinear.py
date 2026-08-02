"""
Part 2 -- Non-linear portfolio (options across asset classes).
==============================================================
Long options on interest-rate, equity, FX and commodity underlyings. IMPLIED
VOLATILITY is modelled as a SEPARATE stochastic RISK FACTOR (mean-reverting
log-OU) used to reprice the options -- distinct from the REALIZED volatility that
diffuses the underlyings.

Three moneyness cases (options deep-ITM / ATM / OTM). Two shock families:
    * spot risk-factor shock   (methodology B: realized vol fixed;  A: realized vol shocked)
    * implied-vol risk-factor shock  (realized vol fixed)
We compare ΔEEPE against ΔCE and ΔMTM.

Outputs (../figures):
    part2_sim_envelopes.png  -- spot AND implied-vol path envelopes + portfolio, base vs shocked.
    part2_tracking.png       -- Δ metrics vs shock: spot shock (top) & implied-vol shock (bottom).
    part2_scatter.png        -- ΔEEPE vs ΔCE / ΔMTM with slopes (delta tracks, vega under-tracks).
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from ccrlib import Asset, Trade, Scenario, Simulator, price_portfolio, exposure_metrics, R_DISC

FIG = os.path.join(os.path.dirname(__file__), "..", "figures"); os.makedirs(FIG, exist_ok=True)

# ----- risk factors: spots + stochastic implied vols -------------------------
ASSETS = [
    Asset('rate', 'ou', 0.03, 0.010, kappa=0.10, theta=0.03, imp0=0.010, iv_kappa=3.0, iv_nu=0.6, iv_rho=-0.1),
    Asset('eq', 'gbm', 100.0, 0.20, q=0.02, imp0=0.20, iv_kappa=3.0, iv_nu=0.8, iv_rho=-0.7),
    Asset('fx', 'gbm', 1.10, 0.10, q=0.01, imp0=0.10, iv_kappa=3.0, iv_nu=0.5, iv_rho=-0.2),
    Asset('cm', 'gbm', 80.0, 0.25, q=0.01, imp0=0.25, iv_kappa=3.0, iv_nu=0.6, iv_rho=-0.3),
]
IMP0 = {a.name: a.imp0 for a in ASSETS}
CORR = np.array([[1, -.2, .1, .05], [-.2, 1, .25, .30], [.1, .25, 1, .40], [.05, .30, .40, 1.]])
TIMES = np.array([k / 12 for k in range(13)])          # monthly, 1y (options mature at 1y)
SIM = Simulator(ASSETS, CORR, TIMES, n_paths=24000, seed=2)
T = 1.0

# forwards (for strike placement)
FWD = {'eq': 100 * np.exp((R_DISC - 0.02) * T), 'fx': 1.10 * np.exp((R_DISC - 0.01) * T),
       'cm': 80 * np.exp((R_DISC - 0.01) * T), 'rate': 0.03}

def book(moneyness):
    mult = {'deep_itm': 0.85, 'atm': 1.00, 'otm': 1.15}[moneyness]
    return [Trade('rate', 'rate_call', FWD['rate'] * mult, 100.0, T, duration=4.5),
            Trade('eq', 'call', FWD['eq'] * mult, 0.50, T),
            Trade('fx', 'call', FWD['fx'] * mult, 50.0, T),
            Trade('cm', 'call', FWD['cm'] * mult, 0.30, T)]

BETA_VOL = 1.5
def sc_spot(s, shock_vol):
    return Scenario(spot_mult={'eq': 1 + s, 'fx': 1 + 0.5 * s, 'cm': 1 + s},
                    rate_add={'rate': 0.008 * s},
                    vol_mult=(1 + BETA_VOL * s) if shock_vol else 1.0)
def sc_iv(div_rel):
    return Scenario(iv_add={a: IMP0[a] * div_rel for a in IMP0}, vol_mult=1.0)

def metrics_for(trades, sc):
    paths, ivp = SIM.simulate(sc)
    V = price_portfolio(trades, ASSETS, paths, ivp, TIMES)
    return exposure_metrics(V, TIMES, eepe_window=1.0)

# ============================================================================ #
#  Figure 1 -- simulations after shock: spots + IMPLIED VOL (risk factor) + book
# ============================================================================ #
def fig_envelopes():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    base_p, base_iv = SIM.simulate(Scenario())
    sp_p, sp_iv = SIM.simulate(sc_spot(0.15, True))        # spot + realized-vol shock
    iv_p, iv_iv = SIM.simulate(sc_iv(0.30))                # implied-vol +30%
    trades = book('atm')
    Vb = price_portfolio(trades, ASSETS, base_p, base_iv, TIMES)
    Vs = price_portfolio(trades, ASSETS, sp_p, sp_iv, TIMES)
    Vi = price_portfolio(trades, ASSETS, iv_p, iv_iv, TIMES)

    def env(P): return np.percentile(P, 5, 1), np.percentile(P, 95, 1)
    def band(ax, P, c, lab):
        lo, hi = env(P); ax.fill_between(TIMES, lo, hi, color=c, alpha=0.12)
        ax.plot(TIMES, hi, color=c, lw=1.3); ax.plot(TIMES, lo, color=c, lw=1.3, label=lab)

    fig, axs = plt.subplots(2, 3, figsize=(15, 8.5))
    # equity spot + equity implied vol
    band(axs[0, 0], base_p['eq'], 'grey', 'base'); band(axs[0, 0], sp_p['eq'], 'C3', 'spot+vol shock')
    axs[0, 0].set(title='Equity SPOT envelope', xlabel='years')
    band(axs[0, 1], base_iv['eq'] * 100, 'grey', 'base'); band(axs[0, 1], iv_iv['eq'] * 100, 'C4', 'implied +30%')
    axs[0, 1].set(title='Equity IMPLIED VOL envelope (a risk factor!)', ylabel='vol %', xlabel='years')
    # commodity spot + implied vol
    band(axs[0, 2], base_p['cm'], 'grey', 'base'); band(axs[0, 2], sp_p['cm'], 'C3', 'spot+vol shock')
    axs[0, 2].set(title='Commodity SPOT envelope', xlabel='years')
    band(axs[1, 0], base_iv['cm'] * 100, 'grey', 'base'); band(axs[1, 0], iv_iv['cm'] * 100, 'C4', 'implied +30%')
    axs[1, 0].set(title='Commodity IMPLIED VOL envelope', ylabel='vol %', xlabel='years')
    # rate
    band(axs[1, 1], base_p['rate'] * 100, 'grey', 'base'); band(axs[1, 1], sp_p['rate'] * 100, 'C3', 'spot+vol shock')
    axs[1, 1].set(title='Interest-rate envelope (%)', xlabel='years')
    # portfolio
    band(axs[1, 2], Vb, 'grey', 'base'); band(axs[1, 2], Vs, 'C3', 'spot+vol shock'); band(axs[1, 2], Vi, 'C4', 'implied +30%')
    axs[1, 2].set(title='Portfolio value V (ATM options)', xlabel='years')
    for ax in axs.ravel(): ax.grid(alpha=.3); ax.legend(fontsize=7)
    fig.suptitle('Part 2 — Simulated risk factors incl. IMPLIED VOL as its own stochastic factor', fontsize=13)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "part2_sim_envelopes.png"), dpi=125); plt.close(fig)

# ============================================================================ #
#  Figures 2 & 3 -- tracking under spot shock and implied-vol shock
# ============================================================================ #
def fig_tracking():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    cases = [('deep_itm', 'ITM options'), ('atm', 'ATM options'), ('otm', 'OTM options')]
    Sspot = np.linspace(-0.12, 0.12, 13)
    Div = np.linspace(-0.30, 0.30, 13)
    res = {}
    for key, _ in cases:
        tr = book(key); base = metrics_for(tr, Scenario())
        d = dict(baseMTM=base['MTM'], baseCE=base['CE'], baseEEPE=base['EEPE'])
        # spot shock, methodologies B and A
        d['sMTM'] = []; d['sCE'] = []; d['sEEPE_B'] = []; d['sEEPE_A'] = []
        for s in Sspot:
            mB = metrics_for(tr, sc_spot(s, False)); mA = metrics_for(tr, sc_spot(s, True))
            d['sMTM'].append(mB['MTM'] - base['MTM']); d['sCE'].append(mB['CE'] - base['CE'])
            d['sEEPE_B'].append(mB['EEPE'] - base['EEPE']); d['sEEPE_A'].append(mA['EEPE'] - base['EEPE'])
        # implied-vol shock
        d['vMTM'] = []; d['vCE'] = []; d['vEEPE'] = []
        for dv in Div:
            m = metrics_for(tr, sc_iv(dv))
            d['vMTM'].append(m['MTM'] - base['MTM']); d['vCE'].append(m['CE'] - base['CE'])
            d['vEEPE'].append(m['EEPE'] - base['EEPE'])
        res[key] = {k: (np.array(v) if isinstance(v, list) else v) for k, v in d.items()}

    # Figure 2: metrics vs shock. Row 0 = spot shock, Row 1 = implied-vol shock
    fig, axs = plt.subplots(2, 3, figsize=(15, 9))
    for col, (key, title) in enumerate(cases):
        r = res[key]
        ax = axs[0, col]
        ax.plot(Sspot, r['sMTM'], 's-', color='C7', label='ΔMTM')
        ax.plot(Sspot, r['sCE'], 'd-', color='C1', label='ΔCE')
        ax.plot(Sspot, r['sEEPE_B'], 'o-', color='C0', label='ΔEEPE (B: realized vol fixed)')
        ax.plot(Sspot, r['sEEPE_A'], 'o-', color='C3', label='ΔEEPE (A: realized vol shocked)')
        ax.axhline(0, color='k', lw=.6); ax.axvline(0, color='k', lw=.6)
        ax.set(title=f'SPOT shock — {title} (baseMTM={r["baseMTM"]:.1f})', xlabel='spot shock s'); ax.grid(alpha=.3); ax.legend(fontsize=6.5)
        ax = axs[1, col]
        ax.plot(Div * 100, r['vMTM'], 's-', color='C7', label='ΔMTM')
        ax.plot(Div * 100, r['vCE'], 'd-', color='C1', label='ΔCE (full t0 vega)')
        ax.plot(Div * 100, r['vEEPE'], 'o-', color='C4', label='ΔEEPE (damped vega)')
        ax.axhline(0, color='k', lw=.6); ax.axvline(0, color='k', lw=.6)
        ax.set(title=f'IMPLIED-VOL shock — {title}', xlabel='implied vol shift (%)'); ax.grid(alpha=.3); ax.legend(fontsize=6.5)
    fig.suptitle('Part 2 — Δ metrics vs shock. Spot shock (top): delta tracks. Implied-vol shock (bottom): CE over-states (vega term structure)', fontsize=12)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "part2_tracking.png"), dpi=125); plt.close(fig)

    # Figure 3: scatter ΔEEPE vs ΔCE for spot and implied shocks
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.8))
    for ax, (key, title) in zip(axs, cases):
        r = res[key]
        def slope(x, y): return np.polyfit(x, y, 1)[0] if np.ptp(x) > 1e-6 else float('nan')
        ax.scatter(r['sCE'], r['sEEPE_B'], c='C0', s=20, label=f'spot shock (slope {slope(r["sCE"], r["sEEPE_B"]):.2f})')
        ax.scatter(r['vCE'], r['vEEPE'], c='C4', s=20, label=f'implied-vol shock (slope {slope(r["vCE"], r["vEEPE"]):.2f})')
        allx = np.concatenate([r['sCE'], r['vCE']])
        xs = np.array([allx.min(), allx.max()]); ax.plot(xs, xs, 'k--', lw=1, label='45°')
        ax.set(title=f'{title}: ΔEEPE vs ΔCE', xlabel='ΔCE (mm)', ylabel='ΔEEPE (mm)'); ax.grid(alpha=.3); ax.legend(fontsize=7)
    fig.suptitle('Part 2 — ΔEEPE vs ΔCE: spot (delta, ~1) vs implied-vol (vega, <1) shocks', fontsize=13)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "part2_scatter.png"), dpi=125); plt.close(fig)

    print("Part 2 summary:")
    for key, title in cases:
        r = res[key]
        def sl(x, y): return np.polyfit(x, y, 1)[0] if np.ptp(x) > 1e-6 else float('nan')
        print(f"  {title:12s} baseMTM={r['baseMTM']:6.2f} baseCE={r['baseCE']:6.2f} baseEEPE={r['baseEEPE']:6.2f}"
              f" | spot slope(EEPE/CE)={sl(r['sCE'], r['sEEPE_B']):.2f}"
              f"  implied slope(EEPE/CE)={sl(r['vCE'], r['vEEPE']):.2f}")


def fig_cloud():
    """Many different shocks (spot per factor, implied-vol, combined, random):
    spot family lands on the 45° line (delta), implied family falls below (vega)
    -> ΔCE is a cloud w.r.t. ΔEEPE once both dimensions are present."""
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    rng = np.random.default_rng(0)
    fams = ['equity spot', 'commodity spot', 'rate spot', 'implied vol (all)', 'combined spot', 'random spot+implied']
    col = dict(zip(fams, ['C0', 'C2', 'C8', 'C4', 'k', 'C7']))

    def shocks():
        L = []
        for s in np.linspace(-0.12, 0.12, 9):
            if abs(s) < 1e-9: continue
            L += [(Scenario(spot_mult={'eq': 1 + s}), 'equity spot'),
                  (Scenario(spot_mult={'cm': 1 + s}), 'commodity spot'),
                  (Scenario(rate_add={'rate': 0.008 * s}), 'rate spot'),
                  (sc_spot(s, False), 'combined spot')]
        for dv in np.linspace(-0.30, 0.30, 9):
            if abs(dv) < 1e-9: continue
            L.append((sc_iv(dv), 'implied vol (all)'))
        for _ in range(35):
            L.append((Scenario(spot_mult={'eq': 1 + rng.uniform(-.12, .12), 'cm': 1 + rng.uniform(-.12, .12)},
                               rate_add={'rate': rng.uniform(-.004, .004)},
                               iv_add={a: IMP0[a] * rng.uniform(-.25, .25) for a in IMP0}), 'random spot+implied'))
        return L

    cases = [('atm', 'ATM options'), ('otm', 'OTM options')]
    fig, axs = plt.subplots(2, 2, figsize=(13, 10))
    for ci, (key, title) in enumerate(cases):
        tr = book(key); base = metrics_for(tr, Scenario())
        data = {f: ([], [], []) for f in fams}          # dMTM, dCE, dEEPE
        for sc, fam in shocks():
            m = metrics_for(tr, sc)
            data[fam][0].append(m['MTM'] - base['MTM'])
            data[fam][1].append(m['CE'] - base['CE'])
            data[fam][2].append(m['EEPE'] - base['EEPE'])
        allE = np.concatenate([np.array(data[f][2]) for f in fams])
        for row, (xlab, idx) in enumerate([('ΔMTM', 0), ('ΔCE', 1)]):
            ax = axs[row, ci]; xall = np.concatenate([np.array(data[f][idx]) for f in fams])
            for f in fams:
                xs = np.array(data[f][idx])
                if len(xs): ax.scatter(xs, np.array(data[f][2]), s=15, color=col[f], alpha=0.8, label=f)
            a, b = np.polyfit(xall, allE, 1); resid = allE - (a * xall + b)
            lim = np.array([xall.min(), xall.max()]); ax.plot(lim, lim, 'k--', lw=1, label='45°')
            ax.set_title(f'{title}: ΔEEPE vs {xlab}  (slope {a:.2f}, σ = {resid.std():.2f} mm)')
            ax.set(xlabel=f'{xlab} (mm)', ylabel='ΔEEPE (mm)'); ax.grid(alpha=.3); ax.legend(fontsize=6)
    fig.suptitle('Part 2 — ΔEEPE vs ΔMTM (top) and ΔCE (bottom): spot tracks (45°), implied-vol falls below (vega)', fontsize=12)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "part2_cloud.png"), dpi=125); plt.close(fig)
    print("Part 2 cloud written.")


if __name__ == "__main__":
    fig_envelopes()
    fig_tracking()
    fig_cloud()
    print("figures written to", os.path.abspath(FIG))
