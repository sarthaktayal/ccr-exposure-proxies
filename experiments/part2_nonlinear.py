"""
Part 2 -- Non-linear portfolio (options across asset classes).
==============================================================
Long options on interest-rate, equity, FX and commodity underlyings. IMPLIED
VOLATILITY is a SEPARATE stochastic RISK FACTOR (mean-reverting log-OU) used to
reprice the options -- distinct from the REALIZED volatility that diffuses the
underlyings. As a risk factor, implied vol has BOTH a LEVEL (the base implied vol)
and a DIFFUSION vol (its vol-of-vol).

Shock methodology (mirrors Part 1, now with implied vol included):
  A LEVEL shock moves ALL risk-factor levels -- spot levels AND the base implied
  vol level. Two cases:
    (B) do NOT shock any DIFFUSION vol  -> realized spot vol fixed AND vol-of-vol fixed
    (A) DO shock the diffusion vols     -> realized spot vol AND vol-of-vol shocked
  CE and MTM (t0 marks) are IDENTICAL under A and B -- they see the level shock
  (spot delta + implied-vol vega) but are blind to the diffusion/width change.
  Only EEPE responds to (A) vs (B).

Outputs (../figures):
    part2_sim_envelopes.png  -- spot AND implied-vol envelopes + portfolio, base / B / A.
    part2_tracking.png       -- Δ metrics vs shock, per moneyness (ΔMTM, ΔCE, ΔEEPE B/A).
    part2_scatter.png        -- ΔEEPE vs ΔMTM (top) and ΔCE (bottom), B and A.
    part2_cloud.png          -- many different shocks: line vs cloud.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from ccrlib import Asset, Trade, Scenario, Simulator, price_portfolio, exposure_metrics, R_DISC

FIG = os.path.join(os.path.dirname(__file__), "..", "figures"); os.makedirs(FIG, exist_ok=True)

ASSETS = [
    Asset('rate', 'ou', 0.03, 0.010, kappa=0.10, theta=0.03, imp0=0.010, iv_kappa=3.0, iv_nu=0.6, iv_rho=-0.1),
    Asset('eq', 'gbm', 100.0, 0.20, q=0.02, imp0=0.20, iv_kappa=3.0, iv_nu=0.8, iv_rho=-0.7),
    Asset('fx', 'gbm', 1.10, 0.10, q=0.01, imp0=0.10, iv_kappa=3.0, iv_nu=0.5, iv_rho=-0.2),
    Asset('cm', 'gbm', 80.0, 0.25, q=0.01, imp0=0.25, iv_kappa=3.0, iv_nu=0.6, iv_rho=-0.3),
]
IMP0 = {a.name: a.imp0 for a in ASSETS}
CORR = np.array([[1, -.2, .1, .05], [-.2, 1, .25, .30], [.1, .25, 1, .40], [.05, .30, .40, 1.]])
TIMES = np.array([k / 12 for k in range(13)])          # monthly, 1y
SIM = Simulator(ASSETS, CORR, TIMES, n_paths=24000, seed=2)
T = 1.0
FWD = {'eq': 100 * np.exp((R_DISC - 0.02) * T), 'fx': 1.10 * np.exp((R_DISC - 0.01) * T),
       'cm': 80 * np.exp((R_DISC - 0.01) * T), 'rate': 0.03}

def book(moneyness):
    mult = {'deep_itm': 0.85, 'atm': 1.00, 'otm': 1.15}[moneyness]
    return [Trade('rate', 'rate_call', FWD['rate'] * mult, 100.0, T, duration=4.5),
            Trade('eq', 'call', FWD['eq'] * mult, 0.50, T),
            Trade('fx', 'call', FWD['fx'] * mult, 50.0, T),
            Trade('cm', 'call', FWD['cm'] * mult, 0.30, T)]

BETA_VOL = 1.5        # diffusion-vol shock coupling (methodology A): realized vol AND vol-of-vol
IV_LEVEL_BETA = 1.5   # base implied-vol LEVEL moves with the level shock

def sc_all(s, shock_diff_vol):
    """Shock ALL levels (spot + base implied vol). If shock_diff_vol, ALSO shock all
    diffusion vols (realized spot vol AND implied-vol's vol-of-vol)."""
    m = (1 + BETA_VOL * s) if shock_diff_vol else 1.0
    return Scenario(spot_mult={'eq': 1 + s, 'fx': 1 + 0.5 * s, 'cm': 1 + s},
                    rate_add={'rate': 0.008 * s},
                    iv_add={a: IMP0[a] * (IV_LEVEL_BETA * s) for a in IMP0},
                    vol_mult=m, iv_vov_mult=m)

def metrics_for(trades, sc):
    paths, ivp = SIM.simulate(sc)
    V = price_portfolio(trades, ASSETS, paths, ivp, TIMES)
    return exposure_metrics(V, TIMES, eepe_window=1.0)

# ============================================================================ #
#  Figure 1 -- simulations after shock: spots + implied vol, base / B / A
# ============================================================================ #
def fig_envelopes():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    s_ill = 0.15
    base, _ = SIM.simulate(Scenario())
    Bp, Biv = SIM.simulate(sc_all(s_ill, False))     # levels (incl implied) shocked; diffusion vols FIXED
    Ap, Aiv = SIM.simulate(sc_all(s_ill, True))      # + diffusion vols shocked (realized & vol-of-vol)
    _, base_iv = SIM.simulate(Scenario())
    trades = book('atm')
    Vb = price_portfolio(trades, ASSETS, base, base_iv, TIMES)
    VB = price_portfolio(trades, ASSETS, Bp, Biv, TIMES)
    VA = price_portfolio(trades, ASSETS, Ap, Aiv, TIMES)
    STY = {'base': ('grey', '-', 1.8, 7), 'B': ('C0', '-', 2.6, 11), 'A': ('C3', '--', 1.6, 5)}

    def draw(ax, series, scale=1.0):
        for P, key, lab in series:
            c, ls, lw, ms = STY[key]; lo, md, hi = (np.percentile(P, q, 1) * scale for q in (5, 50, 95))
            ax.fill_between(TIMES, lo, hi, color=c, alpha=0.10)
            ax.plot(TIMES, hi, color=c, lw=0.9, ls=ls); ax.plot(TIMES, lo, color=c, lw=0.9, ls=ls)
            ax.plot(TIMES, md, color=c, lw=lw, ls=ls, label=f'{lab} (median)')
            ax.plot(0, md[0], 'o', color=c, ms=ms, zorder=5)
        ax.grid(alpha=.3); ax.legend(fontsize=6.5)

    fig, axs = plt.subplots(2, 3, figsize=(15, 8.5))
    draw(axs[0, 0], [(base['eq'], 'base', 'base'), (Bp['eq'], 'B', '(B) levels'), (Ap['eq'], 'A', '(A) levels+diff vol')])
    axs[0, 0].set(title='Equity SPOT', xlabel='years')
    draw(axs[0, 1], [(base_iv['eq'], 'base', 'base'), (Biv['eq'], 'B', '(B) levels'), (Aiv['eq'], 'A', '(A) levels+vov')], scale=100)
    axs[0, 1].set(title='Equity IMPLIED VOL (%) — level shocked in B & A; vol-of-vol only in A', ylabel='vol %', xlabel='years')
    draw(axs[0, 2], [(base['cm'], 'base', 'base'), (Bp['cm'], 'B', '(B) levels'), (Ap['cm'], 'A', '(A) levels+diff vol')])
    axs[0, 2].set(title='Commodity SPOT', xlabel='years')
    draw(axs[1, 0], [(base_iv['cm'], 'base', 'base'), (Biv['cm'], 'B', '(B) levels'), (Aiv['cm'], 'A', '(A) levels+vov')], scale=100)
    axs[1, 0].set(title='Commodity IMPLIED VOL (%)', ylabel='vol %', xlabel='years')
    draw(axs[1, 1], [(base['rate'], 'base', 'base'), (Bp['rate'], 'B', '(B) levels'), (Ap['rate'], 'A', '(A) levels+diff vol')], scale=100)
    axs[1, 1].set(title='Interest rate (%)', xlabel='years')
    draw(axs[1, 2], [(Vb, 'base', 'base'), (VB, 'B', '(B) levels'), (VA, 'A', '(A) levels+diff vol')])
    axs[1, 2].set(title='Portfolio value V (ATM options)', xlabel='years')
    fig.suptitle('Part 2 — All levels shocked (spot + base implied vol). B: diffusion vols fixed. A: diffusion vols (realized & vol-of-vol) also shocked.', fontsize=11.5)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "part2_sim_envelopes.png"), dpi=125); plt.close(fig)

# ============================================================================ #
#  Figures 2 & 3 -- tracking and scatter under the unified shock
# ============================================================================ #
def fig_tracking():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    cases = [('deep_itm', 'ITM options'), ('atm', 'ATM options'), ('otm', 'OTM options')]
    S = np.linspace(-0.12, 0.12, 13); res = {}
    for key, _ in cases:
        tr = book(key); base = metrics_for(tr, Scenario())
        r = {'dMTM': [], 'dCE': [], 'dEEPE_B': [], 'dEEPE_A': []}
        for s in S:
            mB = metrics_for(tr, sc_all(s, False)); mA = metrics_for(tr, sc_all(s, True))
            r['dMTM'].append(mB['MTM'] - base['MTM']); r['dCE'].append(mB['CE'] - base['CE'])
            r['dEEPE_B'].append(mB['EEPE'] - base['EEPE']); r['dEEPE_A'].append(mA['EEPE'] - base['EEPE'])
        res[key] = {k: np.array(v) for k, v in r.items()}; res[key]['baseMTM'] = base['MTM']

    fig, axs = plt.subplots(1, 3, figsize=(15, 4.8))
    for ax, (key, title) in zip(axs, cases):
        r = res[key]
        ax.plot(S, r['dMTM'], 's-', color='C7', label='ΔMTM')
        ax.plot(S, r['dCE'], 'd-', color='C1', label='ΔCE')
        ax.plot(S, r['dEEPE_B'], 'o-', color='C0', label='ΔEEPE (B: diffusion vols fixed)')
        ax.plot(S, r['dEEPE_A'], 'o-', color='C3', label='ΔEEPE (A: diffusion vols shocked)')
        ax.axhline(0, color='k', lw=.6); ax.axvline(0, color='k', lw=.6)
        ax.set(title=f'{title} (baseMTM={r["baseMTM"]:.1f})', xlabel='level shock s'); ax.grid(alpha=.3); ax.legend(fontsize=6.5)
    fig.suptitle('Part 2 — Δ metrics vs a shock to ALL levels (spot + base implied vol). CE=MTM here (long options); A vs B differ only by the diffusion-vol shock.', fontsize=11)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "part2_tracking.png"), dpi=125); plt.close(fig)

    fig, axs = plt.subplots(2, 3, figsize=(15, 9))
    for col, (key, title) in enumerate(cases):
        r = res[key]
        for row, (xk, xlab) in enumerate([('dMTM', 'ΔMTM'), ('dCE', 'ΔCE')]):
            ax = axs[row, col]; x = r[xk]
            for dee, c, lab in [(r['dEEPE_B'], 'C0', 'B: diff vols fixed'), (r['dEEPE_A'], 'C3', 'A: diff vols shocked')]:
                sl = np.polyfit(x, dee, 1)[0] if np.ptp(x) > 1e-6 else float('nan')
                ax.scatter(x, dee, c=c, s=20, label=f'{lab} (slope {sl:.2f})')
            if np.ptp(x) > 1e-6:
                xs = np.array([x.min(), x.max()]); ax.plot(xs, xs, 'k--', lw=1, label='45°')
            ax.set(title=f'{title}: ΔEEPE vs {xlab}', xlabel=f'{xlab} (mm)', ylabel='ΔEEPE (mm)'); ax.grid(alpha=.3); ax.legend(fontsize=6.5)
    fig.suptitle('Part 2 — ΔEEPE vs ΔMTM (top) and ΔCE (bottom): B vs A (diffusion-vol shock)', fontsize=12)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "part2_scatter.png"), dpi=125); plt.close(fig)

    print("Part 2 summary (slope dEEPE/dCE):")
    for key, title in cases:
        r = res[key]
        def sl(y): return np.polyfit(r['dCE'], y, 1)[0] if np.ptp(r['dCE']) > 1e-6 else float('nan')
        print(f"  {title:12s} baseMTM={r['baseMTM']:6.2f} | B={sl(r['dEEPE_B']):.3f}  A={sl(r['dEEPE_A']):.3f}")

# ============================================================================ #
#  Figure 4 -- many different shocks: line or cloud?
# ============================================================================ #
def fig_cloud():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    rng = np.random.default_rng(0)
    fams = ['spot levels only', 'implied-vol level only', 'all levels (B)', 'all levels + diff vols (A)', 'random']
    COL = dict(zip(fams, ['C0', 'C4', 'k', 'C3', 'C7']))

    def shocks():
        L = []
        for s in np.linspace(-0.12, 0.12, 9):
            if abs(s) < 1e-9: continue
            L += [(Scenario(spot_mult={'eq': 1 + s, 'fx': 1 + 0.5 * s, 'cm': 1 + s}, rate_add={'rate': 0.008 * s}), 'spot levels only'),
                  (Scenario(iv_add={a: IMP0[a] * (IV_LEVEL_BETA * s) for a in IMP0}), 'implied-vol level only'),
                  (sc_all(s, False), 'all levels (B)'),
                  (sc_all(s, True), 'all levels + diff vols (A)')]
        for _ in range(35):
            s = rng.uniform(-.12, .12); dv = rng.uniform(0.6, 1.6)
            L.append((Scenario(spot_mult={'eq': 1 + rng.uniform(-.12, .12), 'cm': 1 + rng.uniform(-.12, .12)},
                               rate_add={'rate': rng.uniform(-.004, .004)},
                               iv_add={a: IMP0[a] * rng.uniform(-.25, .25) for a in IMP0},
                               vol_mult=dv, iv_vov_mult=dv), 'random'))
        return L

    cases = [('atm', 'ATM options'), ('otm', 'OTM options')]
    fig, axs = plt.subplots(2, 2, figsize=(13, 10))
    for ci, (key, title) in enumerate(cases):
        tr = book(key); base = metrics_for(tr, Scenario())
        data = {f: ([], [], []) for f in fams}
        for sc, fam in shocks():
            m = metrics_for(tr, sc)
            data[fam][0].append(m['MTM'] - base['MTM']); data[fam][1].append(m['CE'] - base['CE']); data[fam][2].append(m['EEPE'] - base['EEPE'])
        allE = np.concatenate([np.array(data[f][2]) for f in fams])
        for row, (xlab, idx) in enumerate([('ΔMTM', 0), ('ΔCE', 1)]):
            ax = axs[row, ci]; xall = np.concatenate([np.array(data[f][idx]) for f in fams])
            for f in fams:
                xs = np.array(data[f][idx])
                if len(xs): ax.scatter(xs, np.array(data[f][2]), s=15, color=COL[f], alpha=0.8, label=f)
            a, b = np.polyfit(xall, allE, 1); resid = allE - (a * xall + b)
            lim = np.array([xall.min(), xall.max()]); ax.plot(lim, lim, 'k--', lw=1, label='45°')
            ax.set_title(f'{title}: ΔEEPE vs {xlab}  (slope {a:.2f}, σ = {resid.std():.2f} mm)')
            ax.set(xlabel=f'{xlab} (mm)', ylabel='ΔEEPE (mm)'); ax.grid(alpha=.3); ax.legend(fontsize=6)
    fig.suptitle('Part 2 — many shocks: implied-vol level (purple) and diffusion-vol (A) shocks pull ΔEEPE off the ΔCE/ΔMTM line', fontsize=11.5)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "part2_cloud.png"), dpi=125); plt.close(fig)
    print("Part 2 cloud written.")


if __name__ == "__main__":
    fig_envelopes()
    fig_tracking()
    fig_cloud()
    print("figures written to", os.path.abspath(FIG))
