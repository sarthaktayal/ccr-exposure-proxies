"""
Part 1 -- Linear portfolio.
===========================
Mock linear multi-asset netting set (IR swap + FX/equity/commodity forwards).
For three MONEYNESS cases (deep-ITM, ATM, OTM) we shock the risk factors and
compare the change in EEPE against the change in CE and MTM, under TWO shocking
methodologies:
    (B) do NOT shock the simulation (realized) volatility  -> envelope width fixed
    (A) DO shock the simulation volatility with the scenario -> envelope width moves

Outputs (../figures):
    part1_sim_envelopes.png  -- how the simulated paths change after a shock, per
                                risk factor and for the portfolio value, A vs B.
    part1_tracking.png       -- dEEPE (A & B) vs dCE and dMTM across the shock, per
                                moneyness case.
    part1_scatter.png        -- dEEPE vs dCE scatter with fitted slopes.
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from ccrlib import Asset, Trade, Scenario, Simulator, price_portfolio, exposure_metrics

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(FIG, exist_ok=True)

# ----- risk factors (shared) -------------------------------------------------
ASSETS = [Asset('rate', 'ou', 0.03, 0.010, kappa=0.10, theta=0.03),
          Asset('eq', 'gbm', 100.0, 0.20, q=0.02),
          Asset('fx', 'gbm', 1.10, 0.10, q=0.01),
          Asset('cm', 'gbm', 80.0, 0.25, q=0.01)]
CORR = np.array([[1, -.2, .1, .05], [-.2, 1, .25, .30], [.1, .25, 1, .40], [.05, .30, .40, 1.]])
TIMES = np.array([k / 12 for k in range(61)])          # monthly, 5y
SIM = Simulator(ASSETS, CORR, TIMES, n_paths=24000, seed=1)

# ----- linear books at three moneyness levels --------------------------------
def book(moneyness):
    #                         swap fixed   fxK    eqK    cmK
    K = {'deep_itm': (0.018, 1.05, 88.0, 68.0),
         'atm':      (0.0305, 1.168, 103.0, 86.6),
         'otm':      (0.045, 1.30, 120.0, 103.0)}[moneyness]
    return [Trade('rate', 'swap', K[0], 100.0, 5.0, duration=4.5),
            Trade('fx', 'fwd', K[1], 50.0, 3.0),
            Trade('eq', 'fwd', K[2], 0.50, 3.0),
            Trade('cm', 'fwd', K[3], 0.40, 4.0)]

# ----- shock definition ------------------------------------------------------
BETA_VOL = 1.5   # methodology A couples realized-vol shock to the level shock

def scenario(s, shock_vol):
    return Scenario(
        spot_mult={'eq': 1 + 1.0 * s, 'fx': 1 + 0.5 * s, 'cm': 1 + 1.0 * s},
        rate_add={'rate': 0.008 * s},
        vol_mult=(1 + BETA_VOL * s) if shock_vol else 1.0,
    )

def metrics_for(trades, s, shock_vol):
    paths, ivp = SIM.simulate(scenario(s, shock_vol))
    V = price_portfolio(trades, ASSETS, paths, ivp, TIMES)
    return exposure_metrics(V, TIMES, eepe_window=1.0)

# ============================================================================ #
#  Figure 1 -- how the simulations change after a shock (per factor + portfolio)
# ============================================================================ #
def fig_envelopes():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    s_ill = 0.18
    base_p, _ = SIM.simulate(scenario(0.0, False))
    B_p, _ = SIM.simulate(scenario(s_ill, False))    # level only
    A_p, _ = SIM.simulate(scenario(s_ill, True))     # level + vol
    trades = book('deep_itm')
    base_V = price_portfolio(trades, ASSETS, base_p, {}, TIMES)
    B_V = price_portfolio(trades, ASSETS, B_p, {}, TIMES)
    A_V = price_portfolio(trades, ASSETS, A_p, {}, TIMES)

    def env(P):  # 5/50/95 percentile bands
        return np.percentile(P, 5, 1), np.percentile(P, 50, 1), np.percentile(P, 95, 1)

    fig, axs = plt.subplots(2, 3, figsize=(15, 8))
    panels = [('rate', 'Interest rate'), ('eq', 'Equity'), ('fx', 'FX'), ('cm', 'Commodity')]
    for ax, (nm, title) in zip(axs.ravel(), panels):
        for P, c, lab in [(base_p[nm], 'grey', 'base'), (B_p[nm], 'C0', '(B) level only'),
                          (A_p[nm], 'C3', '(A) level + vol')]:
            lo, md, hi = env(P)
            ax.fill_between(TIMES, lo, hi, color=c, alpha=0.12)
            ax.plot(TIMES, hi, color=c, lw=1.0); ax.plot(TIMES, lo, color=c, lw=1.0)
            ax.plot(TIMES, md, color=c, lw=1.8, label=f'{lab} (median)')   # median = the LEVEL
        ax.set(title=f'{title} — 5–95% envelope + median', xlabel='years'); ax.grid(alpha=.3)
        ax.legend(fontsize=7)
    # portfolio value envelope
    ax = axs.ravel()[4]
    for V, c, lab in [(base_V, 'grey', 'base'), (B_V, 'C0', '(B) level only'), (A_V, 'C3', '(A) level+vol')]:
        lo, md, hi = env(V)
        ax.fill_between(TIMES, lo, hi, color=c, alpha=0.12)
        ax.plot(TIMES, hi, color=c, lw=1.0); ax.plot(TIMES, lo, color=c, lw=1.0)
        ax.plot(TIMES, md, color=c, lw=1.8, label=f'{lab} (median)')
    ax.set(title='Portfolio value V (deep-ITM book)', xlabel='years'); ax.grid(alpha=.3); ax.legend(fontsize=7)
    axs.ravel()[5].axis('off')
    axs.ravel()[5].text(0.02, 0.5,
        "Shock s = +0.18.\n\nThe LEVEL shock is IDENTICAL under A\nand B: the median lines (bold) of A (red)\n"
        "and B (blue) COINCIDE.\n\n(B) leaves the envelope width unchanged.\n\n"
        "(A) keeps the same level but WIDENS\n     the envelope (realized-vol up).\n\n"
        "CE and MTM (t=0 marks) are identical\nunder A and B — blind to the width\n"
        "change. Only EEPE sees it.", fontsize=10, va='center')
    fig.suptitle('Part 1 — How the simulated risk factors change after a shock (A: vol shocked vs B: not)',
                 fontsize=13)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "part1_sim_envelopes.png"), dpi=125)
    plt.close(fig)

# ============================================================================ #
#  Figures 2 & 3 -- tracking of dEEPE vs dCE / dMTM across moneyness & method
# ============================================================================ #
def fig_tracking():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    cases = [('deep_itm', 'Deep ITM'), ('atm', 'ATM'), ('otm', 'OTM')]
    S = np.linspace(-0.12, 0.12, 13)
    results = {}
    for key, _ in cases:
        tr = book(key)
        base = metrics_for(tr, 0.0, False)
        rows = {'s': S, 'dMTM': [], 'dCE': [], 'dEEPE_B': [], 'dEEPE_A': []}
        for s in S:
            mB = metrics_for(tr, s, False); mA = metrics_for(tr, s, True)
            rows['dMTM'].append(mB['MTM'] - base['MTM'])
            rows['dCE'].append(mB['CE'] - base['CE'])
            rows['dEEPE_B'].append(mB['EEPE'] - base['EEPE'])
            rows['dEEPE_A'].append(mA['EEPE'] - base['EEPE'])
        results[key] = {k: np.array(v) for k, v in rows.items()}
        results[key]['baseMTM'] = base['MTM']

    # Figure 2: metrics vs shock, per moneyness
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.8))
    for ax, (key, title) in zip(axs, cases):
        r = results[key]
        ax.plot(S, r['dMTM'], 's-', color='C7', label='ΔMTM')
        ax.plot(S, r['dCE'], 'd-', color='C1', label='ΔCE')
        ax.plot(S, r['dEEPE_B'], 'o-', color='C0', label='ΔEEPE (B: vol fixed)')
        ax.plot(S, r['dEEPE_A'], 'o-', color='C3', label='ΔEEPE (A: vol shocked)')
        ax.axhline(0, color='k', lw=.6); ax.axvline(0, color='k', lw=.6)
        ax.set(title=f'{title}  (base MTM={r["baseMTM"]:.1f})', xlabel='risk-factor shock s'); ax.grid(alpha=.3)
        ax.legend(fontsize=7)
    fig.suptitle('Part 1 — Δ metrics vs shock: ΔEEPE (A/B) compared with ΔCE and ΔMTM', fontsize=13)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "part1_tracking.png"), dpi=125); plt.close(fig)

    # Figure 3: dEEPE vs dMTM (top) and vs dCE (bottom), per moneyness
    fig, axs = plt.subplots(2, 3, figsize=(15, 9))
    for col, (key, title) in enumerate(cases):
        r = results[key]
        for row, (xk, xlab) in enumerate([('dMTM', 'ΔMTM'), ('dCE', 'ΔCE')]):
            ax = axs[row, col]; x = r[xk]
            for dee, c, lab in [(r['dEEPE_B'], 'C0', 'B: vol fixed'), (r['dEEPE_A'], 'C3', 'A: vol shocked')]:
                m = np.polyfit(x, dee, 1)[0] if np.ptp(x) > 1e-6 else np.nan
                ax.scatter(x, dee, c=c, s=22, label=f'{lab} (slope {m:.2f})')
            if np.ptp(x) > 1e-6:
                xs = np.array([x.min(), x.max()]); ax.plot(xs, xs, 'k--', lw=1, label='45°')
            else:
                ax.axvline(x[0], color='C3', ls=':', label='ΔCE frozen (OTM)')
            ax.set(title=f'{title}: ΔEEPE vs {xlab}', xlabel=f'{xlab} (mm)', ylabel='ΔEEPE (mm)')
            ax.grid(alpha=.3); ax.legend(fontsize=7)
    fig.suptitle('Part 1 — ΔEEPE vs ΔMTM (top) and ΔCE (bottom): tracking by moneyness & method', fontsize=13)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "part1_scatter.png"), dpi=125); plt.close(fig)

    # console summary
    print("Part 1 summary (slope dEEPE/dCE):")
    for key, title in cases:
        r = results[key]
        def slope(dee):
            return np.polyfit(r['dCE'], dee, 1)[0] if np.ptp(r['dCE']) > 1e-6 else float('nan')
        print(f"  {title:9s} baseMTM={r['baseMTM']:7.2f} | B(vol fixed)={slope(r['dEEPE_B']):.3f}"
              f"  A(vol shocked)={slope(r['dEEPE_A']):.3f}")


if __name__ == "__main__":
    fig_envelopes()
    fig_tracking()
    print("figures written to", os.path.abspath(FIG))
