# CCR Exposure Proxies — when does ΔCE track ΔEEPE?

A small, reproducible research repo that studies **how well the change in Current Exposure (ΔCE) and Mark-to-Market (ΔMTM) proxy the change in Effective EPE (ΔEEPE)** for counterparty-credit portfolios under market shocks — across moneyness regimes, for linear and non-linear (options) books, and treating **implied volatility as a risk factor distinct from the simulation volatility**.

> **Read the full write-up:** [`article/article.md`](article/article.md)

## TL;DR findings

| Regime / shock | ΔCE vs ΔEEPE | Why |
|---|---|---|
| Linear, deep-ITM, level shock | tracks (slope ≈ 1) | max inactive, no vega |
| Linear, ATM / OTM | degrades → **fails** | hockey-stick kink; CE floors at 0 |
| **Shock the simulation volatility** | **fails** (CE flat) | linear book has no vega |
| Options, **spot** shock | tracks (slope ≈ 1) | CE carries delta/gamma |
| Options, **implied-vol** shock | **over-states** (slope ≈ 2/3) | vega term structure |

## Repository layout

```
ccr-exposure-proxies/
├── ccrlib/                 # the Monte-Carlo engine
│   └── core.py             #   risk-factor sim, pricing, EE/EEPE/CE/MTM metrics
├── experiments/
│   ├── part1_linear.py     # linear book: moneyness × {shock vol / don't}
│   └── part2_nonlinear.py  # options: implied vol as a separate risk factor
├── article/article.md      # the academic write-up (methodology + figures)
├── figures/                # generated PNGs
├── run_all.py              # regenerate every figure
└── requirements.txt
```

## Quick start

```bash
pip install -r requirements.txt
python run_all.py            # regenerates all figures into figures/
```

or run the two studies individually:

```bash
python experiments/part1_linear.py
python experiments/part2_nonlinear.py
```

## What each part does

**Part 1 — Linear portfolio.** An IR swap + equity/FX/commodity forwards, in three moneyness cases (deep-ITM / ATM / OTM). We shock the risk factors two ways — **(B)** holding the simulation volatility fixed, and **(A)** shocking it — and compare ΔEEPE against ΔCE and ΔMTM. Figures show how the *simulated paths* change under each methodology (per risk factor and for the portfolio), and how the tracking degrades with moneyness.

**Part 2 — Non-linear portfolio.** Long options on rate/equity/FX/commodity underlyings. **Implied volatility is a stochastic risk factor** (mean-reverting log-OU with a leverage correlation) kept *separate* from the realized volatility that diffuses the underlyings. We shock spot and implied vol independently and show that ΔCE tracks a spot (delta) shock but **over-states an implied-vol (vega) shock** by the vega-term-structure factor of ≈ 2/3.

## Key figures

| | |
|---|---|
| Linear — sims change under shock | `figures/part1_sim_envelopes.png` |
| Linear — ΔEEPE vs ΔCE/ΔMTM by moneyness | `figures/part1_tracking.png`, `part1_scatter.png` |
| Non-linear — implied vol as a risk factor | `figures/part2_sim_envelopes.png` |
| Non-linear — spot vs implied-vol shocks | `figures/part2_tracking.png`, `part2_scatter.png` |

## Caveats

This is a **teaching / research** engine, not a booking system: single netting set, flat discount curve, non-discounted EE, parallel implied-vol surface (no skew/term-structure dynamics). See the article's discussion for how each simplification would change the conclusions.

## License

MIT.
