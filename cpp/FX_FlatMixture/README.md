# VolatilitySmileTwoModels — balanced one-point boundary

This Visual Studio/C++ project compares:

1. **Model 1** — arbitrage-free log-normal mixture with exactly one boundary
   artificial quote, chosen to balance/minimize the two neighboring butterfly
   masses, tiny artificial calibration weight `1e-10`, and
   `sigma_j > 10% * sigma_ATM`.
2. **Model 2** — natural cubic interpolation of total implied variance versus
   log-moneyness.

All prior recursive/midpoint artificial-quote experiments are removed.

## Build

Open `VolatilitySmileTwoModels.sln`, choose `Release | x64`, rebuild, and run.
A portable CMake build is also included.

## Bundled datasets

Existing tests: EURO STOXX 50, SPX, EUR/USD.

Additional equity tests: **AAPL, NVDA and TSLA** from small Cboe delayed-chain
snapshots.  See `data/README.md` for the exact normalization convention and
provenance.

## New-equity numerical results

| Dataset | Maturity | points | Model 1 max IV error (bp) | Model 2 min density | boundary | balanced boundary weight |
|---|---|---:|---:|---:|---|---:|
| AAPL | 2026-09-11 | 5 | 0.0076 | 2.6977 | right | 0.139170 |
| NVDA | 2026-09-14 | 6 | 0.0000 | 1.4587 | left | 0.172353 |
| TSLA | 2026-09-14 | 6 | 0.4935 | 1.3697 | right | 0.135729 |

For every tested slice, the artificial butterfly mass and its adjacent market
butterfly mass agree to numerical precision, as intended by the minimax/L2
construction.

## Main outputs

- `output/all_market_vs_models.csv`
- `output/all_model1_components.csv`
- `output/all_model1_artificial_quote_diagnostics.csv`
- `output/all_continuous_smiles.csv`
- `output/all_densities.csv`
- `output/all_multi_maturity_comparison.csv`
- `output/summary.csv`
- `VALIDATION_BALANCED_BOUNDARY.txt`
- `BALANCED_BOUNDARY_NOTES.md`

## Plots

Run `python python/plot_two_models.py` to regenerate:

- `output/plots/market_vs_models/`
- `output/plots/continuous_smiles_delta_range/`
- `output/plots/densities_delta_range/`
- `output/plots/model1_component_kernels/`

The component plots show the individual weighted log-normal kernels and mark
the artificial boundary component with a dashed line.
