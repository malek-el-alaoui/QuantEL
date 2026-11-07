# Rich equity-chain volatility-smile comparison

This Visual Studio/C++ project applies the latest two-model equity methodology to public option chains with 22–30 retained market quotes per maturity.

## Models

**Model 1 — bounded log-normal mixture**
- one balanced boundary artificial quote;
- artificial quote calibration weight = 1e-10;
- 0.10 sigma_ATM <= sigma_j <= 3.00 sigma_ATM;
- component-weighted quadratic right-wing decay only outside the last market quote;
- direct static-arbitrage safety checks.

**Model 2 — natural cubic total-variance spline**
- x = log(K/F);
- w(x) = sigma_imp(x)^2 T;
- exact interpolation at retained market knots.

## Rich public data

- SPY 2026-10-16: 30 retained quotes
- TSLA 2026-10-16: 30 retained quotes
- NVDA 2026-10-16: 22 retained quotes
- AAPL 2026-10-16: 22 retained quotes

See `RICH_EQUITY_DATA_NOTES.md` for URLs, forward proxies and quote filtering.

## Build

Open `VolatilitySmileTwoModels.sln` in Visual Studio, choose x64/Release, rebuild, and run.

CMake is also supported:

```text
cmake -S . -B build
cmake --build build --config Release
```

Run all slices:

```text
VolatilitySmileTwoModels all
```

Regenerate plots:

```text
python python/plot_two_models.py
```

## Main outputs

- `output/all_market_vs_models.csv`
- `output/all_model1_components.csv`
- `output/all_model1_artificial_quote_diagnostics.csv`
- `output/all_continuous_smiles.csv`
- `output/all_densities.csv`
- `output/all_wing_extrapolation.csv`
- `output/summary.csv`
- `output/plots/`

See `RESULTS_SUMMARY.md` and `VALIDATION_RICH_EQUITY_CHAINS.txt`.
