# Results summary — rich equity chains

The project was rerun on four public 2026-10-16 option-chain slices with 22–30 retained quotes each.

| Slice | Quotes | Model 1 max IV error (bp) | Model 2 min density on 2-delta grid | Model 2 negative-density grid points | Max kernel / ATM | Kernels at 3x cap | Right-wing alpha | Model 1 min wing density |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| SPY_RICH | 30 | 3.062 | -8.17599 | 6 | 1.303x | 0 | 232.433 | 0.002486 |
| TSLA_RICH | 30 | 20.890 | -3.61461 | 25 | 2.584x | 0 | 14.1114 | 0.003801 |
| NVDA_RICH | 22 | 20.428 | -0.026178 | 6 | 1.086x | 0 | 77.5414 | 0.000687 |
| AAPL_RICH | 22 | 21.066 | -1.42553 | 117 | 3.000x | 1 | 3.53942e-12 | 0.027088 |

## Main observations

1. All four input slices pass the discrete static-arbitrage check before calibration.
2. Model 1 remains positive on the tested wing safety grids and passes the direct call monotonicity/convexity test for all four slices.
3. With many market knots, an exact natural cubic total-variance interpolant can oscillate between knots: all four dense-chain examples show at least one negative Model-2 density point on the tested grid.
4. SPY is the cleanest dense benchmark for Model 1: 30 retained quotes and about 3.1 bp maximum market-IV error.
5. TSLA, NVDA and AAPL are harder fits under the 10%-to-300%-ATM kernel bounds; their maximum quote errors are about 20–21 bp.
6. The right-wing control is applied only beyond the last retained market strike. If the observed chain itself turns upward on high strikes, the model must reproduce that market shape before the extrapolation rule begins.
