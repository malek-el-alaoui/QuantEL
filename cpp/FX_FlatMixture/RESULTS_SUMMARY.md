# Results summary — balanced one-point mixture vs cubic variance spline

## Core checks

- Exactly one boundary artificial quote per slice (or none if already balanced).
- Artificial quote calibration objective weight: `1e-10`.
- Model-1 kernel-vol floor: `sigma_j > 0.10 * sigma_ATM`.
- Boundary/artificial butterfly and adjacent market butterfly are equalized whenever the unconstrained equal split is feasible.
- All 18 bundled slices are discretely arbitrage-free on input.
- Model 1 density is nonnegative on every generated output grid point.

## All slices

| Dataset | Maturity | N | M1 max IV err bp | M2 min density | M2 negative? | boundary | w_boundary | w_neighbor |
|---|---|---:|---:|---:|---|---|---:|---:|
| eurostoxx50 | 2014-10-17 | 10 | 0.0000 | -0.11507 | yes | left | 0.0171258 | 0.0171258 |
| eurostoxx50 | 2014-11-21 | 11 | 0.0000 | 0.323743 | no | left | 0.0522012 | 0.0522012 |
| eurostoxx50 | 2014-12-19 | 11 | 0.4364 | 0.373318 | no | left | 0.0766648 | 0.0766648 |
| eurostoxx50 | 2015-03-20 | 11 | 1.0083 | 0.00838348 | no | left | 0.119029 | 0.119029 |
| eurostoxx50 | 2015-09-18 | 11 | 3.7657 | 0.021444 | no | left | 0.185578 | 0.185578 |
| spx | 2026-07-31 | 5 | 24.6468 | -44.397 | yes | left | 0.192466 | 0.192466 |
| eurusd | 1WK | 5 | 0.0000 | 5.71677 | no | right | 0.0751951 | 0.0751951 |
| eurusd | 1MO | 5 | 0.0000 | 2.56624 | no | left | 0.0790668 | 0.0790668 |
| eurusd | 2MO | 5 | 0.0000 | 1.85371 | no | left | 0.0753998 | 0.0753998 |
| eurusd | 3MO | 5 | 0.0000 | 1.41353 | no | left | 0.0732375 | 0.0732375 |
| eurusd | 6MO | 5 | 0.0000 | 0.887839 | no | left | 0.0716042 | 0.0716042 |
| eurusd | 9MO | 5 | 0.0000 | 0.682597 | no | left | 0.0690993 | 0.0690993 |
| eurusd | 1YR | 5 | 0.0000 | 0.535047 | no | left | 0.0670375 | 0.0670375 |
| eurusd | 18MO | 5 | 0.0000 | 0.389745 | no | left | 0.069933 | 0.069933 |
| eurusd | 2YR | 5 | 0.0000 | 0.302397 | no | left | 0.0707234 | 0.0707234 |
| aapl | 2026-09-11 | 5 | 0.0076 | 2.69772 | no | right | 0.13917 | 0.13917 |
| nvda | 2026-09-14 | 6 | 0.0000 | 1.4587 | no | left | 0.172353 | 0.172353 |
| tsla | 2026-09-14 | 6 | 0.4935 | 1.36974 | no | right | 0.135729 | 0.135729 |

## Interpretation

The balanced boundary construction minimizes the largest of the two local endpoint butterfly masses because their sum is fixed by the neighboring market slope.  It therefore removes the arbitrary choice of a very small boundary mass and replaces it with the unique minimax/L2 split.

It does **not** guarantee a hump-free continuous mixture density by itself: the calibrated kernel widths also matter.  This is why the project includes individual component-kernel plots.  In the added very short-dated equity slices, Model 1 remains positive and fits closely, but narrow peaks can still be visible.

The cubic variance spline interpolates the market IV knots exactly, but it is not constrained to be butterfly-arbitrage-free; negative spline density remains visible where it occurs.
