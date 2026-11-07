# Imaginary Local Volatility — C++ project

Repackaging of **two previously developed C++ numerical experiments**, not a new implementation of the full November 2026 article.

## Modules

- `src/dual_arbitrage_margin.cpp` — numerical arbitrage-margin experiment from `Dual_Dupire_Arbitrage_Margin_LaTeX_Project_v5.zip`: synthetic volatility shocks, beta margin, positive/regularized local-variance PDEs, vanilla repricing and autocall finite-bump vega stability. Generates `autocall_vega_step_stability_cpp.csv` in its *current working directory*.
- `src/extended_dupire_comparison.cpp` — earlier `extended_dupire_classic_vs_extended.zip` benchmark for synthetic EURO STOXX-like arbitrageable surfaces, calibration and validation against repaired classical Dupire. Its original README is also provided for the data/output conventions.

The C++ source of both modules has been preserved as recovered from their original packages. The combined build wrapper is new.

## Compile (C++17, CMake)

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

On Windows with Visual Studio, build with `cmake -S . -B build` followed by `cmake --build build --config Release`.

## Run

```bash
# From the project root, so CSV output is created here:
./build/dual_arbitrage_margin
./build/extended_dupire_comparison
```

Both experiments may take time due to finite-difference grids and calibration.

## Scope and limitations

These are exploratory numerical prototypes previously used for the dual-Dupire research and synthetic benchmark. They **do not** constitute a full implementation of the manuscript's complex payoff/measure, uniqueness of Jordan decomposition, or the full continuously recalibrated dynamic delta-hedging scheme.

Synthetic outputs and diagnostics are demonstration results, not market-data verification. See the original extended-Dupire README for experiment assumptions.
