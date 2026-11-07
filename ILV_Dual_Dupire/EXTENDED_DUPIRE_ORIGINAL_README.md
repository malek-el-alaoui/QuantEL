# Classical vs Extended Dupire on Arbitrage-Contaminated EURO STOXX-like Surfaces

This project compares:

1. **Classical Dupire**: raw finite-difference local variance from the market surface. Because the synthetic surface deliberately violates no-arbitrage, some raw variances are negative/singular. To obtain a genuine real diffusion, the executable benchmark clips variance to `[1e-6, 1.0]` and reprices with a fully implicit forward Dupire PDE.
2. **Extended two-clock Dupire**: the market surface is represented as
   `C = (1+beta) C_plus - beta C_minus`, where `C_plus` and `C_minus` each come from a positive-variance diffusion. The signed projection reproduces the inconsistent market surface.

## Main numerical results

### Moderate arbitrage surface
- 184 quotes
- 34 non-positive/invalid market-gamma nodes
- 46 calendar violations
- Extended RMSE: `8.29e-12`
- Classical admissible Dupire RMSE: `5.8018`
- Classical max absolute error: `18.9926`

### Violent butterfly-arbitrage surface
- 164 quotes
- worst discrete gamma: about `-3.604e-3`
- Extended RMSE: `2.28e-12`
- Classical admissible Dupire RMSE: `560.022`
- Classical max absolute error: `1495.505`

The violent case is intentionally synthetic and extreme. It is included to make the structural limitation of a positive-density one-factor Dupire diffusion visually obvious.

## Build and run

Linux/macOS:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/extended_dupire
python3 plot_results.py
```

Windows with Visual Studio / CMake:

```bat
cmake -S . -B build
cmake --build build --config Release
build\Release\extended_dupire.exe
python plot_results.py
```

## Folder layout

- `src/main.cpp` - complete C++17 implementation
- `data/` - moderate and violent synthetic market surfaces
- `output/regular/` - classical and extended diagnostics for moderate arbitrage
- `output/violent/` - classical and extended diagnostics for violent arbitrage
- `plots/regular/` - separate price, gamma, and variance graphs for every maturity
- `plots/violent/` - separate price, gamma, and variance graphs for every maturity
- `paper/extended_dupire_comparison.pdf` - research-paper-style writeup
- `paper/extended_dupire_comparison.tex` - LaTeX source
- `plot_results.py` - graph generation

## Important modeling caveat

The extended signed pricing functional is **not** a single positive risk-neutral expectation and therefore is not an arbitrage-free model for a traded underlying. It is a representation of an inconsistent/contaminated market quote surface. The two hidden branches themselves remain ordinary positive-variance diffusions.
