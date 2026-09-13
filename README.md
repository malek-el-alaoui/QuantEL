# Log-normal mixture article — full revision package

This package contains the revised LaTeX article **and the C++/Visual Studio implementations used for the numerical experiments**.

## Folders

- `article/` — LaTeX source, bibliography embedded in `main.tex`, generated figures, and the compiled revised PDF.
- `cpp/FX_FlatMixture/` — FX implementation for EUR/USD and USD/HKD. The mixture kernel widths remain flat outside the quoted region; **no equity right-wing control is applied to FX**.
- `cpp/Equity_RightWing/` — dense equity/index implementation (SPY, TSLA, NVDA, AAPL, EURO STOXX 50) with the bounded-kernel calibration and the optional right-wing control used in the article's extrapolation study.

Both C++ projects contain a Visual Studio `.sln`, `.vcxproj`, CMake build files, CSV inputs, Python plotting scripts, and CSV output routines.

## Kernel floor eta

Both projects now accept the kernel floor as the second command-line argument:

```text
VolatilitySmileTwoModels.exe <dataset> <eta>
```

Examples:

```text
VolatilitySmileTwoModels.exe eurusd 0.10
VolatilitySmileTwoModels.exe eurusd 0.02
VolatilitySmileTwoModels.exe usdhkd 0.10
VolatilitySmileTwoModels.exe usdhkd 0.02
```

For equities:

```text
VolatilitySmileTwoModels.exe spy_rich 0.10
VolatilitySmileTwoModels.exe spy_rich 0.02
VolatilitySmileTwoModels.exe tsla_rich 0.10
VolatilitySmileTwoModels.exe eurostoxx50 0.10
```

The supplied `run_fx_examples.bat` and `run_equity_examples.bat` automate the 10% / 2% runs and preserve each output set in a separate folder.

## USD/HKD

`cpp/FX_FlatMixture/data/usdhkd_20250521_delta.csv` contains the four maturities used in the revised article (1M, 3M, 6M, 1Y), with the five forward-delta quotes 10P, 25P, ATM, 25C, 10C.

## Building with Visual Studio

Open the `.sln` in the desired C++ project, retarget the Platform Toolset if Visual Studio asks, select `Release | x64`, and build.

Alternatively with CMake:

```text
cmake -S . -B build
cmake --build build --config Release
```

## Plotting

Create a Python environment and install:

```text
pip install pandas numpy matplotlib
```

Then run the plotting script in the project's `python/` folder after the C++ executable has generated the CSV files.

## Article

Compile from `article/` with:

```text
latexmk -pdf main.tex
```

The already compiled PDF is included as `article/revised_article_v4.pdf`.
