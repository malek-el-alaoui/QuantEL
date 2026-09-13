MAIN FILES
==========
all_market_vs_models.csv       Market quotes vs Model 1 and Model 2.
all_model1_artificial_quote_diagnostics.csv  Boundary quote / minimum-weight diagnostics.
all_continuous_smiles.csv       Delta-range dense smile grids incl. extrapolation: Model 1 vs Model 2.
all_densities.csv               2-delta density grids: flat vs wing-controlled mixture vs cubic spline.
all_wing_extrapolation.csv       Extended 0.5-delta wing smile/density diagnostics.
all_multi_maturity_comparison.csv Both models at quoted and intermediate maturities.
summary.csv                     Calibration and arbitrage diagnostics.

MODEL DEFINITIONS
=================
Model 1 = one-point paper log-normal mixture with a balanced boundary mass split.
The artificial strike is chosen to minimize the larger of the boundary and adjacent-market butterfly masses; equivalently it minimizes their squared sum.
Because their sum is fixed, the unconstrained optimum equalizes the two local butterfly masses.
The artificial quote has only 1e-10 objective weight in calibration; every real market quote has weight 1.
All Model 1 kernel total vols satisfy 0.10*sigma_ATM <= sigma_i <= 3.0*sigma_ATM inside the constrained linear solve.
This equity-only cap prevents very broad kernels from dominating the remote right tail.
The left wing remains exactly the calibrated flat mixture. Only the right wing is modified.
For x=log(K/F)>x_R, sigma_j(x)=sigma_j_star*exp[-alpha_R*(sigma_j_star/sigma_ATM)^2*(x-x_R)^2].
Inside the market sigma_j(x)=sigma_j_star exactly. The quadratic attachment is C1 at x_R.
alpha_R is chosen as the smallest statically admissible value making the aggregate right-wing IV non-increasing after a tiny attachment buffer;
if this is impossible before the static-arbitrage boundary, 95% of the critical alpha is used and the residual overshoot is reported.
Model 2 = natural cubic spline of total implied variance w(x)=IV^2*T vs x=log(K/F).
The Model 2 density is computed analytically from w, w' and w'' using the standard butterfly factor g(x).
Unlike Model 1, an unconstrained cubic spline is NOT guaranteed to be butterfly-arbitrage-free.
Negative Model 2 densities are intentionally left visible in the CSV instead of being clipped.
Dense outputs use a market-proportional 2-delta range based on the nearest observed wing volatilities.
Change OUTPUT_DELTA_CUTOFF in src/main.cpp from 0.02 to 0.05 for a 5-delta range.
Rows outside the observed quote interval are marked is_extrapolation=1.
If cubic extrapolation makes raw total variance non-positive, Model 2 IV is left blank and model2_variance_valid=0.

PLOTTING
========
Run python/plot_two_models.py after installing pandas and matplotlib.
