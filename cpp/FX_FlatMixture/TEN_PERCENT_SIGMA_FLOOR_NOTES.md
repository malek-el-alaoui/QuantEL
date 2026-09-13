# Model 1: one boundary artificial quote + 10% ATM kernel-volatility floor

This project deliberately removes the recursive midpoint artificial-quote refinement.

For each maturity, Model 1 uses only the single boundary artificial quote from the practical one-point construction. No midpoint quotes are inserted between market strikes.

The nonlinear calibration retains the hard lower bound

    sigma_i > 0.10 * sigma_ATM

inside the bound-constrained linearized Levenberg-Marquardt solve. The default eta is 0.10 and effectiveEta = max(0.10, eta).

The boundary artificial quote keeps its very small objective weight from the one-point project. Model 2 remains the natural cubic spline of total implied variance versus log-moneyness. The market-proportional output range remains the 2-delta range by default.
