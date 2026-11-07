# Equity-only bounded-kernel mixture with right-wing control

## Model 1 calibration

This project returns to the clean one-artificial-quote log-normal mixture. The synthetic quote has objective weight `1e-10`; all real market quotes have weight 1.

The kernel total volatilities are solved under the box constraint

`0.10 * sigma_ATM <= sigma_j <= 3.0 * sigma_ATM`.

The upper bound is imposed **inside the bound-constrained linearized LM solve**, not clipped only after solving.

## Equity right-wing rule

The left wing is left unchanged. Inside the observed market interval every kernel remains exactly at its calibrated value.

For `x = log(K/F) > x_R`,

`sigma_j(x) = sigma_j* exp[-alpha_R * (sigma_j*/sigma_ATM)^2 * (x-x_R)^2]`.

Thus broader kernels are damped more aggressively. This is deliberate: very broad kernels were the main cause of the ascending TSLA/SPX right tails.

The code first finds the largest `alpha_R` compatible with direct normalized-call monotonicity/convexity, nonnegative density and the implied-total-variance butterfly test. It then searches for the **smallest** admissible alpha that makes the aggregate right-wing IV non-increasing after a very small C1 attachment buffer. If that target is impossible before the static-arbitrage limit, it uses 95% of the critical alpha and reports the residual overshoot.

Because the quadratic attachment has zero first derivative exactly at `x_R`, an observed last-market upturn cannot always be reversed instantaneously while keeping the market fit exactly unchanged. TSLA is the main example: the long ascending extrapolation is removed, but a small shoulder remains before the wing turns down.

## Key validation

- lower sigma-bound violations: 0
- upper sigma-bound violations: 0
- max flat-vs-winged IV change inside observed market range: 0.000e+00 bp
- negative Model-1 density grid points on extended wing outputs: 0
- minimum tested Model-1 density: 0
- minimum tested butterfly factor g: 0.00864475
- all direct call static-arbitrage checks pass: yes
- artificial quote objective weight: <= 1.0e-10
- max |artificial boundary mass - adjacent market mass|: 1.290e-13
