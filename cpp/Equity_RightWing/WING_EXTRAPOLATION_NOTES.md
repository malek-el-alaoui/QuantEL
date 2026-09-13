# All-kernel exponential wing parametrization

## Purpose

The calibrated one-boundary-quote log-normal mixture is left completely unchanged on the observed market interval.  Outside that interval, **every** calibrated kernel volatility is made strike-dependent.  This differs from the preceding experiment, where only the leftmost/rightmost kernels were changed and the aggregate smile sometimes barely moved because interior kernels dominated the tails.

Let `x = log(K/F)`, `xL = log(k_market_min)`, `xR = log(k_market_max)`, and let `sigma_j*` be the calibrated total volatility of component `j`.

For every component `j`:

```
sigma_j(x) = sigma_j*                                      xL <= x <= xR
sigma_j(x) = sigma_j* exp[-alpha_L (xL-x)^3]               x < xL
sigma_j(x) = sigma_j* exp[-alpha_R (x-xR)^3]               x > xR
```

Thus all kernels receive the same side-specific multiplicative factor. Their calibrated relative widths are preserved.  The cubic distance power makes the multiplier, its first derivative, and its second derivative join the flat interior continuously at the market boundary.

## Choice of alpha

Static no-arbitrage does not determine a unique decay speed.  For each side the code first searches for the largest one-sided `alpha_crit` that passes a dense safety-grid test of:

1. normalized call bounds, monotonicity and convexity directly in strike;
2. no negative density grid points; and
3. the implied-total-variance butterfly factor `g(x)` remaining above a small finite-difference tolerance.

A transparent shape target is then used: at an approximate 0.5-delta diagnostic wing point (`delta = 0.005`), target a common kernel-sigma multiplier of 0.90.  The selected alpha is

```
alpha = min(alpha_target, 0.90 * alpha_crit).
```

Finally both wings are turned on together and retested.  If their interaction were to fail the static tests, both alphas would be scaled down together.  In the supplied runs that extra scale is 1.0.

## Important consequence

Because the multiplier is applied to all components, it now materially changes the aggregate equity extrapolation.  At the ~0.5-delta endpoints in the supplied tests:

- AAPL: left IV 36.59% -> 33.64%; right IV 46.91% -> 44.44%.
- NVDA: left IV 15.81% -> 14.42%; right IV 28.22% -> 26.55%.
- TSLA: left IV 64.09% -> 63.13%; right IV 62.16% -> 60.88%.
- SPX: left IV 44.14% -> 43.54%; right IV 42.96% -> 42.51%.

The TSLA/SPX multipliers cannot reach the 0.90 target on the diagnostic endpoints because the no-arbitrage critical alpha is binding.  This is intentional: the code prioritizes the static-arbitrage constraints over the shape target.

## Market fit

For every market strike, the multiplier is exactly 1.  Therefore the calibrated call prices and implied volatilities at all market quotes are unchanged by the wing parametrization.
