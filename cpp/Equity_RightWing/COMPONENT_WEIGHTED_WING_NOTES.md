# Component-weighted all-kernel wing extrapolation

This revision fixes the main weakness of the previous all-kernel experiment.  The previous version multiplied every kernel sigma by the SAME wing factor, so it mostly preserved the mixture geometry and only rescaled all widths together.

For x=log(K/F) outside the observed market interval, component j now uses

    sigma_j(x) = sigma_j* exp(-alpha_side * rho_j,side * d(x)^3)

where rho is fixed by the component center.  On the left wing, components centered farther to the right receive stronger damping; on the right wing, components centered farther to the left receive stronger damping.  Explicitly, with xL=log(kmin), xR=log(kmax):

    rho_j,L = 1 + 3 clamp((x_j-xL)/(xR-xL),0,1)
    rho_j,R = 1 + 3 clamp((xR-x_j)/(xR-xL),0,1).

Thus alpha_L and alpha_R remain the only free wing parameters.  All sigmas remain exactly equal to the calibrated values for x in [xL,xR].  The exp(-alpha*d^3) join is C2 at the market boundaries.

The reference wing target is now a 30% reduction (ratio 0.70) of the wing-nearest component around the 0.5-delta diagnostic point, capped at 85% of the numerically determined no-arbitrage critical alpha.  The full mixture is then retested for call monotonicity/convexity, nonnegative density and total-variance butterfly g.

## Key validation

- Market-region wing-vs-flat IV difference: 0.000e+00 bp.
- Negative winged-density points on the diagnostic grids: 0.
- Minimum butterfly g on the output grids: 0.00849221.

See output/component_weighted_wing_summary.csv for per-slice changes.
