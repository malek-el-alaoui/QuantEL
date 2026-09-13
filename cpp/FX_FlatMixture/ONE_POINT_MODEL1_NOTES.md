# Model 1 notes - one artificial quote, low objective weight

## Boundary construction

For each maturity:

- `p1 = c1 + k1 - 1`;
- `cN` is the rightmost normalized market call.

The practical one-point construction uses:

- if `p1 > cN`: one left artificial quote satisfying `p0 = cN`;
- if `cN > p1`: one right artificial quote satisfying `cN+1 = p1`;
- if they coincide numerically: no artificial quote.

The artificial **mixture mass** is kept small. On the left the code begins close to the finite theoretical infimum; on the right the search starts at `wN+1 = 1e-8` because the mathematical infimum zero is reached only at infinite strike.

## Artificial quote calibration weight

The nonlinear fit now gives the synthetic quote a separate objective weight:

```text
market quote objective weight     = 1
artificial quote objective weight = 1e-8
```

Thus the artificial point controls the mixture construction but has negligible influence compared with actual market observations when solving the kernel volatilities.

The value is written to `*_model1_artificial_quote_diagnostics.csv`.

## Constrained kernel vols

The Gauss-Newton / LM subproblem is solved under

```text
sigma_i > 0.10 * sigma_ATM
```

with the lower bounds active inside the quadratic solve.

## Dense extrapolation range

Continuous smiles and densities are generated only to a market-proportional delta wing, default 2-delta. Change `OUTPUT_DELTA_CUTOFF` from `0.02` to `0.05` in `src/main.cpp` to use 5-delta.
