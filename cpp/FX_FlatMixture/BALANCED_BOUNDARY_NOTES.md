# Balanced one-artificial-quote Model 1

This version removes all recursive/interior artificial-point experiments.  The
only synthetic quote is the **single boundary quote** used by the practical
one-point construction.

## One-point rule

For market calls `(k_i,c_i)`, define `p1=c1+k1-1` and let `cN` be the last
market call.

- If `p1 > cN`, add one left quote satisfying `p0=cN`, i.e.
  `c0=1-k0+cN`.
- If `cN > p1`, add one right quote satisfying `c_{N+1}=p1`.

This is the same one-point branch used previously; only the artificial strike
selection is changed.

## New artificial-strike choice: minimize the two boundary butterfly masses

### Left branch

Let `d1` be the artificial-to-first-market slope and `d2` the first
market-to-market slope.  The two local butterfly masses are

    w0 = 1 + d1,
    w1 = d2 - d1.

Their sum is fixed:

    w0 + w1 = 1 + d2.

Therefore it is impossible to make both arbitrarily small without changing the
market quotes.  The unique minimizer of either

    max(w0,w1)

or

    w0^2 + w1^2

is the equal split

    w0 = w1 = (1+d2)/2,
    d1* = (d2-1)/2.

Using `p0=cN`,

    d1 = -1 + (p1-cN)/(k1-k0),

so the artificial strike follows explicitly.  If the equal-split point would
fall outside the strictly admissible interval, the slope is projected onto the
feasible interval.

### Right branch

Let `dN` be the last market-market slope and `dR` the last-market-to-artificial
slope.  Then

    wN   = dR-dN,
    wN+1 = -dR,
    wN+wN+1 = -dN.

The same minimax/L2 argument gives

    dR = dN/2,
    wN = wN+1 = -dN/2.

The artificial strike is then obtained from `c_{N+1}=p1`.

## Calibration weight and sigma floor

- Real market quote objective weight: `1.0`.
- Single artificial quote objective weight: `1e-10`.
- Kernel-volatility constraint: `sigma_j > 0.10*sigma_ATM` inside the
  bound-constrained LM subproblem.

The tiny artificial residual weight means the boundary quote determines the
mixture probability/forward construction but has essentially no direct pull on
the market fit.

## Important limitation

Balancing the two boundary masses is the best possible **local mass split** for
fixed market quotes, but it does not mathematically guarantee that the final
continuous density is hump-free.  Kernel widths are calibrated variables and
short-dated smiles can still create narrow peaks.  The component-kernel plots
in `output/plots/model1_component_kernels/` are included specifically to make
this visible rather than hide it.
