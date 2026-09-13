# Rich public equity-chain data

This project replaces the earlier sparse 5–6 point equity examples with public single-expiry option-chain snapshots containing 22–30 retained market quotes per slice.

## Public source

The raw option-chain values were taken from WhaleQuant single-expiry chain pages for the 2026-10-16 expiry, retrieved on 2026-09-13:

- SPY: https://whalequant.io/en/stocks/SPY/options/2026-10-16
- TSLA: https://whalequant.io/en/stocks/TSLA/options/2026-10-16
- NVDA: https://whalequant.io/en/stocks/NVDA/options/2026-10-16
- AAPL: https://whalequant.io/en/stocks/AAPL/options/2026-10-16

The source pages publish strike, bid, ask, volume, open interest, and implied volatility for calls and puts.

## Final retained slices

- SPY_RICH: 30 call-IV quotes.
- TSLA_RICH: 30 call-IV quotes.
- NVDA_RICH: 22 call-IV quotes.
- AAPL_RICH: 22 call-IV quotes.

We use one quoted-IV convention within each slice (call IV) instead of mixing call-IV and put-IV conventions. Source bid/ask, volume, open interest and source IV are retained in the input CSVs.

## Forward proxies

The calibration uses normalized moneyness K/F. Forward proxies are estimated from near-ATM call-put parity using quotes from the same expiry and DF approximately equal to one over the short 33-day horizon:

- SPY: F = 765.30
- TSLA: F = 366.21
- NVDA: F = 219.108
- AAPL: F = 333.558

These are practical proxies for the experiment, not exchange-published forwards.

## Static-arbitrage cleaning

Vendor IV grids can contain local quote noise. Before calibration each candidate IV is converted to a normalized Black call and the retained set is required to satisfy the discrete decreasing/convex call conditions.

Only offending nodes are removed; retained source IVs are not altered.

The final four input slices all pass the project's discrete static-arbitrage test.

## Model comparison

Model 1 is the latest equity-only mixture specification:

- one balanced boundary artificial quote,
- artificial-quote objective weight 1e-10,
- 0.10 sigma_ATM <= sigma_j <= 3.0 sigma_ATM,
- component-weighted quadratic right-wing decay outside the last observed market strike,
- direct call monotonicity/convexity and nonnegative-density safety checks.

Model 2 is the natural cubic spline of total implied variance versus log-moneyness.

The wing rule never changes any market quote: it switches on only beyond the last observed market strike.


EURO STOXX 50 rich chain: this project uses data/eurostoxx50_rich_20140930.csv, a 31-strike per maturity dense chain obtained by arbitrage-preserving piecewise-linear interpolation of the original observed call mid quotes, then reconverted to implied vols.
