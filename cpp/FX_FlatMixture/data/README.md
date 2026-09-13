# Bundled market snapshots

These CSVs are **offline test snapshots**, not live feeds and not investment data.

## EURO STOXX 50 — 2014-09-30

`eurostoxx50_20140930.csv` contains a selected, convex set of historical EURO STOXX 50 European call quotes for five expiries. The public Yves Hilpisch CQF / *Derivatives Analytics with Python* material labels these as real EURO STOXX 50 option quotes and attributes the option/index data to Eurex and STOXX.

Source material:
- https://gist.github.com/yhilpisch/e2af84bc2ff1cab88d6be019459f5672
- Original source attribution in that material: www.eurexchange.com and www.stoxx.com

The file stores the source call mids and implied vols derived from those prices with the historical sample convention used by the source code: index level 3225.93, ACT/365 maturity and r = -5%. The calibration itself consumes normalized K/F and annualized implied vol.

## SPX — 2026-07-30 delayed snapshot

`spx_20260730.csv` contains a small convex subset of Cboe SPXW calls for the 2026-07-31 expiry observed on Cboe's delayed quote page on 2026-07-30. The CSV preserves the actual bid/ask midpoint and the source page's rounded IV. `market_iv` is recomputed from the midpoint for a more precise fit comparison, using T=1/365 and a documented 4% rate/no-dividend demo assumption.

Source:
- https://www.cboe.com/delayed_quotes/spx/quote_table

This is deliberately a tiny near-ATM slice. Replace it with your own production chain if you need a full SPX smile.

## EUR/USD — 2024-05-15

`eurusd_20240515_delta.csv` contains published EUR/USD ATM, 25-delta call/put, and 10-delta call/put vol observations across 1W, 1M, 2M, 3M, 6M, 9M, 1Y, 18M and 2Y.

Source:
- https://finpricing.com/lib/FxVolIntroduction.html

The source itself explains that FX quotes depend on the ATM and delta convention. For a self-contained demo, `loadFxDeltaSurfaceCsv` converts the quotes to K/F using the **forward-delta** convention and sets ATM to K/F=1. This conversion assumption is not claimed to be the original broker convention. For production, replace it with the exact pair/tenor convention used by your market data source.


## Additional US equity tests — Cboe delayed chains, September 2026

The files `aapl_cboe_20260910.csv`, `nvda_cboe_20260910.csv`, and
`tsla_cboe_20260911.csv` are small near-the-money test slices transcribed from
Cboe delayed quote pages.  They are included only to stress-test the numerical
construction on additional equity smiles.

Sources:
- AAPL: https://www.cboe.com/delayed_quotes/AAPL/quote_table/
- NVDA: https://www.cboe.com/delayed_quotes/nvda/quote_table/
- TSLA: https://www.cboe.com/delayed_quotes/tsla/quote_table/

For these three small tests the model does **not** assume a spot/rate/dividend
forward.  Instead it reconstructs normalized moneyness directly from the
Cboe-displayed Black IV and call delta under the forward-delta identity

    d1 = N^{-1}(Delta_F),
    K/F = exp(0.5*v^2 - d1*v),  v = IV*sqrt(T).

The CSV therefore stores this normalized K/F in the `strike` column with
`forward=1`.  The physical exchange strike and reported delta are retained in
`pillar`; displayed call mid and IV are retained as source metadata.  This
makes the test self-contained and avoids silently assuming an equity forward.
