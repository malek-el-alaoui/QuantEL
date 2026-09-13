Revised article package (v4)
============================

Files
- main.tex: LaTeX source.
- revised_article_v4.pdf: compiled article.
- figures/: all PDF/PNG figures used by the article.
- make_figures_v4.py: script used to regenerate the numerical figures from the calibration outputs available in the working environment.

Revision highlights
- Preserves the organization and mathematical construction of the original article.
- Calendar-spread compact interval corrected from [0,k_max] to [eta,k_max], eta>0.
- Cousot (2007) added and explicitly distinguished from the strike-wise chordal interpolation used in the article.
- FX results include EUR/USD and all four USD/HKD maturities from the cited BofA table.
- Both eta=10% and eta=2% kernel-width floors are reported for FX and equity/index examples.
- Equity wing extrapolation remains a separate numerical study; SPX includes mild, medium and full right-wing control strengths.
- Market quotes are plotted with small square markers.

Compile with:
  latexmk -pdf main.tex
