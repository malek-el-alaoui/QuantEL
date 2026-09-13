from pathlib import Path
import re
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "output"
PLOTS = OUTPUT / "plots"
MARKET_DIR = PLOTS / "market_vs_models"
SMILE_DIR = PLOTS / "continuous_smiles_2delta"
DENSITY_DIR = PLOTS / "densities_2delta"
COMPONENT_DIR = PLOTS / "model1_component_kernels_winged"
WING_SMILE_DIR = PLOTS / "wing_extrapolation_0p5delta_smiles"
WING_DENSITY_DIR = PLOTS / "wing_extrapolation_0p5delta_densities"
for folder in (MARKET_DIR, SMILE_DIR, DENSITY_DIR, COMPONENT_DIR, WING_SMILE_DIR, WING_DENSITY_DIR):
    folder.mkdir(parents=True, exist_ok=True)

def safe_name(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value))

def market_window(market_lo, market_hi, frac_left=0.05, frac_right=0.05):
    width = max(market_hi - market_lo, 1e-8)
    return float(market_lo - frac_left * width), float(market_hi + frac_right * width)


def smile_window(market_lo, market_hi, instrument=None):
    # Tighter zoom for standard interpolation plots, especially for dense EURO STOXX smiles.
    if instrument is not None and 'EUROSTOXX50' in str(instrument).upper():
        return market_window(market_lo, market_hi, 0.04, 0.04)
    return market_window(market_lo, market_hi, 0.06, 0.06)


def wing_window(k_values, market_lo, market_hi, instrument=None):
    vals = np.asarray(k_values, dtype=float)
    vals = vals[np.isfinite(vals) & (vals > 0)]
    if vals.size == 0:
        vals = np.array([market_lo, market_hi], dtype=float)
    width = max(market_hi - market_lo, 1e-8)
    # balanced but slightly right-biased window for extrapolation plots
    left_extra = min(max(market_lo - vals.min(), 0.0), 0.18 * width)
    right_extra = min(max(vals.max() - market_hi, 0.0), 0.28 * width)
    lower = market_lo - left_extra - 0.02 * width
    upper = market_hi + right_extra + 0.03 * width
    return float(lower), float(upper)


def add_market_markers(ax, market_lo, market_hi):
    ax.axvline(market_lo, linewidth=0.9, linestyle=':', alpha=0.9)
    ax.axvline(market_hi, linewidth=0.9, linestyle=':', alpha=0.9)


market = pd.read_csv(OUTPUT / "all_market_vs_models.csv")
smiles = pd.read_csv(OUTPUT / "all_continuous_smiles.csv")
densities = pd.read_csv(OUTPUT / "all_densities.csv")
wings = pd.read_csv(OUTPUT / "all_wing_extrapolation.csv")

# Market fit: wing parametrization is exactly flat on the observed range, so
# it preserves the calibrated market values.
for (instrument, maturity), g in market.groupby(["instrument", "maturity"]):
    g = g.sort_values("k")
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(g["k"], 100 * g["market_iv"], marker="o", linestyle="None", label="Market")
    ax.plot(g["k"], 100 * g["model1_lnmix_iv"], marker="x", linestyle="--", label="Model 1 - wing-controlled mixture")
    ax.plot(g["k"], 100 * g["model2_cubic_iv"], marker="+", linestyle=":", label="Model 2 - cubic variance spline")
    ax.set_xlim(*smile_window(float(g["k"].min()), float(g["k"].max()), instrument))
    ax.set_xlabel("Moneyness K/F")
    ax.set_ylabel("Implied volatility (%)")
    ax.set_title(f"{instrument} - {maturity}: market fit")
    ax.grid(True, alpha=0.3); ax.legend(); fig.tight_layout()
    fig.savefig(MARKET_DIR / f"{safe_name(instrument)}_{safe_name(maturity)}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

# Standard 2-delta smile comparison: baseline flat-kernel mixture is shown as
# a diagnostic so the effect of the new wing parametrization is explicit.
for (instrument, maturity), g in smiles.groupby(["instrument", "maturity"]):
    g = g.sort_values("k")
    q = market[(market["instrument"] == instrument) & (market["maturity"] == maturity)].sort_values("k")
    k_lo = float(g["market_k_min"].iloc[0]); k_hi = float(g["market_k_max"].iloc[0])
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(g["k"], 100 * g["model1_flat_iv"], linestyle=":", label="Model 1 baseline - flat kernel sigmas")
    ax.plot(g["k"], 100 * g["model1_winged_iv"], label="Model 1 - bounded-kernel right-wing control")
    ax.plot(g["k"], 100 * g["model2_cubic_iv"], linestyle="--", label="Model 2 - cubic variance spline")
    ax.scatter(q["k"], 100 * q["market_iv"], marker="o", s=28, label="Market quotes", zorder=5)
    add_market_markers(ax, k_lo, k_hi)
    ax.set_xlim(*smile_window(k_lo, k_hi, instrument))
    ax.set_xlabel("Moneyness K/F"); ax.set_ylabel("Implied volatility (%)")
    ax.set_title(f"{instrument} - {maturity}: 2-delta wing comparison")
    ax.grid(True, alpha=0.3); ax.legend(); fig.tight_layout()
    fig.savefig(SMILE_DIR / f"{safe_name(instrument)}_{safe_name(maturity)}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

for (instrument, maturity), g in densities.groupby(["instrument", "maturity"]):
    g = g.sort_values("k")
    k_lo = float(g["market_k_min"].iloc[0]); k_hi = float(g["market_k_max"].iloc[0])
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(g["k"], g["model1_flat_density"], linestyle=":", label="Model 1 baseline density")
    ax.plot(g["k"], g["model1_winged_density"], label="Model 1 bounded-kernel right-wing density")
    ax.plot(g["k"], g["model2_density"], linestyle="--", label="Model 2 cubic density")
    ax.axhline(0.0, linewidth=1.0)
    add_market_markers(ax, k_lo, k_hi)
    ax.set_xlim(*smile_window(k_lo, k_hi, instrument))
    ax.set_xlabel("Moneyness K/F"); ax.set_ylabel("Normalized risk-neutral density")
    ax.set_title(f"{instrument} - {maturity}: density to 2-delta wings")
    ax.grid(True, alpha=0.3); ax.legend(); fig.tight_layout()
    fig.savefig(DENSITY_DIR / f"{safe_name(instrument)}_{safe_name(maturity)}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

# Wider 0.5-delta diagnostic view, useful for seeing the extrapolation itself.
for (instrument, maturity), g in wings.groupby(["instrument", "maturity"]):
    g = g.sort_values("k")
    k_lo = float(g["market_k_min"].iloc[0]); k_hi = float(g["market_k_max"].iloc[0])
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(g["k"], 100 * g["model1_flat_iv"], linestyle=":", label="Flat-kernel mixture")
    ax.plot(g["k"], 100 * g["model1_winged_iv"], label="Bounded-kernel right-wing mixture")
    ax.plot(g["k"], 100 * g["model2_cubic_iv"], linestyle="--", label="Cubic variance spline")
    add_market_markers(ax, k_lo, k_hi)
    ax.set_xlim(*wing_window(g["k"], k_lo, k_hi, instrument))
    ax.set_xlabel("Moneyness K/F"); ax.set_ylabel("Implied volatility (%)")
    ax.set_title(f"{instrument} - {maturity}: extended 0.5-delta wing extrapolation")
    ax.grid(True, alpha=0.3); ax.legend(); fig.tight_layout()
    fig.savefig(WING_SMILE_DIR / f"{safe_name(instrument)}_{safe_name(maturity)}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(g["k"], g["model1_flat_density"], linestyle=":", label="Flat-kernel mixture")
    ax.plot(g["k"], g["model1_winged_density"], label="Bounded-kernel right-wing mixture")
    ax.plot(g["k"], g["model2_density"], linestyle="--", label="Cubic variance spline")
    ax.axhline(0.0, linewidth=1.0)
    add_market_markers(ax, k_lo, k_hi)
    ax.set_xlim(*wing_window(g["k"], k_lo, k_hi, instrument))
    ax.set_xlabel("Moneyness K/F"); ax.set_ylabel("Normalized risk-neutral density")
    ax.set_title(f"{instrument} - {maturity}: extended 0.5-delta density")
    ax.grid(True, alpha=0.3); ax.legend(); fig.tight_layout()
    fig.savefig(WING_DENSITY_DIR / f"{safe_name(instrument)}_{safe_name(maturity)}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

# Individual kernels with component-dependent right-wing damping.
# The damping strength scales as (sigma_j/sigma_ATM)^2; the left wing remains flat.
parts=[]
for path in OUTPUT.glob("*_model1_components.csv"):
    if path.name.startswith("all_"): continue
    d=pd.read_csv(path)
    if not d.empty: parts.append(d)
if parts:
    comp=pd.concat(parts, ignore_index=True)
    for (instrument,maturity), g in comp.groupby(["instrument","maturity"]):
        dg=wings[(wings.instrument==instrument)&(wings.maturity==maturity)].sort_values("k")
        if dg.empty: continue
        xgrid=dg["k"].to_numpy(float)
        total=np.zeros_like(xgrid); flat_total=np.zeros_like(xgrid)
        kmin=float(g.market_k_min.iloc[0]); kmax=float(g.market_k_max.iloc[0])
        aL=float(g.left_wing_alpha.iloc[0]); aR=float(g.right_wing_alpha.iloc[0]); pwr=float(g.wing_power.iloc[0])
        atm_total=float(g.atm_total_vol.iloc[0]) if "atm_total_vol" in g.columns else float(np.median(g.kernel_total_vol))
        n=len(g)
        fig,ax=plt.subplots(figsize=(10,6))
        for idx,(_,row) in enumerate(g.sort_values("component").iterrows()):
            center=float(row.center_k); weight=float(row.weight); base=float(row.kernel_total_vol)
            sigma=np.full_like(xgrid,base)
            # Left wing intentionally remains flat for equity skew.
            mask=xgrid>kmax
            if np.any(mask):
                d=np.log(xgrid[mask]/kmax)
                factor=np.clip((base/max(atm_total,1e-12))**2,0.01,9.0)
                sigma[mask]=np.maximum(1e-10,base*np.exp(-aR*factor*d**pwr))
            z=np.log(xgrid/center)/sigma
            y=weight*np.exp(-0.5*z*z)/(xgrid*sigma*math.sqrt(2*math.pi))
            z0=np.log(xgrid/center)/base
            y0=weight*np.exp(-0.5*z0*z0)/(xgrid*base*math.sqrt(2*math.pi))
            total+=y; flat_total+=y0
            ax.plot(xgrid,y,linestyle="-",linewidth=0.75,alpha=0.45)
        ax.plot(xgrid,flat_total,linestyle=":",linewidth=1.7,label="Flat-kernel component sum")
        # The simple sum of variable-sigma lognormal PDFs is NOT the exact
        # density once sigma depends on strike, so show the exact density from
        # the CSV as the authoritative wing-controlled curve.
        ax.plot(dg["k"],dg["model1_winged_density"],linewidth=2.0,label="Exact wing-controlled density")
        add_market_markers(ax, kmin, kmax)
        ax.set_xlim(*wing_window(xgrid, kmin, kmax, instrument))
        ax.set_xlabel("Moneyness K/F"); ax.set_ylabel("Density / component contribution")
        ax.set_title(f"{instrument} - {maturity}: bounded-kernel right-wing control")
        ax.grid(True,alpha=0.3); ax.legend(); fig.tight_layout()
        fig.savefig(COMPONENT_DIR/f"{safe_name(instrument)}_{safe_name(maturity)}.png",dpi=180,bbox_inches="tight")
        plt.close(fig)

print("Plots saved under:", PLOTS)
