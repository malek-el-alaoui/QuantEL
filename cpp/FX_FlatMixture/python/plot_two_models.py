from pathlib import Path
import re

import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "output"
PLOTS = OUTPUT / "plots"
MARKET_DIR = PLOTS / "market_vs_models"
SMILE_DIR = PLOTS / "continuous_smiles_delta_range"
DENSITY_DIR = PLOTS / "densities_delta_range"
COMPONENT_DIR = PLOTS / "model1_component_kernels"

for folder in (MARKET_DIR, SMILE_DIR, DENSITY_DIR, COMPONENT_DIR):
    folder.mkdir(parents=True, exist_ok=True)


def safe_name(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value))


market = pd.read_csv(OUTPUT / "all_market_vs_models.csv")
smiles = pd.read_csv(OUTPUT / "all_continuous_smiles.csv")
densities = pd.read_csv(OUTPUT / "all_densities.csv")

# 1) Market quote fit.
for (instrument, maturity), g in market.groupby(["instrument", "maturity"]):
    g = g.sort_values("k")
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(g["k"], 100 * g["market_iv"], marker="o", linestyle="None", label="Market")
    ax.plot(g["k"], 100 * g["model1_lnmix_iv"], marker="x", linestyle="--", label="Model 1 - balanced-boundary mixture")
    ax.plot(g["k"], 100 * g["model2_cubic_iv"], marker="+", linestyle=":", label="Model 2 - cubic variance spline")
    ax.set_xlabel("Moneyness K/F")
    ax.set_ylabel("Implied volatility (%)")
    ax.set_title(f"{instrument} - {maturity}: market vs models")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(MARKET_DIR / f"{safe_name(instrument)}_{safe_name(maturity)}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

# 2) Delta-range continuous smiles. The full x-domain is generated from a
# market-wing forward-delta cutoff (2% by default); the shaded strip is the
# actually observed quote range.
for (instrument, maturity), g in smiles.groupby(["instrument", "maturity"]):
    g = g.sort_values("k")
    q = market[(market["instrument"] == instrument) & (market["maturity"] == maturity)].sort_values("k")
    k_lo = float(g["market_k_min"].iloc[0])
    k_hi = float(g["market_k_max"].iloc[0])

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(g["k"], 100 * g["model1_lnmix_iv"], label="Model 1 - balanced-boundary mixture")
    ax.plot(g["k"], 100 * g["model2_cubic_iv"], linestyle="--", label="Model 2 - cubic variance spline")
    ax.scatter(q["k"], 100 * q["market_iv"], marker="o", s=30, label="Market quotes", zorder=5)
    ax.axvspan(k_lo, k_hi, alpha=0.08, label="Observed market range")
    ax.axvline(k_lo, linewidth=0.8, linestyle=":")
    ax.axvline(k_hi, linewidth=0.8, linestyle=":")
    ax.set_xlabel("Moneyness K/F")
    ax.set_ylabel("Implied volatility (%)")
    delta = 100 * float(g["range_delta_cutoff"].iloc[0])
    ax.set_title(f"{instrument} - {maturity}: smile to ~{delta:g}-delta wings")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(SMILE_DIR / f"{safe_name(instrument)}_{safe_name(maturity)}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)

# 3) Densities over the same market-proportional delta range. Negative cubic
# densities remain visible intentionally.
for (instrument, maturity), g in densities.groupby(["instrument", "maturity"]):
    g = g.sort_values("k")
    k_lo = float(g["market_k_min"].iloc[0])
    k_hi = float(g["market_k_max"].iloc[0])

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(g["k"], g["model1_density"], label="Model 1 - balanced-boundary mixture")
    ax.plot(g["k"], g["model2_density"], linestyle="--", label="Model 2 - cubic variance spline")
    ax.axhline(0.0, linewidth=1.0)
    ax.axvspan(k_lo, k_hi, alpha=0.08, label="Observed market range")
    ax.axvline(k_lo, linewidth=0.8, linestyle=":")
    ax.axvline(k_hi, linewidth=0.8, linestyle=":")
    ax.set_xlabel("Moneyness K/F")
    ax.set_ylabel("Normalized risk-neutral density")
    delta = 100 * float(g["range_delta_cutoff"].iloc[0])
    ax.set_title(f"{instrument} - {maturity}: densities to ~{delta:g}-delta wings")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(DENSITY_DIR / f"{safe_name(instrument)}_{safe_name(maturity)}.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


# 4) Model-1 individual weighted log-normal kernels.  This is useful for
# diagnosing endpoint humps: the boundary component is dashed while market
# centered components are solid.  Their sum is the Model-1 density.
components = []
for path in OUTPUT.glob("*_model1_components.csv"):
    if path.name.startswith("all_"):
        continue
    dfc = pd.read_csv(path)
    if not dfc.empty:
        components.append(dfc)
if components:
    import numpy as np
    import math
    comp = pd.concat(components, ignore_index=True)
    for (instrument, maturity), g in comp.groupby(["instrument", "maturity"]):
        dg = densities[(densities["instrument"] == instrument) & (densities["maturity"] == maturity)].sort_values("k")
        if dg.empty:
            continue
        xgrid = dg["k"].to_numpy(float)
        fig, ax = plt.subplots(figsize=(10, 6))
        total = np.zeros_like(xgrid)
        for _, row in g.iterrows():
            center = float(row["center_k"])
            weight = float(row["weight"])
            sigma = float(row["kernel_total_vol"])
            if weight <= 0 or sigma <= 0:
                continue
            z = np.log(xgrid / center) / sigma
            y = weight * np.exp(-0.5 * z * z) / (xgrid * sigma * math.sqrt(2.0 * math.pi))
            total += y
            ls = "--" if int(row["is_artificial_component"]) else "-"
            ax.plot(xgrid, y, linestyle=ls, linewidth=0.8, alpha=0.55)
        ax.plot(xgrid, total, linewidth=2.0, label="Total Model 1 density")
        ax.axhline(0.0, linewidth=1.0)
        ax.set_xlabel("Moneyness K/F")
        ax.set_ylabel("Weighted density contribution")
        ax.set_title(f"{instrument} - {maturity}: Model 1 component kernels")
        ax.grid(True, alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(COMPONENT_DIR / f"{safe_name(instrument)}_{safe_name(maturity)}.png", dpi=180, bbox_inches="tight")
        plt.close(fig)

print("Plots saved under:", PLOTS)
