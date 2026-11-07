
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent

def fmt_T(T):
    return f"{T:.2f}".replace(".", "_")

def make_scenario_plots(name, title_prefix):
    outdir = BASE / "plots" / name
    outdir.mkdir(parents=True, exist_ok=True)
    cmp = pd.read_csv(BASE / "output" / name / "model_comparison.csv")
    raw = pd.read_csv(BASE / "output" / name / "classic_dupire_local_variance.csv")
    ext = pd.read_csv(BASE / "output" / name / "extended_complex_dupire.csv")

    # Separate price graph for every maturity.
    for T in sorted(cmp["maturity"].unique()):
        g = cmp[np.isclose(cmp["maturity"], T)].copy()
        plt.figure(figsize=(8.2, 5.2))
        plt.scatter(g["strike"], g["market_price"], s=24, label="Market")
        plt.plot(g["strike"], g["extended_price"], linewidth=2.0, label="Extended Dupire")
        plt.plot(g["strike"], g["classic_price"], linewidth=1.8, label="Classic Dupire (admissible)")
        plt.xlabel("Strike K")
        plt.ylabel("Call price")
        plt.title(f"{title_prefix} - T={T:g}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(outdir / f"smile_prices_T{fmt_T(T)}.png", dpi=190)
        plt.close()

        # Separate local-variance diagnostic for the same maturity.
        r = raw[np.isclose(raw["maturity"], T)].copy()
        e = ext[np.isclose(ext["maturity"], T)].copy()
        plt.figure(figsize=(8.2, 5.2))
        plt.plot(r["strike"], r["raw_dupire_variance"], label="Raw classic Dupire variance")
        plt.plot(r["strike"], r["regularized_variance"], label="Admissible classic variance")
        plt.plot(e["strike"], e["sigma_ext2_real"], label="Re extended variance")
        plt.plot(e["strike"], e["sigma_ext2_imag"], label="Im extended variance")
        plt.axhline(0.0, linewidth=1)
        plt.xlabel("Strike K")
        plt.ylabel("Variance")
        plt.title(f"Local-variance diagnostics - {title_prefix}, T={T:g}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(outdir / f"local_variance_T{fmt_T(T)}.png", dpi=190)
        plt.close()

        # Market gamma: highlights butterfly violations.
        plt.figure(figsize=(8.2, 5.2))
        plt.plot(r["strike"], r["gamma_fd"], marker="o", markersize=3)
        plt.axhline(0.0, linewidth=1)
        plt.xlabel("Strike K")
        plt.ylabel("Discrete C_KK")
        plt.title(f"Market strike convexity - {title_prefix}, T={T:g}")
        plt.tight_layout()
        plt.savefig(outdir / f"market_gamma_T{fmt_T(T)}.png", dpi=190)
        plt.close()

    # All quote errors.
    plt.figure(figsize=(8.2, 5.2))
    for T in sorted(cmp["maturity"].unique()):
        g = cmp[np.isclose(cmp["maturity"], T)]
        plt.plot(g["strike"], g["classic_abs_error"], label=f"Classic T={T:g}")
    plt.xlabel("Strike K")
    plt.ylabel("Absolute price error")
    plt.title(f"Classic Dupire pricing errors - {title_prefix}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "classic_errors_all_maturities.png", dpi=190)
    plt.close()

    plt.figure(figsize=(8.2, 5.2))
    for T in sorted(cmp["maturity"].unique()):
        g = cmp[np.isclose(cmp["maturity"], T)]
        plt.plot(g["strike"], g["extended_abs_error"], label=f"Extended T={T:g}")
    plt.yscale("log")
    plt.xlabel("Strike K")
    plt.ylabel("Absolute price error")
    plt.title(f"Extended Dupire reproduction errors - {title_prefix}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "extended_errors_all_maturities.png", dpi=190)
    plt.close()

make_scenario_plots("regular", "Moderate arbitrage surface")
make_scenario_plots("violent", "Violent butterfly-arbitrage surface")

# Cross-scenario RMSE comparison.
rows = []
for name, label in [("regular", "Moderate"), ("violent", "Violent")]:
    d = pd.read_csv(BASE / "output" / name / "model_comparison.csv")
    rows.append({
        "scenario": label,
        "classic": float(np.sqrt(np.mean(d["classic_error"]**2))),
        "extended": float(np.sqrt(np.mean(d["extended_error"]**2))),
    })
m = pd.DataFrame(rows)
m.to_csv(BASE / "output" / "rmse_comparison.csv", index=False)

x = np.arange(len(m))
w = 0.35
plt.figure(figsize=(7.4, 4.8))
plt.bar(x - w/2, m["classic"], width=w, label="Classic Dupire")
plt.bar(x + w/2, m["extended"], width=w, label="Extended Dupire")
plt.xticks(x, m["scenario"])
plt.yscale("log")
plt.ylabel("RMSE (log scale)")
plt.title("Pricing-fit comparison")
plt.legend()
plt.tight_layout()
plt.savefig(BASE / "plots" / "rmse_comparison.png", dpi=190)
plt.close()

# Emphasized violent T=0.50 figure.
d = pd.read_csv(BASE / "output" / "violent" / "model_comparison.csv")
g = d[np.isclose(d["maturity"], 0.50)]
plt.figure(figsize=(9.0, 5.8))
plt.scatter(g["strike"], g["market_price"], s=32, label="Market")
plt.plot(g["strike"], g["extended_price"], linewidth=2.2, label="Extended Dupire")
plt.plot(g["strike"], g["classic_price"], linewidth=2.0, label="Classic Dupire (admissible)")
plt.xlabel("Strike K")
plt.ylabel("Call price")
plt.title("Violent arbitrage test at T=0.50")
plt.legend()
plt.tight_layout()
plt.savefig(BASE / "plots" / "violent_smile_T0_50_highlight.png", dpi=210)
plt.close()
