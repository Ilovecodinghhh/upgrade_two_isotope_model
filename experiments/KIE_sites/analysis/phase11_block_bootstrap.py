#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 11: year-block bootstrap diagnostics for KIE_sites."""

from pathlib import Path
import json

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


EXPT_DIR = Path(__file__).resolve().parent.parent
PHASE1_DIR = EXPT_DIR / "results" / "phase1_data"
OUT_DIR = EXPT_DIR / "results" / "phase11_block_bootstrap"
FIG_DIR = EXPT_DIR / "figures"
OUT_JSON = OUT_DIR / "block_bootstrap_results.json"
OUT_FIG = FIG_DIR / "fig20_block_bootstrap_alpha.png"

SITE_CODES = ["ALT", "ZEP", "BRW", "CBA", "MHD", "KUM", "CGO", "SPO"]
ALPHA_D_OH = 1.294


def block_bootstrap_years(years, n_boot=1000, seed=123):
    """Draw year blocks with replacement; return Python lists for stable tests."""
    unique_years = np.array(sorted(set(int(y) for y in years)), dtype=int)
    rng = np.random.default_rng(seed)
    return [
        [int(y) for y in rng.choice(unique_years, size=len(unique_years), replace=True)]
        for _ in range(n_boot)
    ]


def fit_harmonic(t, y):
    """Fit annual harmonic plus intercept to a resampled monthly series."""
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    x = 2 * np.pi * t
    design = np.column_stack([np.ones_like(t), np.sin(x), np.cos(x)])
    coeffs, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
    _, B, C = coeffs
    return {"B": float(B), "C": float(C), "amplitude": float(np.hypot(B, C))}


def fit_resampled_site(df, sampled_years):
    """Fit one site after resampling whole years with replacement."""
    pieces = []
    for i, year in enumerate(sampled_years):
        block = df[df["year"] == int(year)].copy()
        if block.empty:
            continue
        block["decimal_year"] = block["decimal_year"] + i * 100.0
        pieces.append(block)
    if not pieces:
        return {"ratio": np.nan, "alpha_13c_oh": np.nan, "n_months": 0}
    sample = pd.concat(pieces, ignore_index=True)
    if len(sample) < 8:
        return {"ratio": np.nan, "alpha_13c_oh": np.nan, "n_months": int(len(sample))}
    f13 = fit_harmonic(sample["decimal_year"], sample["d13C_mean"])
    fD = fit_harmonic(sample["decimal_year"], sample["dD_mean"])
    ratio = f13["amplitude"] / fD["amplitude"] if fD["amplitude"] > 0 else np.nan
    return {
        "ratio": float(ratio),
        "alpha_13c_oh": float(1.0 + ratio * (ALPHA_D_OH - 1.0)),
        "n_months": int(len(sample)),
    }


def summarize_samples(values):
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if len(arr) == 0:
        return {"median": np.nan, "ci95": [np.nan, np.nan], "n": 0}
    return {
        "median": float(np.median(arr)),
        "ci95": [float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))],
        "n": int(len(arr)),
    }


def run_site_bootstrap(code, n_boot=500, seed=123):
    path = PHASE1_DIR / f"site_monthly_{code}.csv"
    df = pd.read_csv(path)
    years = sorted(df["year"].unique())
    draws = block_bootstrap_years(years, n_boot=n_boot, seed=seed)
    fits = [fit_resampled_site(df, draw) for draw in draws]
    return {
        "years": [int(y) for y in years],
        "n_boot": n_boot,
        "ratio": summarize_samples([f["ratio"] for f in fits]),
        "alpha_13c_oh": summarize_samples([f["alpha_13c_oh"] for f in fits]),
    }


def plot_results(results):
    codes = list(results["sites"].keys())
    med = [results["sites"][c]["alpha_13c_oh"]["median"] for c in codes]
    lo = [results["sites"][c]["alpha_13c_oh"]["ci95"][0] for c in codes]
    hi = [results["sites"][c]["alpha_13c_oh"]["ci95"][1] for c in codes]
    x = np.arange(len(codes))
    yerr = np.vstack([np.array(med) - np.array(lo), np.array(hi) - np.array(med)])

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.errorbar(x, med, yerr=yerr, fmt="o", capsize=3)
    ax.axhline(1.0039, color="C2", ls="--", lw=1, label="Saueressig")
    ax.axhline(1.0054, color="C3", ls="--", lw=1, label="Cantrell")
    ax.set_xticks(x)
    ax.set_xticklabels(codes)
    ax.set_ylabel("Year-block bootstrap alpha13C_OH")
    ax.set_title("Year-block bootstrap sensitivity")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=180)
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    sites = {code: run_site_bootstrap(code) for code in SITE_CODES if (PHASE1_DIR / f"site_monthly_{code}.csv").exists()}
    sh_values = []
    for code in ["CGO", "SPO"]:
        if code in sites:
            sh_values.append(sites[code]["alpha_13c_oh"]["median"])
    output = {
        "metadata": {"method": "Year-block bootstrap of Phase 1 monthly paired data"},
        "sites": sites,
        "sh_only_median_alpha": float(np.nanmedian(sh_values)) if sh_values else np.nan,
    }
    OUT_JSON.write_text(json.dumps(output, indent=2), encoding="utf-8")
    plot_results(output)
    print(f"Wrote {OUT_JSON}")
    print(f"Wrote {OUT_FIG}")


if __name__ == "__main__":
    main()
