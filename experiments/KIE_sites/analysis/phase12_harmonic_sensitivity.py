#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Harmonic-model sensitivity checks for KIE_sites."""

from pathlib import Path
import json

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


EXPT_DIR = Path(__file__).resolve().parent.parent
PHASE1_DIR = EXPT_DIR / "results" / "phase1_data"
OUT_DIR = EXPT_DIR / "results" / "phase12_harmonic_sensitivity"
FIG_DIR = EXPT_DIR / "figures"
OUT_JSON = OUT_DIR / "harmonic_sensitivity_results.json"
OUT_FIG = FIG_DIR / "fig21_harmonic_model_comparison.png"

SITE_CODES = ["ALT", "ZEP", "BRW", "CBA", "MHD", "KUM", "CGO", "SPO"]


def _fit_design(t, y, design):
    coeffs, _, _, _ = np.linalg.lstsq(design, np.asarray(y, dtype=float), rcond=None)
    return coeffs


def fit_annual(t, y):
    """Fit intercept + annual harmonic."""
    t = np.asarray(t, dtype=float)
    x1 = 2 * np.pi * t
    design = np.column_stack([np.ones_like(t), np.sin(x1), np.cos(x1)])
    c0, B1, C1 = _fit_design(t, y, design)
    return {
        "intercept": float(c0),
        "B": float(B1),
        "C": float(C1),
        "amplitude": float(np.hypot(B1, C1)),
    }


def fit_annual_plus_semiannual(t, y):
    """Fit intercept + annual + semiannual harmonics."""
    t = np.asarray(t, dtype=float)
    x1 = 2 * np.pi * t
    x2 = 4 * np.pi * t
    design = np.column_stack(
        [np.ones_like(t), np.sin(x1), np.cos(x1), np.sin(x2), np.cos(x2)]
    )
    c0, B1, C1, B2, C2 = _fit_design(t, y, design)
    return {
        "intercept": float(c0),
        "B1": float(B1),
        "C1": float(C1),
        "B2": float(B2),
        "C2": float(C2),
        "annual_amplitude": float(np.hypot(B1, C1)),
        "semiannual_amplitude": float(np.hypot(B2, C2)),
    }


def fit_monthly_fixed_effect(months, values):
    """Estimate seasonal amplitude from monthly means."""
    months = np.asarray(months, dtype=int)
    values = np.asarray(values, dtype=float)
    monthly = []
    for month in range(1, 13):
        vals = values[months == month]
        monthly.append(float(np.nanmean(vals)) if len(vals) else np.nan)
    arr = np.asarray(monthly, dtype=float)
    amplitude = 0.5 * (np.nanmax(arr) - np.nanmin(arr))
    return {
        "monthly_means": [float(v) if np.isfinite(v) else np.nan for v in arr],
        "amplitude": float(amplitude),
    }


def leave_one_year_out_ratios(records):
    """Return annual-harmonic ratios after excluding each year once.

    Records are tuples of (year, decimal_year, d13C, dD).
    """
    arr = list(records)
    years = sorted(set(int(r[0]) for r in arr))
    results = {}
    for held_out in years:
        kept = [r for r in arr if int(r[0]) != held_out]
        if len(kept) < 8:
            results[held_out] = np.nan
            continue
        t = np.array([r[1] for r in kept], dtype=float)
        y13 = np.array([r[2] for r in kept], dtype=float)
        yD = np.array([r[3] for r in kept], dtype=float)
        fit13 = fit_annual(t, y13)
        fitD = fit_annual(t, yD)
        results[held_out] = (
            float(fit13["amplitude"] / fitD["amplitude"])
            if fitD["amplitude"] > 0
            else np.nan
        )
    return results


def json_safe(value):
    """Recursively convert non-finite floats to None for strict JSON output."""
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    if isinstance(value, tuple):
        return [json_safe(item) for item in value]
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    if isinstance(value, np.floating):
        item = float(value)
        return item if np.isfinite(item) else None
    return value


def analyze_site(code):
    """Run all sensitivity models for one site."""
    df = pd.read_csv(PHASE1_DIR / f"site_monthly_{code}.csv")
    t = df["decimal_year"].to_numpy(dtype=float)
    months = df["month"].to_numpy(dtype=int)
    y13 = df["d13C_mean"].to_numpy(dtype=float)
    yD = df["dD_mean"].to_numpy(dtype=float)

    annual_13 = fit_annual(t, y13)
    annual_D = fit_annual(t, yD)
    semi_13 = fit_annual_plus_semiannual(t, y13)
    semi_D = fit_annual_plus_semiannual(t, yD)
    fixed_13 = fit_monthly_fixed_effect(months, y13)
    fixed_D = fit_monthly_fixed_effect(months, yD)
    records = list(zip(df["year"], df["decimal_year"], df["d13C_mean"], df["dD_mean"]))

    def ratio(a13, aD):
        return float(a13 / aD) if aD > 0 else np.nan

    loo = leave_one_year_out_ratios(records)
    finite_loo = np.array([v for v in loo.values() if np.isfinite(v)], dtype=float)
    return {
        "annual": {
            "R": ratio(annual_13["amplitude"], annual_D["amplitude"]),
            "A13": annual_13["amplitude"],
            "AD": annual_D["amplitude"],
        },
        "annual_plus_semiannual": {
            "R": ratio(semi_13["annual_amplitude"], semi_D["annual_amplitude"]),
            "A13": semi_13["annual_amplitude"],
            "AD": semi_D["annual_amplitude"],
            "semiannual_A13": semi_13["semiannual_amplitude"],
            "semiannual_AD": semi_D["semiannual_amplitude"],
        },
        "monthly_fixed_effect": {
            "R": ratio(fixed_13["amplitude"], fixed_D["amplitude"]),
            "A13": fixed_13["amplitude"],
            "AD": fixed_D["amplitude"],
        },
        "leave_one_year_out": {
            "ratios": {str(k): float(v) if np.isfinite(v) else np.nan for k, v in loo.items()},
            "median": float(np.nanmedian(finite_loo)) if len(finite_loo) else np.nan,
            "range": [
                float(np.nanmin(finite_loo)) if len(finite_loo) else np.nan,
                float(np.nanmax(finite_loo)) if len(finite_loo) else np.nan,
            ],
        },
    }


def plot_results(results):
    codes = list(results["sites"].keys())
    x = np.arange(len(codes))
    annual = [results["sites"][c]["annual"]["R"] for c in codes]
    semi = [results["sites"][c]["annual_plus_semiannual"]["R"] for c in codes]
    monthly = [results["sites"][c]["monthly_fixed_effect"]["R"] for c in codes]
    loo_median = [results["sites"][c]["leave_one_year_out"]["median"] for c in codes]
    loo_range = [results["sites"][c]["leave_one_year_out"]["range"] for c in codes]
    loo_low = np.array(
        [
            median - bounds[0] if np.isfinite(median) and np.isfinite(bounds[0]) else np.nan
            for median, bounds in zip(loo_median, loo_range)
        ],
        dtype=float,
    )
    loo_high = np.array(
        [
            bounds[1] - median if np.isfinite(median) and np.isfinite(bounds[1]) else np.nan
            for median, bounds in zip(loo_median, loo_range)
        ],
        dtype=float,
    )

    fig, ax = plt.subplots(figsize=(9.2, 4.4))
    offsets = {
        "annual": -0.27,
        "annual+semiannual": -0.09,
        "monthly fixed effect": 0.09,
        "leave-one-year-out": 0.27,
    }
    ax.scatter(x + offsets["annual"], annual, marker="o", s=42, label="Annual harmonic")
    ax.scatter(
        x + offsets["annual+semiannual"],
        semi,
        marker="s",
        s=40,
        label="Annual component of annual+semiannual fit",
    )
    ax.scatter(
        x + offsets["monthly fixed effect"],
        monthly,
        marker="^",
        s=44,
        label="Monthly fixed effect",
    )
    ax.errorbar(
        x + offsets["leave-one-year-out"],
        loo_median,
        yerr=[loo_low, loo_high],
        fmt="D",
        ms=4.5,
        capsize=3,
        label="Leave-one-year-out median and range",
    )
    ax.set_xticks(x)
    ax.set_xticklabels(codes)
    ax.set_ylabel("R = A(d13C) / A(dD)")
    ax.set_title("Harmonic model sensitivity")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=7, ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=180)
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    sites = {
        code: analyze_site(code)
        for code in SITE_CODES
        if (PHASE1_DIR / f"site_monthly_{code}.csv").exists()
    }
    output = {
        "metadata": {
            "method": "Compare annual harmonic, annual+semiannual, monthly fixed effects, and leave-one-year-out ratios"
        },
        "sites": sites,
    }
    OUT_JSON.write_text(json.dumps(json_safe(output), indent=2, allow_nan=False), encoding="utf-8")
    plot_results(output)
    print(f"Wrote {OUT_JSON}")
    print(f"Wrote {OUT_FIG}")


if __name__ == "__main__":
    main()
