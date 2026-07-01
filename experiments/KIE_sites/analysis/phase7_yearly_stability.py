#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Individual-year seasonal fits for KIE_sites.

This robustness check asks whether annual amplitude ratios are stable across
the short 2005-2010 co-located isotope window. Unstable year-to-year ratios
would imply that the full-record harmonic uncertainty may be too optimistic.
"""

from pathlib import Path
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


EXPT_DIR = Path(__file__).resolve().parent.parent
PHASE1_DIR = EXPT_DIR / "results" / "phase1_data"
PHASE2_JSON = EXPT_DIR / "results" / "phase2_harmonics" / "harmonic_fits.json"
OUT_DIR = EXPT_DIR / "results" / "phase7_yearly_stability"
FIG_DIR = EXPT_DIR / "figures"

OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

OUT_JSON = OUT_DIR / "yearly_stability_results.json"
OUT_FIG = FIG_DIR / "fig12_yearly_stability.png"

SITE_ORDER = ["ALT", "ZEP", "BRW", "CBA", "MHD", "AZR",
              "MLO", "KUM", "ASC", "SMO", "CGO", "SPO"]
SITE_LATS = {
    "ALT": 82.45, "ZEP": 78.91, "BRW": 71.32, "CBA": 55.21,
    "MHD": 53.33, "AZR": 38.77, "MLO": 19.54, "KUM": 19.56,
    "ASC": -7.97, "SMO": -14.25, "CGO": -40.68, "SPO": -89.98,
}
MIN_MONTHS = 8
MIN_USABLE_YEARS_FOR_STABILITY = 2


def fit_harmonic(t, y):
    """Fit annual harmonic plus linear trend."""
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    t_ref = np.mean(t)
    omega = 2 * np.pi
    X = np.column_stack([
        np.ones_like(t),
        t - t_ref,
        np.sin(omega * t),
        np.cos(omega * t),
    ])
    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    c0, c1, B, C = coeffs
    amplitude = float(np.hypot(B, C))
    phase_rad = float(np.arctan2(C, B))
    peak_frac = ((np.pi / 2 - phase_rad) / omega) % 1.0
    peak_month = float(peak_frac * 12.0 + 1.0)
    if peak_month > 12.5:
        peak_month -= 12.0
    return {
        "amplitude": amplitude,
        "phase_rad": phase_rad,
        "peak_month": peak_month,
        "trend": float(c1),
        "intercept": float(c0),
        "B": float(B),
        "C": float(C),
    }


def phase_diff_months(phase_13c, phase_dD):
    """Return wrapped phase difference in months."""
    diff = phase_13c - phase_dD
    diff = (diff + np.pi) % (2 * np.pi) - np.pi
    return float(diff / (2 * np.pi) * 12.0)


def fit_site_year(site_code, year, df, min_months=MIN_MONTHS):
    """Fit one site-year if enough paired monthly observations exist."""
    yearly = df[df["year"] == year].sort_values("month")
    if len(yearly) < min_months:
        return None

    t = yearly["decimal_year"].to_numpy()
    fit_13c = fit_harmonic(t, yearly["d13C_mean"].to_numpy())
    fit_dD = fit_harmonic(t, yearly["dD_mean"].to_numpy())
    amp_dD = fit_dD["amplitude"]
    ratio = fit_13c["amplitude"] / amp_dD if amp_dD > 0 else np.nan

    return {
        "site": site_code,
        "year": int(year),
        "n_months": int(len(yearly)),
        "ratio": float(ratio),
        "phase_diff_months": phase_diff_months(fit_13c["phase_rad"], fit_dD["phase_rad"]),
        "d13C": {k: fit_13c[k] for k in ("amplitude", "peak_month", "B", "C")},
        "dD": {k: fit_dD[k] for k in ("amplitude", "peak_month", "B", "C")},
    }


def ci_to_sigma(ci):
    return (ci[1] - ci[0]) / (2 * 1.96)


def summarize_site(code, yearly_records, phase2):
    """Summarize yearly scatter relative to full-period harmonic ratio."""
    full_ratio = phase2[code]["ratio"]["value"]
    full_ci = phase2[code]["ratio"]["ci95"]
    full_sigma = ci_to_sigma(full_ci)

    ratios = np.array([r["ratio"] for r in yearly_records], dtype=float)
    ratios = ratios[np.isfinite(ratios)]
    n = len(ratios)
    if n == 0:
        return {
            "n_years": 0,
            "full_period_ratio": float(full_ratio),
            "full_period_ratio_ci95": full_ci,
        }

    mean = float(np.mean(ratios))
    sd = float(np.std(ratios, ddof=1)) if n >= 2 else np.nan
    rms_vs_full = float(np.sqrt(np.mean((ratios - full_ratio) ** 2)))
    scatter_to_phase2_sigma = float(sd / full_sigma) if n >= 2 and full_sigma > 0 else np.nan
    within_ci = float(np.mean((ratios >= full_ci[0]) & (ratios <= full_ci[1])))

    return {
        "n_years": int(n),
        "full_period_ratio": float(full_ratio),
        "full_period_ratio_ci95": [float(full_ci[0]), float(full_ci[1])],
        "yearly_ratio_mean": mean,
        "yearly_ratio_sd": sd,
        "rms_difference_from_full_period": rms_vs_full,
        "scatter_to_phase2_sigma": scatter_to_phase2_sigma,
        "fraction_years_within_phase2_ci": within_ci,
        "interpretation": classify_stability(n, scatter_to_phase2_sigma, within_ci),
    }


def classify_stability(n_years, scatter_to_phase2_sigma, within_ci):
    if n_years < MIN_USABLE_YEARS_FOR_STABILITY:
        return "insufficient_years"
    if np.isfinite(scatter_to_phase2_sigma) and scatter_to_phase2_sigma > 1.5:
        return "yearly_scatter_exceeds_phase2_uncertainty"
    if within_ci < 0.5:
        return "most_years_outside_phase2_ci"
    return "yearly_ratios_broadly_consistent"


def run_yearly_stability(min_months=MIN_MONTHS):
    with open(PHASE2_JSON, encoding="utf-8") as f:
        phase2 = json.load(f)

    results = {
        "metadata": {
            "method": "Individual-year annual harmonic fits",
            "min_months_per_year": min_months,
            "min_usable_years_for_stability": MIN_USABLE_YEARS_FOR_STABILITY,
            "note": "Diagnostic only: many site-years have sparse monthly coverage.",
        },
        "sites": {},
        "summary": {},
    }

    for csv_path in sorted(PHASE1_DIR.glob("site_monthly_*.csv")):
        code = csv_path.stem.replace("site_monthly_", "")
        if code not in phase2:
            continue
        df = pd.read_csv(csv_path)
        yearly_records = []
        for year in sorted(df["year"].unique()):
            rec = fit_site_year(code, int(year), df, min_months=min_months)
            if rec is not None:
                yearly_records.append(rec)

        results["sites"][code] = {
            "latitude": SITE_LATS.get(code),
            "yearly_fits": yearly_records,
            "stability": summarize_site(code, yearly_records, phase2),
        }

    summaries = [s["stability"] for s in results["sites"].values()]
    usable = [s for s in summaries if s.get("n_years", 0) >= 2]
    flagged = [
        code for code, site in results["sites"].items()
        if site["stability"].get("interpretation") in (
            "yearly_scatter_exceeds_phase2_uncertainty",
            "most_years_outside_phase2_ci",
        )
    ]
    results["summary"] = {
        "sites_with_at_least_two_usable_years": len(usable),
        "sites_flagged_unstable": flagged,
        "median_scatter_to_phase2_sigma": (
            float(np.nanmedian([s["scatter_to_phase2_sigma"] for s in usable]))
            if usable else np.nan
        ),
    }
    return results


def plot_yearly_stability(results):
    """Create yearly R stability by site."""
    fig, axes = plt.subplots(4, 3, figsize=(13, 11), sharex=True)
    axes = axes.flatten()
    sc = None

    for ax, code in zip(axes, SITE_ORDER):
        site = results["sites"].get(code)
        if not site:
            ax.axis("off")
            continue
        stability = site["stability"]
        yearly = site["yearly_fits"]
        full = stability["full_period_ratio"]
        ci = stability["full_period_ratio_ci95"]

        ax.axhspan(ci[0], ci[1], color="0.7", alpha=0.22, lw=0)
        ax.axhline(full, color="0.25", lw=1.3, ls="--")

        if yearly:
            years = [r["year"] for r in yearly]
            ratios = [r["ratio"] for r in yearly]
            n_months = [r["n_months"] for r in yearly]
            sc = ax.scatter(
                years, ratios, c=n_months, cmap="viridis", vmin=8, vmax=12,
                s=45, edgecolor="black", linewidth=0.4, zorder=4,
            )
            ax.plot(years, ratios, color="C0", alpha=0.45, lw=1.0)

        interp = stability.get("interpretation", "no_fit")
        color = "C3" if "exceeds" in interp or "outside" in interp else "black"
        ax.set_title(
            f"{code} ({SITE_LATS.get(code, 0):+.0f})  n={stability.get('n_years', 0)}",
            fontsize=9,
            color=color,
        )
        ax.set_xlim(2004.6, 2010.4)
        ax.set_xticks(range(2005, 2011))
        ax.tick_params(axis="x", labelrotation=45, labelsize=7)
        ax.tick_params(axis="y", labelsize=7)
        ax.grid(alpha=0.2)

    for ax in axes[::3]:
        ax.set_ylabel("Yearly R = A(d13C)/A(dD)", fontsize=8)

    fig.subplots_adjust(left=0.07, right=0.90, top=0.91, bottom=0.14,
                        hspace=0.42, wspace=0.22)
    cax = fig.add_axes([0.925, 0.18, 0.018, 0.62])
    cbar = fig.colorbar(sc, cax=cax)
    cbar.set_label("Paired months in plotted year (>=8)", fontsize=9)

    fig.suptitle(
        "Individual-year seasonal amplitude ratios (2005-2010)",
        fontsize=13,
        fontweight="bold",
    )
    fig.text(
        0.5, 0.055,
        "Only site-years with >=8 paired months are plotted; >=2 plotted years are required for stability classification.",
        ha="center",
        fontsize=9,
    )
    fig.text(
        0.5, 0.030,
        "Dashed line = full-period annual-harmonic ratio; grey band = full-period 95% CI. "
        "Red titles flag yearly scatter larger than full-period uncertainty.",
        ha="center",
        fontsize=9,
    )
    fig.savefig(OUT_FIG, dpi=300)
    plt.close(fig)


def main():
    results = run_yearly_stability()
    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    plot_yearly_stability(results)

    print(f"Saved: {OUT_JSON}")
    print(f"Saved: {OUT_FIG}")
    print("Sites with >=2 usable years:", results["summary"]["sites_with_at_least_two_usable_years"])
    print("Flagged unstable sites:", results["summary"]["sites_flagged_unstable"])
    print("Median yearly scatter / full-period sigma:",
          f"{results['summary']['median_scatter_to_phase2_sigma']:.2f}")


if __name__ == "__main__":
    main()
