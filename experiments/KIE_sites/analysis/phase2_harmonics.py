#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
phase2_harmonics.py — Seasonal harmonic fitting for the KIE_sites experiment
=============================================================================

For each co-located site, fits an annual harmonic + linear trend to the
monthly-mean δ¹³C and δD time series from Phase 1:

    δ(t) = c₀ + c₁·(t − t_ref) + B·sin(2πt) + C·cos(2πt)

Extracts:
    - Amplitude  A = √(B² + C²)
    - Phase      φ = atan2(C, B)   →  peak month = (π/2 − φ) / (2π) × 12
    - Amplitude ratio  R = A_d13C / A_dD
    - Phase difference Δφ = φ_d13C − φ_dD  (months)

Bootstrap (N=2000) provides 95% CIs on all quantities.

Also loads CH₄ ppb event data and fits its seasonal cycle for reference
(the CH₄ cycle reflects the total source−sink balance without fractionation).

Output:
    results/phase2_harmonics/
        harmonic_fits.json     — per-site fit results + bootstrap CIs
    figures/
        fig2_seasonal_cycles.png     — detrended seasonal overlays per site
        fig2_harmonic_summary.png    — amplitude ratio + phase diff vs latitude
"""

from pathlib import Path
import json

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================================
# PATHS
# ============================================================================
REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
PHASE1_DIR = Path(__file__).resolve().parent.parent / "results" / "phase1_data"
RESULTS_DIR = Path(__file__).resolve().parent.parent / "results" / "phase2_harmonics"
FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"
CH4_EVENT_DIR = REPO_ROOT / "sitesdata" / "methane_ppb" / "noaa_gml_2025_event"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

N_BOOTSTRAP = 2000
RNG_SEED = 42

# Site metadata (must match Phase 1 order)
SITE_CODES = ["ALT", "ZEP", "BRW", "CBA", "MHD", "AZR",
              "MLO", "KUM", "ASC", "SMO", "CGO", "SPO"]
SITE_LAT = {"ALT": 82.45, "ZEP": 78.91, "BRW": 71.32, "CBA": 55.21,
            "MHD": 53.33, "AZR": 38.77, "MLO": 19.54, "KUM": 19.56,
            "ASC": -7.97, "SMO": -14.25, "CGO": -40.68, "SPO": -89.98}

# CH₄ ppb file code (lowercase, matching NOAA naming convention)
CH4_FILE_CODE = {c: c.lower() for c in SITE_CODES}


# ============================================================================
# HARMONIC FITTING
# ============================================================================

def fit_harmonic(t: np.ndarray, y: np.ndarray):
    """Fit annual harmonic + linear trend to a time series.

    Model: y(t) = c0 + c1*(t - t_ref) + B*sin(2πt) + C*cos(2πt)

    Parameters
    ----------
    t : decimal year array (e.g. 2005.04, 2005.12, ...)
    y : isotope values (‰) or CH₄ (ppb)

    Returns
    -------
    dict with keys:
        amplitude : √(B² + C²)
        phase_rad : atan2(C, B) — phase in radians
        peak_month : month of year when harmonic peaks (1–12, fractional)
        trend : c1 (slope per year)
        intercept : c0
        B, C : harmonic coefficients
        residuals : y - y_fit
        y_fit : fitted values
    """
    t_ref = np.mean(t)
    omega = 2.0 * np.pi  # annual frequency (period = 1 year)

    # Design matrix: [1, (t - t_ref), sin(ωt), cos(ωt)]
    X = np.column_stack([
        np.ones_like(t),
        t - t_ref,
        np.sin(omega * t),
        np.cos(omega * t),
    ])

    # Ordinary least squares
    coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    c0, c1, B, C = coeffs

    amplitude = np.sqrt(B**2 + C**2)
    phase_rad = np.arctan2(C, B)

    # Peak occurs where sin(ωt + φ) = 1, i.e. ωt + φ = π/2
    # → t_peak (fractional year) = (π/2 − φ) / ω
    t_peak_frac = (np.pi / 2 - phase_rad) / omega
    # Wrap to [0, 1) and convert to month (1 = Jan, 12 = Dec)
    t_peak_frac = t_peak_frac % 1.0
    peak_month = t_peak_frac * 12.0 + 1.0  # 1-indexed month
    if peak_month > 12.5:
        peak_month -= 12.0

    y_fit = X @ coeffs
    residuals = y - y_fit

    return {
        "amplitude": float(amplitude),
        "phase_rad": float(phase_rad),
        "peak_month": float(peak_month),
        "trend": float(c1),
        "intercept": float(c0),
        "B": float(B),
        "C": float(C),
        "residuals": residuals,
        "y_fit": y_fit,
    }


def bootstrap_harmonic(t: np.ndarray, y: np.ndarray, n_boot: int,
                       rng: np.random.Generator):
    """Bootstrap the harmonic fit to get confidence intervals.

    Resamples (t, y) pairs with replacement and refits each time.

    Returns
    -------
    dict with arrays of length n_boot:
        amplitude, phase_rad, peak_month, trend
    """
    n = len(t)
    amps = np.zeros(n_boot)
    phases = np.zeros(n_boot)
    peaks = np.zeros(n_boot)
    trends = np.zeros(n_boot)

    for b in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        fit = fit_harmonic(t[idx], y[idx])
        amps[b] = fit["amplitude"]
        phases[b] = fit["phase_rad"]
        peaks[b] = fit["peak_month"]
        trends[b] = fit["trend"]

    return {"amplitude": amps, "phase_rad": phases,
            "peak_month": peaks, "trend": trends}


def ci_95(arr: np.ndarray):
    """Return (median, lower 2.5%, upper 97.5%) from bootstrap samples."""
    return (float(np.nanmedian(arr)),
            float(np.nanpercentile(arr, 2.5)),
            float(np.nanpercentile(arr, 97.5)))


def month_name(m: float) -> str:
    """Convert fractional month (1–12) to abbreviated name."""
    names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
             "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    idx = int(round(m - 1)) % 12
    return names[idx]


# ============================================================================
# CH₄ CONCENTRATION LOADING
# ============================================================================

def load_ch4_event(site_code_lower: str, year_min: int, year_max: int) -> pd.DataFrame:
    """Load NOAA GML CH₄ ppb event data, QC-filtered and trimmed to year range.

    Returns DataFrame with columns: decimal_year, year, month, ch4_ppb.
    Returns empty DataFrame if file not found.
    """
    pattern = f"ch4_{site_code_lower}_*_event.txt"
    matches = list(CH4_EVENT_DIR.glob(pattern))
    if not matches:
        return pd.DataFrame()
    fpath = matches[0]

    with open(fpath) as f:
        first_line = f.readline()
    n_header = int(first_line.split(":")[-1].strip())

    df = pd.read_csv(
        fpath, sep=r"\s+", skiprows=n_header,
        names=["site", "year", "month", "day", "hour", "minute", "second",
               "datetime", "decimal_year", "flask_id", "value", "value_unc",
               "lat", "lon", "alt", "elev", "intake_ht", "method",
               "event_num", "instrument", "analysis_dt", "qcflag"],
        dtype={"qcflag": str},
        na_values=["-999.999", "-999.99", "-9.999"],
    )
    df = df[df["qcflag"].str[0] == "."].copy()
    df = df.rename(columns={"value": "ch4_ppb"})
    df = df[(df["year"] >= year_min) & (df["year"] <= year_max)]
    df = df[["decimal_year", "year", "month", "ch4_ppb"]].dropna(subset=["ch4_ppb"])

    # Monthly means
    monthly = (df.groupby(["year", "month"])
               .agg(ch4_mean=("ch4_ppb", "mean"))
               .reset_index())
    monthly["decimal_year"] = monthly["year"] + (monthly["month"] - 0.5) / 12.0
    return monthly.sort_values("decimal_year").reset_index(drop=True)


# ============================================================================
# PLOTTING
# ============================================================================

def plot_seasonal_cycles(all_results: dict) -> None:
    """Figure: Detrended seasonal overlays for δ¹³C, δD, and CH₄ at each site.

    For each site, subtracts the linear trend and plots residuals vs month-of-year,
    with the fitted harmonic curve overlaid.
    """
    codes = [c for c in SITE_CODES if c in all_results]
    n = len(codes)
    fig, axes = plt.subplots(n, 3, figsize=(15, 2.2 * n), squeeze=False)

    col_titles = ["δ¹³C-CH₄ (‰ VPDB)", "δD-CH₄ (‰ VSMOW)", "CH₄ (ppb)"]
    for j, title in enumerate(col_titles):
        axes[0, j].annotate(title, xy=(0.5, 1.35), xycoords="axes fraction",
                            ha="center", fontsize=11, fontweight="bold")

    # Distinct colors for each data type: δ¹³C = blue, δD = orange, CH₄ = green
    DATA_COLORS = {"d13C": "C0", "dD": "C1", "ch4": "C2"}
    FIT_COLORS = {"d13C": "darkblue", "dD": "darkred", "ch4": "darkgreen"}

    for i, code in enumerate(codes):
        res = all_results[code]
        lat = SITE_LAT[code]

        for j, key in enumerate(["d13C", "dD", "ch4"]):
            ax = axes[i, j]
            if key not in res:
                ax.set_visible(False)
                continue

            r = res[key]
            t = np.array(r["_t"])
            y = np.array(r["_y"])

            # Detrend: subtract linear component, keep harmonic + residual
            y_detrended = y - (r["intercept"] + r["trend_per_yr"] * (t - np.mean(t)))
            # Fractional month for x-axis (1 = Jan, 12 = Dec)
            month_frac = (t % 1.0) * 12.0 + 1.0

            ax.scatter(month_frac, y_detrended, s=12, alpha=0.6,
                       color=DATA_COLORS[key], zorder=2)

            # Overlay fitted harmonic curve
            t_smooth = np.linspace(0, 1, 200)
            omega = 2.0 * np.pi
            y_harmonic = r["B"] * np.sin(omega * t_smooth) + \
                         r["C"] * np.cos(omega * t_smooth)
            ax.plot(t_smooth * 12 + 1, y_harmonic, "-", lw=1.5,
                    color=FIT_COLORS[key], zorder=3)

            ax.axhline(0, color="gray", lw=0.5, ls="--")
            ax.set_xlim(0.5, 12.5)
            ax.set_xticks(range(1, 13))
            ax.tick_params(labelsize=6)

            if j == 0:
                ax.set_ylabel(f"{code}\n({lat:+.0f}°)", fontsize=8)

            if i == n - 1:
                ax.set_xticklabels(["J", "F", "M", "A", "M", "J",
                                    "J", "A", "S", "O", "N", "D"], fontsize=6)
                ax.set_xlabel("Month", fontsize=8)
            else:
                ax.set_xticklabels([])

            # Annotate amplitude and peak month
            amp_str = f"A={r['amplitude']:.3f}" if key != "ch4" else f"A={r['amplitude']:.1f}"
            ax.annotate(f"{amp_str}, pk={month_name(r['peak_month'])}",
                        xy=(0.98, 0.92), xycoords="axes fraction",
                        ha="right", fontsize=6, color=FIT_COLORS[key])

    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(FIGURES_DIR / "fig2_seasonal_cycles.png", dpi=300)
    plt.close(fig)
    print(f"✓ Figure saved: {FIGURES_DIR / 'fig2_seasonal_cycles.png'}")


def plot_harmonic_summary(all_results: dict) -> None:
    """Figure: Amplitude ratio and phase difference vs latitude, with 95% CIs."""
    codes = [c for c in SITE_CODES
             if c in all_results and "ratio" in all_results[c]]

    lats = [SITE_LAT[c] for c in codes]
    ratios = [all_results[c]["ratio"]["value"] for c in codes]
    ratio_lo = [all_results[c]["ratio"]["ci95"][0] for c in codes]
    ratio_hi = [all_results[c]["ratio"]["ci95"][1] for c in codes]
    phase_diffs = [all_results[c]["phase_diff_months"]["value"] for c in codes]
    phase_lo = [all_results[c]["phase_diff_months"]["ci95"][0] for c in codes]
    phase_hi = [all_results[c]["phase_diff_months"]["ci95"][1] for c in codes]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # ── Panel (a): Amplitude ratio vs latitude ──
    # Asymmetric errorbars from bootstrap CIs (use abs to guard against edge cases)
    ratio_err_lo = np.abs(np.array(ratios) - np.array(ratio_lo))
    ratio_err_hi = np.abs(np.array(ratio_hi) - np.array(ratios))
    ax1.errorbar(lats, ratios,
                 yerr=[ratio_err_lo, ratio_err_hi],
                 fmt="o", capsize=4, color="C0", ms=6)
    for c, lat, r in zip(codes, lats, ratios):
        ax1.annotate(c, (lat, r), textcoords="offset points",
                     xytext=(5, 5), fontsize=7)

    # Predicted ratio bands from OH KIE
    # Saueressig: (1.0039−1)/(1.31−1) = 0.0126;  Cantrell: (1.0054−1)/(1.31−1) = 0.0174
    ax1.axhspan(0.0126, 0.0174, alpha=0.15, color="green",
                label="Pure OH prediction\n(Saueressig–Cantrell)")
    ax1.axhline(0.0126, color="green", ls="--", lw=0.8, alpha=0.5)
    ax1.axhline(0.0174, color="green", ls="--", lw=0.8, alpha=0.5)

    ax1.set_xlabel("Latitude (°)", fontsize=10)
    ax1.set_ylabel("Amplitude ratio  A(δ¹³C) / A(δD)", fontsize=10)
    ax1.set_title("(a) Seasonal amplitude ratio vs latitude", fontsize=11)
    ax1.legend(fontsize=8, loc="upper right")

    # ── Panel (b): Phase difference vs latitude ──
    phase_err_lo = np.abs(np.array(phase_diffs) - np.array(phase_lo))
    phase_err_hi = np.abs(np.array(phase_hi) - np.array(phase_diffs))
    ax2.errorbar(lats, phase_diffs,
                 yerr=[phase_err_lo, phase_err_hi],
                 fmt="s", capsize=4, color="C1", ms=6)
    for c, lat, p in zip(codes, lats, phase_diffs):
        ax2.annotate(c, (lat, p), textcoords="offset points",
                     xytext=(5, 5), fontsize=7)

    ax2.axhline(0, color="gray", ls="--", lw=1, label="Perfect alignment")
    ax2.axhspan(-1, 1, alpha=0.1, color="gray", label="±1 month")

    ax2.set_xlabel("Latitude (°)", fontsize=10)
    ax2.set_ylabel("Phase difference  δ¹³C − δD (months)", fontsize=10)
    ax2.set_title("(b) Phase difference vs latitude", fontsize=11)
    ax2.legend(fontsize=8, loc="upper right")

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "fig2_harmonic_summary.png", dpi=300)
    plt.close(fig)
    print(f"✓ Figure saved: {FIGURES_DIR / 'fig2_harmonic_summary.png'}")


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("Phase 2: Seasonal harmonic fitting")
    print("=" * 70)

    rng = np.random.default_rng(RNG_SEED)
    all_results = {}

    for code in SITE_CODES:
        csv_path = PHASE1_DIR / f"site_monthly_{code}.csv"
        if not csv_path.exists():
            print(f"\n--- {code}: no Phase 1 data, skipping ---")
            continue

        monthly = pd.read_csv(csv_path)
        if len(monthly) < 12:
            print(f"\n--- {code}: only {len(monthly)} months, need ≥12 for harmonic fit ---")
            continue

        lat = SITE_LAT[code]
        print(f"\n--- {code} ({lat:+.1f}°, {len(monthly)} months) ---")

        results = {}

        # ── Fit δ¹³C and δD ──
        for key, col in [("d13C", "d13C_mean"), ("dD", "dD_mean")]:
            t = monthly["decimal_year"].values
            y = monthly[col].values

            fit = fit_harmonic(t, y)
            boot = bootstrap_harmonic(t, y, N_BOOTSTRAP, rng)

            amp_med, amp_lo, amp_hi = ci_95(boot["amplitude"])
            peak_med, peak_lo, peak_hi = ci_95(boot["peak_month"])

            results[key] = {
                "amplitude": fit["amplitude"],
                "amplitude_ci95": [amp_lo, amp_hi],
                "peak_month": fit["peak_month"],
                "peak_month_ci95": [peak_lo, peak_hi],
                "trend_per_yr": fit["trend"],
                "phase_rad": fit["phase_rad"],
                "B": fit["B"],
                "C": fit["C"],
                "intercept": fit["intercept"],
                # Store raw data for plotting (not saved to JSON)
                "_t": t.tolist(),
                "_y": y.tolist(),
            }

            unit = "‰" if key != "ch4" else "ppb"
            print(f"  {key}: A = {fit['amplitude']:.4f}{unit}  "
                  f"[{amp_lo:.4f}, {amp_hi:.4f}],  "
                  f"peak = {month_name(fit['peak_month'])} "
                  f"({fit['peak_month']:.1f})")

        # ── Fit CH₄ ppb (for reference) ──
        year_min = int(monthly["year"].min())
        year_max = int(monthly["year"].max())
        ch4_monthly = load_ch4_event(CH4_FILE_CODE[code], year_min, year_max)

        if len(ch4_monthly) >= 12:
            t_ch4 = ch4_monthly["decimal_year"].values
            y_ch4 = ch4_monthly["ch4_mean"].values
            fit_ch4 = fit_harmonic(t_ch4, y_ch4)
            boot_ch4 = bootstrap_harmonic(t_ch4, y_ch4, N_BOOTSTRAP, rng)
            amp_med_ch4, amp_lo_ch4, amp_hi_ch4 = ci_95(boot_ch4["amplitude"])

            results["ch4"] = {
                "amplitude": fit_ch4["amplitude"],
                "amplitude_ci95": [amp_lo_ch4, amp_hi_ch4],
                "peak_month": fit_ch4["peak_month"],
                "trend_per_yr": fit_ch4["trend"],
                "phase_rad": fit_ch4["phase_rad"],
                "B": fit_ch4["B"],
                "C": fit_ch4["C"],
                "intercept": fit_ch4["intercept"],
                "_t": t_ch4.tolist(),
                "_y": y_ch4.tolist(),
            }
            print(f"  CH₄: A = {fit_ch4['amplitude']:.1f} ppb,  "
                  f"peak = {month_name(fit_ch4['peak_month'])}")
        else:
            print(f"  CH₄: insufficient data ({len(ch4_monthly)} months)")

        # ── Compute amplitude ratio and phase difference ──
        A_c13 = results["d13C"]["amplitude"]
        A_dD = results["dD"]["amplitude"]
        ratio = A_c13 / A_dD if A_dD > 0 else np.nan

        # Phase difference in months: (φ_d13C − φ_dD) converted from radians
        # One full cycle = 2π rad = 12 months
        phase_diff_rad = results["d13C"]["phase_rad"] - results["dD"]["phase_rad"]
        # Wrap to [−π, π]
        phase_diff_rad = (phase_diff_rad + np.pi) % (2 * np.pi) - np.pi
        phase_diff_months = phase_diff_rad / (2 * np.pi) * 12.0

        # Bootstrap the ratio and phase difference using PAIRED resampling.
        # Both isotopes are resampled with the SAME random indices in each
        # bootstrap iteration.  This preserves the month-to-month correlation
        # between δ¹³C and δD (measured from the same flasks at INSTAAR sites),
        # which is essential for a correct CI on their ratio and phase difference.
        # Independent resampling would overestimate the uncertainty because it
        # breaks the natural pairing structure.
        t_iso = monthly["decimal_year"].values
        y_c13 = monthly["d13C_mean"].values
        y_dD = monthly["dD_mean"].values
        n_pts = len(t_iso)

        boot_ratios = np.zeros(N_BOOTSTRAP)
        boot_phase_diffs = np.zeros(N_BOOTSTRAP)
        for b in range(N_BOOTSTRAP):
            idx = rng.choice(n_pts, size=n_pts, replace=True)
            fit_c = fit_harmonic(t_iso[idx], y_c13[idx])
            fit_d = fit_harmonic(t_iso[idx], y_dD[idx])
            boot_ratios[b] = fit_c["amplitude"] / fit_d["amplitude"] if fit_d["amplitude"] > 0 else np.nan
            pd_rad = fit_c["phase_rad"] - fit_d["phase_rad"]
            pd_rad = (pd_rad + np.pi) % (2 * np.pi) - np.pi
            boot_phase_diffs[b] = pd_rad / (2 * np.pi) * 12.0

        ratio_med, ratio_lo, ratio_hi = ci_95(boot_ratios)
        pd_med, pd_lo, pd_hi = ci_95(boot_phase_diffs)

        results["ratio"] = {
            "value": ratio,
            "ci95": [ratio_lo, ratio_hi],
        }
        results["phase_diff_months"] = {
            "value": phase_diff_months,
            "ci95": [pd_lo, pd_hi],
        }

        print(f"  ── Ratio A(δ¹³C)/A(δD) = {ratio:.4f}  [{ratio_lo:.4f}, {ratio_hi:.4f}]")
        print(f"  ── Phase diff (δ¹³C − δD) = {phase_diff_months:+.1f} months  "
              f"[{pd_lo:+.1f}, {pd_hi:+.1f}]")

        # Comparison with pure-OH prediction
        oh_saueressig = (1.0039 - 1) / (1.31 - 1)  # 0.0126
        oh_cantrell = (1.0054 - 1) / (1.31 - 1)     # 0.0174
        if ratio < oh_saueressig:
            print(f"       → BELOW pure-OH range ({oh_saueressig:.4f}–{oh_cantrell:.4f})")
        elif ratio > oh_cantrell:
            print(f"       → ABOVE pure-OH range ({oh_saueressig:.4f}–{oh_cantrell:.4f})")
        else:
            print(f"       → WITHIN pure-OH range ({oh_saueressig:.4f}–{oh_cantrell:.4f})")

        all_results[code] = results

    # ── Save results (strip raw data arrays for JSON) ──
    json_results = {}
    for code, res in all_results.items():
        jr = {}
        for key, val in res.items():
            if isinstance(val, dict):
                jr[key] = {k: v for k, v in val.items() if not k.startswith("_")}
            else:
                jr[key] = val
        json_results[code] = jr

    with open(RESULTS_DIR / "harmonic_fits.json", "w") as f:
        json.dump(json_results, f, indent=2)
    print(f"\n✓ Results saved to {RESULTS_DIR / 'harmonic_fits.json'}")

    # ── Figures ──
    plot_seasonal_cycles(all_results)
    plot_harmonic_summary(all_results)

    # ── Summary table ──
    print("\n" + "=" * 70)
    print("HARMONIC FIT SUMMARY")
    print("=" * 70)
    print(f"{'Site':<5} {'Lat':>6} {'A(δ¹³C)/‰':>10} {'pk':>4} "
          f"{'A(δD)/‰':>9} {'pk':>4} {'Ratio':>8} {'Δφ(mo)':>8} {'vs OH?':>10}")
    print("-" * 74)
    oh_lo = (1.0039 - 1) / (1.31 - 1)
    oh_hi = (1.0054 - 1) / (1.31 - 1)
    for code in SITE_CODES:
        if code not in all_results or "ratio" not in all_results[code]:
            continue
        r = all_results[code]
        ratio = r["ratio"]["value"]
        if ratio < oh_lo:
            verdict = "BELOW"
        elif ratio > oh_hi:
            verdict = "ABOVE"
        else:
            verdict = "WITHIN"
        print(f"{code:<5} {SITE_LAT[code]:>+6.1f} "
              f"{r['d13C']['amplitude']:>10.4f} {month_name(r['d13C']['peak_month']):>4} "
              f"{r['dD']['amplitude']:>9.3f} {month_name(r['dD']['peak_month']):>4} "
              f"{ratio:>8.4f} {r['phase_diff_months']['value']:>+7.1f} "
              f"{verdict:>10}")


if __name__ == "__main__":
    main()
