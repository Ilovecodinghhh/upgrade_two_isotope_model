#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
phase3_synthesis.py — Cross-site synthesis for the KIE_sites experiment
========================================================================

Synthesizes the per-site harmonic fit results from Phase 2 to:

1. Classify sites as "clean" (small phase offset, remote MBL, adequate SNR)
   vs "contaminated" (large phase offset, continental source influence, or
   weak seasonal signal).

2. Compute a weighted-mean amplitude ratio R = A(δ¹³C)/A(δD) from the clean
   subset, using inverse-variance weights derived from bootstrap 95% CIs.

3. Generate publication-quality figures:
   - fig3_ratio_vs_latitude.png:  R vs latitude with clean/contaminated markers
   - fig3_phase_diagnostic.png:   phase difference vs latitude as a QC diagnostic
   - fig3_site_classification.png: overview of classification criteria

4. Interpret the observed ratios relative to pure-OH predictions.

Output:
    results/phase3_synthesis/
        synthesis_results.json   — classification, weighted mean, interpretation
    figures/
        fig3_ratio_vs_latitude.png
        fig3_phase_diagnostic.png
        fig3_site_classification.png
"""

from pathlib import Path
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================================
# PATHS
# ============================================================================
PHASE2_DIR = Path(__file__).resolve().parent.parent / "results" / "phase2_harmonics"
RESULTS_DIR = Path(__file__).resolve().parent.parent / "results" / "phase3_synthesis"
FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"
PHASE1_DIR = Path(__file__).resolve().parent.parent / "results" / "phase1_data"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# PURE-OH PREDICTIONS
# ============================================================================
# α_D for OH + CH₃D reaction (Saueressig et al. 2001, well-constrained)
ALPHA_D_OH = 1.31

# Contested α_13C for OH + ¹³CH₄ reaction:
ALPHA_13C_SAUERESSIG = 1.0039  # Saueressig et al. (2001)
ALPHA_13C_CANTRELL = 1.0054    # Cantrell et al. (1990)

# Predicted amplitude ratios under pure-OH seasonal control:
#   R_pure_OH = (α_13C − 1) / (α_D − 1)
R_OH_SAUERESSIG = (ALPHA_13C_SAUERESSIG - 1) / (ALPHA_D_OH - 1)  # ≈ 0.0126
R_OH_CANTRELL = (ALPHA_13C_CANTRELL - 1) / (ALPHA_D_OH - 1)      # ≈ 0.0174

# ============================================================================
# CLASSIFICATION CRITERIA
# ============================================================================
# A site is classified as "clean" if ALL of the following are met:
#   1. Phase offset |Δφ| ≤ PHASE_THRESHOLD months
#      → ensures δ¹³C and δD respond to the same seasonal driver
#   2. δ¹³C amplitude > MIN_D13C_AMPLITUDE ‰
#      → ensures adequate signal-to-noise (measurement precision ~0.05‰)
#   3. MBL (marine boundary layer) site = True
#      → reduces continental source contamination

PHASE_THRESHOLD = 2.0   # months — allows ±2 month tolerance
MIN_D13C_AMPLITUDE = 0.04  # ‰ — need at least ~4× measurement precision
MBL_REQUIRED = True

# MBL site flags from Phase 1 plan
MBL_SITES = {"ALT", "ZEP", "BRW", "CBA", "MHD", "KUM", "ASC", "SMO", "CGO", "SPO"}
NON_MBL_SITES = {"AZR", "MLO"}  # AZR marginal, MLO elevated non-MBL


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def ci_to_sigma(ci_lo: float, ci_hi: float) -> float:
    """Approximate 1σ uncertainty from a 95% CI, assuming normal distribution.

    For a 95% CI: width = 2 × 1.96 × σ  →  σ ≈ (ci_hi − ci_lo) / (2 × 1.96).
    """
    return (ci_hi - ci_lo) / (2 * 1.96)


def weighted_mean_with_uncertainty(values: np.ndarray,
                                   sigmas: np.ndarray) -> dict:
    """Compute inverse-variance weighted mean and its uncertainty.

    Parameters
    ----------
    values : array of point estimates
    sigmas : array of 1σ uncertainties

    Returns
    -------
    dict with: mean, sigma, chi2, n_sites, chi2_reduced
    """
    weights = 1.0 / sigmas**2
    w_sum = np.sum(weights)
    w_mean = np.sum(weights * values) / w_sum
    w_sigma = np.sqrt(1.0 / w_sum)

    # Chi-squared goodness of fit: do all sites agree within their errors?
    chi2 = np.sum(weights * (values - w_mean)**2)
    n = len(values)
    chi2_red = chi2 / (n - 1) if n > 1 else np.nan

    return {
        "mean": float(w_mean),
        "sigma": float(w_sigma),
        "chi2": float(chi2),
        "n_sites": int(n),
        "chi2_reduced": float(chi2_red),
    }


# ============================================================================
# SITE CLASSIFICATION
# ============================================================================

def classify_sites(fits: dict) -> dict:
    """Classify each site as clean or contaminated based on Phase 2 results.

    Returns dict: code → {classification, reasons, ...}
    """
    classifications = {}

    for code, data in fits.items():
        if "ratio" not in data:
            continue

        ratio = data["ratio"]["value"]
        ratio_ci = data["ratio"]["ci95"]
        phase_diff = data["phase_diff_months"]["value"]
        phase_ci = data["phase_diff_months"]["ci95"]
        amp_d13c = data["d13C"]["amplitude"]
        amp_dD = data["dD"]["amplitude"]

        is_mbl = code in MBL_SITES
        phase_ok = abs(phase_diff) <= PHASE_THRESHOLD
        snr_ok = amp_d13c >= MIN_D13C_AMPLITUDE

        reasons = []
        if not is_mbl:
            reasons.append(f"non-MBL site")
        if not phase_ok:
            reasons.append(f"|Δφ| = {abs(phase_diff):.1f} mo > {PHASE_THRESHOLD} mo")
        if not snr_ok:
            reasons.append(f"A(δ¹³C) = {amp_d13c:.4f}‰ < {MIN_D13C_AMPLITUDE}‰")

        is_clean = phase_ok and snr_ok and (is_mbl or not MBL_REQUIRED)

        classifications[code] = {
            "classification": "clean" if is_clean else "contaminated",
            "is_mbl": is_mbl,
            "phase_ok": phase_ok,
            "snr_ok": snr_ok,
            "reasons": reasons if reasons else ["all criteria met"],
            "ratio": ratio,
            "ratio_ci95": ratio_ci,
            "ratio_sigma": ci_to_sigma(ratio_ci[0], ratio_ci[1]),
            "phase_diff_months": phase_diff,
            "phase_diff_ci95": phase_ci,
            "amp_d13c": amp_d13c,
            "amp_dD": amp_dD,
            "latitude": _get_latitude(code, data),
        }

    return classifications


def _get_latitude(code: str, data: dict) -> float:
    """Extract latitude from Phase 1 summary or hardcoded table."""
    lat_table = {
        "ALT": 82.45, "ZEP": 78.91, "BRW": 71.32, "CBA": 55.21,
        "MHD": 53.33, "AZR": 38.77, "MLO": 19.54, "KUM": 19.56,
        "ASC": -7.97, "SMO": -14.25, "CGO": -40.68, "SPO": -89.98,
    }
    return lat_table.get(code, 0.0)


# ============================================================================
# INTERPRETATION
# ============================================================================

def interpret_ratios(classifications: dict, w_mean: dict) -> dict:
    """Provide scientific interpretation of the observed amplitude ratios.

    The key question: are the ratios consistent with pure-OH seasonal control?
    If not, what does the excess ratio imply?
    """
    r_obs = w_mean["mean"]
    r_obs_sigma = w_mean["sigma"]

    # How many σ above the pure-OH predictions?
    excess_saueressig = (r_obs - R_OH_SAUERESSIG) / r_obs_sigma
    excess_cantrell = (r_obs - R_OH_CANTRELL) / r_obs_sigma

    # The excess ratio R_obs >> R_OH implies an additional source of δ¹³C
    # seasonal variability beyond OH fractionation. The most likely candidate
    # is microbial source seasonality:
    #   - Wetlands and rice paddies emit CH₄ that is isotopically lighter
    #     (more negative δ¹³C ≈ −60‰, more negative δD ≈ −310‰)
    #   - These sources peak in summer (same season as OH maximum)
    #   - For δ¹³C: source depletion (~60‰ lighter than atmosphere) OPPOSES
    #     the OH enrichment, but the δ¹³C amplitude is still ENHANCED because
    #     the total seasonal swing involves both effects
    #   - For δD: source depletion (~230‰ lighter) has a relatively smaller
    #     FRACTIONAL effect on the seasonal amplitude compared to δ¹³C
    #
    # This asymmetry arises because:
    #   ε_OH_13C / (δ_source_13C − δ_atm_13C) ≠ ε_OH_D / (δ_source_D − δ_atm_D)
    #
    # The ratio of seasonal amplitudes thus reflects both KIE and source terms.

    interpretation = {
        "observed_weighted_mean_ratio": r_obs,
        "observed_sigma": r_obs_sigma,
        "OH_prediction_saueressig": R_OH_SAUERESSIG,
        "OH_prediction_cantrell": R_OH_CANTRELL,
        "sigma_above_saueressig": round(excess_saueressig, 1),
        "sigma_above_cantrell": round(excess_cantrell, 1),
        "conclusion": "",
    }

    if r_obs > R_OH_CANTRELL + 2 * r_obs_sigma:
        interpretation["conclusion"] = (
            "The observed amplitude ratio is significantly ABOVE the pure-OH "
            "prediction range, indicating that source seasonality (likely "
            "microbial emissions) contributes substantially to the δ¹³C seasonal "
            "cycle. The simple ratio method CANNOT directly constrain the OH KIE "
            "without source deconvolution (Phase 4)."
        )
    elif r_obs > R_OH_SAUERESSIG - 2 * r_obs_sigma:
        interpretation["conclusion"] = (
            "The observed ratio overlaps with or is near the pure-OH prediction "
            "range. Source deconvolution may still improve the constraint."
        )
    else:
        interpretation["conclusion"] = (
            "The observed ratio is BELOW the pure-OH prediction, suggesting "
            "anomalous suppression of the δ¹³C seasonal cycle — possibly due to "
            "source isotope signatures or interhemispheric transport effects."
        )

    return interpretation


# ============================================================================
# PLOTTING
# ============================================================================

def plot_ratio_vs_latitude(classifications: dict, w_mean_clean: dict,
                           w_mean_all: dict) -> None:
    """Figure 3a: Amplitude ratio vs latitude with clean/contaminated markers.

    Clean sites shown as filled circles, contaminated as open triangles.
    Horizontal bands show pure-OH predictions and weighted-mean from clean sites.
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    # Separate clean and contaminated sites
    for tag, marker, facecolor, label in [
        ("clean", "o", "C0", "Clean sites"),
        ("contaminated", "^", "none", "Excluded sites"),
    ]:
        codes = [c for c, v in classifications.items()
                 if v["classification"] == tag]
        if not codes:
            continue
        lats = [classifications[c]["latitude"] for c in codes]
        rats = [classifications[c]["ratio"] for c in codes]
        errs_lo = [abs(classifications[c]["ratio"] - classifications[c]["ratio_ci95"][0])
                   for c in codes]
        errs_hi = [abs(classifications[c]["ratio_ci95"][1] - classifications[c]["ratio"])
                   for c in codes]

        edgecolor = "C0" if tag == "clean" else "gray"
        ax.errorbar(lats, rats, yerr=[errs_lo, errs_hi],
                    fmt=marker, capsize=4, color=edgecolor, ms=8,
                    markerfacecolor=facecolor, markeredgecolor=edgecolor,
                    label=label, zorder=4)

        for c, lat, r in zip(codes, lats, rats):
            ax.annotate(c, (lat, r), textcoords="offset points",
                        xytext=(6, 6), fontsize=7,
                        color="black" if tag == "clean" else "gray")

    # Pure-OH prediction band
    ax.axhspan(R_OH_SAUERESSIG, R_OH_CANTRELL, alpha=0.15, color="green",
               label=f"Pure OH ({R_OH_SAUERESSIG:.4f}–{R_OH_CANTRELL:.4f})")
    ax.axhline(R_OH_SAUERESSIG, color="green", ls="--", lw=0.8, alpha=0.6)
    ax.axhline(R_OH_CANTRELL, color="green", ls="--", lw=0.8, alpha=0.6)

    # Weighted mean from clean sites
    if w_mean_clean["n_sites"] > 0:
        wm = w_mean_clean["mean"]
        ws = w_mean_clean["sigma"]
        ax.axhline(wm, color="C3", ls="-", lw=1.5, alpha=0.8,
                   label=f"Clean-site mean: {wm:.4f} ± {ws:.4f}")
        ax.axhspan(wm - 1.96*ws, wm + 1.96*ws, alpha=0.1, color="C3")

    ax.set_xlabel("Latitude (°)", fontsize=11)
    ax.set_ylabel("Amplitude ratio  R = A(δ¹³C) / A(δD)", fontsize=11)
    ax.set_title("Seasonal amplitude ratio vs latitude — site classification",
                 fontsize=12)
    ax.legend(fontsize=9, loc="upper left")
    ax.set_xlim(-100, 100)

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "fig3_ratio_vs_latitude.png", dpi=300)
    plt.close(fig)
    print(f"✓ Figure saved: {FIGURES_DIR / 'fig3_ratio_vs_latitude.png'}")


def plot_phase_diagnostic(classifications: dict) -> None:
    """Figure 3b: Phase difference vs latitude as a QC diagnostic.

    The phase threshold (±2 months) is shown as gray bands. Sites outside
    this band are flagged as contaminated.
    """
    fig, ax = plt.subplots(figsize=(10, 5))

    codes = sorted(classifications.keys(),
                   key=lambda c: -classifications[c]["latitude"])
    lats = [classifications[c]["latitude"] for c in codes]
    phases = [classifications[c]["phase_diff_months"] for c in codes]
    phase_lo = [abs(classifications[c]["phase_diff_months"] -
                    classifications[c]["phase_diff_ci95"][0]) for c in codes]
    phase_hi = [abs(classifications[c]["phase_diff_ci95"][1] -
                    classifications[c]["phase_diff_months"]) for c in codes]
    colors = ["C0" if classifications[c]["classification"] == "clean" else "gray"
              for c in codes]

    ax.errorbar(lats, phases, yerr=[phase_lo, phase_hi],
                fmt="none", capsize=3, ecolor="gray", zorder=1)
    ax.scatter(lats, phases, c=colors, s=60, zorder=3, edgecolors="k", lw=0.5)

    for c, lat, p in zip(codes, lats, phases):
        ax.annotate(c, (lat, p), textcoords="offset points",
                    xytext=(6, 6), fontsize=7)

    # Threshold bands
    ax.axhspan(-PHASE_THRESHOLD, PHASE_THRESHOLD, alpha=0.1, color="green",
               label=f"±{PHASE_THRESHOLD} month threshold")
    ax.axhline(0, color="gray", ls="--", lw=1)

    ax.set_xlabel("Latitude (°)", fontsize=11)
    ax.set_ylabel("Phase difference  δ¹³C − δD (months)", fontsize=11)
    ax.set_title("Phase alignment diagnostic — sites with large |Δφ| flagged",
                 fontsize=12)
    ax.legend(fontsize=9)
    ax.set_xlim(-100, 100)

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "fig3_phase_diagnostic.png", dpi=300)
    plt.close(fig)
    print(f"✓ Figure saved: {FIGURES_DIR / 'fig3_phase_diagnostic.png'}")


def plot_site_classification(classifications: dict) -> None:
    """Figure 3c: Overview table of classification criteria per site.

    A horizontal bar-style figure showing which criteria each site passes.
    """
    codes = sorted(classifications.keys(),
                   key=lambda c: -classifications[c]["latitude"])
    n = len(codes)

    fig, ax = plt.subplots(figsize=(10, max(4, 0.5 * n)))

    criteria_labels = ["MBL site", "|Δφ| ≤ 2 mo", "A(δ¹³C) ≥ 0.04‰"]
    cmap_pass = "#4CAF50"  # green
    cmap_fail = "#F44336"  # red

    for i, code in enumerate(codes):
        cl = classifications[code]
        checks = [cl["is_mbl"], cl["phase_ok"], cl["snr_ok"]]

        for j, passed in enumerate(checks):
            color = cmap_pass if passed else cmap_fail
            ax.barh(i, 0.8, left=j, height=0.6, color=color, edgecolor="white",
                    lw=0.5)
            symbol = "✓" if passed else "✗"
            ax.text(j + 0.4, i, symbol, ha="center", va="center",
                    fontsize=10, fontweight="bold",
                    color="white")

        # Site label with classification tag
        tag = cl["classification"].upper()
        tag_color = "darkgreen" if tag == "CLEAN" else "darkred"
        ax.text(-0.2, i, f"{code} ({cl['latitude']:+.0f}°)",
                ha="right", va="center", fontsize=9)
        ax.text(len(criteria_labels) + 0.2, i, tag,
                ha="left", va="center", fontsize=8, fontweight="bold",
                color=tag_color)

    ax.set_yticks([])
    ax.set_xticks([j + 0.4 for j in range(len(criteria_labels))])
    ax.set_xticklabels(criteria_labels, fontsize=9)
    ax.set_xlim(-0.5, len(criteria_labels) + 1.5)
    ax.set_ylim(-0.5, n - 0.5)
    ax.invert_yaxis()
    ax.set_title("Site classification criteria (green = pass, red = fail)",
                 fontsize=11)

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "fig3_site_classification.png", dpi=300)
    plt.close(fig)
    print(f"✓ Figure saved: {FIGURES_DIR / 'fig3_site_classification.png'}")


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("Phase 3: Cross-site synthesis")
    print("=" * 70)

    # ── Load Phase 2 results ──
    with open(PHASE2_DIR / "harmonic_fits.json") as f:
        fits = json.load(f)
    print(f"Loaded harmonic fits for {len(fits)} sites")

    # ── Classify sites ──
    classifications = classify_sites(fits)

    clean_codes = [c for c, v in classifications.items()
                   if v["classification"] == "clean"]
    contam_codes = [c for c, v in classifications.items()
                    if v["classification"] == "contaminated"]

    print(f"\n{'='*50}")
    print("SITE CLASSIFICATION")
    print(f"{'='*50}")
    print(f"\nClean sites ({len(clean_codes)}):")
    for c in sorted(clean_codes, key=lambda x: -classifications[x]["latitude"]):
        cl = classifications[c]
        print(f"  {c:>4} ({cl['latitude']:>+6.1f}°): R = {cl['ratio']:.4f} "
              f"[{cl['ratio_ci95'][0]:.4f}, {cl['ratio_ci95'][1]:.4f}], "
              f"Δφ = {cl['phase_diff_months']:+.1f} mo")

    print(f"\nExcluded sites ({len(contam_codes)}):")
    for c in sorted(contam_codes, key=lambda x: -classifications[x]["latitude"]):
        cl = classifications[c]
        print(f"  {c:>4} ({cl['latitude']:>+6.1f}°): R = {cl['ratio']:.4f}, "
              f"Δφ = {cl['phase_diff_months']:+.1f} mo  "
              f"— Reason: {'; '.join(cl['reasons'])}")

    # ── Compute weighted means ──
    # Clean sites only
    if clean_codes:
        clean_ratios = np.array([classifications[c]["ratio"] for c in clean_codes])
        clean_sigmas = np.array([classifications[c]["ratio_sigma"]
                                 for c in clean_codes])
        w_mean_clean = weighted_mean_with_uncertainty(clean_ratios, clean_sigmas)
    else:
        w_mean_clean = {"mean": np.nan, "sigma": np.nan, "chi2": np.nan,
                        "n_sites": 0, "chi2_reduced": np.nan}

    # All sites (for comparison)
    all_codes = list(classifications.keys())
    all_ratios = np.array([classifications[c]["ratio"] for c in all_codes])
    all_sigmas = np.array([classifications[c]["ratio_sigma"] for c in all_codes])
    w_mean_all = weighted_mean_with_uncertainty(all_ratios, all_sigmas)

    print(f"\n{'='*50}")
    print("WEIGHTED MEAN AMPLITUDE RATIO")
    print(f"{'='*50}")
    print(f"\n  Clean sites only ({w_mean_clean['n_sites']} sites):")
    print(f"    R = {w_mean_clean['mean']:.4f} ± {w_mean_clean['sigma']:.4f}")
    print(f"    χ²_red = {w_mean_clean['chi2_reduced']:.2f} "
          f"(χ² = {w_mean_clean['chi2']:.1f}, dof = {w_mean_clean['n_sites']-1})")
    print(f"\n  All sites ({w_mean_all['n_sites']} sites):")
    print(f"    R = {w_mean_all['mean']:.4f} ± {w_mean_all['sigma']:.4f}")
    print(f"    χ²_red = {w_mean_all['chi2_reduced']:.2f}")

    print(f"\n  Pure-OH predictions:")
    print(f"    Saueressig (α=1.0039): R = {R_OH_SAUERESSIG:.4f}")
    print(f"    Cantrell   (α=1.0054): R = {R_OH_CANTRELL:.4f}")

    # ── Interpret ──
    interpretation = interpret_ratios(classifications, w_mean_clean)

    print(f"\n{'='*50}")
    print("INTERPRETATION")
    print(f"{'='*50}")
    print(f"\n  Observed clean-site mean: {interpretation['observed_weighted_mean_ratio']:.4f}"
          f" ± {interpretation['observed_sigma']:.4f}")
    print(f"  {interpretation['sigma_above_saueressig']:.1f}σ above Saueressig prediction")
    print(f"  {interpretation['sigma_above_cantrell']:.1f}σ above Cantrell prediction")
    print(f"\n  {interpretation['conclusion']}")

    # ── Save results ──
    output = {
        "classification_criteria": {
            "phase_threshold_months": PHASE_THRESHOLD,
            "min_d13C_amplitude_permil": MIN_D13C_AMPLITUDE,
            "mbl_required": MBL_REQUIRED,
        },
        "site_classifications": classifications,
        "weighted_mean_clean": w_mean_clean,
        "weighted_mean_all": w_mean_all,
        "oh_predictions": {
            "saueressig_ratio": R_OH_SAUERESSIG,
            "cantrell_ratio": R_OH_CANTRELL,
            "alpha_D_OH": ALPHA_D_OH,
            "alpha_13C_saueressig": ALPHA_13C_SAUERESSIG,
            "alpha_13C_cantrell": ALPHA_13C_CANTRELL,
        },
        "interpretation": interpretation,
    }

    with open(RESULTS_DIR / "synthesis_results.json", "w") as f:
        json.dump(output, f, indent=2)
    print(f"\n✓ Results saved to {RESULTS_DIR / 'synthesis_results.json'}")

    # ── Figures ──
    plot_ratio_vs_latitude(classifications, w_mean_clean, w_mean_all)
    plot_phase_diagnostic(classifications)
    plot_site_classification(classifications)


if __name__ == "__main__":
    main()
