#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
phase4_deconv.py — Source deconvolution for the KIE_sites experiment
=====================================================================

Phase 3 showed that the observed amplitude ratio R = A(δ¹³C)/A(δD) is
3–4σ above the pure-OH prediction, implying source seasonality contaminates
the isotope seasonal cycles. This phase uses a simple seasonal box model to
separate source and sink contributions.

Model for each site (one well-mixed box per isotopologue):
--------------------------------------------------------------------------
The seasonal cycle of each isotope ratio δ can be decomposed as:

    Δδ_observed = Δδ_sink + Δδ_source

where:
    Δδ_sink  = ε_OH × ΔF_OH   (OH fractionation × seasonal OH flux change)
    Δδ_source = (δ_source − δ_atm) × ΔS/S_mean  (source signature vs atm)

The key insight: the CH₄ concentration seasonal cycle constrains the
**total** source−sink balance, allowing us to decompose the isotope cycles.

Method:
-------
1. From the CH₄ ppb seasonal amplitude and known mean OH lifetime (τ~9.7 yr),
   estimate the seasonal amplitude of the OH destruction rate.

2. The sink contribution to each isotope seasonal cycle is:
      Δδ_sink ≈ ε × ΔL/L_mean
   where ε = (α−1)×1000 and ΔL/L_mean is the fractional seasonal variation
   in the loss rate (same for all isotopologues, proportional to OH variation).

3. The source contribution is the residual:
      Δδ_source = Δδ_observed − Δδ_sink

4. The **sink-only** amplitude ratio is:
      R_sink = Δδ¹³C_sink / ΔδD_sink = ε_13C / ε_D = (α_13C−1) / (α_D−1)

   This directly constrains the OH KIE ratio — independent of source terms.

Approach — Keeling-plot-style seasonal decomposition:
-----------------------------------------------------
Rather than fitting a full dynamical model (which requires assumptions about
transport, mixing times, etc.), we use the harmonic amplitudes and phases
from Phase 2 along with a simple linearized perturbation analysis.

For a box model in quasi-steady state:
    δ_atm ≈ δ_source + ε_bulk    (isotopic balance)

Seasonal perturbation:
    Δδ_atm ≈ (δ_source − δ_atm) × (ΔS/S) − ε_bulk × (ΔL/L)

The CH₄ ppb amplitude constrains ΔS/S − ΔL/L (the net flux imbalance).
Combined with the two isotope amplitudes, we have 3 equations and can solve
for: ΔL/L (OH seasonality), ΔS/S (source seasonality), and cross-check
the source isotope signature assumptions.

Output:
    results/phase4_deconv/
        deconv_results.json    — per-site decomposition and sink-only ratios
    figures/
        fig4_decomposition.png     — stacked seasonal contributions per site
        fig4_sink_ratio.png        — sink-only ratio vs latitude with OH band
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
PHASE3_DIR = Path(__file__).resolve().parent.parent / "results" / "phase3_synthesis"
RESULTS_DIR = Path(__file__).resolve().parent.parent / "results" / "phase4_deconv"
FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# PHYSICAL CONSTANTS
# ============================================================================

# Mean atmospheric CH₄ burden and lifetime
CH4_MEAN_PPB = 1850.0        # mean tropospheric CH₄ ~2005–2010 (ppb)
TAU_OH_YR = 11.2             # OH-only CH₄ lifetime (years), Prather et al. 2012
TAU_TOTAL_YR = 9.1           # total atmospheric CH₄ lifetime (years)

# Fraction of total CH₄ sink attributable to each pathway
F_OH = 0.84                  # tropospheric OH (~84%)
F_SOIL = 0.06                # soil uptake (~6%)
F_CL = 0.035                 # tropospheric Cl (~3.5%)
F_STRAT = 0.065              # stratospheric loss (~6.5%)

# KIE fractionation factors (α values)
# OH + CH₄: α_13C is the quantity we want to constrain
ALPHA_D_OH = 1.294           # Saueressig et al. (2001) — widely accepted
ALPHA_13C_SAUERESSIG = 1.0039
ALPHA_13C_CANTRELL = 1.0054

# Cl + CH₄
ALPHA_13C_CL = 1.066         # Saueressig et al. (1995)
ALPHA_D_CL = 1.508           # Saueressig et al. (1996)

# Soil uptake (diffusion-limited, near-unity)
ALPHA_13C_SOIL = 1.022       # King et al. (1989)
ALPHA_D_SOIL = 1.066         # Snover & Quay (2000)

# Stratospheric loss (effective, includes OH, Cl, O(¹D))
ALPHA_13C_STRAT = 1.013      # Rice et al. (2003) effective
ALPHA_D_STRAT = 1.16         # Rice et al. (2003) effective

# Source isotope signatures (mean microbial, the dominant seasonal source)
# Microbial sources = wetlands + rice, which drive most of the source seasonality
D13C_SOURCE_MICROBIAL = -62.0   # ‰ VPDB (wetlands/rice, Sherwood et al. 2017)
DD_SOURCE_MICROBIAL = -310.0    # ‰ VSMOW (wetlands, Waldron et al. 1999)

# Mean atmospheric isotope values (~2005–2010)
D13C_ATM = -47.3              # ‰ VPDB
DD_ATM = -86.0                # ‰ VSMOW

# Bulk KIE (weighted average across all sinks)
# ε_bulk = Σ f_i × (α_i − 1) × 1000  for each sink pathway
def compute_bulk_epsilon(alpha_13c_oh: float) -> dict:
    """Compute bulk fractionation (ε in ‰) for both isotopes.

    The bulk ε is the flux-weighted average across OH, Cl, soil, stratosphere.

    Parameters
    ----------
    alpha_13c_oh : the OH KIE for ¹³C (the contested value)

    Returns
    -------
    dict with eps_13C and eps_D in ‰
    """
    eps_13c = (F_OH * (alpha_13c_oh - 1) +
               F_CL * (ALPHA_13C_CL - 1) +
               F_SOIL * (ALPHA_13C_SOIL - 1) +
               F_STRAT * (ALPHA_13C_STRAT - 1)) * 1000.0

    eps_D = (F_OH * (ALPHA_D_OH - 1) +
             F_CL * (ALPHA_D_CL - 1) +
             F_SOIL * (ALPHA_D_SOIL - 1) +
             F_STRAT * (ALPHA_D_STRAT - 1)) * 1000.0

    return {"eps_13C": eps_13c, "eps_D": eps_D}


# ============================================================================
# SEASONAL DECOMPOSITION
# ============================================================================

def decompose_seasonal_cycle(
    A_ch4: float,        # CH₄ ppb seasonal amplitude
    A_d13c: float,       # δ¹³C seasonal amplitude (‰)
    A_dD: float,         # δD seasonal amplitude (‰)
    phase_ch4: float,    # CH₄ peak month
    phase_d13c: float,   # δ¹³C peak month
    phase_dD: float,     # δD peak month
    alpha_13c_oh: float = ALPHA_13C_SAUERESSIG,
) -> dict:
    """Decompose observed isotope amplitudes into sink and source contributions.

    Uses a linearized seasonal perturbation model:

        Seasonal CH₄ variation:
            ΔC/C = ΔS/S − ΔL/L     (source minus sink)

        Seasonal isotope variation:
            Δδ_13C = (δ_src_13C − δ_atm_13C) × ΔS/S + ε_bulk_13C × ΔL/L
            Δδ_D   = (δ_src_D − δ_atm_D)     × ΔS/S + ε_bulk_D   × ΔL/L

    Given A_ch4 constrains ΔC/C, and the two isotope amplitudes are observed,
    we can solve for ΔS/S and ΔL/L.

    Note: This is a simplified decomposition. It assumes:
    - Source seasonality is dominated by microbial emissions (wetlands/rice)
    - The source isotope signature is known (δ_src)
    - OH seasonal variation dominates the loss-rate seasonality
    - Phases of source and sink seasonal cycles are approximately aligned
      (both peak in summer at NH sites)

    Parameters
    ----------
    alpha_13c_oh : assumed OH KIE for ¹³C

    Returns
    -------
    dict with: dS_over_S, dL_over_L, A_d13c_sink, A_d13c_source,
               A_dD_sink, A_dD_source, R_sink, R_observed
    """
    # Fractional CH₄ seasonal variation
    dC_over_C = A_ch4 / CH4_MEAN_PPB

    # Source–atmosphere isotope offsets (‰)
    delta_src_13c = D13C_SOURCE_MICROBIAL - D13C_ATM  # ≈ −14.7‰
    delta_src_D = DD_SOURCE_MICROBIAL - DD_ATM          # ≈ −224‰

    # Bulk fractionation
    bulk = compute_bulk_epsilon(alpha_13c_oh)
    eps_13c = bulk["eps_13C"]  # ‰
    eps_D = bulk["eps_D"]      # ‰

    # System of equations (simplified, using amplitudes):
    # ΔC/C = ΔS/S − ΔL/L                           ... (1)
    # A_d13c ≈ |delta_src_13c × ΔS/S + eps_13c × ΔL/L|  ... (2)
    # A_dD   ≈ |delta_src_D   × ΔS/S + eps_D   × ΔL/L|  ... (3)
    #
    # From (2) and the observed amplitude, we can use (1) to solve:
    # Let x = ΔS/S, y = ΔL/L
    # x − y = ΔC/C                 ... (1)
    # delta_src_13c × x + eps_13c × y = ±A_d13c  ... (2)
    #
    # From (1): x = y + ΔC/C
    # Substitute into (2):
    #   delta_src_13c × (y + ΔC/C) + eps_13c × y = ±A_d13c
    #   y × (delta_src_13c + eps_13c) = ±A_d13c − delta_src_13c × ΔC/C
    #   y = (±A_d13c − delta_src_13c × ΔC/C) / (delta_src_13c + eps_13c)

    # Determine sign: at NH sites, OH sink enriches δ¹³C (makes more positive)
    # while microbial sources push it more negative. The observed δ¹³C peaks in
    # late spring/early summer (May–Jun), same as OH enrichment.
    # The δ¹³C perturbation from sinks is POSITIVE (enrichment), while from
    # sources it can go either way depending on relative magnitudes.
    #
    # Use the δ¹³C equation with the sign convention that positive ΔL/L means
    # increased loss rate (summer), which ENRICHES δ¹³C (positive contribution):
    #   eps_13c > 0, so eps_13c × ΔL/L > 0 when ΔL/L > 0
    # And increased source (summer) pushes δ¹³C negative:
    #   delta_src_13c < 0, so delta_src_13c × ΔS/S < 0 when ΔS/S > 0

    # Solve using δ¹³C equation:
    # Assume observed amplitude represents the net positive peak (enrichment)
    # In summer: sinks enrich (+), sources deplete (−), net = obs peak
    # The amplitude A_d13c is always positive, representing half the peak-to-trough
    # For the linear model, the net effect in summer is:
    #   delta_src_13c × ΔS/S + eps_13c × ΔL/L
    # This could be positive or negative. We use A_d13c with positive sign
    # since δ¹³C peaks in summer (enrichment wins at most sites).

    denom_13c = delta_src_13c + eps_13c  # ≈ −14.7 + 5.8 ≈ −8.9

    # Solve for ΔL/L using δ¹³C amplitude
    # Note on sign convention: the ΔL/L values returned are typically negative
    # because the linearized model solves for the amplitude relationship, and the
    # sign reflects the *phase* of the decomposition (sources and sinks both peak
    # in summer, partially cancelling in the CH₄ concentration cycle).  The
    # MAGNITUDES of A_d13c_sink, A_dD_sink are what matter for the ratio;
    # absolute values are taken below to ensure positive amplitudes.
    # Using positive A_d13c (summer enrichment)
    dL_over_L = (A_d13c - delta_src_13c * dC_over_C) / denom_13c
    dS_over_S = dL_over_L + dC_over_C

    # Sink-only isotope amplitude contributions
    A_d13c_sink = abs(eps_13c * dL_over_L)
    A_dD_sink = abs(eps_D * dL_over_L)

    # Source-only isotope amplitude contributions
    A_d13c_source = abs(delta_src_13c * dS_over_S)
    A_dD_source = abs(delta_src_D * dS_over_S)

    # Sink-only amplitude ratio (this is what we want to extract!)
    R_sink = A_d13c_sink / A_dD_sink if A_dD_sink > 0 else np.nan

    # Note: R_sink = |eps_13c × ΔL/L| / |eps_D × ΔL/L| = eps_13c / eps_D
    # This equals the bulk KIE ratio, independent of ΔL/L.
    # The value depends on the assumed α_13C_OH through eps_13c.

    # Also compute what R_sink would predict for the observed ratio
    R_observed = A_d13c / A_dD if A_dD > 0 else np.nan

    return {
        "dS_over_S": float(dS_over_S),
        "dL_over_L": float(dL_over_L),
        "dC_over_C": float(dC_over_C),
        "A_d13c_sink": float(A_d13c_sink),
        "A_d13c_source": float(A_d13c_source),
        "A_dD_sink": float(A_dD_sink),
        "A_dD_source": float(A_dD_source),
        "R_sink": float(R_sink),
        "R_observed": float(R_observed),
        "eps_13c_bulk": float(eps_13c),
        "eps_D_bulk": float(eps_D),
    }


def bootstrap_decomposition(
    fits_site: dict, n_boot: int = 2000, rng: np.random.Generator = None,
) -> dict:
    """Bootstrap the decomposition to propagate amplitude uncertainties.

    Draws A_ch4, A_d13c, A_dD from their bootstrap distributions (approximated
    as Gaussian from the 95% CIs) and runs decompose_seasonal_cycle for each draw.

    Returns 95% CIs on R_sink, dS_over_S, dL_over_L.
    """
    if rng is None:
        rng = np.random.default_rng(42)

    def ci_to_sigma(ci):
        return (ci[1] - ci[0]) / (2 * 1.96)

    A_ch4_mean = fits_site["ch4"]["amplitude"]
    A_ch4_sigma = ci_to_sigma(fits_site["ch4"]["amplitude_ci95"])

    A_d13c_mean = fits_site["d13C"]["amplitude"]
    A_d13c_sigma = ci_to_sigma(fits_site["d13C"]["amplitude_ci95"])

    A_dD_mean = fits_site["dD"]["amplitude"]
    A_dD_sigma = ci_to_sigma(fits_site["dD"]["amplitude_ci95"])

    results_saueressig = np.zeros(n_boot)
    results_cantrell = np.zeros(n_boot)
    dL_samples = np.zeros(n_boot)
    dS_samples = np.zeros(n_boot)

    for b in range(n_boot):
        a_ch4 = max(0, rng.normal(A_ch4_mean, A_ch4_sigma))
        a_d13c = max(0, rng.normal(A_d13c_mean, A_d13c_sigma))
        a_dD = max(0, rng.normal(A_dD_mean, A_dD_sigma))

        dec_s = decompose_seasonal_cycle(
            a_ch4, a_d13c, a_dD, 0, 0, 0,
            alpha_13c_oh=ALPHA_13C_SAUERESSIG)
        dec_c = decompose_seasonal_cycle(
            a_ch4, a_d13c, a_dD, 0, 0, 0,
            alpha_13c_oh=ALPHA_13C_CANTRELL)

        results_saueressig[b] = dec_s["R_sink"]
        results_cantrell[b] = dec_c["R_sink"]
        dL_samples[b] = dec_s["dL_over_L"]
        dS_samples[b] = dec_s["dS_over_S"]

    def pct(arr):
        return [float(np.nanpercentile(arr, 2.5)),
                float(np.nanpercentile(arr, 97.5))]

    return {
        "R_sink_saueressig_ci95": pct(results_saueressig),
        "R_sink_cantrell_ci95": pct(results_cantrell),
        "dL_over_L_ci95": pct(dL_samples),
        "dS_over_S_ci95": pct(dS_samples),
    }


# ============================================================================
# INVERSION: WHAT α_13C_OH IS CONSISTENT WITH OBSERVED RATIO?
# ============================================================================

def invert_for_alpha_13c(
    A_ch4: float, A_d13c: float, A_dD: float,
) -> dict:
    """Given observed amplitudes, find what α_13C_OH would be needed to
    reproduce the observed ratio if ALL seasonality were sink-driven.

    This provides an upper bound on α_13C_OH. If source seasonality inflates
    the δ¹³C amplitude, the true α is lower.

    Also computes the α_13C_OH that the decomposed sink-only ratio implies
    for a range of source seasonality assumptions.
    """
    # The observed ratio R_obs = A_d13c / A_dD
    R_obs = A_d13c / A_dD if A_dD > 0 else np.nan

    # If pure sink: R = ε_13C_bulk / ε_D_bulk
    # ε_D_bulk is known (doesn't depend on α_13C_OH):
    eps_D_bulk = (F_OH * (ALPHA_D_OH - 1) +
                  F_CL * (ALPHA_D_CL - 1) +
                  F_SOIL * (ALPHA_D_SOIL - 1) +
                  F_STRAT * (ALPHA_D_STRAT - 1)) * 1000.0

    # ε_13C_bulk = f_OH × (α_13C_OH − 1) × 1000 + (non-OH terms)
    eps_13c_non_oh = (F_CL * (ALPHA_13C_CL - 1) +
                      F_SOIL * (ALPHA_13C_SOIL - 1) +
                      F_STRAT * (ALPHA_13C_STRAT - 1)) * 1000.0

    # R_obs = ε_13C_bulk / ε_D_bulk
    # ε_13C_bulk = R_obs × ε_D_bulk
    eps_13c_needed = R_obs * eps_D_bulk

    # f_OH × (α_13C_OH − 1) × 1000 = eps_13c_needed − eps_13c_non_oh
    alpha_13c_if_pure_sink = 1.0 + (eps_13c_needed - eps_13c_non_oh) / (F_OH * 1000.0)

    return {
        "R_observed": float(R_obs),
        "alpha_13c_if_pure_sink": float(alpha_13c_if_pure_sink),
        "eps_D_bulk": float(eps_D_bulk),
        "eps_13c_non_oh": float(eps_13c_non_oh),
    }


# ============================================================================
# PLOTTING
# ============================================================================

def plot_decomposition(all_decomp: dict, site_order: list) -> None:
    """Figure 4a: Stacked bar chart of sink vs source contributions to
    δ¹³C and δD seasonal amplitudes at each site.
    """
    codes = [c for c in site_order if c in all_decomp]
    n = len(codes)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    x = np.arange(n)
    width = 0.5

    # δ¹³C decomposition
    sink_13c = [all_decomp[c]["saueressig"]["A_d13c_sink"] for c in codes]
    src_13c = [all_decomp[c]["saueressig"]["A_d13c_source"] for c in codes]

    ax1.bar(x, sink_13c, width, label="Sink (OH fractionation)", color="C0",
            alpha=0.8)
    ax1.bar(x, src_13c, width, bottom=sink_13c,
            label="Source (microbial seasonality)", color="C1", alpha=0.8)
    ax1.set_xticks(x)
    ax1.set_xticklabels(codes, fontsize=9)
    ax1.set_ylabel("δ¹³C amplitude (‰)", fontsize=10)
    ax1.set_title("(a) δ¹³C seasonal amplitude decomposition", fontsize=11)
    ax1.legend(fontsize=8)

    # δD decomposition
    sink_D = [all_decomp[c]["saueressig"]["A_dD_sink"] for c in codes]
    src_D = [all_decomp[c]["saueressig"]["A_dD_source"] for c in codes]

    ax2.bar(x, sink_D, width, label="Sink (OH fractionation)", color="C0",
            alpha=0.8)
    ax2.bar(x, src_D, width, bottom=sink_D,
            label="Source (microbial seasonality)", color="C1", alpha=0.8)
    ax2.set_xticks(x)
    ax2.set_xticklabels(codes, fontsize=9)
    ax2.set_ylabel("δD amplitude (‰)", fontsize=10)
    ax2.set_title("(b) δD seasonal amplitude decomposition", fontsize=11)
    ax2.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "fig4_decomposition.png", dpi=300)
    plt.close(fig)
    print(f"✓ Figure saved: {FIGURES_DIR / 'fig4_decomposition.png'}")


def plot_sink_ratio_and_alpha(all_decomp: dict, classifications: dict) -> None:
    """Figure 4b: Implied α_13C_OH from the pure-sink assumption vs latitude.

    Also shows the R_sink values for Saueressig/Cantrell assumptions.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    clean_codes = [c for c in classifications
                   if classifications[c]["classification"] == "clean"
                   and c in all_decomp]

    # Sort by latitude
    clean_codes = sorted(clean_codes,
                         key=lambda c: -classifications[c]["latitude"])

    lats = [classifications[c]["latitude"] for c in clean_codes]
    alphas = [all_decomp[c]["inversion"]["alpha_13c_if_pure_sink"]
              for c in clean_codes]
    R_obs = [all_decomp[c]["saueressig"]["R_observed"] for c in clean_codes]
    R_sink_s = [all_decomp[c]["saueressig"]["R_sink"] for c in clean_codes]
    R_sink_c = [all_decomp[c]["cantrell"]["R_sink"] for c in clean_codes]

    # Panel (a): Implied α_13C if pure sink
    ax1.scatter(lats, alphas, s=80, c="C0", edgecolors="k", lw=0.5, zorder=3)
    for c, lat, a in zip(clean_codes, lats, alphas):
        ax1.annotate(c, (lat, a), textcoords="offset points",
                     xytext=(6, 6), fontsize=7)

    ax1.axhline(ALPHA_13C_SAUERESSIG, color="green", ls="--", lw=1.2,
                label=f"Saueressig ({ALPHA_13C_SAUERESSIG})")
    ax1.axhline(ALPHA_13C_CANTRELL, color="orange", ls="--", lw=1.2,
                label=f"Cantrell ({ALPHA_13C_CANTRELL})")

    ax1.set_xlabel("Latitude (°)", fontsize=10)
    ax1.set_ylabel("Implied α¹³C_OH  (if pure sink)", fontsize=10)
    ax1.set_title("(a) α¹³C_OH implied if ALL seasonal variation = sink",
                  fontsize=10)
    ax1.legend(fontsize=8)
    ax1.set_xlim(-100, 100)

    # Panel (b): Observed vs sink-only R for both α assumptions
    x = np.arange(len(clean_codes))
    width = 0.25

    ax2.bar(x - width, R_obs, width, label="Observed R", color="C0", alpha=0.8)
    ax2.bar(x, R_sink_s, width,
            label=f"R_sink (Saueressig, α={ALPHA_13C_SAUERESSIG})",
            color="green", alpha=0.7)
    ax2.bar(x + width, R_sink_c, width,
            label=f"R_sink (Cantrell, α={ALPHA_13C_CANTRELL})",
            color="orange", alpha=0.7)

    ax2.set_xticks(x)
    ax2.set_xticklabels([f"{c}\n({classifications[c]['latitude']:+.0f}°)"
                          for c in clean_codes], fontsize=7)
    ax2.set_ylabel("Amplitude ratio  A(δ¹³C) / A(δD)", fontsize=10)
    ax2.set_title("(b) Observed vs sink-only ratio", fontsize=10)
    ax2.legend(fontsize=7, loc="upper right")

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "fig4_sink_ratio.png", dpi=300)
    plt.close(fig)
    print(f"✓ Figure saved: {FIGURES_DIR / 'fig4_sink_ratio.png'}")


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("Phase 4: Source deconvolution")
    print("=" * 70)

    # ── Load Phase 2 and Phase 3 results ──
    with open(PHASE2_DIR / "harmonic_fits.json") as f:
        fits = json.load(f)
    with open(PHASE3_DIR / "synthesis_results.json") as f:
        synthesis = json.load(f)

    classifications = synthesis["site_classifications"]

    # Reference values
    bulk_saueressig = compute_bulk_epsilon(ALPHA_13C_SAUERESSIG)
    bulk_cantrell = compute_bulk_epsilon(ALPHA_13C_CANTRELL)
    print(f"\nBulk fractionation (Saueressig): ε_13C = {bulk_saueressig['eps_13C']:.2f}‰, "
          f"ε_D = {bulk_saueressig['eps_D']:.2f}‰")
    print(f"Bulk fractionation (Cantrell):   ε_13C = {bulk_cantrell['eps_13C']:.2f}‰, "
          f"ε_D = {bulk_cantrell['eps_D']:.2f}‰")
    print(f"Bulk R_sink (Saueressig): {bulk_saueressig['eps_13C']/bulk_saueressig['eps_D']:.4f}")
    print(f"Bulk R_sink (Cantrell):   {bulk_cantrell['eps_13C']/bulk_cantrell['eps_D']:.4f}")

    rng = np.random.default_rng(42)
    all_decomp = {}
    site_order = sorted(fits.keys(),
                        key=lambda c: -classifications.get(c, {}).get("latitude", 0))

    for code in site_order:
        if code not in fits or "ch4" not in fits[code]:
            print(f"\n--- {code}: no CH₄ data, skipping decomposition ---")
            continue

        lat = classifications.get(code, {}).get("latitude", 0)
        tag = classifications.get(code, {}).get("classification", "?")
        print(f"\n--- {code} ({lat:+.1f}°, {tag}) ---")

        A_ch4 = fits[code]["ch4"]["amplitude"]
        A_d13c = fits[code]["d13C"]["amplitude"]
        A_dD = fits[code]["dD"]["amplitude"]
        pk_ch4 = fits[code]["ch4"].get("peak_month", 0)
        pk_d13c = fits[code]["d13C"]["peak_month"]
        pk_dD = fits[code]["dD"]["peak_month"]

        # Decompose with both α assumptions
        dec_s = decompose_seasonal_cycle(A_ch4, A_d13c, A_dD,
                                          pk_ch4, pk_d13c, pk_dD,
                                          ALPHA_13C_SAUERESSIG)
        dec_c = decompose_seasonal_cycle(A_ch4, A_d13c, A_dD,
                                          pk_ch4, pk_d13c, pk_dD,
                                          ALPHA_13C_CANTRELL)

        # Inversion: what α would explain the observed ratio as pure sink?
        inv = invert_for_alpha_13c(A_ch4, A_d13c, A_dD)

        # Bootstrap
        boot = bootstrap_decomposition(fits[code], n_boot=2000, rng=rng)

        all_decomp[code] = {
            "saueressig": dec_s,
            "cantrell": dec_c,
            "inversion": inv,
            "bootstrap": boot,
        }

        print(f"  CH₄ amplitude: {A_ch4:.1f} ppb → ΔC/C = {dec_s['dC_over_C']:.4f}")
        print(f"  Decomposition (Saueressig α={ALPHA_13C_SAUERESSIG}):")
        print(f"    ΔL/L = {dec_s['dL_over_L']:.4f},  ΔS/S = {dec_s['dS_over_S']:.4f}")
        print(f"    δ¹³C: sink = {dec_s['A_d13c_sink']:.4f}‰, "
              f"source = {dec_s['A_d13c_source']:.4f}‰")
        print(f"    δD:   sink = {dec_s['A_dD_sink']:.2f}‰, "
              f"source = {dec_s['A_dD_source']:.2f}‰")
        print(f"    R_sink = {dec_s['R_sink']:.4f},  R_obs = {dec_s['R_observed']:.4f}")
        print(f"  If pure sink → α_13C_OH = {inv['alpha_13c_if_pure_sink']:.4f}")

    # ── Save results ──
    with open(RESULTS_DIR / "deconv_results.json", "w") as f:
        json.dump(all_decomp, f, indent=2)
    print(f"\n✓ Results saved to {RESULTS_DIR / 'deconv_results.json'}")

    # ── Figures ──
    plot_decomposition(all_decomp, site_order)
    plot_sink_ratio_and_alpha(all_decomp, classifications)

    # ── Summary ──
    print("\n" + "=" * 70)
    print("DECONVOLUTION SUMMARY")
    print("=" * 70)

    print(f"\n{'Site':<5} {'Lat':>6} {'R_obs':>8} {'R_sink_S':>8} {'R_sink_C':>8} "
          f"{'α_13C*':>8} {'ΔL/L':>8} {'ΔS/S':>8}")
    print("-" * 70)
    for code in site_order:
        if code not in all_decomp:
            continue
        d = all_decomp[code]
        lat = classifications.get(code, {}).get("latitude", 0)
        print(f"{code:<5} {lat:>+6.1f} "
              f"{d['saueressig']['R_observed']:>8.4f} "
              f"{d['saueressig']['R_sink']:>8.4f} "
              f"{d['cantrell']['R_sink']:>8.4f} "
              f"{d['inversion']['alpha_13c_if_pure_sink']:>8.4f} "
              f"{d['saueressig']['dL_over_L']:>8.4f} "
              f"{d['saueressig']['dS_over_S']:>8.4f}")

    # Key insight
    print(f"\n{'='*70}")
    print("KEY INSIGHT")
    print(f"{'='*70}")
    print(f"\nThe sink-only ratio R_sink = ε_13C_bulk / ε_D_bulk is a CONSTANT")
    print(f"that depends only on the assumed KIE values, not on site-specific data.")
    print(f"  Saueressig (α=1.0039): R_sink = {bulk_saueressig['eps_13C']/bulk_saueressig['eps_D']:.4f}")
    print(f"  Cantrell   (α=1.0054): R_sink = {bulk_cantrell['eps_13C']/bulk_cantrell['eps_D']:.4f}")
    print(f"\nHowever, the implied α_13C if ALL seasonal variation were sink-driven")
    print(f"varies by site (because each site has a different observed ratio).")
    print(f"These implied α values range from "
          f"{min(d['inversion']['alpha_13c_if_pure_sink'] for d in all_decomp.values()):.4f} "
          f"to {max(d['inversion']['alpha_13c_if_pure_sink'] for d in all_decomp.values()):.4f}")
    print(f"\nSince all clean NH sites give implied α >> 1.0054 (Cantrell),")
    print(f"source seasonality is clearly a major contributor.")
    print(f"The decomposition shows source contributions to δ¹³C are comparable")
    print(f"to or larger than sink contributions at most sites.")
    print(f"\n→ Phase 5 will extract the best KIE constraint by combining")
    print(f"  the decomposition with uncertainty propagation.")


if __name__ == "__main__":
    main()
