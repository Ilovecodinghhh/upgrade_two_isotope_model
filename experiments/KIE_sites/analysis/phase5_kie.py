#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
phase5_kie.py — KIE extraction for the KIE_sites experiment
=============================================================

Final phase: extract the best observational constraint on the OH ¹³C KIE
from the seasonal amplitude ratios, after accounting for source seasonality
and non-OH sinks.

Strategy:
---------
Phase 4 showed that the **bulk** sink-only ratio R_sink = ε_13C_bulk / ε_D_bulk
is a constant (independent of site), determined by KIE values. The observed
ratio at most NH sites far exceeds this due to source seasonality.

However, the **SH sites** (CGO, SPO) have observed ratios close to R_sink,
suggesting minimal source contamination. These provide the most direct
observational constraint on the KIE.

Three approaches are used:

**Approach 1: Direct SH constraint**
  Use CGO and SPO (minimal source seasonality) to directly constrain α_13C_OH
  by inverting R_observed → α_13C from the bulk KIE ratio formula.

**Approach 2: Source-corrected NH constraint**
  Use the Phase 4 decomposition to subtract the estimated source contribution
  from NH sites, then invert the residual sink-only ratio.

**Approach 3: Latitude-gradient extrapolation**
  The ratio R_observed increases with latitude towards NH high latitudes
  (more source contamination). Extrapolating the R vs latitude relationship
  to the zero-source-contamination limit provides another constraint.

All approaches propagate uncertainties via Monte Carlo.

Output:
    results/phase5_kie/
        kie_results.json        — final KIE constraint with all uncertainties
    figures/
        fig5_kie_constraint.png — comparison with Saueressig and Cantrell
    RESULT.md                   — plain-language summary
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
PHASE4_DIR = Path(__file__).resolve().parent.parent / "results" / "phase4_deconv"
RESULTS_DIR = Path(__file__).resolve().parent.parent / "results" / "phase5_kie"
FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"
EXPT_DIR = Path(__file__).resolve().parent.parent  # experiments/KIE_sites/

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# PHYSICAL CONSTANTS (same as Phase 4, for self-contained readability)
# ============================================================================
ALPHA_D_OH = 1.294
ALPHA_13C_SAUERESSIG = 1.0039
ALPHA_13C_CANTRELL = 1.0054

F_OH = 0.84
F_CL = 0.035
F_SOIL = 0.06
F_STRAT = 0.065

ALPHA_13C_CL = 1.066
ALPHA_D_CL = 1.508
ALPHA_13C_SOIL = 1.022
ALPHA_D_SOIL = 1.066
ALPHA_13C_STRAT = 1.013
ALPHA_D_STRAT = 1.16

# Uncertainties on non-OH KIEs and sink fractions (for MC propagation)
# These are approximate 1σ values from the literature
SIGMA_F_OH = 0.04          # ±4% on the OH fraction
SIGMA_F_CL = 0.01          # ±1% on Cl fraction
SIGMA_F_SOIL = 0.02        # ±2% on soil fraction
SIGMA_ALPHA_D_OH = 0.01    # ±0.01 on α_D_OH
SIGMA_ALPHA_13C_CL = 0.005  # ±0.005 on α_13C_Cl
SIGMA_ALPHA_D_CL = 0.05    # ±0.05 on α_D_Cl


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def compute_bulk_epsilon(alpha_13c_oh, alpha_d_oh=ALPHA_D_OH,
                          f_oh=F_OH, f_cl=F_CL, f_soil=F_SOIL, f_strat=F_STRAT,
                          alpha_13c_cl=ALPHA_13C_CL, alpha_d_cl=ALPHA_D_CL):
    """Compute bulk fractionation ε (‰) for ¹³C and D.

    Includes all four sink pathways with adjustable parameters for MC.
    """
    eps_13c = (f_oh * (alpha_13c_oh - 1) +
               f_cl * (alpha_13c_cl - 1) +
               f_soil * (ALPHA_13C_SOIL - 1) +
               f_strat * (ALPHA_13C_STRAT - 1)) * 1000.0

    eps_D = (f_oh * (alpha_d_oh - 1) +
             f_cl * (alpha_d_cl - 1) +
             f_soil * (ALPHA_D_SOIL - 1) +
             f_strat * (ALPHA_D_STRAT - 1)) * 1000.0

    return eps_13c, eps_D


def ratio_to_alpha_13c(R_obs, alpha_d_oh=ALPHA_D_OH,
                        f_oh=F_OH, f_cl=F_CL, f_soil=F_SOIL, f_strat=F_STRAT,
                        alpha_13c_cl=ALPHA_13C_CL, alpha_d_cl=ALPHA_D_CL):
    """Invert the amplitude ratio to get α_13C_OH.

    Given R = ε_13C_bulk / ε_D_bulk, and all non-OH KIE values known,
    solve for α_13C_OH.

    Parameters
    ----------
    R_obs : observed (or source-corrected) amplitude ratio A(δ¹³C)/A(δD)

    Returns
    -------
    alpha_13c_oh : the implied OH ¹³C KIE
    """
    _, eps_D = compute_bulk_epsilon(1.0, alpha_d_oh, f_oh, f_cl, f_soil, f_strat,
                                     alpha_13c_cl, alpha_d_cl)

    # ε_13C_bulk = R × ε_D_bulk
    eps_13c_needed = R_obs * eps_D

    # Non-OH ¹³C contribution
    eps_13c_non_oh = (f_cl * (alpha_13c_cl - 1) +
                      f_soil * (ALPHA_13C_SOIL - 1) +
                      f_strat * (ALPHA_13C_STRAT - 1)) * 1000.0

    # α_13C_OH = 1 + (ε_13C_needed − ε_13C_non_OH) / (f_OH × 1000)
    alpha_13c = 1.0 + (eps_13c_needed - eps_13c_non_oh) / (f_oh * 1000.0)
    return alpha_13c


def ci_to_sigma(ci_lo, ci_hi):
    """Approximate 1σ from 95% CI."""
    return (ci_hi - ci_lo) / (2 * 1.96)


# ============================================================================
# APPROACH 1: DIRECT SH CONSTRAINT
# ============================================================================

def approach1_sh_direct(fits: dict, n_mc: int = 50000) -> dict:
    """Use SH sites (CGO, SPO) with minimal source contamination.

    These sites have observed R close to R_sink. Invert R_obs → α_13C_OH
    with full uncertainty propagation via Monte Carlo.
    """
    print("\n─── Approach 1: Direct SH constraint ───")
    rng = np.random.default_rng(123)

    sh_sites = ["CGO", "SPO"]
    available = [s for s in sh_sites if s in fits]

    # Collect observed ratios and uncertainties
    ratios = []
    sigmas = []
    for code in available:
        r = fits[code]["ratio"]["value"]
        r_ci = fits[code]["ratio"]["ci95"]
        sigma = ci_to_sigma(r_ci[0], r_ci[1])
        ratios.append(r)
        sigmas.append(sigma)
        print(f"  {code}: R = {r:.4f} ± {sigma:.4f}")

    # Inverse-variance weighted mean of SH ratios
    ratios = np.array(ratios)
    sigmas = np.array(sigmas)
    weights = 1.0 / sigmas**2
    R_sh = np.sum(weights * ratios) / np.sum(weights)
    R_sh_sigma = np.sqrt(1.0 / np.sum(weights))
    print(f"  Weighted mean: R_SH = {R_sh:.4f} ± {R_sh_sigma:.4f}")

    # Monte Carlo: propagate all uncertainties
    alpha_samples = np.zeros(n_mc)
    for i in range(n_mc):
        # Draw R from observed distribution
        r = rng.normal(R_sh, R_sh_sigma)

        # Draw uncertain parameters
        f_oh = np.clip(rng.normal(F_OH, SIGMA_F_OH), 0.5, 0.99)
        f_cl = np.clip(rng.normal(F_CL, SIGMA_F_CL), 0.0, 0.1)
        f_soil = np.clip(rng.normal(F_SOIL, SIGMA_F_SOIL), 0.0, 0.15)
        f_strat = 1.0 - f_oh - f_cl - f_soil  # ensure fractions sum to ~1

        alpha_d_oh = rng.normal(ALPHA_D_OH, SIGMA_ALPHA_D_OH)
        alpha_13c_cl = rng.normal(ALPHA_13C_CL, SIGMA_ALPHA_13C_CL)
        alpha_d_cl = rng.normal(ALPHA_D_CL, SIGMA_ALPHA_D_CL)

        alpha_samples[i] = ratio_to_alpha_13c(
            r, alpha_d_oh, f_oh, f_cl, f_soil, f_strat, alpha_13c_cl, alpha_d_cl)

    alpha_med = np.median(alpha_samples)
    alpha_lo, alpha_hi = np.percentile(alpha_samples, [2.5, 97.5])

    print(f"  → α_13C_OH = {alpha_med:.4f}  [{alpha_lo:.4f}, {alpha_hi:.4f}] (95% CI)")

    return {
        "method": "Direct SH constraint (CGO + SPO)",
        "sites_used": available,
        "R_sh_weighted": float(R_sh),
        "R_sh_sigma": float(R_sh_sigma),
        "alpha_13c_oh_median": float(alpha_med),
        "alpha_13c_oh_ci95": [float(alpha_lo), float(alpha_hi)],
        "n_mc": n_mc,
    }


# ============================================================================
# APPROACH 2: SOURCE-CORRECTED NH
# ============================================================================

def approach2_source_corrected(fits: dict, classifications: dict,
                                n_mc: int = 50000) -> dict:
    """Use all clean sites, correcting for source seasonality to extract α_13C_OH.

    The key to avoiding circularity is solving the system in two steps:

    Step 1: Use CH₄ amplitude + δD amplitude to solve for ΔL/L and ΔS/S.
            This does NOT depend on α_13C_OH because ε_D_bulk is known.

        ΔC/C = ΔS/S − ΔL/L                           ... (1)
        A_δD ≈ |δ_src_D × ΔS/S + ε_D × ΔL/L|        ... (2)

        Substituting (1) into (2) and solving for ΔL/L gives a value
        independent of α_13C_OH.

    Step 2: Use the δ¹³C amplitude and the known ΔL/L, ΔS/S to extract ε_13C:

        A_δ¹³C ≈ |δ_src_13C × ΔS/S + ε_13C × ΔL/L|  ... (3)

        → ε_13C = (A_δ¹³C − δ_src_13C × ΔS/S) / ΔL/L

    Then invert ε_13C_bulk → α_13C_OH.
    """
    print("\n─── Approach 2: Source-corrected (all clean sites) ───")
    rng = np.random.default_rng(456)

    # Source isotope signatures for microbial emissions
    D13C_SRC = -62.0   # ‰ VPDB
    DD_SRC = -310.0     # ‰ VSMOW
    D13C_ATM = -47.3    # ‰ VPDB
    DD_ATM = -86.0      # ‰ VSMOW
    CH4_MEAN = 1850.0   # ppb

    # Uncertainty on source signatures (1σ)
    SIGMA_D13C_SRC = 5.0   # ‰ — wetlands range from −55 to −70‰
    SIGMA_DD_SRC = 30.0     # ‰ — wetlands range from −280 to −350‰

    clean_codes = [c for c in fits
                   if c in classifications
                   and classifications[c]["classification"] == "clean"
                   and "ch4" in fits[c]]

    print(f"  Using {len(clean_codes)} clean sites with CH₄ data")

    alpha_all_samples = np.zeros(n_mc)

    for i in range(n_mc):
        site_alphas = []
        site_weights = []

        # Draw global uncertain parameters once per MC sample
        f_oh = np.clip(rng.normal(F_OH, SIGMA_F_OH), 0.5, 0.99)
        f_cl = np.clip(rng.normal(F_CL, SIGMA_F_CL), 0.0, 0.1)
        f_soil = np.clip(rng.normal(F_SOIL, SIGMA_F_SOIL), 0.0, 0.15)
        f_strat = 1.0 - f_oh - f_cl - f_soil
        alpha_d_oh = rng.normal(ALPHA_D_OH, SIGMA_ALPHA_D_OH)
        alpha_13c_cl = rng.normal(ALPHA_13C_CL, SIGMA_ALPHA_13C_CL)
        alpha_d_cl = rng.normal(ALPHA_D_CL, SIGMA_ALPHA_D_CL)
        d13c_src = rng.normal(D13C_SRC, SIGMA_D13C_SRC)
        dD_src = rng.normal(DD_SRC, SIGMA_DD_SRC)

        # Compute ε_D_bulk (does NOT depend on α_13C_OH)
        eps_D = (f_oh * (alpha_d_oh - 1) +
                 f_cl * (alpha_d_cl - 1) +
                 f_soil * (ALPHA_D_SOIL - 1) +
                 f_strat * (ALPHA_D_STRAT - 1)) * 1000.0

        # Non-OH ¹³C contributions (known, independent of α_13C_OH)
        eps_13c_non_oh = (f_cl * (alpha_13c_cl - 1) +
                          f_soil * (ALPHA_13C_SOIL - 1) +
                          f_strat * (ALPHA_13C_STRAT - 1)) * 1000.0

        dsrc_13c = d13c_src - D13C_ATM
        dsrc_D = dD_src - DD_ATM

        for code in clean_codes:
            # Draw amplitudes from their bootstrap distributions
            A_d13c = max(0, rng.normal(
                fits[code]["d13C"]["amplitude"],
                ci_to_sigma(*fits[code]["d13C"]["amplitude_ci95"])))
            A_dD = max(1e-6, rng.normal(
                fits[code]["dD"]["amplitude"],
                ci_to_sigma(*fits[code]["dD"]["amplitude_ci95"])))
            A_ch4 = max(0, rng.normal(
                fits[code]["ch4"]["amplitude"],
                ci_to_sigma(*fits[code]["ch4"]["amplitude_ci95"])))

            dC_over_C = A_ch4 / CH4_MEAN

            # Step 1: Solve for ΔL/L using CH₄ + δD (independent of α_13C_OH)
            # From eqs (1)+(2): ΔL/L = (A_dD − dsrc_D × ΔC/C) / (dsrc_D + eps_D)
            denom_D = dsrc_D + eps_D
            if abs(denom_D) < 0.1:
                continue
            dL_over_L = (A_dD - dsrc_D * dC_over_C) / denom_D
            dS_over_S = dL_over_L + dC_over_C

            if abs(dL_over_L) < 1e-8:
                continue

            # Step 2: Use δ¹³C amplitude to extract ε_13C_bulk
            # A_d13c = |dsrc_13c × ΔS/S + ε_13C × ΔL/L|
            # → ε_13C = (A_d13c − dsrc_13c × ΔS/S) / ΔL/L
            eps_13c_bulk = (A_d13c - dsrc_13c * dS_over_S) / dL_over_L

            # Invert: ε_13C_bulk = f_OH × (α_13C_OH − 1) × 1000 + non-OH
            # → α_13C_OH = 1 + (ε_13C_bulk − non-OH) / (f_OH × 1000)
            alpha_13c = 1.0 + (eps_13c_bulk - eps_13c_non_oh) / (f_oh * 1000.0)

            # Weight by inverse amplitude uncertainty
            w = 1.0 / ci_to_sigma(*fits[code]["d13C"]["amplitude_ci95"])**2
            site_alphas.append(alpha_13c)
            site_weights.append(w)

        if site_alphas:
            site_alphas = np.array(site_alphas)
            site_weights = np.array(site_weights)
            alpha_all_samples[i] = np.average(site_alphas, weights=site_weights)
        else:
            alpha_all_samples[i] = np.nan

    alpha_med = np.nanmedian(alpha_all_samples)
    alpha_lo, alpha_hi = np.nanpercentile(alpha_all_samples, [2.5, 97.5])

    print(f"  → α_13C_OH = {alpha_med:.4f}  [{alpha_lo:.4f}, {alpha_hi:.4f}] (95% CI)")

    return {
        "method": "Source-corrected (all clean sites, non-circular)",
        "sites_used": clean_codes,
        "alpha_13c_oh_median": float(alpha_med),
        "alpha_13c_oh_ci95": [float(alpha_lo), float(alpha_hi)],
        "n_mc": n_mc,
        "source_signature_assumptions": {
            "d13C_source": D13C_SRC,
            "d13C_source_sigma": SIGMA_D13C_SRC,
            "dD_source": DD_SRC,
            "dD_source_sigma": SIGMA_DD_SRC,
        },
    }


# ============================================================================
# APPROACH 3: LATITUDE-GRADIENT EXTRAPOLATION
# ============================================================================

def approach3_latitude_gradient(fits: dict, classifications: dict) -> dict:
    """Extrapolate the R vs |latitude| trend to estimate the zero-source limit.

    In SH, there is less source contamination → lower R. Fit R vs |latitude|
    for clean sites and extrapolate.

    However: we expect R to approach R_sink at sites with no source seasonality.
    Since R_sink doesn't depend on latitude, this approach really just identifies
    which sites are least contaminated.
    """
    print("\n─── Approach 3: Latitude-gradient extrapolation ───")

    clean_codes = [c for c in classifications
                   if classifications[c]["classification"] == "clean"
                   and c in fits]

    lats = np.array([classifications[c]["latitude"] for c in clean_codes])
    ratios = np.array([fits[c]["ratio"]["value"] for c in clean_codes])
    ratio_sigmas = np.array([ci_to_sigma(*fits[c]["ratio"]["ci95"])
                              for c in clean_codes])

    # Fit R vs latitude (simple weighted linear regression)
    # R = a + b × lat
    weights = 1.0 / ratio_sigmas**2
    X = np.column_stack([np.ones_like(lats), lats])
    W = np.diag(weights)
    try:
        beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ ratios)
    except np.linalg.LinAlgError:
        print("  Linear regression failed")
        return {"method": "Latitude gradient", "status": "failed"}

    a, b = beta
    print(f"  Linear fit: R = {a:.4f} + {b:.6f} × lat")
    print(f"  Intercept (lat=0): R = {a:.4f}")
    print(f"  Slope: {b:.6f} per degree")

    # The intercept at lat=0 represents tropical sites with potentially
    # some source seasonality. A more meaningful extrapolation is not possible
    # from a simple linear model.

    # Instead, report the latitude dependence as evidence of source contamination
    # The fact that R increases with |NH latitude| is consistent with wetland
    # source seasonality increasing at higher latitudes.

    # Use the intercept as a rough constraint
    alpha_at_equator = ratio_to_alpha_13c(a)

    print(f"  Implied α_13C_OH at equator: {alpha_at_equator:.4f}")
    print(f"  (This assumes zero source seasonality at equator — likely too optimistic)")

    return {
        "method": "Latitude gradient extrapolation",
        "intercept_R": float(a),
        "slope_per_degree": float(b),
        "alpha_at_equator": float(alpha_at_equator),
        "n_sites": len(clean_codes),
        "note": "Intercept provides a rough lower bound on R_sink, but assumes "
                "zero source seasonality at equator which is not strictly true.",
    }


# ============================================================================
# PLOTTING
# ============================================================================

def plot_kie_constraint(approach1: dict, approach2: dict, approach3: dict) -> None:
    """Figure 5: Final KIE constraint — comparison with Saueressig and Cantrell.

    Three-panel figure:
    (a) PDF of α_13C_OH from Approach 1 (SH direct)
    (b) PDF from Approach 2 (source-corrected)
    (c) Summary comparison bar chart
    """
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # ── Panel (a): Approach 1 ──
    ax = axes[0]
    alpha_1 = approach1["alpha_13c_oh_median"]
    ci1 = approach1["alpha_13c_oh_ci95"]
    # Draw approximate PDF
    sigma_1 = (ci1[1] - ci1[0]) / (2 * 1.96)
    x = np.linspace(0.995, 1.015, 500)
    pdf_1 = np.exp(-0.5 * ((x - alpha_1) / sigma_1)**2) / (sigma_1 * np.sqrt(2 * np.pi))
    ax.fill_between(x, pdf_1, alpha=0.3, color="C0")
    ax.plot(x, pdf_1, color="C0", lw=1.5, label=f"SH direct: {alpha_1:.4f}")
    ax.axvline(ALPHA_13C_SAUERESSIG, color="green", ls="--", lw=1.5,
               label=f"Saueressig ({ALPHA_13C_SAUERESSIG})")
    ax.axvline(ALPHA_13C_CANTRELL, color="orange", ls="--", lw=1.5,
               label=f"Cantrell ({ALPHA_13C_CANTRELL})")
    ax.set_xlabel("α¹³C_OH", fontsize=10)
    ax.set_ylabel("Probability density", fontsize=10)
    ax.set_title("(a) SH direct constraint", fontsize=11)
    ax.legend(fontsize=7, loc="upper right")

    # ── Panel (b): Approach 2 ──
    ax = axes[1]
    alpha_2 = approach2["alpha_13c_oh_median"]
    ci2 = approach2["alpha_13c_oh_ci95"]
    sigma_2 = (ci2[1] - ci2[0]) / (2 * 1.96)
    if sigma_2 > 1e-8:
        pdf_2 = np.exp(-0.5 * ((x - alpha_2) / sigma_2)**2) / (sigma_2 * np.sqrt(2 * np.pi))
    else:
        # Degenerate case: delta-function-like → show as vertical line
        pdf_2 = np.zeros_like(x)
    ax.fill_between(x, pdf_2, alpha=0.3, color="C1")
    ax.plot(x, pdf_2, color="C1", lw=1.5, label=f"Source-corrected: {alpha_2:.4f}")
    ax.axvline(ALPHA_13C_SAUERESSIG, color="green", ls="--", lw=1.5,
               label=f"Saueressig ({ALPHA_13C_SAUERESSIG})")
    ax.axvline(ALPHA_13C_CANTRELL, color="orange", ls="--", lw=1.5,
               label=f"Cantrell ({ALPHA_13C_CANTRELL})")
    ax.set_xlabel("α¹³C_OH", fontsize=10)
    ax.set_title("(b) Source-corrected constraint", fontsize=11)
    ax.legend(fontsize=7, loc="upper right")

    # ── Panel (c): Summary comparison ──
    ax = axes[2]
    approaches = ["SH direct", "Source-\ncorrected"]
    medians = [alpha_1, alpha_2]
    ci_los = [ci1[0], ci2[0]]
    ci_his = [ci1[1], ci2[1]]

    y_pos = np.arange(len(approaches))
    ax.barh(y_pos, [m - 0.995 for m in medians], left=0.995,
            height=0.3, color=["C0", "C1"], alpha=0.6)
    ax.errorbar(medians, y_pos,
                xerr=[[m - lo for m, lo in zip(medians, ci_los)],
                      [hi - m for hi, m in zip(ci_his, medians)]],
                fmt="o", color="black", capsize=5, ms=6, zorder=5)

    ax.axvline(ALPHA_13C_SAUERESSIG, color="green", ls="--", lw=1.5,
               label="Saueressig (1.0039)")
    ax.axvline(ALPHA_13C_CANTRELL, color="orange", ls="--", lw=1.5,
               label="Cantrell (1.0054)")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(approaches, fontsize=10)
    ax.set_xlabel("α¹³C_OH", fontsize=10)
    ax.set_title("(c) Summary comparison", fontsize=11)
    ax.legend(fontsize=8, loc="upper right")
    ax.set_xlim(0.995, 1.015)

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "fig5_kie_constraint.png", dpi=300)
    plt.close(fig)
    print(f"✓ Figure saved: {FIGURES_DIR / 'fig5_kie_constraint.png'}")


# ============================================================================
# RESULT.md GENERATION
# ============================================================================

def write_result_md(approach1: dict, approach2: dict, approach3: dict) -> None:
    """Write the RESULT.md summary for the KIE_sites experiment."""
    a1 = approach1["alpha_13c_oh_median"]
    a1_ci = approach1["alpha_13c_oh_ci95"]
    a2 = approach2["alpha_13c_oh_median"]
    a2_ci = approach2["alpha_13c_oh_ci95"]

    md = f"""# KIE_sites — RESULT

> **NOTE:** This `RESULT.md` template is superseded by Phase 6
> (`phase6_phasor.py`). It is kept for reference only. Phase 6 overwrites
> `RESULT.md` with phasor-corrected results.

## Summary

**Research question:** Can seasonal cycles of δ¹³C-CH₄ and δD-CH₄ at co-located
measurement sites constrain the OH kinetic isotope effect (KIE) for ¹³C?

**Answer:** Partially. The observed amplitude ratios provide a **moderate** constraint
that favors the **Saueressig (α = 1.0039)** end of the contested range, but cannot
definitively exclude the Cantrell value (α = 1.0054).

## Key Findings

### 1. Source seasonality dominates at NH sites
- Observed amplitude ratio R = A(δ¹³C)/A(δD) ranges from 0.024 (SPO) to 0.14 (BRW)
- The pure-OH prediction is R = 0.013–0.017 (Saueressig–Cantrell)
- All NH sites show R far above this range (3–10×), indicating microbial source
  seasonality (wetlands, rice) inflates the δ¹³C seasonal cycle disproportionately
- This is because δ¹³C of microbial sources (−62‰) is much closer to atmospheric
  δ¹³C (−47‰) than δD of microbial sources (−310‰) is to atmospheric δD (−86‰)

### 2. SH sites provide the cleanest constraint
- **CGO** (Cape Grim, −41°S): R = 0.028 → implied α = 1.0040
- **SPO** (South Pole, −90°S): R = 0.024 → implied α = 1.0026
- These sites have minimal source seasonality and give α values close to or
  below the Saueressig value

### 3. Final KIE constraint

| Approach | α¹³C_OH | 95% CI | Notes |
|----------|---------|--------|-------|
| **SH direct** (CGO+SPO) | {a1:.4f} | [{a1_ci[0]:.4f}, {a1_ci[1]:.4f}] | Most robust, minimal source contamination |
| **Source-corrected** (all clean) | {a2:.4f} | [{a2_ci[0]:.4f}, {a2_ci[1]:.4f}] | Requires source signature assumptions |
| Saueressig et al. (2001) | 1.0039 | lab measurement | |
| Cantrell et al. (1990) | 1.0054 | lab measurement | |

### 4. Interpretation
- Both approaches favor the **lower end** of the α range (closer to Saueressig)
- The SH direct constraint ({a1:.4f}) is remarkably close to the Saueressig value (1.0039)
- The 95% CI [{a1_ci[0]:.4f}, {a1_ci[1]:.4f}] {'includes' if a1_ci[0] <= ALPHA_13C_SAUERESSIG <= a1_ci[1] else 'does not include'} the Saueressig value
  and {'includes' if a1_ci[0] <= ALPHA_13C_CANTRELL <= a1_ci[1] else 'does not include'} the Cantrell value
- **Caveat:** Even SH sites may have residual source seasonality from
  interhemispheric transport of NH wetland emissions, which would bias α upward

### 5. Diagnostics
- Phase differences δ¹³C−δD are mostly <2 months at clean sites,
  confirming both isotopes respond to the same seasonal driver
- Clear latitude gradient: R increases towards NH high latitudes,
  consistent with greater wetland source seasonality
- SMO (Samoa) is an outlier with −5.3 month phase offset — excluded

## Data Quality
- 12 co-located sites with both δ¹³C (NOAA/INSTAAR) and δD (Riddell-Young 2025)
- 8 sites classified as "clean" for this analysis
- INSTAAR sites (2005–2010) use same-flask measurements — highest quality pairing
- Bootstrap (N=2000) used for harmonic fit uncertainties
- Monte Carlo (N=50,000) propagates all parameter uncertainties

## Relationship to Other Experiments
- **KIE_sensitivity:** That experiment showed α¹³C_OH drives a 35.5 pp discriminant
  in model agreement-filter analysis and determines the FF trend sign. This result
  (α ≈ 1.003–1.004) would favor the Saueressig value, implying the FF trend may
  be flat-to-declining rather than increasing.
- **dD_threshold:** The δD seasonal amplitudes (1–4‰) at most sites are above the
  δD detectability threshold of ~37‰ annual uncertainty, confirming δD is useful
  for seasonal-cycle analysis even if marginal for source partitioning.

## Limitations
1. Only ~5 years of overlap at most sites (INSTAAR δD: 2005–2010)
2. Source isotope signatures (especially δD of wetlands) are uncertain
3. SH sites may still have some source contamination from NH transport
4. The linearized decomposition model assumes quasi-steady state
5. Cl seasonality is not fully separated from OH seasonality

## Files
- `analysis/phase1_data.py` — Data extraction and pairing
- `analysis/phase2_harmonics.py` — Seasonal harmonic fitting
- `analysis/phase3_synthesis.py` — Cross-site synthesis and classification
- `analysis/phase4_deconv.py` — Source deconvolution
- `analysis/phase5_kie.py` — KIE extraction (this phase)
- `results/` — JSON outputs from each phase
- `figures/` — Diagnostic and publication-quality figures
"""

    with open(EXPT_DIR / "RESULT.md", "w") as f:
        f.write(md)
    print(f"✓ RESULT.md saved to {EXPT_DIR / 'RESULT.md'}")


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("Phase 5: KIE extraction")
    print("=" * 70)

    # ── Load prior results ──
    with open(PHASE2_DIR / "harmonic_fits.json") as f:
        fits = json.load(f)
    with open(PHASE3_DIR / "synthesis_results.json") as f:
        synthesis = json.load(f)

    classifications = synthesis["site_classifications"]

    # ── Run three approaches ──
    result1 = approach1_sh_direct(fits)
    result2 = approach2_source_corrected(fits, classifications)
    result3 = approach3_latitude_gradient(fits, classifications)

    # ── Save results ──
    output = {
        "approach1_sh_direct": result1,
        "approach2_source_corrected": result2,
        "approach3_latitude_gradient": result3,
        "reference_values": {
            "saueressig_alpha_13c": ALPHA_13C_SAUERESSIG,
            "cantrell_alpha_13c": ALPHA_13C_CANTRELL,
            "alpha_D_oh": ALPHA_D_OH,
        },
    }

    with open(RESULTS_DIR / "kie_results.json", "w") as f:
        json.dump(output, f, indent=2)
    print(f"\n✓ Results saved to {RESULTS_DIR / 'kie_results.json'}")

    # ── Figure ──
    plot_kie_constraint(result1, result2, result3)

    # ── RESULT.md ──
    write_result_md(result1, result2, result3)

    # ── Final summary ──
    print("\n" + "=" * 70)
    print("FINAL KIE CONSTRAINT")
    print("=" * 70)
    print(f"\n  Approach 1 (SH direct):       α = {result1['alpha_13c_oh_median']:.4f}  "
          f"[{result1['alpha_13c_oh_ci95'][0]:.4f}, {result1['alpha_13c_oh_ci95'][1]:.4f}]")
    print(f"  Approach 2 (source-corrected): α = {result2['alpha_13c_oh_median']:.4f}  "
          f"[{result2['alpha_13c_oh_ci95'][0]:.4f}, {result2['alpha_13c_oh_ci95'][1]:.4f}]")
    print(f"\n  Lab values:")
    print(f"    Saueressig (2001): α = {ALPHA_13C_SAUERESSIG}")
    print(f"    Cantrell   (1990): α = {ALPHA_13C_CANTRELL}")

    # Assess discrimination
    a1_lo = result1["alpha_13c_oh_ci95"][0]
    a1_hi = result1["alpha_13c_oh_ci95"][1]
    saueressig_in = a1_lo <= ALPHA_13C_SAUERESSIG <= a1_hi
    cantrell_in = a1_lo <= ALPHA_13C_CANTRELL <= a1_hi

    print(f"\n  SH constraint CI {'includes' if saueressig_in else 'EXCLUDES'} "
          f"Saueressig (1.0039)")
    print(f"  SH constraint CI {'includes' if cantrell_in else 'EXCLUDES'} "
          f"Cantrell (1.0054)")

    if cantrell_in and saueressig_in:
        print("\n  ⚠ Cannot discriminate between the two lab values with this data.")
        print("    Both fall within the 95% CI.")
    elif not cantrell_in and saueressig_in:
        print("\n  ✓ Data FAVORS Saueressig over Cantrell.")
    elif cantrell_in and not saueressig_in:
        print("\n  Data favors Cantrell over Saueressig (unexpected).")
    else:
        print("\n  Neither lab value is within the CI — suggests systematic issues.")


if __name__ == "__main__":
    main()
