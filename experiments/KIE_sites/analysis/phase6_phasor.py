#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
phase6_phasor.py — Phasor source correction for the KIE_sites experiment
=========================================================================

Removes wetland source contamination from the observed seasonal amplitude
ratio R = A(δ¹³C)/A(δD) using vector (phasor) subtraction, then inverts
the corrected ratio to constrain the OH ¹³C KIE (α_13C_OH).

Background
----------
The observed seasonal cycle of each isotope is a vector sum of sink-driven
enrichment (OH + other sinks) and source-driven depletion (wetlands):

    Z_obs = Z_sink + Z_source         (complex phasor, ‰)

Because OH (~July peak) and wetlands (~July–Aug peak) have different phases,
simple scalar subtraction is incorrect. We use phasor (vector) decomposition
in the complex B + iC plane, where B and C are the sin/cos harmonic
coefficients from Phase 2.

Algorithm
---------
For each clean site and each isotope (δ¹³C, δD):

  1. Z_obs   = B_obs + i·C_obs                            [from Phase 2]
  2. Z_frac  = (B_Q + i·C_Q) / Q_total                    [wetland fractional seasonality]
  3. Z_src   = (δ_source − δ_atm) × Z_frac                [source phasor, ‰]
  4. Z_sink  = Z_obs − Z_src                               [vector subtraction]

  R_corrected = |Z_sink(δ¹³C)| / |Z_sink(δD)|
  → α_13C_OH via ratio_to_alpha_13c()

Consistency check: arg(Z_sink) should be ≈ same month for δ¹³C and δD
                   (both driven by OH, which peaks ~July in NH).

Source region assignment
------------------------
NH sites see nearby-band wetland emissions directly. SH sites (CGO, SPO)
see a mix of local SH wetlands (tiny) and NH emissions attenuated by
interhemispheric transport. A 2-box model gives:
  - NH→SH attenuation: ~1/(1 + τ_mix·ω)² ≈ 0.11 for τ_mix = 1.3 yr
  - Phase shift: ~arctan(τ_mix·ω) ≈ 2.8 months
Rather than model this imprecisely, we assign SH sites to the SH_extra
band (local sources only). This is conservative: any transported NH signal
adds ≲10% of the NH amplitude, and its ~3-month phase shift means it
partially cancels the local SH signal. The effect on R is <0.002.

Inputs
------
  results/phase2_harmonics/harmonic_fits.json   — observed B, C per site
  data/wetland_seasonality.json                 — wetland B_Q, C_Q per site
  data/dD_source_database.json                  — site-specific δD_wetland

Outputs
-------
  results/phase6_phasor/phasor_results.json     — corrected R, α_13C_OH per site
  figures/fig8_phasor_decomposition.png          — vector diagrams
  figures/fig9_corrected_ratio.png               — R_corr vs latitude
  figures/fig10_alpha_constraint.png             — final KIE constraint

References
----------
  Li et al. (2026) ESSD — wetland emission seasonality
  Douglas et al. (2021) Biogeosciences — δD source signatures
  Saueressig et al. (2001) — α_13C_OH = 1.0039
  Cantrell et al. (1990) — α_13C_OH = 1.0054
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
EXPT_DIR    = Path(__file__).resolve().parent.parent
PHASE2_JSON = EXPT_DIR / "results" / "phase2_harmonics" / "harmonic_fits.json"
WETLAND_JSON = EXPT_DIR / "data" / "wetland_seasonality.json"
DD_SRC_JSON  = EXPT_DIR / "data" / "dD_source_database.json"
SYNTH_JSON   = EXPT_DIR / "results" / "phase3_synthesis" / "synthesis_results.json"
OUT_DIR      = EXPT_DIR / "results" / "phase6_phasor"
FIG_DIR      = EXPT_DIR / "figures"

OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

OUT_JSON = OUT_DIR / "phasor_results.json"

# ============================================================================
# PHYSICAL CONSTANTS  (identical to Phase 5 for consistency)
# ============================================================================
# Atmospheric isotope composition
D13C_ATM = -47.3   # ‰ VPDB
DD_ATM   = -86.0   # ‰ VSMOW

# Wetland source δ¹³C (global; all sites are C₃-dominated)
D13C_WETLAND       = -62.0    # ‰
D13C_WETLAND_SIGMA =   5.0    # ‰

# Total global CH₄ source (Saunois 2025)
Q_TOTAL_TG_YR       = 580.0    # Tg/yr
Q_TOTAL_TG_YR_SIGMA =  50.0    # Tg/yr
Q_TOTAL_TG_MONTH    = Q_TOTAL_TG_YR / 12.0   # 48.3 Tg/month

# Wetland emission uncertainty (Li2026 ensemble spread)
WETLAND_BC_FRAC_SIGMA = 0.20   # ±20% on B_Q, C_Q

# Sink fractionation and fractions (same as Phase 5)
ALPHA_D_OH     = 1.294
ALPHA_13C_CL   = 1.066;   ALPHA_D_CL   = 1.508
ALPHA_13C_SOIL = 1.022;   ALPHA_D_SOIL = 1.066
ALPHA_13C_STRAT= 1.013;   ALPHA_D_STRAT= 1.16

F_OH   = 0.84;    SIGMA_F_OH   = 0.04
F_CL   = 0.035;   SIGMA_F_CL   = 0.01
F_SOIL = 0.06;    SIGMA_F_SOIL = 0.02
F_STRAT= 0.065

SIGMA_ALPHA_D_OH    = 0.01
SIGMA_ALPHA_13C_CL  = 0.005
SIGMA_ALPHA_D_CL    = 0.05

# Saueressig / Cantrell reference values
ALPHA_13C_SAUERESSIG = 1.0039
ALPHA_13C_CANTRELL   = 1.0054

# Clean sites (from Phase 3)
CLEAN_SITES = ["ALT", "ZEP", "BRW", "CBA", "MHD", "KUM", "CGO", "SPO"]

# ============================================================================
# HARMONIC CONVENTION CONVERSION
# ============================================================================
# Phase 2 uses:   δ(t) = c0 + c1·(t−t_ref) + B·sin(2πt) + C·cos(2πt)
#                 where t is in fractional years (t = 0.0 = Jan 1, t = 0.5 = Jul 1)
#
# Wetland file uses: Q(m) = Q̄ + B_Q·sin(2πm/12) + C_Q·cos(2πm/12)
#                    where m = 0..11 (Jan=0, Dec=11)
#
# The bases have the same annual frequency, but Phase 1 stores monthly isotope
# means at month midpoints:
# (month_index + 0.5) / 12, while wetland climatology coefficients are fit
# against month_index / 12. Rotate wetland B_Q,C_Q into the Phase 2 midpoint
# basis before phasor subtraction.


# ============================================================================
# HELPER: ratio_to_alpha_13c  (reproduced from Phase 5 for self-containment)
# ============================================================================
def ratio_to_alpha_13c(R_obs, alpha_d_oh=ALPHA_D_OH,
                        f_oh=F_OH, f_cl=F_CL, f_soil=F_SOIL, f_strat=F_STRAT,
                        alpha_13c_cl=ALPHA_13C_CL, alpha_d_cl=ALPHA_D_CL):
    """Invert amplitude ratio R = ε_13C_bulk / ε_D_bulk → α_13C_OH.

    Given all non-OH KIE values, solve:
        ε_13C_bulk = R × ε_D_bulk
        α_13C_OH   = 1 + (ε_13C_bulk − ε_13C_non_OH) / (f_OH × 1000)
    """
    # Bulk ε_D (does not depend on α_13C_OH)
    eps_D = (f_oh * (alpha_d_oh - 1) +
             f_cl * (alpha_d_cl - 1) +
             f_soil * (ALPHA_D_SOIL - 1) +
             f_strat * (ALPHA_D_STRAT - 1)) * 1000.0

    eps_13c_needed = R_obs * eps_D

    eps_13c_non_oh = (f_cl * (alpha_13c_cl - 1) +
                      f_soil * (ALPHA_13C_SOIL - 1) +
                      f_strat * (ALPHA_13C_STRAT - 1)) * 1000.0

    return 1.0 + (eps_13c_needed - eps_13c_non_oh) / (f_oh * 1000.0)


def ci_to_sigma(ci_lo, ci_hi):
    """Approximate 1σ from 95% CI (assuming normal distribution)."""
    return (ci_hi - ci_lo) / (2 * 1.96)


def summarize_alpha_filters(alpha_samples):
    """Compare legacy narrow and wider alpha filters for pooled MC samples."""
    finite_alpha = np.asarray(alpha_samples)
    finite_alpha = finite_alpha[np.isfinite(finite_alpha)]

    def summarize_window(lo, hi):
        filtered = finite_alpha[(finite_alpha > lo) & (finite_alpha < hi)]
        if len(filtered) == 0:
            return {
                "filter": [lo, hi],
                "alpha_13c_oh_median": np.nan,
                "alpha_13c_oh_ci95": [np.nan, np.nan],
                "n_samples": 0,
                "n_excluded": int(len(finite_alpha)),
            }
        return {
            "filter": [lo, hi],
            "alpha_13c_oh_median": float(np.median(filtered)),
            "alpha_13c_oh_ci95": [
                float(np.percentile(filtered, 2.5)),
                float(np.percentile(filtered, 97.5)),
            ],
            "n_samples": int(len(filtered)),
            "n_excluded": int(len(finite_alpha) - len(filtered)),
        }

    narrow = summarize_window(0.99, 1.02)
    wide = summarize_window(0.98, 1.05)
    return {
        "narrow": narrow,
        "wide": wide,
        "impact_wide_minus_narrow": {
            "median": float(wide["alpha_13c_oh_median"] - narrow["alpha_13c_oh_median"]),
            "ci95": [
                float(wide["alpha_13c_oh_ci95"][0] - narrow["alpha_13c_oh_ci95"][0]),
                float(wide["alpha_13c_oh_ci95"][1] - narrow["alpha_13c_oh_ci95"][1]),
            ],
            "n_samples": int(wide["n_samples"] - narrow["n_samples"]),
        },
    }


def phasor_peak_month(B, C):
    """Convert harmonic B, C to peak month (0 = Jan, 6 = Jul).

    For signal = B·sin(2πt) + C·cos(2πt), peak at t where derivative = 0:
        peak_t = atan2(C, B) / (2π)  →  then mod 1.0
    But sin(2πt) peaks at t = 0.25, so:
        peak_t = (0.25 − atan2(C, B)/(2π)) mod 1.0
    → month = peak_t × 12
    """
    phase = np.arctan2(C, B)
    peak_t = (0.25 - phase / (2 * np.pi)) % 1.0
    return peak_t * 12.0


def convert_wetland_to_phase2_phasor(B_Q, C_Q):
    """Rotate wetland month-index coefficients into Phase 2 midpoint time.

    Wetland coefficients are fit against month_index / 12. Phase 2 monthly
    isotope means are fit at (month_index + 0.5) / 12. For

        Q = B_Q*sin(w*x) + C_Q*cos(w*x)

    with x = month_index/12 and t = x + 0.5/12, the equivalent coefficients
    in the Phase 2 basis are R(-pi/12) applied to (B_Q, C_Q).
    """
    delta = 2 * np.pi * 0.5 / 12.0
    cos_d = np.cos(delta)
    sin_d = np.sin(delta)
    B_mid = B_Q * cos_d + C_Q * sin_d
    C_mid = -B_Q * sin_d + C_Q * cos_d
    return B_mid, C_mid


# ============================================================================
# LOAD DATA
# ============================================================================
def load_all():
    """Load Phase 2 harmonics, wetland seasonality, δD source database."""
    with open(PHASE2_JSON) as f:
        fits = json.load(f)
    with open(WETLAND_JSON) as f:
        wetland = json.load(f)
    with open(DD_SRC_JSON) as f:
        dd_db = json.load(f)
    with open(SYNTH_JSON) as f:
        synth = json.load(f)
    return fits, wetland, dd_db, synth


# ============================================================================
# PHASOR DECOMPOSITION (single site, deterministic)
# ============================================================================
def phasor_decompose(B_obs_13c, C_obs_13c,
                     B_obs_dD, C_obs_dD,
                     B_Q, C_Q, Q_mean,
                     Q_total,
                     d13C_wetland, dD_wetland):
    """Perform phasor source subtraction for one site.

    Parameters
    ----------
    B_obs_13c, C_obs_13c : harmonic coefficients for observed δ¹³C (‰)
    B_obs_dD, C_obs_dD   : harmonic coefficients for observed δD (‰)
    B_Q, C_Q             : wetland emission harmonic (Tg/month)
    Q_mean               : mean wetland emission for this band (Tg/month)
    Q_total              : total global CH₄ source (Tg/month)
    d13C_wetland         : δ¹³C of wetland CH₄ (‰)
    dD_wetland           : δD of wetland CH₄ (‰)

    Returns
    -------
    dict with corrected amplitudes, phases, ratio
    """
    # Observed phasors  (complex: B + iC)
    Z_obs_13c = complex(B_obs_13c, C_obs_13c)
    Z_obs_dD  = complex(B_obs_dD,  C_obs_dD)

    # Fractional wetland seasonality phasor, converted to Phase 2 midpoint time.
    # Z_frac = (B_Q + i·C_Q) / Q_total
    B_Q_mid, C_Q_mid = convert_wetland_to_phase2_phasor(B_Q, C_Q)
    Z_frac = complex(B_Q_mid, C_Q_mid) / Q_total

    # Source phasors  (‰)
    # Source pulls δ_atm toward δ_source → seasonal Δδ = (δ_src − δ_atm) × ΔS/S
    gap_13c = d13C_wetland - D13C_ATM    # negative (−62 − (−47.3) = −14.7)
    gap_dD  = dD_wetland   - DD_ATM      # negative (e.g., −374 − (−86) = −288)

    Z_src_13c = gap_13c * Z_frac
    Z_src_dD  = gap_dD  * Z_frac

    # Vector subtraction: sink = observed − source
    Z_sink_13c = Z_obs_13c - Z_src_13c
    Z_sink_dD  = Z_obs_dD  - Z_src_dD

    # Corrected amplitudes and ratio
    A_sink_13c = abs(Z_sink_13c)
    A_sink_dD  = abs(Z_sink_dD)
    R_corrected = A_sink_13c / A_sink_dD if A_sink_dD > 0 else np.nan

    # Peak months
    peak_obs_13c  = phasor_peak_month(B_obs_13c, C_obs_13c)
    peak_obs_dD   = phasor_peak_month(B_obs_dD,  C_obs_dD)
    peak_sink_13c = phasor_peak_month(Z_sink_13c.real, Z_sink_13c.imag)
    peak_sink_dD  = phasor_peak_month(Z_sink_dD.real, Z_sink_dD.imag)
    peak_src_13c  = phasor_peak_month(Z_src_13c.real, Z_src_13c.imag)
    peak_src_dD   = phasor_peak_month(Z_src_dD.real, Z_src_dD.imag)

    return {
        # Observed
        "Z_obs_13c": [Z_obs_13c.real, Z_obs_13c.imag],
        "Z_obs_dD":  [Z_obs_dD.real,  Z_obs_dD.imag],
        "A_obs_13c": abs(Z_obs_13c),
        "A_obs_dD":  abs(Z_obs_dD),
        "R_obs":     abs(Z_obs_13c) / abs(Z_obs_dD) if abs(Z_obs_dD) > 0 else np.nan,
        "peak_obs_13c": peak_obs_13c,
        "peak_obs_dD":  peak_obs_dD,
        # Source
        "Z_src_13c": [Z_src_13c.real, Z_src_13c.imag],
        "Z_src_dD":  [Z_src_dD.real,  Z_src_dD.imag],
        "A_src_13c": abs(Z_src_13c),
        "A_src_dD":  abs(Z_src_dD),
        "peak_src_13c": peak_src_13c,
        "peak_src_dD":  peak_src_dD,
        "gap_13c":   gap_13c,
        "gap_dD":    gap_dD,
        # Sink (corrected)
        "Z_sink_13c": [Z_sink_13c.real, Z_sink_13c.imag],
        "Z_sink_dD":  [Z_sink_dD.real,  Z_sink_dD.imag],
        "A_sink_13c": A_sink_13c,
        "A_sink_dD":  A_sink_dD,
        "R_corrected": R_corrected,
        "peak_sink_13c": peak_sink_13c,
        "peak_sink_dD":  peak_sink_dD,
        "sink_phase_diff_months": float(min(abs(peak_sink_13c - peak_sink_dD),
                                             12 - abs(peak_sink_13c - peak_sink_dD))),
    }


# ============================================================================
# MONTE CARLO UNCERTAINTY
# ============================================================================
def mc_phasor(B_obs_13c, C_obs_13c, amp_ci_13c, peak_ci_13c,
              B_obs_dD, C_obs_dD, amp_ci_dD, peak_ci_dD,
              B_Q, C_Q, Q_mean,
              dD_wetland, dD_sigma,
              n_mc=50000, rng=None):
    """Monte Carlo uncertainty propagation for phasor correction at one site.

    Draws from:
      - Observed amplitude ± CI  and  phase ± CI  (independently)
      - δD_wetland ± sigma
      - δ¹³C_wetland ± 5 ‰
      - Q_total ± 50 Tg/yr
      - Wetland B_Q, C_Q ± 20%
      - Non-OH KIE parameters (same as Phase 5)

    Returns array of R_corrected samples and α_13C_OH samples.
    """
    if rng is None:
        rng = np.random.default_rng(42)

    # Observed amplitude and phase uncertainties
    A_obs_13c = np.sqrt(B_obs_13c**2 + C_obs_13c**2)
    A_obs_dD  = np.sqrt(B_obs_dD**2  + C_obs_dD**2)
    sigma_amp_13c = ci_to_sigma(amp_ci_13c[0], amp_ci_13c[1])
    sigma_amp_dD  = ci_to_sigma(amp_ci_dD[0],  amp_ci_dD[1])

    # Phase uncertainty: from peak_month CI → convert to radians (1 month ≈ 2π/12)
    phase_obs_13c = np.arctan2(C_obs_13c, B_obs_13c)
    phase_obs_dD  = np.arctan2(C_obs_dD,  B_obs_dD)
    sigma_phase_13c = ci_to_sigma(peak_ci_13c[0], peak_ci_13c[1]) * (2 * np.pi / 12)
    sigma_phase_dD  = ci_to_sigma(peak_ci_dD[0],  peak_ci_dD[1])  * (2 * np.pi / 12)

    # Draw observed B, C with independent amplitude and phase perturbation.
    a13 = np.maximum(rng.normal(A_obs_13c, sigma_amp_13c, n_mc), 1e-6)
    p13 = phase_obs_13c + rng.normal(0, sigma_phase_13c, n_mc)
    b13 = a13 * np.cos(p13)
    c13 = a13 * np.sin(p13)

    aD = np.maximum(rng.normal(A_obs_dD, sigma_amp_dD, n_mc), 1e-6)
    pD = phase_obs_dD + rng.normal(0, sigma_phase_dD, n_mc)
    bD = aD * np.cos(pD)
    cD = aD * np.sin(pD)

    # Draw wetland parameters.
    d13c_w = rng.normal(D13C_WETLAND, D13C_WETLAND_SIGMA, n_mc)
    dD_w = rng.normal(dD_wetland, dD_sigma, n_mc)
    q_tot = np.maximum(rng.normal(Q_TOTAL_TG_YR, Q_TOTAL_TG_YR_SIGMA, n_mc), 300.0) / 12.0
    bq = B_Q * (1 + rng.normal(0, WETLAND_BC_FRAC_SIGMA, n_mc))
    cq = C_Q * (1 + rng.normal(0, WETLAND_BC_FRAC_SIGMA, n_mc))
    bq_mid, cq_mid = convert_wetland_to_phase2_phasor(bq, cq)

    # Phasor decomposition.
    Z_obs_13c = b13 + 1j * c13
    Z_obs_dD = bD + 1j * cD
    Z_frac = (bq_mid + 1j * cq_mid) / q_tot
    Z_src_13c = (d13c_w - D13C_ATM) * Z_frac
    Z_src_dD = (dD_w - DD_ATM) * Z_frac
    Z_sink_13c = Z_obs_13c - Z_src_13c
    Z_sink_dD = Z_obs_dD - Z_src_dD

    A_sink_13c = np.abs(Z_sink_13c)
    A_sink_dD = np.abs(Z_sink_dD)

    R_samples = np.divide(
        A_sink_13c,
        A_sink_dD,
        out=np.full(n_mc, np.nan),
        where=A_sink_dD > 1e-6,
    )

    # Draw non-OH KIE parameters.
    f_oh = np.clip(rng.normal(F_OH, SIGMA_F_OH, n_mc), 0.5, 0.99)
    f_cl = np.clip(rng.normal(F_CL, SIGMA_F_CL, n_mc), 0.0, 0.1)
    f_soil = np.clip(rng.normal(F_SOIL, SIGMA_F_SOIL, n_mc), 0.0, 0.15)
    f_strat = 1.0 - f_oh - f_cl - f_soil
    alpha_d_oh = rng.normal(ALPHA_D_OH, SIGMA_ALPHA_D_OH, n_mc)
    alpha_13c_cl = rng.normal(ALPHA_13C_CL, SIGMA_ALPHA_13C_CL, n_mc)
    alpha_d_cl = rng.normal(ALPHA_D_CL, SIGMA_ALPHA_D_CL, n_mc)

    alpha_samples = ratio_to_alpha_13c(
        R_samples, alpha_d_oh, f_oh, f_cl, f_soil, f_strat, alpha_13c_cl, alpha_d_cl)

    return R_samples, alpha_samples

# ============================================================================
# MAIN ANALYSIS
# ============================================================================
def main():
    print("=" * 70)
    print("Phase 6: Phasor source correction")
    print("=" * 70)

    fits, wetland, dd_db, synth = load_all()

    results = {
        "metadata": {
            "method": "Phasor (vector) subtraction of wetland source seasonality",
            "inputs": ["phase2_harmonics", "wetland_seasonality", "dD_source_database"],
            "convention": "Z = B + iC, same as Phase 2 (B sin + C cos, fractional year)",
            "n_mc": 50000,
            "date": "2026-05-20",
        },
        "sites": {},
        "multi_site_result": {},
    }

    # ── Per-site analysis ──
    print(f"\n{'Site':<5} {'R_obs':>7} {'R_corr':>7} {'Δ':>7} "
          f"{'pk_obs13':>8} {'pk_obD':>7} {'pk_sk13':>7} {'pk_skD':>7} "
          f"{'A_src13':>8} {'A_srcD':>7}")
    print("-" * 85)

    valid_R_corr = []
    valid_alpha  = []
    valid_sigmas = []
    valid_codes  = []

    for code in CLEAN_SITES:
        if code not in fits or code not in wetland["site_assignment"]:
            continue

        site_fit = fits[code]
        site_wet = wetland["site_assignment"][code]
        site_dd  = dd_db["sites"].get(code, {})

        # Observed B, C
        B_obs_13c = site_fit["d13C"]["B"]
        C_obs_13c = site_fit["d13C"]["C"]
        B_obs_dD  = site_fit["dD"]["B"]
        C_obs_dD  = site_fit["dD"]["C"]

        # Wetland harmonic coefficients (Tg/month)
        B_Q    = site_wet["B_Q_Tg_month"]
        C_Q    = site_wet["C_Q_Tg_month"]
        Q_mean = site_wet["Q_mean_Tg_month"]

        # δD source signature
        dD_wetland = site_dd["recommended"]["dD_CH4"]
        dD_sigma   = site_dd["recommended"]["sigma"]

        # ── Deterministic decomposition ──
        decomp = phasor_decompose(
            B_obs_13c, C_obs_13c,
            B_obs_dD, C_obs_dD,
            B_Q, C_Q, Q_mean,
            Q_TOTAL_TG_MONTH,
            D13C_WETLAND, dD_wetland,
        )

        # ── Monte Carlo ──
        amp_ci_13c  = site_fit["d13C"]["amplitude_ci95"]
        amp_ci_dD   = site_fit["dD"]["amplitude_ci95"]
        peak_ci_13c = site_fit["d13C"].get("peak_month_ci95", [site_fit["d13C"]["peak_month"] - 1.0,
                                                                 site_fit["d13C"]["peak_month"] + 1.0])
        peak_ci_dD  = site_fit["dD"].get("peak_month_ci95",  [site_fit["dD"]["peak_month"] - 1.0,
                                                                site_fit["dD"]["peak_month"] + 1.0])

        R_mc, alpha_mc = mc_phasor(
            B_obs_13c, C_obs_13c, amp_ci_13c, peak_ci_13c,
            B_obs_dD, C_obs_dD, amp_ci_dD, peak_ci_dD,
            B_Q, C_Q, Q_mean,
            dD_wetland, dD_sigma,
        )

        # Filter valid MC samples (finite, positive ratio)
        mask = np.isfinite(R_mc) & np.isfinite(alpha_mc) & (R_mc > 0) & (R_mc < 1)
        R_mc_valid     = R_mc[mask]
        alpha_mc_valid = alpha_mc[mask]

        R_med    = float(np.median(R_mc_valid)) if len(R_mc_valid) > 100 else np.nan
        R_lo, R_hi = (float(np.percentile(R_mc_valid, 2.5)),
                       float(np.percentile(R_mc_valid, 97.5))) if len(R_mc_valid) > 100 else (np.nan, np.nan)
        a_med    = float(np.median(alpha_mc_valid)) if len(alpha_mc_valid) > 100 else np.nan
        a_lo, a_hi = (float(np.percentile(alpha_mc_valid, 2.5)),
                       float(np.percentile(alpha_mc_valid, 97.5))) if len(alpha_mc_valid) > 100 else (np.nan, np.nan)

        R_sigma = ci_to_sigma(R_lo, R_hi) if np.isfinite(R_lo) else np.nan

        print(f"{code:<5} {decomp['R_obs']:>7.4f} {decomp['R_corrected']:>7.4f} "
              f"{decomp['R_corrected'] - decomp['R_obs']:>+7.4f} "
              f"{decomp['peak_obs_13c']:>8.1f} {decomp['peak_obs_dD']:>7.1f} "
              f"{decomp['peak_sink_13c']:>7.1f} {decomp['peak_sink_dD']:>7.1f} "
              f"{decomp['A_src_13c']:>8.4f} {decomp['A_src_dD']:>7.3f}")

        # Store results
        site_result = {
            "source_band": site_wet["source_band"],
            "dD_wetland": dD_wetland,
            "dD_sigma": dD_sigma,
            **{k: (round(v, 6) if isinstance(v, float) else v)
               for k, v in decomp.items()},
            "mc": {
                "R_corrected_median": R_med,
                "R_corrected_ci95": [R_lo, R_hi],
                "alpha_13c_oh_median": a_med,
                "alpha_13c_oh_ci95": [a_lo, a_hi],
                "n_valid": int(np.sum(mask)),
            },
        }
        results["sites"][code] = site_result

        # Collect for multi-site weighted mean
        if np.isfinite(R_med) and np.isfinite(R_sigma) and R_sigma > 0:
            valid_R_corr.append(R_med)
            valid_alpha.append(alpha_mc_valid)
            valid_sigmas.append(R_sigma)
            valid_codes.append(code)

    # ── Multi-site weighted mean ──
    print("\n" + "=" * 70)
    print("MULTI-SITE PHASOR-CORRECTED RESULT")
    print("=" * 70)

    if len(valid_R_corr) >= 2:
        R_arr = np.array(valid_R_corr)
        sig_arr = np.array(valid_sigmas)
        w = 1.0 / sig_arr**2
        R_wm  = np.sum(w * R_arr) / np.sum(w)
        R_wm_sig = np.sqrt(1.0 / np.sum(w))

        # Pool all alpha MC samples for combined constraint. Report the legacy
        # narrow filter and a wider conservative filter so tail clipping is visible.
        alpha_filter_summary = summarize_alpha_filters(np.concatenate(valid_alpha))
        alpha_wide = alpha_filter_summary["wide"]
        alpha_narrow = alpha_filter_summary["narrow"]
        alpha_impact = alpha_filter_summary["impact_wide_minus_narrow"]

        a_wm_med = alpha_wide["alpha_13c_oh_median"]
        a_wm_lo, a_wm_hi = alpha_wide["alpha_13c_oh_ci95"]

        print(f"  Sites used: {valid_codes}")
        print(f"  Weighted mean R_corrected = {R_wm:.4f} ± {R_wm_sig:.4f}")
        print(f"  Pooled α_13C_OH (wide filter 0.98-1.05) = "
              f"{a_wm_med:.4f}  [{a_wm_lo:.4f}, {a_wm_hi:.4f}] (95% CI)")
        print(f"  Legacy narrow filter 0.99-1.02 = "
              f"{alpha_narrow['alpha_13c_oh_median']:.4f}  "
              f"[{alpha_narrow['alpha_13c_oh_ci95'][0]:.4f}, "
              f"{alpha_narrow['alpha_13c_oh_ci95'][1]:.4f}]")
        print(f"  Wide-minus-narrow impact: median {alpha_impact['median']:+.4f}, "
              f"CI [{alpha_impact['ci95'][0]:+.4f}, {alpha_impact['ci95'][1]:+.4f}], "
              f"+{alpha_impact['n_samples']} samples")
        print(f"  Saueressig = {ALPHA_13C_SAUERESSIG:.4f}")
        print(f"  Cantrell   = {ALPHA_13C_CANTRELL:.4f}")

        results["multi_site_result"] = {
            "sites_used": valid_codes,
            "R_weighted_mean": round(float(R_wm), 6),
            "R_weighted_sigma": round(float(R_wm_sig), 6),
            "alpha_13c_oh_median": round(a_wm_med, 6),
            "alpha_13c_oh_ci95": [round(a_wm_lo, 6), round(a_wm_hi, 6)],
            "n_pooled_samples": int(alpha_wide["n_samples"]),
            "alpha_filter_used": "wide_0.98_1.05",
            "alpha_filter_sensitivity": alpha_filter_summary,
        }

        # ── SH-only subset ──
        sh_codes = [c for c in valid_codes if c in ("CGO", "SPO")]
        if len(sh_codes) >= 1:
            sh_filter_summary = summarize_alpha_filters(
                np.concatenate([valid_alpha[valid_codes.index(c)] for c in sh_codes])
            )
            sh_wide = sh_filter_summary["wide"]
            if sh_wide["n_samples"] > 100:
                results["multi_site_result"]["sh_only"] = {
                    "sites": sh_codes,
                    "alpha_13c_oh_median": round(sh_wide["alpha_13c_oh_median"], 6),
                    "alpha_13c_oh_ci95": [round(sh_wide["alpha_13c_oh_ci95"][0], 6),
                                           round(sh_wide["alpha_13c_oh_ci95"][1], 6)],
                    "alpha_filter_used": "wide_0.98_1.05",
                    "alpha_filter_sensitivity": sh_filter_summary,
                }
                print(f"\n  SH-only ({sh_codes}, wide filter): α = "
                      f"{sh_wide['alpha_13c_oh_median']:.4f} "
                      f"[{sh_wide['alpha_13c_oh_ci95'][0]:.4f}, "
                      f"{sh_wide['alpha_13c_oh_ci95'][1]:.4f}]")

    # ── Save ──
    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2, default=lambda x: float(x) if isinstance(x, (np.floating, np.integer)) else x)
    print(f"\n✓ Saved: {OUT_JSON}")

    # ── Figures ──
    plot_phasor_diagrams(results)
    plot_corrected_ratio(results, fits)
    plot_alpha_constraint(results)


# ============================================================================
# FIGURE 8: PHASOR VECTOR DIAGRAMS
# ============================================================================
def plot_phasor_diagrams(results):
    """Fig 8: Phasor decomposition for 4 representative sites.

    Each panel shows three vectors from origin:
      green  = Z_obs   (observed)
      red    = Z_src   (wetland source)
      blue   = Z_sink  (corrected = obs − src), also shown from src tip to obs tip
    Axes are auto-scaled per panel to fit all vectors.
    """
    show_sites = ["BRW", "CBA", "CGO", "SPO"]
    available = [s for s in show_sites if s in results["sites"]]
    n = len(available)
    if n == 0:
        return

    fig, axes = plt.subplots(2, n, figsize=(4.5 * n, 8))
    if n == 1:
        axes = axes.reshape(2, 1)

    for j, code in enumerate(available):
        s = results["sites"][code]

        for row, iso, label in [(0, "13c", "δ¹³C"), (1, "dD", "δD")]:
            ax = axes[row, j]

            Z_obs  = s[f"Z_obs_{iso}"]
            Z_src  = s[f"Z_src_{iso}"]
            Z_sink = s[f"Z_sink_{iso}"]

            # Collect all endpoints to determine axis limits
            all_x = [0, Z_obs[0], Z_src[0], Z_sink[0]]
            all_y = [0, Z_obs[1], Z_src[1], Z_sink[1]]
            pad_x = max(0.01, (max(all_x) - min(all_x)) * 0.25)
            pad_y = max(0.01, (max(all_y) - min(all_y)) * 0.25)

            # ── Draw vectors as arrows ──
            hw = 0.04   # relative head width

            # Observed (green): origin → obs
            ax.annotate("", xy=Z_obs, xytext=(0, 0),
                        arrowprops=dict(arrowstyle="-|>", color="green", lw=2.5,
                                        mutation_scale=15))
            # Source (red dashed): origin → src
            ax.annotate("", xy=Z_src, xytext=(0, 0),
                        arrowprops=dict(arrowstyle="-|>", color="red", lw=2,
                                        mutation_scale=15, linestyle="dashed"))
            # Sink from origin (blue dotted): origin → sink
            ax.annotate("", xy=Z_sink, xytext=(0, 0),
                        arrowprops=dict(arrowstyle="-|>", color="C0", lw=1.5,
                                        mutation_scale=12, linestyle="dotted"))
            # Sink as vector addition (blue solid): src tip → obs tip
            # This shows Z_obs = Z_src + Z_sink visually
            ax.annotate("", xy=Z_obs, xytext=Z_src,
                        arrowprops=dict(arrowstyle="-|>", color="C0", lw=2.5,
                                        mutation_scale=15))

            # Dot labels at vector tips
            ax.plot(*Z_obs,  "go", ms=5, zorder=5)
            ax.plot(*Z_src,  "r^", ms=5, zorder=5)
            ax.plot(*Z_sink, "bs", ms=5, zorder=5)
            ax.plot(0, 0, "ko", ms=5, zorder=5)

            # Annotate amplitudes
            A_obs  = np.sqrt(Z_obs[0]**2  + Z_obs[1]**2)
            A_src  = np.sqrt(Z_src[0]**2  + Z_src[1]**2)
            A_sink = np.sqrt(Z_sink[0]**2 + Z_sink[1]**2)

            ax.annotate(f"|obs|={A_obs:.3f}", xy=Z_obs, fontsize=6.5,
                        textcoords="offset points", xytext=(5, 5), color="green")
            if A_src > 1e-4:
                ax.annotate(f"|src|={A_src:.3f}", xy=Z_src, fontsize=6.5,
                            textcoords="offset points", xytext=(5, -10), color="red")
            ax.annotate(f"|sink|={A_sink:.3f}", xy=Z_sink, fontsize=6.5,
                        textcoords="offset points", xytext=(5, 5), color="C0")

            # Axis setup — scale to data, NOT equal aspect
            ax.set_xlim(min(all_x) - pad_x, max(all_x) + pad_x)
            ax.set_ylim(min(all_y) - pad_y, max(all_y) + pad_y)
            ax.axhline(0, color="gray", lw=0.5, zorder=0)
            ax.axvline(0, color="gray", lw=0.5, zorder=0)
            ax.set_xlabel("B  (sin coeff, ‰)", fontsize=8)
            ax.set_ylabel("C  (cos coeff, ‰)", fontsize=8)
            ax.set_title(f"{code} — {label}\n"
                         f"R_obs={s['R_obs']:.4f} → R_corr={s['R_corrected']:.4f}",
                         fontsize=9)
            ax.grid(True, alpha=0.2)

            # Legend (top-left panel only)
            if j == 0 and row == 0:
                ax.plot([], [], "g-",  lw=2.5, label="Z_obs (observed)")
                ax.plot([], [], "r--", lw=2,   label="Z_src (wetland source)")
                ax.plot([], [], "b-",  lw=2.5, label="Z_sink (corrected)")
                ax.legend(fontsize=7, loc="best")

    fig.suptitle("Fig 8: Phasor decomposition — vector subtraction of wetland source",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = FIG_DIR / "fig8_phasor_decomposition.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"✓ Figure saved: {out}")


# ============================================================================
# FIGURE 9: CORRECTED RATIO VS LATITUDE
# ============================================================================
def plot_corrected_ratio(results, fits):
    """Fig 9: R_corrected vs R_observed vs latitude."""
    fig, ax = plt.subplots(figsize=(10, 5.5))

    codes_all = sorted(results["sites"].keys(),
                       key=lambda c: -results["sites"][c].get("peak_obs_13c", 0))
    # Sort by latitude
    codes_all = sorted(codes_all, key=lambda c: fits[c]["d13C"]["B"])  # proxy; use actual lat
    # Actually use the site_fit lat info - approximate from known values
    site_lats = {"ALT": 82.5, "ZEP": 78.9, "BRW": 71.3, "CBA": 55.2, "MHD": 53.3,
                 "AZR": 38.8, "MLO": 19.5, "KUM": 19.6, "ASC": -8.0, "SMO": -14.2,
                 "CGO": -40.7, "SPO": -90.0}
    codes_all = sorted(results["sites"].keys(), key=lambda c: -site_lats.get(c, 0))

    lats = [site_lats.get(c, 0) for c in codes_all]
    R_obs  = [results["sites"][c]["R_obs"] for c in codes_all]
    R_corr = [results["sites"][c]["R_corrected"] for c in codes_all]

    # MC CIs
    R_corr_lo = [results["sites"][c]["mc"]["R_corrected_ci95"][0] for c in codes_all]
    R_corr_hi = [results["sites"][c]["mc"]["R_corrected_ci95"][1] for c in codes_all]
    R_err_lo = [r - lo for r, lo in zip(R_corr, R_corr_lo)]
    R_err_hi = [hi - r for r, hi in zip(R_corr, R_corr_hi)]

    ax.scatter(lats, R_obs, s=80, c="gray", marker="o", alpha=0.5, label="R observed (uncorrected)", zorder=3)
    ax.errorbar(lats, R_corr, yerr=[R_err_lo, R_err_hi],
                fmt="s", color="C0", ms=8, capsize=4, lw=1.5,
                label="R corrected (phasor)", zorder=4)

    # OH-only predictions (inline; avoids import from phase5_kie)
    def _bulk_eps(alpha_13c_oh):
        e13 = (F_OH*(alpha_13c_oh-1) + F_CL*(ALPHA_13C_CL-1)
               + F_SOIL*(ALPHA_13C_SOIL-1) + F_STRAT*(ALPHA_13C_STRAT-1)) * 1000
        eD  = (F_OH*(ALPHA_D_OH-1) + F_CL*(ALPHA_D_CL-1)
               + F_SOIL*(ALPHA_D_SOIL-1) + F_STRAT*(ALPHA_D_STRAT-1)) * 1000
        return e13, eD
    eps13_s, epsD_s = _bulk_eps(ALPHA_13C_SAUERESSIG)
    eps13_c, epsD_c = _bulk_eps(ALPHA_13C_CANTRELL)
    R_saueressig = eps13_s / epsD_s
    R_cantrell   = eps13_c / epsD_c

    ax.axhline(R_saueressig, color="C2", ls="--", lw=1.5, label=f"Saueressig (R={R_saueressig:.4f})")
    ax.axhline(R_cantrell,   color="C3", ls="--", lw=1.5, label=f"Cantrell (R={R_cantrell:.4f})")

    for c, lat, ro, rc in zip(codes_all, lats, R_obs, R_corr):
        ax.annotate(c, (lat, rc), textcoords="offset points", xytext=(6, 6), fontsize=8)
        # Arrow from observed to corrected
        ax.annotate("", xy=(lat, rc), xytext=(lat, ro),
                    arrowprops=dict(arrowstyle="->", color="gray", lw=0.8, alpha=0.5))

    ax.set_xlabel("Latitude (°)", fontsize=11)
    ax.set_ylabel("Amplitude ratio R = A(δ¹³C) / A(δD)", fontsize=11)
    ax.set_title("Fig 9: Phasor-corrected amplitude ratio vs latitude", fontsize=12)
    ax.legend(fontsize=8, loc="upper right")
    ax.set_ylim(bottom=-0.01)

    fig.tight_layout()
    out = FIG_DIR / "fig9_corrected_ratio.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"✓ Figure saved: {out}")


# ============================================================================
# FIGURE 10: FINAL α_13C_OH CONSTRAINT
# ============================================================================
def plot_alpha_constraint(results):
    """Fig 10: Final α_13C_OH constraint vs Saueressig / Cantrell."""
    fig, ax = plt.subplots(figsize=(10, 5))

    site_lats = {"ALT": 82.5, "ZEP": 78.9, "BRW": 71.3, "CBA": 55.2, "MHD": 53.3,
                 "KUM": 19.6, "CGO": -40.7, "SPO": -90.0}

    codes = sorted(results["sites"].keys(), key=lambda c: -site_lats.get(c, 0))
    y_pos = np.arange(len(codes))

    for i, code in enumerate(codes):
        s = results["sites"][code]["mc"]
        med = s["alpha_13c_oh_median"]
        lo, hi = s["alpha_13c_oh_ci95"]
        ax.errorbar(med, i, xerr=[[med - lo], [hi - med]], fmt="o", color="C0",
                    ms=8, capsize=4, lw=2)

    # Multi-site
    ms = results["multi_site_result"]
    if ms:
        y_multi = len(codes) + 1
        ax.errorbar(ms["alpha_13c_oh_median"], y_multi,
                    xerr=[[ms["alpha_13c_oh_median"] - ms["alpha_13c_oh_ci95"][0]],
                          [ms["alpha_13c_oh_ci95"][1] - ms["alpha_13c_oh_median"]]],
                    fmt="D", color="black", ms=10, capsize=5, lw=2.5)
        codes_labels = codes + ["", "ALL SITES\n(phasor-corrected)"]

        if "sh_only" in ms:
            y_sh = len(codes) + 2
            sh = ms["sh_only"]
            ax.errorbar(sh["alpha_13c_oh_median"], y_sh,
                        xerr=[[sh["alpha_13c_oh_median"] - sh["alpha_13c_oh_ci95"][0]],
                              [sh["alpha_13c_oh_ci95"][1] - sh["alpha_13c_oh_median"]]],
                        fmt="D", color="C1", ms=10, capsize=5, lw=2.5)
            codes_labels.append("SH ONLY\n(phasor-corrected)")
    else:
        codes_labels = list(codes)

    # Reference lines
    ax.axvline(ALPHA_13C_SAUERESSIG, color="C2", ls="--", lw=2,
               label=f"Saueressig (α = {ALPHA_13C_SAUERESSIG})")
    ax.axvline(ALPHA_13C_CANTRELL, color="C3", ls="--", lw=2,
               label=f"Cantrell (α = {ALPHA_13C_CANTRELL})")
    ax.axvspan(ALPHA_13C_SAUERESSIG, ALPHA_13C_CANTRELL, alpha=0.1, color="gray",
               label="Contested range")

    ax.set_yticks(range(len(codes_labels)))
    ax.set_yticklabels(codes_labels, fontsize=9)
    ax.set_xlabel("α¹³C_OH", fontsize=11)
    ax.set_title("Fig 10: OH ¹³C KIE constraint from phasor-corrected seasonal amplitudes",
                 fontsize=11)
    ax.legend(fontsize=8, loc="upper right")
    ax.set_xlim(0.995, 1.015)

    fig.tight_layout()
    out = FIG_DIR / "fig10_alpha_constraint.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"✓ Figure saved: {out}")


def plot_phasor_diagrams(results):
    """Fig 8: Polar phasor clocks for representative source correction sites."""
    show_sites = ["BRW", "CBA", "CGO", "SPO"]
    available = [s for s in show_sites if s in results["sites"]]
    if not available:
        return

    month_names = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

    def vector_amp_month(z):
        B, C = z
        return float(np.hypot(B, C)), phasor_peak_month(B, C)

    def setup_clock(ax, rmax, title=None, show_rlabels=True):
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_xticks(np.linspace(0, 2 * np.pi, 12, endpoint=False))
        ax.set_xticklabels(month_names, fontsize=7)
        ax.set_ylim(0, rmax)
        ax.set_rlabel_position(225)
        if show_rlabels:
            ax.tick_params(axis="y", labelsize=6, pad=0)
        else:
            ax.set_yticklabels([])
        ax.grid(True, alpha=0.22, lw=0.7)
        if title:
            ax.set_title(title, fontsize=9.5, fontweight="bold", pad=16)

    def draw_arrow(ax, month, amp, color, linestyle, linewidth, alpha=1.0):
        theta = 2 * np.pi * month / 12.0
        ax.annotate(
            "",
            xy=(theta, amp),
            xytext=(theta, 0),
            arrowprops=dict(
                arrowstyle="-|>",
                color=color,
                lw=linewidth,
                linestyle=linestyle,
                mutation_scale=13,
                alpha=alpha,
                shrinkA=0,
                shrinkB=0,
            ),
        )
        ax.plot(theta, amp, "o", ms=4.5, color=color, alpha=alpha, zorder=5)

    row_max = {}
    for iso in ["13c", "dD"]:
        amps = []
        for code in available:
            site = results["sites"][code]
            for kind in ["obs", "src", "sink"]:
                amp, _ = vector_amp_month(site[f"Z_{kind}_{iso}"])
                amps.append(amp)
        row_max[iso] = max(amps) * 1.18

    fig, axes = plt.subplots(
        2, len(available), figsize=(4.35 * len(available), 8.2),
        subplot_kw={"projection": "polar"},
    )
    if len(available) == 1:
        axes = axes.reshape(2, 1)

    colors = {"obs": "#2ca02c", "src": "#d62728", "sink": "#1f77b4"}
    styles = {"obs": "-", "src": "--", "sink": "-"}
    widths = {"obs": 2.2, "src": 2.0, "sink": 2.6}

    for j, code in enumerate(available):
        site = results["sites"][code]
        for row, iso, iso_label in [(0, "13c", "d13C"), (1, "dD", "dD")]:
            ax = axes[row, j]
            vectors = {}
            for kind in ["obs", "src", "sink"]:
                amp, month = vector_amp_month(site[f"Z_{kind}_{iso}"])
                vectors[kind] = {"amp": amp, "month": month}

            title = (
                f"{code} - {iso_label}\n"
                f"R {site['R_obs']:.3f} -> {site['R_corrected']:.3f}"
                if row == 0 else iso_label
            )
            setup_clock(ax, row_max[iso], title=title, show_rlabels=(j == 0))

            for kind in ["obs", "src", "sink"]:
                draw_arrow(
                    ax,
                    vectors[kind]["month"],
                    vectors[kind]["amp"],
                    colors[kind],
                    styles[kind],
                    widths[kind],
                    alpha=0.9 if kind == "src" else 1.0,
                )

            amp_text = "\n".join(
                f"{kind} {vectors[kind]['amp']:.2f}" for kind in ["obs", "src", "sink"]
            )
            ax.text(
                0.04, 0.04, amp_text,
                transform=ax.transAxes,
                fontsize=7,
                ha="left",
                va="bottom",
                bbox=dict(facecolor="white", edgecolor="0.85", alpha=0.86, pad=2),
            )

            if code in ("CGO", "SPO"):
                inset = ax.inset_axes([0.58, 0.02, 0.40, 0.40], projection="polar")
                local_max = max(v["amp"] for v in vectors.values()) * 1.18
                setup_clock(inset, local_max, show_rlabels=False)
                inset.set_xticklabels([])
                inset.set_title("zoom", fontsize=6.5, pad=1)
                for kind in ["obs", "src", "sink"]:
                    draw_arrow(
                        inset,
                        vectors[kind]["month"],
                        vectors[kind]["amp"],
                        colors[kind],
                        styles[kind],
                        1.4,
                        alpha=0.9 if kind == "src" else 1.0,
                    )

    legend_elements = [
        plt.Line2D([0], [0], color=colors["obs"], lw=2.2, label="Z_obs observed"),
        plt.Line2D([0], [0], color=colors["src"], lw=2.0, ls="--", label="Z_src wetland"),
        plt.Line2D([0], [0], color=colors["sink"], lw=2.6, label="Z_sink corrected"),
    ]
    fig.legend(
        handles=legend_elements,
        loc="lower center",
        ncol=3,
        fontsize=9,
        framealpha=0.92,
        bbox_to_anchor=(0.5, 0.015),
    )

    fig.suptitle(
        "Fig 8: Phasor clocks - wetland source subtraction in amplitude/phase space",
        fontsize=13,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0.06, 1, 0.93])
    out = FIG_DIR / "fig8_phasor_decomposition.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"鉁?Figure saved: {out}")


if __name__ == "__main__":
    main()
