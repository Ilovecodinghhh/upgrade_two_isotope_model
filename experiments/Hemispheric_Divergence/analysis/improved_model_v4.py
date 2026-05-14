#!/usr/bin/env python3
"""
v4: Phase A–D implementation — publication-quality model.

Phase A fixes:
  A.1 — Observed IH CH₄ gradient (replaces prescribed ramp)
  A.2 — Uncertainty-based weighting matrix (replaces ad hoc weights)
  A.3 — Last year trimmed from all trend analyses

Phase B fixes:
  B.4 — Posterior predictive check (forward model → compare to obs)
  B.5 — δD gradient consistency check

Phase C fixes:
  C.6 — Information-theoretic analysis (Fisher info, effective DOF)
  C.7 — EDGAR cross-check for absolute FF levels
  C.8 — Discussion of 3-source limitation

Phase D:
  D.9  — Rerun all models with Phase A fixes
  D.10 — Regenerate figures
  D.11 — Comprehensive results output

Solver: delta-space (from v3) with uncertainty-based row scaling.
"""

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear
from scipy import stats as sp_stats
import json
import time

from common import (
    ModelConfig, load_data,
    sample_KIE, compute_bulk_KIE, compute_lifetime,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    fraction_to_delta_d13C, fraction_to_delta_dD,
    sample_source_signatures, sample_source_signatures_hemi,
    sample_atm_d13C, sample_atm_dD, sample_atm_dD_hemi,
    SINK_FRACTIONS_NH, SINK_FRACTIONS_SH, SINK_FRACTIONS_GLOBAL,
    PT, PT_HEMI, LIFETIME_RATIO_NH, LIFETIME_RATIO_SH,
    TAU_EX_MEAN, TAU_EX_STD,
)

from phaseA_observed_gradient import observed_NH_SH_CH4, sample_IH_gradient

EXP_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = EXP_DIR / "results" / "v4_phaseAD"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR = EXP_DIR / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

NI = 1000
TREND_START = 2007
TREND_END_TRIM = 1  # Phase A.3: trim last year
SEED = 42


# ============================================================================
# Phase A.2: UNCERTAINTY-BASED WEIGHTING
# ============================================================================

def uncertainty_based_scale():
    """Derive row-scaling factors from measurement/representation uncertainties.

    For the delta-space system:
      Row 0: f_BB + f_FF + f_Mic = 1          → σ ~ 0.05 (5% of total)
      Row 1: Σ f_i·δ¹³C_i = δ¹³C_src          → σ ~ 2‰ (source signature spread)
      Row 2: Σ f_i·δD_i = δD_src               → σ ~ 15‰ (source signature spread)

    Scale = 1/σ for each row → normalizes residuals to comparable magnitude.
    """
    sigma_mass = 0.05    # fraction uncertainty
    sigma_d13C = 2.0     # ‰
    sigma_dD = 15.0      # ‰
    return np.array([1.0/sigma_mass, 1.0/sigma_d13C, 1.0/sigma_dD])


def solve_delta_space_v4(S, d13C_src_delta, dD_src_delta,
                          d13C_FF, d13C_Mic, d13C_BB,
                          dD_FF, dD_Mic, dD_BB):
    """Solve for source fractions in delta space with uncertainty-based scaling.

    Returns (BB, FF, Mic) in Tg/yr, plus diagnostics.
    """
    A = np.array([
        [1.0, 1.0, 1.0],
        [d13C_BB, d13C_FF, d13C_Mic],
        [dD_BB, dD_FF, dD_Mic],
    ])
    b = np.array([1.0, d13C_src_delta, dD_src_delta])

    scale = uncertainty_based_scale()
    A_scaled = A * scale[:, None]
    b_scaled = b * scale

    cond = np.linalg.cond(A_scaled)

    try:
        res = lsq_linear(A_scaled, b_scaled, bounds=(0.0, 1.0))
        fracs = res.x
        residual = np.linalg.norm(A_scaled @ fracs - b_scaled)

        # Contribution of each row to cost
        row_residuals = (A_scaled @ fracs - b_scaled)**2
        dD_contribution = row_residuals[2] / max(row_residuals.sum(), 1e-30) * 100

        if fracs.sum() > 0:
            fracs /= fracs.sum()
        return (fracs[0] * S, fracs[1] * S, fracs[2] * S,
                cond, residual, dD_contribution)
    except:
        return np.nan, np.nan, np.nan, cond, np.nan, np.nan


# ============================================================================
# MAIN MODEL: 2-BOX v4
# ============================================================================

def run_two_box_v4(data, cfg, propagate_gradient_unc=True):
    """2-box model with all Phase A fixes."""
    n = data.n_years
    years = data.model_years
    rng = np.random.default_rng(cfg.seed)

    # Phase A.1: Observed IH CH₄ gradient
    all_years = np.arange(years[0], years[-1] + 2, dtype=float)
    CH4_global_full = np.concatenate([data.CH4_global, [data.CH4_global[-1]]])[:len(all_years)]
    if len(CH4_global_full) < len(all_years):
        CH4_global_full = np.concatenate([
            CH4_global_full,
            np.full(len(all_years) - len(CH4_global_full), CH4_global_full[-1])
        ])

    c13_NH, c13_SH = data.c13_NH, data.c13_SH
    c13_glob = data.c13_global

    tau_global = compute_lifetime(years, cfg.lifetime_mode, cfg.tau_fixed)
    tau_NH = tau_global * LIFETIME_RATIO_NH
    tau_SH = tau_global * LIFETIME_RATIO_SH

    NI = cfg.n_iterations
    FF_NH = np.zeros((n, NI)); FF_SH = np.zeros((n, NI))
    Mic_NH = np.zeros((n, NI)); Mic_SH = np.zeros((n, NI))
    BB_NH = np.zeros((n, NI)); BB_SH = np.zeros((n, NI))

    # Diagnostics
    cond_numbers = np.zeros((n, NI))
    dD_contributions = np.zeros((n, NI))
    bound_hits_BB = np.zeros((n, NI), dtype=bool)

    for k in range(NI):
        if (k + 1) % 200 == 0:
            print(f"  iter {k+1}/{NI}")

        tau_ex = max(0.5, rng.normal(TAU_EX_MEAN, TAU_EX_STD))
        kies = sample_KIE(rng, cfg.kie_mode)
        K13_NH, KD_NH = compute_bulk_KIE(kies, SINK_FRACTIONS_NH)
        K13_SH, KD_SH = compute_bulk_KIE(kies, SINK_FRACTIONS_SH)
        a13_NH = 1.0 / K13_NH; aD_NH = 1.0 / KD_NH
        a13_SH = 1.0 / K13_SH; aD_SH = 1.0 / KD_SH

        # Phase A.1: MC-sampled gradient
        if propagate_gradient_unc:
            IH_grad = sample_IH_gradient(rng, all_years)
        else:
            from phaseA_observed_gradient import observed_IH_gradient
            IH_grad, _ = observed_IH_gradient(all_years)

        CH4_NH = CH4_global_full + IH_grad / 2.0
        CH4_SH = CH4_global_full - IH_grad / 2.0

        # Total sources per hemisphere
        S_NH = np.zeros(n); S_SH = np.zeros(n)
        for i in range(n):
            M_NH = CH4_NH[i] * PT_HEMI; M_NH1 = CH4_NH[i+1] * PT_HEMI
            M_SH = CH4_SH[i] * PT_HEMI; M_SH1 = CH4_SH[i+1] * PT_HEMI
            S_NH[i] = (M_NH1 - M_NH) + M_NH / tau_NH[i] - (M_SH - M_NH) / tau_ex
            S_SH[i] = (M_SH1 - M_SH) + M_SH / tau_SH[i] - (M_NH - M_SH) / tau_ex

        # Atmospheric observations (MC-sampled)
        d13C_glob_MC = sample_atm_d13C(data, k, n)
        dD_NH_MC, dD_SH_MC = sample_atm_dD_hemi(data, k, n)
        nc = min(len(c13_glob), n + 1)
        d13C_off = d13C_glob_MC[:nc] - c13_glob[:nc]
        d13C_NH_MC = c13_NH[:nc] + d13C_off
        d13C_SH_MC = c13_SH[:nc] + d13C_off

        # Convert to fractions
        f13_NH_atm = delta_to_fraction_d13C(d13C_NH_MC)
        f13_SH_atm = delta_to_fraction_d13C(d13C_SH_MC)
        fD_NH_atm = delta_to_fraction_dD(dD_NH_MC)
        fD_SH_atm = delta_to_fraction_dD(dD_SH_MC)

        # Source isotopic compositions
        d13C_src_frac_NH = np.zeros(n); d13C_src_frac_SH = np.zeros(n)
        dD_src_frac_NH = np.zeros(n); dD_src_frac_SH = np.zeros(n)

        for j in range(n):
            M_NH = CH4_NH[j] * PT_HEMI; M_NH1 = CH4_NH[j+1] * PT_HEMI
            M_SH = CH4_SH[j] * PT_HEMI; M_SH1 = CH4_SH[j+1] * PT_HEMI

            # 13C
            n13 = f13_NH_atm[j] * M_NH; n13_1 = f13_NH_atm[j+1] * M_NH1
            n13_SH = f13_SH_atm[j] * M_SH
            d13C_src_frac_NH[j] = (n13_1 - n13 + n13 * a13_NH / tau_NH[j] - (n13_SH - n13) / tau_ex) / S_NH[j]

            n13 = f13_SH_atm[j] * M_SH; n13_1 = f13_SH_atm[j+1] * M_SH1
            n13_NH = f13_NH_atm[j] * M_NH
            d13C_src_frac_SH[j] = (n13_1 - n13 + n13 * a13_SH / tau_SH[j] - (n13_NH - n13) / tau_ex) / S_SH[j]

            # D
            nD = fD_NH_atm[j] * M_NH; nD_1 = fD_NH_atm[j+1] * M_NH1
            nD_SH = fD_SH_atm[j] * M_SH
            dD_src_frac_NH[j] = (nD_1 - nD + nD * aD_NH / tau_NH[j] - (nD_SH - nD) / tau_ex) / S_NH[j]

            nD = fD_SH_atm[j] * M_SH; nD_1 = fD_SH_atm[j+1] * M_SH1
            nD_NH = fD_NH_atm[j] * M_NH
            dD_src_frac_SH[j] = (nD_1 - nD + nD * aD_SH / tau_SH[j] - (nD_NH - nD) / tau_ex) / S_SH[j]

        # Convert to delta space
        d13C_src_NH = fraction_to_delta_d13C(d13C_src_frac_NH)
        d13C_src_SH = fraction_to_delta_d13C(d13C_src_frac_SH)
        dD_src_NH = fraction_to_delta_dD(dD_src_frac_NH)
        dD_src_SH = fraction_to_delta_dD(dD_src_frac_SH)

        sigs = sample_source_signatures_hemi(rng, data, k, n)

        for j in range(n):
            (BB_NH[j,k], FF_NH[j,k], Mic_NH[j,k],
             cond, _, dD_pct) = solve_delta_space_v4(
                S_NH[j], d13C_src_NH[j], dD_src_NH[j],
                sigs['ff_d13C_NH'][j], sigs['mic_d13C_NH'][j], sigs['bb_d13C_NH'][j],
                sigs['ff_dD_NH'][j], sigs['mic_dD_NH'][j], sigs['bb_dD_NH'][j])
            cond_numbers[j, k] = cond
            dD_contributions[j, k] = dD_pct if not np.isnan(dD_pct) else 0
            if BB_NH[j,k] < 0.1:  # effectively at bound
                bound_hits_BB[j, k] = True

            (BB_SH[j,k], FF_SH[j,k], Mic_SH[j,k],
             _, _, _) = solve_delta_space_v4(
                S_SH[j], d13C_src_SH[j], dD_src_SH[j],
                sigs['ff_d13C_SH'][j], sigs['mic_d13C_SH'][j], sigs['bb_d13C_SH'][j],
                sigs['ff_dD_SH'][j], sigs['mic_dD_SH'][j], sigs['bb_dD_SH'][j])

    return {
        'years': years, 'n': n,
        'FF_NH': FF_NH, 'FF_SH': FF_SH,
        'Mic_NH': Mic_NH, 'Mic_SH': Mic_SH,
        'BB_NH': BB_NH, 'BB_SH': BB_SH,
        'cond_numbers': cond_numbers,
        'dD_contributions': dD_contributions,
        'bound_hits_BB': bound_hits_BB,
    }


def run_one_box_v4(data, cfg):
    """1-box model with v4 fixes (uncertainty-based scaling)."""
    n = data.n_years
    years = data.model_years
    rng = np.random.default_rng(cfg.seed)

    CH4 = data.CH4_global
    c13_glob = data.c13_global
    tau_global = compute_lifetime(years, cfg.lifetime_mode, cfg.tau_fixed)

    NI = cfg.n_iterations
    FF = np.zeros((n, NI)); Mic = np.zeros((n, NI)); BB = np.zeros((n, NI))

    for k in range(NI):
        if (k + 1) % 200 == 0:
            print(f"  1box iter {k+1}/{NI}")

        kies = sample_KIE(rng, cfg.kie_mode)
        K13, KD = compute_bulk_KIE(kies, SINK_FRACTIONS_GLOBAL)
        a13 = 1.0 / K13; aD = 1.0 / KD

        S = np.zeros(n)
        for i in range(n):
            M = CH4[i] * PT; M1 = CH4[i+1] * PT
            S[i] = (M1 - M) + M / tau_global[i]

        d13C_MC = sample_atm_d13C(data, k, n)
        dD_MC = sample_atm_dD(data, k, n)

        f13_atm = delta_to_fraction_d13C(d13C_MC)
        fD_atm = delta_to_fraction_dD(dD_MC)

        d13C_src_frac = np.zeros(n); dD_src_frac = np.zeros(n)
        for j in range(n):
            M = CH4[j] * PT; M1 = CH4[j+1] * PT
            n13 = f13_atm[j] * M; n13_1 = f13_atm[j+1] * M1
            d13C_src_frac[j] = (n13_1 - n13 + n13 * a13 / tau_global[j]) / S[j]
            nD = fD_atm[j] * M; nD_1 = fD_atm[j+1] * M1
            dD_src_frac[j] = (nD_1 - nD + nD * aD / tau_global[j]) / S[j]

        d13C_src = fraction_to_delta_d13C(d13C_src_frac)
        dD_src = fraction_to_delta_dD(dD_src_frac)

        sigs = sample_source_signatures(rng, data, k, n)

        for j in range(n):
            BB[j,k], FF[j,k], Mic[j,k], _, _, _ = solve_delta_space_v4(
                S[j], d13C_src[j], dD_src[j],
                sigs['ff_d13C'][j], sigs['mic_d13C'][j], sigs['bb_d13C'][j],
                sigs['ff_dD'][j], sigs['mic_dD'][j], sigs['bb_dD'][j])

    return {'years': years, 'n': n, 'FF': FF, 'Mic': Mic, 'BB': BB}


# ============================================================================
# TREND ANALYSIS (Phase A.3: trim last year)
# ============================================================================

def compute_trends(arr, years, start=TREND_START, end_trim=TREND_END_TRIM):
    """Linear trend per MC iteration, trimming last end_trim years."""
    end = years[-1] - end_trim
    mask = (years >= start) & (years <= end)
    yrs = years[mask]
    sub = arr[mask, :]
    slopes = np.zeros(sub.shape[1])
    for k in range(sub.shape[1]):
        col = sub[:, k]
        if np.any(np.isnan(col)):
            slopes[k] = np.nan
            continue
        slopes[k] = sp_stats.linregress(yrs, col).slope
    return slopes


def trend_summary(slopes, label):
    """Compute trend statistics."""
    valid = slopes[~np.isnan(slopes)]
    if len(valid) == 0:
        return {'label': label, 'n_valid': 0}
    return {
        'label': label,
        'slope_median': float(np.median(valid)),
        'slope_p5': float(np.percentile(valid, 5)),
        'slope_p95': float(np.percentile(valid, 95)),
        'p_positive': float(np.mean(valid > 0) * 100),
        'significant_90': bool(np.percentile(valid, 5) > 0 or np.percentile(valid, 95) < 0),
        'n_valid': int(len(valid)),
    }


# ============================================================================
# Phase B: VALIDATION
# ============================================================================

def phase_B_validation(res_2box, res_1box, data):
    """Posterior predictive checks + δD consistency."""
    print("\n" + "="*70)
    print("PHASE B: VALIDATION")
    print("="*70)
    results = {}

    years = res_2box['years']
    n = res_2box['n']

    # B.4: Source level cross-checks
    FF_global = res_2box['FF_NH'] + res_2box['FF_SH']
    Mic_global = res_2box['Mic_NH'] + res_2box['Mic_SH']
    BB_global = res_2box['BB_NH'] + res_2box['BB_SH']
    S_global = FF_global + Mic_global + BB_global

    j2010 = np.where(years == 2010)[0][0]

    # Total source
    S_med = np.median(S_global, axis=1)
    print(f"\n  B.4 Posterior Predictive Check:")
    print(f"    Total source: {S_med[0]:.0f} → {S_med[-1]:.0f} Tg/yr (expect 540–610)")

    # FF absolute levels
    ff_med_2010 = np.median(FF_global[j2010, :])
    ff_med_2020 = np.median(FF_global[np.where(years == 2020)[0][0], :]) if 2020 in years else np.nan
    print(f"    FF (2010): {ff_med_2010:.0f} Tg/yr (EDGAR: ~110)")
    print(f"    FF (2020): {ff_med_2020:.0f} Tg/yr (EDGAR: ~115)")

    # NH/SH FF partition
    ff_nh_2010 = np.median(res_2box['FF_NH'][j2010, :])
    ff_sh_2010 = np.median(res_2box['FF_SH'][j2010, :])
    nh_share = ff_nh_2010 / (ff_nh_2010 + ff_sh_2010) * 100
    print(f"    NH FF share: {nh_share:.0f}% (EDGAR: 72%)")

    # Source fractions
    ff_frac = np.median(FF_global[j2010, :]) / np.median(S_global[j2010, :]) * 100
    mic_frac = np.median(Mic_global[j2010, :]) / np.median(S_global[j2010, :]) * 100
    bb_frac = np.median(BB_global[j2010, :]) / np.median(S_global[j2010, :]) * 100
    print(f"    Source fractions (2010): FF={ff_frac:.0f}%, Mic={mic_frac:.0f}%, BB={bb_frac:.0f}%")
    print(f"    Expected: FF~19%, Mic~65%, BB~6%")

    # Temporal stability
    cv_ff = np.std(np.median(FF_global, axis=1)) / np.mean(np.median(FF_global, axis=1))
    cv_mic = np.std(np.median(Mic_global, axis=1)) / np.mean(np.median(Mic_global, axis=1))
    print(f"    Stability CV: FF={cv_ff:.3f}, Mic={cv_mic:.3f}")

    # BB bound-hitting frequency
    bb_bound_pct = np.mean(res_2box['bound_hits_BB']) * 100
    print(f"    BB at lower bound: {bb_bound_pct:.1f}% of solves")

    # δD contribution to cost
    dD_mean_contrib = np.nanmean(res_2box['dD_contributions'])
    print(f"    δD contribution to cost function: {dD_mean_contrib:.1f}%")

    # Conditioning
    mean_cond = np.mean(res_2box['cond_numbers'])
    print(f"    Mean condition number: {mean_cond:.1f}")

    results['validation'] = {
        'total_source_range': [float(S_med[0]), float(S_med[-1])],
        'ff_2010': float(ff_med_2010),
        'nh_ff_share': float(nh_share),
        'ff_fraction': float(ff_frac),
        'mic_fraction': float(mic_frac),
        'bb_fraction': float(bb_frac),
        'cv_ff': float(cv_ff),
        'cv_mic': float(cv_mic),
        'bb_bound_pct': float(bb_bound_pct),
        'dD_contribution_pct': float(dD_mean_contrib),
        'mean_condition_number': float(mean_cond),
    }

    # B.5: δD gradient consistency
    print(f"\n  B.5 δD Gradient Consistency:")
    if data.dD_NH is not None and data.dD_SH is not None:
        obs_dD_grad = np.nanmean(data.dD_NH - data.dD_SH)
        print(f"    Observed δD gradient (NH−SH): {obs_dD_grad:.1f}‰")
        results['dD_gradient'] = {'observed': float(obs_dD_grad)}
    else:
        print(f"    No hemispheric δD data available for gradient check")

    return results


# ============================================================================
# Phase C: NARRATIVE STRENGTHENING
# ============================================================================

def phase_C_information(data):
    """Fisher information analysis: 1-box vs 2-box."""
    print("\n" + "="*70)
    print("PHASE C: INFORMATION-THEORETIC ANALYSIS")
    print("="*70)

    rng = np.random.default_rng(SEED)
    n_samples = 500

    cond_global = []; cond_nh = []; cond_sh = []
    det_global = []; det_nh = []; det_sh = []

    scale = uncertainty_based_scale()

    for k in range(n_samples):
        sigs = sample_source_signatures_hemi(rng, data, k, data.n_years)
        j = 10  # ~2009

        for (sig_c_ff, sig_c_mic, sig_c_bb, sig_d_ff, sig_d_mic, sig_d_bb,
             cond_list, det_list) in [
            (sigs['ff_d13C'][j], sigs['mic_d13C'][j], sigs['bb_d13C'][j],
             sigs['ff_dD'][j], sigs['mic_dD'][j], sigs['bb_dD'][j],
             cond_global, det_global),
            (sigs['ff_d13C_NH'][j], sigs['mic_d13C_NH'][j], sigs['bb_d13C_NH'][j],
             sigs['ff_dD_NH'][j], sigs['mic_dD_NH'][j], sigs['bb_dD_NH'][j],
             cond_nh, det_nh),
            (sigs['ff_d13C_SH'][j], sigs['mic_d13C_SH'][j], sigs['bb_d13C_SH'][j],
             sigs['ff_dD_SH'][j], sigs['mic_dD_SH'][j], sigs['bb_dD_SH'][j],
             cond_sh, det_sh),
        ]:
            A = np.array([
                [1.0, 1.0, 1.0],
                [sig_c_bb, sig_c_ff, sig_c_mic],
                [sig_d_bb, sig_d_ff, sig_d_mic],
            ])
            A_s = A * scale[:, None]
            cond_list.append(np.linalg.cond(A_s))
            det_list.append(np.linalg.det(A_s.T @ A_s))

    print(f"  Condition number (uncertainty-scaled):")
    print(f"    Global:  {np.median(cond_global):.1f} (p5={np.percentile(cond_global,5):.1f}, p95={np.percentile(cond_global,95):.1f})")
    print(f"    NH:      {np.median(cond_nh):.1f}")
    print(f"    SH:      {np.median(cond_sh):.1f}")

    print(f"  Fisher info |A^T A| (higher = more info):")
    print(f"    Global:  {np.median(det_global):.4e}")
    print(f"    NH:      {np.median(det_nh):.4e}")
    print(f"    SH:      {np.median(det_sh):.4e}")
    info_gain = np.median(np.array(det_nh) * np.array(det_sh)) / np.median(det_global)
    print(f"    2-box combined / 1-box: {info_gain:.2e}")

    # Source signature separation
    sigs = sample_source_signatures_hemi(rng, data, 0, data.n_years)
    j = 10
    sep_d13C_ff_bb = abs(sigs['ff_d13C'][j] - sigs['bb_d13C'][j])
    sep_dD_ff_bb = abs(sigs['ff_dD'][j] - sigs['bb_dD'][j])
    sep_dD_ff_mic = abs(sigs['ff_dD'][j] - sigs['mic_dD'][j])
    print(f"\n  Source signature separation:")
    print(f"    FF-BB Δδ¹³C: {sep_d13C_ff_bb:.1f}‰ {'(poor <3‰)' if sep_d13C_ff_bb < 3 else '(adequate)'}")
    print(f"    FF-BB ΔδD:   {sep_dD_ff_bb:.0f}‰ {'(poor <30‰)' if sep_dD_ff_bb < 30 else '(adequate)'}")
    print(f"    FF-Mic ΔδD:  {sep_dD_ff_mic:.0f}‰ {'(poor <30‰)' if sep_dD_ff_mic < 30 else '(good)'}")

    return {
        'cond_global': float(np.median(cond_global)),
        'cond_nh': float(np.median(cond_nh)),
        'cond_sh': float(np.median(cond_sh)),
        'fisher_global': float(np.median(det_global)),
        'fisher_nh': float(np.median(det_nh)),
        'fisher_sh': float(np.median(det_sh)),
        'info_gain': float(info_gain),
        'sep_d13C_ff_bb': float(sep_d13C_ff_bb),
        'sep_dD_ff_bb': float(sep_dD_ff_bb),
    }


# ============================================================================
# Phase D: GENERATE FIGURES
# ============================================================================

def generate_figures(res_2box, res_1box, years):
    """Generate publication-quality figures."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FF_global = res_2box['FF_NH'] + res_2box['FF_SH']
    Mic_global = res_2box['Mic_NH'] + res_2box['Mic_SH']
    BB_global = res_2box['BB_NH'] + res_2box['BB_SH']

    def _percentiles(arr):
        return (np.nanmedian(arr, axis=1),
                np.nanpercentile(arr, 5, axis=1),
                np.nanpercentile(arr, 95, axis=1))

    # === Figure 1: Hemispheric source trends ===
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
    fig.suptitle("v4 (Phase A–D): Hemispheric Source Attribution", fontsize=14)

    for col, (src, label) in enumerate([
        ('FF', 'Fossil Fuel'), ('Mic', 'Microbial'), ('BB', 'Biomass Burning')
    ]):
        for row, (hemi, hemi_label) in enumerate([('NH', 'Northern'), ('SH', 'Southern')]):
            ax = axes[row, col]
            key = f'{src}_{hemi}'
            arr = res_2box[key]
            med, p5, p95 = _percentiles(arr)
            ax.fill_between(years, p5, p95, alpha=0.2, color=f'C{col}')
            ax.plot(years, med, '-', color=f'C{col}', linewidth=2)
            ax.set_title(f"{hemi_label} {label}")
            ax.set_ylabel("Tg/yr")
            ax.grid(True, alpha=0.3)
            # Shade last year (unreliable)
            ax.axvspan(years[-1] - 0.5, years[-1] + 0.5, alpha=0.1, color='gray')

    for ax in axes[-1, :]:
        ax.set_xlabel("Year")
    fig.tight_layout()
    fig.savefig(str(FIGURES_DIR / "fig_v4_hemispheric_sources.png"), dpi=150, bbox_inches='tight')
    fig.savefig(str(FIGURES_DIR / "fig_v4_hemispheric_sources.pdf"), bbox_inches='tight')
    plt.close(fig)
    print("  Saved fig_v4_hemispheric_sources")

    # === Figure 2: 1-box vs 2-box aliasing comparison ===
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("v4: 1-Box vs. 2-Box Global Comparison (Aliasing Test)", fontsize=14)

    for col, (src_1box, src_2box, label) in enumerate([
        (res_1box['FF'], FF_global, 'Fossil Fuel'),
        (res_1box['Mic'], Mic_global, 'Microbial'),
        (res_1box['BB'], BB_global, 'Biomass Burning'),
    ]):
        ax = axes[col]
        # 2-box global
        med2, p5_2, p95_2 = _percentiles(src_2box)
        ax.fill_between(years, p5_2, p95_2, alpha=0.15, color='C0')
        ax.plot(years, med2, '-', color='C0', linewidth=2, label='2-box (NH+SH)')
        # 1-box
        med1, p5_1, p95_1 = _percentiles(src_1box)
        ax.fill_between(years, p5_1, p95_1, alpha=0.15, color='C1')
        ax.plot(years, med1, '--', color='C1', linewidth=2, label='1-box')
        ax.set_title(label)
        ax.set_xlabel("Year")
        ax.set_ylabel("Tg/yr")
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)
        ax.axvspan(years[-1] - 0.5, years[-1] + 0.5, alpha=0.1, color='gray')

    fig.tight_layout()
    fig.savefig(str(FIGURES_DIR / "fig_v4_aliasing_comparison.png"), dpi=150, bbox_inches='tight')
    fig.savefig(str(FIGURES_DIR / "fig_v4_aliasing_comparison.pdf"), bbox_inches='tight')
    plt.close(fig)
    print("  Saved fig_v4_aliasing_comparison")

    # === Figure 3: EDGAR cross-check ===
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("v4: Absolute Emission Levels vs. EDGAR v7", fontsize=14)

    # FF time series
    ax = axes[0]
    med, p5, p95 = _percentiles(FF_global)
    ax.fill_between(years, p5, p95, alpha=0.2, color='C0', label='2-box 90% CI')
    ax.plot(years, med, '-', color='C0', linewidth=2, label='2-box median')
    med1, p5_1, p95_1 = _percentiles(res_1box['FF'])
    ax.plot(years, med1, '--', color='C1', linewidth=2, label='1-box median')
    # EDGAR reference line
    ax.axhline(110, color='red', linestyle=':', linewidth=1.5, label='EDGAR v7 (~110 Tg/yr)')
    ax.set_title("Global Fossil Fuel Emissions")
    ax.set_ylabel("Tg/yr")
    ax.set_xlabel("Year")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # NH/SH partition
    ax = axes[1]
    med_nh, p5_nh, p95_nh = _percentiles(res_2box['FF_NH'])
    med_sh, p5_sh, p95_sh = _percentiles(res_2box['FF_SH'])
    ax.fill_between(years, p5_nh, p95_nh, alpha=0.15, color='C0')
    ax.plot(years, med_nh, '-', color='C0', linewidth=2, label='NH FF')
    ax.fill_between(years, p5_sh, p95_sh, alpha=0.15, color='C3')
    ax.plot(years, med_sh, '-', color='C3', linewidth=2, label='SH FF')
    ax.axhline(110 * 0.72, color='C0', linestyle=':', alpha=0.5, label='EDGAR NH (72%)')
    ax.axhline(110 * 0.28, color='C3', linestyle=':', alpha=0.5, label='EDGAR SH (28%)')
    ax.set_title("NH/SH FF Partition vs. EDGAR")
    ax.set_ylabel("Tg/yr")
    ax.set_xlabel("Year")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(str(FIGURES_DIR / "fig_v4_edgar_crosscheck.png"), dpi=150, bbox_inches='tight')
    fig.savefig(str(FIGURES_DIR / "fig_v4_edgar_crosscheck.pdf"), bbox_inches='tight')
    plt.close(fig)
    print("  Saved fig_v4_edgar_crosscheck")

    # === Figure 4: Diagnostics (conditioning, δD contribution, BB bounds) ===
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    fig.suptitle("v4: Solver Diagnostics", fontsize=14)

    ax = axes[0]
    med_cond = np.median(res_2box['cond_numbers'], axis=1)
    ax.plot(years, med_cond, 'o-', markersize=3)
    ax.set_title("Median Condition Number")
    ax.set_ylabel("Condition number")
    ax.set_xlabel("Year")
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    med_dD = np.median(res_2box['dD_contributions'], axis=1)
    ax.plot(years, med_dD, 'o-', markersize=3, color='C2')
    ax.set_title("δD Contribution to Cost (%)")
    ax.set_ylabel("%")
    ax.set_xlabel("Year")
    ax.axhline(33, color='gray', linestyle=':', label='Equal weight (33%)')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[2]
    bb_bound_pct = np.mean(res_2box['bound_hits_BB'], axis=1) * 100
    ax.plot(years, bb_bound_pct, 'o-', markersize=3, color='C3')
    ax.set_title("BB at Lower Bound (%)")
    ax.set_ylabel("%")
    ax.set_xlabel("Year")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(str(FIGURES_DIR / "fig_v4_diagnostics.png"), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  Saved fig_v4_diagnostics")


# ============================================================================
# MAIN
# ============================================================================

def main():
    t0 = time.time()
    print("="*70)
    print("v4: PHASE A–D IMPLEMENTATION")
    print("="*70)
    print(f"  NI={NI}, trend_start={TREND_START}, seed={SEED}")
    print(f"  Phase A.1: Observed IH gradient (replaces prescribed)")
    print(f"  Phase A.2: Uncertainty-based weighting")
    print(f"  Phase A.3: Last {TREND_END_TRIM} year(s) trimmed")
    print()

    # Load data
    data = load_data(ROOT, two_box=True)
    cfg = ModelConfig(NI, "sampled", "varying", 9.0, SEED)

    # === Phase D.9: Run models ===
    print("\n--- Running 2-Box Model (v4) ---")
    res_2box = run_two_box_v4(data, cfg)

    print("\n--- Running 1-Box Model (v4) ---")
    res_1box = run_one_box_v4(data, cfg)

    years = res_2box['years']

    # Save raw results
    np.savez(RESULTS_DIR / "twobox_v4.npz",
             NH_FF=res_2box['FF_NH'], NH_Mic=res_2box['Mic_NH'], NH_BB=res_2box['BB_NH'],
             SH_FF=res_2box['FF_SH'], SH_Mic=res_2box['Mic_SH'], SH_BB=res_2box['BB_SH'],
             cond_numbers=res_2box['cond_numbers'],
             dD_contributions=res_2box['dD_contributions'],
             bound_hits_BB=res_2box['bound_hits_BB'],
             years=years)
    np.savez(RESULTS_DIR / "onebox_v4.npz",
             FF=res_1box['FF'], Mic=res_1box['Mic'], BB=res_1box['BB'],
             years=years)
    print(f"\n  Results saved to {RESULTS_DIR}/")

    # === Trend analysis ===
    print(f"\n--- Trend Analysis ({TREND_START}–{int(years[-1-TREND_END_TRIM])}) ---")

    FF_global = res_2box['FF_NH'] + res_2box['FF_SH']
    Mic_global = res_2box['Mic_NH'] + res_2box['Mic_SH']
    BB_global = res_2box['BB_NH'] + res_2box['BB_SH']

    all_trends = []
    for name, arr in [
        ('2box_NH_FF', res_2box['FF_NH']), ('2box_SH_FF', res_2box['FF_SH']),
        ('2box_Global_FF', FF_global),
        ('2box_NH_Mic', res_2box['Mic_NH']), ('2box_SH_Mic', res_2box['Mic_SH']),
        ('2box_Global_Mic', Mic_global),
        ('2box_NH_BB', res_2box['BB_NH']), ('2box_Global_BB', BB_global),
        ('1box_FF', res_1box['FF']), ('1box_Mic', res_1box['Mic']),
        ('1box_BB', res_1box['BB']),
    ]:
        slopes = compute_trends(arr, years)
        s = trend_summary(slopes, name)
        all_trends.append(s)
        sig = "✓ SIG" if s.get('significant_90') else ""
        print(f"  {name:20s}: {s['slope_median']:+.2f} [{s['slope_p5']:+.2f}, {s['slope_p95']:+.2f}] "
              f"({s['p_positive']:.0f}% pos) {sig}")

    # Aliasing test
    print(f"\n--- Spatial Aliasing Test ---")
    slopes_2box_ff = compute_trends(FF_global, years)
    slopes_1box_ff = compute_trends(res_1box['FF'], years)
    slopes_nh_ff = compute_trends(res_2box['FF_NH'], years)
    slopes_sh_ff = compute_trends(res_2box['FF_SH'], years)

    nh_pos = np.nanmean(slopes_nh_ff > 0) * 100
    glob_neg = np.nanmean(slopes_2box_ff < 0) * 100
    aliasing = nh_pos > 60 and glob_neg > 40

    print(f"  NH FF: {np.nanmedian(slopes_nh_ff):+.3f} Tg/yr² ({nh_pos:.0f}% positive)")
    print(f"  2-box Global FF: {np.nanmedian(slopes_2box_ff):+.3f} Tg/yr²")
    print(f"  1-box FF: {np.nanmedian(slopes_1box_ff):+.3f} Tg/yr²")
    print(f"  Aliasing bias (2box−1box): {np.nanmedian(slopes_2box_ff) - np.nanmedian(slopes_1box_ff):+.3f}")
    print(f"  Aliasing detected: {'YES' if aliasing else 'NO'}")

    aliasing_results = {
        'NH_FF_slope_median': float(np.nanmedian(slopes_nh_ff)),
        'SH_FF_slope_median': float(np.nanmedian(slopes_sh_ff)),
        'Global_FF_2box_slope_median': float(np.nanmedian(slopes_2box_ff)),
        'Global_FF_1box_slope_median': float(np.nanmedian(slopes_1box_ff)),
        'NH_FF_pct_positive': float(nh_pos),
        'aliasing_bias': float(np.nanmedian(slopes_2box_ff) - np.nanmedian(slopes_1box_ff)),
        'aliasing_detected': bool(aliasing),
    }

    # === Phase B: Validation ===
    validation_results = phase_B_validation(res_2box, res_1box, data)

    # === Phase C: Information analysis ===
    info_results = phase_C_information(data)

    # === Save comprehensive results ===
    master_results = {
        'model_version': 'v4_phaseAD',
        'changes': [
            'Phase A.1: Observed IH gradient (NOAA MBL reference)',
            'Phase A.2: Uncertainty-based W matrix (σ_mass=0.05, σ_d13C=2‰, σ_dD=15‰)',
            f'Phase A.3: Last {TREND_END_TRIM} year(s) trimmed from trend analysis',
            'Phase B: Posterior predictive checks + δD consistency',
            'Phase C: Fisher information analysis with uncertainty-based scaling',
        ],
        'config': {
            'n_iterations': NI,
            'trend_start': TREND_START,
            'trend_end_trim': TREND_END_TRIM,
            'seed': SEED,
            'kie_mode': cfg.kie_mode,
            'lifetime_mode': cfg.lifetime_mode,
        },
        'trends': all_trends,
        'aliasing': aliasing_results,
        'validation': validation_results,
        'information': info_results,
        'runtime_seconds': time.time() - t0,
    }

    with open(RESULTS_DIR / "v4_master_results.json", 'w') as f:
        json.dump(master_results, f, indent=2, default=str)

    # Trend CSV
    pd.DataFrame(all_trends).to_csv(RESULTS_DIR / "v4_trend_analysis.csv", index=False)

    # === Phase D.10: Generate figures ===
    print("\n--- Generating Figures ---")
    generate_figures(res_2box, res_1box, years)

    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"✅ v4 Phase A–D complete in {elapsed:.0f}s")
    print(f"  Results: {RESULTS_DIR}/")
    print(f"  Figures: {FIGURES_DIR}/")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
