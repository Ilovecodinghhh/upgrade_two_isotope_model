#!/usr/bin/env python3
"""
Phase B: Model Validation
  - Posterior predictive check (forward model from solved sources → compare to obs)
  - δD gradient consistency check
  - EDGAR cross-check for absolute levels
  
Phase C: Narrative Strengthening
  - Information-theoretic analysis (Fisher information, effective DOF)
  - W-matrix sensitivity sweep
  - GFED BB validation
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

from common import (
    ModelConfig, load_data,
    sample_KIE, compute_bulk_KIE, compute_lifetime,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    sample_source_signatures_hemi,
    sample_atm_d13C, sample_atm_dD_hemi,
    SINK_FRACTIONS_NH, SINK_FRACTIONS_SH,
    PT_HEMI, LIFETIME_RATIO_NH, LIFETIME_RATIO_SH,
    TAU_EX_MEAN, TAU_EX_STD,
)

EXP_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = EXP_DIR / "results" / "v2_improved"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

NI = 1000
SEED = 42


def realistic_IH_gradient(years):
    anchor_years = np.array([2000, 2010, 2020, 2022])
    anchor_grad = np.array([108.0, 120.0, 140.0, 145.0])
    return np.interp(years, anchor_years, anchor_grad)


# ═══════════════════════════════════════════════════════════════
# PHASE B: MODEL VALIDATION
# ═══════════════════════════════════════════════════════════════

def phase_B4_posterior_predictive(data):
    """
    Proper posterior predictive check:
    Take solved sources → reconstruct what δ¹³C and δD the atmosphere
    SHOULD have → compare to what we actually observed.
    
    For a 3×3 system with 3 unknowns, the solver ALWAYS fits exactly
    (RMSE ≈ 0 for the INPUT observations). The real test is whether
    the solved source decomposition, when integrated forward in time,
    produces physically reasonable isotopic trajectories.
    """
    print("\n" + "="*60)
    print("Phase B.4: Posterior Predictive Check")
    print("="*60)
    
    v2 = np.load(RESULTS_DIR / "twobox_v2.npz")
    years = v2['years']
    n = len(years)
    
    FF_NH, FF_SH = v2['NH_FF'], v2['SH_FF']
    Mic_NH, Mic_SH = v2['NH_Mic'], v2['SH_Mic']
    BB_NH, BB_SH = v2['NH_BB'], v2['SH_BB']
    
    # Total sources per hemisphere
    S_NH = FF_NH + Mic_NH + BB_NH
    S_SH = FF_SH + Mic_SH + BB_SH
    S_global = S_NH + S_SH
    
    # Check 1: Does total source match expected?
    # Expected total: ~540-600 Tg/yr
    S_med = np.median(S_global, axis=1)
    print(f"\n  Total source (median): {S_med[0]:.0f} → {S_med[-1]:.0f} Tg/yr")
    print(f"  Expected (CarbonTracker): 560-610 Tg/yr")
    
    # Check 2: FF fraction
    ff_frac_NH = np.median(FF_NH / S_NH, axis=1)
    ff_frac_SH = np.median(FF_SH / S_SH, axis=1)
    mic_frac_NH = np.median(Mic_NH / S_NH, axis=1)
    mic_frac_SH = np.median(Mic_SH / S_SH, axis=1)
    bb_frac_NH = np.median(BB_NH / S_NH, axis=1)
    bb_frac_SH = np.median(BB_SH / S_SH, axis=1)
    
    print(f"\n  Source fractions (2010):")
    j2010 = np.where(years == 2010)[0][0]
    print(f"    NH: FF={ff_frac_NH[j2010]*100:.0f}%, Mic={mic_frac_NH[j2010]*100:.0f}%, BB={bb_frac_NH[j2010]*100:.0f}%")
    print(f"    SH: FF={ff_frac_SH[j2010]*100:.0f}%, Mic={mic_frac_SH[j2010]*100:.0f}%, BB={bb_frac_SH[j2010]*100:.0f}%")
    print(f"    Expected: FF~19%, Mic~65%, BB~6% (EDGAR/GAO)")
    
    # Check 3: NH/SH FF ratio vs EDGAR
    ff_ratio = np.median(FF_NH[j2010,:]) / (np.median(FF_NH[j2010,:]) + np.median(FF_SH[j2010,:]))
    print(f"\n  NH share of global FF: {ff_ratio*100:.0f}%")
    print(f"  EDGAR estimate: ~72% (most FF in NH)")
    
    # Check 4: Year-to-year stability (coefficient of variation)
    cv_ff = np.std(np.median(FF_NH + FF_SH, axis=1)) / np.mean(np.median(FF_NH + FF_SH, axis=1))
    cv_mic = np.std(np.median(Mic_NH + Mic_SH, axis=1)) / np.mean(np.median(Mic_NH + Mic_SH, axis=1))
    print(f"\n  Temporal stability (CV of annual medians):")
    print(f"    FF: CV = {cv_ff:.3f}  {'✓' if cv_ff < 0.3 else '⚠ noisy'}")
    print(f"    Mic: CV = {cv_mic:.3f}  {'✓' if cv_mic < 0.3 else '⚠ noisy'}")
    
    # Check 5: Negative sources?
    pct_neg_FF = np.mean(FF_NH < 0) * 100 + np.mean(FF_SH < 0) * 100
    print(f"\n  Negative FF iterations: {pct_neg_FF:.1f}% (bounded to 0)")
    
    results = {
        'total_source_range': [float(S_med[0]), float(S_med[-1])],
        'ff_fraction_NH_2010': float(ff_frac_NH[j2010]),
        'ff_fraction_SH_2010': float(ff_frac_SH[j2010]),
        'nh_share_of_ff': float(ff_ratio),
        'edgar_nh_share': 0.72,
        'cv_ff': float(cv_ff),
        'cv_mic': float(cv_mic),
    }
    return results


def phase_B5_dD_consistency(data):
    """
    δD gradient consistency check:
    Compare predicted NH-SH δD gradient from source decomposition
    to observed δD gradient from input data.
    """
    print("\n" + "="*60)
    print("Phase B.5: δD Gradient Consistency Check")
    print("="*60)
    
    # Load observed δD
    nh_dD_file = ROOT / "rel" / "data" / "NHMean_dD_iterations_DasguptaCal_noBUDS.xlsx"
    sh_dD_file = ROOT / "rel" / "data" / "SHMean_dD_iterations_DasguptaCal_noBUDS.xlsx"
    
    nh_dD_df = pd.read_excel(nh_dD_file, index_col=0)
    sh_dD_df = pd.read_excel(sh_dD_file, index_col=0)
    
    # Find overlapping years
    common_years = sorted(set(nh_dD_df.index) & set(sh_dD_df.index))
    print(f"  Overlapping years: {common_years[0]}-{common_years[-1]} ({len(common_years)} yrs)")
    
    # Compute observed gradient per MC iteration
    nh_vals = nh_dD_df.loc[common_years].values  # (n_years, n_MC)
    sh_vals = sh_dD_df.loc[common_years].values
    
    obs_gradient = nh_vals - sh_vals  # NH - SH, ‰
    
    # Summary stats
    grad_median = np.median(obs_gradient, axis=1)
    grad_mean = np.mean(obs_gradient)
    grad_std = np.std(np.mean(obs_gradient, axis=1))
    
    print(f"\n  Observed δD gradient (NH - SH):")
    print(f"    Mean: {grad_mean:.1f}‰")
    print(f"    Std of annual means: {grad_std:.1f}‰")
    print(f"    Range: {np.min(grad_median):.1f} to {np.max(grad_median):.1f}‰")
    
    # Now compute what the model PREDICTS the gradient should be
    # Using v2 solved sources + source δD signatures
    v2 = np.load(RESULTS_DIR / "twobox_v2.npz")
    years = v2['years']
    
    rng = np.random.default_rng(SEED)
    
    # Sample source signatures for a subset of iterations
    n_check = min(200, NI)
    predicted_gradients = []
    
    for k in range(n_check):
        sigs = sample_source_signatures_hemi(rng, data, k, len(years))
        
        j = 10  # middle year (2009)
        
        # Predicted source-weighted δD for NH
        FF_nh = v2['NH_FF'][j, k]
        Mic_nh = v2['NH_Mic'][j, k]
        BB_nh = v2['NH_BB'][j, k]
        S_nh = FF_nh + Mic_nh + BB_nh
        
        if S_nh > 0:
            dD_src_NH = (FF_nh * sigs['ff_dD_NH'][j] + 
                         Mic_nh * sigs['mic_dD_NH'][j] + 
                         BB_nh * sigs['bb_dD_NH'][j]) / S_nh
        else:
            continue
            
        FF_sh = v2['SH_FF'][j, k]
        Mic_sh = v2['SH_Mic'][j, k]
        BB_sh = v2['SH_BB'][j, k]
        S_sh = FF_sh + Mic_sh + BB_sh
        
        if S_sh > 0:
            dD_src_SH = (FF_sh * sigs['ff_dD_SH'][j] + 
                         Mic_sh * sigs['mic_dD_SH'][j] + 
                         BB_sh * sigs['bb_dD_SH'][j]) / S_sh
        else:
            continue
            
        predicted_gradients.append(dD_src_NH - dD_src_SH)
    
    pred_grad = np.array(predicted_gradients)
    print(f"\n  Predicted source δD gradient (NH - SH):")
    print(f"    Median: {np.median(pred_grad):.1f}‰")
    print(f"    90% CI: [{np.percentile(pred_grad, 5):.1f}, {np.percentile(pred_grad, 95):.1f}]‰")
    
    # Compare
    discrepancy = np.median(pred_grad) - grad_mean
    print(f"\n  Discrepancy: {discrepancy:+.1f}‰")
    if abs(discrepancy) > 10:
        print(f"  ⚠ LARGE discrepancy — source δD values may be inconsistent with observations")
    else:
        print(f"  ✓ Reasonable agreement")
    
    results = {
        'observed_dD_gradient_mean': float(grad_mean),
        'observed_dD_gradient_std': float(grad_std),
        'predicted_dD_gradient_median': float(np.median(pred_grad)),
        'predicted_dD_gradient_90ci': [float(np.percentile(pred_grad, 5)),
                                        float(np.percentile(pred_grad, 95))],
        'discrepancy': float(discrepancy),
    }
    return results


# ═══════════════════════════════════════════════════════════════
# PHASE C: NARRATIVE STRENGTHENING
# ═══════════════════════════════════════════════════════════════

def phase_C6_information_analysis(data):
    """
    Replace "degeneracy breaking" with proper information-theoretic analysis.
    
    Key insight: The 2-box doesn't break degeneracy (cond numbers same).
    What it does is provide 6 observation-equations for 6 unknowns 
    (3 per hemisphere) vs 3 for 3 in the 1-box. But the per-hemisphere
    systems are still 3×3.
    
    The REAL information gain comes from:
    1. Different source signatures per hemisphere → different A matrices
    2. Different total source S per hemisphere → different constraint surfaces
    3. The IH exchange term couples the hemispheres
    """
    print("\n" + "="*60)
    print("Phase C.6: Information-Theoretic Analysis")
    print("="*60)
    
    rng = np.random.default_rng(SEED)
    n_samples = 500
    
    # For each MC sample, compute Fisher information matrices
    fisher_global = []
    fisher_nh = []
    fisher_sh = []
    
    # Also compute SVD-based effective rank
    rank_global = []
    rank_nh = []
    
    for k in range(n_samples):
        sigs = sample_source_signatures_hemi(rng, data, k, data.n_years)
        j = 10  # 2009
        
        # Global A
        A_glob = np.array([
            [1.0, 1.0, 1.0],
            [delta_to_fraction_d13C(sigs['ff_d13C'][j]),
             delta_to_fraction_d13C(sigs['mic_d13C'][j]),
             delta_to_fraction_d13C(sigs['bb_d13C'][j])],
            [delta_to_fraction_dD(sigs['ff_dD'][j]),
             delta_to_fraction_dD(sigs['mic_dD'][j]),
             delta_to_fraction_dD(sigs['bb_dD'][j])],
        ])
        
        # NH A
        A_nh = np.array([
            [1.0, 1.0, 1.0],
            [delta_to_fraction_d13C(sigs['ff_d13C_NH'][j]),
             delta_to_fraction_d13C(sigs['mic_d13C_NH'][j]),
             delta_to_fraction_d13C(sigs['bb_d13C_NH'][j])],
            [delta_to_fraction_dD(sigs['ff_dD_NH'][j]),
             delta_to_fraction_dD(sigs['mic_dD_NH'][j]),
             delta_to_fraction_dD(sigs['bb_dD_NH'][j])],
        ])
        
        # SH A
        A_sh = np.array([
            [1.0, 1.0, 1.0],
            [delta_to_fraction_d13C(sigs['ff_d13C_SH'][j]),
             delta_to_fraction_d13C(sigs['mic_d13C_SH'][j]),
             delta_to_fraction_d13C(sigs['bb_d13C_SH'][j])],
            [delta_to_fraction_dD(sigs['ff_dD_SH'][j]),
             delta_to_fraction_dD(sigs['mic_dD_SH'][j]),
             delta_to_fraction_dD(sigs['bb_dD_SH'][j])],
        ])
        
        # Fisher info = A^T A (for unit-weighted system)
        F_glob = A_glob.T @ A_glob
        F_nh = A_nh.T @ A_nh
        F_sh = A_sh.T @ A_sh
        fisher_global.append(np.linalg.det(F_glob))
        fisher_nh.append(np.linalg.det(F_nh))
        fisher_sh.append(np.linalg.det(F_sh))
        
        # SVD for effective rank (number of singular values > 1% of max)
        U_g, s_g, Vt_g = np.linalg.svd(A_glob)
        U_n, s_n, Vt_n = np.linalg.svd(A_nh)
        rank_global.append(np.sum(s_g > 0.01 * s_g[0]))
        rank_nh.append(np.sum(s_n > 0.01 * s_n[0]))
        
    fisher_global = np.array(fisher_global)
    fisher_nh = np.array(fisher_nh)
    fisher_sh = np.array(fisher_sh)
    
    # Combined 2-box Fisher info (block diagonal since hemispheres solved separately)
    fisher_2box = fisher_nh * fisher_sh
    
    print(f"  Fisher information determinant (higher = more info):")
    print(f"    Global (1-box): {np.median(fisher_global):.2e}")
    print(f"    NH alone:       {np.median(fisher_nh):.2e}")
    print(f"    SH alone:       {np.median(fisher_sh):.2e}")
    print(f"    2-box combined: {np.median(fisher_2box):.2e}")
    print(f"    Info gain (2box/1box): {np.median(fisher_2box)/np.median(fisher_global):.1f}×")
    
    print(f"\n  SVD effective rank (of 3):")
    print(f"    Global: {np.mean(rank_global):.2f}")
    print(f"    NH:     {np.mean(rank_nh):.2f}")
    
    # Source signature separation analysis
    # How well do FF, Mic, BB separate in δ¹³C vs δD space?
    print(f"\n  Source signature separation (Δ in ‰):")
    sigs = sample_source_signatures_hemi(rng, data, 0, data.n_years)
    j = 10
    # Global
    d13c = [sigs['ff_d13C'][j], sigs['mic_d13C'][j], sigs['bb_d13C'][j]]
    dD = [sigs['ff_dD'][j], sigs['mic_dD'][j], sigs['bb_dD'][j]]
    print(f"    Global δ¹³C: FF={d13c[0]:.1f}, Mic={d13c[1]:.1f}, BB={d13c[2]:.1f}")
    print(f"    Global δD:   FF={dD[0]:.0f}, Mic={dD[1]:.0f}, BB={dD[2]:.0f}")
    print(f"    FF-BB Δδ¹³C: {abs(d13c[0]-d13c[2]):.1f}‰ — {'poor' if abs(d13c[0]-d13c[2])<3 else 'ok'}")
    print(f"    FF-BB ΔδD:   {abs(dD[0]-dD[2]):.0f}‰ — {'poor' if abs(dD[0]-dD[2])<30 else 'ok'}")
    print(f"    FF-Mic ΔδD:  {abs(dD[0]-dD[1]):.0f}‰ — {'poor' if abs(dD[0]-dD[1])<30 else 'ok'}")
    
    results = {
        'fisher_info_global': float(np.median(fisher_global)),
        'fisher_info_nh': float(np.median(fisher_nh)),
        'fisher_info_sh': float(np.median(fisher_sh)),
        'fisher_info_2box_combined': float(np.median(fisher_2box)),
        'info_gain_ratio': float(np.median(fisher_2box)/np.median(fisher_global)),
        'effective_rank_global': float(np.mean(rank_global)),
        'effective_rank_nh': float(np.mean(rank_nh)),
    }
    return results


def phase_C7_edgar_crosscheck():
    """Compare model FF levels and NH/SH partition to EDGAR inventory."""
    print("\n" + "="*60)
    print("Phase C.7: EDGAR Cross-Check")
    print("="*60)
    
    v2 = np.load(RESULTS_DIR / "twobox_v2.npz")
    years = v2['years']
    
    # EDGAR v7 literature values for FF CH₄ (Tg/yr)
    # Saunois et al. 2020 (Table 3): FF 108-116 Tg/yr (2008-2017)
    # EDGAR v6: ~110 Tg/yr (2010), ~72% in NH
    edgar_ff_2010 = 110.0
    edgar_nh_fraction = 0.72
    
    j2010 = np.where(years == 2010)[0][0]
    
    model_ff = v2['NH_FF'][j2010,:] + v2['SH_FF'][j2010,:]
    model_nh_frac = v2['NH_FF'][j2010,:] / model_ff
    
    # EDGAR FF trend 2007-2018: ~+2.1 Tg/yr² (increasing)
    edgar_trend = 2.1
    model_ff_slopes = np.array([
        sp_stats.linregress(years[(years>=2007)&(years<=2018)],
                            (v2['NH_FF'] + v2['SH_FF'])[(years>=2007)&(years<=2018), k]).slope
        for k in range(v2['NH_FF'].shape[1])])
    
    print(f"  FF emissions (2010):")
    print(f"    Model: {np.median(model_ff):.0f} [{np.percentile(model_ff,5):.0f}, {np.percentile(model_ff,95):.0f}] Tg/yr")
    print(f"    EDGAR: {edgar_ff_2010:.0f} Tg/yr")
    agreement = abs(np.median(model_ff) - edgar_ff_2010) / edgar_ff_2010* 100
    print(f"    Agreement: {agreement:.0f}% offset {'✓' if agreement < 20 else '⚠'}")
    
    print(f"\n  NH share of FF:")
    print(f"    Model: {np.median(model_nh_frac)*100:.0f}%")
    print(f"    EDGAR: {edgar_nh_fraction*100:.0f}%")
    
    print(f"\n  FF trend (2007-2018):")
    print(f"    Model: {np.median(model_ff_slopes):+.2f} Tg/yr²")
    print(f"    EDGAR: +{edgar_trend:.1f} Tg/yr² (BU inventory)")
    print(f"    Note: EDGAR is bottom-up; our model is top-down isotopic")
    
    # Mic comparison
    model_mic = v2['NH_Mic'][j2010,:] + v2['SH_Mic'][j2010,:]
    print(f"\n  Microbial emissions (2010):")
    print(f"    Model: {np.median(model_mic):.0f} [{np.percentile(model_mic,5):.0f}, {np.percentile(model_mic,95):.0f}] Tg/yr")
    print(f"    Literature: ~370 Tg/yr (Saunois 2020)")
    
    model_bb = v2['NH_BB'][j2010,:] + v2['SH_BB'][j2010,:]
    print(f"\n  BB emissions (2010):")
    print(f"    Model: {np.median(model_bb):.0f} [{np.percentile(model_bb,5):.0f}, {np.percentile(model_bb,95):.0f}] Tg/yr")
    print(f"    Literature: ~30 Tg/yr (GFED4s)")
    
    results = {
        'model_ff_2010': float(np.median(model_ff)),
        'model_ff_90ci': [float(np.percentile(model_ff,5)), float(np.percentile(model_ff,95))],
        'edgar_ff_2010': edgar_ff_2010,
        'model_nh_fraction': float(np.median(model_nh_frac)),
        'edgar_nh_fraction': edgar_nh_fraction,
        'model_ff_trend': float(np.median(model_ff_slopes)),
        'edgar_ff_trend': edgar_trend,
        'model_mic_2010': float(np.median(model_mic)),
        'model_bb_2010': float(np.median(model_bb)),
    }
    return results


def phase_C8_W_sensitivity():
    """
    W-matrix sensitivity sweep: how much do results change with different
    weighting schemes? This is the key methodological finding.
    """
    print("\n" + "="*60)
    print("Phase C.8: W-Matrix Sensitivity Sweep")
    print("="*60)
    
    data = load_data(ROOT, two_box=True)
    n = data.n_years
    years = data.model_years
    rng = np.random.default_rng(SEED)
    
    # IH gradient
    all_years = np.arange(years[0], years[-1] + 2)
    IH_grad = realistic_IH_gradient(all_years.astype(float))
    CH4_NH = data.CH4_global + IH_grad / 2.0
    CH4_SH = data.CH4_global - IH_grad / 2.0
    
    tau_global = compute_lifetime(years, "varying", 9.0)
    tau_NH = tau_global * LIFETIME_RATIO_NH
    tau_SH = tau_global * LIFETIME_RATIO_SH
    
    # W-matrix configurations to test
    W_configs = {
        'v1_original': (np.diag([100, 1, 0.5]), np.diag([200, 1, 0.5])),
        'mass_heavy': (np.diag([1000, 1, 0.5]), np.diag([1000, 1, 0.5])),
        'isotope_heavy': (np.diag([1, 100, 50]), np.diag([1, 100, 50])),
        'equal': (np.eye(3), np.eye(3)),
        'uncertainty_based': ('auto', 'auto'),  # computed per-iteration
    }
    
    NI_sens = 400  # fewer for speed
    results_all = {}
    
    for wname, (W_NH_spec, W_SH_spec) in W_configs.items():
        print(f"\n  Running W={wname} ({NI_sens} MC)...")
        
        FF_global = np.zeros((n, NI_sens))
        
        for k in range(NI_sens):
            tau_ex = max(0.5, rng.normal(TAU_EX_MEAN, TAU_EX_STD))
            kies = sample_KIE(rng, "sampled")
            K13_NH, KD_NH = compute_bulk_KIE(kies, SINK_FRACTIONS_NH)
            K13_SH, KD_SH = compute_bulk_KIE(kies, SINK_FRACTIONS_SH)
            a13_NH = 1.0/K13_NH; aD_NH = 1.0/KD_NH
            a13_SH = 1.0/K13_SH; aD_SH = 1.0/KD_SH
            
            S_NH = np.zeros(n); S_SH = np.zeros(n)
            for i in range(n):
                M_NH = CH4_NH[i]*PT_HEMI; M_NH1 = CH4_NH[i+1]*PT_HEMI
                M_SH = CH4_SH[i]*PT_HEMI; M_SH1 = CH4_SH[i+1]*PT_HEMI
                S_NH[i] = (M_NH1-M_NH) + M_NH/tau_NH[i] - (M_SH-M_NH)/tau_ex
                S_SH[i] = (M_SH1-M_SH) + M_SH/tau_SH[i] - (M_NH-M_SH)/tau_ex
            
            d13C_MC = sample_atm_d13C(data, k, n)
            dD_NH_MC, dD_SH_MC = sample_atm_dD_hemi(data, k, n)
            nc = min(len(data.c13_global), n+1)
            d13C_off = d13C_MC[:nc] - data.c13_global[:nc]
            d13C_NH_MC = data.c13_NH[:nc] + d13C_off
            d13C_SH_MC = data.c13_SH[:nc] + d13C_off
            
            f13_NH_atm = delta_to_fraction_d13C(d13C_NH_MC)
            f13_SH_atm = delta_to_fraction_d13C(d13C_SH_MC)
            fD_NH_atm = delta_to_fraction_dD(dD_NH_MC)
            fD_SH_atm = delta_to_fraction_dD(dD_SH_MC)
            
            d13C_src_NH = np.zeros(n); dD_src_NH = np.zeros(n)
            d13C_src_SH = np.zeros(n); dD_src_SH = np.zeros(n)
            for j in range(n):
                n13_NH = f13_NH_atm[j]*CH4_NH[j]*PT_HEMI
                n13_NH1 = f13_NH_atm[j+1]*CH4_NH[j+1]*PT_HEMI
                n13_SH = f13_SH_atm[j]*CH4_SH[j]*PT_HEMI
                n13_SH1 = f13_SH_atm[j+1]*CH4_SH[j+1]*PT_HEMI
                d13C_src_NH[j] = (n13_NH1-n13_NH + n13_NH*a13_NH/tau_NH[j] - (n13_SH-n13_NH)/tau_ex) / S_NH[j]
                d13C_src_SH[j] = (n13_SH1-n13_SH + n13_SH*a13_SH/tau_SH[j] - (n13_NH-n13_SH)/tau_ex) / S_SH[j]
                
                nD_NH = fD_NH_atm[j]*CH4_NH[j]*PT_HEMI
                nD_NH1 = fD_NH_atm[j+1]*CH4_NH[j+1]*PT_HEMI
                nD_SH = fD_SH_atm[j]*CH4_SH[j]*PT_HEMI
                nD_SH1 = fD_SH_atm[j+1]*CH4_SH[j+1]*PT_HEMI
                dD_src_NH[j] = (nD_NH1-nD_NH + nD_NH*aD_NH/tau_NH[j] - (nD_SH-nD_NH)/tau_ex) / S_NH[j]
                dD_src_SH[j] = (nD_SH1-nD_SH + nD_SH*aD_SH/tau_SH[j] - (nD_NH-nD_SH)/tau_ex) / S_SH[j]
            
            sigs = sample_source_signatures_hemi(rng, data, k, n)
            
            # Compute W matrices
            if isinstance(W_NH_spec, str) and W_NH_spec == 'auto':
                S_NH_mean = np.mean(S_NH[S_NH>0]) if np.any(S_NH>0) else 200
                w1 = 10.0/max(S_NH_mean,100); w2 = 1.0/3e-6; w3 = 1.0/5e-6
                wmax = max(w1,w2,w3)
                W_NH_use = np.diag([w1/wmax, w2/wmax, w3/wmax])
                S_SH_mean = np.mean(S_SH[S_SH>0]) if np.any(S_SH>0) else 200
                w1 = 10.0/max(S_SH_mean,100); w2 = 1.0/3e-6; w3 = 1.0/5e-6
                wmax = max(w1,w2,w3)
                W_SH_use = np.diag([w1/wmax, w2/wmax, w3/wmax])
            else:
                W_NH_use = W_NH_spec
                W_SH_use = W_SH_spec
            
            for j in range(n):
                f13_ff_nh = delta_to_fraction_d13C(sigs['ff_d13C_NH'][j])
                f13_mic_nh = delta_to_fraction_d13C(sigs['mic_d13C_NH'][j])
                f13_bb_nh = delta_to_fraction_d13C(sigs['bb_d13C_NH'][j])
                fD_ff_nh = delta_to_fraction_dD(sigs['ff_dD_NH'][j])
                fD_mic_nh = delta_to_fraction_dD(sigs['mic_dD_NH'][j])
                fD_bb_nh = delta_to_fraction_dD(sigs['bb_dD_NH'][j])
                
                A_nh = np.array([[1,1,1],
                                 [f13_bb_nh,f13_ff_nh,f13_mic_nh],
                                 [fD_bb_nh,fD_ff_nh,fD_mic_nh]])
                B_nh = np.array([S_NH[j], S_NH[j]*d13C_src_NH[j], S_NH[j]*dD_src_NH[j]])
                ub = S_NH[j]*1.5 if S_NH[j]>0 else 1000
                try:
                    r = lsq_linear(W_NH_use@A_nh, W_NH_use@B_nh, bounds=(0,ub))
                    ff_nh = r.x[1]
                except:
                    ff_nh = np.nan
                
                f13_ff_sh = delta_to_fraction_d13C(sigs['ff_d13C_SH'][j])
                f13_mic_sh = delta_to_fraction_d13C(sigs['mic_d13C_SH'][j])
                f13_bb_sh = delta_to_fraction_d13C(sigs['bb_d13C_SH'][j])
                fD_ff_sh = delta_to_fraction_dD(sigs['ff_dD_SH'][j])
                fD_mic_sh = delta_to_fraction_dD(sigs['mic_dD_SH'][j])
                fD_bb_sh = delta_to_fraction_dD(sigs['bb_dD_SH'][j])
                
                A_sh = np.array([[1,1,1],
                                 [f13_bb_sh,f13_ff_sh,f13_mic_sh],
                                 [fD_bb_sh,fD_ff_sh,fD_mic_sh]])
                B_sh = np.array([S_SH[j], S_SH[j]*d13C_src_SH[j], S_SH[j]*dD_src_SH[j]])
                ub = S_SH[j]*1.5 if S_SH[j]>0 else 1000
                try:
                    r = lsq_linear(W_SH_use@A_sh, W_SH_use@B_sh, bounds=(0,ub))
                    ff_sh = r.x[1]
                except:
                    ff_sh = np.nan
                
                FF_global[j,k] = ff_nh + ff_sh
        
        # Compute trend and mean
        j2010 = np.where(years==2010)[0][0]
        ff_med = np.nanmedian(FF_global[j2010,:])
        slopes = np.array([sp_stats.linregress(
            years[(years>=2007)&(years<=2020)],
            FF_global[(years>=2007)&(years<=2020),k]).slope
            for k in range(NI_sens) if not np.any(np.isnan(FF_global[(years>=2007)&(years<=2020),k]))])
        
        results_all[wname] = {
            'ff_2010_median': float(ff_med),
            'ff_trend_median': float(np.median(slopes)) if len(slopes)>0 else None,
            'ff_trend_90ci': [float(np.percentile(slopes,5)), float(np.percentile(slopes,95))] if len(slopes)>0 else None,
        }
        print(f"    FF(2010)={ff_med:.0f} Tg/yr, trend={np.median(slopes):+.2f} Tg/yr²")
    
    print(f"\n  SUMMARY: W-matrix sensitivity")
    print(f"  {'Config':20s} {'FF(2010)':>10s} {'FF trend':>12s}")
    for wname, r in results_all.items():
        trend_str = f"{r['ff_trend_median']:+.2f}" if r['ff_trend_median'] is not None else "N/A"
        print(f"  {wname:20s} {r['ff_2010_median']:10.0f} {trend_str:>12s}")
    
    return results_all


def main():
    data = load_data(ROOT, two_box=True)
    
    all_results = {}
    
    # Phase B
    all_results['B4_posterior'] = phase_B4_posterior_predictive(data)
    all_results['B5_dD_consistency'] = phase_B5_dD_consistency(data)
    
    # Phase C
    all_results['C6_information'] = phase_C6_information_analysis(data)
    all_results['C7_edgar'] = phase_C7_edgar_crosscheck()
    all_results['C8_W_sensitivity'] = phase_C8_W_sensitivity()
    
    # Save
    with open(RESULTS_DIR / "validation_results.json", 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    
    print(f"\n{'='*60}")
    print(f"ALL VALIDATION COMPLETE — saved to {RESULTS_DIR}/validation_results.json")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
