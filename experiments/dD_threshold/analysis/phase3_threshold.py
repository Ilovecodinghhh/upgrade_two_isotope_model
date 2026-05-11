#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
phase3_threshold.py — The δD Threshold Experiment (Core Result)
================================================================

Systematically inflates the microbial δD source-signature uncertainty and
measures when δD stops providing useful constraint on FF emissions.

This is the central experiment of Title 1: finding the critical threshold
above which δD no longer helps discriminate sources.
"""

import sys
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))
from common import (
    load_data, sample_KIE, compute_bulk_KIE, compute_lifetime,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    smooth_5yr, pad_to_length,
    SINK_FRACTIONS_GLOBAL, SINK_FRACTIONS_NH, SINK_FRACTIONS_SH,
    PT, PT_HEMI, LIFETIME_RATIO_NH, LIFETIME_RATIO_SH,
    DD_IH_OFFSET, TAU_EX_MEAN, TAU_EX_STD,
    BB_NH_FRACTION, BB_SH_FRACTION,
)

OUT_DIR = Path(__file__).resolve().parent.parent / "results" / "phase3_threshold"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def inflate_dD_uncertainty(sigs, rng, multiplier, year_idx):
    """
    Inflate δD signatures by a multiplier around their mean.
    
    Takes an already-sampled signature draw and adds additional spread
    equivalent to (multiplier - 1) × original_std.
    
    The original MC draw has σ ≈ 7-8‰ for Mic_dD; we inflate to simulate
    the larger uncertainties used by Thanwerdas et al. (±128‰).
    """
    if multiplier <= 1.0:
        return sigs
    
    # Additional spread (in ‰) scales roughly with multiplier
    # Baseline σ values from Phase 2: Mic_dD ±8.25‰, FF_dD ±0.70‰, BB_dD ±7.09‰
    extra_mic = rng.normal() * 8.25 * (multiplier - 1)
    extra_ff = rng.normal() * 0.70 * (multiplier - 1)  # Small (FF_dD well-constrained)
    extra_bb = rng.normal() * 7.09 * (multiplier - 1)
    
    sigs_new = dict(sigs)
    sigs_new['mic_dD'] = sigs['mic_dD'] + extra_mic
    sigs_new['ff_dD'] = sigs['ff_dD'] + extra_ff
    sigs_new['bb_dD'] = sigs['bb_dD'] + extra_bb
    return sigs_new


def run_twobox_with_inflation(data, multiplier, n_iter=500, seed=42, mode="dual"):
    """
    Run 2-box model with inflated δD uncertainty.
    
    Returns compiled FF arrays (NH, SH, Global).
    """
    from common import sample_source_signatures, sample_atm_d13C, sample_atm_dD
    
    rng = np.random.default_rng(seed)
    n = data.n_years
    years = data.model_years
    
    CH4_NH, CH4_SH = data.CH4_NH, data.CH4_SH
    c13_NH, c13_SH = data.c13_NH, data.c13_SH
    c13_glob = data.c13_global
    
    tau_global = compute_lifetime(years, "varying", 9.0)
    tau_NH = tau_global * LIFETIME_RATIO_NH
    tau_SH = tau_global * LIFETIME_RATIO_SH
    
    BB_hemi_NH = data.BB_global_mean * BB_NH_FRACTION
    BB_hemi_SH = data.BB_global_mean * BB_SH_FRACTION
    
    FF_NH = np.zeros((n, n_iter)); FF_SH = np.zeros((n, n_iter))
    Mic_NH = np.zeros((n, n_iter)); Mic_SH = np.zeros((n, n_iter))
    
    W_NH = np.diag([100.0, 1.0, 0.5])
    W_SH = np.diag([200.0, 1.0, 0.5])
    
    for k in range(n_iter):
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
            ex_NH = (M_SH-M_NH)/tau_ex; ex_SH = (M_NH-M_SH)/tau_ex
            S_NH[i] = (M_NH1-M_NH) + M_NH/tau_NH[i] - ex_NH
            S_SH[i] = (M_SH1-M_SH) + M_SH/tau_SH[i] - ex_SH
        
        d13C_glob = sample_atm_d13C(data, k, n)
        nc = min(len(c13_glob), n+1)
        d13C_off = d13C_glob[:nc] - c13_glob[:nc]
        d13C_NH_MC = c13_NH[:nc] + d13C_off
        d13C_SH_MC = c13_SH[:nc] + d13C_off
        
        f13_NH_atm = delta_to_fraction_d13C(d13C_NH_MC)
        f13_SH_atm = delta_to_fraction_d13C(d13C_SH_MC)
        
        d13C_src_NH = np.zeros(n); d13C_src_SH = np.zeros(n)
        for j in range(n):
            n13_NH = f13_NH_atm[j]*CH4_NH[j]*PT_HEMI
            n13_NH1 = f13_NH_atm[j+1]*CH4_NH[j+1]*PT_HEMI
            n13_SH = f13_SH_atm[j]*CH4_SH[j]*PT_HEMI
            n13_SH1 = f13_SH_atm[j+1]*CH4_SH[j+1]*PT_HEMI
            ex13_NH = (n13_SH-n13_NH)/tau_ex; ex13_SH = (n13_NH-n13_SH)/tau_ex
            d13C_src_NH[j] = (n13_NH1-n13_NH + n13_NH*a13_NH/tau_NH[j] - ex13_NH)/S_NH[j]
            d13C_src_SH[j] = (n13_SH1-n13_SH + n13_SH*a13_SH/tau_SH[j] - ex13_SH)/S_SH[j]
        
        sigs = sample_source_signatures(rng, data, k, n)
        
        # INFLATE δD uncertainty
        if multiplier > 1.0 and mode == "dual":
            sigs = inflate_dD_uncertainty(sigs, rng, multiplier, year_idx=16)
        
        f13_bb = delta_to_fraction_d13C(sigs['bb_d13C'])
        f13_ff = delta_to_fraction_d13C(sigs['ff_d13C'])
        f13_mic = delta_to_fraction_d13C(sigs['mic_d13C'])
        
        if mode == "dual":
            dD_glob = sample_atm_dD(data, k, n)
            dD_NH_MC = dD_glob - DD_IH_OFFSET
            dD_SH_MC = dD_glob + DD_IH_OFFSET
            fD_NH_atm = delta_to_fraction_dD(dD_NH_MC)
            fD_SH_atm = delta_to_fraction_dD(dD_SH_MC)
            
            dD_src_NH = np.zeros(n); dD_src_SH = np.zeros(n)
            for j in range(n):
                nD_NH = fD_NH_atm[j]*CH4_NH[j]*PT_HEMI
                nD_NH1 = fD_NH_atm[j+1]*CH4_NH[j+1]*PT_HEMI
                nD_SH = fD_SH_atm[j]*CH4_SH[j]*PT_HEMI
                nD_SH1 = fD_SH_atm[j+1]*CH4_SH[j+1]*PT_HEMI
                exD_NH = (nD_SH-nD_NH)/tau_ex; exD_SH = (nD_NH-nD_SH)/tau_ex
                dD_src_NH[j] = (nD_NH1-nD_NH + nD_NH*aD_NH/tau_NH[j] - exD_NH)/S_NH[j]
                dD_src_SH[j] = (nD_SH1-nD_SH + nD_SH*aD_SH/tau_SH[j] - exD_SH)/S_SH[j]
            
            fD_bb = delta_to_fraction_dD(sigs['bb_dD'])
            fD_ff = delta_to_fraction_dD(sigs['ff_dD'])
            fD_mic = delta_to_fraction_dD(sigs['mic_dD'])
            
            for j in range(n):
                A = np.array([
                    [1.0, 1.0, 1.0],
                    [f13_bb[j], f13_ff[j], f13_mic[j]],
                    [fD_bb[j], fD_ff[j], fD_mic[j]],
                ])
                B_nh = np.array([S_NH[j], S_NH[j]*d13C_src_NH[j], S_NH[j]*dD_src_NH[j]])
                try:
                    res = lsq_linear(W_NH @ A, W_NH @ B_nh, bounds=(0, S_NH[j]*1.5))
                    x = res.x
                except: x = np.array([np.nan]*3)
                FF_NH[j,k] = x[1]; Mic_NH[j,k] = x[2]
                
                B_sh = np.array([S_SH[j], S_SH[j]*d13C_src_SH[j], S_SH[j]*dD_src_SH[j]])
                try:
                    res = lsq_linear(W_SH @ A, W_SH @ B_sh, bounds=(0, S_SH[j]*1.5))
                    x = res.x
                except: x = np.array([np.nan]*3)
                FF_SH[j,k] = x[1]; Mic_SH[j,k] = x[2]
        
        else:  # d13C_only
            for j in range(n):
                denom = f13_ff[j] - f13_mic[j]
                if abs(denom) < 1e-15:
                    FF_NH[j,k] = np.nan; Mic_NH[j,k] = np.nan
                    FF_SH[j,k] = np.nan; Mic_SH[j,k] = np.nan
                    continue
                # NH
                S_rem = S_NH[j] - BB_hemi_NH
                rhs = S_NH[j]*d13C_src_NH[j] - BB_hemi_NH*f13_bb[j]
                FF_NH[j,k] = (rhs - S_rem*f13_mic[j]) / denom
                Mic_NH[j,k] = S_rem - FF_NH[j,k]
                # SH
                S_rem = S_SH[j] - BB_hemi_SH
                rhs = S_SH[j]*d13C_src_SH[j] - BB_hemi_SH*f13_bb[j]
                FF_SH[j,k] = (rhs - S_rem*f13_mic[j]) / denom
                Mic_SH[j,k] = S_rem - FF_SH[j,k]
    
    return {
        'FF_NH': FF_NH, 'FF_SH': FF_SH, 'FF_G': FF_NH + FF_SH,
        'Mic_NH': Mic_NH, 'Mic_SH': Mic_SH, 'Mic_G': Mic_NH + Mic_SH,
    }


def compute_metrics(results, label):
    """Compute summary metrics: CI width, std, % non-physical."""
    FF_s = smooth_5yr(results['FF_G'])
    # CI width over post-2007 period (indices 8+)
    ci_5 = np.nanpercentile(FF_s[8:], 5, axis=1)
    ci_95 = np.nanpercentile(FF_s[8:], 95, axis=1)
    ci_width = ci_95 - ci_5
    
    return {
        'label': label,
        'FF_mean': float(np.nanmean(FF_s[8:])),
        'FF_std_mean': float(np.nanmean(np.nanstd(FF_s[8:], axis=1))),
        'FF_CI95_width_mean': float(np.nanmean(ci_width)),
        'FF_CI95_width_2015': float(ci_width[16-8] if len(ci_width) > 8 else np.nan),
        'pct_negative': float(np.sum(FF_s < 0) / FF_s.size * 100),
    }


def main():
    print("=" * 70)
    print("PHASE 3: δD THRESHOLD EXPERIMENT")
    print("=" * 70)
    
    data = load_data(REPO_ROOT, two_box=True)
    
    # Multipliers to test — from below baseline (0.5) through extreme (15)
    # At multiplier=1: Mic_dD σ ≈ 8‰
    # At multiplier=15: Mic_dD σ ≈ 120‰ (≈ Thanwerdas prior)
    multipliers = [0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 16.0]
    N_ITER = 500  # Reduced for speed in sensitivity sweep
    
    all_results = {}
    
    # Reference: δ¹³C-only (no δD — shouldn't change with multiplier)
    print("\n[Baseline] 2-box δ¹³C-only...")
    ref_c = run_twobox_with_inflation(data, 1.0, N_ITER, 42, mode="d13C_only")
    ref_c_metrics = compute_metrics(ref_c, "d13C_only_ref")
    all_results['d13C_only_ref'] = ref_c_metrics
    print(f"  FF CI width: {ref_c_metrics['FF_CI95_width_mean']:.1f} Tg/yr")
    
    # Sweep multipliers with dual isotopes
    for mult in multipliers:
        print(f"\n[Multiplier = {mult:.1f}×] 2-box dual-isotope with inflated δD...")
        res = run_twobox_with_inflation(data, mult, N_ITER, 42, mode="dual")
        m = compute_metrics(res, f"dual_mult_{mult}")
        all_results[f'dual_mult_{mult}'] = m
        print(f"  FF CI width: {m['FF_CI95_width_mean']:.1f} Tg/yr  "
              f"(approx Mic δD σ ≈ {8.25*mult:.0f}‰)")
    
    # === ANALYSIS ===
    # Constraint improvement = (d13C_width - dual_width) / d13C_width
    ref_width = ref_c_metrics['FF_CI95_width_mean']
    
    improvements = {}
    threshold_passed = None
    for mult in multipliers:
        dual_width = all_results[f'dual_mult_{mult}']['FF_CI95_width_mean']
        improvement_pct = (ref_width - dual_width) / ref_width * 100
        mic_dD_sigma = 8.25 * mult  # approximate
        improvements[mult] = {
            'multiplier': mult,
            'mic_dD_sigma_permil': mic_dD_sigma,
            'dual_CI_width': dual_width,
            'd13C_CI_width': ref_width,
            'improvement_pct': improvement_pct,
            'dD_adds_value': improvement_pct >= 10.0,  # Criterion A
        }
        if threshold_passed is None and improvement_pct < 10.0:
            threshold_passed = mult
    
    # Print summary
    print(f"\n{'='*70}")
    print("THRESHOLD RESULTS")
    print(f"{'='*70}")
    print(f"\n{'Mult':>6} {'σ(Mic δD)':>12} {'Dual CI':>10} {'δ¹³C CI':>10} {'Improve%':>10} {'Helps?':>8}")
    for mult in multipliers:
        i = improvements[mult]
        flag = "✓" if i['dD_adds_value'] else "✗"
        print(f"{mult:>6.1f} {i['mic_dD_sigma_permil']:>10.1f}‰   "
              f"{i['dual_CI_width']:>8.1f}  {i['d13C_CI_width']:>8.1f}  "
              f"{i['improvement_pct']:>8.1f}%  {flag:>6}")
    
    if threshold_passed:
        print(f"\n  🎯 THRESHOLD: δD stops adding value at multiplier ≈ {threshold_passed:.1f}×")
        print(f"     (Mic δD σ ≈ {8.25*threshold_passed:.0f}‰)")
    else:
        print(f"\n  δD still adds value at all tested multipliers (up to {multipliers[-1]}×)")
    
    # Save
    summary = {
        'multipliers': multipliers,
        'baseline_d13C_CI': ref_width,
        'improvements': improvements,
        'threshold_multiplier': threshold_passed,
        'threshold_mic_dD_sigma_permil': (8.25 * threshold_passed) if threshold_passed else None,
        'thanwerdas_reference_sigma_permil': 128.0,
        'all_metrics': all_results,
    }
    with open(OUT_DIR / "threshold_results.json", 'w') as f:
        json.dump(summary, f, indent=2, default=float)
    
    # CSV
    df = pd.DataFrame([improvements[m] for m in multipliers])
    df.to_csv(OUT_DIR / "threshold_sweep.csv", index=False)
    
    print(f"\nSaved: {OUT_DIR}/threshold_results.json")
    print(f"       {OUT_DIR}/threshold_sweep.csv")


if __name__ == "__main__":
    main()
