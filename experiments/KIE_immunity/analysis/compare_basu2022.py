#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compare_basu2022.py — Compare KIE spread to Basu 2022 (ACP)
=============================================================
v2: Tests three δD configurations (none / offset / real hemispheric)
    and a new finding: whether real hemispheric δD changes which KIE
    (Saueressig vs Cantrell) is preferred.

Basu et al. (2022, ACP): Two TM5-4DVAR inversions differing only in 
OH-¹³C KIE (Saueressig 1.0039 vs Cantrell 1.0054).
Post-2007 FF spread: ~13 Tg/yr.
"""

import sys
import json
from pathlib import Path

import numpy as np
from scipy.optimize import lsq_linear

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))
from common import (
    load_data, compute_bulk_KIE, compute_lifetime,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    sample_source_signatures, sample_source_signatures_hemi,
    sample_atm_d13C, sample_atm_dD, sample_atm_dD_hemi,
    smooth_5yr,
    SINK_FRACTIONS_NH, SINK_FRACTIONS_SH,
    PT_HEMI, LIFETIME_RATIO_NH, LIFETIME_RATIO_SH,
    DD_IH_OFFSET, TAU_EX_MEAN, TAU_EX_STD,
    BB_NH_FRACTION, BB_SH_FRACTION,
)

OUT_DIR = Path(__file__).resolve().parent.parent / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

BASU_KIE_SPREAD = 13.0  # Tg/yr


def run_at_fixed_kie(data, mode, oh_13c_kie, n_iter=400, seed=42,
                     use_real_hemi_dD=False):
    """Run 2-box at FIXED OH-13C KIE. Return FF time series MC array."""
    rng = np.random.default_rng(seed)
    n = data.n_years
    years = data.model_years
    CH4_NH, CH4_SH = data.CH4_NH, data.CH4_SH
    c13_NH, c13_SH = data.c13_NH, data.c13_SH
    c13_glob = data.c13_global
    
    tau_global = compute_lifetime(years, "varying", 9.0)
    tau_NH = tau_global * LIFETIME_RATIO_NH
    tau_SH = tau_global * LIFETIME_RATIO_SH
    BB_NH = data.BB_global_mean * BB_NH_FRACTION
    BB_SH = data.BB_global_mean * BB_SH_FRACTION
    
    kies = {
        'OH_13C': oh_13c_kie,
        'OH_D': 0.5*(1.294+1.327),
        'Cl_13C': 1.066, 'Cl_D': 1.520,
        'Strat_13C': 1.003, 'Strat_D': 1.050,
        'Soil_13C': 1.0201, 'Soil_D': 1.103,
    }
    K13_NH, KD_NH = compute_bulk_KIE(kies, SINK_FRACTIONS_NH)
    K13_SH, KD_SH = compute_bulk_KIE(kies, SINK_FRACTIONS_SH)
    a13_NH, aD_NH = 1.0/K13_NH, 1.0/KD_NH
    a13_SH, aD_SH = 1.0/K13_SH, 1.0/KD_SH
    
    FF_G = np.zeros((n, n_iter))
    W = np.diag([100.0, 1.0, 0.5])
    
    for k in range(n_iter):
        tau_ex = max(0.5, rng.normal(TAU_EX_MEAN, TAU_EX_STD))
        
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
        
        if use_real_hemi_dD:
            sigs = sample_source_signatures_hemi(rng, data, k, n)
        else:
            sigs = sample_source_signatures(rng, data, k, n)

        # δ¹³C source signatures: hemispheric when available, else global for both
        if use_real_hemi_dD and 'ff_d13C_NH' in sigs:
            f13_bb_NH  = delta_to_fraction_d13C(sigs['bb_d13C_NH'])
            f13_ff_NH  = delta_to_fraction_d13C(sigs['ff_d13C_NH'])
            f13_mic_NH = delta_to_fraction_d13C(sigs['mic_d13C_NH'])
            f13_bb_SH  = delta_to_fraction_d13C(sigs['bb_d13C_SH'])
            f13_ff_SH  = delta_to_fraction_d13C(sigs['ff_d13C_SH'])
            f13_mic_SH = delta_to_fraction_d13C(sigs['mic_d13C_SH'])
        else:
            f13_bb  = delta_to_fraction_d13C(sigs['bb_d13C'])
            f13_ff  = delta_to_fraction_d13C(sigs['ff_d13C'])
            f13_mic = delta_to_fraction_d13C(sigs['mic_d13C'])
            f13_bb_NH = f13_bb_SH = f13_bb
            f13_ff_NH = f13_ff_SH = f13_ff
            f13_mic_NH = f13_mic_SH = f13_mic
        
        if mode == "dual":
            if use_real_hemi_dD:
                dD_NH_MC, dD_SH_MC = sample_atm_dD_hemi(data, k, n)
            else:
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
            
            if use_real_hemi_dD:
                fD_bb_NH = delta_to_fraction_dD(sigs['bb_dD_NH'])
                fD_ff_NH = delta_to_fraction_dD(sigs['ff_dD_NH'])
                fD_mic_NH = delta_to_fraction_dD(sigs['mic_dD_NH'])
                fD_bb_SH = delta_to_fraction_dD(sigs['bb_dD_SH'])
                fD_ff_SH = delta_to_fraction_dD(sigs['ff_dD_SH'])
                fD_mic_SH = delta_to_fraction_dD(sigs['mic_dD_SH'])
            else:
                fD_bb = delta_to_fraction_dD(sigs['bb_dD'])
                fD_ff = delta_to_fraction_dD(sigs['ff_dD'])
                fD_mic = delta_to_fraction_dD(sigs['mic_dD'])
                fD_bb_NH = fD_bb_SH = fD_bb
                fD_ff_NH = fD_ff_SH = fD_ff
                fD_mic_NH = fD_mic_SH = fD_mic
            
            for j in range(n):
                for S, d13C_src, dD_src, f13_bb_h, f13_ff_h, f13_mic_h, fD_bb_h, fD_ff_h, fD_mic_h in [
                    (S_NH[j], d13C_src_NH[j], dD_src_NH[j],
                     f13_bb_NH[j], f13_ff_NH[j], f13_mic_NH[j],
                     fD_bb_NH[j], fD_ff_NH[j], fD_mic_NH[j]),
                    (S_SH[j], d13C_src_SH[j], dD_src_SH[j],
                     f13_bb_SH[j], f13_ff_SH[j], f13_mic_SH[j],
                     fD_bb_SH[j], fD_ff_SH[j], fD_mic_SH[j]),
                ]:
                    A = np.array([
                        [1.0, 1.0, 1.0],
                        [f13_bb_h, f13_ff_h, f13_mic_h],
                        [fD_bb_h, fD_ff_h, fD_mic_h],
                    ])
                    B = np.array([S, S*d13C_src, S*dD_src])
                    try:
                        res = lsq_linear(W@A, W@B, bounds=(0, S*1.5))
                        FF_G[j,k] += res.x[1]
                    except: pass
        else:
            for j in range(n):
                denom_NH = f13_ff_NH[j] - f13_mic_NH[j]
                denom_SH = f13_ff_SH[j] - f13_mic_SH[j]
                for S, d13C_src, BB, denom, f13_bb_h, f13_mic_h in [
                    (S_NH[j], d13C_src_NH[j], BB_NH, denom_NH, f13_bb_NH[j], f13_mic_NH[j]),
                    (S_SH[j], d13C_src_SH[j], BB_SH, denom_SH, f13_bb_SH[j], f13_mic_SH[j]),
                ]:
                    if abs(denom) < 1e-15: continue
                    S_rem = S - BB
                    rhs = S*d13C_src - BB*f13_bb_h
                    ff = (rhs - S_rem*f13_mic_h) / denom
                    FF_G[j,k] += max(0, min(ff, S*1.5))
    
    return FF_G


def trend_post_vs_pre2007(FF, years):
    """Mean(2010-2018) - Mean(2000-2006)."""
    FF_s = smooth_5yr(FF)
    yrs_arr = np.array(years)
    pre_idx = np.where((yrs_arr >= 2000) & (yrs_arr <= 2006))[0]
    post_idx = np.where((yrs_arr >= 2010) & (yrs_arr <= 2018))[0]
    return np.nanmean(FF_s[post_idx]) - np.nanmean(FF_s[pre_idx])


def compute_residuals(data, mode, oh_13c_kie, n_iter=400, seed=42,
                      use_real_hemi_dD=False):
    """Compute per-year residual RMS from 3×3 system (goodness of fit).
    Lower = better fit = that KIE is more consistent with observations."""
    rng = np.random.default_rng(seed)
    n = data.n_years
    years = data.model_years
    CH4_NH, CH4_SH = data.CH4_NH, data.CH4_SH
    c13_NH, c13_SH = data.c13_NH, data.c13_SH
    c13_glob = data.c13_global
    
    tau_global = compute_lifetime(years, "varying", 9.0)
    tau_NH = tau_global * LIFETIME_RATIO_NH
    tau_SH = tau_global * LIFETIME_RATIO_SH
    
    kies = {
        'OH_13C': oh_13c_kie,
        'OH_D': 0.5*(1.294+1.327),
        'Cl_13C': 1.066, 'Cl_D': 1.520,
        'Strat_13C': 1.003, 'Strat_D': 1.050,
        'Soil_13C': 1.0201, 'Soil_D': 1.103,
    }
    K13_NH, KD_NH = compute_bulk_KIE(kies, SINK_FRACTIONS_NH)
    K13_SH, KD_SH = compute_bulk_KIE(kies, SINK_FRACTIONS_SH)
    a13_NH, aD_NH = 1.0/K13_NH, 1.0/KD_NH
    a13_SH, aD_SH = 1.0/K13_SH, 1.0/KD_SH
    
    W = np.diag([100.0, 1.0, 0.5])
    residuals = np.zeros((n, n_iter))
    nonphysical = np.zeros((n, n_iter))
    
    for k in range(n_iter):
        tau_ex = max(0.5, rng.normal(TAU_EX_MEAN, TAU_EX_STD))
        
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
            ex13_NH = (n13_SH-n13_NH)/tau_ex
            ex13_SH = (n13_NH-n13_SH)/tau_ex
            d13C_src_NH[j] = (n13_NH1-n13_NH + n13_NH*a13_NH/tau_NH[j] - ex13_NH)/S_NH[j]
            d13C_src_SH[j] = (n13_SH1-n13_SH + n13_SH*a13_SH/tau_SH[j] - ex13_SH)/S_SH[j]
        
        if use_real_hemi_dD:
            sigs = sample_source_signatures_hemi(rng, data, k, n)
            dD_NH_MC, dD_SH_MC = sample_atm_dD_hemi(data, k, n)
        else:
            sigs = sample_source_signatures(rng, data, k, n)
            dD_glob = sample_atm_dD(data, k, n)
            dD_NH_MC = dD_glob - DD_IH_OFFSET
            dD_SH_MC = dD_glob + DD_IH_OFFSET

        # δ¹³C source signatures: hemispheric when available
        if use_real_hemi_dD and 'ff_d13C_NH' in sigs:
            f13_bb_NH  = delta_to_fraction_d13C(sigs['bb_d13C_NH'])
            f13_ff_NH  = delta_to_fraction_d13C(sigs['ff_d13C_NH'])
            f13_mic_NH = delta_to_fraction_d13C(sigs['mic_d13C_NH'])
            f13_bb_SH  = delta_to_fraction_d13C(sigs['bb_d13C_SH'])
            f13_ff_SH  = delta_to_fraction_d13C(sigs['ff_d13C_SH'])
            f13_mic_SH = delta_to_fraction_d13C(sigs['mic_d13C_SH'])
        else:
            f13_bb  = delta_to_fraction_d13C(sigs['bb_d13C'])
            f13_ff  = delta_to_fraction_d13C(sigs['ff_d13C'])
            f13_mic = delta_to_fraction_d13C(sigs['mic_d13C'])
            f13_bb_NH = f13_bb_SH = f13_bb
            f13_ff_NH = f13_ff_SH = f13_ff
            f13_mic_NH = f13_mic_SH = f13_mic
        
        fD_NH_atm = delta_to_fraction_dD(dD_NH_MC)
        fD_SH_atm = delta_to_fraction_dD(dD_SH_MC)
        
        dD_src_NH = np.zeros(n); dD_src_SH = np.zeros(n)
        for j in range(n):
            nD_NH = fD_NH_atm[j]*CH4_NH[j]*PT_HEMI
            nD_NH1 = fD_NH_atm[j+1]*CH4_NH[j+1]*PT_HEMI
            nD_SH = fD_SH_atm[j]*CH4_SH[j]*PT_HEMI
            nD_SH1 = fD_SH_atm[j+1]*CH4_SH[j+1]*PT_HEMI
            exD_NH = (nD_SH-nD_NH)/tau_ex
            exD_SH = (nD_NH-nD_SH)/tau_ex
            dD_src_NH[j] = (nD_NH1-nD_NH + nD_NH*aD_NH/tau_NH[j] - exD_NH)/S_NH[j]
            dD_src_SH[j] = (nD_SH1-nD_SH + nD_SH*aD_SH/tau_SH[j] - exD_SH)/S_SH[j]
        
        if use_real_hemi_dD:
            fD_bb_NH = delta_to_fraction_dD(sigs['bb_dD_NH'])
            fD_ff_NH = delta_to_fraction_dD(sigs['ff_dD_NH'])
            fD_mic_NH = delta_to_fraction_dD(sigs['mic_dD_NH'])
            fD_bb_SH = delta_to_fraction_dD(sigs['bb_dD_SH'])
            fD_ff_SH = delta_to_fraction_dD(sigs['ff_dD_SH'])
            fD_mic_SH = delta_to_fraction_dD(sigs['mic_dD_SH'])
        else:
            fD_bb = delta_to_fraction_dD(sigs['bb_dD'])
            fD_ff = delta_to_fraction_dD(sigs['ff_dD'])
            fD_mic = delta_to_fraction_dD(sigs['mic_dD'])
            fD_bb_NH = fD_bb_SH = fD_bb
            fD_ff_NH = fD_ff_SH = fD_ff
            fD_mic_NH = fD_mic_SH = fD_mic
        
        for j in range(n):
            yr_resid = 0.0
            yr_np = 0
            for S, d13C_src, dD_src, f13_bb_h, f13_ff_h, f13_mic_h, fD_bb_h, fD_ff_h, fD_mic_h in [
                (S_NH[j], d13C_src_NH[j], dD_src_NH[j],
                 f13_bb_NH[j], f13_ff_NH[j], f13_mic_NH[j],
                 fD_bb_NH[j], fD_ff_NH[j], fD_mic_NH[j]),
                (S_SH[j], d13C_src_SH[j], dD_src_SH[j],
                 f13_bb_SH[j], f13_ff_SH[j], f13_mic_SH[j],
                 fD_bb_SH[j], fD_ff_SH[j], fD_mic_SH[j]),
            ]:
                A = np.array([
                    [1.0, 1.0, 1.0],
                    [f13_bb_h, f13_ff_h, f13_mic_h],
                    [fD_bb_h, fD_ff_h, fD_mic_h],
                ])
                B = np.array([S, S*d13C_src, S*dD_src])
                try:
                    res = lsq_linear(W@A, W@B, bounds=(0, S*1.5))
                    yr_resid += res.cost
                    # Check if solution hits bounds
                    if any(res.x <= 1e-6) or any(res.x >= S*1.5 - 1e-6):
                        yr_np += 1
                except:
                    yr_resid += 1e6
                    yr_np += 1
            residuals[j, k] = yr_resid
            nonphysical[j, k] = yr_np
    
    return residuals, nonphysical


def main():
    print("="*70)
    print("BASU 2022 COMPARISON v2: With real hemispheric δD")
    print("="*70)
    
    data = load_data(REPO_ROOT, two_box=True)
    N = 400
    
    results = {'basu_spread': BASU_KIE_SPREAD}
    
    configs = [
        ("d13C_only", "d13C_only", False),
        ("dual_offset", "dual", False),
        ("dual_real_hemi", "dual", True),
    ]
    
    for label, mode, use_real in configs:
        print(f"\n  Config: {label}")
        FF_S = run_at_fixed_kie(data, mode, 1.0039, N, 42, use_real)
        FF_C = run_at_fixed_kie(data, mode, 1.0054, N, 42, use_real)
        
        trend_S = np.median([trend_post_vs_pre2007(FF_S[:, k:k+1], data.model_years) 
                              for k in range(N)])
        trend_C = np.median([trend_post_vs_pre2007(FF_C[:, k:k+1], data.model_years) 
                              for k in range(N)])
        spread = abs(trend_C - trend_S)
        
        print(f"    Saueressig (1.0039): ΔFF = {trend_S:+.1f} Tg/yr")
        print(f"    Cantrell   (1.0054): ΔFF = {trend_C:+.1f} Tg/yr")
        print(f"    KIE spread:          {spread:.1f} Tg/yr")
        
        results[label] = {
            'trend_saueressig': float(trend_S),
            'trend_cantrell': float(trend_C),
            'kie_spread': float(spread),
        }
    
    # NEW: Residual analysis — which KIE fits better with real hemispheric δD?
    print(f"\n{'='*70}")
    print("RESIDUAL ANALYSIS: Which KIE is more consistent with observations?")
    print(f"{'='*70}")
    
    for label, use_real in [("dual_offset", False), ("dual_real_hemi", True)]:
        print(f"\n  Config: {label}")
        res_S, np_S = compute_residuals(data, "dual", 1.0039, N, 42, use_real)
        res_C, np_C = compute_residuals(data, "dual", 1.0054, N, 42, use_real)
        
        mean_res_S = np.nanmean(res_S)
        mean_res_C = np.nanmean(res_C)
        mean_np_S = np.nanmean(np_S)
        mean_np_C = np.nanmean(np_C)
        
        print(f"    Saueressig: mean residual = {mean_res_S:.2f}, mean bound-hits = {mean_np_S:.2f}")
        print(f"    Cantrell:   mean residual = {mean_res_C:.2f}, mean bound-hits = {mean_np_C:.2f}")
        
        if mean_res_C < mean_res_S:
            print(f"    → Cantrell fits {(mean_res_S/mean_res_C - 1)*100:.0f}% better")
            preferred = "Cantrell"
        else:
            print(f"    → Saueressig fits {(mean_res_C/mean_res_S - 1)*100:.0f}% better")
            preferred = "Saueressig"
        
        results[f'{label}_residuals'] = {
            'saueressig_mean_residual': float(mean_res_S),
            'cantrell_mean_residual': float(mean_res_C),
            'saueressig_mean_boundhits': float(mean_np_S),
            'cantrell_mean_boundhits': float(mean_np_C),
            'preferred': preferred,
        }
    
    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")
    print(f"  Basu 2022 (3D, δ¹³C-only):        KIE spread = {BASU_KIE_SPREAD:.1f} Tg/yr")
    for label in ['d13C_only', 'dual_offset', 'dual_real_hemi']:
        r = results[label]
        print(f"  Our 2-box ({label:<18s}):  KIE spread = {r['kie_spread']:.1f} Tg/yr")
    
    # Key question: does real δD change the preferred KIE?
    if 'dual_offset_residuals' in results and 'dual_real_hemi_residuals' in results:
        old_pref = results['dual_offset_residuals']['preferred']
        new_pref = results['dual_real_hemi_residuals']['preferred']
        print(f"\n  KIE preference (residual analysis):")
        print(f"    Old (offset δD):     {old_pref}")
        print(f"    New (real hemi δD):  {new_pref}")
        if old_pref != new_pref:
            print(f"    ⚡ FINDING: Real hemispheric δD REVERSES the KIE preference!")
        else:
            print(f"    ✓ Both configurations prefer {new_pref}")
    
    with open(OUT_DIR / "basu_comparison_v2.json", 'w') as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\n  Saved: {OUT_DIR}/basu_comparison_v2.json")


if __name__ == "__main__":
    main()
