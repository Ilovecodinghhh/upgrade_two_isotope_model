#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
compare_basu2022.py — Compare our KIE spread to Basu 2022 (ACP)
================================================================

Basu et al. (2022, ACP) ran two TM5-4DVAR inversions of CH₄+δ¹³C
that differed only in OH KIE choice (Saueressig vs Cantrell).
Their post-2007 FF emission attributions differed by ~13 Tg/yr.

Our model reproduces this when running δ¹³C-only.
With dual isotopes, the spread should collapse.

References for Basu 2022 numbers (extracted from MASTER_DATA_INVENTORY):
  - Saueressig posterior:  FF increase ≈ +12 Tg/yr (post-2007 vs pre-2007)
  - Cantrell posterior:    FF increase ≈ +25 Tg/yr
  - KIE-driven ambiguity:  ~13 Tg/yr
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
    sample_source_signatures, sample_atm_d13C, sample_atm_dD,
    smooth_5yr,
    SINK_FRACTIONS_NH, SINK_FRACTIONS_SH,
    PT_HEMI, LIFETIME_RATIO_NH, LIFETIME_RATIO_SH,
    DD_IH_OFFSET, TAU_EX_MEAN, TAU_EX_STD,
    BB_NH_FRACTION, BB_SH_FRACTION,
)

OUT_DIR = Path(__file__).resolve().parent.parent / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Basu 2022 reported values (Saueressig vs Cantrell experiment)
BASU_SAUERESSIG_FF_TREND = 12.0  # Tg/yr (approximate, post-2007)
BASU_CANTRELL_FF_TREND = 25.0
BASU_KIE_SPREAD = abs(BASU_CANTRELL_FF_TREND - BASU_SAUERESSIG_FF_TREND)  # ≈13


def run_at_fixed_kie(data, mode, oh_13c_kie, n_iter=400, seed=42):
    """Run 2-box at a FIXED OH-13C KIE value. Return FF time series MC."""
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
    
    # Fixed KIE config (only OH_13C varies between Basu's two runs)
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
        
        sigs = sample_source_signatures(rng, data, k, n)
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
                for S, d13C_src, dD_src in [
                    (S_NH[j], d13C_src_NH[j], dD_src_NH[j]),
                    (S_SH[j], d13C_src_SH[j], dD_src_SH[j]),
                ]:
                    B = np.array([S, S*d13C_src, S*dD_src])
                    try:
                        res = lsq_linear(W@A, W@B, bounds=(0, S*1.5))
                        FF_G[j,k] += res.x[1]
                    except: pass
        else:
            for j in range(n):
                denom = f13_ff[j] - f13_mic[j]
                if abs(denom) < 1e-15: continue
                for S, d13C_src, BB in [
                    (S_NH[j], d13C_src_NH[j], BB_NH),
                    (S_SH[j], d13C_src_SH[j], BB_SH),
                ]:
                    S_rem = S - BB
                    rhs = S*d13C_src - BB*f13_bb[j]
                    ff = (rhs - S_rem*f13_mic[j]) / denom
                    FF_G[j,k] += max(0, min(ff, S*1.5))
    
    return FF_G


def trend_post_vs_pre2007(FF, years):
    """Compute FF increase: mean(2010-2018) - mean(2000-2006)."""
    FF_s = smooth_5yr(FF)
    yrs_arr = np.array(years)
    pre_idx = np.where((yrs_arr >= 2000) & (yrs_arr <= 2006))[0]
    post_idx = np.where((yrs_arr >= 2010) & (yrs_arr <= 2018))[0]
    pre = np.nanmean(FF_s[pre_idx])
    post = np.nanmean(FF_s[post_idx])
    return post - pre


def main():
    print("="*70)
    print("BASU 2022 COMPARISON: KIE spread in FF emission attribution")
    print("="*70)
    print(f"\nBasu 2022 (ACP) reports:")
    print(f"  Saueressig (KIE=1.0039): ΔFF post-2007 ≈ +{BASU_SAUERESSIG_FF_TREND} Tg/yr")
    print(f"  Cantrell    (KIE=1.0054): ΔFF post-2007 ≈ +{BASU_CANTRELL_FF_TREND} Tg/yr")
    print(f"  Spread: {BASU_KIE_SPREAD} Tg/yr")
    
    data = load_data(REPO_ROOT, two_box=True)
    N = 400
    
    results = {'basu_spread': BASU_KIE_SPREAD}
    
    for mode in ['d13C_only', 'dual']:
        print(f"\n  Running our 2-box ({mode})...")
        FF_S = run_at_fixed_kie(data, mode, 1.0039, n_iter=N, seed=42)
        FF_C = run_at_fixed_kie(data, mode, 1.0054, n_iter=N, seed=42)
        
        # Median trend across MC
        trend_S = np.median([trend_post_vs_pre2007(FF_S[:, k:k+1], data.model_years) 
                              for k in range(N)])
        trend_C = np.median([trend_post_vs_pre2007(FF_C[:, k:k+1], data.model_years) 
                              for k in range(N)])
        spread = abs(trend_C - trend_S)
        
        print(f"    Saueressig: ΔFF = {trend_S:+.1f} Tg/yr")
        print(f"    Cantrell:   ΔFF = {trend_C:+.1f} Tg/yr")
        print(f"    Spread:     {spread:.1f} Tg/yr")
        
        results[f'our_{mode}'] = {
            'trend_saueressig': float(trend_S),
            'trend_cantrell': float(trend_C),
            'kie_spread': float(spread),
        }
    
    # Comparison
    print(f"\n{'='*70}")
    print("KEY RESULT")
    print(f"{'='*70}")
    print(f"  Basu 2022 (δ¹³C-only, 3D):        KIE spread = {BASU_KIE_SPREAD:.1f} Tg/yr")
    print(f"  Our 2-box (δ¹³C-only):           KIE spread = {results['our_d13C_only']['kie_spread']:.1f} Tg/yr")
    print(f"  Our 2-box (dual-isotope):        KIE spread = {results['our_dual']['kie_spread']:.1f} Tg/yr")
    
    reduction = (1 - results['our_dual']['kie_spread'] /
                  max(results['our_d13C_only']['kie_spread'], 0.1)) * 100
    print(f"\n  → Dual isotopes reduce KIE spread by {reduction:.0f}%")
    print(f"\n  IMPLICATION: If Basu 2022 had included δD, their KIE-driven uncertainty")
    print(f"  would have collapsed from ~13 Tg/yr to ~{BASU_KIE_SPREAD*(1-reduction/100):.1f} Tg/yr.")
    
    with open(OUT_DIR / "basu_comparison.json", 'w') as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\n  Saved: {OUT_DIR}/basu_comparison.json")


if __name__ == "__main__":
    main()
