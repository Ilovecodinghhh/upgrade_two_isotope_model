#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
phase3b_thanwerdas.py — Replicate Thanwerdas 2024 uncertainty specification
===========================================================================

Test: Is it the SPATIAL FRAMEWORK or the UNCERTAINTY SPECIFICATION that kills δD?

Thanwerdas 2024 (3D CTM, LMDz-SACS) found δD adds only "minor influence"
because they used very large prior uncertainties:
  - Mic(WET) δD: ±128‰  (40% of -320‰)
  - Mic(AGW) δD: ±93‰   (30% of -310‰)
  - FF δD: ±37‰         (20% of -183‰)
  - BB δD: ±70‰         (35% of -200‰)

We test:
  1. Our 2-box with OUR tight uncertainties → strong δD constraint
  2. Our 2-box with THANWERDAS uncertainties → δD should become useless
  3. If (2) shows δD useless: it's the uncertainty that kills it, not the model framework
"""

import sys
import json
from pathlib import Path

import numpy as np
from scipy.optimize import lsq_linear

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))
from common import (
    load_data, sample_KIE, compute_bulk_KIE, compute_lifetime,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    sample_source_signatures, sample_source_signatures_hemi,
    sample_atm_d13C, sample_atm_dD, sample_atm_dD_hemi,
    smooth_5yr,
    SINK_FRACTIONS_NH, SINK_FRACTIONS_SH,
    PT_HEMI, LIFETIME_RATIO_NH, LIFETIME_RATIO_SH,
    DD_IH_OFFSET, TAU_EX_MEAN, TAU_EX_STD,
    BB_NH_FRACTION, BB_SH_FRACTION,
)

OUT_DIR = Path(__file__).resolve().parent.parent / "results" / "phase3b_thanwerdas"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def apply_thanwerdas_uncertainties(sigs, rng):
    """
    Replace the MC-sampled δD source signatures with draws from
    Thanwerdas 2024's much larger prior distributions.
    
    Our baseline MC has: Mic σ≈8‰, FF σ≈0.7‰, BB σ≈7‰
    Thanwerdas has:       Mic σ≈110‰, FF σ≈37‰, BB σ≈70‰
    
    Handles both global and hemispheric keys.
    """
    sigs_new = dict(sigs)
    
    # Mean values from our data, but inflate to Thanwerdas spread
    mic_dD_mean = -305.0  # typical microbial δD
    ff_dD_mean = -193.0   # typical FF δD
    bb_dD_mean = -227.0   # typical BB δD
    
    # Thanwerdas uncertainties (1σ, approximate from their Table 2)
    mic_dD_sigma = 110.0  # average of WET(128) and AGW(93)
    ff_dD_sigma = 37.0
    bb_dD_sigma = 70.0
    
    # Draw single perturbations (coherent across time)
    mic_pert = rng.normal(0, mic_dD_sigma)
    ff_pert = rng.normal(0, ff_dD_sigma)
    bb_pert = rng.normal(0, bb_dD_sigma)
    
    # Apply to global keys
    n = len(sigs['mic_dD'])
    sigs_new['mic_dD'] = np.full(n, mic_dD_mean) + mic_pert
    sigs_new['ff_dD'] = np.full(n, ff_dD_mean) + ff_pert
    sigs_new['bb_dD'] = np.full(n, bb_dD_mean) + bb_pert
    
    # Apply same perturbation to hemispheric keys if present
    for hemi in ('NH', 'SH'):
        for src, mean_val, pert in [('mic_dD', mic_dD_mean, mic_pert),
                                     ('ff_dD', ff_dD_mean, ff_pert),
                                     ('bb_dD', bb_dD_mean, bb_pert)]:
            key = f'{src}_{hemi}'
            if key in sigs:
                sigs_new[key] = np.full(n, mean_val) + pert
    
    return sigs_new


def run_twobox(data, mode="dual", uncertainty="ours", n_iter=500, seed=42):
    """Run 2-box model with specified isotope mode and uncertainty."""
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
    
    FF_G = np.zeros((n, n_iter))
    Mic_G = np.zeros((n, n_iter))
    
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
        
        sigs = sample_source_signatures_hemi(rng, data, k, n)
        
        # Apply Thanwerdas uncertainties if requested
        if uncertainty == "thanwerdas":
            sigs = apply_thanwerdas_uncertainties(sigs, rng)
        
        f13_bb = delta_to_fraction_d13C(sigs['bb_d13C'])
        f13_ff = delta_to_fraction_d13C(sigs['ff_d13C'])
        f13_mic = delta_to_fraction_d13C(sigs['mic_d13C'])
        
        if mode == "dual":
            dD_NH_MC, dD_SH_MC = sample_atm_dD_hemi(data, k, n)
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
            
            # Hemispheric δD source signatures
            fD_bb_NH = delta_to_fraction_dD(sigs['bb_dD_NH'])
            fD_ff_NH = delta_to_fraction_dD(sigs['ff_dD_NH'])
            fD_mic_NH = delta_to_fraction_dD(sigs['mic_dD_NH'])
            fD_bb_SH = delta_to_fraction_dD(sigs['bb_dD_SH'])
            fD_ff_SH = delta_to_fraction_dD(sigs['ff_dD_SH'])
            fD_mic_SH = delta_to_fraction_dD(sigs['mic_dD_SH'])
            
            for j in range(n):
                A_nh = np.array([
                    [1.0, 1.0, 1.0],
                    [f13_bb[j], f13_ff[j], f13_mic[j]],
                    [fD_bb_NH[j], fD_ff_NH[j], fD_mic_NH[j]],
                ])
                B_nh = np.array([S_NH[j], S_NH[j]*d13C_src_NH[j], S_NH[j]*dD_src_NH[j]])
                try:
                    res = lsq_linear(W_NH@A_nh, W_NH@B_nh, bounds=(0, S_NH[j]*1.5))
                    FF_G[j,k] += res.x[1]
                    Mic_G[j,k] += res.x[2]
                except: pass
                
                A_sh = np.array([
                    [1.0, 1.0, 1.0],
                    [f13_bb[j], f13_ff[j], f13_mic[j]],
                    [fD_bb_SH[j], fD_ff_SH[j], fD_mic_SH[j]],
                ])
                B_sh = np.array([S_SH[j], S_SH[j]*d13C_src_SH[j], S_SH[j]*dD_src_SH[j]])
                try:
                    res = lsq_linear(W_SH@A_sh, W_SH@B_sh, bounds=(0, S_SH[j]*1.5))
                    FF_G[j,k] += res.x[1]
                    Mic_G[j,k] += res.x[2]
                except: pass
        
        else:  # d13C_only
            for j in range(n):
                denom = f13_ff[j] - f13_mic[j]
                if abs(denom) < 1e-15: continue
                # NH
                S_rem = S_NH[j] - BB_hemi_NH
                rhs = S_NH[j]*d13C_src_NH[j] - BB_hemi_NH*f13_bb[j]
                FF_G[j,k] += (rhs - S_rem*f13_mic[j]) / denom
                Mic_G[j,k] += S_rem - FF_G[j,k]  # wrong — need to track
                # Fix: redo properly
                ff_nh = (rhs - S_rem*f13_mic[j]) / denom
                mic_nh = S_rem - ff_nh
                # SH
                S_rem = S_SH[j] - BB_hemi_SH
                rhs = S_SH[j]*d13C_src_SH[j] - BB_hemi_SH*f13_bb[j]
                ff_sh = (rhs - S_rem*f13_mic[j]) / denom
                mic_sh = S_rem - ff_sh
                FF_G[j,k] = ff_nh + ff_sh
                Mic_G[j,k] = mic_nh + mic_sh
    
    return years, FF_G, Mic_G


def compute_ci(arr):
    """Compute 90% CI width from smoothed MC array."""
    s = smooth_5yr(arr)
    ci_5 = np.nanpercentile(s[8:], 5, axis=1)
    ci_95 = np.nanpercentile(s[8:], 95, axis=1)
    return float(np.nanmean(ci_95 - ci_5))


def main():
    print("="*70)
    print("PHASE 3b: THANWERDAS REPLICATION")
    print("="*70)
    
    data = load_data(REPO_ROOT, two_box=True)
    N = 500
    
    configs = [
        ("2-box, dual, OUR uncertainties", "dual", "ours"),
        ("2-box, dual, THANWERDAS uncertainties", "dual", "thanwerdas"),
        ("2-box, δ¹³C-only (reference)", "d13C_only", "ours"),
    ]
    
    results = {}
    for label, mode, uncert in configs:
        print(f"\n  Running: {label}...")
        yrs, FF, Mic = run_twobox(data, mode, uncert, N, seed=42)
        ci = compute_ci(FF)
        results[label] = {
            'FF_CI_width': ci,
            'FF_mean': float(np.nanmean(smooth_5yr(FF)[8:])),
        }
        print(f"    FF 90% CI width: {ci:.1f} Tg/yr")
    
    # Analysis
    ci_ours = results["2-box, dual, OUR uncertainties"]['FF_CI_width']
    ci_than = results["2-box, dual, THANWERDAS uncertainties"]['FF_CI_width']
    ci_ref = results["2-box, δ¹³C-only (reference)"]['FF_CI_width']
    
    improve_ours = (ci_ref - ci_ours) / ci_ref * 100
    improve_than = (ci_ref - ci_than) / ci_ref * 100
    
    print(f"\n{'='*70}")
    print("DIAGNOSIS: What kills δD?")
    print(f"{'='*70}")
    print(f"\n  Baseline (δ¹³C-only):    CI = {ci_ref:.1f} Tg/yr")
    print(f"  + δD (our σ ≈ 8‰):       CI = {ci_ours:.1f} Tg/yr  → {improve_ours:+.1f}% improvement")
    print(f"  + δD (Thanwerdas σ≈110‰): CI = {ci_than:.1f} Tg/yr  → {improve_than:+.1f}% improvement")
    
    if improve_than < 10:
        print(f"\n  ✓ CONFIRMED: It's the UNCERTAINTY SPECIFICATION that kills δD,")
        print(f"    not the model framework. Even our simple 2-box model reproduces")
        print(f"    Thanwerdas's finding when given their source-signature priors.")
        print(f"\n  The take-home: δD IS informative with current measurement precision")
        print(f"    (σ ≈ 8‰) but NOT with the overly conservative priors used in 3D CTM.")
    else:
        print(f"\n  ✗ Unexpected: δD still helps under Thanwerdas uncertainties?!")
    
    summary = {
        'our_uncertainties': {'CI': ci_ours, 'improvement_pct': improve_ours},
        'thanwerdas_uncertainties': {'CI': ci_than, 'improvement_pct': improve_than},
        'd13C_reference': {'CI': ci_ref},
        'conclusion': 'uncertainty_kills_dD' if improve_than < 10 else 'framework_matters',
    }
    
    with open(OUT_DIR / "thanwerdas_comparison.json", 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n  Saved: {OUT_DIR}/thanwerdas_comparison.json")


if __name__ == "__main__":
    main()
