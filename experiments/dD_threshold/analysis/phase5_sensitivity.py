#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
phase5_sensitivity.py — KIE and Lifetime Sensitivity at the Threshold
======================================================================

Tests whether the δD threshold (~25‰ for Mic δD) is robust to:
  A. KIE choice: Saueressig vs. Cantrell vs. sampled
  B. Lifetime assumption: fixed 9.0yr vs. varying vs. short (8.5yr)
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

OUT_DIR = Path(__file__).resolve().parent.parent / "results" / "phase5_sensitivity"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def run_threshold_sweep(data, multipliers, kie_mode="sampled", 
                         lifetime_mode="varying", tau_fixed=9.0, 
                         n_iter=300, seed=42):
    """
    Run a threshold sweep for a given KIE/lifetime configuration.
    Returns dict of multiplier → CI width.
    """
    results = {}
    
    for mult in multipliers:
        rng = np.random.default_rng(seed)
        n = data.n_years
        years = data.model_years
        
        CH4_NH, CH4_SH = data.CH4_NH, data.CH4_SH
        c13_NH, c13_SH = data.c13_NH, data.c13_SH
        c13_glob = data.c13_global
        
        tau_global = compute_lifetime(years, lifetime_mode, tau_fixed)
        tau_NH = tau_global * LIFETIME_RATIO_NH
        tau_SH = tau_global * LIFETIME_RATIO_SH
        
        BB_hemi_NH = data.BB_global_mean * BB_NH_FRACTION
        BB_hemi_SH = data.BB_global_mean * BB_SH_FRACTION
        
        FF_G = np.zeros((n, n_iter))
        W_NH = np.diag([100.0, 1.0, 0.5])
        W_SH = np.diag([200.0, 1.0, 0.5])
        
        for k in range(n_iter):
            tau_ex = max(0.5, rng.normal(TAU_EX_MEAN, TAU_EX_STD))
            kies = sample_KIE(rng, kie_mode)
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
            
            # Inflate δD
            if mult > 1.0:
                extra_mic = rng.normal() * 8.25 * (mult - 1)
                extra_ff = rng.normal() * 0.70 * (mult - 1)
                extra_bb = rng.normal() * 7.09 * (mult - 1)
                sigs = dict(sigs)
                sigs['mic_dD'] = sigs['mic_dD'] + extra_mic
                sigs['ff_dD'] = sigs['ff_dD'] + extra_ff
                sigs['bb_dD'] = sigs['bb_dD'] + extra_bb
                # Also inflate hemispheric keys
                for hemi in ('NH', 'SH'):
                    for src, extra in [('mic_dD', extra_mic), ('ff_dD', extra_ff), ('bb_dD', extra_bb)]:
                        key = f'{src}_{hemi}'
                        if key in sigs:
                            sigs[key] = sigs[key] + extra
            
            f13_bb = delta_to_fraction_d13C(sigs['bb_d13C'])
            f13_ff = delta_to_fraction_d13C(sigs['ff_d13C'])
            f13_mic = delta_to_fraction_d13C(sigs['mic_d13C'])
            
            dD_NH_MC, dD_SH_MC = sample_atm_dD_hemi(data, k, target_length=n)
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
                except: pass
        
        # Compute CI width
        FF_s = smooth_5yr(FF_G)
        ci = np.nanpercentile(FF_s[8:], 95, axis=1) - np.nanpercentile(FF_s[8:], 5, axis=1)
        results[mult] = float(np.nanmean(ci))
    
    return results


def main():
    print("="*70)
    print("PHASE 5: SENSITIVITY ANALYSIS")
    print("="*70)
    
    data = load_data(REPO_ROOT, two_box=True)
    multipliers = [1.0, 2.0, 3.0, 5.0, 8.0]
    N = 300  # Faster for sensitivity
    
    # Also get baseline (d13C-only) for reference
    # (reuse from Phase 3 — known value ~101 Tg/yr)
    d13C_ref = 101.3
    
    all_results = {}
    
    # === A. KIE Sensitivity ===
    print("\n" + "─"*50)
    print("A. KIE SENSITIVITY")
    print("─"*50)
    
    kie_configs = [
        ("saueressig", "saueressig"),
        ("cantrell", "cantrell"),
        ("sampled", "sampled"),
    ]
    
    for label, kie_mode in kie_configs:
        print(f"\n  KIE = {label}...")
        res = run_threshold_sweep(data, multipliers, kie_mode=kie_mode, n_iter=N)
        all_results[f'KIE_{label}'] = res
        # Find threshold
        for m in multipliers:
            improvement = (d13C_ref - res[m]) / d13C_ref * 100
            if improvement < 10:
                print(f"    Threshold: mult={m:.1f}× (σ≈{8.25*m:.0f}‰)")
                break
        else:
            print(f"    Threshold: > {multipliers[-1]}× (σ≈{8.25*multipliers[-1]:.0f}‰)")
    
    # === B. Lifetime Sensitivity ===
    print("\n" + "─"*50)
    print("B. LIFETIME SENSITIVITY")
    print("─"*50)
    
    lifetime_configs = [
        ("fixed_9.0", "fixed", 9.0),
        ("varying", "varying", 9.0),
        ("fixed_8.5", "fixed", 8.5),
    ]
    
    for label, mode, tau in lifetime_configs:
        print(f"\n  Lifetime = {label}...")
        res = run_threshold_sweep(data, multipliers, lifetime_mode=mode, 
                                   tau_fixed=tau, n_iter=N)
        all_results[f'tau_{label}'] = res
        for m in multipliers:
            improvement = (d13C_ref - res[m]) / d13C_ref * 100
            if improvement < 10:
                print(f"    Threshold: mult={m:.1f}× (σ≈{8.25*m:.0f}‰)")
                break
        else:
            print(f"    Threshold: > {multipliers[-1]}× (σ≈{8.25*multipliers[-1]:.0f}‰)")
    
    # === Summary ===
    print(f"\n{'='*70}")
    print("SENSITIVITY SUMMARY")
    print(f"{'='*70}")
    
    print(f"\n{'Config':<20} {'mult=1':<10} {'mult=2':<10} {'mult=3':<10} {'mult=5':<10} {'mult=8':<10}")
    for key, res in all_results.items():
        vals = [f"{res[m]:.0f}" for m in multipliers]
        print(f"{key:<20} {'  '.join(vals)}")
    
    print(f"\n  d13C-only reference: {d13C_ref:.1f} Tg/yr")
    print(f"  Threshold criterion: improvement < 10%")
    
    # Save
    with open(OUT_DIR / "sensitivity_results.json", 'w') as f:
        json.dump(all_results, f, indent=2, default=float)
    
    print(f"\n  Saved: {OUT_DIR}/sensitivity_results.json")


if __name__ == "__main__":
    main()
