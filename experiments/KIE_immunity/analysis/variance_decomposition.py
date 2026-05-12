#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
variance_decomposition.py — Decompose FF uncertainty into its sources
======================================================================
v2: Uses REAL hemispheric δD observations and source signatures
    instead of the DD_IH_OFFSET = ±6‰ hack.

Compares three configurations:
  1. δ¹³C-only (2-box, BB-fixed)
  2. Dual-isotope with OLD δD hemispheric data (±6‰ offset hack)
  3. Dual-isotope with REAL hemispheric δD (station-level MC + gridded sigs)

For each config, decomposes FF variance by:
  - Fixing KIE → measure KIE contribution
  - Fixing source signatures → measure sig contribution
  - Fixing lifetime → measure τ contribution
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

OUT_DIR = Path(__file__).resolve().parent.parent / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def run_2box(data, mode, n_iter, seed,
             fix_kie=False, fix_sigs=False, fix_tau=False,
             use_real_hemi_dD=False):
    """
    Run 2-box model with selectable randomness sources.
    
    mode: "d13C_only" or "dual"
    use_real_hemi_dD: if True and mode=="dual", use real hemispheric δD
                       instead of global ± offset
    """
    rng = np.random.default_rng(seed)
    n = data.n_years
    years = data.model_years
    CH4_NH, CH4_SH = data.CH4_NH, data.CH4_SH
    c13_NH, c13_SH = data.c13_NH, data.c13_SH
    c13_glob = data.c13_global
    
    # Lifetime
    if fix_tau:
        tau_global = compute_lifetime(years, "fixed", 9.0)
    else:
        tau_global = compute_lifetime(years, "varying", 9.0)
    tau_NH = tau_global * LIFETIME_RATIO_NH
    tau_SH = tau_global * LIFETIME_RATIO_SH
    
    BB_NH = data.BB_global_mean * BB_NH_FRACTION
    BB_SH = data.BB_global_mean * BB_SH_FRACTION
    
    FF_G = np.zeros((n, n_iter))
    W = np.diag([100.0, 1.0, 0.5]) if mode == "dual" else None
    
    # Pre-compute fixed values if needed
    if fix_kie:
        kie_fixed = {'OH_13C': 0.5*(1.0039+1.0054),
                     'OH_D': 0.5*(1.294+1.327),
                     'Cl_13C': 1.066, 'Cl_D': 1.520,
                     'Strat_13C': 1.003, 'Strat_D': 1.050,
                     'Soil_13C': 1.0201, 'Soil_D': 1.103}
    if fix_sigs:
        rng_tmp = np.random.default_rng(0)
        if use_real_hemi_dD:
            sigs_fixed = sample_source_signatures_hemi(rng_tmp, data, 0, n)
        else:
            sigs_fixed = sample_source_signatures(rng_tmp, data, 0, n)
    
    for k in range(n_iter):
        tau_ex = TAU_EX_MEAN if fix_tau else max(0.5, rng.normal(TAU_EX_MEAN, TAU_EX_STD))
        
        if fix_kie:
            kies = kie_fixed
        else:
            kies = sample_KIE(rng, "sampled")
        
        K13_NH, KD_NH = compute_bulk_KIE(kies, SINK_FRACTIONS_NH)
        K13_SH, KD_SH = compute_bulk_KIE(kies, SINK_FRACTIONS_SH)
        a13_NH, aD_NH = 1.0/K13_NH, 1.0/KD_NH
        a13_SH, aD_SH = 1.0/K13_SH, 1.0/KD_SH
        
        # Total source strength
        S_NH = np.zeros(n); S_SH = np.zeros(n)
        for i in range(n):
            M_NH = CH4_NH[i]*PT_HEMI; M_NH1 = CH4_NH[i+1]*PT_HEMI
            M_SH = CH4_SH[i]*PT_HEMI; M_SH1 = CH4_SH[i+1]*PT_HEMI
            ex_NH = (M_SH-M_NH)/tau_ex; ex_SH = (M_NH-M_SH)/tau_ex
            S_NH[i] = (M_NH1-M_NH) + M_NH/tau_NH[i] - ex_NH
            S_SH[i] = (M_SH1-M_SH) + M_SH/tau_SH[i] - ex_SH
        
        # δ¹³C source composition
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
        
        # Source signatures
        if fix_sigs:
            sigs = sigs_fixed
        elif use_real_hemi_dD:
            sigs = sample_source_signatures_hemi(rng, data, k, n)
        else:
            sigs = sample_source_signatures(rng, data, k, n)
        
        f13_bb = delta_to_fraction_d13C(sigs['bb_d13C'])
        f13_ff = delta_to_fraction_d13C(sigs['ff_d13C'])
        f13_mic = delta_to_fraction_d13C(sigs['mic_d13C'])
        
        if mode == "dual":
            # δD atmospheric observations
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
            
            # δD source signatures — hemispheric or global
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
                for hemi_idx, (S, d13C_src, dD_src, fD_bb_h, fD_ff_h, fD_mic_h) in enumerate([
                    (S_NH[j], d13C_src_NH[j], dD_src_NH[j], fD_bb_NH[j], fD_ff_NH[j], fD_mic_NH[j]),
                    (S_SH[j], d13C_src_SH[j], dD_src_SH[j], fD_bb_SH[j], fD_ff_SH[j], fD_mic_SH[j]),
                ]):
                    A = np.array([
                        [1.0, 1.0, 1.0],
                        [f13_bb[j], f13_ff[j], f13_mic[j]],
                        [fD_bb_h, fD_ff_h, fD_mic_h],
                    ])
                    B = np.array([S, S*d13C_src, S*dD_src])
                    try:
                        res = lsq_linear(W@A, W@B, bounds=(0, S*1.5))
                        FF_G[j,k] += res.x[1]
                    except: pass
        else:
            # δ¹³C-only: BB fixed, solve for FF and Mic
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


def compute_var(data, mode, n_iter, seed, use_real_hemi_dD=False, **flags):
    """Run and return mean post-2007 variance."""
    FF = run_2box(data, mode, n_iter, seed,
                  use_real_hemi_dD=use_real_hemi_dD, **flags)
    FF_s = smooth_5yr(FF)
    return float(np.nanmean(np.nanvar(FF_s[8:], axis=1)))


def main():
    print("="*70)
    print("VARIANCE DECOMPOSITION v2: Real hemispheric δD")
    print("="*70)
    
    data = load_data(REPO_ROOT, two_box=True)
    N = 400
    SEED = 42
    
    results = {}
    
    configs = [
        ("d13C_only", False),
        ("dual_offset", False),     # Old: global δD ± 6‰
        ("dual_real_hemi", True),   # New: real hemispheric δD
    ]
    
    for label, use_real in configs:
        mode = "d13C_only" if label == "d13C_only" else "dual"
        print(f"\n  Config: {label} (mode={mode}, real_hemi={use_real})")
        
        var_total = compute_var(data, mode, N, SEED, use_real_hemi_dD=use_real)
        var_no_kie = compute_var(data, mode, N, SEED, use_real_hemi_dD=use_real, fix_kie=True)
        var_no_sigs = compute_var(data, mode, N, SEED, use_real_hemi_dD=use_real, fix_sigs=True)
        var_no_tau = compute_var(data, mode, N, SEED, use_real_hemi_dD=use_real, fix_tau=True)
        
        kie_contrib = max(0, var_total - var_no_kie)
        sigs_contrib = max(0, var_total - var_no_sigs)
        tau_contrib = max(0, var_total - var_no_tau)
        residual = max(0, var_total - kie_contrib - sigs_contrib - tau_contrib)
        
        total = kie_contrib + sigs_contrib + tau_contrib + residual
        if total < 1e-9:
            total = 1.0
        
        sigma = np.sqrt(var_total)
        print(f"    Total variance:    {var_total:8.1f} (Tg/yr)² → σ = {sigma:.1f} Tg/yr")
        print(f"    Fix KIE  → var =   {var_no_kie:8.1f}  (KIE contrib = {kie_contrib/total*100:.1f}%)")
        print(f"    Fix sigs → var =   {var_no_sigs:8.1f}  (Sig contrib = {sigs_contrib/total*100:.1f}%)")
        print(f"    Fix tau  → var =   {var_no_tau:8.1f}  (Tau contrib = {tau_contrib/total*100:.1f}%)")
        print(f"    Residual (atm + interaction): {residual/total*100:.1f}%")
        
        results[label] = {
            'var_total': var_total,
            'sigma': sigma,
            'kie_pct': kie_contrib/total*100,
            'sigs_pct': sigs_contrib/total*100,
            'tau_pct': tau_contrib/total*100,
            'residual_pct': residual/total*100,
        }
    
    # Comparison table
    print(f"\n{'='*70}")
    print("COMPARISON TABLE")
    print(f"{'='*70}")
    print(f"{'Config':<20} {'σ(FF)':>8} {'KIE%':>7} {'Sig%':>7} {'τ%':>7} {'Resid%':>7}")
    print("-"*60)
    for label in ['d13C_only', 'dual_offset', 'dual_real_hemi']:
        r = results[label]
        print(f"{label:<20} {r['sigma']:>7.1f}  {r['kie_pct']:>6.1f} {r['sigs_pct']:>7.1f} {r['tau_pct']:>6.1f} {r['residual_pct']:>7.1f}")
    
    # Key scientific findings
    print(f"\n{'='*70}")
    print("KEY FINDINGS")
    print(f"{'='*70}")
    
    d13c = results['d13C_only']
    old = results['dual_offset']
    new = results['dual_real_hemi']
    
    print(f"\n1. δ¹³C-only → dual (offset):  σ reduction = {d13c['sigma']:.1f} → {old['sigma']:.1f} Tg/yr ({(1-old['sigma']/d13c['sigma'])*100:.0f}%)")
    print(f"2. δ¹³C-only → dual (real):    σ reduction = {d13c['sigma']:.1f} → {new['sigma']:.1f} Tg/yr ({(1-new['sigma']/d13c['sigma'])*100:.0f}%)")
    print(f"3. dual (offset) → dual (real): σ change = {old['sigma']:.1f} → {new['sigma']:.1f} Tg/yr ({(new['sigma']/old['sigma']-1)*100:+.0f}%)")
    
    kie_reduction_offset = d13c['kie_pct'] - old['kie_pct']
    kie_reduction_real = d13c['kie_pct'] - new['kie_pct']
    print(f"\n4. KIE contribution: δ¹³C-only = {d13c['kie_pct']:.1f}% → offset = {old['kie_pct']:.1f}% → real = {new['kie_pct']:.1f}%")
    
    if new['sigma'] < old['sigma']:
        print(f"\n✓ REAL hemispheric δD FURTHER reduces FF uncertainty!")
        print(f"  This means the ±6‰ offset was too crude — real NH/SH δD")
        print(f"  gradient (~15‰) provides stronger constraint.")
    elif new['sigma'] > old['sigma'] * 1.1:
        print(f"\n⚠ REAL hemispheric δD INCREASES FF uncertainty!")
        print(f"  This means the ±6‰ offset was artificially constraining —")
        print(f"  real NH/SH δD data introduce realistic scatter that the")
        print(f"  crude offset suppressed. This is actually MORE HONEST.")
    else:
        print(f"\n≈ Real vs offset δD yield similar variance")
        print(f"  The ±6‰ offset was a reasonable first-order approximation.")
    
    with open(OUT_DIR / "variance_decomposition_v2.json", 'w') as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\n  Saved: {OUT_DIR}/variance_decomposition_v2.json")


if __name__ == "__main__":
    main()
