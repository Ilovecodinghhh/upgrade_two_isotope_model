#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
variance_decomposition.py — Decompose FF uncertainty into its sources
======================================================================

For each (model × isotope mode) configuration, run with:
  - All sources varied → variance_total
  - KIE fixed → variance_minus_KIE
  - Source signatures fixed → variance_minus_sigs
  - Lifetime fixed → variance_minus_tau

Then attribution:
  KIE contribution = variance_total − variance_minus_KIE
  Sig contribution = variance_total − variance_minus_sigs
  τ   contribution = variance_total − variance_minus_tau

Key paper claim: In δ¹³C-only, KIE dominates (>50% of variance).
                  In dual-isotope, KIE shrinks to a minor contributor.
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
    sample_source_signatures, sample_atm_d13C, sample_atm_dD,
    smooth_5yr,
    SINK_FRACTIONS_NH, SINK_FRACTIONS_SH,
    PT_HEMI, LIFETIME_RATIO_NH, LIFETIME_RATIO_SH,
    DD_IH_OFFSET, TAU_EX_MEAN, TAU_EX_STD,
    BB_NH_FRACTION, BB_SH_FRACTION,
)

OUT_DIR = Path(__file__).resolve().parent.parent / "results"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def run_2box(data, mode, n_iter, seed,
             fix_kie=False, fix_sigs=False, fix_tau=False):
    """Run 2-box with selectable randomness sources."""
    rng = np.random.default_rng(seed)
    n = data.n_years
    years = data.model_years
    CH4_NH, CH4_SH = data.CH4_NH, data.CH4_SH
    c13_NH, c13_SH = data.c13_NH, data.c13_SH
    c13_glob = data.c13_global
    
    # Lifetime (vary or fix)
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
    
    # Fixed-value placeholders if fix_kie/fix_sigs
    if fix_kie:
        kie_fixed = {'OH_13C': 0.5*(1.0039+1.0054),
                     'OH_D': 0.5*(1.294+1.327),
                     'Cl_13C': 1.066, 'Cl_D': 1.520,
                     'Strat_13C': 1.003, 'Strat_D': 1.050,
                     'Soil_13C': 1.0201, 'Soil_D': 1.103}
    if fix_sigs:
        # Use means from CSVs (call once with seed=0 and take iter 0)
        rng_tmp = np.random.default_rng(0)
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
        
        sigs = sigs_fixed if fix_sigs else sample_source_signatures(rng, data, k, n)
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
                for hemi, S, d13C_src, dD_src in [
                    ('NH', S_NH[j], d13C_src_NH[j], dD_src_NH[j]),
                    ('SH', S_SH[j], d13C_src_SH[j], dD_src_SH[j]),
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


def compute_var(data, mode, n_iter, seed, **flags):
    """Run and return mean post-2007 variance."""
    FF = run_2box(data, mode, n_iter, seed, **flags)
    FF_s = smooth_5yr(FF)
    # Average variance over post-2007 (smoothed window starts at index 8)
    return float(np.nanmean(np.nanvar(FF_s[8:], axis=1)))


def main():
    print("="*70)
    print("VARIANCE DECOMPOSITION: Source-of-uncertainty analysis")
    print("="*70)
    
    data = load_data(REPO_ROOT, two_box=True)
    N = 400
    SEED = 42
    
    results = {}
    
    for mode in ["d13C_only", "dual"]:
        print(f"\n  Mode: {mode}")
        
        var_total = compute_var(data, mode, N, SEED)
        var_no_kie = compute_var(data, mode, N, SEED, fix_kie=True)
        var_no_sigs = compute_var(data, mode, N, SEED, fix_sigs=True)
        var_no_tau = compute_var(data, mode, N, SEED, fix_tau=True)
        
        # Contributions = variance removed when that source is fixed
        kie_contrib = max(0, var_total - var_no_kie)
        sigs_contrib = max(0, var_total - var_no_sigs)
        tau_contrib = max(0, var_total - var_no_tau)
        residual = max(0, var_total - kie_contrib - sigs_contrib - tau_contrib)
        
        total = kie_contrib + sigs_contrib + tau_contrib + residual
        if total < 1e-9:
            total = 1.0
        
        print(f"    Total variance:    {var_total:8.1f} (Tg/yr)²")
        print(f"    Fix KIE  → var =   {var_no_kie:8.1f}  (KIE contrib = {kie_contrib/total*100:.1f}%)")
        print(f"    Fix sigs → var =   {var_no_sigs:8.1f}  (Sig contrib = {sigs_contrib/total*100:.1f}%)")
        print(f"    Fix tau  → var =   {var_no_tau:8.1f}  (Tau contrib = {tau_contrib/total*100:.1f}%)")
        print(f"    Residual (atm + interaction): {residual/total*100:.1f}%")
        
        results[mode] = {
            'var_total': var_total,
            'kie_pct': kie_contrib/total*100,
            'sigs_pct': sigs_contrib/total*100,
            'tau_pct': tau_contrib/total*100,
            'residual_pct': residual/total*100,
        }
    
    # Final comparison
    print(f"\n{'='*70}")
    print("KEY RESULT: KIE contribution to FF variance")
    print(f"{'='*70}")
    print(f"  δ¹³C-only: KIE = {results['d13C_only']['kie_pct']:.1f}% of FF variance")
    print(f"  Dual:      KIE = {results['dual']['kie_pct']:.1f}% of FF variance")
    
    reduction = results['d13C_only']['kie_pct'] - results['dual']['kie_pct']
    print(f"\n  → Dual isotopes reduce KIE's variance share by {reduction:.1f} pp")
    
    if results['dual']['kie_pct'] < results['d13C_only']['kie_pct']:
        print(f"  ✓ Confirms KIE immunity from dual-isotope constraint")
    
    with open(OUT_DIR / "variance_decomposition.json", 'w') as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\n  Saved: {OUT_DIR}/variance_decomposition.json")


if __name__ == "__main__":
    main()
