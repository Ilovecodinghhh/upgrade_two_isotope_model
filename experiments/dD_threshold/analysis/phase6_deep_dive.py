#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
phase6_deep_dive.py — Extended Analysis (plan.md items)
========================================================

Covers:
A. Fine-grid threshold interpolation (exact crossover point)
B. Hemispheric CI breakdown (which hemisphere drives the threshold?)
C. Bootstrap CI on CI (uncertainty on uncertainty estimates)
D. Year-range sensitivity (does pre-2005 padding bias results?)
E. Bound-hit diagnostics (how often does LSQ clamp?)
"""

import sys
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear
from scipy.interpolate import interp1d

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
    sample_source_signatures, sample_source_signatures_hemi,
    sample_atm_d13C, sample_atm_dD, sample_atm_dD_hemi,
)

OUT_DIR = Path(__file__).resolve().parent.parent / "results" / "phase6_deep_dive"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# SHARED: 2-box runner with diagnostics
# ============================================================================

def inflate_dD_uncertainty(sigs, rng, multiplier):
    if multiplier <= 1.0:
        return sigs
    extra_mic = rng.normal() * 8.25 * (multiplier - 1)
    extra_ff  = rng.normal() * 0.70 * (multiplier - 1)
    extra_bb  = rng.normal() * 7.09 * (multiplier - 1)
    sigs_new = dict(sigs)
    sigs_new['mic_dD'] = sigs['mic_dD'] + extra_mic
    sigs_new['ff_dD']  = sigs['ff_dD']  + extra_ff
    sigs_new['bb_dD']  = sigs['bb_dD']  + extra_bb
    for hemi in ('NH', 'SH'):
        for src, extra in [('mic_dD', extra_mic), ('ff_dD', extra_ff), ('bb_dD', extra_bb)]:
            key = f'{src}_{hemi}'
            if key in sigs_new:
                sigs_new[key] = sigs_new[key] + extra
    return sigs_new


def run_twobox(data, multiplier, n_iter, seed, mode="dual",
               year_start_idx=0, track_bounds=False):
    """Run 2-box with optional year clipping and bound-hit tracking."""
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
    bound_hits = np.zeros((n, n_iter), dtype=int) if track_bounds else None

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
            ex_NH = (M_SH - M_NH)/tau_ex; ex_SH = (M_NH - M_SH)/tau_ex
            S_NH[i] = (M_NH1 - M_NH) + M_NH/tau_NH[i] - ex_NH
            S_SH[i] = (M_SH1 - M_SH) + M_SH/tau_SH[i] - ex_SH

        d13C_glob_mc = sample_atm_d13C(data, k, n)
        nc = min(len(c13_glob), n+1)
        d13C_off = d13C_glob_mc[:nc] - c13_glob[:nc]
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
            ex13_NH = (n13_SH - n13_NH)/tau_ex; ex13_SH = (n13_NH - n13_SH)/tau_ex
            d13C_src_NH[j] = (n13_NH1 - n13_NH + n13_NH*a13_NH/tau_NH[j] - ex13_NH)/S_NH[j]
            d13C_src_SH[j] = (n13_SH1 - n13_SH + n13_SH*a13_SH/tau_SH[j] - ex13_SH)/S_SH[j]

        sigs = sample_source_signatures_hemi(rng, data, k, n)
        if multiplier > 1.0 and mode == "dual":
            sigs = inflate_dD_uncertainty(sigs, rng, multiplier)

        f13_bb  = delta_to_fraction_d13C(sigs['bb_d13C'])
        f13_ff  = delta_to_fraction_d13C(sigs['ff_d13C'])
        f13_mic = delta_to_fraction_d13C(sigs['mic_d13C'])
        # Hemispheric d13C source signatures
        f13_bb_NH  = delta_to_fraction_d13C(sigs['bb_d13C_NH'])
        f13_ff_NH  = delta_to_fraction_d13C(sigs['ff_d13C_NH'])
        f13_mic_NH = delta_to_fraction_d13C(sigs['mic_d13C_NH'])
        f13_bb_SH  = delta_to_fraction_d13C(sigs['bb_d13C_SH'])
        f13_ff_SH  = delta_to_fraction_d13C(sigs['ff_d13C_SH'])
        f13_mic_SH = delta_to_fraction_d13C(sigs['mic_d13C_SH'])

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
                exD_NH = (nD_SH - nD_NH)/tau_ex; exD_SH = (nD_NH - nD_SH)/tau_ex
                dD_src_NH[j] = (nD_NH1 - nD_NH + nD_NH*aD_NH/tau_NH[j] - exD_NH)/S_NH[j]
                dD_src_SH[j] = (nD_SH1 - nD_SH + nD_SH*aD_SH/tau_SH[j] - exD_SH)/S_SH[j]

            fD_bb_NH  = delta_to_fraction_dD(sigs['bb_dD_NH'])
            fD_ff_NH  = delta_to_fraction_dD(sigs['ff_dD_NH'])
            fD_mic_NH = delta_to_fraction_dD(sigs['mic_dD_NH'])
            fD_bb_SH  = delta_to_fraction_dD(sigs['bb_dD_SH'])
            fD_ff_SH  = delta_to_fraction_dD(sigs['ff_dD_SH'])
            fD_mic_SH = delta_to_fraction_dD(sigs['mic_dD_SH'])

            for j in range(n):
                A_nh = np.array([[1,1,1],
                    [f13_bb_NH[j], f13_ff_NH[j], f13_mic_NH[j]],
                    [fD_bb_NH[j], fD_ff_NH[j], fD_mic_NH[j]]])
                B_nh = np.array([S_NH[j], S_NH[j]*d13C_src_NH[j], S_NH[j]*dD_src_NH[j]])
                try:
                    res = lsq_linear(W_NH@A_nh, W_NH@B_nh, bounds=(0, S_NH[j]*1.5))
                    FF_NH[j,k] = res.x[1]
                    if track_bounds:
                        lb = np.any(res.x < 1e-6)
                        ub = np.any(res.x > S_NH[j]*1.5 - 1e-6)
                        if lb or ub: bound_hits[j,k] |= 1
                except: pass

                A_sh = np.array([[1,1,1],
                    [f13_bb_SH[j], f13_ff_SH[j], f13_mic_SH[j]],
                    [fD_bb_SH[j], fD_ff_SH[j], fD_mic_SH[j]]])
                B_sh = np.array([S_SH[j], S_SH[j]*d13C_src_SH[j], S_SH[j]*dD_src_SH[j]])
                try:
                    res = lsq_linear(W_SH@A_sh, W_SH@B_sh, bounds=(0, S_SH[j]*1.5))
                    FF_SH[j,k] = res.x[1]
                    if track_bounds:
                        lb = np.any(res.x < 1e-6)
                        ub = np.any(res.x > S_SH[j]*1.5 - 1e-6)
                        if lb or ub: bound_hits[j,k] |= 2
                except: pass

        else:  # d13C_only
            for j in range(n):
                # NH
                denom_nh = f13_ff_NH[j] - f13_mic_NH[j]
                if abs(denom_nh) < 1e-15:
                    FF_NH[j,k] = np.nan
                else:
                    S_rem = S_NH[j] - BB_hemi_NH
                    rhs = S_NH[j]*d13C_src_NH[j] - BB_hemi_NH*f13_bb_NH[j]
                    FF_NH[j,k] = (rhs - S_rem*f13_mic_NH[j]) / denom_nh
                # SH
                denom_sh = f13_ff_SH[j] - f13_mic_SH[j]
                if abs(denom_sh) < 1e-15:
                    FF_SH[j,k] = np.nan
                else:
                    S_rem = S_SH[j] - BB_hemi_SH
                    rhs = S_SH[j]*d13C_src_SH[j] - BB_hemi_SH*f13_bb_SH[j]
                    FF_SH[j,k] = (rhs - S_rem*f13_mic_SH[j]) / denom_sh

    return {
        'FF_NH': FF_NH, 'FF_SH': FF_SH, 'FF_G': FF_NH + FF_SH,
        'bound_hits': bound_hits,
    }


def ci_width(arr, start_idx=8):
    """90% CI width over years starting at start_idx."""
    s = smooth_5yr(arr)
    ci5  = np.nanpercentile(s[start_idx:], 5, axis=1)
    ci95 = np.nanpercentile(s[start_idx:], 95, axis=1)
    return float(np.nanmean(ci95 - ci5))


def ci_width_hemi(arr, start_idx=8):
    """90% CI width for a single hemisphere's FF array."""
    s = smooth_5yr(arr)
    ci5  = np.nanpercentile(s[start_idx:], 5, axis=1)
    ci95 = np.nanpercentile(s[start_idx:], 95, axis=1)
    return float(np.nanmean(ci95 - ci5))


# ============================================================================
# A: FINE-GRID THRESHOLD
# ============================================================================

def fine_grid_threshold(data):
    print("\n" + "="*70)
    print("A. FINE-GRID THRESHOLD INTERPOLATION")
    print("="*70)

    # Coarse grid first, then fine around crossover
    coarse = [0.5, 1.0, 2.0, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0]
    N = 500
    
    ref = run_twobox(data, 1.0, N, 42, mode="d13C_only")
    ref_ci = ci_width(ref['FF_G'])
    print(f"  d13C-only reference CI: {ref_ci:.1f} Tg/yr")

    results = {}
    for mult in coarse:
        res = run_twobox(data, mult, N, 42, mode="dual")
        w = ci_width(res['FF_G'])
        imp = (ref_ci - w) / ref_ci * 100
        results[mult] = {'ci': w, 'improvement': imp}
        print(f"  {mult:5.1f}x: CI={w:6.1f}  improvement={imp:+.1f}%")

    # Interpolate to find exact crossover (improvement = 0)
    mults = sorted(results.keys())
    imps = [results[m]['improvement'] for m in mults]
    
    # Find where improvement crosses 0 and 10
    f_interp = interp1d(mults, imps, kind='linear', fill_value='extrapolate')
    
    # Search for 0% crossover
    from scipy.optimize import brentq
    try:
        cross_0 = brentq(lambda x: float(f_interp(x)), 3.0, 8.0)
        sigma_0 = 8.25 * cross_0
        print(f"\n  Exact crossover (improvement = 0%): {cross_0:.2f}x (sigma = {sigma_0:.1f} permil)")
    except:
        cross_0 = None; sigma_0 = None
        print("  Could not find 0% crossover in range")

    # Search for 10% crossover
    try:
        cross_10 = brentq(lambda x: float(f_interp(x)) - 10.0, 2.0, 8.0)
        sigma_10 = 8.25 * cross_10
        print(f"  10% improvement threshold:           {cross_10:.2f}x (sigma = {sigma_10:.1f} permil)")
    except:
        cross_10 = None; sigma_10 = None
        print("  Could not find 10% crossover in range")

    summary = {
        'ref_ci': ref_ci,
        'sweep': {str(m): results[m] for m in mults},
        'crossover_0pct': {'multiplier': cross_0, 'sigma_permil': sigma_0},
        'crossover_10pct': {'multiplier': cross_10, 'sigma_permil': sigma_10},
    }
    return summary


# ============================================================================
# B: HEMISPHERIC CI BREAKDOWN
# ============================================================================

def hemispheric_breakdown(data):
    print("\n" + "="*70)
    print("B. HEMISPHERIC CI BREAKDOWN")
    print("="*70)
    
    mults = [1.0, 3.0, 5.0, 8.0]
    N = 500
    
    ref = run_twobox(data, 1.0, N, 42, mode="d13C_only")
    ref_nh = ci_width_hemi(ref['FF_NH'])
    ref_sh = ci_width_hemi(ref['FF_SH'])
    ref_g  = ci_width(ref['FF_G'])
    print(f"  d13C-only:  NH={ref_nh:.1f}  SH={ref_sh:.1f}  Global={ref_g:.1f}")
    
    results = {'d13C_only': {'NH': ref_nh, 'SH': ref_sh, 'Global': ref_g}}
    
    for mult in mults:
        res = run_twobox(data, mult, N, 42, mode="dual")
        nh_ci = ci_width_hemi(res['FF_NH'])
        sh_ci = ci_width_hemi(res['FF_SH'])
        g_ci  = ci_width(res['FF_G'])
        imp_nh = (ref_nh - nh_ci) / ref_nh * 100
        imp_sh = (ref_sh - sh_ci) / ref_sh * 100
        imp_g  = (ref_g  - g_ci)  / ref_g  * 100
        results[f'dual_{mult}x'] = {
            'NH': nh_ci, 'SH': sh_ci, 'Global': g_ci,
            'imp_NH': imp_nh, 'imp_SH': imp_sh, 'imp_G': imp_g,
        }
        print(f"  {mult:.1f}x dual:  NH={nh_ci:.1f} ({imp_nh:+.1f}%)  "
              f"SH={sh_ci:.1f} ({imp_sh:+.1f}%)  Global={g_ci:.1f} ({imp_g:+.1f}%)")
    
    return results


# ============================================================================
# C: BOOTSTRAP CI ON CI
# ============================================================================

def bootstrap_ci_on_ci(data):
    print("\n" + "="*70)
    print("C. BOOTSTRAP CI ON CI")
    print("="*70)
    
    N = 1000  # More iters for bootstrapping
    n_bootstrap = 200
    
    # Run once with more iterations
    ref = run_twobox(data, 1.0, N, 42, mode="d13C_only")
    dual = run_twobox(data, 1.0, N, 42, mode="dual")
    
    FF_ref = smooth_5yr(ref['FF_G'])[8:]   # post-2007
    FF_dual = smooth_5yr(dual['FF_G'])[8:]
    
    rng = np.random.default_rng(123)
    
    boot_ref_ci = []
    boot_dual_ci = []
    boot_improvement = []
    
    for _ in range(n_bootstrap):
        idx = rng.choice(N, N, replace=True)
        ref_boot = FF_ref[:, idx]
        dual_boot = FF_dual[:, idx]
        
        ref_w = float(np.nanmean(
            np.nanpercentile(ref_boot, 95, axis=1) - 
            np.nanpercentile(ref_boot, 5, axis=1)))
        dual_w = float(np.nanmean(
            np.nanpercentile(dual_boot, 95, axis=1) - 
            np.nanpercentile(dual_boot, 5, axis=1)))
        
        boot_ref_ci.append(ref_w)
        boot_dual_ci.append(dual_w)
        boot_improvement.append((ref_w - dual_w) / ref_w * 100)
    
    boot_ref_ci = np.array(boot_ref_ci)
    boot_dual_ci = np.array(boot_dual_ci)
    boot_improvement = np.array(boot_improvement)
    
    print(f"  d13C-only CI: {np.mean(boot_ref_ci):.1f} +/- {np.std(boot_ref_ci):.1f} Tg/yr "
          f"[{np.percentile(boot_ref_ci, 2.5):.1f}, {np.percentile(boot_ref_ci, 97.5):.1f}]")
    print(f"  Dual CI:      {np.mean(boot_dual_ci):.1f} +/- {np.std(boot_dual_ci):.1f} Tg/yr "
          f"[{np.percentile(boot_dual_ci, 2.5):.1f}, {np.percentile(boot_dual_ci, 97.5):.1f}]")
    print(f"  Improvement:  {np.mean(boot_improvement):.1f} +/- {np.std(boot_improvement):.1f}% "
          f"[{np.percentile(boot_improvement, 2.5):.1f}, {np.percentile(boot_improvement, 97.5):.1f}]")
    print(f"  P(improvement > 0): {np.mean(boot_improvement > 0)*100:.1f}%")
    print(f"  P(improvement > 30): {np.mean(boot_improvement > 30)*100:.1f}%")
    
    return {
        'd13C_CI': {
            'mean': float(np.mean(boot_ref_ci)),
            'std': float(np.std(boot_ref_ci)),
            'CI95': [float(np.percentile(boot_ref_ci, 2.5)),
                     float(np.percentile(boot_ref_ci, 97.5))],
        },
        'dual_CI': {
            'mean': float(np.mean(boot_dual_ci)),
            'std': float(np.std(boot_dual_ci)),
            'CI95': [float(np.percentile(boot_dual_ci, 2.5)),
                     float(np.percentile(boot_dual_ci, 97.5))],
        },
        'improvement_pct': {
            'mean': float(np.mean(boot_improvement)),
            'std': float(np.std(boot_improvement)),
            'CI95': [float(np.percentile(boot_improvement, 2.5)),
                     float(np.percentile(boot_improvement, 97.5))],
            'prob_positive': float(np.mean(boot_improvement > 0)),
            'prob_gt30': float(np.mean(boot_improvement > 30)),
        },
        'n_mc': N, 'n_bootstrap': n_bootstrap,
    }


# ============================================================================
# D: YEAR-RANGE SENSITIVITY
# ============================================================================

def year_range_sensitivity(data):
    print("\n" + "="*70)
    print("D. YEAR-RANGE SENSITIVITY")
    print("="*70)
    
    N = 500
    n = data.n_years
    years = data.model_years
    
    # Full range (1999–2021, idx 0–22) — using all years for CI
    # δD-valid range (2005–2021, idx 6–22) — only years with real δD data
    
    configs = [
        ("Full (1999-2021)", 0),   # start CI from year index 0
        ("Post-padding (2005-2021)", 6),  # start CI from year index 6
        ("Post-2007 (2007-2021)", 8),  # standard analysis window
    ]
    
    results = {}
    for label, start_idx in configs:
        ref = run_twobox(data, 1.0, N, 42, mode="d13C_only")
        dual = run_twobox(data, 1.0, N, 42, mode="dual")
        
        ref_ci = ci_width(ref['FF_G'], start_idx=start_idx)
        dual_ci = ci_width(dual['FF_G'], start_idx=start_idx)
        imp = (ref_ci - dual_ci) / ref_ci * 100
        
        results[label] = {
            'd13C_CI': ref_ci, 'dual_CI': dual_ci, 'improvement': imp,
        }
        print(f"  {label:30s}: d13C={ref_ci:.1f}  dual={dual_ci:.1f}  improvement={imp:+.1f}%")
    
    return results


# ============================================================================
# E: BOUND-HIT DIAGNOSTICS
# ============================================================================

def bound_diagnostics(data):
    print("\n" + "="*70)
    print("E. BOUND-HIT DIAGNOSTICS")
    print("="*70)
    
    N = 500
    mults = [1.0, 3.0, 5.0, 8.0]
    
    results = {}
    for mult in mults:
        res = run_twobox(data, mult, N, 42, mode="dual", track_bounds=True)
        bh = res['bound_hits']
        nh_hits = np.sum(bh & 1) / bh.size * 100
        sh_hits = np.sum(bh & 2) / bh.size * 100
        any_hits = np.sum(bh > 0) / bh.size * 100
        
        # Per-year breakdown (post-2007)
        bh_post = bh[8:]
        any_post = np.sum(bh_post > 0) / bh_post.size * 100
        
        results[f'{mult}x'] = {
            'NH_bound_pct': nh_hits,
            'SH_bound_pct': sh_hits,
            'any_bound_pct': any_hits,
            'post2007_bound_pct': any_post,
        }
        print(f"  {mult:.1f}x: NH={nh_hits:.1f}% SH={sh_hits:.1f}% "
              f"any={any_hits:.1f}% (post-2007: {any_post:.1f}%)")
    
    return results


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("PHASE 6: DEEP DIVE (plan.md items)")
    print("=" * 70)
    
    data = load_data(REPO_ROOT, two_box=True)
    
    all_results = {}
    
    # A: Fine grid
    all_results['A_fine_grid'] = fine_grid_threshold(data)
    
    # B: Hemispheric breakdown
    all_results['B_hemispheric'] = hemispheric_breakdown(data)
    
    # C: Bootstrap
    all_results['C_bootstrap'] = bootstrap_ci_on_ci(data)
    
    # D: Year range
    all_results['D_year_range'] = year_range_sensitivity(data)
    
    # E: Bound hits
    all_results['E_bound_hits'] = bound_diagnostics(data)
    
    # Save
    with open(OUT_DIR / "deep_dive_results.json", 'w') as f:
        json.dump(all_results, f, indent=2, default=float)
    
    print(f"\n{'='*70}")
    print("SAVED:", OUT_DIR / "deep_dive_results.json")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
