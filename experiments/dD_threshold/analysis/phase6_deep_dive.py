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

from core import run_twobox, inflate_dD_uncertainty, ci_width, ci_width_hemi


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
