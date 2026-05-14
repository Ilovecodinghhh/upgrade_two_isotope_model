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
    sample_source_signatures, sample_source_signatures_hemi,
    sample_atm_d13C, sample_atm_dD, sample_atm_dD_hemi,
)

OUT_DIR = Path(__file__).resolve().parent.parent / "results" / "phase3_threshold"
OUT_DIR.mkdir(parents=True, exist_ok=True)


from core import run_twobox as _core_run_twobox, inflate_dD_uncertainty, ci_width


def run_twobox_with_inflation(data, multiplier, n_iter=500, seed=42, mode="dual"):
    """Run 2-box model with inflated δD uncertainty."""
    return _core_run_twobox(data, multiplier, n_iter, seed, mode=mode)

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
