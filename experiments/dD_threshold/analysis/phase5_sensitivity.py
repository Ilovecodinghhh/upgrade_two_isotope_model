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


from core import run_twobox as _core_run_twobox, ci_width


def run_threshold_sweep(data, multipliers, kie_mode="sampled",
                         lifetime_mode="varying", tau_fixed=9.0,
                         n_iter=300, seed=42):
    """Run a threshold sweep for a given KIE/lifetime configuration."""
    results = {}
    for mult in multipliers:
        res = _core_run_twobox(data, mult, n_iter, seed, mode="dual",
                               kie_mode=kie_mode, lifetime_mode=lifetime_mode,
                               tau_fixed=tau_fixed)
        results[mult] = ci_width(res['FF_G'])
    return results

def main():
    print("="*70)
    print("PHASE 5: SENSITIVITY ANALYSIS")
    print("="*70)
    
    data = load_data(REPO_ROOT, two_box=True)
    multipliers = [1.0, 2.0, 3.0, 5.0, 8.0]
    N = 300  # Faster for sensitivity
    
    # Compute d13C-only reference dynamically
    ref = _core_run_twobox(data, 1.0, N, 42, mode="d13C_only")
    d13C_ref = ci_width(ref['FF_G'])
    
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
