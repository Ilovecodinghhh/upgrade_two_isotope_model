#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
phase2_dfs.py — Degrees of Freedom for Signal (DFS) Calculator
================================================================

Computes the theoretical information content of δ¹³C-only vs dual-isotope
observations for methane source attribution, at 1-box and 2-box resolution.

DFS = trace(HBHᵀ(HBHᵀ + R)⁻¹)

Where:
  H = Jacobian (∂observations / ∂sources)
  B = Prior covariance of sources (from MC spread of signatures)
  R = Observation error covariance
"""

import sys
import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))
from common import (
    load_data, sample_KIE, compute_bulk_KIE,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    SINK_FRACTIONS_GLOBAL, SINK_FRACTIONS_NH, SINK_FRACTIONS_SH,
    KIE_FIXED, PT, PT_HEMI,
)

OUT_DIR = Path(__file__).resolve().parent.parent / "results" / "phase2_dfs"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def compute_jacobian_1box(f13_bb, f13_ff, f13_mic, fD_bb=None, fD_ff=None, fD_mic=None, mode="dual"):
    """
    Compute Jacobian H for the 1-box model.
    
    For a linear system Ax = b where:
      Row 0: mass balance (∂total_obs / ∂BB = 1, ∂/∂FF = 1, ∂/∂Mic = 1)
      Row 1: δ¹³C (∂δ¹³C_obs / ∂BB = f13_bb, etc.)
      Row 2: δD (∂δD_obs / ∂BB = fD_bb, etc.)  [only in dual mode]
    
    The Jacobian H maps from source space (3 sources: BB, FF, Mic) to 
    observation space (2 or 3 observations).
    """
    if mode == "dual":
        H = np.array([
            [1.0, 1.0, 1.0],
            [f13_bb, f13_ff, f13_mic],
            [fD_bb, fD_ff, fD_mic],
        ])
    else:
        H = np.array([
            [1.0, 1.0, 1.0],
            [f13_bb, f13_ff, f13_mic],
        ])
    return H


def compute_prior_covariance(data, n_samples=500, seed=42):
    """
    Estimate prior covariance B of source emissions from MC signature spread.
    
    Uses the spread of MC-sampled source signatures to estimate how much
    prior uncertainty there is in the source partition.
    
    Returns B (3×3) and source signature stats.
    """
    rng = np.random.default_rng(seed)
    
    # Sample source-signature variability
    from common import sample_source_signatures, compute_lifetime, pad_to_length
    
    # Use year 2015 (index ~16) as representative
    yr_idx = 16
    
    ff_d13C_samples = []
    mic_d13C_samples = []
    bb_d13C_samples = []
    ff_dD_samples = []
    mic_dD_samples = []
    bb_dD_samples = []
    
    for k in range(min(n_samples, 1000)):
        sigs = sample_source_signatures(rng, data, k, data.n_years)
        ff_d13C_samples.append(sigs['ff_d13C'][yr_idx])
        mic_d13C_samples.append(sigs['mic_d13C'][yr_idx])
        bb_d13C_samples.append(sigs['bb_d13C'][yr_idx])
        ff_dD_samples.append(sigs['ff_dD'][yr_idx])
        mic_dD_samples.append(sigs['mic_dD'][yr_idx])
        bb_dD_samples.append(sigs['bb_dD'][yr_idx])
    
    stats = {
        'ff_d13C': {'mean': np.mean(ff_d13C_samples), 'std': np.std(ff_d13C_samples)},
        'mic_d13C': {'mean': np.mean(mic_d13C_samples), 'std': np.std(mic_d13C_samples)},
        'bb_d13C': {'mean': np.mean(bb_d13C_samples), 'std': np.std(bb_d13C_samples)},
        'ff_dD': {'mean': np.mean(ff_dD_samples), 'std': np.std(ff_dD_samples)},
        'mic_dD': {'mean': np.mean(mic_dD_samples), 'std': np.std(mic_dD_samples)},
        'bb_dD': {'mean': np.mean(bb_dD_samples), 'std': np.std(bb_dD_samples)},
    }
    
    # Prior emission uncertainties (Tg/yr, 1σ)
    # Based on literature spread: FF ~100±30, Mic ~230±50, BB ~30±15
    B = np.diag([15.0**2, 30.0**2, 50.0**2])  # BB, FF, Mic
    
    return B, stats


def compute_obs_covariance(mode="dual"):
    """
    Observation error covariance R.
    
    For mass balance: total source uncertainty ~ 20 Tg/yr (from lifetime + CH₄ growth)
    For δ¹³C: measurement uncertainty 0.04‰ → convert to fraction space
    For δD: measurement uncertainty 3‰ → convert to fraction space
    """
    # Convert δ uncertainties to fraction-space uncertainties
    # d(fraction)/d(delta) ≈ std / 1000 * C13_STD for d13C
    # At δ¹³C ≈ -47‰: df/ddelta ≈ C13_STD / 1000 ≈ 1.1e-5 per ‰
    # At δD ≈ -86‰: df/ddelta ≈ D_STD / 1000 ≈ 1.56e-7 per ‰
    
    sigma_mass = 20.0  # Tg/yr
    sigma_d13C_frac = 0.04 * 0.011113 / 1000.0  # 0.04‰ in fraction space
    sigma_dD_frac = 3.0 * 0.00015576 / 1000.0   # 3‰ in fraction space
    
    # Scale to emission-weighted observation: multiply by total source (~560 Tg/yr)
    S_total = 560.0
    sigma_d13C_em = sigma_d13C_frac * S_total
    sigma_dD_em = sigma_dD_frac * S_total
    
    if mode == "dual":
        R = np.diag([sigma_mass**2, sigma_d13C_em**2, sigma_dD_em**2])
    else:
        R = np.diag([sigma_mass**2, sigma_d13C_em**2])
    return R


def compute_DFS(H, B, R):
    """
    DFS = trace(HBHᵀ (HBHᵀ + R)⁻¹)
    
    This measures how many independent pieces of information the observations
    provide about the sources.
    """
    HBHT = H @ B @ H.T
    M = HBHT + R
    try:
        DFS = np.trace(HBHT @ np.linalg.inv(M))
    except np.linalg.LinAlgError:
        DFS = np.nan
    return DFS


def main():
    print("=" * 60)
    print("PHASE 2: DFS CALCULATOR")
    print("=" * 60)
    
    data = load_data(REPO_ROOT, two_box=False)
    B, sig_stats = compute_prior_covariance(data)
    
    print(f"\nSource signature statistics (year 2015):")
    for k, v in sig_stats.items():
        print(f"  {k}: {v['mean']:.2f} ± {v['std']:.2f} ‰")
    
    # Compute representative Jacobians using mean signatures
    f13_bb = delta_to_fraction_d13C(sig_stats['bb_d13C']['mean'])
    f13_ff = delta_to_fraction_d13C(sig_stats['ff_d13C']['mean'])
    f13_mic = delta_to_fraction_d13C(sig_stats['mic_d13C']['mean'])
    fD_bb = delta_to_fraction_dD(sig_stats['bb_dD']['mean'])
    fD_ff = delta_to_fraction_dD(sig_stats['ff_dD']['mean'])
    fD_mic = delta_to_fraction_dD(sig_stats['mic_dD']['mean'])
    
    results = {}
    
    # === 1-BOX ===
    # d13C-only
    H_1box_c = compute_jacobian_1box(f13_bb, f13_ff, f13_mic, mode="d13C_only")
    R_1box_c = compute_obs_covariance("d13C_only")
    dfs_1box_c = compute_DFS(H_1box_c, B, R_1box_c)
    
    # dual
    H_1box_d = compute_jacobian_1box(f13_bb, f13_ff, f13_mic, fD_bb, fD_ff, fD_mic, mode="dual")
    R_1box_d = compute_obs_covariance("dual")
    dfs_1box_d = compute_DFS(H_1box_d, B, R_1box_d)
    
    results['onebox_d13C_only'] = {'DFS': float(dfs_1box_c)}
    results['onebox_dual'] = {'DFS': float(dfs_1box_d)}
    results['onebox_delta_DFS'] = float(dfs_1box_d - dfs_1box_c)
    
    # === 2-BOX ===
    # For 2-box: 6 sources (BB_NH, FF_NH, Mic_NH, BB_SH, FF_SH, Mic_SH)
    # Observations: 2 or 3 per hemisphere = 4 or 6 total
    B_2box = np.diag([
        7.5**2, 15.0**2, 25.0**2,   # NH: BB, FF, Mic (half of global uncertainties)
        7.5**2, 15.0**2, 25.0**2,   # SH
    ])
    
    # d13C-only: 4 observations (mass_NH, d13C_NH, mass_SH, d13C_SH)
    H_2box_c = np.array([
        [1, 1, 1, 0, 0, 0],              # mass_NH
        [f13_bb, f13_ff, f13_mic, 0, 0, 0],  # d13C_NH
        [0, 0, 0, 1, 1, 1],              # mass_SH
        [0, 0, 0, f13_bb, f13_ff, f13_mic],  # d13C_SH
    ])
    S_hemi = 280.0
    sigma_mass_h = 10.0
    sigma_d13C_em_h = 0.04 * 0.011113 / 1000.0 * S_hemi
    R_2box_c = np.diag([sigma_mass_h**2, sigma_d13C_em_h**2,
                         sigma_mass_h**2, sigma_d13C_em_h**2])
    dfs_2box_c = compute_DFS(H_2box_c, B_2box, R_2box_c)
    
    # dual: 6 observations (mass_NH, d13C_NH, dD_NH, mass_SH, d13C_SH, dD_SH)
    sigma_dD_em_h = 3.0 * 0.00015576 / 1000.0 * S_hemi
    H_2box_d = np.array([
        [1, 1, 1, 0, 0, 0],
        [f13_bb, f13_ff, f13_mic, 0, 0, 0],
        [fD_bb, fD_ff, fD_mic, 0, 0, 0],
        [0, 0, 0, 1, 1, 1],
        [0, 0, 0, f13_bb, f13_ff, f13_mic],
        [0, 0, 0, fD_bb, fD_ff, fD_mic],
    ])
    R_2box_d = np.diag([sigma_mass_h**2, sigma_d13C_em_h**2, sigma_dD_em_h**2,
                         sigma_mass_h**2, sigma_d13C_em_h**2, sigma_dD_em_h**2])
    dfs_2box_d = compute_DFS(H_2box_d, B_2box, R_2box_d)
    
    results['twobox_d13C_only'] = {'DFS': float(dfs_2box_c)}
    results['twobox_dual'] = {'DFS': float(dfs_2box_d)}
    results['twobox_delta_DFS'] = float(dfs_2box_d - dfs_2box_c)
    
    # === Condition numbers ===
    results['onebox_cond_d13C'] = float(np.linalg.cond(H_1box_c))
    results['onebox_cond_dual'] = float(np.linalg.cond(H_1box_d))
    results['twobox_cond_d13C'] = float(np.linalg.cond(H_2box_c))
    results['twobox_cond_dual'] = float(np.linalg.cond(H_2box_d))
    
    # === Source signature stats for reference ===
    results['source_signatures'] = {k: {kk: float(vv) for kk, vv in v.items()} 
                                     for k, v in sig_stats.items()}
    
    # Print results
    print(f"\n{'='*60}")
    print("DFS RESULTS")
    print(f"{'='*60}")
    print(f"\n  1-Box:")
    print(f"    δ¹³C-only DFS: {dfs_1box_c:.3f}")
    print(f"    Dual DFS:      {dfs_1box_d:.3f}")
    print(f"    ΔDFS:          {dfs_1box_d - dfs_1box_c:.3f}")
    print(f"    Cond(d13C):    {np.linalg.cond(H_1box_c):.1f}")
    print(f"    Cond(dual):    {np.linalg.cond(H_1box_d):.1f}")
    print(f"\n  2-Box:")
    print(f"    δ¹³C-only DFS: {dfs_2box_c:.3f}")
    print(f"    Dual DFS:      {dfs_2box_d:.3f}")
    print(f"    ΔDFS:          {dfs_2box_d - dfs_2box_c:.3f}")
    print(f"    Cond(d13C):    {np.linalg.cond(H_2box_c):.1f}")
    print(f"    Cond(dual):    {np.linalg.cond(H_2box_d):.1f}")
    
    # Save
    with open(OUT_DIR / "dfs_results.json", 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n  Saved to: {OUT_DIR}/dfs_results.json")


if __name__ == "__main__":
    main()
