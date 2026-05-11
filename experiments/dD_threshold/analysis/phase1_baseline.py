#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
phase1_baseline.py — Run 1-box and 2-box models in δ¹³C-only vs dual-isotope mode
==================================================================================

Phase 1 of the δD Threshold experiment:
  - Runs 3x3_one equivalent in δ¹³C-only mode (2×2 system: mass + δ¹³C, BB free)
  - Runs 3x3_one in dual-isotope mode (existing 3×3)
  - Runs 3x3_two equivalent in δ¹³C-only mode (2×2 per hemisphere)
  - Runs 3x3_two in dual-isotope mode (existing 3×3 per hemisphere)

Saves all iteration data as .npz for downstream analysis.
"""

import sys
import json
from pathlib import Path

import numpy as np
from scipy.optimize import lsq_linear

# Add repo root to path (experiments/dD_threshold/analysis/ → 3 levels up to workspace)
REPO_ROOT = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, str(REPO_ROOT))
from common import (
    ModelConfig, QualityMonitor, load_data,
    sample_KIE, compute_bulk_KIE, compute_lifetime,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    sample_source_signatures, sample_atm_d13C, sample_atm_dD,
    smooth_5yr, pad_to_length,
    SINK_FRACTIONS_GLOBAL, SINK_FRACTIONS_NH, SINK_FRACTIONS_SH,
    PT, PT_HEMI, LIFETIME_RATIO_NH, LIFETIME_RATIO_SH,
    DD_IH_OFFSET, TAU_EX_MEAN, TAU_EX_STD,
)

BASE_DIR = REPO_ROOT
OUT_DIR = Path(__file__).resolve().parent.parent / "results" / "phase1_baseline"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def run_one_box(mode="dual", n_iter=1000, seed=42):
    """Run global 1-box model in 'dual' or 'd13C_only' mode.
    
    In d13C_only mode: solve a 2×2 system (mass + δ¹³C) for FF and Mic,
    with BB fixed from CarbonTracker mean.
    
    In dual mode: solve full 3×3 (mass + δ¹³C + δD) for BB, FF, Mic.
    """
    print(f"\n{'='*60}")
    print(f"ONE-BOX MODEL — mode={mode}, N={n_iter}")
    print(f"{'='*60}")
    
    cfg = ModelConfig(n_iterations=n_iter, kie_mode="sampled",
                      lifetime_mode="varying", seed=seed)
    data = load_data(BASE_DIR, two_box=False)
    n = data.n_years
    years = data.model_years
    CH4 = data.CH4_global
    rng = np.random.default_rng(seed)
    
    tau = compute_lifetime(years, cfg.lifetime_mode, cfg.tau_fixed)
    
    # Total source
    SumSource = np.zeros(n)
    for i in range(n):
        SumSource[i] = CH4[i+1] * PT - CH4[i] * PT + CH4[i] * PT / tau[i]
    
    NI = n_iter
    BB_comp = np.zeros((n, NI))
    FF_comp = np.zeros((n, NI))
    Mic_comp = np.zeros((n, NI))
    
    for k in range(NI):
        if (k+1) % 200 == 0:
            print(f"  iter {k+1}/{NI}")
        
        kies = sample_KIE(rng, cfg.kie_mode)
        KIE_13C, KIE_D = compute_bulk_KIE(kies, SINK_FRACTIONS_GLOBAL)
        alpha_13C = 1.0 / KIE_13C
        alpha_D = 1.0 / KIE_D
        
        d13C_atm = sample_atm_d13C(data, k, n)
        f13 = delta_to_fraction_d13C(d13C_atm)
        n13 = f13 * CH4 * PT
        
        d13C_src = np.zeros(n)
        for j in range(n):
            d13C_src[j] = (n13[j+1] - n13[j] + n13[j] * alpha_13C / tau[j]) / SumSource[j]
        
        sigs = sample_source_signatures(rng, data, k, n)
        f13_bb = delta_to_fraction_d13C(sigs['bb_d13C'])
        f13_ff = delta_to_fraction_d13C(sigs['ff_d13C'])
        f13_mic = delta_to_fraction_d13C(sigs['mic_d13C'])
        
        if mode == "dual":
            # Full 3×3
            dD_atm = sample_atm_dD(data, k, n)
            fD = delta_to_fraction_dD(dD_atm)
            nD = fD * CH4 * PT
            dD_src = np.zeros(n)
            for j in range(n):
                dD_src[j] = (nD[j+1] - nD[j] + nD[j] * alpha_D / tau[j]) / SumSource[j]
            
            fD_bb = delta_to_fraction_dD(sigs['bb_dD'])
            fD_ff = delta_to_fraction_dD(sigs['ff_dD'])
            fD_mic = delta_to_fraction_dD(sigs['mic_dD'])
            
            for j in range(n):
                A = np.array([
                    [1.0,       1.0,       1.0],
                    [f13_bb[j], f13_ff[j], f13_mic[j]],
                    [fD_bb[j],  fD_ff[j],  fD_mic[j]],
                ])
                B = np.array([SumSource[j], SumSource[j]*d13C_src[j], SumSource[j]*dD_src[j]])
                try:
                    x = np.linalg.solve(A, B)
                except np.linalg.LinAlgError:
                    x = np.array([np.nan, np.nan, np.nan])
                BB_comp[j,k] = x[0]; FF_comp[j,k] = x[1]; Mic_comp[j,k] = x[2]
        
        else:  # d13C_only
            # 2×2 system: BB fixed, solve for FF and Mic
            BB_fixed = data.BB_global_mean  # CarbonTracker mean
            for j in range(n):
                S_remaining = SumSource[j] - BB_fixed
                # δ¹³C constraint: S*f13_src = BB*f13_BB + FF*f13_FF + Mic*f13_Mic
                # With FF + Mic = S_remaining:
                #   FF*f13_FF + (S_remaining - FF)*f13_Mic = S*f13_src - BB*f13_BB
                # Solve for FF:
                rhs = SumSource[j] * d13C_src[j] - BB_fixed * f13_bb[j]
                denom = f13_ff[j] - f13_mic[j]
                if abs(denom) > 1e-15:
                    FF_val = (rhs - S_remaining * f13_mic[j]) / denom
                    Mic_val = S_remaining - FF_val
                else:
                    FF_val = np.nan; Mic_val = np.nan
                
                BB_comp[j,k] = BB_fixed
                FF_comp[j,k] = FF_val
                Mic_comp[j,k] = Mic_val
    
    print("  Done!")
    return years, BB_comp, FF_comp, Mic_comp, SumSource, tau


def run_two_box(mode="dual", n_iter=1000, seed=42):
    """Run 2-box (NH/SH) model in 'dual' or 'd13C_only' mode.
    
    In d13C_only: solve 2×2 per hemisphere (mass + δ¹³C, BB fixed).
    In dual: solve 3×3 per hemisphere (mass + δ¹³C + δD, BB free, bounded LS).
    """
    print(f"\n{'='*60}")
    print(f"TWO-BOX MODEL — mode={mode}, N={n_iter}")
    print(f"{'='*60}")
    
    cfg = ModelConfig(n_iterations=n_iter, kie_mode="sampled",
                      lifetime_mode="varying", seed=seed)
    data = load_data(BASE_DIR, two_box=True)
    n = data.n_years
    years = data.model_years
    rng = np.random.default_rng(seed)
    
    CH4_NH, CH4_SH = data.CH4_NH, data.CH4_SH
    c13_NH, c13_SH = data.c13_NH, data.c13_SH
    c13_glob = data.c13_global
    
    tau_global = compute_lifetime(years, cfg.lifetime_mode, cfg.tau_fixed)
    tau_NH = tau_global * LIFETIME_RATIO_NH
    tau_SH = tau_global * LIFETIME_RATIO_SH
    
    NI = n_iter
    FF_NH_comp = np.zeros((n, NI)); Mic_NH_comp = np.zeros((n, NI)); BB_NH_comp = np.zeros((n, NI))
    FF_SH_comp = np.zeros((n, NI)); Mic_SH_comp = np.zeros((n, NI)); BB_SH_comp = np.zeros((n, NI))
    
    # BB split
    from common import BB_NH_FRACTION, BB_SH_FRACTION
    BB_hemi_NH = data.BB_global_mean * BB_NH_FRACTION
    BB_hemi_SH = data.BB_global_mean * BB_SH_FRACTION
    
    W_NH = np.diag([100.0, 1.0, 0.5])
    W_SH = np.diag([200.0, 1.0, 0.5])
    
    for k in range(NI):
        if (k+1) % 200 == 0:
            print(f"  iter {k+1}/{NI}")
        
        tau_ex = max(0.5, rng.normal(TAU_EX_MEAN, TAU_EX_STD))
        kies = sample_KIE(rng, cfg.kie_mode)
        K13_NH, KD_NH = compute_bulk_KIE(kies, SINK_FRACTIONS_NH)
        K13_SH, KD_SH = compute_bulk_KIE(kies, SINK_FRACTIONS_SH)
        a13_NH = 1.0 / K13_NH; aD_NH = 1.0 / KD_NH
        a13_SH = 1.0 / K13_SH; aD_SH = 1.0 / KD_SH
        
        # Total source per hemisphere
        S_NH = np.zeros(n); S_SH = np.zeros(n)
        for i in range(n):
            M_NH = CH4_NH[i] * PT_HEMI; M_NH1 = CH4_NH[i+1] * PT_HEMI
            M_SH = CH4_SH[i] * PT_HEMI; M_SH1 = CH4_SH[i+1] * PT_HEMI
            ex_NH = (M_SH - M_NH) / tau_ex
            ex_SH = (M_NH - M_SH) / tau_ex
            S_NH[i] = (M_NH1 - M_NH) + M_NH / tau_NH[i] - ex_NH
            S_SH[i] = (M_SH1 - M_SH) + M_SH / tau_SH[i] - ex_SH
        
        # Atmospheric δ¹³C
        d13C_glob_MC = sample_atm_d13C(data, k, n)
        nc = min(len(c13_glob), n+1)
        d13C_off = d13C_glob_MC[:nc] - c13_glob[:nc]
        d13C_NH_MC = c13_NH[:nc] + d13C_off
        d13C_SH_MC = c13_SH[:nc] + d13C_off
        
        f13_NH_atm = delta_to_fraction_d13C(d13C_NH_MC)
        f13_SH_atm = delta_to_fraction_d13C(d13C_SH_MC)
        
        # δ¹³C source fractions
        d13C_src_NH = np.zeros(n); d13C_src_SH = np.zeros(n)
        for j in range(n):
            n13_NH = f13_NH_atm[j] * CH4_NH[j] * PT_HEMI
            n13_NH1 = f13_NH_atm[j+1] * CH4_NH[j+1] * PT_HEMI
            n13_SH = f13_SH_atm[j] * CH4_SH[j] * PT_HEMI
            n13_SH1 = f13_SH_atm[j+1] * CH4_SH[j+1] * PT_HEMI
            ex13_NH = (n13_SH - n13_NH) / tau_ex
            ex13_SH = (n13_NH - n13_SH) / tau_ex
            d13C_src_NH[j] = (n13_NH1 - n13_NH + n13_NH * a13_NH / tau_NH[j] - ex13_NH) / S_NH[j]
            d13C_src_SH[j] = (n13_SH1 - n13_SH + n13_SH * a13_SH / tau_SH[j] - ex13_SH) / S_SH[j]
        
        sigs = sample_source_signatures(rng, data, k, n)
        f13_bb = delta_to_fraction_d13C(sigs['bb_d13C'])
        f13_ff = delta_to_fraction_d13C(sigs['ff_d13C'])
        f13_mic = delta_to_fraction_d13C(sigs['mic_d13C'])
        
        if mode == "dual":
            # Also compute δD
            dD_glob_MC = sample_atm_dD(data, k, n)
            dD_NH_MC = dD_glob_MC - DD_IH_OFFSET
            dD_SH_MC = dD_glob_MC + DD_IH_OFFSET
            fD_NH_atm = delta_to_fraction_dD(dD_NH_MC)
            fD_SH_atm = delta_to_fraction_dD(dD_SH_MC)
            
            dD_src_NH = np.zeros(n); dD_src_SH = np.zeros(n)
            for j in range(n):
                nD_NH = fD_NH_atm[j] * CH4_NH[j] * PT_HEMI
                nD_NH1 = fD_NH_atm[j+1] * CH4_NH[j+1] * PT_HEMI
                nD_SH = fD_SH_atm[j] * CH4_SH[j] * PT_HEMI
                nD_SH1 = fD_SH_atm[j+1] * CH4_SH[j+1] * PT_HEMI
                exD_NH = (nD_SH - nD_NH) / tau_ex
                exD_SH = (nD_NH - nD_SH) / tau_ex
                dD_src_NH[j] = (nD_NH1 - nD_NH + nD_NH * aD_NH / tau_NH[j] - exD_NH) / S_NH[j]
                dD_src_SH[j] = (nD_SH1 - nD_SH + nD_SH * aD_SH / tau_SH[j] - exD_SH) / S_SH[j]
            
            fD_bb = delta_to_fraction_dD(sigs['bb_dD'])
            fD_ff = delta_to_fraction_dD(sigs['ff_dD'])
            fD_mic = delta_to_fraction_dD(sigs['mic_dD'])
            
            for j in range(n):
                A = np.array([
                    [1.0, 1.0, 1.0],
                    [f13_bb[j], f13_ff[j], f13_mic[j]],
                    [fD_bb[j], fD_ff[j], fD_mic[j]],
                ])
                # NH
                B_nh = np.array([S_NH[j], S_NH[j]*d13C_src_NH[j], S_NH[j]*dD_src_NH[j]])
                A_w = W_NH @ A; B_w = W_NH @ B_nh
                try:
                    res = lsq_linear(A_w, B_w, bounds=(0, S_NH[j]*1.5))
                    x = res.x
                except:
                    x = np.array([np.nan, np.nan, np.nan])
                BB_NH_comp[j,k] = x[0]; FF_NH_comp[j,k] = x[1]; Mic_NH_comp[j,k] = x[2]
                
                # SH
                B_sh = np.array([S_SH[j], S_SH[j]*d13C_src_SH[j], S_SH[j]*dD_src_SH[j]])
                A_w = W_SH @ A; B_w = W_SH @ B_sh
                try:
                    res = lsq_linear(A_w, B_w, bounds=(0, S_SH[j]*1.5))
                    x = res.x
                except:
                    x = np.array([np.nan, np.nan, np.nan])
                BB_SH_comp[j,k] = x[0]; FF_SH_comp[j,k] = x[1]; Mic_SH_comp[j,k] = x[2]
        
        else:  # d13C_only — 2×2 per hemisphere with BB fixed
            for j in range(n):
                # NH
                S_rem_nh = S_NH[j] - BB_hemi_NH
                rhs_nh = S_NH[j] * d13C_src_NH[j] - BB_hemi_NH * f13_bb[j]
                denom = f13_ff[j] - f13_mic[j]
                if abs(denom) > 1e-15:
                    FF_val = (rhs_nh - S_rem_nh * f13_mic[j]) / denom
                    Mic_val = S_rem_nh - FF_val
                else:
                    FF_val = np.nan; Mic_val = np.nan
                BB_NH_comp[j,k] = BB_hemi_NH; FF_NH_comp[j,k] = FF_val; Mic_NH_comp[j,k] = Mic_val
                
                # SH
                S_rem_sh = S_SH[j] - BB_hemi_SH
                rhs_sh = S_SH[j] * d13C_src_SH[j] - BB_hemi_SH * f13_bb[j]
                if abs(denom) > 1e-15:
                    FF_val = (rhs_sh - S_rem_sh * f13_mic[j]) / denom
                    Mic_val = S_rem_sh - FF_val
                else:
                    FF_val = np.nan; Mic_val = np.nan
                BB_SH_comp[j,k] = BB_hemi_SH; FF_SH_comp[j,k] = FF_val; Mic_SH_comp[j,k] = Mic_val
    
    print("  Done!")
    # Global = NH + SH
    FF_G = FF_NH_comp + FF_SH_comp
    Mic_G = Mic_NH_comp + Mic_SH_comp
    BB_G = BB_NH_comp + BB_SH_comp
    
    return years, {
        'FF_NH': FF_NH_comp, 'Mic_NH': Mic_NH_comp, 'BB_NH': BB_NH_comp,
        'FF_SH': FF_SH_comp, 'Mic_SH': Mic_SH_comp, 'BB_SH': BB_SH_comp,
        'FF_G': FF_G, 'Mic_G': Mic_G, 'BB_G': BB_G,
        'S_NH': S_NH, 'S_SH': S_SH,
    }


def main():
    N = 1000
    SEED = 42
    
    # === ONE-BOX ===
    print("\n" + "█"*70)
    print("  PHASE 1: BASELINE RUNS")
    print("█"*70)
    
    # 1-box dual
    yrs, BB_1d, FF_1d, Mic_1d, SS, tau = run_one_box("dual", N, SEED)
    np.savez(OUT_DIR / "onebox_dual.npz",
             years=yrs, BB=BB_1d, FF=FF_1d, Mic=Mic_1d,
             SumSource=SS, tau=tau)
    
    # 1-box d13C-only
    yrs, BB_1c, FF_1c, Mic_1c, SS, tau = run_one_box("d13C_only", N, SEED)
    np.savez(OUT_DIR / "onebox_d13C_only.npz",
             years=yrs, BB=BB_1c, FF=FF_1c, Mic=Mic_1c,
             SumSource=SS, tau=tau)
    
    # === TWO-BOX ===
    # 2-box dual
    yrs, results_2d = run_two_box("dual", N, SEED)
    np.savez(OUT_DIR / "twobox_dual.npz", years=yrs, **results_2d)
    
    # 2-box d13C-only
    yrs, results_2c = run_two_box("d13C_only", N, SEED)
    np.savez(OUT_DIR / "twobox_d13C_only.npz", years=yrs, **results_2c)
    
    # === SUMMARY STATISTICS ===
    summary = {}
    for name, FF_arr in [("onebox_dual", FF_1d), ("onebox_d13C_only", FF_1c),
                          ("twobox_dual", results_2d['FF_G']), ("twobox_d13C_only", results_2c['FF_G'])]:
        FF_s = smooth_5yr(FF_arr)
        ci_width = np.nanpercentile(FF_s, 95, axis=1) - np.nanpercentile(FF_s, 5, axis=1)
        summary[name] = {
            'FF_mean': float(np.nanmean(FF_s)),
            'FF_95CI_width_mean': float(np.nanmean(ci_width)),
            'FF_95CI_width_2007_2022': float(np.nanmean(ci_width[8:])),  # post-2007
        }
    
    # Compute constraint improvement
    for box in ['onebox', 'twobox']:
        dual_w = summary[f'{box}_dual']['FF_95CI_width_mean']
        c13_w = summary[f'{box}_d13C_only']['FF_95CI_width_mean']
        improvement = (c13_w - dual_w) / c13_w * 100
        summary[f'{box}_improvement_pct'] = round(improvement, 1)
    
    with open(OUT_DIR / "summary.json", 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\n{'='*60}")
    print("PHASE 1 SUMMARY")
    print(f"{'='*60}")
    for key, val in summary.items():
        print(f"  {key}: {val}")
    print(f"\nResults saved to: {OUT_DIR}")


if __name__ == "__main__":
    main()
