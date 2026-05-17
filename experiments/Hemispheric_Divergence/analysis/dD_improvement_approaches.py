#!/usr/bin/env python3
"""
Test 5 δD improvement approaches on the v3 delta-space model.

A. Constrain Mic δD with source-water δD (GNIP/OIPC proxy)
B. EDGAR subcategory-weighted FF δD
C. C3/C4-dependent BB δD (extend Luo 2024 approach)
D. Bayesian inversion with informative priors from A+B+C
E. Add observed NH-SH δD gradient as 7th constraint

Each approach modifies how source δD signatures are sampled,
then runs the full 2-box + 1-box v3 model (1000 MC).
"""

import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear, minimize
from scipy import stats as sp_stats
import json

from common import (
    ModelConfig, load_data,
    sample_KIE, compute_bulk_KIE, compute_lifetime,
    delta_to_fraction_d13C, delta_to_fraction_dD,
    sample_source_signatures, sample_source_signatures_hemi,
    sample_atm_d13C, sample_atm_dD, sample_atm_dD_hemi,
    SINK_FRACTIONS_NH, SINK_FRACTIONS_SH, SINK_FRACTIONS_GLOBAL,
    PT, PT_HEMI, LIFETIME_RATIO_NH, LIFETIME_RATIO_SH,
    TAU_EX_MEAN, TAU_EX_STD,
)

EXP_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = EXP_DIR / "results" / "dD_improvements"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR = ROOT / "rel" / "data"

NI = 1000
TREND_START = 2007
SEED = 42
N_YEARS = 23  # 1999-2021


def realistic_IH_gradient(years):
    anchor_years = np.array([2000, 2010, 2020, 2022])
    anchor_grad = np.array([108.0, 120.0, 140.0, 145.0])
    return np.interp(years, anchor_years, anchor_grad)


def fraction_to_delta_d13C(f):
    R_std = 0.0112372
    R = f / (1.0 - f)
    return (R / R_std - 1.0) * 1000.0


def fraction_to_delta_dD(f):
    R_std = 0.00015576
    R = f / (1.0 - f)
    return (R / R_std - 1.0) * 1000.0


def solve_delta_space(S, d13C_src_delta, dD_src_delta,
                      d13C_FF, d13C_Mic, d13C_BB,
                      dD_FF, dD_Mic, dD_BB):
    A = np.array([
        [1.0, 1.0, 1.0],
        [d13C_BB, d13C_FF, d13C_Mic],
        [dD_BB, dD_FF, dD_Mic],
    ])
    b = np.array([1.0, d13C_src_delta, dD_src_delta])
    scale = np.array([1.0, 1.0/50.0, 1.0/250.0])
    A_s = A * scale[:, None]
    b_s = b * scale
    try:
        res = lsq_linear(A_s, b_s, bounds=(0.0, 1.0))
        fracs = res.x
        if fracs.sum() > 0:
            fracs /= fracs.sum()
        return fracs[0]*S, fracs[1]*S, fracs[2]*S  # BB, FF, Mic
    except:
        return np.nan, np.nan, np.nan


def solve_delta_space_with_gradient_constraint(
        S_NH, d13C_src_NH, dD_src_NH,
        S_SH, d13C_src_SH, dD_src_SH,
        d13C_FF_NH, d13C_Mic_NH, d13C_BB_NH,
        dD_FF_NH, dD_Mic_NH, dD_BB_NH,
        d13C_FF_SH, d13C_Mic_SH, d13C_BB_SH,
        dD_FF_SH, dD_Mic_SH, dD_BB_SH,
        obs_dD_gradient=-14.5, gradient_weight=0.5):
    """
    Approach E: Solve NH+SH jointly with observed δD gradient as 7th equation.
    
    6 unknowns: f_BB_NH, f_FF_NH, f_Mic_NH, f_BB_SH, f_FF_SH, f_Mic_SH
    7 equations:
      NH: f_BB+f_FF+f_Mic = 1, δ¹³C mixing, δD mixing  (3)
      SH: f_BB+f_FF+f_Mic = 1, δ¹³C mixing, δD mixing  (3)
      Gradient: predicted_atm_dD_NH - predicted_atm_dD_SH ≈ obs  (1)
    
    The 7th equation uses the steady-state relationship:
      Δ(δD_atm) ≈ transfer_ratio × Δ(δD_src_weighted)
    where transfer_ratio ≈ 0.71 (from τ_ex, τ_sink, α_D)
    """
    TRANSFER_RATIO = 0.71  # from our steady-state calculation
    
    # Build combined system
    # x = [f_BB_NH, f_FF_NH, f_Mic_NH, f_BB_SH, f_FF_SH, f_Mic_SH]
    A = np.zeros((7, 6))
    b = np.zeros(7)
    
    # NH mass: row 0
    A[0, 0:3] = [1, 1, 1]
    b[0] = 1.0
    # NH d13C: row 1
    A[1, 0:3] = [d13C_BB_NH, d13C_FF_NH, d13C_Mic_NH]
    b[1] = d13C_src_NH
    # NH dD: row 2
    A[2, 0:3] = [dD_BB_NH, dD_FF_NH, dD_Mic_NH]
    b[2] = dD_src_NH
    # SH mass: row 3
    A[3, 3:6] = [1, 1, 1]
    b[3] = 1.0
    # SH d13C: row 4
    A[4, 3:6] = [d13C_BB_SH, d13C_FF_SH, d13C_Mic_SH]
    b[4] = d13C_src_SH
    # SH dD: row 5
    A[5, 3:6] = [dD_BB_SH, dD_FF_SH, dD_Mic_SH]
    b[5] = dD_src_SH
    # Gradient: row 6
    # predicted_atm_gradient = TRANSFER_RATIO * (src_dD_NH - src_dD_SH)
    # src_dD_NH = sum(f_i * dD_i_NH), src_dD_SH = sum(f_j * dD_j_SH)
    A[6, 0:3] = [dD_BB_NH, dD_FF_NH, dD_Mic_NH]
    A[6, 3:6] = [-dD_BB_SH, -dD_FF_SH, -dD_Mic_SH]
    A[6, :] *= TRANSFER_RATIO
    b[6] = obs_dD_gradient
    
    # Scale rows
    scale = np.array([1.0, 1/50, 1/250, 1.0, 1/50, 1/250, gradient_weight/15])
    A_s = A * scale[:, None]
    b_s = b * scale
    
    try:
        res = lsq_linear(A_s, b_s, bounds=(0.0, 1.0))
        x = res.x
        f_NH = x[0:3]; f_SH = x[3:6]
        if f_NH.sum() > 0: f_NH /= f_NH.sum()
        if f_SH.sum() > 0: f_SH /= f_SH.sum()
        return (f_NH[0]*S_NH, f_NH[1]*S_NH, f_NH[2]*S_NH,
                f_SH[0]*S_SH, f_SH[1]*S_SH, f_SH[2]*S_SH)
    except:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan


# ═══════════════════════════════════════════════════════════
# δD SIGNATURE MODIFICATION APPROACHES
# ═══════════════════════════════════════════════════════════

def load_original_sigs():
    """Load the original MC source δD signatures."""
    sigs = {}
    for src in ['FF', 'Mic', 'BB']:
        for hemi in ['NH', 'SH']:
            df = pd.read_csv(DATA_DIR / f'{src}_dD_{hemi}_MC.csv')
            sigs[f'{src}_dD_{hemi}'] = df.values  # (24, 1001)
    return sigs


def approach_A_source_water(rng, orig_sigs, k, n):
    """
    A. Constrain Mic δD with source-water δD maps.
    
    Microbial CH₄ δD is controlled by:
    δD_CH4 ≈ δD_water - α_methanogenesis (where α ≈ 1.15-1.35, i.e. -150 to -350‰ offset)
    
    Source water δD (GNIP/OIPC data):
      NH mean: ~-50‰ (area-weighted, wetland regions)
      SH mean: ~-25‰ (warmer, less continental)
      
    Methanogenesis fractionation: -270 ± 20‰ (Waldron 1999, Whiticar 1999)
    
    So: Mic δD ≈ water_δD - 270 ± 20‰
      NH: -50 - 270 ± 20 = -320 ± 20‰ (was -317 ± 74)
      SH: -25 - 270 ± 20 = -295 ± 20‰ (was -305 ± 74)
    """
    sigs = {}
    for src in ['FF', 'BB']:
        for hemi in ['NH', 'SH']:
            key = f'{src}_dD_{hemi}'
            sigs[key] = orig_sigs[key][:n, k % orig_sigs[key].shape[1]]
    
    # Constrained Mic δD
    sigs['Mic_dD_NH'] = rng.normal(-320.0, 20.0, n)
    sigs['Mic_dD_SH'] = rng.normal(-295.0, 20.0, n)
    
    return sigs


def approach_B_edgar_ff(rng, orig_sigs, k, n):
    """
    B. EDGAR subcategory-weighted FF δD.
    
    FF subcategories and their δD (Sherwood et al. 2017):
      Coal:            -160 ± 15‰  (NH: 45%, SH: 30%)
      Conventional gas: -175 ± 15‰  (NH: 30%, SH: 35%)  
      Shale/tight gas: -210 ± 20‰  (NH: 15%, SH: 5%)
      Oil production:  -185 ± 20‰  (NH: 10%, SH: 30%)
    
    EDGAR v7 subcategory fractions (approximate):
    """
    sigs = {}
    for src in ['Mic', 'BB']:
        for hemi in ['NH', 'SH']:
            key = f'{src}_dD_{hemi}'
            sigs[key] = orig_sigs[key][:n, k % orig_sigs[key].shape[1]]
    
    # EDGAR-weighted FF δD
    # NH: heavy coal + shale
    coal_nh = 0.45; gas_nh = 0.30; shale_nh = 0.15; oil_nh = 0.10
    ff_dD_nh_mean = coal_nh*(-160) + gas_nh*(-175) + shale_nh*(-210) + oil_nh*(-185)  # = -175.5
    ff_dD_nh_std = np.sqrt(coal_nh*15**2 + gas_nh*15**2 + shale_nh*20**2 + oil_nh*20**2)  # ~16
    
    # SH: more oil, less shale
    coal_sh = 0.30; gas_sh = 0.35; shale_sh = 0.05; oil_sh = 0.30
    ff_dD_sh_mean = coal_sh*(-160) + gas_sh*(-175) + shale_sh*(-210) + oil_sh*(-185)  # = -176.8
    ff_dD_sh_std = np.sqrt(coal_sh*15**2 + gas_sh*15**2 + shale_sh*20**2 + oil_sh*20**2)  # ~16
    
    sigs['FF_dD_NH'] = rng.normal(ff_dD_nh_mean, ff_dD_nh_std, n)
    sigs['FF_dD_SH'] = rng.normal(ff_dD_sh_mean, ff_dD_sh_std, n)
    
    return sigs


def approach_C_c3c4_bb(rng, orig_sigs, k, n):
    """
    C. C3/C4-dependent BB δD.
    
    BB δD depends on burned biomass type (Snover & Quay 2000):
      C3 (forests):  -215 ± 15‰
      C4 (savanna):  -185 ± 15‰
    
    Luo 2024 C4 fractions (already in your d13C BB data):
      NH: ~30% C4 (tropical savannas)
      SH: ~55% C4 (African/Australian savannas dominate)
    
    These fractions are time-varying (fires shift location year to year).
    """
    sigs = {}
    for src in ['FF', 'Mic']:
        for hemi in ['NH', 'SH']:
            key = f'{src}_dD_{hemi}'
            sigs[key] = orig_sigs[key][:n, k % orig_sigs[key].shape[1]]
    
    # C3/C4 weighted BB δD
    # NH: ~30% C4
    c4_frac_nh = rng.normal(0.30, 0.05, n).clip(0.05, 0.95)
    bb_nh = c4_frac_nh * (-185) + (1 - c4_frac_nh) * (-215) + rng.normal(0, 10, n)
    
    # SH: ~55% C4
    c4_frac_sh = rng.normal(0.55, 0.08, n).clip(0.05, 0.95)
    bb_sh = c4_frac_sh * (-185) + (1 - c4_frac_sh) * (-215) + rng.normal(0, 10, n)
    
    sigs['BB_dD_NH'] = bb_nh
    sigs['BB_dD_SH'] = bb_sh
    
    return sigs


def approach_D_bayesian(rng, orig_sigs, k, n):
    """
    D. Bayesian: combine A+B+C (all informative priors).
    """
    sigs = {}
    
    # Mic from approach A (source-water constrained)
    sigs['Mic_dD_NH'] = rng.normal(-320.0, 20.0, n)
    sigs['Mic_dD_SH'] = rng.normal(-295.0, 20.0, n)
    
    # FF from approach B (EDGAR-weighted)
    sigs['FF_dD_NH'] = rng.normal(-175.5, 16.0, n)
    sigs['FF_dD_SH'] = rng.normal(-176.8, 16.0, n)
    
    # BB from approach C (C3/C4-weighted)
    c4_nh = rng.normal(0.30, 0.05, n).clip(0.05, 0.95)
    sigs['BB_dD_NH'] = c4_nh*(-185) + (1-c4_nh)*(-215) + rng.normal(0, 10, n)
    c4_sh = rng.normal(0.55, 0.08, n).clip(0.05, 0.95)
    sigs['BB_dD_SH'] = c4_sh*(-185) + (1-c4_sh)*(-215) + rng.normal(0, 10, n)
    
    return sigs


# ═══════════════════════════════════════════════════════════
# MODEL RUNNER
# ═══════════════════════════════════════════════════════════

def run_model(data, approach_name, sig_modifier=None, use_gradient_constraint=False):
    """Run 2-box + 1-box with specified δD approach."""
    print(f"\n{'='*60}")
    print(f"  Approach {approach_name}")
    print(f"{'='*60}")
    
    n = data.n_years
    years = data.model_years
    rng = np.random.default_rng(SEED)
    
    all_years = np.arange(years[0], years[-1] + 2)
    IH_grad = realistic_IH_gradient(all_years.astype(float))
    CH4_NH = data.CH4_global + IH_grad / 2.0
    CH4_SH = data.CH4_global - IH_grad / 2.0
    
    tau_global = compute_lifetime(years, "varying", 9.0)
    tau_NH = tau_global * LIFETIME_RATIO_NH
    tau_SH = tau_global * LIFETIME_RATIO_SH
    
    orig_sigs = load_original_sigs()
    
    FF_NH = np.zeros((n, NI)); FF_SH = np.zeros((n, NI))
    Mic_NH = np.zeros((n, NI)); Mic_SH = np.zeros((n, NI))
    BB_NH = np.zeros((n, NI)); BB_SH = np.zeros((n, NI))
    
    for k in range(NI):
        if (k+1) % 250 == 0:
            print(f"    iter {k+1}/{NI}")
        
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
            S_NH[i] = (M_NH1-M_NH) + M_NH/tau_NH[i] - (M_SH-M_NH)/tau_ex
            S_SH[i] = (M_SH1-M_SH) + M_SH/tau_SH[i] - (M_NH-M_SH)/tau_ex
        
        d13C_MC = sample_atm_d13C(data, k, n)
        dD_NH_MC, dD_SH_MC = sample_atm_dD_hemi(data, k, n)
        nc = min(len(data.c13_global), n+1)
        d13C_off = d13C_MC[:nc] - data.c13_global[:nc]
        d13C_NH_MC = data.c13_NH[:nc] + d13C_off
        d13C_SH_MC = data.c13_SH[:nc] + d13C_off
        
        f13_NH = delta_to_fraction_d13C(d13C_NH_MC)
        f13_SH = delta_to_fraction_d13C(d13C_SH_MC)
        fD_NH = delta_to_fraction_dD(dD_NH_MC)
        fD_SH = delta_to_fraction_dD(dD_SH_MC)
        
        # Source isotopic compositions
        d13C_src_frac_NH = np.zeros(n); d13C_src_frac_SH = np.zeros(n)
        dD_src_frac_NH = np.zeros(n); dD_src_frac_SH = np.zeros(n)
        
        for j in range(n):
            M_NH_ = CH4_NH[j]*PT_HEMI; M_NH1_ = CH4_NH[j+1]*PT_HEMI
            M_SH_ = CH4_SH[j]*PT_HEMI; M_SH1_ = CH4_SH[j+1]*PT_HEMI
            
            n13 = f13_NH[j]*M_NH_; n13_1 = f13_NH[j+1]*M_NH1_
            n13_s = f13_SH[j]*M_SH_
            d13C_src_frac_NH[j] = (n13_1-n13 + n13*a13_NH/tau_NH[j] - (n13_s-n13)/tau_ex) / S_NH[j]
            
            n13 = f13_SH[j]*M_SH_; n13_1 = f13_SH[j+1]*M_SH1_
            n13_n = f13_NH[j]*M_NH_
            d13C_src_frac_SH[j] = (n13_1-n13 + n13*a13_SH/tau_SH[j] - (n13_n-n13)/tau_ex) / S_SH[j]
            
            nD = fD_NH[j]*M_NH_; nD_1 = fD_NH[j+1]*M_NH1_
            nD_s = fD_SH[j]*M_SH_
            dD_src_frac_NH[j] = (nD_1-nD + nD*aD_NH/tau_NH[j] - (nD_s-nD)/tau_ex) / S_NH[j]
            
            nD = fD_SH[j]*M_SH_; nD_1 = fD_SH[j+1]*M_SH1_
            nD_n = fD_NH[j]*M_NH_
            dD_src_frac_SH[j] = (nD_1-nD + nD*aD_SH/tau_SH[j] - (nD_n-nD)/tau_ex) / S_SH[j]
        
        d13C_src_NH = fraction_to_delta_d13C(d13C_src_frac_NH)
        d13C_src_SH = fraction_to_delta_d13C(d13C_src_frac_SH)
        dD_src_NH = fraction_to_delta_dD(dD_src_frac_NH)
        dD_src_SH = fraction_to_delta_dD(dD_src_frac_SH)
        
        # Get source signatures (original for d13C, modified for dD)
        sigs_orig = sample_source_signatures_hemi(rng, data, k, n)
        
        if sig_modifier is not None:
            mod_sigs = sig_modifier(rng, orig_sigs, k, n)
        else:
            mod_sigs = None
        
        for j in range(n):
            # d13C always from original
            d13C_ff_nh = sigs_orig['ff_d13C_NH'][j]
            d13C_mic_nh = sigs_orig['mic_d13C_NH'][j]
            d13C_bb_nh = sigs_orig['bb_d13C_NH'][j]
            d13C_ff_sh = sigs_orig['ff_d13C_SH'][j]
            d13C_mic_sh = sigs_orig['mic_d13C_SH'][j]
            d13C_bb_sh = sigs_orig['bb_d13C_SH'][j]
            
            # dD: modified or original
            if mod_sigs is not None:
                dD_ff_nh = mod_sigs['FF_dD_NH'][j]
                dD_mic_nh = mod_sigs['Mic_dD_NH'][j]
                dD_bb_nh = mod_sigs['BB_dD_NH'][j]
                dD_ff_sh = mod_sigs['FF_dD_SH'][j]
                dD_mic_sh = mod_sigs['Mic_dD_SH'][j]
                dD_bb_sh = mod_sigs['BB_dD_SH'][j]
            else:
                dD_ff_nh = sigs_orig['ff_dD_NH'][j]
                dD_mic_nh = sigs_orig['mic_dD_NH'][j]
                dD_bb_nh = sigs_orig['bb_dD_NH'][j]
                dD_ff_sh = sigs_orig['ff_dD_SH'][j]
                dD_mic_sh = sigs_orig['mic_dD_SH'][j]
                dD_bb_sh = sigs_orig['bb_dD_SH'][j]
            
            if use_gradient_constraint:
                # Approach E: joint solve with gradient constraint
                result = solve_delta_space_with_gradient_constraint(
                    S_NH[j], d13C_src_NH[j], dD_src_NH[j],
                    S_SH[j], d13C_src_SH[j], dD_src_SH[j],
                    d13C_ff_nh, d13C_mic_nh, d13C_bb_nh,
                    dD_ff_nh, dD_mic_nh, dD_bb_nh,
                    d13C_ff_sh, d13C_mic_sh, d13C_bb_sh,
                    dD_ff_sh, dD_mic_sh, dD_bb_sh)
                BB_NH[j,k], FF_NH[j,k], Mic_NH[j,k] = result[0], result[1], result[2]
                BB_SH[j,k], FF_SH[j,k], Mic_SH[j,k] = result[3], result[4], result[5]
            else:
                BB_NH[j,k], FF_NH[j,k], Mic_NH[j,k] = solve_delta_space(
                    S_NH[j], d13C_src_NH[j], dD_src_NH[j],
                    d13C_ff_nh, d13C_mic_nh, d13C_bb_nh,
                    dD_ff_nh, dD_mic_nh, dD_bb_nh)
                BB_SH[j,k], FF_SH[j,k], Mic_SH[j,k] = solve_delta_space(
                    S_SH[j], d13C_src_SH[j], dD_src_SH[j],
                    d13C_ff_sh, d13C_mic_sh, d13C_bb_sh,
                    dD_ff_sh, dD_mic_sh, dD_bb_sh)
    
    return {
        'years': years,
        'FF_NH': FF_NH, 'FF_SH': FF_SH,
        'Mic_NH': Mic_NH, 'Mic_SH': Mic_SH,
        'BB_NH': BB_NH, 'BB_SH': BB_SH,
    }


def analyze(res, label):
    """Compute summary statistics for a model run."""
    years = res['years']
    n = len(years)
    
    FF_g = res['FF_NH'] + res['FF_SH']
    Mic_g = res['Mic_NH'] + res['Mic_SH']
    BB_g = res['BB_NH'] + res['BB_SH']
    
    j2010 = np.where(years == 2010)[0][0]
    
    # Medians at 2010
    ff = np.nanmedian(FF_g[j2010,:])
    mic = np.nanmedian(Mic_g[j2010,:])
    bb = np.nanmedian(BB_g[j2010,:])
    total = ff + mic + bb
    
    # Trends (trim last year)
    def trend(arr):
        end = years[-1] - 1
        mask = (years >= TREND_START) & (years <= end)
        yrs = years[mask]
        sub = arr[mask, :]
        slopes = np.array([sp_stats.linregress(yrs, sub[:,k]).slope 
                           for k in range(sub.shape[1]) 
                           if not np.any(np.isnan(sub[:,k]))])
        return slopes
    
    ff_slopes = trend(FF_g)
    mic_slopes = trend(Mic_g)
    bb_slopes = trend(BB_g)
    
    # NH FF share
    nh_share = np.nanmedian(res['FF_NH'][j2010,:]) / ff if ff > 0 else 0
    
    # Uncertainty width (90% CI width)
    ff_ci = np.nanpercentile(FF_g[j2010,:], 95) - np.nanpercentile(FF_g[j2010,:], 5)
    bb_ci = np.nanpercentile(BB_g[j2010,:], 95) - np.nanpercentile(BB_g[j2010,:], 5)
    
    result = {
        'label': label,
        'ff_2010': float(ff),
        'mic_2010': float(mic),
        'bb_2010': float(bb),
        'total_2010': float(total),
        'ff_pct': float(ff/total*100) if total > 0 else 0,
        'bb_pct': float(bb/total*100) if total > 0 else 0,
        'nh_ff_share': float(nh_share),
        'ff_trend': float(np.median(ff_slopes)) if len(ff_slopes) > 0 else None,
        'ff_trend_90ci': [float(np.percentile(ff_slopes, 5)), float(np.percentile(ff_slopes, 95))] if len(ff_slopes) > 0 else None,
        'mic_trend': float(np.median(mic_slopes)) if len(mic_slopes) > 0 else None,
        'bb_trend': float(np.median(bb_slopes)) if len(bb_slopes) > 0 else None,
        'ff_90ci_width': float(ff_ci),
        'bb_90ci_width': float(bb_ci),
        'ff_trend_significant': bool(np.percentile(ff_slopes, 5) > 0 or np.percentile(ff_slopes, 95) < 0) if len(ff_slopes) > 0 else False,
        'mic_trend_significant': bool(np.percentile(mic_slopes, 5) > 0 or np.percentile(mic_slopes, 95) < 0) if len(mic_slopes) > 0 else False,
    }
    
    print(f"\n  {label}:")
    print(f"    FF={ff:.0f} Tg/yr ({ff/total*100:.0f}%), Mic={mic:.0f}, BB={bb:.0f}")
    print(f"    FF trend: {result['ff_trend']:+.2f} [{result['ff_trend_90ci'][0]:+.2f}, {result['ff_trend_90ci'][1]:+.2f}]")
    print(f"    Mic trend: {result['mic_trend']:+.2f}{'✓' if result['mic_trend_significant'] else ''}")
    print(f"    BB trend: {result['bb_trend']:+.2f}")
    print(f"    NH FF share: {nh_share*100:.0f}%")
    print(f"    FF 90% CI width: {ff_ci:.0f} Tg/yr")
    print(f"    BB 90% CI width: {bb_ci:.0f} Tg/yr")
    
    return result


def main():
    print("=" * 70)
    print("δD IMPROVEMENT APPROACHES: Testing A through E")
    print("=" * 70)
    
    data = load_data(ROOT, two_box=True)
    
    approaches = [
        ("Baseline (v3)", None, False),
        ("A: Source-water Mic δD", approach_A_source_water, False),
        ("B: EDGAR-weighted FF δD", approach_B_edgar_ff, False),
        ("C: C3/C4 BB δD", approach_C_c3c4_bb, False),
        ("D: Bayesian (A+B+C combined)", approach_D_bayesian, False),
        ("E: δD gradient constraint", None, True),
    ]
    
    all_results = []
    
    for name, modifier, use_grad in approaches:
        res = run_model(data, name, sig_modifier=modifier, use_gradient_constraint=use_grad)
        stats = analyze(res, name)
        all_results.append(stats)
        
        # Save per-approach results
        safe_name = name.split(":")[0].strip().replace(" ", "_").replace("(", "").replace(")", "")
        np.savez(RESULTS_DIR / f"{safe_name}.npz",
                 NH_FF=res['FF_NH'], SH_FF=res['FF_SH'],
                 NH_Mic=res['Mic_NH'], SH_Mic=res['Mic_SH'],
                 NH_BB=res['BB_NH'], SH_BB=res['BB_SH'],
                 years=res['years'])
    
    # Summary table
    print("\n" + "=" * 100)
    print("SUMMARY COMPARISON")
    print("=" * 100)
    print(f"{'Approach':35s} {'FF':>6s} {'Mic':>6s} {'BB':>6s} {'FF%':>5s} {'BB%':>5s} {'NH_FF%':>7s} {'FF_trend':>10s} {'FF_CI':>7s} {'BB_CI':>7s}")
    print("-" * 100)
    for r in all_results:
        ft = f"{r['ff_trend']:+.2f}" if r['ff_trend'] is not None else "N/A"
        print(f"{r['label']:35s} {r['ff_2010']:6.0f} {r['mic_2010']:6.0f} {r['bb_2010']:6.0f} "
              f"{r['ff_pct']:5.0f} {r['bb_pct']:5.0f} {r['nh_ff_share']*100:7.0f} "
              f"{ft:>10s} {r['ff_90ci_width']:7.0f} {r['bb_90ci_width']:7.0f}")
    
    print(f"\n  EDGAR reference: FF=110, Mic=370, BB=30, FF%=19, NH_FF=72%")
    
    # Save summary
    with open(RESULTS_DIR / "comparison_summary.json", 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    
    print(f"\n  All results saved to {RESULTS_DIR}/")


if __name__ == "__main__":
    main()
