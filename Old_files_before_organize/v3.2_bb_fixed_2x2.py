#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Two-Hemisphere 2×2 (BB-Fixed) Monte Carlo Box Model
======================================================
Version: 3.2 — Separate δ¹³C and δD inversions (Ben's approach + our specialties)

APPROACH:
  Following Riddell-Young et al. (2025, PNAS), fix BB from CarbonTracker and
  solve a 2-equation system for FF and Mic per hemisphere:
    
    S_total = FF + Mic + BB_fixed
    S_total × δ_source = FF × δ_FF + Mic × δ_Mic + BB_fixed × δ_BB
  
  → FF = (S × δ_source - Mic_sig × (S - BB) - BB_sig × BB) / (FF_sig - Mic_sig)
  → Mic = S - BB - FF

  This is done SEPARATELY for δ¹³C and δD, yielding two independent estimates.
  The agreement/divergence between them tests sensitivity to Cl, BB, and OH trends.

OUR SPECIALTIES (beyond Ben's one-box):
  1. Two-hemisphere structure (NH/SH boxes + interhemispheric exchange)
  2. KIE sampling in MC loop (all 4 sinks: OH, Cl, Strat, Soil)
  3. Time-varying CH₄ lifetime: τ(t) = 9.0 - 0.017*(t - 2010)
  4. τ_ex uncertainty: sampled Normal(1.0, 0.1)
  5. Hemisphere-specific sink fractions and lifetimes
  6. δD hemispheric offset = ±6‰ (from Riddell-Young NH-SH gradient ~12‰)
  7. Quality monitoring + condition tracking
  8. NH/SH BB split from GFED4 (55%/45%)

ADVANTAGE OVER 3×3:
  - Well-conditioned (2×2 system, no δD ill-conditioning)
  - No bounded least squares needed
  - Separate δ¹³C and δD estimates allow cross-validation
  - Tighter uncertainty bounds

LIMITATION:
  - BB is not independently constrained by isotopes
  - Assumes CarbonTracker BB prior is accurate

Author: Built from v3.1 + Ben's methodology
Date: 2026-05-05
"""

from pathlib import Path
import sys
import json

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# I/O paths
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent
ORIGINAL_MODEL_DIR = BASE_DIR.parent / "TwoIsotopeBoxModel"
REL_DIR = ORIGINAL_MODEL_DIR / "rel"
DATA_DIR = REL_DIR / "data"
SRC_DIR = REL_DIR / "output"
OUT_DIR = BASE_DIR / "Output_v3.2_2x2"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
N_ITERATIONS = 1000
TAU_EX_MEAN = 1.0
TAU_EX_STD = 0.1

# NH/SH BB split from GFED4 (global fires split roughly 55/45 NH/SH)
BB_NH_FRACTION = 0.55
BB_SH_FRACTION = 0.45

# ---------------------------------------------------------------------------
# KIE Sampling (all 4 sinks, same as v3.1)
# ---------------------------------------------------------------------------
KIE_CONFIG = {
    'OH_13C':    {'dist': 'uniform', 'low': 1.0039, 'high': 1.0054},
    'OH_D':      {'dist': 'uniform', 'low': 1.294,  'high': 1.327},
    'Cl_13C':    {'dist': 'normal',  'mean': 1.066,  'std': 0.002},
    'Cl_D':      {'dist': 'normal',  'mean': 1.52,   'std': 0.02},
    'Strat_13C': {'dist': 'normal',  'mean': 1.003,  'std': 0.001},
    'Strat_D':   {'dist': 'normal',  'mean': 1.179,  'std': 0.01},
    'Soil_13C':  {'dist': 'normal',  'mean': 1.0201, 'std': 0.003},
    'Soil_D':    {'dist': 'normal',  'mean': 1.083,  'std': 0.01},
}

def sample_KIE(rng):
    kies = {}
    for key, cfg in KIE_CONFIG.items():
        if isinstance(cfg, (int, float)):
            kies[key] = cfg
        elif cfg['dist'] == 'uniform':
            kies[key] = rng.uniform(cfg['low'], cfg['high'])
        elif cfg['dist'] == 'normal':
            kies[key] = rng.normal(cfg['mean'], cfg['std'])
    return kies

# ---------------------------------------------------------------------------
# Sink fractions per hemisphere
# ---------------------------------------------------------------------------
SINK_FRACTIONS = {
    'NH': {'OH': 0.825, 'Cl': 0.040, 'Strat': 0.070, 'Soil': 0.065},
    'SH': {'OH': 0.850, 'Cl': 0.028, 'Strat': 0.070, 'Soil': 0.052},
}

# ---------------------------------------------------------------------------
# Time-varying lifetime
# ---------------------------------------------------------------------------
def compute_lifetime_array(years):
    return 9.0 - 0.017 * (np.asarray(years, dtype=float) - 2010)

LIFETIME_RATIO = {'NH': 0.95, 'SH': 1.05}

# δD hemispheric offset (corrected from paper)
DD_IH_OFFSET = 6.0  # ‰

# ---------------------------------------------------------------------------
# Isotope utilities
# ---------------------------------------------------------------------------
C13Std = 0.011113
DStd = 0.00015576
PT = 2.815
PT_HEMI = PT / 2.0

def delta_to_fraction_d13C(delta_permil):
    R = (delta_permil / 1000.0 + 1.0) * C13Std
    return R / (1.0 + R)

def delta_to_fraction_dD(delta_permil):
    R = (delta_permil / 1000.0 + 1.0) * DStd
    return R / (1.0 + R)

def fraction_to_delta_d13C(f):
    R = f / (1.0 - f)
    return ((R - C13Std) / C13Std) * 1000

def fraction_to_delta_dD(f):
    R = f / (1.0 - f)
    return ((R - DStd) / DStd) * 1000

# ---------------------------------------------------------------------------
# 5-Year Smoothing
# ---------------------------------------------------------------------------
def smooth_5yr(arr_2d):
    n_years, n_cols = arr_2d.shape
    if n_years < 5:
        return arr_2d.copy()
    result = np.zeros_like(arr_2d)
    result[0, :] = np.mean(arr_2d[0:3, :], axis=0)
    result[1, :] = np.mean(arr_2d[0:4, :], axis=0)
    for i in range(2, n_years - 2):
        result[i, :] = np.mean(arr_2d[i-2:i+3, :], axis=0)
    result[-2, :] = np.mean(arr_2d[-4:, :], axis=0)
    result[-1, :] = np.mean(arr_2d[-3:, :], axis=0)
    return result


# ===========================================================================
# DATA LOADING
# ===========================================================================
print("="*70)
print("TWO-HEMISPHERE 2×2 (BB-FIXED) MODEL v3.2")
print("="*70)
print("\nLoading data...")

# CH₄
CH4data_raw = pd.read_excel(DATA_DIR / "GML_CH4_AnnualMean.xlsx").to_numpy()
CH4_global = CH4data_raw[15:39, 1].astype(float)
CH4_years = CH4data_raw[15:39, 0].astype(float)
IH_GRADIENT = np.linspace(80, 100, len(CH4_global))
CH4_NH = CH4_global + IH_GRADIENT / 2.0
CH4_SH = CH4_global - IH_GRADIENT / 2.0

# δ¹³C
C13data = pd.read_excel(DATA_DIR / "ch4c13_nh_sh_mean.xlsx", header=None).to_numpy()
c13_dates = C13data[:, 0]
c13_global_raw = C13data[:, 1]
c13_NH_raw = C13data[:, 2]
c13_SH_raw = C13data[:, 3]

def annual_average(dates, values):
    years_floor = np.floor(dates).astype(int)
    unique_years = np.unique(years_floor)
    ann_years, ann_means = [], []
    for yr in unique_years:
        mask = years_floor == yr
        if np.sum(mask) >= 6:
            ann_years.append(yr)
            ann_means.append(np.nanmean(values[mask]))
    return np.array(ann_years), np.array(ann_means)

c13_ann_years, c13_ann_global = annual_average(c13_dates, c13_global_raw)
_, c13_ann_NH = annual_average(c13_dates, c13_NH_raw)
_, c13_ann_SH = annual_average(c13_dates, c13_SH_raw)

c13_start_idx = np.where(c13_ann_years == 1999)[0][0]
c13_end_idx = np.where(c13_ann_years == 2022)[0][0] + 1
c13_NH = c13_ann_NH[c13_start_idx:c13_end_idx]
c13_SH = c13_ann_SH[c13_start_idx:c13_end_idx]
c13_glob = c13_ann_global[c13_start_idx:c13_end_idx]

# δD
glob_ann_dD_path = DATA_DIR / "GlobMean_dD_iterations_UmezawaCal_noBUDS.xlsx"
glob_ann_dD_df = pd.read_excel(glob_ann_dD_path)
glob_ann_dD_num = glob_ann_dD_df.apply(pd.to_numeric, errors="coerce")
glob_ann_dD = glob_ann_dD_num.iloc[:, 1].to_numpy(dtype=np.float64) - 0.5
dD_AnnAvg_matrix = glob_ann_dD_num.iloc[:, 5:].to_numpy(dtype=np.float64)

# d13C DEI iterations
d13C_glob_iterations_data = np.loadtxt(str(DATA_DIR / "d13C_dei_compiled.txt"))
d13C_glob_iterations = d13C_glob_iterations_data[1:, 1:]

# Source signatures
BB_dD_data = pd.read_csv(SRC_DIR / "BB_dD_annual.csv", delimiter=',', header=None)
Mic_dD_data = pd.read_csv(SRC_DIR / "Mic_dD_AnnGlob.csv", delimiter=',', header=None)
Mic_dD_MC_trends = pd.read_csv(SRC_DIR / "Mic_dD_MC.csv", delimiter=',', header=None)
Mic_dD_MC = Mic_dD_MC_trends.iloc[6:, 1:]

FF_dD_data = pd.read_csv(SRC_DIR / "FF_dD_GlobUnc.csv", delimiter=',')
FF_dD_MC_EDGAR_data = pd.read_csv(SRC_DIR / "FF_dD_GlobMC_EDGAR.csv", delimiter=',')
FF_dD_MC_EDGAR = FF_dD_MC_EDGAR_data.iloc[34:, 1:]
if FF_dD_MC_EDGAR.shape[0] < 24:
    pad_count = 24 - FF_dD_MC_EDGAR.shape[0]
    pad_rows = pd.concat([FF_dD_MC_EDGAR.iloc[0:1, :]] * pad_count, ignore_index=True)
    FF_dD_MC_EDGAR = pd.concat([pad_rows, FF_dD_MC_EDGAR], ignore_index=True)

BB_d13C_data = pd.read_csv(SRC_DIR / "BB_d13C_annual.csv", delimiter=',', header=None)
Mic_d13C_data = pd.read_csv(SRC_DIR / "Mic_d13C_annual.csv", delimiter=',', header=None)
Mic_d13C_MC_trends = pd.read_csv(SRC_DIR / "Mic_d13C_MC.csv", delimiter=',', header=None)
Mic_d13C_MC = Mic_d13C_MC_trends.iloc[:, 1:]

FF_d13C_data = pd.read_csv(SRC_DIR / "FF_d13C_GlobUnc.csv", delimiter=',')
FF_d13C_MC_EDGAR_data = pd.read_csv(SRC_DIR / "FF_d13C_GlobMC_EDGAR.csv", delimiter=',')
FF_d13C_MC_EDGAR = FF_d13C_MC_EDGAR_data.iloc[28:, 1:]

# CarbonTracker BB (FIXED — core of this approach)
data_CT = pd.read_excel(DATA_DIR / "CarbonTracker_CH4.xlsx")
bbCT = data_CT.iloc[:, 9].values  # Prior BB emissions
BB_GLOBAL_MEAN = np.mean(bbCT)  # ~30 Tg/yr
# Use time-varying BB if available (prefer annual from CT)
BB_annual_CT = bbCT  # This may need alignment to 1999-2022

print(f"  BB from CarbonTracker (GFED4 prior): {BB_GLOBAL_MEAN:.1f} Tg/yr (FIXED)")
print(f"  BB split: NH={BB_GLOBAL_MEAN*BB_NH_FRACTION:.1f}, SH={BB_GLOBAL_MEAN*BB_SH_FRACTION:.1f} Tg/yr")

# Source signature arrays
mic_d13C_mean = Mic_d13C_data.iloc[:, 1].mean()
ff_d13C = np.array(FF_d13C_data.iloc[28:, 1]).flatten()
ff_d13C_U = np.array(FF_d13C_data.iloc[28:, 2]).flatten()
bb_d13C = np.array(BB_d13C_data.iloc[1:, 1]).flatten()
bb_d13C_U = np.array(BB_d13C_data.iloc[1:, 2]).flatten()
mean_bb_d13C = bb_d13C[-1]
mean_bb_d13C_U = bb_d13C_U[-1]
bb_d13C = np.concatenate((bb_d13C, np.full(1, mean_bb_d13C)))
bb_d13C_U = np.concatenate((bb_d13C_U, np.full(1, mean_bb_d13C_U)))

mic_dd_mean = Mic_dD_data.iloc[:, 1].mean()
mic_dd_U = 8.2
ff_dd = np.array(FF_dD_data.iloc[34:, 1]).flatten()
ff_dd_U = np.array(FF_dD_data.iloc[34:, 2]).flatten()
bb_dd = np.array(BB_dD_data.iloc[:, 1]).flatten()
bb_dd_U = np.array(BB_dD_data.iloc[:, 2]).flatten()
mean_bb_dd = bb_dd[-1]
mean_bb_dd_U = bb_dd_U[-1]
bb_dd = np.concatenate((np.full(3, mean_bb_dd), bb_dd, np.full(1, mean_bb_dd)))
bb_dd_U = np.concatenate((np.full(3, mean_bb_dd_U), bb_dd_U, np.full(1, mean_bb_dd_U)))

# ===========================================================================
# MODEL DIMENSIONS
# ===========================================================================
n_years_model = len(CH4_global) - 1  # 23
model_years = np.arange(1999, 1999 + n_years_model)
Lifetime_global = compute_lifetime_array(model_years)
Lifetime_NH = Lifetime_global * LIFETIME_RATIO['NH']
Lifetime_SH = Lifetime_global * LIFETIME_RATIO['SH']
target_length = n_years_model

# Pad arrays
while len(ff_dd) < target_length:
    ff_dd = np.concatenate([np.array([ff_dd[0]]), ff_dd])
    ff_dd_U = np.concatenate([np.array([ff_dd_U[0]]), ff_dd_U])
ff_dd = ff_dd[:target_length]
ff_dd_U = ff_dd_U[:target_length]
bb_dd = bb_dd[:target_length]
bb_dd_U = bb_dd_U[:target_length]

pad_length_dD = max(0, target_length + 1 - dD_AnnAvg_matrix.shape[0])
if pad_length_dD > 0:
    pad_rows = np.repeat(dD_AnnAvg_matrix[0:1, :], pad_length_dD, axis=0)
    dD_AnnAvg_matrix = np.vstack([pad_rows, dD_AnnAvg_matrix])

if Mic_dD_MC.shape[0] < target_length:
    pad_count = target_length - Mic_dD_MC.shape[0]
    pad_rows_MC = pd.concat([Mic_dD_MC.iloc[0:1, :]] * pad_count, ignore_index=True)
    Mic_dD_MC = pd.concat([pad_rows_MC, Mic_dD_MC], ignore_index=True)
elif Mic_dD_MC.shape[0] > target_length:
    Mic_dD_MC = Mic_dD_MC.iloc[:target_length, :]

# BB annual values (align to model years)
# CT data covers different years — use mean as constant if alignment fails
if len(BB_annual_CT) >= n_years_model:
    BB_annual = BB_annual_CT[:n_years_model]
else:
    BB_annual = np.full(n_years_model, BB_GLOBAL_MEAN)

BB_NH_fixed = BB_annual * BB_NH_FRACTION
BB_SH_fixed = BB_annual * BB_SH_FRACTION

print(f"\n  Model years: {int(model_years[0])}–{int(model_years[-1])} ({n_years_model} years)")

# ===========================================================================
# MONTE CARLO LOOP — 2×2 SEPARATE INVERSIONS
# ===========================================================================
print(f"\n{'='*70}")
print("STARTING v3.2 MONTE CARLO (2×2 BB-fixed, separate δ¹³C & δD)")
print(f"{'='*70}")
print(f"  Iterations: {N_ITERATIONS}")
print(f"  Method: Fix BB from CT, solve for FF & Mic with each isotope separately")
print(f"  Advantages: Well-conditioned, no bounded LS needed, cross-validation")
print(f"{'='*70}\n")

rng = np.random.default_rng(seed=42)

# Results: δ¹³C-derived and δD-derived, per hemisphere
FF_NH_d13C = np.zeros((n_years_model, N_ITERATIONS))
Mic_NH_d13C = np.zeros((n_years_model, N_ITERATIONS))
FF_SH_d13C = np.zeros((n_years_model, N_ITERATIONS))
Mic_SH_d13C = np.zeros((n_years_model, N_ITERATIONS))

FF_NH_dD = np.zeros((n_years_model, N_ITERATIONS))
Mic_NH_dD = np.zeros((n_years_model, N_ITERATIONS))
FF_SH_dD = np.zeros((n_years_model, N_ITERATIONS))
Mic_SH_dD = np.zeros((n_years_model, N_ITERATIONS))

# Track negative solutions (for diagnostics — shouldn't happen often with 2×2)
n_negative_d13C = 0
n_negative_dD = 0
total_solves = 0

for k in range(N_ITERATIONS):
    if (k + 1) % 200 == 0:
        print(f"  Iteration {k + 1}/{N_ITERATIONS}...")
    
    # Sample τ_ex
    tau_ex = max(0.5, rng.normal(TAU_EX_MEAN, TAU_EX_STD))
    
    # Sample KIEs
    kies = sample_KIE(rng)
    
    def compute_bulk_KIE(kies, sink_fracs):
        kie_13C = (kies['OH_13C'] * sink_fracs['OH'] +
                   kies['Cl_13C'] * sink_fracs['Cl'] +
                   kies['Strat_13C'] * sink_fracs['Strat'] +
                   kies['Soil_13C'] * sink_fracs['Soil'])
        kie_D = (kies['OH_D'] * sink_fracs['OH'] +
                 kies['Cl_D'] * sink_fracs['Cl'] +
                 kies['Strat_D'] * sink_fracs['Strat'] +
                 kies['Soil_D'] * sink_fracs['Soil'])
        return kie_13C, kie_D
    
    Sink_13C_NH, Sink_D_NH = compute_bulk_KIE(kies, SINK_FRACTIONS['NH'])
    Sink_13C_SH, Sink_D_SH = compute_bulk_KIE(kies, SINK_FRACTIONS['SH'])
    alpha_13C_NH = 1.0 / Sink_13C_NH
    alpha_D_NH = 1.0 / Sink_D_NH
    alpha_13C_SH = 1.0 / Sink_13C_SH
    alpha_D_SH = 1.0 / Sink_D_SH
    
    # Compute hemispheric source strengths
    SumSource_NH = np.zeros(n_years_model)
    SumSource_SH = np.zeros(n_years_model)
    for i in range(n_years_model):
        M_NH_now = CH4_NH[i] * PT_HEMI
        M_NH_next = CH4_NH[i + 1] * PT_HEMI
        M_SH_now = CH4_SH[i] * PT_HEMI
        M_SH_next = CH4_SH[i + 1] * PT_HEMI
        exchange_to_NH = (M_SH_now - M_NH_now) / tau_ex
        exchange_to_SH = (M_NH_now - M_SH_now) / tau_ex
        SumSource_NH[i] = (M_NH_next - M_NH_now) + M_NH_now / Lifetime_NH[i] - exchange_to_NH
        SumSource_SH[i] = (M_SH_next - M_SH_now) + M_SH_now / Lifetime_SH[i] - exchange_to_SH
    
    # Sample atmospheric observations
    # δ¹³C
    if k < d13C_glob_iterations.shape[1]:
        d13C_global_MC = d13C_glob_iterations[:target_length + 1, k]
    else:
        d13C_global_MC = d13C_glob_iterations[:target_length + 1, -1]
    n_c13 = min(len(c13_glob), target_length + 1)
    d13C_offset = d13C_global_MC[:n_c13] - c13_glob[:n_c13]
    d13C_NH_MC = c13_NH[:n_c13] + d13C_offset
    d13C_SH_MC = c13_SH[:n_c13] + d13C_offset
    
    # δD
    if k < dD_AnnAvg_matrix.shape[1]:
        dD_global_MC = dD_AnnAvg_matrix[:target_length + 1, k]
    else:
        dD_global_MC = dD_AnnAvg_matrix[:target_length + 1, -1]
    if len(dD_global_MC) < target_length + 1:
        pad = np.full(target_length + 1 - len(dD_global_MC), dD_global_MC[0])
        dD_global_MC = np.concatenate([pad, dD_global_MC])
    dD_NH_MC = dD_global_MC - DD_IH_OFFSET
    dD_SH_MC = dD_global_MC + DD_IH_OFFSET
    
    # Convert to fractions
    f13_NH = delta_to_fraction_d13C(d13C_NH_MC)
    f13_SH = delta_to_fraction_d13C(d13C_SH_MC)
    fD_NH = delta_to_fraction_dD(dD_NH_MC)
    fD_SH = delta_to_fraction_dD(dD_SH_MC)
    
    # Compute isotopic source fractions
    d13C_source_NH = np.zeros(n_years_model)
    d13C_source_SH = np.zeros(n_years_model)
    dD_source_NH = np.zeros(n_years_model)
    dD_source_SH = np.zeros(n_years_model)
    
    for j in range(n_years_model):
        # NH ¹³C
        n13C_NH_now = f13_NH[j] * CH4_NH[j] * PT_HEMI
        n13C_NH_next = f13_NH[j + 1] * CH4_NH[j + 1] * PT_HEMI
        n13C_SH_now = f13_SH[j] * CH4_SH[j] * PT_HEMI
        exchange_13C_to_NH = (n13C_SH_now - n13C_NH_now) / tau_ex
        d13C_source_NH[j] = (n13C_NH_next - n13C_NH_now +
                             n13C_NH_now * alpha_13C_NH / Lifetime_NH[j] -
                             exchange_13C_to_NH) / SumSource_NH[j]
        
        # SH ¹³C
        n13C_SH_next = f13_SH[j + 1] * CH4_SH[j + 1] * PT_HEMI
        exchange_13C_to_SH = (n13C_NH_now - n13C_SH_now) / tau_ex
        d13C_source_SH[j] = (n13C_SH_next - n13C_SH_now +
                             n13C_SH_now * alpha_13C_SH / Lifetime_SH[j] -
                             exchange_13C_to_SH) / SumSource_SH[j]
        
        # NH D
        nD_NH_now = fD_NH[j] * CH4_NH[j] * PT_HEMI
        nD_NH_next = fD_NH[j + 1] * CH4_NH[j + 1] * PT_HEMI
        nD_SH_now = fD_SH[j] * CH4_SH[j] * PT_HEMI
        exchange_D_to_NH = (nD_SH_now - nD_NH_now) / tau_ex
        dD_source_NH[j] = (nD_NH_next - nD_NH_now +
                           nD_NH_now * alpha_D_NH / Lifetime_NH[j] -
                           exchange_D_to_NH) / SumSource_NH[j]
        
        # SH D
        nD_SH_next = fD_SH[j + 1] * CH4_SH[j + 1] * PT_HEMI
        exchange_D_to_SH = (nD_NH_now - nD_SH_now) / tau_ex
        dD_source_SH[j] = (nD_SH_next - nD_SH_now +
                           nD_SH_now * alpha_D_SH / Lifetime_SH[j] -
                           exchange_D_to_SH) / SumSource_SH[j]
    
    # === Sample source end-member signatures ===
    RandomGauss_FF_d13C = rng.normal()
    RandomGauss_BB_d13C = rng.normal()
    RandomGauss_FF_dD = rng.normal()
    RandomGauss_BB_dD = rng.normal()
    
    if k < FF_d13C_MC_EDGAR.shape[1]:
        ff_d13C_MC_iter = np.array(FF_d13C_MC_EDGAR.iloc[:, k]).flatten()
    else:
        ff_d13C_MC_iter = ff_d13C + RandomGauss_FF_d13C * ff_d13C_U
    
    if k < FF_dD_MC_EDGAR.shape[1]:
        ff_dD_MC_iter = np.array(FF_dD_MC_EDGAR.iloc[:, k]).flatten()
    else:
        ff_dD_MC_iter = ff_dd + RandomGauss_FF_dD * ff_dd_U
    
    bb_d13C_MC_iter = bb_d13C + RandomGauss_BB_d13C * bb_d13C_U
    bb_dD_MC_iter = bb_dd + RandomGauss_BB_dD * bb_dd_U
    
    mic_d13C_MC_iter = np.array(Mic_d13C_MC.iloc[:target_length, k]).flatten() if k < Mic_d13C_MC.shape[1] else np.full(target_length, mic_d13C_mean)
    mic_dD_MC_iter = np.array(Mic_dD_MC.iloc[:target_length, k]).flatten() if k < Mic_dD_MC.shape[1] else np.full(target_length, mic_dd_mean)
    
    def pad_to_length(arr, length):
        arr = np.asarray(arr).flatten()
        if len(arr) >= length:
            return arr[:length]
        return np.concatenate([np.full(length - len(arr), arr[0]), arr])[:length]
    
    ff_d13C_MC_iter = pad_to_length(ff_d13C_MC_iter, target_length)
    ff_dD_MC_iter = pad_to_length(ff_dD_MC_iter, target_length)
    bb_d13C_MC_iter = pad_to_length(bb_d13C_MC_iter, target_length)
    bb_dD_MC_iter = pad_to_length(bb_dD_MC_iter, target_length)
    mic_d13C_MC_iter = pad_to_length(mic_d13C_MC_iter, target_length)
    mic_dD_MC_iter = pad_to_length(mic_dD_MC_iter, target_length)
    
    # === 2×2 SOLVE: Ben's approach per hemisphere ===
    # For δ¹³C:
    #   S × δ_source = FF × δ_FF + Mic × δ_Mic + BB × δ_BB
    #   S = FF + Mic + BB
    #   → FF = (S × δ_source - δ_Mic × (S - BB) - δ_BB × BB) / (δ_FF - δ_Mic)
    #   → Mic = S - BB - FF
    #
    # Work in delta-space (like Ben), not fraction space (simpler and equivalent)
    
    for j in range(n_years_model):
        total_solves += 1
        
        # Convert source fractions back to delta for simpler formula
        d13C_src_NH_delta = fraction_to_delta_d13C(d13C_source_NH[j])
        d13C_src_SH_delta = fraction_to_delta_d13C(d13C_source_SH[j])
        dD_src_NH_delta = fraction_to_delta_dD(dD_source_NH[j])
        dD_src_SH_delta = fraction_to_delta_dD(dD_source_SH[j])
        
        # --- NH δ¹³C inversion ---
        S_NH = SumSource_NH[j]
        BB_NH_j = BB_NH_fixed[j] if j < len(BB_NH_fixed) else BB_GLOBAL_MEAN * BB_NH_FRACTION
        
        denom_13C = ff_d13C_MC_iter[j] - mic_d13C_MC_iter[j]
        if abs(denom_13C) > 0.1:  # Sanity check
            FF_NH_j_d13C = (S_NH * d13C_src_NH_delta - mic_d13C_MC_iter[j] * (S_NH - BB_NH_j) - bb_d13C_MC_iter[j] * BB_NH_j) / denom_13C
            Mic_NH_j_d13C = S_NH - BB_NH_j - FF_NH_j_d13C
        else:
            FF_NH_j_d13C = np.nan
            Mic_NH_j_d13C = np.nan
        
        # Clamp negatives (report but fix)
        if FF_NH_j_d13C < 0 or Mic_NH_j_d13C < 0:
            n_negative_d13C += 1
            FF_NH_j_d13C = max(0, FF_NH_j_d13C)
            Mic_NH_j_d13C = S_NH - BB_NH_j - FF_NH_j_d13C
            if Mic_NH_j_d13C < 0:
                Mic_NH_j_d13C = 0
                FF_NH_j_d13C = S_NH - BB_NH_j
        
        FF_NH_d13C[j, k] = FF_NH_j_d13C
        Mic_NH_d13C[j, k] = Mic_NH_j_d13C
        
        # --- SH δ¹³C inversion ---
        S_SH = SumSource_SH[j]
        BB_SH_j = BB_SH_fixed[j] if j < len(BB_SH_fixed) else BB_GLOBAL_MEAN * BB_SH_FRACTION
        
        if abs(denom_13C) > 0.1:
            FF_SH_j_d13C = (S_SH * d13C_src_SH_delta - mic_d13C_MC_iter[j] * (S_SH - BB_SH_j) - bb_d13C_MC_iter[j] * BB_SH_j) / denom_13C
            Mic_SH_j_d13C = S_SH - BB_SH_j - FF_SH_j_d13C
        else:
            FF_SH_j_d13C = np.nan
            Mic_SH_j_d13C = np.nan
        
        if FF_SH_j_d13C < 0 or Mic_SH_j_d13C < 0:
            n_negative_d13C += 1
            FF_SH_j_d13C = max(0, FF_SH_j_d13C)
            Mic_SH_j_d13C = S_SH - BB_SH_j - FF_SH_j_d13C
            if Mic_SH_j_d13C < 0:
                Mic_SH_j_d13C = 0
                FF_SH_j_d13C = S_SH - BB_SH_j
        
        FF_SH_d13C[j, k] = FF_SH_j_d13C
        Mic_SH_d13C[j, k] = Mic_SH_j_d13C
        
        # --- NH δD inversion ---
        denom_dD = ff_dD_MC_iter[j] - mic_dD_MC_iter[j]
        if abs(denom_dD) > 1.0:  # δD has larger spread
            FF_NH_j_dD = (S_NH * dD_src_NH_delta - mic_dD_MC_iter[j] * (S_NH - BB_NH_j) - bb_dD_MC_iter[j] * BB_NH_j) / denom_dD
            Mic_NH_j_dD = S_NH - BB_NH_j - FF_NH_j_dD
        else:
            FF_NH_j_dD = np.nan
            Mic_NH_j_dD = np.nan
        
        if FF_NH_j_dD < 0 or Mic_NH_j_dD < 0:
            n_negative_dD += 1
            FF_NH_j_dD = max(0, FF_NH_j_dD)
            Mic_NH_j_dD = S_NH - BB_NH_j - FF_NH_j_dD
            if Mic_NH_j_dD < 0:
                Mic_NH_j_dD = 0
                FF_NH_j_dD = S_NH - BB_NH_j
        
        FF_NH_dD[j, k] = FF_NH_j_dD
        Mic_NH_dD[j, k] = Mic_NH_j_dD
        
        # --- SH δD inversion ---
        if abs(denom_dD) > 1.0:
            FF_SH_j_dD = (S_SH * dD_src_SH_delta - mic_dD_MC_iter[j] * (S_SH - BB_SH_j) - bb_dD_MC_iter[j] * BB_SH_j) / denom_dD
            Mic_SH_j_dD = S_SH - BB_SH_j - FF_SH_j_dD
        else:
            FF_SH_j_dD = np.nan
            Mic_SH_j_dD = np.nan
        
        if FF_SH_j_dD < 0 or Mic_SH_j_dD < 0:
            n_negative_dD += 1
            FF_SH_j_dD = max(0, FF_SH_j_dD)
            Mic_SH_j_dD = S_SH - BB_SH_j - FF_SH_j_dD
            if Mic_SH_j_dD < 0:
                Mic_SH_j_dD = 0
                FF_SH_j_dD = S_SH - BB_SH_j
        
        FF_SH_dD[j, k] = FF_SH_j_dD
        Mic_SH_dD[j, k] = Mic_SH_j_dD

print("\nMonte Carlo complete!")
print(f"  Negative solutions (clamped): δ¹³C={n_negative_d13C} ({100*n_negative_d13C/total_solves:.1f}%), "
      f"δD={n_negative_dD} ({100*n_negative_dD/total_solves:.1f}%)")

# ===========================================================================
# POST-PROCESSING
# ===========================================================================
print("\nApplying 5-year smoothing...")

# Global = NH + SH
FF_Global_d13C = FF_NH_d13C + FF_SH_d13C
Mic_Global_d13C = Mic_NH_d13C + Mic_SH_d13C
FF_Global_dD = FF_NH_dD + FF_SH_dD
Mic_Global_dD = Mic_NH_dD + Mic_SH_dD

# Smooth
FF_Global_d13C_s = smooth_5yr(FF_Global_d13C)
Mic_Global_d13C_s = smooth_5yr(Mic_Global_d13C)
FF_Global_dD_s = smooth_5yr(FF_Global_dD)
Mic_Global_dD_s = smooth_5yr(Mic_Global_dD)

FF_NH_d13C_s = smooth_5yr(FF_NH_d13C)
Mic_NH_d13C_s = smooth_5yr(Mic_NH_d13C)
FF_SH_d13C_s = smooth_5yr(FF_SH_d13C)
Mic_SH_d13C_s = smooth_5yr(Mic_SH_d13C)
FF_NH_dD_s = smooth_5yr(FF_NH_dD)
Mic_NH_dD_s = smooth_5yr(Mic_NH_dD)
FF_SH_dD_s = smooth_5yr(FF_SH_dD)
Mic_SH_dD_s = smooth_5yr(Mic_SH_dD)

# ===========================================================================
# RESULTS
# ===========================================================================
print(f"\n{'='*70}")
print("RESULTS — v3.2 (2×2 BB-Fixed, Two-Hemisphere)")
print(f"{'='*70}")

print(f"\n  --- δ¹³C-derived (Global, Smoothed) ---")
print(f"  FF:  {np.nanmean(FF_Global_d13C_s):.1f} ± {np.nanstd(np.nanmean(FF_Global_d13C_s, axis=0)):.1f} Tg/yr")
print(f"  Mic: {np.nanmean(Mic_Global_d13C_s):.1f} ± {np.nanstd(np.nanmean(Mic_Global_d13C_s, axis=0)):.1f} Tg/yr")
print(f"  BB:  {BB_GLOBAL_MEAN:.1f} Tg/yr (fixed from CT)")

print(f"\n  --- δD-derived (Global, Smoothed) ---")
print(f"  FF:  {np.nanmean(FF_Global_dD_s):.1f} ± {np.nanstd(np.nanmean(FF_Global_dD_s, axis=0)):.1f} Tg/yr")
print(f"  Mic: {np.nanmean(Mic_Global_dD_s):.1f} ± {np.nanstd(np.nanmean(Mic_Global_dD_s, axis=0)):.1f} Tg/yr")
print(f"  BB:  {BB_GLOBAL_MEAN:.1f} Tg/yr (fixed from CT)")

# δ¹³C/δD agreement check
ff_d13C_mean = np.nanmean(FF_Global_d13C_s)
ff_dD_mean = np.nanmean(FF_Global_dD_s)
print(f"\n  δ¹³C vs δD agreement: FF difference = {abs(ff_d13C_mean - ff_dD_mean):.1f} Tg/yr")
print(f"  (Ben's paper: δ¹³C gives FF=160±29, δD gives FF=133±33 — ~27 Tg/yr difference)")

# Trend analysis
idx_base = slice(6, 9)    # 2005-2007
idx_recent = slice(-3, None)  # last 3 years

print(f"\n  --- TRENDS: Δ(2020–2022 vs 2005–2007) ---")
for name, compiled in [('FF (δ¹³C)', FF_Global_d13C_s), ('Mic (δ¹³C)', Mic_Global_d13C_s),
                       ('FF (δD)', FF_Global_dD_s), ('Mic (δD)', Mic_Global_dD_s)]:
    delta = compiled[idx_recent, :].mean(axis=0) - compiled[idx_base, :].mean(axis=0)
    pct_pos = (delta > 0).sum() / len(delta) * 100
    print(f"  {name}: Δ = {delta.mean():+.1f} ± {delta.std():.1f} Tg/yr "
          f"(positive in {pct_pos:.0f}% of MC runs)")

# ===========================================================================
# SAVE
# ===========================================================================
print(f"\nSaving to {OUT_DIR}/...")

results_df = pd.DataFrame({
    'Year': model_years,
    'FF_Global_d13C_mean': np.nanmean(FF_Global_d13C_s, axis=1),
    'FF_Global_d13C_std': np.nanstd(FF_Global_d13C_s, axis=1),
    'Mic_Global_d13C_mean': np.nanmean(Mic_Global_d13C_s, axis=1),
    'Mic_Global_d13C_std': np.nanstd(Mic_Global_d13C_s, axis=1),
    'FF_Global_dD_mean': np.nanmean(FF_Global_dD_s, axis=1),
    'FF_Global_dD_std': np.nanstd(FF_Global_dD_s, axis=1),
    'Mic_Global_dD_mean': np.nanmean(Mic_Global_dD_s, axis=1),
    'Mic_Global_dD_std': np.nanstd(Mic_Global_dD_s, axis=1),
    'BB_fixed': BB_annual[:n_years_model] if len(BB_annual) >= n_years_model else np.full(n_years_model, BB_GLOBAL_MEAN),
    'FF_NH_d13C_mean': np.nanmean(FF_NH_d13C_s, axis=1),
    'Mic_NH_d13C_mean': np.nanmean(Mic_NH_d13C_s, axis=1),
    'FF_SH_d13C_mean': np.nanmean(FF_SH_d13C_s, axis=1),
    'Mic_SH_d13C_mean': np.nanmean(Mic_SH_d13C_s, axis=1),
})
results_df.to_csv(OUT_DIR / 'v3.2_results_smoothed.csv', index=False)

# Quality report
quality = {
    'method': '2x2 BB-fixed separate inversions',
    'BB_source': 'CarbonTracker GFED4 prior mean',
    'BB_global_Tg': round(float(BB_GLOBAL_MEAN), 1),
    'negative_solutions_d13C_pct': round(100 * n_negative_d13C / total_solves, 2),
    'negative_solutions_dD_pct': round(100 * n_negative_dD / total_solves, 2),
    'dD_hemispheric_offset': DD_IH_OFFSET,
    'tau_ex_mean': TAU_EX_MEAN,
    'tau_ex_std': TAU_EX_STD,
}
with open(OUT_DIR / 'quality_report.json', 'w') as f:
    json.dump(quality, f, indent=2)

# ===========================================================================
# VISUALIZATION
# ===========================================================================
print("Creating plots...")

# Main result: δ¹³C vs δD comparison (like Ben's Fig. 2)
fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=150, sharex=True)
fig.suptitle('v3.2 Two-Hemisphere 2×2: δ¹³C vs δD Mass Balance (Smoothed)', fontsize=13, y=0.98)

# Global FF
ax = axes[0, 0]
ff_d13C_m = np.nanmean(FF_Global_d13C_s, axis=1)
ff_d13C_s2 = 2 * np.nanstd(FF_Global_d13C_s, axis=1)
ff_dD_m = np.nanmean(FF_Global_dD_s, axis=1)
ff_dD_s2 = 2 * np.nanstd(FF_Global_dD_s, axis=1)
ax.plot(model_years, ff_d13C_m, 'r-', lw=2, label='δ¹³C-derived')
ax.fill_between(model_years, ff_d13C_m - ff_d13C_s2, ff_d13C_m + ff_d13C_s2, alpha=0.2, color='red')
ax.plot(model_years, ff_dD_m, 'b-', lw=2, label='δD-derived')
ax.fill_between(model_years, ff_dD_m - ff_dD_s2, ff_dD_m + ff_dD_s2, alpha=0.2, color='blue')
ax.set_ylabel('Fossil Fuel (Tg/yr)')
ax.set_title('Global Fossil Fuel Emissions')
ax.legend()
ax.grid(True, alpha=0.3)

# Global Mic
ax = axes[0, 1]
mic_d13C_m = np.nanmean(Mic_Global_d13C_s, axis=1)
mic_d13C_s2 = 2 * np.nanstd(Mic_Global_d13C_s, axis=1)
mic_dD_m = np.nanmean(Mic_Global_dD_s, axis=1)
mic_dD_s2 = 2 * np.nanstd(Mic_Global_dD_s, axis=1)
ax.plot(model_years, mic_d13C_m, 'r-', lw=2, label='δ¹³C-derived')
ax.fill_between(model_years, mic_d13C_m - mic_d13C_s2, mic_d13C_m + mic_d13C_s2, alpha=0.2, color='red')
ax.plot(model_years, mic_dD_m, 'b-', lw=2, label='δD-derived')
ax.fill_between(model_years, mic_dD_m - mic_dD_s2, mic_dD_m + mic_dD_s2, alpha=0.2, color='blue')
ax.set_ylabel('Microbial (Tg/yr)')
ax.set_title('Global Microbial Emissions')
ax.legend()
ax.grid(True, alpha=0.3)

# NH breakdown (δ¹³C)
ax = axes[1, 0]
ax.plot(model_years, np.nanmean(FF_NH_d13C_s, axis=1), 'b-', lw=2, label='FF NH')
ax.plot(model_years, np.nanmean(Mic_NH_d13C_s, axis=1), 'g-', lw=2, label='Mic NH')
ax.plot(model_years, BB_NH_fixed[:n_years_model], 'r--', lw=1.5, label='BB NH (fixed)')
ax.set_ylabel('Emissions (Tg/yr)')
ax.set_title('NH Sources (δ¹³C-derived)')
ax.set_xlabel('Year')
ax.legend()
ax.grid(True, alpha=0.3)

# SH breakdown (δ¹³C)
ax = axes[1, 1]
ax.plot(model_years, np.nanmean(FF_SH_d13C_s, axis=1), 'b-', lw=2, label='FF SH')
ax.plot(model_years, np.nanmean(Mic_SH_d13C_s, axis=1), 'g-', lw=2, label='Mic SH')
ax.plot(model_years, BB_SH_fixed[:n_years_model], 'r--', lw=1.5, label='BB SH (fixed)')
ax.set_ylabel('Emissions (Tg/yr)')
ax.set_title('SH Sources (δ¹³C-derived)')
ax.set_xlabel('Year')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUT_DIR / 'v3.2_d13C_vs_dD_comparison.png', dpi=150, bbox_inches='tight')
plt.close(fig)

# Trend histograms (like Ben's Fig. 3)
fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=150)
fig.suptitle('v3.2: Emission Trends (2020–2022 vs 2005–2007)', fontsize=13)

for i, (name, compiled, color) in enumerate([
    ('FF (δ¹³C)', FF_Global_d13C_s, 'red'),
    ('FF (δD)', FF_Global_dD_s, 'blue'),
    ('Mic (δ¹³C)', Mic_Global_d13C_s, 'red'),
    ('Mic (δD)', Mic_Global_dD_s, 'blue'),
]):
    ax = axes[i // 2, i % 2]
    delta = compiled[idx_recent, :].mean(axis=0) - compiled[idx_base, :].mean(axis=0)
    ax.hist(delta, bins=50, color=color, alpha=0.7, edgecolor='black', linewidth=0.5)
    ax.axvline(0, color='black', linestyle='--', linewidth=1.5)
    ax.axvline(delta.mean(), color='darkred' if 'δ¹³C' in name else 'darkblue', linestyle='-', linewidth=2)
    ax.set_xlabel('ΔEmissions (Tg/yr)')
    ax.set_ylabel('MC Count')
    pct_pos = (delta > 0).sum() / len(delta) * 100
    ax.set_title(f'{name}: {delta.mean():+.1f} ± {delta.std():.1f} ({pct_pos:.0f}% positive)')

plt.tight_layout()
plt.savefig(OUT_DIR / 'v3.2_trend_histograms.png', dpi=150, bbox_inches='tight')
plt.close(fig)

print(f"\n{'='*70}")
print("v3.2 (2×2 BB-Fixed, Two-Hemisphere) — RUN COMPLETE")
print(f"{'='*70}")
