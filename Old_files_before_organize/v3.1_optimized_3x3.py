#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Two-Hemisphere Dual-Isotope (δ¹³C + δD) Monte Carlo Box Model
================================================================
Version: 3.1 — Optimized 3×3 (Ben-informed improvements)

CHANGES FROM v3.0:
  1. δD hemispheric offset corrected: ±6‰ (paper shows NH is ~12‰ lower than SH)
  2. All 4 sinks in KIE: OH + Cl + Stratosphere + Soil (all sampled or fixed)
  3. 5-year smoothing post-processing (Ben's approach for robust trend extraction)
  4. Mean/Trend uncertainty separation in MC (Ben's key innovation)
  5. Updated source signature central values from Riddell-Young (2025, PNAS)
  6. τ_ex uncertainty: sampled as Normal(1.0, 0.1) per MC iteration
  7. Net sink KIE computed per hemisphere per iteration (was already done in v3.0)

SCIENTIFIC REFERENCE:
  Riddell-Young et al. (2025, PNAS): "Microbial driver of 2006–2023 CH₄ growth"
  - NH δD–CH₄ is ~12‰ lower than SH δD–CH₄
  - Net sink KIE_13C = 1.0082, KIE_D = 1.281
  - Source signatures: Table 1 (with trend uncertainties)

Author: Upgraded from v3.0 with Ben-model insights
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
from scipy.optimize import lsq_linear

# ---------------------------------------------------------------------------
# I/O paths
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent
ORIGINAL_MODEL_DIR = BASE_DIR.parent / "TwoIsotopeBoxModel"
REL_DIR = ORIGINAL_MODEL_DIR / "rel"
DATA_DIR = REL_DIR / "data"
SRC_DIR = REL_DIR / "output"
OUT_DIR = BASE_DIR / "Output_v3.1_3x3"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
MB_DEBUG = "--debug" in sys.argv
N_ITERATIONS = 1000
TAU_EX_MEAN = 1.0     # Interhemispheric exchange time (years)
TAU_EX_STD = 0.1      # Uncertainty in τ_ex (NEW in v3.1)

# NH/SH emission fractions (same as v3.0)
EMISSION_FRACTIONS = {
    'FF': {'NH': 0.85, 'SH': 0.15},
    'Mic': {'NH': 0.65, 'SH': 0.35},
    'BB': {'NH': 0.55, 'SH': 0.45},
}

# ---------------------------------------------------------------------------
# KIE Sampling Configuration — now includes Strat and Soil as sampable
# Updated from Riddell-Young (2025) SI Table S3
# ---------------------------------------------------------------------------
KIE_CONFIG = {
    'OH_13C':    {'dist': 'uniform', 'low': 1.0039, 'high': 1.0054},
    'OH_D':      {'dist': 'uniform', 'low': 1.294,  'high': 1.327},
    'Cl_13C':    {'dist': 'normal',  'mean': 1.066,  'std': 0.002},
    'Cl_D':      {'dist': 'normal',  'mean': 1.52,   'std': 0.02},
    'Strat_13C': {'dist': 'normal',  'mean': 1.003,  'std': 0.001},  # Dyonisius 2020
    'Strat_D':   {'dist': 'normal',  'mean': 1.179,  'std': 0.01},   # Beck 2018
    'Soil_13C':  {'dist': 'normal',  'mean': 1.0201, 'std': 0.003},  # Snover & Quay avg
    'Soil_D':    {'dist': 'normal',  'mean': 1.083,  'std': 0.01},   # Snover & Quay
}

def sample_KIE(rng):
    """Draw KIE values from configured distributions."""
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
# Sink fractions per hemisphere (unchanged from v3.0)
# ---------------------------------------------------------------------------
SINK_FRACTIONS = {
    'NH': {'OH': 0.825, 'Cl': 0.040, 'Strat': 0.070, 'Soil': 0.065},
    'SH': {'OH': 0.850, 'Cl': 0.028, 'Strat': 0.070, 'Soil': 0.052},
}

# ---------------------------------------------------------------------------
# Time-varying lifetime
# ---------------------------------------------------------------------------
def compute_lifetime_array(years):
    """τ(t) = 9.0 - 0.017*(t - 2010)"""
    return 9.0 - 0.017 * (np.asarray(years, dtype=float) - 2010)

LIFETIME_RATIO = {'NH': 0.95, 'SH': 1.05}

# ---------------------------------------------------------------------------
# δD hemispheric offset — CRITICAL FIX from v3.0
# Riddell-Young (2025): "NH δD–CH₄ is ~12‰ lower than SH δD–CH₄"
# Seasonal range: 8–15‰. Annual mean offset ≈ 12‰.
# So: δD_NH ≈ global - 6‰, δD_SH ≈ global + 6‰
# ---------------------------------------------------------------------------
DD_IH_OFFSET = 6.0  # ‰, NH is 6‰ more negative than global mean (was 1.5 in v3.0!)

# ---------------------------------------------------------------------------
# Isotope utility functions
# ---------------------------------------------------------------------------
C13Std = 0.011113
DStd = 0.00015576
PT = 2.815
PT_HEMI = PT / 2.0

def delta_to_R_d13C(delta_permil):
    return (delta_permil / 1000.0 + 1.0) * C13Std

def delta_to_R_dD(delta_permil):
    return (delta_permil / 1000.0 + 1.0) * DStd

def R_to_fraction(R):
    return R / (1.0 + R)

def fraction_to_R(f):
    return f / (1.0 - f)

def delta_to_fraction_d13C(delta_permil):
    return R_to_fraction(delta_to_R_d13C(delta_permil))

def delta_to_fraction_dD(delta_permil):
    return R_to_fraction(delta_to_R_dD(delta_permil))

def fraction_to_delta_d13C(f):
    return ((fraction_to_R(f) - C13Std) / C13Std) * 1000

def fraction_to_delta_dD(f):
    return ((fraction_to_R(f) - DStd) / DStd) * 1000

# ---------------------------------------------------------------------------
# Solution Quality Monitor
# ---------------------------------------------------------------------------
class SolutionQualityMonitor:
    def __init__(self, n_years, n_iterations, label=""):
        self.label = label
        self.n_years = n_years
        self.n_iterations = n_iterations
        self.condition_numbers = np.zeros((n_years, n_iterations))
        self.is_nonphysical = np.zeros((n_years, n_iterations), dtype=bool)
        self.is_nan = np.zeros((n_years, n_iterations), dtype=bool)
        
    def record(self, year_idx, iter_idx, A, x):
        cond = np.linalg.cond(A)
        self.condition_numbers[year_idx, iter_idx] = cond
        if np.any(~np.isfinite(x)):
            self.is_nan[year_idx, iter_idx] = True
            self.is_nonphysical[year_idx, iter_idx] = True
        elif np.any(x < 0):
            self.is_nonphysical[year_idx, iter_idx] = True
    
    def summary(self):
        total = self.n_years * self.n_iterations
        n_nonphys = np.sum(self.is_nonphysical)
        mean_cond = np.mean(self.condition_numbers)
        max_cond = np.max(self.condition_numbers)
        pct_nonphys = 100.0 * n_nonphys / total
        report = {
            'hemisphere': self.label,
            'total_solves': total,
            'nonphysical_count': int(n_nonphys),
            'nonphysical_pct': round(pct_nonphys, 2),
            'nan_count': int(np.sum(self.is_nan)),
            'mean_condition_number': round(float(mean_cond), 1),
            'max_condition_number': round(float(max_cond), 1),
        }
        print(f"  [{self.label}] {pct_nonphys:.1f}% non-physical, "
              f"mean cond = {mean_cond:.1f} (max = {max_cond:.1f})")
        return report

# ---------------------------------------------------------------------------
# 5-Year Smoothing (from Ben's approach)
# ---------------------------------------------------------------------------
def smooth_5yr(arr_2d):
    """Apply 5-year moving average to each column. Input: [years × iterations]."""
    n_years, n_cols = arr_2d.shape
    if n_years < 5:
        return arr_2d.copy()
    
    result = np.zeros_like(arr_2d)
    # First 2 points: shorter window
    result[0, :] = np.mean(arr_2d[0:3, :], axis=0)
    result[1, :] = np.mean(arr_2d[0:4, :], axis=0)
    # Middle: full 5-year window
    for i in range(2, n_years - 2):
        result[i, :] = np.mean(arr_2d[i-2:i+3, :], axis=0)
    # Last 2 points
    result[-2, :] = np.mean(arr_2d[-4:, :], axis=0)
    result[-1, :] = np.mean(arr_2d[-3:, :], axis=0)
    return result


# ===========================================================================
# DATA LOADING
# ===========================================================================
print("="*70)
print("TWO-HEMISPHERE DUAL-ISOTOPE BOX MODEL v3.1 (Optimized 3×3)")
print("="*70)
print("\nLoading data...")

# --- CH₄ concentrations ---
CH4data_raw = pd.read_excel(DATA_DIR / "GML_CH4_AnnualMean.xlsx").to_numpy()
CH4_global = CH4data_raw[15:39, 1].astype(float)  # 1999-2022
CH4_years = CH4data_raw[15:39, 0].astype(float)

# IH gradient (same as v3.0)
IH_GRADIENT = np.linspace(80, 100, len(CH4_global))
CH4_NH = CH4_global + IH_GRADIENT / 2.0
CH4_SH = CH4_global - IH_GRADIENT / 2.0

print(f"  CH₄ (1999): Global={CH4_global[0]:.1f}, NH={CH4_NH[0]:.1f}, SH={CH4_SH[0]:.1f} ppb")
print(f"  CH₄ (2022): Global={CH4_global[-1]:.1f}, NH={CH4_NH[-1]:.1f}, SH={CH4_SH[-1]:.1f} ppb")

# --- δ¹³C-CH₄: NH and SH ---
C13data = pd.read_excel(DATA_DIR / "ch4c13_nh_sh_mean.xlsx", header=None).to_numpy()
c13_dates = C13data[:, 0]
c13_global = C13data[:, 1]
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

c13_ann_years, c13_ann_global = annual_average(c13_dates, c13_global)
_, c13_ann_NH = annual_average(c13_dates, c13_NH_raw)
_, c13_ann_SH = annual_average(c13_dates, c13_SH_raw)

c13_start_idx = np.where(c13_ann_years == 1999)[0][0]
c13_end_idx = np.where(c13_ann_years == 2022)[0][0] + 1
c13_NH = c13_ann_NH[c13_start_idx:c13_end_idx]
c13_SH = c13_ann_SH[c13_start_idx:c13_end_idx]
c13_glob = c13_ann_global[c13_start_idx:c13_end_idx]

print(f"  δ¹³C (1999): NH={c13_NH[0]:.3f}‰, SH={c13_SH[0]:.3f}‰")
print(f"  δD hemispheric offset: ±{DD_IH_OFFSET}‰ (corrected from ±1.5‰ in v3.0)")

# --- δD-CH₄: Global MC iterations ---
glob_ann_dD_path = DATA_DIR / "GlobMean_dD_iterations_UmezawaCal_noBUDS.xlsx"
glob_ann_dD_df = pd.read_excel(glob_ann_dD_path)
glob_ann_dD_num = glob_ann_dD_df.apply(pd.to_numeric, errors="coerce")
glob_ann_dD_years = glob_ann_dD_num.iloc[:, 0].to_numpy(dtype=np.float64)
glob_ann_dD = glob_ann_dD_num.iloc[:, 1].to_numpy(dtype=np.float64) - 0.5
dD_AnnAvg_matrix = glob_ann_dD_num.iloc[:, 5:].to_numpy(dtype=np.float64)

# --- d13C DEI iterations ---
d13C_glob_iterations_data = np.loadtxt(str(DATA_DIR / "d13C_dei_compiled.txt"))
d13C_glob_iterations = d13C_glob_iterations_data[1:, 1:]

# --- Source signatures ---
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

# CarbonTracker BB
data_CT = pd.read_excel(DATA_DIR / "CarbonTracker_CH4.xlsx")
bbCT = data_CT.iloc[:, 9].values
BB_global = np.mean(bbCT)

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
mic_dd_U = 8.2  # Updated from Ben (was 7.0)
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

# Pad arrays to target_length
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

print(f"\n  Model years: {int(model_years[0])}–{int(model_years[-1])} ({n_years_model} years)")
print(f"  τ_ex = {TAU_EX_MEAN} ± {TAU_EX_STD} yr (sampled per iteration)")
print(f"  Lifetime NH: {Lifetime_NH[0]:.3f}–{Lifetime_NH[-1]:.3f} yr")
print(f"  Lifetime SH: {Lifetime_SH[0]:.3f}–{Lifetime_SH[-1]:.3f} yr")

# ===========================================================================
# MONTE CARLO LOOP
# ===========================================================================
print(f"\n{'='*70}")
print("STARTING v3.1 MONTE CARLO ANALYSIS (3×3 with Ben optimizations)")
print(f"{'='*70}")
print(f"  Iterations: {N_ITERATIONS}")
print(f"  Key changes: δD offset=±{DD_IH_OFFSET}‰, τ_ex sampled, all 4 sinks sampled")
print(f"{'='*70}\n")

rng = np.random.default_rng(seed=42)

# Result arrays
BB_NH_compiled = np.zeros((n_years_model, N_ITERATIONS))
FF_NH_compiled = np.zeros((n_years_model, N_ITERATIONS))
Mic_NH_compiled = np.zeros((n_years_model, N_ITERATIONS))
BB_SH_compiled = np.zeros((n_years_model, N_ITERATIONS))
FF_SH_compiled = np.zeros((n_years_model, N_ITERATIONS))
Mic_SH_compiled = np.zeros((n_years_model, N_ITERATIONS))

qm_NH = SolutionQualityMonitor(n_years_model, N_ITERATIONS, "NH")
qm_SH = SolutionQualityMonitor(n_years_model, N_ITERATIONS, "SH")

# Track τ_ex samples
tau_ex_samples = np.zeros(N_ITERATIONS)

for k in range(N_ITERATIONS):
    if (k + 1) % 200 == 0:
        print(f"  Iteration {k + 1}/{N_ITERATIONS}...")
    
    # === Sample τ_ex for this iteration (NEW in v3.1) ===
    tau_ex = max(0.5, rng.normal(TAU_EX_MEAN, TAU_EX_STD))  # Clamp ≥ 0.5
    tau_ex_samples[k] = tau_ex
    
    # === Sample KIE values ===
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
    
    # === Compute hemispheric source strengths (depends on τ_ex) ===
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
    
    # === Sample atmospheric observations ===
    # δ¹³C
    if k < d13C_glob_iterations.shape[1]:
        d13C_global_MC = d13C_glob_iterations[:target_length + 1, k]
    else:
        d13C_global_MC = d13C_glob_iterations[:target_length + 1, -1]
    
    n_c13 = min(len(c13_glob), target_length + 1)
    d13C_offset = d13C_global_MC[:n_c13] - c13_glob[:n_c13]
    d13C_NH_MC = c13_NH[:n_c13] + d13C_offset
    d13C_SH_MC = c13_SH[:n_c13] + d13C_offset
    
    # δD with CORRECTED offset
    if k < dD_AnnAvg_matrix.shape[1]:
        dD_global_MC = dD_AnnAvg_matrix[:target_length + 1, k]
    else:
        dD_global_MC = dD_AnnAvg_matrix[:target_length + 1, -1]
    if len(dD_global_MC) < target_length + 1:
        pad = np.full(target_length + 1 - len(dD_global_MC), dD_global_MC[0])
        dD_global_MC = np.concatenate([pad, dD_global_MC])
    
    dD_NH_MC = dD_global_MC - DD_IH_OFFSET   # NH is MORE depleted in D
    dD_SH_MC = dD_global_MC + DD_IH_OFFSET   # SH is LESS depleted in D
    
    # Convert to fractions
    f13_NH = delta_to_fraction_d13C(d13C_NH_MC)
    f13_SH = delta_to_fraction_d13C(d13C_SH_MC)
    fD_NH = delta_to_fraction_dD(dD_NH_MC)
    fD_SH = delta_to_fraction_dD(dD_SH_MC)
    
    # === Compute isotopic source signatures per hemisphere ===
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
    
    # Convert to fractions
    f13_bb = delta_to_fraction_d13C(bb_d13C_MC_iter)
    f13_ff = delta_to_fraction_d13C(ff_d13C_MC_iter)
    f13_mic = delta_to_fraction_d13C(mic_d13C_MC_iter)
    fD_bb = delta_to_fraction_dD(bb_dD_MC_iter)
    fD_ff = delta_to_fraction_dD(ff_dD_MC_iter)
    fD_mic = delta_to_fraction_dD(mic_dD_MC_iter)
    
    # === SOLVE 3×3 per hemisphere per year ===
    for j in range(n_years_model):
        # NH
        A_NH = np.array([
            [1.0, 1.0, 1.0],
            [f13_bb[j], f13_ff[j], f13_mic[j]],
            [fD_bb[j], fD_ff[j], fD_mic[j]]
        ])
        B_NH = np.array([
            SumSource_NH[j],
            SumSource_NH[j] * d13C_source_NH[j],
            SumSource_NH[j] * dD_source_NH[j]
        ])
        W_NH = np.diag([100.0, 1.0, 0.5])
        A_w = W_NH @ A_NH
        B_w = W_NH @ B_NH
        ub = SumSource_NH[j] * 1.5
        try:
            result = lsq_linear(A_w, B_w, bounds=(0, ub))
            x_NH = result.x
            qm_NH.record(j, k, A_NH, x_NH)
        except Exception:
            x_NH = np.array([np.nan, np.nan, np.nan])
            qm_NH.is_nan[j, k] = True
            qm_NH.is_nonphysical[j, k] = True
        
        BB_NH_compiled[j, k] = x_NH[0]
        FF_NH_compiled[j, k] = x_NH[1]
        Mic_NH_compiled[j, k] = x_NH[2]
        
        # SH
        A_SH = np.array([
            [1.0, 1.0, 1.0],
            [f13_bb[j], f13_ff[j], f13_mic[j]],
            [fD_bb[j], fD_ff[j], fD_mic[j]]
        ])
        B_SH = np.array([
            SumSource_SH[j],
            SumSource_SH[j] * d13C_source_SH[j],
            SumSource_SH[j] * dD_source_SH[j]
        ])
        W_SH = np.diag([200.0, 1.0, 0.5])
        A_w = W_SH @ A_SH
        B_w = W_SH @ B_SH
        ub = SumSource_SH[j] * 1.5
        try:
            result = lsq_linear(A_w, B_w, bounds=(0, ub))
            x_SH = result.x
            qm_SH.record(j, k, A_SH, x_SH)
        except Exception:
            x_SH = np.array([np.nan, np.nan, np.nan])
            qm_SH.is_nan[j, k] = True
            qm_SH.is_nonphysical[j, k] = True
        
        BB_SH_compiled[j, k] = x_SH[0]
        FF_SH_compiled[j, k] = x_SH[1]
        Mic_SH_compiled[j, k] = x_SH[2]

print("\nMonte Carlo complete!")

# ===========================================================================
# POST-PROCESSING: 5-Year Smoothing
# ===========================================================================
print("\nApplying 5-year smoothing...")

# Global compiled
BB_Global_compiled = BB_NH_compiled + BB_SH_compiled
FF_Global_compiled = FF_NH_compiled + FF_SH_compiled
Mic_Global_compiled = Mic_NH_compiled + Mic_SH_compiled

# Smoothed versions
BB_Global_smooth = smooth_5yr(BB_Global_compiled)
FF_Global_smooth = smooth_5yr(FF_Global_compiled)
Mic_Global_smooth = smooth_5yr(Mic_Global_compiled)
BB_NH_smooth = smooth_5yr(BB_NH_compiled)
FF_NH_smooth = smooth_5yr(FF_NH_compiled)
Mic_NH_smooth = smooth_5yr(Mic_NH_compiled)
BB_SH_smooth = smooth_5yr(BB_SH_compiled)
FF_SH_smooth = smooth_5yr(FF_SH_compiled)
Mic_SH_smooth = smooth_5yr(Mic_SH_compiled)

# ===========================================================================
# QUALITY & STATISTICS
# ===========================================================================
print(f"\n{'='*70}")
print("SOLUTION QUALITY REPORT")
print(f"{'='*70}")
report_NH = qm_NH.summary()
report_SH = qm_SH.summary()

print(f"\n  τ_ex sampling: mean={tau_ex_samples.mean():.3f}, std={tau_ex_samples.std():.3f}")

# Compute stats on SMOOTHED data
def compute_stats(compiled, smooth):
    """Compute mean/std for both raw and smoothed."""
    return {
        'raw_mean': np.nanmean(compiled, axis=1),
        'raw_std': np.nanstd(compiled, axis=1),
        'smooth_mean': np.nanmean(smooth, axis=1),
        'smooth_std': np.nanstd(smooth, axis=1),
    }

stats = {
    'BB_NH': compute_stats(BB_NH_compiled, BB_NH_smooth),
    'FF_NH': compute_stats(FF_NH_compiled, FF_NH_smooth),
    'Mic_NH': compute_stats(Mic_NH_compiled, Mic_NH_smooth),
    'BB_SH': compute_stats(BB_SH_compiled, BB_SH_smooth),
    'FF_SH': compute_stats(FF_SH_compiled, FF_SH_smooth),
    'Mic_SH': compute_stats(Mic_SH_compiled, Mic_SH_smooth),
    'BB_Global': compute_stats(BB_Global_compiled, BB_Global_smooth),
    'FF_Global': compute_stats(FF_Global_compiled, FF_Global_smooth),
    'Mic_Global': compute_stats(Mic_Global_compiled, Mic_Global_smooth),
}

print(f"\n  --- GLOBAL (Smoothed) ---")
print(f"  BB:  {stats['BB_Global']['smooth_mean'].mean():.1f} ± {stats['BB_Global']['smooth_std'].mean():.1f} Tg/yr")
print(f"  FF:  {stats['FF_Global']['smooth_mean'].mean():.1f} ± {stats['FF_Global']['smooth_std'].mean():.1f} Tg/yr")
print(f"  Mic: {stats['Mic_Global']['smooth_mean'].mean():.1f} ± {stats['Mic_Global']['smooth_std'].mean():.1f} Tg/yr")

print(f"\n  --- NORTHERN HEMISPHERE (Smoothed) ---")
print(f"  BB:  {stats['BB_NH']['smooth_mean'].mean():.1f} ± {stats['BB_NH']['smooth_std'].mean():.1f} Tg/yr")
print(f"  FF:  {stats['FF_NH']['smooth_mean'].mean():.1f} ± {stats['FF_NH']['smooth_std'].mean():.1f} Tg/yr")
print(f"  Mic: {stats['Mic_NH']['smooth_mean'].mean():.1f} ± {stats['Mic_NH']['smooth_std'].mean():.1f} Tg/yr")

print(f"\n  --- SOUTHERN HEMISPHERE (Smoothed) ---")
print(f"  BB:  {stats['BB_SH']['smooth_mean'].mean():.1f} ± {stats['BB_SH']['smooth_std'].mean():.1f} Tg/yr")
print(f"  FF:  {stats['FF_SH']['smooth_mean'].mean():.1f} ± {stats['FF_SH']['smooth_std'].mean():.1f} Tg/yr")
print(f"  Mic: {stats['Mic_SH']['smooth_mean'].mean():.1f} ± {stats['Mic_SH']['smooth_std'].mean():.1f} Tg/yr")

# ===========================================================================
# TREND ANALYSIS (Ben's approach: 2005-07 avg vs 2020-22 avg)
# ===========================================================================
print(f"\n{'='*70}")
print("TREND ANALYSIS: Δ(2020–2022 vs 2005–2007)")
print(f"{'='*70}")

# Indices relative to 1999: 2005=idx6, 2007=idx8, 2020=idx21, 2022=idx23(or end)
idx_base = slice(6, 9)   # 2005-2007
idx_recent = slice(-3, None)  # last 3 years

for name, compiled in [('FF_Global', FF_Global_smooth), ('Mic_Global', Mic_Global_smooth),
                       ('BB_Global', BB_Global_smooth)]:
    delta = compiled[idx_recent, :].mean(axis=0) - compiled[idx_base, :].mean(axis=0)
    print(f"  {name}: Δ = {delta.mean():+.1f} ± {delta.std():.1f} Tg/yr "
          f"(positive in {(delta > 0).sum()/len(delta)*100:.0f}% of MC runs)")

# ===========================================================================
# SAVE RESULTS
# ===========================================================================
print(f"\nSaving results to {OUT_DIR}/...")

results_df = pd.DataFrame({
    'Year': model_years,
    'BB_Global_smooth_mean': stats['BB_Global']['smooth_mean'],
    'BB_Global_smooth_std': stats['BB_Global']['smooth_std'],
    'FF_Global_smooth_mean': stats['FF_Global']['smooth_mean'],
    'FF_Global_smooth_std': stats['FF_Global']['smooth_std'],
    'Mic_Global_smooth_mean': stats['Mic_Global']['smooth_mean'],
    'Mic_Global_smooth_std': stats['Mic_Global']['smooth_std'],
    'BB_NH_smooth_mean': stats['BB_NH']['smooth_mean'],
    'FF_NH_smooth_mean': stats['FF_NH']['smooth_mean'],
    'Mic_NH_smooth_mean': stats['Mic_NH']['smooth_mean'],
    'BB_SH_smooth_mean': stats['BB_SH']['smooth_mean'],
    'FF_SH_smooth_mean': stats['FF_SH']['smooth_mean'],
    'Mic_SH_smooth_mean': stats['Mic_SH']['smooth_mean'],
})
results_df.to_csv(OUT_DIR / 'v3.1_results_smoothed.csv', index=False)

with open(OUT_DIR / 'quality_report.json', 'w') as f:
    json.dump({'NH': report_NH, 'SH': report_SH,
               'tau_ex_mean': float(tau_ex_samples.mean()),
               'tau_ex_std': float(tau_ex_samples.std())}, f, indent=2)

# ===========================================================================
# VISUALIZATION
# ===========================================================================
print("Creating plots...")

fig, axes = plt.subplots(3, 2, figsize=(14, 12), dpi=150, sharex=True)
fig.suptitle('v3.1 Two-Hemisphere 3×3 Model (Ben-optimized): Smoothed Sources', fontsize=13, y=0.98)

for row, (name, color) in enumerate([('BB', 'red'), ('FF', 'blue'), ('Mic', 'green')]):
    # NH
    ax = axes[row, 0]
    sm = stats[f'{name}_NH']['smooth_mean']
    ss = stats[f'{name}_NH']['smooth_std']
    ax.plot(model_years, sm, '-', lw=2.5, color=color)
    ax.fill_between(model_years, sm - ss, sm + ss, alpha=0.3, color=color)
    ax.set_ylabel(f'{name} (Tg/yr)')
    if row == 0: ax.set_title('Northern Hemisphere (5yr smoothed)')
    ax.grid(True, alpha=0.3)
    
    # SH
    ax = axes[row, 1]
    sm = stats[f'{name}_SH']['smooth_mean']
    ss = stats[f'{name}_SH']['smooth_std']
    ax.plot(model_years, sm, '-', lw=2.5, color=color)
    ax.fill_between(model_years, sm - ss, sm + ss, alpha=0.3, color=color)
    if row == 0: ax.set_title('Southern Hemisphere (5yr smoothed)')
    ax.grid(True, alpha=0.3)
    
    if row == 2:
        axes[row, 0].set_xlabel('Year')
        axes[row, 1].set_xlabel('Year')

plt.tight_layout()
plt.savefig(OUT_DIR / 'v3.1_hemispheric_sources_smoothed.png', dpi=150, bbox_inches='tight')
plt.close(fig)

# Trend histogram
fig, axes = plt.subplots(1, 3, figsize=(15, 5), dpi=150)
fig.suptitle('v3.1: Emission Change (2020–2022 vs 2005–2007)', fontsize=13)

for i, (name, color, compiled) in enumerate([
    ('Fossil Fuel', 'blue', FF_Global_smooth),
    ('Microbial', 'green', Mic_Global_smooth),
    ('Biomass Burning', 'red', BB_Global_smooth)
]):
    delta = compiled[idx_recent, :].mean(axis=0) - compiled[idx_base, :].mean(axis=0)
    axes[i].hist(delta, bins=50, color=color, alpha=0.7, edgecolor='black', linewidth=0.5)
    axes[i].axvline(0, color='black', linestyle='--', linewidth=1.5)
    axes[i].axvline(delta.mean(), color=color, linestyle='-', linewidth=2)
    axes[i].set_xlabel('ΔEmissions (Tg/yr)')
    axes[i].set_ylabel('MC Count')
    axes[i].set_title(f'{name}: {delta.mean():+.1f} ± {delta.std():.1f}')

plt.tight_layout()
plt.savefig(OUT_DIR / 'v3.1_trend_histograms.png', dpi=150, bbox_inches='tight')
plt.close(fig)

print(f"\n{'='*70}")
print("v3.1 (3×3 Optimized) — RUN COMPLETE")
print(f"{'='*70}")
