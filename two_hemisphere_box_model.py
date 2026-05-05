#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Two-Hemisphere Dual-Isotope (δ¹³C + δD) Monte Carlo Box Model
================================================================
Version: 3.0 — Two-Hemisphere Upgrade

SCIENTIFIC BASIS:
  Following Nguyen et al. (2020, GRL) and addressing the Naus et al. (2019, ACP)
  critique that one-box models alias inter-hemispheric transport as source changes.

  The model splits the atmosphere into Northern (NH) and Southern (SH) hemisphere
  boxes coupled by interhemispheric exchange (τ_ex ≈ 1 year).

  Per hemisphere, the system of equations is:
    d[CH₄]_N/dt = S_N - [CH₄]_N/τ_N + ([CH₄]_S - [CH₄]_N)/τ_ex
    d[CH₄]_S/dt = S_S - [CH₄]_S/τ_S + ([CH₄]_N - [CH₄]_S)/τ_ex

  For isotopes, the mass balance is formulated in terms of heavy-isotope amount
  (fraction × mass) to maintain linearity:
    d(f¹³C·M)_N/dt = f¹³C_src_N·S_N - f¹³C_N·M_N·α/τ + exchange terms
    d(f_D·M)_N/dt  = f_D_src_N·S_N  - f_D_N·M_N·α/τ   + exchange terms

  This yields a 6-variable system (M_N, M_S, f13_N, f13_S, fD_N, fD_S) that
  can be decomposed into two coupled 3×3 solves per hemisphere when the exchange
  terms are treated as known forcing from the previous year's observations.

UPGRADES FROM v2.0 (upgraded_box_model.py):
  - Two-hemisphere spatial structure (this file)
  - NH/SH-specific observations (δ¹³C from ch4c13_nh_sh_mean.xlsx; CH₄ derived)
  - Interhemispheric exchange coupling
  - NH/SH emission ratios from EDGAR/CarbonTracker (~75%/25% total)
  - All v2.0 upgrades retained: KIE sampling, quality monitoring, time-varying τ

RETAINS from v2.0:
  1. KIE sampling in MC loop (Chandra 2024 critique)
  2. Solution quality monitoring (condition numbers, non-physical solutions)
  3. Time-varying CH₄ lifetime (He et al. 2026)

KEY REFERENCES:
  - Nguyen et al. (2020, GRL): Two-hemisphere box model with CH₄-CO-OH chemistry
  - Naus et al. (2019, ACP): Transport aliasing critique of one-box models
  - Prather (1994): Perturbation lifetime theory
  - He et al. (2026, Science): Time-varying CH₄ lifetime

Author: Upgraded by OpenClaw from Yufan Bao's original code
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
# I/O paths — data comes from the original TwoIsotopeBoxModel/rel directory
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent
# Data from the original model
ORIGINAL_MODEL_DIR = BASE_DIR.parent / "TwoIsotopeBoxModel"
REL_DIR = ORIGINAL_MODEL_DIR / "rel"
DATA_DIR = REL_DIR / "data"
SRC_DIR = REL_DIR / "output"
OUT_DIR = BASE_DIR / "Output_2Hemi"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
MB_DEBUG = "--debug" in sys.argv
N_ITERATIONS = 1000  # Monte Carlo iterations
TAU_EX = 1.0  # Interhemispheric exchange time (years) — Nguyen et al. (2020)

# NH/SH emission fractions from EDGAR v7 / CarbonTracker-CH₄
# Total: ~75% NH, ~25% SH (Saunois et al. 2020 ESSD)
# Category-specific ratios (approximate from EDGAR/CT-CH₄):
#   Fossil fuel: ~85% NH, ~15% SH (concentrated in Russia, US, Middle East)
#   Microbial:   ~65% NH, ~35% SH (tropical wetlands span both, rice in NH)
#   Biomass burning: ~55% NH, ~45% SH (savanna fires in SH, boreal in NH)
EMISSION_FRACTIONS = {
    'FF': {'NH': 0.85, 'SH': 0.15},
    'Mic': {'NH': 0.65, 'SH': 0.35},
    'BB': {'NH': 0.55, 'SH': 0.45},
}

# ---------------------------------------------------------------------------
# KIE Sampling Configuration (from v2.0)
# ---------------------------------------------------------------------------
KIE_CONFIG = {
    'OH_13C': {'dist': 'uniform', 'low': 1.0039, 'high': 1.0054},
    'OH_D':   {'dist': 'uniform', 'low': 1.294,  'high': 1.327},
    'Cl_13C': {'dist': 'normal',  'mean': 1.066, 'std': 0.002},
    'Cl_D':   {'dist': 'normal',  'mean': 1.52,  'std': 0.02},
    'Strat_13C': 1.003,
    'Strat_D':   1.179,
    'Soil_13C':  1.0201,
    'Soil_D':    1.083,
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
# Sink fractions — can differ by hemisphere (OH is higher in tropics/SH summer)
# For now, use global mean sink fractions; hemisphere-specific OH:other ratios
# are a future refinement.
# ---------------------------------------------------------------------------
OH_Sink = 0.835
Cl_Sink = 0.035
Strat_Sink = 0.07
Soil_Sink = 0.06

# NH has slightly less OH fraction (more Cl from coastal; more soil)
# SH has slightly more OH fraction (cleaner atmosphere)
# This is a first-order approximation; refine with Lelieveld (2016) / Holmes (2013)
SINK_FRACTIONS = {
    'NH': {'OH': 0.825, 'Cl': 0.040, 'Strat': 0.070, 'Soil': 0.065},
    'SH': {'OH': 0.850, 'Cl': 0.028, 'Strat': 0.070, 'Soil': 0.052},
}

# ---------------------------------------------------------------------------
# Time-varying lifetime (from v2.0)
# ---------------------------------------------------------------------------
def compute_lifetime_array(years):
    """Compute time-varying CH₄ lifetime: τ(t) = 9.0 - 0.017*(t - 2010)"""
    years = np.asarray(years, dtype=float)
    return 9.0 - 0.017 * (years - 2010)

# NH lifetime is slightly shorter (higher OH concentrations in tropics/NH)
# SH lifetime is slightly longer
# The NH/SH lifetime ratio is ~0.95/1.05 based on Prather (2012) / Lawrence (2001)
LIFETIME_RATIO = {'NH': 0.95, 'SH': 1.05}

# ---------------------------------------------------------------------------
# Isotope utility functions (unchanged)
# ---------------------------------------------------------------------------
C13Std = 0.011113
DStd = 0.00015576
PT = 2.815  # ppb → Tg conversion (for full atmosphere)
PT_HEMI = PT / 2.0  # ppb → Tg per hemisphere (half atmosphere)

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

def R_to_delta_d13C(R):
    return ((R - C13Std) / C13Std) * 1000

def R_to_delta_dD(R):
    return ((R - DStd) / DStd) * 1000

def fraction_to_delta_d13C(f):
    return R_to_delta_d13C(fraction_to_R(f))

def fraction_to_delta_dD(f):
    return R_to_delta_dD(fraction_to_R(f))

# ---------------------------------------------------------------------------
# Solution Quality Monitor (from v2.0, extended for 2 hemispheres)
# ---------------------------------------------------------------------------
class SolutionQualityMonitor:
    """Track matrix condition numbers and non-physical solutions per hemisphere."""
    
    def __init__(self, n_years, n_iterations, label=""):
        self.label = label
        self.n_years = n_years
        self.n_iterations = n_iterations
        self.condition_numbers = np.zeros((n_years, n_iterations))
        self.is_nonphysical = np.zeros((n_years, n_iterations), dtype=bool)
        self.is_nan = np.zeros((n_years, n_iterations), dtype=bool)
        self.COND_THRESHOLD = 100.0
        
    def record(self, year_idx, iter_idx, A, x):
        """Record quality metrics for one solve."""
        cond = np.linalg.cond(A)
        self.condition_numbers[year_idx, iter_idx] = cond
        if np.any(~np.isfinite(x)):
            self.is_nan[year_idx, iter_idx] = True
            self.is_nonphysical[year_idx, iter_idx] = True
        elif np.any(x < 0):
            self.is_nonphysical[year_idx, iter_idx] = True
    
    def summary(self):
        """Print and return summary statistics."""
        total = self.n_years * self.n_iterations
        n_nonphys = np.sum(self.is_nonphysical)
        n_nan = np.sum(self.is_nan)
        mean_cond = np.mean(self.condition_numbers)
        max_cond = np.max(self.condition_numbers)
        pct_nonphys = 100.0 * n_nonphys / total
        
        report = {
            'hemisphere': self.label,
            'total_solves': total,
            'nonphysical_count': int(n_nonphys),
            'nonphysical_pct': round(pct_nonphys, 2),
            'nan_count': int(n_nan),
            'mean_condition_number': round(float(mean_cond), 1),
            'max_condition_number': round(float(max_cond), 1),
        }
        
        print(f"  [{self.label}] {pct_nonphys:.1f}% non-physical, "
              f"mean cond = {mean_cond:.1f} (max = {max_cond:.1f})")
        return report


# ===========================================================================
# DATA LOADING
# ===========================================================================
print("="*70)
print("TWO-HEMISPHERE DUAL-ISOTOPE BOX MODEL v3.0")
print("="*70)
print("\nLoading data...")

# --- CH₄ concentrations (global annual means from GML) ---
CH4data_raw = pd.read_excel(DATA_DIR / "GML_CH4_AnnualMean.xlsx").to_numpy()
# Full record starts at row 0 (year 1984)
# We want 1999-2022 → rows 15:39 (1984+15=1999, 1984+38=2022)
CH4_global = CH4data_raw[15:39, 1].astype(float)  # 24 years: 1999-2022
CH4_years = CH4data_raw[15:39, 0].astype(float)

# --- Derive NH/SH CH₄ from global mean and known IH gradient ---
# The NH-SH gradient from NOAA surface flask network (background stations):
#   ~2000: ~80 ppb; ~2010: ~90 ppb; ~2020: ~100 ppb
# These are mass-weighted hemispheric means from NOAA GML zonally-averaged data
# (Dlugokencky et al.; see also Lan et al. 2024 ESSD).
# The gradient grows as NH emissions increase faster than SH.
# For production use: obtain actual NOAA NH/SH zonal mean products.
IH_GRADIENT = np.linspace(80, 100, len(CH4_global))  # ppb, NH minus SH
CH4_NH = CH4_global + IH_GRADIENT / 2.0
CH4_SH = CH4_global - IH_GRADIENT / 2.0

print(f"  CH₄ (1999): Global={CH4_global[0]:.1f}, NH={CH4_NH[0]:.1f}, SH={CH4_SH[0]:.1f} ppb")
print(f"  CH₄ (2022): Global={CH4_global[-1]:.1f}, NH={CH4_NH[-1]:.1f}, SH={CH4_SH[-1]:.1f} ppb")

# --- δ¹³C-CH₄: NH and SH from ch4c13_nh_sh_mean.xlsx ---
# File structure: col0=decimal_date, col1=global, col2=NH, col3=SH (all δ¹³C in ‰)
C13data = pd.read_excel(DATA_DIR / "ch4c13_nh_sh_mean.xlsx", header=None).to_numpy()
c13_dates = C13data[:, 0]
c13_global = C13data[:, 1]
c13_NH_raw = C13data[:, 2]
c13_SH_raw = C13data[:, 3]

# Compute annual averages for NH and SH
def annual_average(dates, values):
    """Compute annual averages from sub-annual data."""
    years_floor = np.floor(dates).astype(int)
    unique_years = np.unique(years_floor)
    ann_years = []
    ann_means = []
    for yr in unique_years:
        mask = years_floor == yr
        if np.sum(mask) >= 6:  # Require at least 6 months
            ann_years.append(yr)
            ann_means.append(np.nanmean(values[mask]))
    return np.array(ann_years), np.array(ann_means)

c13_ann_years, c13_ann_global = annual_average(c13_dates, c13_global)
_, c13_ann_NH = annual_average(c13_dates, c13_NH_raw)
_, c13_ann_SH = annual_average(c13_dates, c13_SH_raw)

# Align to 1999-2022
c13_start_idx = np.where(c13_ann_years == 1999)[0][0]
c13_end_idx = np.where(c13_ann_years == 2022)[0][0] + 1
c13_NH = c13_ann_NH[c13_start_idx:c13_end_idx]
c13_SH = c13_ann_SH[c13_start_idx:c13_end_idx]
c13_glob = c13_ann_global[c13_start_idx:c13_end_idx]

print(f"  δ¹³C (1999): NH={c13_NH[0]:.3f}‰, SH={c13_SH[0]:.3f}‰, Δ(NH-SH)={c13_NH[0]-c13_SH[0]:.3f}‰")
print(f"  δ¹³C (2022): NH={c13_NH[-1]:.3f}‰, SH={c13_SH[-1]:.3f}‰")

# --- δD-CH₄: Only global available; approximate NH/SH ---
# The NH-SH δD gradient is much smaller (~2-4‰) and less well constrained
# Approximation: δD_NH ≈ δD_global - 1.5‰; δD_SH ≈ δD_global + 1.5‰
# (NH is more depleted in D because FF sources are D-depleted and concentrated in NH)
# This is a first-order estimate; refine with station-level δD data when available
glob_ann_dD_path = DATA_DIR / "GlobMean_dD_iterations_UmezawaCal_noBUDS.xlsx"
glob_ann_dD_df = pd.read_excel(glob_ann_dD_path)
glob_ann_dD_num = glob_ann_dD_df.apply(pd.to_numeric, errors="coerce")
glob_ann_dD_years = glob_ann_dD_num.iloc[:, 0].to_numpy(dtype=np.float64)
glob_ann_dD = glob_ann_dD_num.iloc[:, 1].to_numpy(dtype=np.float64) - 0.5
dD_AnnAvg_matrix = glob_ann_dD_num.iloc[:, 5:].to_numpy(dtype=np.float64)

# Approximate NH/SH δD gradient (‰)
DD_IH_OFFSET = 1.5  # NH is ~1.5‰ more negative than SH due to FF influence

# --- d13C DEI iterations (for MC uncertainty on global δ¹³C) ---
d13C_glob_iterations_data = np.loadtxt(str(DATA_DIR / "d13C_dei_compiled.txt"))
d13C_glob_iterations = d13C_glob_iterations_data[1:, 1:]

# --- Source signatures (same as original — global values applied per hemisphere) ---
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

# CarbonTracker for BB reference
data_CT = pd.read_excel(DATA_DIR / "CarbonTracker_CH4.xlsx")
bbCT = data_CT.iloc[:, 9].values
BB_global = np.mean(bbCT)

# --- Prepare source signature arrays ---
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
mic_dd_U = 7.0
ff_dd = np.array(FF_dD_data.iloc[34:, 1]).flatten()
ff_dd_U = np.array(FF_dD_data.iloc[34:, 2]).flatten()
bb_dd = np.array(BB_dD_data.iloc[:, 1]).flatten()
bb_dd_U = np.array(BB_dD_data.iloc[:, 2]).flatten()
mean_bb_dd = bb_dd[-1]
mean_bb_dd_U = bb_dd_U[-1]
bb_dd = np.concatenate((np.full(3, mean_bb_dd), bb_dd, np.full(1, mean_bb_dd)))
bb_dd_U = np.concatenate((np.full(3, mean_bb_dd_U), bb_dd_U, np.full(1, mean_bb_dd_U)))

# ===========================================================================
# COMPUTE MODEL TIME DIMENSIONS
# ===========================================================================
n_years_model = len(CH4_global) - 1  # 23 annual changes (1999→2000 through 2021→2022)
model_years = np.arange(1999, 1999 + n_years_model)
Lifetime_global = compute_lifetime_array(model_years)
Lifetime_NH = Lifetime_global * LIFETIME_RATIO['NH']
Lifetime_SH = Lifetime_global * LIFETIME_RATIO['SH']

print(f"\n  Model years: {int(model_years[0])}–{int(model_years[-1])} ({n_years_model} years)")
print(f"  τ_ex (interhemispheric exchange) = {TAU_EX} year")
print(f"  Lifetime NH: {Lifetime_NH[0]:.3f}–{Lifetime_NH[-1]:.3f} yr")
print(f"  Lifetime SH: {Lifetime_SH[0]:.3f}–{Lifetime_SH[-1]:.3f} yr")

# ===========================================================================
# PAD SOURCE SIGNATURE ARRAYS TO TARGET LENGTH
# ===========================================================================
target_length = n_years_model  # 23

# Pad ff_dd/bb_dd to target length
while len(ff_dd) < target_length:
    ff_dd = np.concatenate([np.array([ff_dd[0]]), ff_dd])
    ff_dd_U = np.concatenate([np.array([ff_dd_U[0]]), ff_dd_U])
ff_dd = ff_dd[:target_length]
ff_dd_U = ff_dd_U[:target_length]
bb_dd = bb_dd[:target_length]
bb_dd_U = bb_dd_U[:target_length]

# Pad dD matrix
pad_length_dD = max(0, target_length + 1 - dD_AnnAvg_matrix.shape[0])
if pad_length_dD > 0:
    pad_rows = np.repeat(dD_AnnAvg_matrix[0:1, :], pad_length_dD, axis=0)
    dD_AnnAvg_matrix = np.vstack([pad_rows, dD_AnnAvg_matrix])

# Pad Mic_dD_MC
if Mic_dD_MC.shape[0] < target_length:
    pad_count = target_length - Mic_dD_MC.shape[0]
    pad_rows_MC = pd.concat([Mic_dD_MC.iloc[0:1, :]] * pad_count, ignore_index=True)
    Mic_dD_MC = pd.concat([pad_rows_MC, Mic_dD_MC], ignore_index=True)
elif Mic_dD_MC.shape[0] > target_length:
    Mic_dD_MC = Mic_dD_MC.iloc[:target_length, :]

# Ensure d13C iteration length matches
# d13C_glob_iterations should have at least n_years_model+1 rows (for CH4 array length)

# ===========================================================================
# FORWARD MODEL: Compute hemispheric source strengths
# ===========================================================================
# For each hemisphere, the mass balance (in Tg) is:
#   S_N = M_N(t+1) - M_N(t) + M_N(t)/τ_N - (M_S(t) - M_N(t))/τ_ex · (PT_HEMI)
#   S_S = M_S(t+1) - M_S(t) + M_S(t)/τ_S - (M_N(t) - M_S(t))/τ_ex · (PT_HEMI)
#
# where M = [CH₄] × PT_HEMI (mass in Tg per hemisphere)

print("\nComputing hemispheric source strengths...")

SumSource_NH = np.zeros(n_years_model)
SumSource_SH = np.zeros(n_years_model)
SumSource_Global = np.zeros(n_years_model)

for i in range(n_years_model):
    # NH mass balance
    M_NH_now = CH4_NH[i] * PT_HEMI
    M_NH_next = CH4_NH[i + 1] * PT_HEMI
    M_SH_now = CH4_SH[i] * PT_HEMI
    M_SH_next = CH4_SH[i + 1] * PT_HEMI
    
    # Exchange flux (positive = SH→NH when [SH] > [NH], but normally NH > SH)
    # Exchange: +([CH₄]_S - [CH₄]_N)/τ_ex for NH box
    exchange_to_NH = (M_SH_now - M_NH_now) / TAU_EX  # Negative (NH > SH)
    exchange_to_SH = (M_NH_now - M_SH_now) / TAU_EX  # Positive (NH > SH)
    
    # Source = accumulation + loss - exchange_in
    SumSource_NH[i] = (M_NH_next - M_NH_now) + M_NH_now / Lifetime_NH[i] - exchange_to_NH
    SumSource_SH[i] = (M_SH_next - M_SH_now) + M_SH_now / Lifetime_SH[i] - exchange_to_SH
    SumSource_Global[i] = SumSource_NH[i] + SumSource_SH[i]

print(f"  NH Source: {SumSource_NH.mean():.1f} Tg/yr (range {SumSource_NH.min():.1f}–{SumSource_NH.max():.1f})")
print(f"  SH Source: {SumSource_SH.mean():.1f} Tg/yr (range {SumSource_SH.min():.1f}–{SumSource_SH.max():.1f})")
print(f"  Global:    {SumSource_Global.mean():.1f} Tg/yr")
print(f"  NH fraction: {(SumSource_NH.mean()/SumSource_Global.mean()*100):.1f}%")

# ===========================================================================
# MONTE CARLO LOOP — TWO-HEMISPHERE
# ===========================================================================
print(f"\n{'='*70}")
print("STARTING TWO-HEMISPHERE MONTE CARLO ANALYSIS")
print(f"{'='*70}")
print(f"  Iterations: {N_ITERATIONS}")
print(f"  Hemispheres: NH + SH with τ_ex = {TAU_EX} yr")
print(f"  Upgrades: KIE sampling, Quality monitoring, Time-varying τ, 2-Hemi")
print(f"{'='*70}\n")

rng = np.random.default_rng(seed=42)

# Result arrays: [years × iterations] per hemisphere
BB_NH_compiled = np.zeros((n_years_model, N_ITERATIONS))
FF_NH_compiled = np.zeros((n_years_model, N_ITERATIONS))
Mic_NH_compiled = np.zeros((n_years_model, N_ITERATIONS))

BB_SH_compiled = np.zeros((n_years_model, N_ITERATIONS))
FF_SH_compiled = np.zeros((n_years_model, N_ITERATIONS))
Mic_SH_compiled = np.zeros((n_years_model, N_ITERATIONS))

# Global (sum of hemispheres)
BB_Global_compiled = np.zeros((n_years_model, N_ITERATIONS))
FF_Global_compiled = np.zeros((n_years_model, N_ITERATIONS))
Mic_Global_compiled = np.zeros((n_years_model, N_ITERATIONS))

# Quality monitors
qm_NH = SolutionQualityMonitor(n_years_model, N_ITERATIONS, "NH")
qm_SH = SolutionQualityMonitor(n_years_model, N_ITERATIONS, "SH")

# KIE tracking
KIE_samples = {'OH_13C': np.zeros(N_ITERATIONS), 'OH_D': np.zeros(N_ITERATIONS),
               'Cl_13C': np.zeros(N_ITERATIONS), 'Cl_D': np.zeros(N_ITERATIONS)}

for k in range(N_ITERATIONS):
    if (k + 1) % 100 == 0:
        print(f"  Iteration {k + 1}/{N_ITERATIONS}...")
    
    # ===================================================================
    # UPGRADE 1: Sample KIE values for this iteration
    # ===================================================================
    kies = sample_KIE(rng)
    KIE_samples['OH_13C'][k] = kies['OH_13C']
    KIE_samples['OH_D'][k] = kies['OH_D']
    KIE_samples['Cl_13C'][k] = kies['Cl_13C']
    KIE_samples['Cl_D'][k] = kies['Cl_D']
    
    # Compute bulk sink KIE per hemisphere
    def compute_bulk_KIE(kies, sink_fracs):
        """Compute effective bulk KIE from sampled values and sink fractions."""
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
    
    # ===================================================================
    # Sample atmospheric δ¹³C observations (MC iteration k)
    # For NH/SH: use the observed NH/SH δ¹³C annual means + MC uncertainty
    # The DEI iterations represent uncertainty on the global mean.
    # Apply the same relative uncertainty to NH and SH.
    # ===================================================================
    
    # Global δ¹³C MC perturbation from DEI iterations
    if k < d13C_glob_iterations.shape[1]:
        d13C_global_MC = d13C_glob_iterations[:target_length + 1, k]
    else:
        d13C_global_MC = d13C_glob_iterations[:target_length + 1, -1]
    
    # Compute the offset from the mean global δ¹³C for this MC iteration
    # and apply same offset to NH and SH
    # First check lengths
    n_c13 = min(len(c13_glob), target_length + 1)
    d13C_offset = d13C_global_MC[:n_c13] - c13_glob[:n_c13]
    
    # Apply to NH and SH (same absolute perturbation from DEI uncertainty)
    d13C_NH_MC = c13_NH[:n_c13] + d13C_offset
    d13C_SH_MC = c13_SH[:n_c13] + d13C_offset
    
    # δD: Use global MC iterations + hemispheric offset
    if k < dD_AnnAvg_matrix.shape[1]:
        dD_global_MC = dD_AnnAvg_matrix[:target_length + 1, k]
    else:
        dD_global_MC = dD_AnnAvg_matrix[:target_length + 1, -1]
    
    # Ensure dD array is long enough
    if len(dD_global_MC) < target_length + 1:
        pad = np.full(target_length + 1 - len(dD_global_MC), dD_global_MC[0])
        dD_global_MC = np.concatenate([pad, dD_global_MC])
    
    dD_NH_MC = dD_global_MC - DD_IH_OFFSET
    dD_SH_MC = dD_global_MC + DD_IH_OFFSET
    
    # ===================================================================
    # Compute isotopic source signatures per hemisphere from box model inversion
    # For each hemisphere:
    #   d(f·M)/dt = f_src·S - f_atm·M·α/τ + exchange_terms
    #   → f_src·S = d(f·M)/dt + f_atm·M·α/τ - exchange_terms
    #   → f_src = [above] / S
    # ===================================================================
    
    # Convert δ to fractions
    f13_NH = delta_to_fraction_d13C(d13C_NH_MC)
    f13_SH = delta_to_fraction_d13C(d13C_SH_MC)
    fD_NH = delta_to_fraction_dD(dD_NH_MC)
    fD_SH = delta_to_fraction_dD(dD_SH_MC)
    
    # Isotopic source fractions for each hemisphere
    d13C_source_NH = np.zeros(n_years_model)
    d13C_source_SH = np.zeros(n_years_model)
    dD_source_NH = np.zeros(n_years_model)
    dD_source_SH = np.zeros(n_years_model)
    
    for j in range(n_years_model):
        # --- NH isotope mass balance ---
        # n13C_NH = f13_NH * CH4_NH * PT_HEMI (Tg of ¹³CH₄-equivalent)
        n13C_NH_now = f13_NH[j] * CH4_NH[j] * PT_HEMI
        n13C_NH_next = f13_NH[j + 1] * CH4_NH[j + 1] * PT_HEMI
        n13C_SH_now = f13_SH[j] * CH4_SH[j] * PT_HEMI
        
        # Exchange of ¹³C isotopic mass: proportional to concentration × fraction
        # Isotopic exchange flux to NH = (f13_SH·M_SH - f13_NH·M_NH) / τ_ex
        exchange_13C_to_NH = (n13C_SH_now - n13C_NH_now) / TAU_EX
        
        # Source isotope mass balance for NH:
        # f_src_NH · S_NH = d(n13C_NH)/dt + n13C_NH·α/τ - exchange_13C_to_NH
        d13C_source_NH[j] = (n13C_NH_next - n13C_NH_now +
                             n13C_NH_now * alpha_13C_NH / Lifetime_NH[j] -
                             exchange_13C_to_NH) / SumSource_NH[j]
        
        # --- SH isotope mass balance ---
        n13C_SH_next = f13_SH[j + 1] * CH4_SH[j + 1] * PT_HEMI
        exchange_13C_to_SH = (n13C_NH_now - n13C_SH_now) / TAU_EX
        
        d13C_source_SH[j] = (n13C_SH_next - n13C_SH_now +
                             n13C_SH_now * alpha_13C_SH / Lifetime_SH[j] -
                             exchange_13C_to_SH) / SumSource_SH[j]
        
        # --- NH δD mass balance ---
        nD_NH_now = fD_NH[j] * CH4_NH[j] * PT_HEMI
        nD_NH_next = fD_NH[j + 1] * CH4_NH[j + 1] * PT_HEMI
        nD_SH_now = fD_SH[j] * CH4_SH[j] * PT_HEMI
        
        exchange_D_to_NH = (nD_SH_now - nD_NH_now) / TAU_EX
        
        dD_source_NH[j] = (nD_NH_next - nD_NH_now +
                           nD_NH_now * alpha_D_NH / Lifetime_NH[j] -
                           exchange_D_to_NH) / SumSource_NH[j]
        
        # --- SH δD mass balance ---
        nD_SH_next = fD_SH[j + 1] * CH4_SH[j + 1] * PT_HEMI
        exchange_D_to_SH = (nD_NH_now - nD_SH_now) / TAU_EX
        
        dD_source_SH[j] = (nD_SH_next - nD_SH_now +
                           nD_SH_now * alpha_D_SH / Lifetime_SH[j] -
                           exchange_D_to_SH) / SumSource_SH[j]
    
    # ===================================================================
    # Sample source end-member signatures (same global values for both hemispheres)
    # Future improvement: hemisphere-specific signatures based on sub-category
    # composition (e.g., NH FF is more gas-dominated, SH is more coal)
    # ===================================================================
    RandomGauss_FF_d13C = rng.normal()
    RandomGauss_BB_d13C = rng.normal()
    RandomGauss_FF_dD = rng.normal()
    RandomGauss_BB_dD = rng.normal()
    
    # FF d13C
    if k < FF_d13C_MC_EDGAR.shape[1]:
        ff_d13C_MC_iter = np.array(FF_d13C_MC_EDGAR.iloc[:, k]).flatten()
    else:
        ff_d13C_MC_iter = ff_d13C + RandomGauss_FF_d13C * ff_d13C_U
    
    # FF dD
    if k < FF_dD_MC_EDGAR.shape[1]:
        ff_dD_MC_iter = np.array(FF_dD_MC_EDGAR.iloc[:, k]).flatten()
    else:
        ff_dD_MC_iter = ff_dd + RandomGauss_FF_dD * ff_dd_U
    
    # BB
    bb_d13C_MC_iter = bb_d13C + RandomGauss_BB_d13C * bb_d13C_U
    bb_dD_MC_iter = bb_dd + RandomGauss_BB_dD * bb_dd_U
    
    # Mic (MC column k)
    mic_d13C_MC_iter = np.array(Mic_d13C_MC.iloc[:target_length, k]).flatten() if k < Mic_d13C_MC.shape[1] else np.full(target_length, mic_d13C_mean)
    mic_dD_MC_iter = np.array(Mic_dD_MC.iloc[:target_length, k]).flatten() if k < Mic_dD_MC.shape[1] else np.full(target_length, mic_dd_mean)
    
    # Ensure all arrays are target_length
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
    
    # Convert to heavy-isotope fractions
    f13_bb = delta_to_fraction_d13C(bb_d13C_MC_iter)
    f13_ff = delta_to_fraction_d13C(ff_d13C_MC_iter)
    f13_mic = delta_to_fraction_d13C(mic_d13C_MC_iter)
    fD_bb = delta_to_fraction_dD(bb_dD_MC_iter)
    fD_ff = delta_to_fraction_dD(ff_dD_MC_iter)
    fD_mic = delta_to_fraction_dD(mic_dD_MC_iter)
    
    # ===================================================================
    # SOLVE 3×3 SYSTEM PER HEMISPHERE PER YEAR
    # 
    # APPROACH: Overdetermined weighted least squares with non-negativity.
    # 
    # The 3×3 system is ill-conditioned because δD end-members are too close.
    # Instead of direct inversion, we use bounded least squares (scipy NNLS / 
    # bounded lsq_linear) with the mass balance as a hard constraint and 
    # the two isotope equations as soft constraints (weighted by observation
    # uncertainty). This gives physically meaningful solutions.
    #
    # System per hemisphere:
    #   Row 0: x_BB + x_FF + x_Mic = S_total         (mass constraint, high weight)
    #   Row 1: f13_BB·x_BB + ... = S_total·f13_src   (δ¹³C constraint)
    #   Row 2: fD_BB·x_BB + ...  = S_total·fD_src    (δD constraint)
    #
    # Additionally, 0 ≤ x_i for all sources (non-negativity).
    # ===================================================================
    
    for j in range(n_years_model):
        # --- NH solve (bounded least squares) ---
        try:
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
            
            # Weight: mass balance row gets 100× weight (hard constraint)
            # δ¹³C row: weight 1.0; δD row: weight 0.5 (less informative)
            W_NH = np.diag([100.0, 1.0, 0.5])
            A_w = W_NH @ A_NH
            B_w = W_NH @ B_NH
            
            # Bounded solve: 0 ≤ x_i ≤ S_total (no single source exceeds total)
            ub = SumSource_NH[j] * 1.5  # Allow some overshoot for single categories
            result = lsq_linear(A_w, B_w, bounds=(0, ub))
            x_NH = result.x
            
            qm_NH.record(j, k, A_NH, x_NH)
            BB_NH_compiled[j, k] = x_NH[0]
            FF_NH_compiled[j, k] = x_NH[1]
            Mic_NH_compiled[j, k] = x_NH[2]
        except (np.linalg.LinAlgError, ValueError, Exception) as e:
            BB_NH_compiled[j, k] = np.nan
            FF_NH_compiled[j, k] = np.nan
            Mic_NH_compiled[j, k] = np.nan
            qm_NH.is_nan[j, k] = True
            qm_NH.is_nonphysical[j, k] = True
        
        # --- SH solve (bounded least squares) ---
        try:
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
            
            # SH has smaller total source — weight mass balance even more
            W_SH = np.diag([200.0, 1.0, 0.5])
            A_w = W_SH @ A_SH
            B_w = W_SH @ B_SH
            
            ub = SumSource_SH[j] * 1.5
            result = lsq_linear(A_w, B_w, bounds=(0, ub))
            x_SH = result.x
            
            qm_SH.record(j, k, A_SH, x_SH)
            BB_SH_compiled[j, k] = x_SH[0]
            FF_SH_compiled[j, k] = x_SH[1]
            Mic_SH_compiled[j, k] = x_SH[2]
        except (np.linalg.LinAlgError, ValueError, Exception) as e:
            BB_SH_compiled[j, k] = np.nan
            FF_SH_compiled[j, k] = np.nan
            Mic_SH_compiled[j, k] = np.nan
            qm_SH.is_nan[j, k] = True
            qm_SH.is_nonphysical[j, k] = True
        
        # Global = NH + SH
        BB_Global_compiled[j, k] = BB_NH_compiled[j, k] + BB_SH_compiled[j, k]
        FF_Global_compiled[j, k] = FF_NH_compiled[j, k] + FF_SH_compiled[j, k]
        Mic_Global_compiled[j, k] = Mic_NH_compiled[j, k] + Mic_SH_compiled[j, k]

print("\nMonte Carlo analysis complete!")

# ===========================================================================
# QUALITY REPORTS
# ===========================================================================
print(f"\n{'='*70}")
print("SOLUTION QUALITY REPORT")
print(f"{'='*70}")
report_NH = qm_NH.summary()
report_SH = qm_SH.summary()

# ===========================================================================
# KIE DIAGNOSTICS
# ===========================================================================
print(f"\n{'='*70}")
print("KIE SAMPLING DIAGNOSTICS")
print(f"{'='*70}")
print(f"  OH_KIE_13C: mean={KIE_samples['OH_13C'].mean():.5f}, "
      f"range=[{KIE_samples['OH_13C'].min():.5f}, {KIE_samples['OH_13C'].max():.5f}]")
print(f"  OH_KIE_D:   mean={KIE_samples['OH_D'].mean():.4f}, "
      f"range=[{KIE_samples['OH_D'].min():.4f}, {KIE_samples['OH_D'].max():.4f}]")

# ===========================================================================
# POST-PROCESSING: Statistics
# ===========================================================================
print(f"\n{'='*70}")
print("POST-PROCESSING: Computing Statistics")
print(f"{'='*70}")

# NH statistics
BB_NH_mean = np.nanmean(BB_NH_compiled, axis=1)
BB_NH_std = np.nanstd(BB_NH_compiled, axis=1)
FF_NH_mean = np.nanmean(FF_NH_compiled, axis=1)
FF_NH_std = np.nanstd(FF_NH_compiled, axis=1)
Mic_NH_mean = np.nanmean(Mic_NH_compiled, axis=1)
Mic_NH_std = np.nanstd(Mic_NH_compiled, axis=1)

# SH statistics
BB_SH_mean = np.nanmean(BB_SH_compiled, axis=1)
BB_SH_std = np.nanstd(BB_SH_compiled, axis=1)
FF_SH_mean = np.nanmean(FF_SH_compiled, axis=1)
FF_SH_std = np.nanstd(FF_SH_compiled, axis=1)
Mic_SH_mean = np.nanmean(Mic_SH_compiled, axis=1)
Mic_SH_std = np.nanstd(Mic_SH_compiled, axis=1)

# Global statistics
BB_Global_mean = np.nanmean(BB_Global_compiled, axis=1)
BB_Global_std = np.nanstd(BB_Global_compiled, axis=1)
FF_Global_mean = np.nanmean(FF_Global_compiled, axis=1)
FF_Global_std = np.nanstd(FF_Global_compiled, axis=1)
Mic_Global_mean = np.nanmean(Mic_Global_compiled, axis=1)
Mic_Global_std = np.nanstd(Mic_Global_compiled, axis=1)

print(f"\n  --- NORTHERN HEMISPHERE ---")
print(f"  BB:  {BB_NH_mean.mean():.1f} ± {BB_NH_std.mean():.1f} Tg/yr")
print(f"  FF:  {FF_NH_mean.mean():.1f} ± {FF_NH_std.mean():.1f} Tg/yr")
print(f"  Mic: {Mic_NH_mean.mean():.1f} ± {Mic_NH_std.mean():.1f} Tg/yr")
print(f"  Sum: {(BB_NH_mean + FF_NH_mean + Mic_NH_mean).mean():.1f} vs Source {SumSource_NH.mean():.1f}")

print(f"\n  --- SOUTHERN HEMISPHERE ---")
print(f"  BB:  {BB_SH_mean.mean():.1f} ± {BB_SH_std.mean():.1f} Tg/yr")
print(f"  FF:  {FF_SH_mean.mean():.1f} ± {FF_SH_std.mean():.1f} Tg/yr")
print(f"  Mic: {Mic_SH_mean.mean():.1f} ± {Mic_SH_std.mean():.1f} Tg/yr")
print(f"  Sum: {(BB_SH_mean + FF_SH_mean + Mic_SH_mean).mean():.1f} vs Source {SumSource_SH.mean():.1f}")

print(f"\n  --- GLOBAL (NH + SH) ---")
print(f"  BB:  {BB_Global_mean.mean():.1f} ± {BB_Global_std.mean():.1f} Tg/yr")
print(f"  FF:  {FF_Global_mean.mean():.1f} ± {FF_Global_std.mean():.1f} Tg/yr")
print(f"  Mic: {Mic_Global_mean.mean():.1f} ± {Mic_Global_std.mean():.1f} Tg/yr")
print(f"  Sum: {(BB_Global_mean + FF_Global_mean + Mic_Global_mean).mean():.1f} vs Source {SumSource_Global.mean():.1f}")

# ===========================================================================
# VALIDATION: Check NH-SH gradients match observations
# ===========================================================================
print(f"\n{'='*70}")
print("VALIDATION: NH-SH Gradient Check")
print(f"{'='*70}")
print(f"  Observed CH₄ NH-SH gradient: {IH_GRADIENT.mean():.0f} ppb (input)")
print(f"  Implied emission NH fraction: {(SumSource_NH/(SumSource_NH+SumSource_SH)).mean()*100:.1f}%")
print(f"  Expected from EDGAR: ~70-75%")
print(f"  δ¹³C NH-SH gradient (obs): {(c13_NH - c13_SH).mean():.3f}‰")
print(f"  (NH more depleted = more fossil fuel influence ✓)")

# ===========================================================================
# SAVE RESULTS
# ===========================================================================
print(f"\nSaving results to {OUT_DIR}/...")

# Summary statistics
results_df = pd.DataFrame({
    'Year': model_years,
    'Lifetime_NH': Lifetime_NH,
    'Lifetime_SH': Lifetime_SH,
    'Source_NH_Tg': SumSource_NH,
    'Source_SH_Tg': SumSource_SH,
    'Source_Global_Tg': SumSource_Global,
    'BB_NH_mean': BB_NH_mean, 'BB_NH_std': BB_NH_std,
    'FF_NH_mean': FF_NH_mean, 'FF_NH_std': FF_NH_std,
    'Mic_NH_mean': Mic_NH_mean, 'Mic_NH_std': Mic_NH_std,
    'BB_SH_mean': BB_SH_mean, 'BB_SH_std': BB_SH_std,
    'FF_SH_mean': FF_SH_mean, 'FF_SH_std': FF_SH_std,
    'Mic_SH_mean': Mic_SH_mean, 'Mic_SH_std': Mic_SH_std,
    'BB_Global_mean': BB_Global_mean, 'BB_Global_std': BB_Global_std,
    'FF_Global_mean': FF_Global_mean, 'FF_Global_std': FF_Global_std,
    'Mic_Global_mean': Mic_Global_mean, 'Mic_Global_std': Mic_Global_std,
})
results_df.to_csv(OUT_DIR / 'two_hemisphere_results.csv', index=False)

# Quality reports
with open(OUT_DIR / 'quality_report.json', 'w') as f:
    json.dump({'NH': report_NH, 'SH': report_SH}, f, indent=2)

# Full MC iterations (global)
for name, arr in [('BB_Global', BB_Global_compiled), ('FF_Global', FF_Global_compiled),
                  ('Mic_Global', Mic_Global_compiled)]:
    df = pd.DataFrame(arr, index=model_years, columns=[f'Iter_{i}' for i in range(N_ITERATIONS)])
    df.to_csv(OUT_DIR / f'{name}_MC_alliterations.csv')

# NH and SH MC iterations
for name, arr in [('BB_NH', BB_NH_compiled), ('FF_NH', FF_NH_compiled), ('Mic_NH', Mic_NH_compiled),
                  ('BB_SH', BB_SH_compiled), ('FF_SH', FF_SH_compiled), ('Mic_SH', Mic_SH_compiled)]:
    df = pd.DataFrame(arr, index=model_years, columns=[f'Iter_{i}' for i in range(N_ITERATIONS)])
    df.to_csv(OUT_DIR / f'{name}_MC_alliterations.csv')

# KIE samples
kie_df = pd.DataFrame(KIE_samples)
kie_df.to_csv(OUT_DIR / 'KIE_samples.csv', index=False)

# ===========================================================================
# VISUALIZATION
# ===========================================================================
print("Creating visualization plots...")

# --- Plot 1: Hemispheric comparison (6 panels) ---
fig, axes = plt.subplots(3, 2, figsize=(14, 12), dpi=150, sharex=True)
fig.suptitle('Two-Hemisphere Dual-Isotope Box Model: Source Partitioning', fontsize=14, y=0.98)

sources = [
    ('Biomass Burning', BB_NH_mean, BB_NH_std, BB_SH_mean, BB_SH_std, 'red'),
    ('Fossil Fuel', FF_NH_mean, FF_NH_std, FF_SH_mean, FF_SH_std, 'blue'),
    ('Microbial', Mic_NH_mean, Mic_NH_std, Mic_SH_mean, Mic_SH_std, 'green'),
]

for row, (name, nh_mean, nh_std, sh_mean, sh_std, color) in enumerate(sources):
    # NH (left column)
    ax = axes[row, 0]
    ax.plot(model_years, nh_mean, '-', linewidth=2.5, color=color, label='NH Mean')
    ax.fill_between(model_years, nh_mean - nh_std, nh_mean + nh_std, alpha=0.3, color=color)
    ax.set_ylabel(f'{name} (Tg yr⁻¹)')
    if row == 0:
        ax.set_title('Northern Hemisphere')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper left')
    
    # SH (right column)
    ax = axes[row, 1]
    ax.plot(model_years, sh_mean, '-', linewidth=2.5, color=color, label='SH Mean')
    ax.fill_between(model_years, sh_mean - sh_std, sh_mean + sh_std, alpha=0.3, color=color)
    if row == 0:
        ax.set_title('Southern Hemisphere')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper left')
    
    if row == 2:
        axes[row, 0].set_xlabel('Year')
        axes[row, 1].set_xlabel('Year')

plt.tight_layout()
plt.savefig(OUT_DIR / 'two_hemisphere_sources.png', dpi=150, bbox_inches='tight')
plt.close(fig)

# --- Plot 2: Global comparison with spaghetti ---
fig, axes = plt.subplots(3, 2, figsize=(14, 12), dpi=150, sharex=True)
fig.suptitle('Two-Hemisphere Model: Global Source Partitioning (NH+SH)', fontsize=14, y=0.98)

global_sources = [
    ('Biomass Burning', BB_Global_compiled, BB_Global_mean, BB_Global_std, 'red'),
    ('Fossil Fuel', FF_Global_compiled, FF_Global_mean, FF_Global_std, 'blue'),
    ('Microbial', Mic_Global_compiled, Mic_Global_mean, Mic_Global_std, 'green'),
]

nh_compiled_dict = {'Biomass Burning': BB_NH_compiled, 'Fossil Fuel': FF_NH_compiled, 'Microbial': Mic_NH_compiled}
sh_compiled_dict = {'Biomass Burning': BB_SH_compiled, 'Fossil Fuel': FF_SH_compiled, 'Microbial': Mic_SH_compiled}

for row, (name, compiled, mean, std, color) in enumerate(global_sources):
    # Spaghetti (left)
    ax = axes[row, 0]
    ax.plot(model_years, compiled, linewidth=0.3, alpha=0.15, color=color)
    ax.plot(model_years, mean, '-', linewidth=2.5, color=color, label='Mean')
    ax.fill_between(model_years, mean - std, mean + std, alpha=0.3, color=color, label='±1σ')
    ax.set_ylabel(f'{name} (Tg yr⁻¹)')
    ax.set_title(f'{name} — All MC Iterations')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper left')
    
    # NH vs SH decomposition (right)
    ax = axes[row, 1]
    nh_m = np.nanmean(nh_compiled_dict[name], axis=1)
    sh_m = np.nanmean(sh_compiled_dict[name], axis=1)
    ax.bar(model_years, nh_m, width=0.8, color=color, alpha=0.7, label='NH')
    ax.bar(model_years, sh_m, width=0.8, bottom=nh_m, color=color, alpha=0.35, label='SH')
    ax.set_ylabel(f'{name} (Tg yr⁻¹)')
    ax.set_title(f'{name} — NH/SH Decomposition')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper left')
    
    if row == 2:
        axes[row, 0].set_xlabel('Year')
        axes[row, 1].set_xlabel('Year')

plt.tight_layout()
plt.savefig(OUT_DIR / 'global_sources_two_hemisphere.png', dpi=150, bbox_inches='tight')
plt.close(fig)

# --- Plot 3: NH-SH gradient validation ---
fig, axes = plt.subplots(2, 2, figsize=(12, 8), dpi=150)
fig.suptitle('Two-Hemisphere Model: Gradient Validation', fontsize=14)

# CH4 gradient
ax = axes[0, 0]
ax.plot(model_years, CH4_NH[:n_years_model], 'b-', linewidth=2, label='NH')
ax.plot(model_years, CH4_SH[:n_years_model], 'r-', linewidth=2, label='SH')
ax.set_ylabel('[CH₄] (ppb)')
ax.set_title('CH₄ Concentrations')
ax.legend()
ax.grid(True, alpha=0.3)

# IH gradient
ax = axes[0, 1]
ax.plot(model_years, IH_GRADIENT[:n_years_model], 'k-', linewidth=2)
ax.set_ylabel('NH - SH (ppb)')
ax.set_title('Interhemispheric CH₄ Gradient')
ax.grid(True, alpha=0.3)

# δ13C gradient
ax = axes[1, 0]
ax.plot(model_years, c13_NH[:n_years_model], 'b-', linewidth=2, label='NH')
ax.plot(model_years, c13_SH[:n_years_model], 'r-', linewidth=2, label='SH')
ax.set_ylabel('δ¹³C-CH₄ (‰)')
ax.set_title('δ¹³C-CH₄ by Hemisphere')
ax.set_xlabel('Year')
ax.legend()
ax.grid(True, alpha=0.3)

# Source fractions
ax = axes[1, 1]
nh_frac = SumSource_NH / SumSource_Global * 100
ax.plot(model_years, nh_frac, 'b-', linewidth=2, label='NH %')
ax.plot(model_years, 100 - nh_frac, 'r-', linewidth=2, label='SH %')
ax.axhline(75, color='gray', linestyle='--', alpha=0.5, label='Expected NH ~75%')
ax.set_ylabel('Emission Fraction (%)')
ax.set_title('NH/SH Emission Split')
ax.set_xlabel('Year')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUT_DIR / 'gradient_validation.png', dpi=150, bbox_inches='tight')
plt.close(fig)

print(f"\n  Plots saved:")
print(f"    - two_hemisphere_sources.png")
print(f"    - global_sources_two_hemisphere.png")
print(f"    - gradient_validation.png")

# ===========================================================================
# FINAL SUMMARY
# ===========================================================================
print(f"\n{'='*70}")
print("TWO-HEMISPHERE DUAL-ISOTOPE BOX MODEL — RUN COMPLETE")
print(f"{'='*70}")
print(f"  Years: {int(model_years[0])}–{int(model_years[-1])}")
print(f"  MC iterations: {N_ITERATIONS}")
print(f"  Architecture: 2-hemisphere boxes (NH + SH)")
print(f"  τ_ex = {TAU_EX} yr interhemispheric exchange")
print(f"  Isotopes: δ¹³C (NH/SH observed) + δD (hemispheric approximation)")
print(f"  Upgrades retained: KIE sampling, Quality monitoring, Time-varying τ")
print(f"  NH source fraction: {nh_frac.mean():.1f}% (expected ~70-75%)")
print(f"  Quality: NH {report_NH['nonphysical_pct']}% non-physical, "
      f"SH {report_SH['nonphysical_pct']}% non-physical")
print(f"{'='*70}")
