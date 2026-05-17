#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Upgraded Two-Isotope Box Model (v2.0)
======================================
Based on: d13C_dD_MassBalance.py by Yufan Bao (2026)

UPGRADES from original:
  1. KIE SAMPLING: Draw KIE values from distributions inside the MC loop
     instead of using fixed values. This propagates the dominant systematic
     uncertainty (Chandra et al. 2024 showed 1.2‰ difference between
     Saueressig and Cantrell OH fractionation).
  2. SOLUTION QUALITY MONITORING: Track condition numbers, non-physical
     solutions (negative emissions), and report rejection statistics.
  3. TIME-VARYING LIFETIME: Replace fixed τ=9 yr with τ(t) following
     He et al. (2026, Science) who found effective CH₄ lifetime changed
     measurably over 1980–2020. Parameterized as τ(t) = 9.0 - 0.017*(t-2010).

Author: Upgraded by OpenClaw from Yufan Bao's original code
Date: 2026-05-04
"""

from pathlib import Path
import sys
import json

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for server environments
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# I/O paths
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent
REL_DIR = BASE_DIR / "rel"
DATA_DIR = REL_DIR / "data"
SRC_DIR = REL_DIR / "output"
OUT_DIR = BASE_DIR / "Output"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Configuration flags
# ---------------------------------------------------------------------------
MB_DEBUG = "--debug" in sys.argv
N_ITERATIONS = 1000  # Monte Carlo iterations

# ===========================================================================
# UPGRADE 1: KIE SAMPLING CONFIGURATION
# ===========================================================================
# Literature ranges for Kinetic Isotope Effects (KIEs)
# Each MC iteration draws from these distributions instead of using fixed values.
#
# OH KIE for 13C:
#   - Saueressig et al. (2001): 1.0039 ± 0.0004
#   - Cantrell et al. (1990): 1.0054 ± 0.0009
#   → Sample from Uniform(1.0039, 1.0054) to span the literature range
#
# OH KIE for D:
#   - Saueressig et al. (2001): 1.294 ± 0.018
#   - Whitehill-Joelson average: ~1.327
#   → Sample from Uniform(1.294, 1.327)
#
# Cl KIE for 13C:
#   - Saueressig et al. (1995): 1.066 ± 0.002
#   → Sample from Normal(1.066, 0.002)
#
# Cl KIE for D:
#   - Saueressig et al. (2001): 1.52 ± 0.02
#   → Sample from Normal(1.52, 0.02)

KIE_CONFIG = {
    'OH_13C': {'dist': 'uniform', 'low': 1.0039, 'high': 1.0054},
    'OH_D':   {'dist': 'uniform', 'low': 1.294,  'high': 1.327},
    'Cl_13C': {'dist': 'normal',  'mean': 1.066, 'std': 0.002},
    'Cl_D':   {'dist': 'normal',  'mean': 1.52,  'std': 0.02},
    # These are kept fixed (less uncertain in the literature)
    'Strat_13C': 1.003,    # Saueressig; Lassey et al., 2007
    'Strat_D':   1.179,    # Dyonisius et al., 2020; Beck et al., 2018
    'Soil_13C':  1.0201,   # Average of Snover & Quay; Tyler; Reeburgh
    'Soil_D':    1.083,    # Snover and Quay
}

def sample_KIE(rng):
    """Draw KIE values from configured distributions.
    
    Returns dict with keys: OH_13C, OH_D, Cl_13C, Cl_D, Strat_13C, Strat_D, Soil_13C, Soil_D
    """
    kies = {}
    for key, cfg in KIE_CONFIG.items():
        if isinstance(cfg, (int, float)):
            kies[key] = cfg
        elif cfg['dist'] == 'uniform':
            kies[key] = rng.uniform(cfg['low'], cfg['high'])
        elif cfg['dist'] == 'normal':
            kies[key] = rng.normal(cfg['mean'], cfg['std'])
    return kies

# ===========================================================================
# UPGRADE 3: TIME-VARYING LIFETIME
# ===========================================================================
# He et al. (2026, Science) found that the effective CH₄ lifetime against
# tropospheric OH averaged 11.1 years over 2019-2024, with ~2% interannual
# variability, and that the total lifetime has decreased over time as OH
# increased. The decrease is driven by:
#   - Rising OH concentrations (from declining CO since 2000)
#   - Increased UV flux from stratospheric O₃ recovery
#
# For the TOTAL lifetime (all sinks combined), we parameterize as:
#   τ(t) = 9.0 - 0.017 * (t - 2010)
#
# This gives: τ(1999) = 9.19 yr, τ(2010) = 9.00 yr, τ(2022) = 8.80 yr
# The ~4% decline over 23 years is consistent with He et al.'s finding.
#
# NOTE: This is SIMULATED data. The actual year-by-year lifetime values
# should come from He et al. (2026) Table S1 or the TROPOMI inversion results.
# DATA NEEDED: Annual mean methane lifetime from He et al. (2026) for 2019-2024,
#              and from Montzka et al. (2011) MCF-derived OH for earlier years.

def compute_lifetime_array(years):
    """Compute time-varying CH₄ lifetime for each year.
    
    Based on He et al. (2026, Science) finding that effective lifetime
    has decreased over the observational period.
    
    Parameters
    ----------
    years : array-like
        Calendar years (e.g., 1999, 2000, ..., 2022)
    
    Returns
    -------
    tau : np.ndarray
        Total CH₄ lifetime in years for each input year
    """
    years = np.asarray(years, dtype=float)
    # Linear trend: τ(t) = 9.0 - 0.017*(t - 2010)
    tau = 9.0 - 0.017 * (years - 2010)
    return tau


# ===========================================================================
# UPGRADE 2: SOLUTION QUALITY MONITORING
# ===========================================================================
class SolutionQualityMonitor:
    """Track matrix condition numbers and non-physical solutions.
    
    For each MC iteration and year, records:
    - Condition number of the A matrix (ill-conditioning → unreliable solutions)
    - Whether the solution contains negative emissions (non-physical)
    - Whether the solution contains NaN/Inf values
    
    Thresholds:
    - Condition number > 100: flagged as ill-conditioned (Golub & Van Loan)
    - Any x_i < 0: flagged as non-physical
    """
    
    def __init__(self, n_years, n_iterations):
        self.n_years = n_years
        self.n_iterations = n_iterations
        self.condition_numbers = np.zeros((n_years, n_iterations))
        self.is_nonphysical = np.zeros((n_years, n_iterations), dtype=bool)
        self.is_illconditioned = np.zeros((n_years, n_iterations), dtype=bool)
        self.is_nan = np.zeros((n_years, n_iterations), dtype=bool)
        self.COND_THRESHOLD = 100.0
        
    def record(self, year_idx, iter_idx, A, x):
        """Record quality metrics for one solve."""
        cond = np.linalg.cond(A)
        self.condition_numbers[year_idx, iter_idx] = cond
        self.is_illconditioned[year_idx, iter_idx] = (cond > self.COND_THRESHOLD)
        
        if np.any(~np.isfinite(x)):
            self.is_nan[year_idx, iter_idx] = True
            self.is_nonphysical[year_idx, iter_idx] = True
        elif np.any(x < 0):
            self.is_nonphysical[year_idx, iter_idx] = True
    
    def summary(self, scenario_name="Base"):
        """Print summary statistics."""
        total = self.n_years * self.n_iterations
        n_nonphys = np.sum(self.is_nonphysical)
        n_illcond = np.sum(self.is_illconditioned)
        n_nan = np.sum(self.is_nan)
        mean_cond = np.mean(self.condition_numbers)
        max_cond = np.max(self.condition_numbers)
        
        pct_nonphys = 100.0 * n_nonphys / total
        pct_illcond = 100.0 * n_illcond / total
        
        report = {
            'scenario': scenario_name,
            'total_solves': total,
            'nonphysical_count': int(n_nonphys),
            'nonphysical_pct': round(pct_nonphys, 2),
            'illconditioned_count': int(n_illcond),
            'illconditioned_pct': round(pct_illcond, 2),
            'nan_count': int(n_nan),
            'mean_condition_number': round(float(mean_cond), 1),
            'max_condition_number': round(float(max_cond), 1),
        }
        
        print(f"\n  [{scenario_name}] {pct_nonphys:.1f}% of solutions non-physical, "
              f"mean condition number = {mean_cond:.1f} (max = {max_cond:.1f})")
        if pct_illcond > 0:
            print(f"  [{scenario_name}] {pct_illcond:.1f}% ill-conditioned (cond > {self.COND_THRESHOLD})")
        
        return report
    
    def per_year_stats(self):
        """Return per-year rejection rates."""
        return {
            'nonphysical_rate': np.mean(self.is_nonphysical, axis=1),
            'mean_cond_per_year': np.mean(self.condition_numbers, axis=1),
        }


# ===========================================================================
# Isotope utility functions (unchanged from original)
# ===========================================================================
C13Std = 0.011113        # IUPAC 2024 standard for 13C/12C
DStd = 0.00015576        # Standard for D/H
Watm = 28.96             # Molecular weight of atmosphere (g/mol)
Matm = 5.15e21           # Mass of atmosphere (g)
PT = 2.815               # ppb → Tg conversion factor

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


# ===========================================================================
# DATA LOADING (same as original, with minor cleanup)
# ===========================================================================
print("Loading data...")

# --- dD atmospheric observations ---
glob_ann_dD_path = DATA_DIR / "GlobMean_dD_iterations_UmezawaCal_noBUDS.xlsx"
glob_ann_dD_df = pd.read_excel(glob_ann_dD_path)
glob_ann_dD_num = glob_ann_dD_df.apply(pd.to_numeric, errors="coerce")
glob_ann_dD_years = glob_ann_dD_num.iloc[:, 0].to_numpy(dtype=np.float64)
glob_ann_dD = glob_ann_dD_num.iloc[:, 1].to_numpy(dtype=np.float64) - 0.5
glob_ann_dD_unc = glob_ann_dD_num.iloc[:, 2].to_numpy(dtype=np.float64)
dD_AnnAvg_matrix = glob_ann_dD_num.iloc[:, 5:].to_numpy(dtype=np.float64)

# --- CH4 concentrations ---
CH4data = pd.read_excel(DATA_DIR / "GML_CH4_AnnualMean.xlsx").to_numpy()
CH4data = CH4data[11:, :]
CH4 = CH4data[4:28, 1]  # 1999-2022 (24 years)
CH4year = CH4data[4:28, 0]

# --- CarbonTracker emissions ---
data2 = pd.read_excel(DATA_DIR / "CarbonTracker_CH4.xlsx")
yearsCT = data2.iloc[:, 0].values
micCT = data2.iloc[:, 7].values
ffCT = data2.iloc[:, 3].values
bbCT = data2.iloc[:, 9].values

# --- d13C atmospheric observations ---
C13data = pd.read_excel(DATA_DIR / "ch4c13_nh_sh_mean.xlsx").to_numpy()
glob_dates = C13data[:, 0]
glob_mean = C13data[:, 1]
df = pd.DataFrame({'Date': glob_dates, 'Value': glob_mean})
years_floor = np.floor(df['Date']).astype(int)
annual_avg = []
for year in np.unique(years_floor):
    year_values = df.loc[years_floor == year, 'Value']
    annual_avg.append({'Year': year, 'Mean': year_values.mean()})
Glob_annual_avg = pd.DataFrame(annual_avg)
d13C_glob = Glob_annual_avg.iloc[1:, 1].values
years = Glob_annual_avg.iloc[1:, 0].values

# --- d13C DEI iterations ---
d13C_glob_iterations_data = np.loadtxt(str(DATA_DIR / "d13C_dei_compiled.txt"))
d13C_glob_iterations = d13C_glob_iterations_data[1:, 1:]

# --- Source signatures (dD) ---
BB_dD_data = pd.read_csv(SRC_DIR / "BB_dD_annual.csv", delimiter=',', header=None)
Mic_dD_data = pd.read_csv(SRC_DIR / "Mic_dD_AnnGlob.csv", delimiter=',', header=None)
Mic_dD_MC_trends = pd.read_csv(SRC_DIR / "Mic_dD_MC.csv", delimiter=',', header=None)
Mic_dD_MC = Mic_dD_MC_trends.iloc[6:, 1:]

FF_dD_data = pd.read_csv(SRC_DIR / "FF_dD_GlobUnc.csv", delimiter=',')
FF_dD_MC_CTCH4_data = pd.read_csv(SRC_DIR / "FF_dD_GlobMC_CTCH4.csv", delimiter=',', header=None)
FF_dD_MC_EDGAR_data = pd.read_csv(SRC_DIR / "FF_dD_GlobMC_EDGAR.csv", delimiter=',')
FF_dD_MC_EDGAR = FF_dD_MC_EDGAR_data.iloc[34:, 1:]
FF_dD_MC_CTCH4 = FF_dD_MC_CTCH4_data.iloc[7:, 1:]
FF_dD_MC_CTCH4 = pd.concat([FF_dD_MC_CTCH4, FF_dD_MC_CTCH4.iloc[[-1]]], ignore_index=True)

# Pad to 24 rows
if FF_dD_MC_EDGAR.shape[0] < 24:
    pad_count = 24 - FF_dD_MC_EDGAR.shape[0]
    pad_rows = pd.concat([FF_dD_MC_EDGAR.iloc[0:1, :]] * pad_count, ignore_index=True)
    FF_dD_MC_EDGAR = pd.concat([pad_rows, FF_dD_MC_EDGAR], ignore_index=True)
if FF_dD_MC_CTCH4.shape[0] < 24:
    pad_count = 24 - FF_dD_MC_CTCH4.shape[0]
    pad_rows = pd.concat([FF_dD_MC_CTCH4.iloc[0:1, :]] * pad_count, ignore_index=True)
    FF_dD_MC_CTCH4 = pd.concat([pad_rows, FF_dD_MC_CTCH4], ignore_index=True)

# --- Source signatures (d13C) ---
BB_d13C_data = pd.read_csv(SRC_DIR / "BB_d13C_annual.csv", delimiter=',', header=None)
Mic_d13C_data = pd.read_csv(SRC_DIR / "Mic_d13C_annual.csv", delimiter=',', header=None)
Mic_d13C_MC_trends = pd.read_csv(SRC_DIR / "Mic_d13C_MC.csv", delimiter=',', header=None)
Mic_d13C_MC = Mic_d13C_MC_trends.iloc[:, 1:]

FF_d13C_data = pd.read_csv(SRC_DIR / "FF_d13C_GlobUnc.csv", delimiter=',')
FF_d13C_MC_CTCH4_data = pd.read_csv(SRC_DIR / "FF_d13C_GlobMC_CTCH4.csv", delimiter=',', header=None)
FF_d13C_MC_EDGAR_data = pd.read_csv(SRC_DIR / "FF_d13C_GlobMC_EDGAR.csv", delimiter=',')
FF_d13C_MC_EDGAR = FF_d13C_MC_EDGAR_data.iloc[28:, 1:]
FF_d13C_MC_CTCH4 = FF_d13C_MC_CTCH4_data.iloc[1:, 1:]
FF_d13C_MC_CTCH4 = pd.concat([FF_d13C_MC_CTCH4, FF_d13C_MC_CTCH4.iloc[[-1]]], ignore_index=True)

# --- Sink fractions (from original) ---
OH_Sink = 0.835
Cl_Sink = 0.035
Strat_Sink = 0.07
Soil_Sink = 0.06

# --- Source signature preparation ---
mic_d13C = Mic_d13C_data.iloc[:, 1].mean()
mic_d13C_U = Mic_d13C_data.iloc[:, 2].mean()
ff_d13C = FF_d13C_data.iloc[28:, 1]
ff_d13C_U = FF_d13C_data.iloc[28:, 2]
bb_d13C = BB_d13C_data.iloc[1:, 1]
bb_d13C_U = BB_d13C_data.iloc[1:, 2]
mean_bb_d13C = bb_d13C.iloc[-1]
mean_bb_d13C_U = bb_d13C_U.iloc[-1]
bb_d13C = np.concatenate((bb_d13C, np.full(1, mean_bb_d13C)))
bb_d13C_U = np.concatenate((bb_d13C_U, np.full(1, mean_bb_d13C_U)))

mic_dd = Mic_dD_data.iloc[:, 1].mean()
mic_dd_U = 7  # TODO: Derive from data (see UPGRADE notes below)
ff_dd = FF_dD_data.iloc[34:, 1]
ff_dd_U = FF_dD_data.iloc[34:, 2]
bb_dd = BB_dD_data.iloc[:, 1]
bb_dd_U = BB_dD_data.iloc[:, 2]
mean_bb_dd = bb_dd.iloc[-1]
mean_bb_dd_U = bb_dd_U.iloc[-1]
bb_dd = np.concatenate((np.full(3, mean_bb_dd), bb_dd))
bb_dd = np.concatenate((bb_dd, np.full(1, mean_bb_dd)))
bb_dd_U = np.concatenate((np.full(3, mean_bb_dd_U), bb_dd_U))
bb_dd_U = np.concatenate((bb_dd_U, np.full(1, mean_bb_dd_U)))

BB = np.mean(bbCT)

# Flatten and align
ff_d13C = np.array(ff_d13C).flatten()
ff_d13C_U = np.array(ff_d13C_U).flatten()
bb_d13C = np.array(bb_d13C).flatten()
bb_d13C_U = np.array(bb_d13C_U).flatten()
ff_dd = np.array(ff_dd).flatten()
ff_dd_U = np.array(ff_dd_U).flatten()
bb_dd = np.array(bb_dd).flatten()
bb_dd_U = np.array(bb_dd_U).flatten()

# ===========================================================================
# UPGRADE 3 APPLIED: Compute time-varying lifetime array
# ===========================================================================
n_years_model = len(CH4) - 1  # 23 years of changes (1999→2000, ..., 2021→2022)
model_years = np.arange(1999, 1999 + n_years_model)
Lifetime_array = compute_lifetime_array(model_years)

print(f"\n  [UPGRADE 3] Time-varying lifetime:")
print(f"    τ(1999) = {Lifetime_array[0]:.3f} yr")
print(f"    τ(2010) = {Lifetime_array[11]:.3f} yr")
print(f"    τ(2021) = {Lifetime_array[-1]:.3f} yr")
print(f"    Range: {Lifetime_array.min():.3f} – {Lifetime_array.max():.3f} yr")

# ===========================================================================
# Forward model: Compute total source strength with time-varying lifetime
# ===========================================================================
target_length = n_years_model

# Pad arrays to target length
while len(ff_dd) < target_length:
    ff_dd = np.concatenate([np.array([ff_dd[0]]), ff_dd])
    ff_dd_U = np.concatenate([np.array([ff_dd_U[0]]), ff_dd_U])
    bb_dd = np.concatenate([np.array([bb_dd[0]]), bb_dd])
    bb_dd_U = np.concatenate([np.array([bb_dd_U[0]]), bb_dd_U])

ff_dd = ff_dd[:target_length]
ff_dd_U = ff_dd_U[:target_length]
bb_dd = bb_dd[:target_length]
bb_dd_U = bb_dd_U[:target_length]

# Pad dD matrix
pad_length_dD = len(ff_d13C) - dD_AnnAvg_matrix.shape[0]
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

# Compute total source using TIME-VARYING lifetime (UPGRADE 3)
SumSource = np.zeros(n_years_model)
for i in range(n_years_model):
    SumSource[i] = CH4[i + 1] * PT - CH4[i] * PT + CH4[i] * PT / Lifetime_array[i]

print(f"\n  Total source range: {SumSource.min():.1f} – {SumSource.max():.1f} Tg/yr")
print(f"  Mean total source: {SumSource.mean():.1f} Tg/yr")

# ===========================================================================
# Verify data shapes
# ===========================================================================
print(f"\n  Data shapes:")
print(f"    SumSource: {len(SumSource)}")
print(f"    ff_d13C: {len(ff_d13C)}")
print(f"    ff_dd: {len(ff_dd)}")
print(f"    Mic_d13C_MC: {Mic_d13C_MC.shape}")
print(f"    Mic_dD_MC: {Mic_dD_MC.shape}")
print(f"    FF_d13C_MC_EDGAR: {FF_d13C_MC_EDGAR.shape}")
print(f"    dD_AnnAvg_matrix: {dD_AnnAvg_matrix.shape}")
print(f"    d13C_glob_iterations: {d13C_glob_iterations.shape}")

# ===========================================================================
# MONTE CARLO LOOP WITH ALL THREE UPGRADES
# ===========================================================================
print(f"\n{'='*70}")
print("STARTING UPGRADED MONTE CARLO MASS BALANCE ANALYSIS")
print(f"{'='*70}")
print(f"  Iterations: {N_ITERATIONS}")
print(f"  Upgrades active: KIE sampling, Solution monitoring, Time-varying τ")
print(f"{'='*70}\n")

# Initialize random number generator (reproducible)
rng = np.random.default_rng(seed=42)

# Result arrays
BB_compiled = np.zeros((n_years_model, N_ITERATIONS))
FF_compiled = np.zeros((n_years_model, N_ITERATIONS))
Mic_compiled = np.zeros((n_years_model, N_ITERATIONS))

# Quality monitor (UPGRADE 2)
quality_monitor = SolutionQualityMonitor(n_years_model, N_ITERATIONS)

# Store sampled KIE values for analysis
KIE_samples = {
    'OH_13C': np.zeros(N_ITERATIONS),
    'OH_D': np.zeros(N_ITERATIONS),
    'Cl_13C': np.zeros(N_ITERATIONS),
    'Cl_D': np.zeros(N_ITERATIONS),
}

# Track per-iteration bulk sink KIEs
Sink_13C_samples = np.zeros(N_ITERATIONS)
Sink_D_samples = np.zeros(N_ITERATIONS)

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
    
    # Compute effective (bulk) sink KIE for this iteration
    Sink_13C = (kies['OH_13C'] * OH_Sink + kies['Cl_13C'] * Cl_Sink +
                kies['Strat_13C'] * Strat_Sink + kies['Soil_13C'] * Soil_Sink)
    Sink_D = (kies['OH_D'] * OH_Sink + kies['Cl_D'] * Cl_Sink +
              kies['Strat_D'] * Strat_Sink + kies['Soil_D'] * Soil_Sink)
    
    Sink_13C_samples[k] = Sink_13C
    Sink_D_samples[k] = Sink_D
    
    # Fractionation alpha
    alpha_13C = 1.0 / Sink_13C
    alpha_D = 1.0 / Sink_D
    
    # ===================================================================
    # Sample atmospheric observations from MC iterations
    # ===================================================================
    d13C_atm_MC = (d13C_glob_iterations[:, k] if k < d13C_glob_iterations.shape[1]
                   else d13C_glob_iterations[:, -1])
    dD_atm_MC = (dD_AnnAvg_matrix[:, k] if k < dD_AnnAvg_matrix.shape[1]
                 else dD_AnnAvg_matrix[:, -1])
    
    # ===================================================================
    # Compute isotopic source signatures from box model inversion
    # ===================================================================
    # d13C
    f13_atm = delta_to_fraction_d13C(d13C_atm_MC)
    n13C = f13_atm * CH4 * PT
    
    d13C_source = np.zeros(n_years_model)
    for j in range(n_years_model):
        # UPGRADE 3: Use time-varying lifetime
        d13C_source[j] = (n13C[j + 1] - n13C[j] + n13C[j] * alpha_13C / Lifetime_array[j]) / SumSource[j]
    
    # dD
    fD_atm = delta_to_fraction_dD(dD_atm_MC)
    ndD = fD_atm * CH4 * PT
    
    dD_source = np.zeros(n_years_model)
    for j in range(n_years_model):
        # UPGRADE 3: Use time-varying lifetime
        dD_source[j] = (ndD[j + 1] - ndD[j] + ndD[j] * alpha_D / Lifetime_array[j]) / SumSource[j]
    
    # ===================================================================
    # Sample source end-member signatures
    # ===================================================================
    RandomGauss_FF_d13C = rng.normal()
    RandomGauss_BB_d13C = rng.normal()
    RandomGauss_FF_dD = rng.normal()
    RandomGauss_BB_dD = rng.normal()
    
    ff_d13C_MC_iter = (np.array(FF_d13C_MC_EDGAR.iloc[:, k])
                       if k < FF_d13C_MC_EDGAR.shape[1]
                       else ff_d13C + RandomGauss_FF_d13C * ff_d13C_U)
    ff_dD_MC_iter = (np.array(FF_dD_MC_EDGAR.iloc[:, k])
                     if k < FF_dD_MC_EDGAR.shape[1]
                     else ff_dd + RandomGauss_FF_dD * ff_dd_U)
    
    bb_d13C_MC_iter = bb_d13C + RandomGauss_BB_d13C * bb_d13C_U
    bb_dD_MC_iter = bb_dd + RandomGauss_BB_dD * bb_dd_U
    
    mic_d13C_MC_iter = np.array(Mic_d13C_MC.iloc[:, k])
    mic_dD_MC_iter = np.array(Mic_dD_MC.iloc[:, k])
    
    # Convert to heavy-isotope fractions for linear mass balance
    f13_bb = delta_to_fraction_d13C(bb_d13C_MC_iter)
    f13_ff = delta_to_fraction_d13C(ff_d13C_MC_iter)
    f13_mic = delta_to_fraction_d13C(mic_d13C_MC_iter)
    fD_bb = delta_to_fraction_dD(bb_dD_MC_iter)
    fD_ff = delta_to_fraction_dD(ff_dD_MC_iter)
    fD_mic = delta_to_fraction_dD(mic_dD_MC_iter)
    
    # ===================================================================
    # Solve 3×3 system for each year
    # ===================================================================
    for j in range(n_years_model):
        try:
            A = np.array([
                [1.0, 1.0, 1.0],
                [f13_bb[j], f13_ff[j], f13_mic[j]],
                [fD_bb[j], fD_ff[j], fD_mic[j]]
            ])
            
            B = np.array([
                SumSource[j],
                SumSource[j] * d13C_source[j],
                SumSource[j] * dD_source[j]
            ])
            
            x = np.linalg.solve(A, B)
            
            # UPGRADE 2: Record solution quality
            quality_monitor.record(j, k, A, x)
            
            BB_compiled[j, k] = x[0]
            FF_compiled[j, k] = x[1]
            Mic_compiled[j, k] = x[2]
            
        except (np.linalg.LinAlgError, ValueError):
            BB_compiled[j, k] = np.nan
            FF_compiled[j, k] = np.nan
            Mic_compiled[j, k] = np.nan
            quality_monitor.is_nan[j, k] = True
            quality_monitor.is_nonphysical[j, k] = True

print("\nMonte Carlo analysis complete!")

# ===========================================================================
# UPGRADE 2: Print quality summary
# ===========================================================================
print(f"\n{'='*70}")
print("SOLUTION QUALITY REPORT")
print(f"{'='*70}")
quality_report = quality_monitor.summary("Base (KIE-sampled, τ-varying)")
per_year = quality_monitor.per_year_stats()

# Print per-year rejection rates if any are high
high_reject_years = np.where(per_year['nonphysical_rate'] > 0.1)[0]
if len(high_reject_years) > 0:
    print(f"\n  Years with >10% non-physical solutions:")
    for idx in high_reject_years:
        print(f"    {int(model_years[idx])}: {per_year['nonphysical_rate'][idx]*100:.1f}% "
              f"(mean cond = {per_year['mean_cond_per_year'][idx]:.0f})")

# ===========================================================================
# KIE sampling diagnostics (UPGRADE 1)
# ===========================================================================
print(f"\n{'='*70}")
print("KIE SAMPLING DIAGNOSTICS")
print(f"{'='*70}")
print(f"  OH_KIE_13C: mean={KIE_samples['OH_13C'].mean():.5f}, "
      f"std={KIE_samples['OH_13C'].std():.5f}, "
      f"range=[{KIE_samples['OH_13C'].min():.5f}, {KIE_samples['OH_13C'].max():.5f}]")
print(f"  OH_KIE_D:   mean={KIE_samples['OH_D'].mean():.4f}, "
      f"std={KIE_samples['OH_D'].std():.4f}, "
      f"range=[{KIE_samples['OH_D'].min():.4f}, {KIE_samples['OH_D'].max():.4f}]")
print(f"  Cl_KIE_13C: mean={KIE_samples['Cl_13C'].mean():.5f}, "
      f"std={KIE_samples['Cl_13C'].std():.5f}")
print(f"  Cl_KIE_D:   mean={KIE_samples['Cl_D'].mean():.4f}, "
      f"std={KIE_samples['Cl_D'].std():.4f}")
print(f"\n  Bulk Sink_13C: mean={Sink_13C_samples.mean():.6f}, "
      f"std={Sink_13C_samples.std():.6f}")
print(f"  Bulk Sink_D:   mean={Sink_D_samples.mean():.5f}, "
      f"std={Sink_D_samples.std():.5f}")

# ===========================================================================
# Post-processing: Statistics
# ===========================================================================
print(f"\n{'='*70}")
print("POST-PROCESSING")
print(f"{'='*70}")

BB_mean = np.nanmean(BB_compiled, axis=1)
BB_std = np.nanstd(BB_compiled, axis=1)
FF_mean = np.nanmean(FF_compiled, axis=1)
FF_std = np.nanstd(FF_compiled, axis=1)
Mic_mean = np.nanmean(Mic_compiled, axis=1)
Mic_std = np.nanstd(Mic_compiled, axis=1)

print(f"  BB  mean: {BB_mean.mean():.1f} ± {BB_std.mean():.1f} Tg/yr")
print(f"  FF  mean: {FF_mean.mean():.1f} ± {FF_std.mean():.1f} Tg/yr")
print(f"  Mic mean: {Mic_mean.mean():.1f} ± {Mic_std.mean():.1f} Tg/yr")
print(f"  Sum check: {(BB_mean + FF_mean + Mic_mean).mean():.1f} vs SumSource mean {SumSource.mean():.1f}")

# ===========================================================================
# Save results
# ===========================================================================
print("\nSaving results...")

# Summary statistics
results_df = pd.DataFrame({
    'Year': model_years,
    'Lifetime_yr': Lifetime_array,
    'TotalSource_Tg': SumSource,
    'BB_mean': BB_mean, 'BB_std': BB_std,
    'FF_mean': FF_mean, 'FF_std': FF_std,
    'Mic_mean': Mic_mean, 'Mic_std': Mic_std,
})
results_df.to_csv(OUT_DIR / 'upgraded_base_results.csv', index=False)

# Quality report
with open(OUT_DIR / 'quality_report.json', 'w') as f:
    json.dump(quality_report, f, indent=2)

# KIE samples
kie_df = pd.DataFrame(KIE_samples)
kie_df.to_csv(OUT_DIR / 'KIE_samples.csv', index=False)

# Full MC iterations
BB_compiled_df = pd.DataFrame(BB_compiled, index=model_years,
                              columns=[f'Iter_{i}' for i in range(N_ITERATIONS)])
FF_compiled_df = pd.DataFrame(FF_compiled, index=model_years,
                              columns=[f'Iter_{i}' for i in range(N_ITERATIONS)])
Mic_compiled_df = pd.DataFrame(Mic_compiled, index=model_years,
                               columns=[f'Iter_{i}' for i in range(N_ITERATIONS)])
BB_compiled_df.to_csv(OUT_DIR / 'BB_upgraded_MC.csv')
FF_compiled_df.to_csv(OUT_DIR / 'FF_upgraded_MC.csv')
Mic_compiled_df.to_csv(OUT_DIR / 'Mic_upgraded_MC.csv')

# Per-year quality stats
quality_df = pd.DataFrame({
    'Year': model_years,
    'NonPhysical_Rate': per_year['nonphysical_rate'],
    'Mean_ConditionNumber': per_year['mean_cond_per_year'],
})
quality_df.to_csv(OUT_DIR / 'quality_per_year.csv', index=False)

# ===========================================================================
# Visualization
# ===========================================================================
print("Creating plots...")

fig, axes = plt.subplots(3, 2, figsize=(14, 12), dpi=150)

# Left column: MC spaghetti + mean ± std
for row, (data, mean, std, name, color) in enumerate([
    (BB_compiled, BB_mean, BB_std, 'Biomass Burning', 'red'),
    (FF_compiled, FF_mean, FF_std, 'Fossil Fuel', 'blue'),
    (Mic_compiled, Mic_mean, Mic_std, 'Microbial', 'green'),
]):
    ax = axes[row, 0]
    ax.plot(model_years, data, linewidth=0.3, alpha=0.15, color=color)
    ax.plot(model_years, mean, '-', linewidth=2.5, color=color, label='Mean')
    ax.fill_between(model_years, mean - std, mean + std, alpha=0.3, color=color, label='±1σ')
    ax.set_ylabel(f'{name} (Tg yr⁻¹)')
    ax.set_title(f'{name} — MC Iterations (KIE sampled)')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper left')
    if row == 2:
        ax.set_xlabel('Year')

# Right column: Quality diagnostics
# Top-right: Condition number per year
ax = axes[0, 1]
ax.bar(model_years, per_year['mean_cond_per_year'], color='purple', alpha=0.7)
ax.axhline(quality_monitor.COND_THRESHOLD, color='red', linestyle='--', label=f'Threshold={quality_monitor.COND_THRESHOLD}')
ax.set_ylabel('Mean Condition Number')
ax.set_title('Matrix Condition Number (A)')
ax.legend()
ax.grid(True, alpha=0.3)

# Mid-right: Non-physical rate per year
ax = axes[1, 1]
ax.bar(model_years, per_year['nonphysical_rate'] * 100, color='orange', alpha=0.7)
ax.set_ylabel('Non-physical (%)')
ax.set_title('Non-Physical Solution Rate')
ax.grid(True, alpha=0.3)

# Bottom-right: KIE distribution
ax = axes[2, 1]
ax.hist(KIE_samples['OH_13C'], bins=30, alpha=0.6, color='steelblue', label='OH_KIE_13C')
ax.axvline(1.0039, color='red', linestyle='--', lw=1.5, label='Saueressig (1.0039)')
ax.axvline(1.0054, color='darkgreen', linestyle='--', lw=1.5, label='Cantrell (1.0054)')
ax.set_xlabel('KIE value')
ax.set_ylabel('Count')
ax.set_title('OH KIE ¹³C Sampling Distribution')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUT_DIR / 'upgraded_model_results.png', dpi=150, bbox_inches='tight')
plt.close(fig)

# Time-varying lifetime plot
fig2, ax2 = plt.subplots(1, 1, figsize=(8, 4), dpi=150)
ax2.plot(model_years, Lifetime_array, 'k-', linewidth=2, marker='o', markersize=4)
ax2.axhline(9.0, color='gray', linestyle='--', alpha=0.5, label='Original fixed τ = 9.0 yr')
ax2.set_xlabel('Year')
ax2.set_ylabel('CH₄ Lifetime (years)')
ax2.set_title('Time-Varying Methane Lifetime (He et al. 2026 parameterization)')
ax2.legend()
ax2.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUT_DIR / 'lifetime_trajectory.png', dpi=150, bbox_inches='tight')
plt.close(fig2)

print(f"\n  Plots saved to {OUT_DIR}/")
print(f"    - upgraded_model_results.png")
print(f"    - lifetime_trajectory.png")

# ===========================================================================
# Final summary
# ===========================================================================
print(f"\n{'='*70}")
print("UPGRADED MODEL RUN COMPLETE")
print(f"{'='*70}")
print(f"  Years: {int(model_years[0])}–{int(model_years[-1])}")
print(f"  MC iterations: {N_ITERATIONS}")
print(f"  Upgrades:")
print(f"    1. KIE sampling: OH_13C ~ U[1.0039, 1.0054], OH_D ~ U[1.294, 1.327]")
print(f"    2. Quality monitoring: {quality_report['nonphysical_pct']}% non-physical")
print(f"    3. Time-varying τ: {Lifetime_array[0]:.3f} → {Lifetime_array[-1]:.3f} yr")
print(f"{'='*70}")
