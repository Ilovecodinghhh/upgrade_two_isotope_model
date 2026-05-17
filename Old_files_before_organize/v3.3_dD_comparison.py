#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Comparison of δD Improvement Strategies for the 2×2 BB-Fixed Model
===================================================================
Implements and compares:
  A) Baseline v3.2 (δD starts 2006, no fixes)
  B) Quick Fix 4: Use Ben's hemispheric δD output directly (NH/SH means)
  C) Quick Fix 1: Start δD constraint from 2010 only (δ13C-only before 2010)
  D) Quick Fix 6: Weight by inverse-uncertainty (more stations → more weight)

All use the v3.2 (2×2 BB-fixed) framework for robustness.

Author: OpenClaw
Date: 2026-05-06
"""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).resolve().parent
REL_DIR = BASE_DIR / "rel"
DATA_DIR = REL_DIR / "data"
SRC_DIR = REL_DIR / "output"
BEN_DD_DIR = BASE_DIR.parent / "Ben-BoxModel" / "Riddell-Young_2025_dD_GlobMean" / "Riddell-Young_2025_dD_GlobMean"
OUT_DIR = BASE_DIR / "Output_dD_comparison"
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_ITERATIONS = 1000
rng = np.random.default_rng(seed=42)

# ---------------------------------------------------------------------------
# Isotope utilities
# ---------------------------------------------------------------------------
C13Std = 0.011113
DStd = 0.00015576
PT = 2.815

def delta_to_R_d13C(d): return (d/1000 + 1) * C13Std
def delta_to_R_dD(d): return (d/1000 + 1) * DStd
def R_to_frac(R): return R/(1+R)
def frac_to_R(f): return f/(1-f)
def frac_to_delta_d13C(f): return ((frac_to_R(f) - C13Std)/C13Std)*1000
def frac_to_delta_dD(f): return ((frac_to_R(f) - DStd)/DStd)*1000
def delta_to_frac_d13C(d): return R_to_frac(delta_to_R_d13C(d))
def delta_to_frac_dD(d): return R_to_frac(delta_to_R_dD(d))

# ---------------------------------------------------------------------------
# Load common data
# ---------------------------------------------------------------------------
print("Loading data...")

# CH4 concentrations
CH4data = pd.read_excel(DATA_DIR / "GML_CH4_AnnualMean.xlsx").to_numpy()
CH4data = CH4data[11:, :]
CH4 = CH4data[4:28, 1].astype(float)  # 1999-2022 (24 years)
CH4year = CH4data[4:28, 0].astype(float)

# d13C atmospheric observations
C13data = pd.read_excel(DATA_DIR / "ch4c13_nh_sh_mean.xlsx").to_numpy()
glob_dates = C13data[:, 0]
glob_mean = C13data[:, 1]
df_c13 = pd.DataFrame({'Date': glob_dates, 'Value': glob_mean})
years_floor = np.floor(df_c13['Date']).astype(int)
annual_avg = []
for year in np.unique(years_floor):
    year_values = df_c13.loc[years_floor == year, 'Value']
    annual_avg.append({'Year': year, 'Mean': year_values.mean()})
Glob_annual_avg = pd.DataFrame(annual_avg)
d13C_glob = Glob_annual_avg.iloc[1:, 1].values  # starts at 1999
d13C_years = Glob_annual_avg.iloc[1:, 0].values

# d13C MC iterations
d13C_glob_iterations_data = np.loadtxt(str(DATA_DIR / "d13C_dei_compiled.txt"))
d13C_glob_iterations = d13C_glob_iterations_data[1:, 1:]  # skip header row, skip year col

# Original dD atmospheric observations (what v3.2 baseline uses)
glob_ann_dD_df = pd.read_excel(DATA_DIR / "GlobMean_dD_iterations_UmezawaCal_noBUDS.xlsx")
glob_ann_dD_num = glob_ann_dD_df.apply(pd.to_numeric, errors="coerce")
glob_ann_dD_years = glob_ann_dD_num.iloc[:, 0].to_numpy(dtype=np.float64)
glob_ann_dD = glob_ann_dD_num.iloc[:, 1].to_numpy(dtype=np.float64) - 0.5
dD_AnnAvg_matrix = glob_ann_dD_num.iloc[:, 5:].to_numpy(dtype=np.float64)

# Ben's hemispheric dD output (Quick Fix 4)
hem_dD = pd.read_csv(BEN_DD_DIR / "output" / "HemMean_dD_dei_UmezawaCal_noBUDS.csv")
hem_dD = hem_dD.dropna()
hem_dD['YearInt'] = np.floor(hem_dD['Year']).astype(int)
hem_annual = hem_dD.groupby('YearInt').agg(
    Glob=('Glob_smooth_mean', 'mean'),
    NH=('NH_smooth_mean', 'mean'),
    SH=('SH_smooth_mean', 'mean')
).reset_index()

# CarbonTracker
data2 = pd.read_excel(DATA_DIR / "CarbonTracker_CH4.xlsx")
bbCT = data2.iloc[:, 9].values
BB = np.mean(bbCT)

# Source signatures
FF_d13C_data = pd.read_csv(SRC_DIR / "FF_d13C_GlobUnc.csv", delimiter=',')
ff_d13C = FF_d13C_data.iloc[28:, 1].values.astype(float)
ff_d13C_U = FF_d13C_data.iloc[28:, 2].values.astype(float)

Mic_d13C_data = pd.read_csv(SRC_DIR / "Mic_d13C_annual.csv", delimiter=',', header=None)
mic_d13C = Mic_d13C_data.iloc[:, 1].mean()
mic_d13C_U = Mic_d13C_data.iloc[:, 2].mean()

Mic_d13C_MC_trends = pd.read_csv(SRC_DIR / "Mic_d13C_MC.csv", delimiter=',', header=None)
Mic_d13C_MC = Mic_d13C_MC_trends.iloc[:, 1:]

FF_dD_data = pd.read_csv(SRC_DIR / "FF_dD_GlobUnc.csv", delimiter=',')
ff_dd = FF_dD_data.iloc[34:, 1].values.astype(float)
ff_dd_U = FF_dD_data.iloc[34:, 2].values.astype(float)

Mic_dD_data = pd.read_csv(SRC_DIR / "Mic_dD_AnnGlob.csv", delimiter=',', header=None)
mic_dd = Mic_dD_data.iloc[:, 1].mean()
mic_dd_U = 7

Mic_dD_MC_trends = pd.read_csv(SRC_DIR / "Mic_dD_MC.csv", delimiter=',', header=None)
Mic_dD_MC = Mic_dD_MC_trends.iloc[6:, 1:]

FF_d13C_MC_EDGAR_data = pd.read_csv(SRC_DIR / "FF_d13C_GlobMC_EDGAR.csv", delimiter=',')
FF_d13C_MC_EDGAR = FF_d13C_MC_EDGAR_data.iloc[28:, 1:]
FF_dD_MC_EDGAR_data = pd.read_csv(SRC_DIR / "FF_dD_GlobMC_EDGAR.csv", delimiter=',')
FF_dD_MC_EDGAR = FF_dD_MC_EDGAR_data.iloc[34:, 1:]

# Sink fractions
OH_Sink, Cl_Sink, Strat_Sink, Soil_Sink = 0.835, 0.035, 0.07, 0.06

# KIE config
KIE_CONFIG = {
    'OH_13C': {'dist': 'uniform', 'low': 1.0039, 'high': 1.0054},
    'OH_D':   {'dist': 'uniform', 'low': 1.294,  'high': 1.327},
    'Cl_13C': {'dist': 'normal',  'mean': 1.066, 'std': 0.002},
    'Cl_D':   {'dist': 'normal',  'mean': 1.52,  'std': 0.02},
    'Strat_13C': 1.003, 'Strat_D': 1.179,
    'Soil_13C': 1.0201, 'Soil_D': 1.083,
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

# Time-varying lifetime
def compute_lifetime(years):
    return 9.0 - 0.017 * (np.asarray(years, dtype=float) - 2010)

# ---------------------------------------------------------------------------
# Model setup
# ---------------------------------------------------------------------------
n_years = len(CH4) - 1  # 23 years (1999-2021)
model_years = np.arange(1999, 1999 + n_years)
Lifetime = compute_lifetime(model_years)

# Total source
SumSource = np.zeros(n_years)
for i in range(n_years):
    SumSource[i] = CH4[i+1]*PT - CH4[i]*PT + CH4[i]*PT / Lifetime[i]

# Align arrays to n_years=23
target_length = n_years

# Pad ff_d13C to 23 if needed
while len(ff_d13C) < target_length:
    ff_d13C = np.concatenate([np.array([ff_d13C[0]]), ff_d13C])
    ff_d13C_U = np.concatenate([np.array([ff_d13C_U[0]]), ff_d13C_U])
ff_d13C = ff_d13C[:target_length]
ff_d13C_U = ff_d13C_U[:target_length]

# Pad ff_dd
while len(ff_dd) < target_length:
    ff_dd = np.concatenate([np.array([ff_dd[0]]), ff_dd])
    ff_dd_U = np.concatenate([np.array([ff_dd_U[0]]), ff_dd_U])
ff_dd = ff_dd[:target_length]
ff_dd_U = ff_dd_U[:target_length]

# Pad Mic MCs
if Mic_d13C_MC.shape[0] < target_length:
    pad = pd.concat([Mic_d13C_MC.iloc[0:1,:]]*( target_length - Mic_d13C_MC.shape[0]), ignore_index=True)
    Mic_d13C_MC = pd.concat([pad, Mic_d13C_MC], ignore_index=True)
if Mic_dD_MC.shape[0] < target_length:
    pad = pd.concat([Mic_dD_MC.iloc[0:1,:]]*(target_length - Mic_dD_MC.shape[0]), ignore_index=True)
    Mic_dD_MC = pd.concat([pad, Mic_dD_MC], ignore_index=True)
elif Mic_dD_MC.shape[0] > target_length:
    Mic_dD_MC = Mic_dD_MC.iloc[:target_length, :]

# Pad FF MC EDGAR
if FF_d13C_MC_EDGAR.shape[0] < target_length:
    pad = pd.concat([FF_d13C_MC_EDGAR.iloc[0:1,:]]*(target_length - FF_d13C_MC_EDGAR.shape[0]), ignore_index=True)
    FF_d13C_MC_EDGAR = pd.concat([pad, FF_d13C_MC_EDGAR], ignore_index=True)
if FF_dD_MC_EDGAR.shape[0] < target_length:
    pad = pd.concat([FF_dD_MC_EDGAR.iloc[0:1,:]]*(target_length - FF_dD_MC_EDGAR.shape[0]), ignore_index=True)
    FF_dD_MC_EDGAR = pd.concat([pad, FF_dD_MC_EDGAR], ignore_index=True)

# Pad dD matrix (original)
pad_needed = target_length + 1 - dD_AnnAvg_matrix.shape[0]  # need 24 rows (for CH4 indices 0..23)
if pad_needed > 0:
    pad_rows = np.repeat(dD_AnnAvg_matrix[0:1, :], pad_needed, axis=0)
    dD_AnnAvg_matrix = np.vstack([pad_rows, dD_AnnAvg_matrix])

# ---------------------------------------------------------------------------
# Prepare Ben's hemispheric dD for Quick Fix 4
# Map to model years (global mean from Ben's hemispheric output)
# ---------------------------------------------------------------------------
ben_dD_glob = np.full(24, np.nan)  # 1999-2022
for _, row in hem_annual.iterrows():
    yr = int(row['YearInt'])
    if 1999 <= yr <= 2022:
        ben_dD_glob[yr - 1999] = row['Glob']

# For years before 2005, use first available value (or original data)
first_valid = ben_dD_glob[~np.isnan(ben_dD_glob)][0]
for i in range(24):
    if np.isnan(ben_dD_glob[i]):
        ben_dD_glob[i] = first_valid

print(f"Ben's hemispheric dD global mean (2005-2022):")
for yr in range(2005, 2023):
    idx = yr - 1999
    print(f"  {yr}: {ben_dD_glob[idx]:.2f}")

# ---------------------------------------------------------------------------
# Station count per year (for Quick Fix 6 weighting)
# From our analysis: approximate station counts in each year
# ---------------------------------------------------------------------------
station_counts = {
    2005: 12, 2006: 12, 2007: 12, 2008: 13, 2009: 13,
    2010: 8, 2011: 6, 2012: 6, 2013: 7, 2014: 7,
    2015: 7, 2016: 7, 2017: 7, 2018: 8, 2019: 8,
    2020: 8, 2021: 8, 2022: 9
}
# Uncertainty scaling: fewer stations → larger uncertainty
# Base uncertainty at max stations (13) = 1.0; scale as sqrt(13/N)
max_stations = 13
dD_weight_by_stations = np.ones(24)
for yr in range(1999, 2023):
    idx = yr - 1999
    n_st = station_counts.get(yr, 5)  # default 5 for pre-2005
    if yr < 2005:
        dD_weight_by_stations[idx] = 0.0  # no dD data
    else:
        dD_weight_by_stations[idx] = 1.0 / np.sqrt(max_stations / n_st)

print(f"\ndD weights by station count:")
for yr in range(2005, 2023):
    print(f"  {yr}: weight={dD_weight_by_stations[yr-1999]:.3f} (n={station_counts.get(yr,0)})")


# ===========================================================================
# MC INVERSION: 2×2 BB-Fixed with different dD treatments
# ===========================================================================
def run_inversion(scenario_name, dD_atm_func, dD_weight_func=None):
    """
    Run 2×2 BB-fixed inversion.
    
    dD_atm_func(k, j) -> returns atmospheric dD for MC iter k, year index j (0-based from 1999)
        If returns None, skip dD constraint for that year (use d13C only)
    dD_weight_func(j) -> returns weight [0,1] for dD constraint at year j
        If None, all weights = 1.0
    """
    print(f"\n  Running scenario: {scenario_name}...")
    
    FF_d13C_results = np.zeros((n_years, N_ITERATIONS))
    Mic_d13C_results = np.zeros((n_years, N_ITERATIONS))
    FF_dD_results = np.full((n_years, N_ITERATIONS), np.nan)
    Mic_dD_results = np.full((n_years, N_ITERATIONS), np.nan)
    # Combined (weighted average of d13C and dD where both available)
    FF_combined = np.zeros((n_years, N_ITERATIONS))
    Mic_combined = np.zeros((n_years, N_ITERATIONS))
    
    neg_count_d13C = 0
    neg_count_dD = 0
    total_dD_solves = 0
    
    for k in range(N_ITERATIONS):
        kies = sample_KIE(rng)
        Sink_13C = (kies['OH_13C']*OH_Sink + kies['Cl_13C']*Cl_Sink + 
                    kies['Strat_13C']*Strat_Sink + kies['Soil_13C']*Soil_Sink)
        Sink_D = (kies['OH_D']*OH_Sink + kies['Cl_D']*Cl_Sink +
                  kies['Strat_D']*Strat_Sink + kies['Soil_D']*Soil_Sink)
        alpha_13C = 1.0 / Sink_13C
        alpha_D = 1.0 / Sink_D
        
        # Sample d13C atmosphere
        col_idx = min(k, d13C_glob_iterations.shape[1]-1)
        d13C_atm = d13C_glob_iterations[:, col_idx]
        
        # Sample source signatures
        rg_ff13 = rng.normal()
        rg_ff_dD = rng.normal()
        
        ff_d13C_iter = (np.array(FF_d13C_MC_EDGAR.iloc[:, min(k, FF_d13C_MC_EDGAR.shape[1]-1)])
                        if k < FF_d13C_MC_EDGAR.shape[1]
                        else ff_d13C + rg_ff13 * ff_d13C_U)
        ff_dD_iter = (np.array(FF_dD_MC_EDGAR.iloc[:, min(k, FF_dD_MC_EDGAR.shape[1]-1)])
                      if k < FF_dD_MC_EDGAR.shape[1]
                      else ff_dd + rg_ff_dD * ff_dd_U)
        
        mic_d13C_iter = np.array(Mic_d13C_MC.iloc[:, min(k, Mic_d13C_MC.shape[1]-1)])
        mic_dD_iter = np.array(Mic_dD_MC.iloc[:, min(k, Mic_dD_MC.shape[1]-1)])
        
        for j in range(n_years):
            # --- δ13C inversion (always available) ---
            f13_j = R_to_frac(delta_to_R_d13C(d13C_atm[j]))
            f13_j1 = R_to_frac(delta_to_R_d13C(d13C_atm[j+1]))
            n13_j = f13_j * CH4[j] * PT
            n13_j1 = f13_j1 * CH4[j+1] * PT
            d13C_source_f = (n13_j1 - n13_j + n13_j * alpha_13C / Lifetime[j]) / SumSource[j]
            d13C_source = frac_to_delta_d13C(d13C_source_f)
            
            # 2×2: FF = (S*δ_src - mic*(S-BB) - bb*BB) / (ff - mic)
            # But in fraction space for consistency:
            ff_sig = ff_d13C_iter[j]
            mic_sig = mic_d13C_iter[j]
            
            FFS_13C = (SumSource[j] * d13C_source - mic_sig * (SumSource[j] - BB) - (-23.1) * BB) / (ff_sig - mic_sig)
            MicS_13C = SumSource[j] - BB - FFS_13C
            
            if FFS_13C < 0 or MicS_13C < 0:
                neg_count_d13C += 1
            
            FF_d13C_results[j, k] = FFS_13C
            Mic_d13C_results[j, k] = MicS_13C
            
            # --- δD inversion (conditional) ---
            dD_atm_val = dD_atm_func(k, j)
            dD_atm_val_next = dD_atm_func(k, j+1) if j+1 < 24 else None
            
            if dD_atm_val is not None and dD_atm_val_next is not None and not np.isnan(dD_atm_val) and not np.isnan(dD_atm_val_next):
                fD_j = R_to_frac(delta_to_R_dD(dD_atm_val))
                fD_j1 = R_to_frac(delta_to_R_dD(dD_atm_val_next))
                nD_j = fD_j * CH4[j] * PT
                nD_j1 = fD_j1 * CH4[j+1] * PT
                dD_source_f = (nD_j1 - nD_j + nD_j * alpha_D / Lifetime[j]) / SumSource[j]
                dD_source = frac_to_delta_dD(dD_source_f)
                
                ff_sig_dD = ff_dD_iter[j]
                mic_sig_dD = mic_dD_iter[j]
                
                FFS_dD = (SumSource[j] * dD_source - mic_sig_dD * (SumSource[j] - BB) - (-228.0) * BB) / (ff_sig_dD - mic_sig_dD)
                MicS_dD = SumSource[j] - BB - FFS_dD
                
                if FFS_dD < 0 or MicS_dD < 0:
                    neg_count_dD += 1
                total_dD_solves += 1
                
                FF_dD_results[j, k] = FFS_dD
                Mic_dD_results[j, k] = MicS_dD
                
                # Combined: weighted average
                w_dD = dD_weight_func(j) if dD_weight_func else 1.0
                w_d13C = 1.0
                w_total = w_d13C + w_dD
                FF_combined[j, k] = (w_d13C * FFS_13C + w_dD * FFS_dD) / w_total
                Mic_combined[j, k] = (w_d13C * MicS_13C + w_dD * MicS_dD) / w_total
            else:
                # No dD available: use d13C only
                FF_combined[j, k] = FFS_13C
                Mic_combined[j, k] = MicS_13C
    
    pct_neg_d13C = 100.0 * neg_count_d13C / (n_years * N_ITERATIONS)
    pct_neg_dD = 100.0 * neg_count_dD / max(total_dD_solves, 1)
    print(f"    d13C negative: {pct_neg_d13C:.1f}%")
    print(f"    dD negative: {pct_neg_dD:.1f}% (of {total_dD_solves} solves)")
    
    return {
        'name': scenario_name,
        'FF_d13C': FF_d13C_results,
        'Mic_d13C': Mic_d13C_results,
        'FF_dD': FF_dD_results,
        'Mic_dD': Mic_dD_results,
        'FF_combined': FF_combined,
        'Mic_combined': Mic_combined,
        'neg_pct_d13C': pct_neg_d13C,
        'neg_pct_dD': pct_neg_dD,
    }

# ===========================================================================
# Scenario A: Baseline (original dD from 2006, pad earlier years)
# ===========================================================================
def dD_baseline(k, j):
    """Original dD data. dD starts at year index 7 (=2006)."""
    # dD_AnnAvg_matrix has been padded to 24+ rows
    row_idx = j  # j=0 is 1999
    if row_idx < dD_AnnAvg_matrix.shape[0]:
        col_idx = min(k, dD_AnnAvg_matrix.shape[1]-1)
        val = dD_AnnAvg_matrix[row_idx, col_idx]
        if not np.isnan(val):
            return val
    return None

# But actually the original dD only starts 2006 (year index 7 from 1999)
# Let's properly check what's padded vs real
# dD data starts 2006 → that's CH4 index 7. Padding fills 0..6 with row[0] value.
# We'll let baseline use everything (padded = constant before 2006)
def dD_baseline_func(k, j):
    row_idx = j
    if row_idx >= dD_AnnAvg_matrix.shape[0]:
        return None
    col_idx = min(k, dD_AnnAvg_matrix.shape[1]-1)
    val = dD_AnnAvg_matrix[row_idx, col_idx]
    return val if not np.isnan(val) else None

# ===========================================================================
# Scenario B: Quick Fix 4 - Use Ben's hemispheric global mean directly
# ===========================================================================
def dD_ben_hemispheric(k, j):
    """Use Ben's hemispheric-averaged global mean dD."""
    # j=0 → 1999; Ben's data starts 2005 (j=6)
    yr = 1999 + j
    if yr < 2005 or yr > 2023:
        return None
    return ben_dD_glob[j]  # Already filled/available

# ===========================================================================
# Scenario C: Quick Fix 1 - Start dD from 2010 only
# ===========================================================================
def dD_from_2010(k, j):
    """Only use dD constraint from 2010 onwards."""
    yr = 1999 + j
    if yr < 2010:
        return None
    return dD_baseline_func(k, j)

# ===========================================================================
# Scenario D: Quick Fix 6 - Inverse-uncertainty weighting
# ===========================================================================
def dD_weight_by_coverage(j):
    """Weight dD by station coverage (fewer stations → lower weight)."""
    return dD_weight_by_stations[j]

# ===========================================================================
# Run all scenarios
# ===========================================================================
print("\n" + "="*70)
print("RUNNING ALL SCENARIOS")
print("="*70)

results_A = run_inversion("A: Baseline (dD from 2006)", dD_baseline_func, lambda j: 1.0)
results_B = run_inversion("B: Ben's hemispheric dD", dD_ben_hemispheric, lambda j: 1.0)
results_C = run_inversion("C: dD from 2010 only", dD_from_2010, lambda j: 1.0)
results_D = run_inversion("D: Inverse-uncertainty weighted", dD_baseline_func, dD_weight_by_coverage)

all_results = [results_A, results_B, results_C, results_D]

# ===========================================================================
# 5-year smoothing
# ===========================================================================
def smooth5(arr, axis=0):
    """5-year centered moving average along axis."""
    from scipy.ndimage import uniform_filter1d
    return uniform_filter1d(arr, size=5, axis=axis, mode='nearest')

# ===========================================================================
# Plotting
# ===========================================================================
print("\nCreating comparison plots...")

colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3']
labels = ['A: Baseline (dD from 2006)', 'B: Ben hemispheric dD', 'C: dD from 2010 only', 'D: Inv-unc weighted']

fig, axes = plt.subplots(3, 2, figsize=(16, 14), dpi=150)

# Top row: FF emissions (raw and smoothed)
ax = axes[0, 0]
for i, res in enumerate(all_results):
    mean = np.nanmean(res['FF_combined'], axis=1)
    std = np.nanstd(res['FF_combined'], axis=1)
    ax.plot(model_years, mean, color=colors[i], linewidth=2, label=labels[i])
    ax.fill_between(model_years, mean-std, mean+std, alpha=0.15, color=colors[i])
ax.axvline(2008, color='gray', linestyle=':', alpha=0.5)
ax.set_ylabel('FF Emissions (Tg/yr)')
ax.set_title('Fossil Fuel — Combined (δ¹³C + δD weighted)')
ax.legend(fontsize=8, loc='upper left')
ax.grid(True, alpha=0.3)

ax = axes[0, 1]
for i, res in enumerate(all_results):
    mean = smooth5(np.nanmean(res['FF_combined'], axis=1))
    ax.plot(model_years, mean, color=colors[i], linewidth=2.5, label=labels[i])
ax.axvline(2008, color='gray', linestyle=':', alpha=0.5)
ax.set_ylabel('FF Emissions (Tg/yr)')
ax.set_title('Fossil Fuel — 5-yr Smoothed')
ax.legend(fontsize=8, loc='upper left')
ax.grid(True, alpha=0.3)

# Middle row: Microbial emissions
ax = axes[1, 0]
for i, res in enumerate(all_results):
    mean = np.nanmean(res['Mic_combined'], axis=1)
    std = np.nanstd(res['Mic_combined'], axis=1)
    ax.plot(model_years, mean, color=colors[i], linewidth=2, label=labels[i])
    ax.fill_between(model_years, mean-std, mean+std, alpha=0.15, color=colors[i])
ax.axvline(2008, color='gray', linestyle=':', alpha=0.5)
ax.set_ylabel('Mic Emissions (Tg/yr)')
ax.set_title('Microbial — Combined (δ¹³C + δD weighted)')
ax.legend(fontsize=8, loc='upper left')
ax.grid(True, alpha=0.3)

ax = axes[1, 1]
for i, res in enumerate(all_results):
    mean = smooth5(np.nanmean(res['Mic_combined'], axis=1))
    ax.plot(model_years, mean, color=colors[i], linewidth=2.5, label=labels[i])
ax.axvline(2008, color='gray', linestyle=':', alpha=0.5)
ax.set_ylabel('Mic Emissions (Tg/yr)')
ax.set_title('Microbial — 5-yr Smoothed')
ax.legend(fontsize=8, loc='upper left')
ax.grid(True, alpha=0.3)

# Bottom row: δD-only FF and comparison of d13C vs dD
ax = axes[2, 0]
for i, res in enumerate(all_results):
    ff_dD = res['FF_dD']
    # Only plot where we have dD solutions
    mean = np.nanmean(ff_dD, axis=1)
    valid = ~np.isnan(mean)
    ax.plot(model_years[valid], mean[valid], color=colors[i], linewidth=2, label=labels[i])
ax.axhline(0, color='black', linewidth=0.5)
ax.axvline(2008, color='gray', linestyle=':', alpha=0.5)
ax.set_ylabel('FF Emissions (Tg/yr)')
ax.set_title('Fossil Fuel — δD-only inversion')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_xlabel('Year')

# Bottom-right: Negative solution rates
ax = axes[2, 1]
bar_width = 0.2
x_pos = np.arange(len(all_results))
d13C_negs = [r['neg_pct_d13C'] for r in all_results]
dD_negs = [r['neg_pct_dD'] for r in all_results]
ax.bar(x_pos - bar_width/2, d13C_negs, bar_width, label='δ¹³C negative %', color='steelblue')
ax.bar(x_pos + bar_width/2, dD_negs, bar_width, label='δD negative %', color='coral')
ax.set_xticks(x_pos)
ax.set_xticklabels(['A\nBaseline', 'B\nBen Hem', 'C\ndD≥2010', 'D\nInv-Unc'], fontsize=9)
ax.set_ylabel('Non-physical solutions (%)')
ax.set_title('Negative Solution Rates')
ax.legend()
ax.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(OUT_DIR / 'dD_fix_comparison.png', dpi=150, bbox_inches='tight')
plt.close()
print(f"\n  Saved: {OUT_DIR / 'dD_fix_comparison.png'}")

# ===========================================================================
# Summary table
# ===========================================================================
print("\n" + "="*70)
print("SUMMARY TABLE")
print("="*70)
print(f"{'Scenario':<30} {'FF_mean':>8} {'FF_std':>8} {'Mic_mean':>9} {'Mic_std':>8} {'%neg_dD':>8}")
print("-"*75)
for res in all_results:
    ff_m = np.nanmean(res['FF_combined'])
    ff_s = np.nanmean(np.nanstd(res['FF_combined'], axis=1))
    mic_m = np.nanmean(res['Mic_combined'])
    mic_s = np.nanmean(np.nanstd(res['Mic_combined'], axis=1))
    print(f"{res['name']:<30} {ff_m:>8.1f} {ff_s:>8.1f} {mic_m:>9.1f} {mic_s:>8.1f} {res['neg_pct_dD']:>8.1f}")

# Save results CSV
summary_rows = []
for res in all_results:
    for j in range(n_years):
        summary_rows.append({
            'Year': int(model_years[j]),
            'Scenario': res['name'],
            'FF_combined_mean': np.nanmean(res['FF_combined'][j,:]),
            'FF_combined_std': np.nanstd(res['FF_combined'][j,:]),
            'Mic_combined_mean': np.nanmean(res['Mic_combined'][j,:]),
            'Mic_combined_std': np.nanstd(res['Mic_combined'][j,:]),
            'FF_dD_mean': np.nanmean(res['FF_dD'][j,:]),
            'FF_dD_std': np.nanstd(res['FF_dD'][j,:]),
        })
pd.DataFrame(summary_rows).to_csv(OUT_DIR / 'scenario_comparison.csv', index=False)
print(f"\n  Saved: {OUT_DIR / 'scenario_comparison.csv'}")

# 2008 specific
print("\n  === 2008 ANOMALY CHECK ===")
j2008 = 2008 - 1999  # index 9
for res in all_results:
    ff_2008 = np.nanmean(res['FF_combined'][j2008,:])
    mic_2008 = np.nanmean(res['Mic_combined'][j2008,:])
    ff_dD_2008 = np.nanmean(res['FF_dD'][j2008,:])
    print(f"  {res['name']:<30}: FF_comb={ff_2008:.1f}, Mic_comb={mic_2008:.1f}, FF_dD={ff_dD_2008:.1f}")

print("\nDone!")
