#!/usr/bin/env python3
"""
v4.0 — Mic vs Non-Mic Box Model
================================

Design philosophy: use δD PRIMARILY for microbial vs non-microbial separation,
and δ¹³C for FF vs BB partitioning within the non-microbial pool.

Architecture:
  Step 1: 1-box forward model → derive total source (S_tot) from CH4 observations
  Step 2: Dual-isotope (δ¹³C + δD) → derive source δ¹³C_src and δD_src  
  Step 3: δD-based 2-category split: Mic vs NonMic
         δD_src = f_mic * δD_mic + (1-f_mic) * δD_nonmic
         → f_mic = (δD_src - δD_nonmic) / (δD_mic - δD_nonmic)
  Step 4: δ¹³C-based sub-split of NonMic into FF + BB
         δ¹³C_nonmic = f_ff_in_nm * δ¹³C_FF + (1-f_ff_in_nm) * δ¹³C_BB
         → Need δ¹³C_nonmic from mass balance:
           δ¹³C_src = f_mic * δ¹³C_mic + (1-f_mic) * δ¹³C_nonmic

Key improvement over Ben 2025:
  - Temperature-dependent D/H KIE (not fixed at lab T)
  - MC uncertainty on effective tropospheric temperature
  - δD used only for Mic/NonMic split (robust: ~100-130‰ gap)
  - FF/BB split uses δ¹³C only (avoids δD's poor FF-BB discrimination)

Author: Generated for methane isotope project
Date: 2026-05-06
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ========================================================================
# CONFIGURATION
# ========================================================================
N_MC = 5000            # Monte Carlo iterations
OUTPUT_DIR = 'Output_dD_comparison'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ========================================================================
# ATMOSPHERIC OBSERVATIONS (NOAA global means)
# ========================================================================
ch4_ppb = {
    1999: 1772.33, 2000: 1773.33, 2001: 1771.22, 2002: 1772.66,
    2003: 1777.33, 2004: 1777.05, 2005: 1774.16, 2006: 1774.96,
    2007: 1781.38, 2008: 1787.01, 2009: 1793.53, 2010: 1798.93,
    2011: 1803.14, 2012: 1808.12, 2013: 1813.41, 2014: 1822.57,
    2015: 1834.26, 2016: 1843.12, 2017: 1849.58, 2018: 1857.33,
    2019: 1866.58, 2020: 1878.93, 2021: 1895.28, 2022: 1911.82,
}

# δ¹³C atmospheric (‰ VPDB) — NOAA INSTAAR global means
d13c_atm = {
    1998: -47.244, 1999: -47.119, 2000: -47.132, 2001: -47.095,
    2002: -47.097, 2003: -47.070, 2004: -47.053, 2005: -47.092,
    2006: -47.102, 2007: -47.089, 2008: -47.061, 2009: -47.135,
    2010: -47.175, 2011: -47.222, 2012: -47.238, 2013: -47.262,
    2014: -47.268, 2015: -47.302, 2016: -47.297, 2017: -47.354,
    2018: -47.403, 2019: -47.485, 2020: -47.478, 2021: -47.578,
    2022: -47.674,
}

# δD atmospheric (‰ VSMOW) — from Ben's global mean (glob_ann_dD.xlsx)
dD_atm = {
    2005: -76.72, 2006: -74.78, 2007: -75.36, 2008: -75.50,
    2009: -75.52, 2010: -76.07, 2011: -76.61, 2012: -77.16,
    2013: -77.71, 2014: -78.25, 2015: -77.38, 2016: -76.83,
    2017: -77.83, 2018: -77.98, 2019: -78.66, 2020: -78.80,
    2021: -80.15, 2022: -81.46,
}

# ========================================================================
# SOURCE SIGNATURES (with uncertainties for MC)
# ========================================================================
# δ¹³C source signatures (time-varying from Ben's database)
# Loading from our processed data
try:
    bb_d13c_df = pd.read_csv('TwoIsotopeBoxModel/BB_d13C_timeseries.csv')
    ff_d13c_df = pd.read_csv('TwoIsotopeBoxModel/FF_d13C_timeseries.csv')
    HAS_TIMESERIES = True
except:
    HAS_TIMESERIES = False
    print("WARNING: Time-series source signatures not found, using fixed values")

# Fixed source signatures (central values and uncertainties)
SIG = {
    # δ¹³C (‰ VPDB)
    'FF_d13C':  {'mean': -44.2, 'std': 0.5},    # Fossil fuel
    'BB_d13C':  {'mean': -23.0, 'std': 3.0},    # Biomass burning (C3/C4 mix)
    'Mic_d13C': {'mean': -61.7, 'std': 1.5},    # Microbial (weighted: livestock+waste+rice+wetlands)
    # δD (‰ VSMOW) — KEY: large separations for Mic vs NonMic
    'Mic_dD':   {'mean': -305, 'std': 10},       # Microbial (livestock -305, waste -312, rice -323, wetland -322)
    'NonMic_dD': {'mean': -190, 'std': 15},      # Non-microbial: FF(-175→-186) + BB(-169→-217) weighted
    'FF_dD':    {'mean': -183, 'std': 8},        # Fossil fuel (gas+oil: -175, coal: -175, thermogenic: -186)
    'BB_dD':    {'mean': -210, 'std': 25},       # Biomass burning (huge range: -169 to -217)
}

# ========================================================================
# SINK PARAMETERS
# ========================================================================
# Saueressig et al. (2001) temperature-dependent formulas
# k(CH4)/k(CH3D) = A * exp(B/T)   [A=1.097, B=49 K]
# k(CH4)/k(CH3D)_Cl = 1.278 * exp(53.31/T)

SINK = {
    'OH_frac':    {'mean': 0.835, 'std': 0.02},
    'Cl_frac':    {'mean': 0.035, 'std': 0.01},
    'Strat_frac': {'mean': 0.070, 'std': 0.01},
    'Soil_frac':  {'mean': 0.060, 'std': 0.01},
    # KIE parameters
    'OH_KIE_13C':   {'mean': 1.0054, 'std': 0.0008},  # Cantrell vs Saueressig
    'Cl_KIE_13C':   {'mean': 1.066, 'std': 0.005},
    'Strat_KIE_13C': {'mean': 1.003, 'std': 0.001},
    'Soil_KIE_13C':  {'mean': 1.0201, 'std': 0.002},
    # D/H KIE — temperature-dependent!
    'T_eff_OH':   {'mean': 272, 'std': 10},  # Effective tropospheric T for OH sink (K)
    'OH_A_D':     1.097,   # Pre-exponential (fixed, from Saueressig)
    'OH_B_D':     49.0,    # Exponent coefficient (K, fixed)
    'Cl_A_D':     1.278,
    'Cl_B_D':     53.31,
    'Strat_KIE_D': {'mean': 1.179, 'std': 0.02},  # Dyonisius 2020
    'Soil_KIE_D':  {'mean': 1.083, 'std': 0.01},
    # Lifetime
    'tau': {'mean': 9.1, 'std': 0.3},
}

ppb_to_Tg = 2.75  # conversion factor

# ========================================================================
# HELPER FUNCTIONS
# ========================================================================
def sample_param(param_dict, n=1):
    """Sample from normal distribution defined by {'mean', 'std'}"""
    return np.random.normal(param_dict['mean'], param_dict['std'], n)

def compute_dD_KIE(T_eff, A_OH, B_OH, A_Cl, B_Cl, strat_kie, soil_kie,
                    f_oh, f_cl, f_strat, f_soil):
    """Compute effective D/H KIE using temperature-dependent Saueressig formulas"""
    OH_KIE_D = A_OH * np.exp(B_OH / T_eff)
    # Cl at tropospheric T (approximate — Cl mostly tropospheric in recent work)
    Cl_KIE_D = A_Cl * np.exp(B_Cl / T_eff)
    net_KIE = OH_KIE_D * f_oh + Cl_KIE_D * f_cl + strat_kie * f_strat + soil_kie * f_soil
    return net_KIE, OH_KIE_D, Cl_KIE_D

def derive_source_isotope(ch4, iso_atm, tau, kie, ppb_to_Tg):
    """Derive source isotopic composition and total source from 1-box model.
    
    kie = k(light)/k(heavy) > 1  (the kinetic isotope effect / fractionation factor)
    For ¹³C: kie = k(¹²CH₄)/k(¹³CH₄) ≈ 1.005
    For D/H: kie = k(CH₄)/k(CH₃D) ≈ 1.30
    
    The heavy isotopologue has a LONGER effective lifetime = kie * tau.
    
    Returns dict: year → {S_tot (Tg/yr), d_src (‰)}
    """
    years = sorted(set(ch4.keys()) & set(iso_atm.keys()))
    results = {}
    for yr in years:
        if yr - 1 not in ch4 or yr + 1 not in ch4:
            continue
        if yr - 1 not in iso_atm or yr + 1 not in iso_atm:
            continue
        C = ch4[yr]
        d = iso_atm[yr]
        dCdt = (ch4[yr + 1] - ch4[yr - 1]) / 2.0
        dddt = (iso_atm[yr + 1] - iso_atm[yr - 1]) / 2.0
        S_ppb = dCdt + C / tau
        S_Tg = S_ppb * ppb_to_Tg
        d_prime = d + 1000.0
        # α = kie > 1: heavy isotopologue removed slower
        # d'_src = d' + [C·dδ'/dt - C·δ'/(τ)·(1 - 1/kie)] / S_ppb
        d_src_prime = d_prime + (C * dddt - C * d_prime / tau * (1 - 1.0 / kie)) / S_ppb
        d_src = d_src_prime - 1000.0
        results[yr] = {'S_tot': S_Tg, 'd_src': d_src, 'dCdt': dCdt}
    return results


# ========================================================================
# MONTE CARLO SIMULATION
# ========================================================================
print("="*80)
print("v4.0 — Mic vs Non-Mic Dual-Isotope Box Model")
print("="*80)
print(f"\nRunning {N_MC} Monte Carlo iterations...")

# Years where we have ALL three observables
common_years = sorted(set(ch4_ppb.keys()) & set(d13c_atm.keys()) & set(dD_atm.keys()))
# Need ±1 year for derivatives
valid_years = [y for y in common_years if y-1 in ch4_ppb and y+1 in ch4_ppb
               and y-1 in d13c_atm and y+1 in d13c_atm
               and y-1 in dD_atm and y+1 in dD_atm]

print(f"Valid years (triple-isotope): {valid_years[0]}–{valid_years[-1]} ({len(valid_years)} years)")

# δ¹³C-only years (longer record)
d13c_years = sorted(set(ch4_ppb.keys()) & set(d13c_atm.keys()))
d13c_valid = [y for y in d13c_years if y-1 in ch4_ppb and y+1 in ch4_ppb
              and y-1 in d13c_atm and y+1 in d13c_atm]
print(f"Valid years (δ¹³C only): {d13c_valid[0]}–{d13c_valid[-1]} ({len(d13c_valid)} years)")

# Storage arrays
results_mc = {yr: {'S_tot': [], 'f_mic': [], 'f_ff': [], 'f_bb': [],
                    'Mic': [], 'FF': [], 'BB': [], 'NonMic': [],
                    'd13C_src': [], 'dD_src': [],
                    'KIE_D': [], 'T_eff': [], 'physical': []}
              for yr in valid_years}

# Extended δ¹³C-only results
results_d13c_only = {yr: {'S_tot': [], 'd13C_src': []} for yr in d13c_valid}

n_nonphysical = 0
n_total = 0

for mc in range(N_MC):
    # --- Sample parameters ---
    tau = max(7.0, sample_param(SINK['tau'])[0])
    
    # Sink fractions (normalize to 1)
    f_oh = max(0.7, sample_param(SINK['OH_frac'])[0])
    f_cl = max(0.005, sample_param(SINK['Cl_frac'])[0])
    f_strat = max(0.03, sample_param(SINK['Strat_frac'])[0])
    f_soil = max(0.02, sample_param(SINK['Soil_frac'])[0])
    f_total = f_oh + f_cl + f_strat + f_soil
    f_oh /= f_total; f_cl /= f_total; f_strat /= f_total; f_soil /= f_total
    
    # ¹³C KIE (weighted) — k(¹²CH₄)/k(¹³CH₄) > 1
    kie_13c = (sample_param(SINK['OH_KIE_13C'])[0] * f_oh +
               sample_param(SINK['Cl_KIE_13C'])[0] * f_cl +
               sample_param(SINK['Strat_KIE_13C'])[0] * f_strat +
               sample_param(SINK['Soil_KIE_13C'])[0] * f_soil)
    
    # D/H KIE — TEMPERATURE DEPENDENT: k(CH₄)/k(CH₃D) > 1
    T_eff = sample_param(SINK['T_eff_OH'])[0]
    strat_kie_d = sample_param(SINK['Strat_KIE_D'])[0]
    soil_kie_d = sample_param(SINK['Soil_KIE_D'])[0]
    net_kie_d, oh_kie_d, cl_kie_d = compute_dD_KIE(
        T_eff, SINK['OH_A_D'], SINK['OH_B_D'],
        SINK['Cl_A_D'], SINK['Cl_B_D'],
        strat_kie_d, soil_kie_d,
        f_oh, f_cl, f_strat, f_soil)
    
    # Source signatures
    d13c_ff = sample_param(SIG['FF_d13C'])[0]
    d13c_bb = sample_param(SIG['BB_d13C'])[0]
    d13c_mic = sample_param(SIG['Mic_d13C'])[0]
    dD_mic = sample_param(SIG['Mic_dD'])[0]
    dD_nonmic = sample_param(SIG['NonMic_dD'])[0]
    
    # --- Derive source isotopic composition ---
    res_13c = derive_source_isotope(ch4_ppb, d13c_atm, tau, kie_13c, ppb_to_Tg)
    res_dD = derive_source_isotope(ch4_ppb, dD_atm, tau, net_kie_d, ppb_to_Tg)
    
    # Store δ¹³C-only results for extended period
    for yr in d13c_valid:
        if yr in res_13c:
            results_d13c_only[yr]['S_tot'].append(res_13c[yr]['S_tot'])
            results_d13c_only[yr]['d13C_src'].append(res_13c[yr]['d_src'])
    
    # --- Step 3: δD-based Mic vs NonMic split ---
    for yr in valid_years:
        if yr not in res_13c or yr not in res_dD:
            continue
        
        n_total += 1
        S_tot = res_13c[yr]['S_tot']  # total source from CH4 balance
        d13c_src = res_13c[yr]['d_src']
        dD_src = res_dD[yr]['d_src']
        
        # Microbial fraction from δD
        denom = dD_mic - dD_nonmic
        if abs(denom) < 5:  # guard against near-zero denominator
            n_nonphysical += 1
            results_mc[yr]['physical'].append(False)
            continue
        
        f_mic = (dD_src - dD_nonmic) / denom
        
        # Check physicality
        if f_mic < 0 or f_mic > 1:
            n_nonphysical += 1
            results_mc[yr]['physical'].append(False)
            # Still record but flag
            f_mic_clipped = np.clip(f_mic, 0, 1)
        else:
            results_mc[yr]['physical'].append(True)
            f_mic_clipped = f_mic
        
        # --- Step 4: δ¹³C-based FF vs BB within NonMic ---
        f_nonmic = 1 - f_mic_clipped
        
        # Derive δ¹³C of non-microbial pool
        if f_nonmic > 0.01:
            d13c_nonmic = (d13c_src - f_mic_clipped * d13c_mic) / f_nonmic
        else:
            d13c_nonmic = d13c_ff  # if almost all microbial, can't determine
        
        # FF fraction within non-microbial
        denom_13c = d13c_ff - d13c_bb
        if abs(denom_13c) < 1:
            f_ff_in_nm = 0.5
        else:
            f_ff_in_nm = (d13c_nonmic - d13c_bb) / denom_13c
            f_ff_in_nm = np.clip(f_ff_in_nm, 0, 1)
        
        # Absolute fluxes
        Mic = S_tot * f_mic_clipped
        NonMic = S_tot * f_nonmic
        FF = NonMic * f_ff_in_nm
        BB = NonMic * (1 - f_ff_in_nm)
        
        # Store
        results_mc[yr]['S_tot'].append(S_tot)
        results_mc[yr]['f_mic'].append(f_mic)  # raw (possibly outside 0-1)
        results_mc[yr]['f_ff'].append(f_ff_in_nm)
        results_mc[yr]['Mic'].append(Mic)
        results_mc[yr]['FF'].append(FF)
        results_mc[yr]['BB'].append(BB)
        results_mc[yr]['NonMic'].append(NonMic)
        results_mc[yr]['d13C_src'].append(d13c_src)
        results_mc[yr]['dD_src'].append(dD_src)
        results_mc[yr]['KIE_D'].append(net_kie_d)
        results_mc[yr]['T_eff'].append(T_eff)

print(f"\nNon-physical solutions: {n_nonphysical}/{n_total} ({100*n_nonphysical/n_total:.1f}%)")

# ========================================================================
# RESULTS SUMMARY
# ========================================================================
print(f"\n{'='*80}")
print("ANNUAL RESULTS (mean ± std of physical solutions)")
print(f"{'='*80}")
print(f"{'Year':>4s} {'S_tot':>8s} {'Mic':>8s} {'FF':>8s} {'BB':>8s} {'f_mic':>7s} {'δD_src':>8s} {'KIE_D':>7s} {'%phys':>6s}")
print(f"{'':>4s} {'Tg/yr':>8s} {'Tg/yr':>8s} {'Tg/yr':>8s} {'Tg/yr':>8s} {'':>7s} {'‰':>8s} {'':>7s} {'':>6s}")
print("-" * 72)

summary_rows = []
for yr in valid_years:
    r = results_mc[yr]
    if len(r['Mic']) == 0:
        continue
    
    phys = np.array(r['physical'])
    n_phys = phys.sum()
    pct_phys = 100 * n_phys / len(phys) if len(phys) > 0 else 0
    
    # Use all solutions (not just physical) for statistics
    Mic_arr = np.array(r['Mic'])
    FF_arr = np.array(r['FF'])
    BB_arr = np.array(r['BB'])
    S_arr = np.array(r['S_tot'])
    fmic_arr = np.array(r['f_mic'])
    dD_arr = np.array(r['dD_src'])
    kie_arr = np.array(r['KIE_D'])
    
    row = {
        'Year': yr,
        'S_tot_mean': np.mean(S_arr), 'S_tot_std': np.std(S_arr),
        'Mic_mean': np.mean(Mic_arr), 'Mic_std': np.std(Mic_arr),
        'FF_mean': np.mean(FF_arr), 'FF_std': np.std(FF_arr),
        'BB_mean': np.mean(BB_arr), 'BB_std': np.std(BB_arr),
        'f_mic_mean': np.mean(fmic_arr),
        'dD_src_mean': np.mean(dD_arr),
        'KIE_D_mean': np.mean(kie_arr),
        'pct_physical': pct_phys,
    }
    summary_rows.append(row)
    
    print(f"{yr:4d} {row['S_tot_mean']:8.1f} {row['Mic_mean']:8.1f} "
          f"{row['FF_mean']:8.1f} {row['BB_mean']:8.1f} "
          f"{row['f_mic_mean']:7.3f} {row['dD_src_mean']:8.1f} "
          f"{row['KIE_D_mean']:7.4f} {row['pct_physical']:5.1f}%")

df_summary = pd.DataFrame(summary_rows)

# ========================================================================
# CSV OUTPUT
# ========================================================================
csv_path = os.path.join(OUTPUT_DIR, 'v4.0_mic_vs_nonmic_results.csv')
df_summary.to_csv(csv_path, index=False, float_format='%.2f')
print(f"\nSaved: {csv_path}")

# Also save MC raw results for comparison
mc_raw_path = os.path.join(OUTPUT_DIR, 'v4.0_mc_raw.csv')
rows = []
for yr in valid_years:
    r = results_mc[yr]
    for i in range(len(r['Mic'])):
        rows.append({
            'Year': yr, 'S_tot': r['S_tot'][i],
            'Mic': r['Mic'][i], 'FF': r['FF'][i], 'BB': r['BB'][i],
            'f_mic': r['f_mic'][i], 'dD_src': r['dD_src'][i],
            'd13C_src': r['d13C_src'][i], 'KIE_D': r['KIE_D'][i],
            'physical': r['physical'][i],
        })
pd.DataFrame(rows).to_csv(mc_raw_path, index=False, float_format='%.4f')

# ========================================================================
# COMPARISON WITH EARLIER MODELS
# ========================================================================
print(f"\n{'='*80}")
print("COMPARISON WITH EARLIER MODEL VERSIONS")
print(f"{'='*80}")

# Load v3.2 (BB-fixed 2×2) results if available
try:
    v32 = pd.read_csv('Output_dD_comparison/v3.2_bb_fixed_results.csv')
    HAS_V32 = True
    print("\nLoaded v3.2 (BB-fixed 2×2) for comparison")
except:
    HAS_V32 = False
    print("\nv3.2 results not found — skipping comparison")

# Load v3.1 (3×3 optimized) results
try:
    v31 = pd.read_csv('Output_dD_comparison/v3.1_optimized_results.csv')
    HAS_V31 = True
    print("Loaded v3.1 (3×3 optimized) for comparison")
except:
    HAS_V31 = False

# Load Ben's original results
try:
    ben_csv = 'Output_dD_comparison/ben_original_results.csv'
    ben = pd.read_csv(ben_csv)
    HAS_BEN = True
    print("Loaded Ben's original results for comparison")
except:
    HAS_BEN = False

# ========================================================================
# MAIN COMPARISON FIGURE
# ========================================================================
fig, axes = plt.subplots(3, 2, figsize=(16, 18))
fig.suptitle('v4.0 Mic vs Non-Mic Model — Comparison with Earlier Versions', 
             fontsize=14, fontweight='bold')

years = df_summary['Year'].values

# --- Panel A: Total Source ---
ax = axes[0, 0]
ax.fill_between(years, df_summary['S_tot_mean'] - df_summary['S_tot_std'],
                df_summary['S_tot_mean'] + df_summary['S_tot_std'],
                alpha=0.3, color='black', label='v4.0 ±1σ')
ax.plot(years, df_summary['S_tot_mean'], 'k-', lw=2, label='v4.0')
if HAS_V32:
    v32y = v32[v32['Year'].isin(years)]
    if 'S_tot_mean' in v32.columns:
        ax.plot(v32y['Year'], v32y['S_tot_mean'], 'b--', lw=1.5, label='v3.2 (2×2)')
if HAS_V31:
    v31y = v31[v31['Year'].isin(years)]
    if 'S_tot_mean' in v31.columns:
        ax.plot(v31y['Year'], v31y['S_tot_mean'], 'r:', lw=1.5, label='v3.1 (3×3)')
ax.set_ylabel('Total Source (Tg/yr)')
ax.set_title('A) Total CH₄ Source')
ax.legend()
ax.grid(True, alpha=0.3)

# --- Panel B: Microbial Source ---
ax = axes[0, 1]
ax.fill_between(years, df_summary['Mic_mean'] - df_summary['Mic_std'],
                df_summary['Mic_mean'] + df_summary['Mic_std'],
                alpha=0.3, color='green')
ax.plot(years, df_summary['Mic_mean'], 'g-', lw=2, label='v4.0 Mic')
if HAS_V32 and 'Mic_mean' in v32.columns:
    v32y = v32[v32['Year'].isin(years)]
    ax.plot(v32y['Year'], v32y['Mic_mean'], 'g--', lw=1.5, label='v3.2 Mic')
ax.set_ylabel('Microbial Source (Tg/yr)')
ax.set_title('B) Microbial Emissions (δD-constrained)')
ax.legend()
ax.grid(True, alpha=0.3)

# --- Panel C: Fossil Fuel ---
ax = axes[1, 0]
ax.fill_between(years, df_summary['FF_mean'] - df_summary['FF_std'],
                df_summary['FF_mean'] + df_summary['FF_std'],
                alpha=0.3, color='red')
ax.plot(years, df_summary['FF_mean'], 'r-', lw=2, label='v4.0 FF')
if HAS_V32 and 'FF_mean' in v32.columns:
    v32y = v32[v32['Year'].isin(years)]
    ax.plot(v32y['Year'], v32y['FF_mean'], 'r--', lw=1.5, label='v3.2 FF')
ax.set_ylabel('FF Source (Tg/yr)')
ax.set_title('C) Fossil Fuel (δ¹³C sub-partitioned)')
ax.legend()
ax.grid(True, alpha=0.3)

# --- Panel D: Biomass Burning ---
ax = axes[1, 1]
ax.fill_between(years, df_summary['BB_mean'] - df_summary['BB_std'],
                df_summary['BB_mean'] + df_summary['BB_std'],
                alpha=0.3, color='orange')
ax.plot(years, df_summary['BB_mean'], '-', color='orange', lw=2, label='v4.0 BB')
if HAS_V32 and 'BB_mean' in v32.columns:
    v32y = v32[v32['Year'].isin(years)]
    ax.plot(v32y['Year'], v32y['BB_mean'], '--', color='orange', lw=1.5, label='v3.2 BB')
ax.set_ylabel('BB Source (Tg/yr)')
ax.set_title('D) Biomass Burning (δ¹³C residual)')
ax.legend()
ax.grid(True, alpha=0.3)

# --- Panel E: Microbial Fraction ---
ax = axes[2, 0]
ax.fill_between(years, 
                [np.percentile(results_mc[yr]['f_mic'], 16) for yr in years],
                [np.percentile(results_mc[yr]['f_mic'], 84) for yr in years],
                alpha=0.3, color='green')
ax.plot(years, df_summary['f_mic_mean'], 'g-', lw=2)
ax.axhline(y=0.73, color='gray', ls=':', label='IPCC AR6 (73%)')
ax.axhline(y=1.0, color='red', ls=':', alpha=0.3)
ax.axhline(y=0.0, color='red', ls=':', alpha=0.3)
ax.set_ylabel('Microbial Fraction')
ax.set_xlabel('Year')
ax.set_title('E) Microbial Fraction (from δD)')
ax.legend()
ax.grid(True, alpha=0.3)

# --- Panel F: Source δD and KIE ---
ax = axes[2, 1]
ax2 = ax.twinx()
dD_means = [np.mean(results_mc[yr]['dD_src']) for yr in years]
dD_stds = [np.std(results_mc[yr]['dD_src']) for yr in years]
ax.fill_between(years, np.array(dD_means) - np.array(dD_stds),
                np.array(dD_means) + np.array(dD_stds),
                alpha=0.2, color='blue')
l1, = ax.plot(years, dD_means, 'b-', lw=2, label='Source δD')
ax.set_ylabel('Source δD (‰)', color='blue')

kie_means = [np.mean(results_mc[yr]['KIE_D']) for yr in years]
l2, = ax2.plot(years, kie_means, 'r--', lw=1.5, label='Net KIE_D')
ax2.set_ylabel('Net KIE_D', color='red')

ax.legend(handles=[l1, l2], loc='lower left')
ax.set_xlabel('Year')
ax.set_title('F) Derived Source δD & D/H KIE')
ax.grid(True, alpha=0.3)

plt.tight_layout()
figpath = os.path.join(OUTPUT_DIR, 'v4.0_mic_vs_nonmic.png')
plt.savefig(figpath, dpi=150, bbox_inches='tight')
print(f"\nSaved: {figpath}")
plt.close()

# ========================================================================
# SENSITIVITY ANALYSIS: KIE temperature
# ========================================================================
print(f"\n{'='*80}")
print("SENSITIVITY ANALYSIS: Effect of T_eff on results")
print(f"{'='*80}")

fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle('Sensitivity of v4.0 to Effective OH Temperature', fontsize=13, fontweight='bold')

T_tests = [260, 270, 280, 290, 296]
colors = plt.cm.coolwarm(np.linspace(0, 1, len(T_tests)))

for i_t, T_test in enumerate(T_tests):
    oh_kie = 1.097 * np.exp(49 / T_test)
    cl_kie = 1.278 * np.exp(53.31 / T_test)
    net_kie_test = oh_kie * 0.835 + cl_kie * 0.035 + 1.179 * 0.07 + 1.083 * 0.06
    tau_test = 9.1
    
    kie_13c_test = 1.0054 * 0.835 + 1.066 * 0.035 + 1.003 * 0.07 + 1.0201 * 0.06
    res_dD_test = derive_source_isotope(ch4_ppb, dD_atm, tau_test, net_kie_test, ppb_to_Tg)
    res_13c_test = derive_source_isotope(ch4_ppb, d13c_atm, tau_test, kie_13c_test, ppb_to_Tg)
    
    test_years = sorted(set(res_dD_test.keys()) & set(res_13c_test.keys()))
    dD_src_list = [res_dD_test[y]['d_src'] for y in test_years]
    d13c_src_list = [res_13c_test[y]['d_src'] for y in test_years]
    
    # Mic fraction
    dD_mic_mean = -305
    dD_nonmic_mean = -190
    f_mic_list = [(d - dD_nonmic_mean) / (dD_mic_mean - dD_nonmic_mean) for d in dD_src_list]
    S_tot_list = [res_13c_test[y]['S_tot'] for y in test_years]
    mic_list = [f * s for f, s in zip(f_mic_list, S_tot_list)]
    
    label = f'T={T_test}K, KIE_D={net_kie_test:.3f}'
    axes[0].plot(test_years, dD_src_list, '-', color=colors[i_t], lw=2, label=label)
    axes[1].plot(test_years, f_mic_list, '-', color=colors[i_t], lw=2, label=f'T={T_test}K')
    axes[2].plot(test_years, mic_list, '-', color=colors[i_t], lw=2, label=f'T={T_test}K')

axes[0].set_ylabel('Source δD (‰)')
axes[0].set_title('Source δD')
axes[0].legend(fontsize=7)
axes[0].grid(True, alpha=0.3)

axes[1].set_ylabel('Microbial Fraction')
axes[1].set_title('Microbial Fraction')
axes[1].axhline(y=1, color='red', ls=':', alpha=0.3)
axes[1].legend(fontsize=8)
axes[1].grid(True, alpha=0.3)

axes[2].set_ylabel('Microbial Source (Tg/yr)')
axes[2].set_title('Microbial Emissions')
axes[2].legend(fontsize=8)
axes[2].grid(True, alpha=0.3)

for ax in axes:
    ax.set_xlabel('Year')

plt.tight_layout()
sens_path = os.path.join(OUTPUT_DIR, 'v4.0_T_sensitivity.png')
plt.savefig(sens_path, dpi=150, bbox_inches='tight')
print(f"Saved: {sens_path}")
plt.close()

# ========================================================================
# DIAGNOSTICS
# ========================================================================
print(f"\n{'='*80}")
print("MODEL DIAGNOSTICS")
print(f"{'='*80}")

# Mean KIE across MC
all_kie = []
all_T = []
for yr in valid_years:
    all_kie.extend(results_mc[yr]['KIE_D'])
    all_T.extend(results_mc[yr]['T_eff'])

print(f"Mean effective T_OH: {np.mean(all_T):.1f} ± {np.std(all_T):.1f} K")
print(f"Mean net KIE_D: {np.mean(all_kie):.4f} ± {np.std(all_kie):.4f}")
print(f"  → α_D = {1/np.mean(all_kie):.4f}")
print(f"  → ε_D = {(1/np.mean(all_kie) - 1)*1000:.1f} ‰")
print(f"  cf. Ben's KIE_D = 1.281, Rice's KIE_D = 1.264")

# Physical solution rate per year
print(f"\nPhysical solution rates:")
for yr in valid_years:
    phys = np.array(results_mc[yr]['physical'])
    print(f"  {yr}: {100*phys.sum()/len(phys):.1f}% physical ({phys.sum()}/{len(phys)})")

# Comparison of trends
if len(df_summary) > 5:
    # Linear trend in microbial
    from numpy.polynomial import polynomial as P
    x = df_summary['Year'].values
    y = df_summary['Mic_mean'].values
    coeffs = np.polyfit(x, y, 1)
    print(f"\nMicrobial emission trend: {coeffs[0]:.1f} Tg/yr per year")
    print(f"  2006: {np.polyval(coeffs, 2006):.0f} Tg/yr")
    print(f"  2021: {np.polyval(coeffs, 2021):.0f} Tg/yr")
    print(f"  Δ(2006-2021): {np.polyval(coeffs, 2021) - np.polyval(coeffs, 2006):.0f} Tg/yr")

print(f"\nDone! All outputs in {OUTPUT_DIR}/")
