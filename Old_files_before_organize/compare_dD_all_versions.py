#!/usr/bin/env python3
"""
Compare ALL box model versions with Ben's δD vs Our Improved δD
================================================================

Uses the ACTUAL derive_source_isotope from v4.0 and the actual model logic
from each version. Runs each model with:
  A) Ben's original δD
  B) Our improved δD (area-weighted, from improved_dD_pipeline.py)

For v2.0, uses actual NH/SH δD from both sources.

Output: Output_improved_dD/
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import os, sys
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

BASE_DIR = Path(__file__).parent
OUTPUT_DIR = BASE_DIR / 'Output_improved_dD'
OUTPUT_DIR.mkdir(exist_ok=True)

PPB_TO_TG = 2.75

# ========================================================================
# ATMOSPHERIC DATA
# ========================================================================
ch4_ppb = {
    1999: 1772.33, 2000: 1773.33, 2001: 1771.22, 2002: 1772.66,
    2003: 1777.33, 2004: 1777.05, 2005: 1774.16, 2006: 1774.96,
    2007: 1781.38, 2008: 1787.01, 2009: 1793.53, 2010: 1798.93,
    2011: 1803.14, 2012: 1808.12, 2013: 1813.41, 2014: 1822.57,
    2015: 1834.26, 2016: 1843.12, 2017: 1849.58, 2018: 1857.33,
    2019: 1866.58, 2020: 1878.93, 2021: 1895.28, 2022: 1911.82,
}

d13c_atm = {
    1998: -47.244, 1999: -47.119, 2000: -47.132, 2001: -47.095,
    2002: -47.097, 2003: -47.070, 2004: -47.053, 2005: -47.092,
    2006: -47.102, 2007: -47.089, 2008: -47.061, 2009: -47.135,
    2010: -47.175, 2011: -47.222, 2012: -47.238, 2013: -47.262,
    2014: -47.268, 2015: -47.302, 2016: -47.297, 2017: -47.354,
    2018: -47.403, 2019: -47.485, 2020: -47.478, 2021: -47.578,
    2022: -47.674,
}

# Ben's original δD (from glob_ann_dD.xlsx, with -0.5‰ adjustment as in models)
ben_dD_atm = {
    2005: -76.72, 2006: -74.78, 2007: -75.36, 2008: -75.50,
    2009: -75.52, 2010: -76.07, 2011: -76.61, 2012: -77.16,
    2013: -77.71, 2014: -78.25, 2015: -77.38, 2016: -76.83,
    2017: -77.83, 2018: -77.98, 2019: -78.66, 2020: -78.80,
    2021: -80.15, 2022: -81.46,
}

# Our improved δD (area-weighted)
our_dD_df = pd.read_csv(BASE_DIR / 'Output_dD_comparison' / 'improved_dD_global_mean.csv')
our_dD_atm = {int(row['Year']): row['dD_glob_area'] if 'dD_glob_area' in our_dD_df.columns 
               else row.iloc[3]  # fallback
               for _, row in our_dD_df.iterrows()}
our_dD_NH = {int(row['Year']): row['dD_NH'] for _, row in our_dD_df.iterrows()}
our_dD_SH = {int(row['Year']): row['dD_SH'] for _, row in our_dD_df.iterrows()}

# Ben's hemispheric δD
BEN_DD_DIR = BASE_DIR.parent / 'Ben-BoxModel' / 'Riddell-Young_2025_dD_GlobMean' / 'Riddell-Young_2025_dD_GlobMean'
hem_df = pd.read_csv(BEN_DD_DIR / 'output' / 'HemMean_dD_dei_UmezawaCal_noBUDS.csv')
hem_df['yr'] = hem_df['Year'].astype(int)
hem_ann = hem_df.groupby('yr').mean(numeric_only=True)
ben_dD_NH = {int(y): row['NH_smooth_mean'] for y, row in hem_ann.iterrows() if not np.isnan(row['NH_smooth_mean'])}
ben_dD_SH = {int(y): row['SH_smooth_mean'] for y, row in hem_ann.iterrows() if not np.isnan(row['SH_smooth_mean'])}

# ========================================================================
# CORE FUNCTION (same as v4.0)
# ========================================================================
def derive_source_isotope(ch4, iso_atm, tau, kie, ppb_to_Tg):
    """Derive source isotopic composition from 1-box mass balance."""
    years = sorted(set(ch4.keys()) & set(iso_atm.keys()))
    results = {}
    for yr in years:
        if yr-1 not in ch4 or yr+1 not in ch4: continue
        if yr-1 not in iso_atm or yr+1 not in iso_atm: continue
        C = ch4[yr]
        d = iso_atm[yr]
        dCdt = (ch4[yr+1] - ch4[yr-1]) / 2.0
        dddt = (iso_atm[yr+1] - iso_atm[yr-1]) / 2.0
        S_ppb = dCdt + C / tau
        S_Tg = S_ppb * ppb_to_Tg
        d_prime = d + 1000.0
        d_src_prime = d_prime + (C * dddt - C * d_prime / tau * (1 - 1.0/kie)) / S_ppb
        d_src = d_src_prime - 1000.0
        results[yr] = {'S_tot': S_Tg, 'd_src': d_src}
    return results

# ========================================================================
# SOURCE SIGNATURES
# ========================================================================
src_13c = {'FF': (-44.0, 0.5), 'BB': (-22.3, 3.7), 'Mic': (-62.0, 1.0)}
src_dD = {'FF': (-183, 8), 'BB': (-210, 25), 'Mic': (-305, 10), 'NonMic': (-190, 15)}

# ========================================================================
# KIE SAMPLING
# ========================================================================
def sample_kie_13c():
    oh = np.random.normal(1.0039, 0.0004)
    cl = np.random.normal(1.066, 0.002)
    return 0.84*oh + 0.03*cl + 0.05*1.018 + 0.08*1.013

def sample_kie_d(T_eff=None):
    if T_eff is None:
        T_eff = np.random.normal(272, 10)
    oh = 1.097 * np.exp(49.0/T_eff)
    cl = 1.278 * np.exp(53.31/T_eff)
    f = np.array([0.84, 0.03, 0.05, 0.08, 0.001])
    k = np.array([oh, cl, 1.083, 1.16, 1.06])
    return np.sum(f*k) / np.sum(f)

# ========================================================================
# MODEL RUNNERS
# ========================================================================
N_MC = 3000

def run_3x3(dD_dict, label, n_mc=N_MC):
    """v1.0/v3.1: 3×3 inversion (FF, BB, Mic) using δ13C + δD"""
    rows = []
    for k in range(n_mc):
        tau = np.random.normal(9.1, 0.9)
        if tau <= 0: continue
        kie_13c = sample_kie_13c()
        kie_d = sample_kie_d()
        
        ff13 = np.random.normal(*src_13c['FF'])
        bb13 = np.random.normal(*src_13c['BB'])
        mic13 = np.random.normal(*src_13c['Mic'])
        ffd = np.random.normal(*src_dD['FF'])
        bbd = np.random.normal(*src_dD['BB'])
        micd = np.random.normal(*src_dD['Mic'])
        
        r13 = derive_source_isotope(ch4_ppb, d13c_atm, tau, kie_13c, PPB_TO_TG)
        rd = derive_source_isotope(ch4_ppb, dD_dict, tau, kie_d, PPB_TO_TG)
        
        common = sorted(set(r13.keys()) & set(rd.keys()))
        for yr in common:
            S = r13[yr]['S_tot']
            d13_s = r13[yr]['d_src']
            dD_s = rd[yr]['d_src']
            
            A = np.array([[1,1,1], [ff13, bb13, mic13], [ffd, bbd, micd]])
            b = np.array([S, d13_s*S, dD_s*S])
            try:
                x = np.linalg.solve(A, b)
                rows.append({'year': yr, 'FF': x[0], 'BB': x[1], 'Mic': x[2],
                           'S_tot': S, 'd13C_src': d13_s, 'dD_src': dD_s,
                           'physical': all(xi >= 0 for xi in x)})
            except: pass
    return pd.DataFrame(rows)

# Load annual GFED4 BB from CarbonTracker (column 9 = "Pyro Prior")
_ct_df = pd.read_excel(BASE_DIR / 'rel' / 'data' / 'CarbonTracker_CH4.xlsx')
BB_GFED = {int(_ct_df.iloc[i, 0]): float(_ct_df.iloc[i, 9]) for i in range(len(_ct_df))}
print(f"  GFED4 BB loaded: {min(BB_GFED.keys())}–{max(BB_GFED.keys())}, "
      f"mean={np.mean(list(BB_GFED.values())):.1f} Tg/yr")

def run_bb_fixed(n_mc=N_MC):
    """v3.2: BB fixed from GFED4 (annual), δ13C only → FF + Mic"""
    rows = []
    for k in range(n_mc):
        tau = np.random.normal(9.1, 0.9)
        if tau <= 0: continue
        kie_13c = sample_kie_13c()
        
        ff13 = np.random.normal(*src_13c['FF'])
        bb13 = np.random.normal(*src_13c['BB'])
        mic13 = np.random.normal(*src_13c['Mic'])
        
        # Add uncertainty to BB: ±20% (CarbonTracker prior uncertainty)
        bb_scale = np.random.normal(1.0, 0.2)
        
        r13 = derive_source_isotope(ch4_ppb, d13c_atm, tau, kie_13c, PPB_TO_TG)
        for yr in sorted(r13.keys()):
            S = r13[yr]['S_tot']
            d13_s = r13[yr]['d_src']
            # Use annual GFED BB, fall back to mean if year not available
            BB = BB_GFED.get(yr, np.mean(list(BB_GFED.values()))) * bb_scale
            if BB < 0: BB = 0
            rem = S - BB
            if rem <= 0: continue
            d13_rem = (d13_s * S - BB * bb13) / rem
            if abs(ff13 - mic13) < 0.1: continue
            f_ff = (d13_rem - mic13) / (ff13 - mic13)
            FF = f_ff * rem
            Mic = (1 - f_ff) * rem
            rows.append({'year': yr, 'FF': FF, 'BB': BB, 'Mic': Mic,
                       'S_tot': S, 'd13C_src': d13_s,
                       'physical': FF >= 0 and Mic >= 0})
    return pd.DataFrame(rows)

def run_mic_nonmic(dD_dict, label, n_mc=N_MC):
    """v4.0: δD → Mic vs NonMic, then δ13C → FF + BB"""
    rows = []
    for k in range(n_mc):
        tau = np.random.normal(9.1, 0.9)
        if tau <= 0: continue
        kie_13c = sample_kie_13c()
        kie_d = sample_kie_d()
        
        micd = np.random.normal(*src_dD['Mic'])
        nonmicd = np.random.normal(*src_dD['NonMic'])
        ff13 = np.random.normal(*src_13c['FF'])
        bb13 = np.random.normal(*src_13c['BB'])
        mic13 = np.random.normal(*src_13c['Mic'])
        
        r13 = derive_source_isotope(ch4_ppb, d13c_atm, tau, kie_13c, PPB_TO_TG)
        rd = derive_source_isotope(ch4_ppb, dD_dict, tau, kie_d, PPB_TO_TG)
        
        common = sorted(set(r13.keys()) & set(rd.keys()))
        for yr in common:
            S = r13[yr]['S_tot']
            d13_s = r13[yr]['d_src']
            dD_s = rd[yr]['d_src']
            
            if abs(micd - nonmicd) < 1: continue
            f_mic = (dD_s - nonmicd) / (micd - nonmicd)
            Mic = f_mic * S
            NonMic = (1 - f_mic) * S
            
            if NonMic <= 0 or abs(ff13 - bb13) < 0.1: continue
            d13_nm = (d13_s * S - Mic * mic13) / NonMic
            f_ff = (d13_nm - bb13) / (ff13 - bb13)
            FF = f_ff * NonMic
            BB = (1 - f_ff) * NonMic
            
            rows.append({'year': yr, 'FF': FF, 'BB': BB, 'Mic': Mic,
                       'S_tot': S, 'd13C_src': d13_s, 'dD_src': dD_s,
                       'f_mic': f_mic,
                       'physical': FF >= 0 and BB >= 0 and Mic >= 0})
    return pd.DataFrame(rows)

def run_2hemi_3x3(dD_NH, dD_SH, label, n_mc=N_MC):
    """v2.0: uses NH/SH δD separately, averages global source"""
    rows = []
    for k in range(n_mc):
        tau = np.random.normal(9.1, 0.9)
        if tau <= 0: continue
        kie_13c = sample_kie_13c()
        kie_d = sample_kie_d()
        
        ff13 = np.random.normal(*src_13c['FF'])
        bb13 = np.random.normal(*src_13c['BB'])
        mic13 = np.random.normal(*src_13c['Mic'])
        ffd = np.random.normal(*src_dD['FF'])
        bbd = np.random.normal(*src_dD['BB'])
        micd = np.random.normal(*src_dD['Mic'])
        
        # Derive source for each hemisphere
        r13 = derive_source_isotope(ch4_ppb, d13c_atm, tau, kie_13c, PPB_TO_TG)
        rd_nh = derive_source_isotope(ch4_ppb, dD_NH, tau, kie_d, PPB_TO_TG)
        rd_sh = derive_source_isotope(ch4_ppb, dD_SH, tau, kie_d, PPB_TO_TG)
        
        common = sorted(set(r13.keys()) & set(rd_nh.keys()) & set(rd_sh.keys()))
        for yr in common:
            S = r13[yr]['S_tot']
            d13_s = r13[yr]['d_src']
            # Average the NH and SH source δD
            dD_s = (rd_nh[yr]['d_src'] + rd_sh[yr]['d_src']) / 2.0
            
            A = np.array([[1,1,1], [ff13, bb13, mic13], [ffd, bbd, micd]])
            b = np.array([S, d13_s*S, dD_s*S])
            try:
                x = np.linalg.solve(A, b)
                rows.append({'year': yr, 'FF': x[0], 'BB': x[1], 'Mic': x[2],
                           'S_tot': S, 'd13C_src': d13_s, 'dD_src': dD_s,
                           'physical': all(xi >= 0 for xi in x)})
            except: pass
    return pd.DataFrame(rows)

def run_3x3_dD_from2010(dD_dict, label, n_mc=N_MC):
    """v3.3: Use δD only from 2010+"""
    dD_2010 = {y: v for y, v in dD_dict.items() if y >= 2010}
    return run_3x3(dD_2010, label, n_mc)

# ========================================================================
# RUN ALL MODELS
# ========================================================================
print("="*80)
print("RUNNING ALL MODEL VERSIONS: Ben's δD vs Our Improved δD")
print("="*80, flush=True)

results = {}

configs = [
    ('v1.0 (3×3)', lambda dD: run_3x3(dD, 'v1.0'), True, False),
    ('v2.0 (2-hemi)', None, True, True),  # special
    ('v3.1 (3×3 opt)', lambda dD: run_3x3(dD, 'v3.1'), True, False),
    ('v3.2 (BB fix)', lambda dD: run_bb_fixed(), False, False),  # no δD
    ('v3.3 (dD≥2010)', lambda dD: run_3x3_dD_from2010(dD, 'v3.3'), True, False),
    ('v4.0 (Mic/NM)', lambda dD: run_mic_nonmic(dD, 'v4.0'), True, False),
]

for name, func, uses_dD, uses_hemi in configs:
    print(f"\n--- {name} ---", flush=True)
    
    if not uses_dD:
        print("  δ13C only — no δD difference")
        res = func(None)
        results[f'{name}_ben'] = res
        results[f'{name}_ours'] = res.copy()
    elif uses_hemi:
        print("  Running with Ben's NH/SH δD...", flush=True)
        res_ben = run_2hemi_3x3(ben_dD_NH, ben_dD_SH, 'v2.0_ben')
        print("  Running with Our NH/SH δD...", flush=True)
        res_ours = run_2hemi_3x3(our_dD_NH, our_dD_SH, 'v2.0_ours')
        results[f'{name}_ben'] = res_ben
        results[f'{name}_ours'] = res_ours
    else:
        print("  Running with Ben's δD...", flush=True)
        res_ben = func(ben_dD_atm)
        print("  Running with Our δD...", flush=True)
        res_ours = func(our_dD_atm)
        results[f'{name}_ben'] = res_ben
        results[f'{name}_ours'] = res_ours
    
    for tag in ['ben', 'ours']:
        r = results[f'{name}_{tag}']
        if len(r) > 0:
            pct = r['physical'].mean() * 100
            print(f"  {tag}: {pct:.1f}% physical ({len(r)} total)", flush=True)

# ========================================================================
# SUMMARY TABLE
# ========================================================================
print(f"\n{'='*80}")
print("SUMMARY: Mean Emissions 2010-2021 (physical solutions, Tg/yr)")
print(f"{'='*80}")
print(f"{'Model':<20s} {'δD':<6s} {'FF':>7s} {'BB':>7s} {'Mic':>7s} {'Tot':>7s} {'Phys%':>7s}")
print("-"*60)

summary_rows = []
for name, *_ in configs:
    for tag in ['ben', 'ours']:
        key = f'{name}_{tag}'
        df = results[key]
        if len(df) == 0: continue
        phys = df[df['physical']]
        sub = phys[(phys['year'] >= 2010) & (phys['year'] <= 2021)]
        if len(sub) == 0: continue
        pct = df['physical'].mean() * 100
        row = {'model': name, 'dD': tag,
               'FF': sub['FF'].mean(), 'BB': sub['BB'].mean(),
               'Mic': sub['Mic'].mean(), 'Tot': sub['S_tot'].mean(),
               'Phys_pct': pct}
        summary_rows.append(row)
        print(f"{name:<20s} {tag:<6s} {row['FF']:7.1f} {row['BB']:7.1f} "
              f"{row['Mic']:7.1f} {row['Tot']:7.1f} {pct:6.1f}%")

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(OUTPUT_DIR / 'summary_comparison.csv', index=False, float_format='%.2f')

# ========================================================================
# DETAILED ANNUAL RESULTS
# ========================================================================
detail_rows = []
for name, *_ in configs:
    for tag in ['ben', 'ours']:
        key = f'{name}_{tag}'
        df = results[key]
        if len(df) == 0: continue
        phys = df[df['physical']]
        if len(phys) == 0: continue
        ann = phys.groupby('year').agg(
            FF_mean=('FF', 'mean'), FF_std=('FF', 'std'),
            BB_mean=('BB', 'mean'), BB_std=('BB', 'std'),
            Mic_mean=('Mic', 'mean'), Mic_std=('Mic', 'std'),
            S_tot=('S_tot', 'mean'), n=('physical', 'size'),
        ).reset_index()
        for _, r in ann.iterrows():
            detail_rows.append({'model': name, 'dD': tag, 'year': int(r['year']),
                              'FF': r['FF_mean'], 'FF_std': r['FF_std'],
                              'BB': r['BB_mean'], 'BB_std': r['BB_std'],
                              'Mic': r['Mic_mean'], 'Mic_std': r['Mic_std'],
                              'S_tot': r['S_tot'], 'n_phys': int(r['n'])})

detail_df = pd.DataFrame(detail_rows)
detail_df.to_csv(OUTPUT_DIR / 'annual_comparison.csv', index=False, float_format='%.2f')

# ========================================================================
# PLOTS
# ========================================================================
print("\nGenerating plots...", flush=True)

model_names = [c[0] for c in configs]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']

# ---- PLOT 1: 6-panel source partitioning ----
fig, axes = plt.subplots(3, 2, figsize=(18, 20), sharex=True)
fig.suptitle('Source Partitioning: Ben\'s δD (solid) vs Our Improved δD (dashed)\nAll Model Versions',
             fontsize=16, fontweight='bold')

for idx, name in enumerate(model_names):
    ax = axes.flat[idx]
    
    for tag, ls, lw_mult in [('ben', '-', 1.0), ('ours', '--', 0.8)]:
        key = f'{name}_{tag}'
        df = results[key]
        if len(df) == 0: continue
        phys = df[df['physical']]
        if len(phys) == 0: continue
        
        ann = phys.groupby('year').agg(
            FF=('FF', 'mean'), BB=('BB', 'mean'), Mic=('Mic', 'mean'),
            FF_lo=('FF', lambda x: np.percentile(x, 16)),
            FF_hi=('FF', lambda x: np.percentile(x, 84)),
            BB_lo=('BB', lambda x: np.percentile(x, 16)),
            BB_hi=('BB', lambda x: np.percentile(x, 84)),
            Mic_lo=('Mic', lambda x: np.percentile(x, 16)),
            Mic_hi=('Mic', lambda x: np.percentile(x, 84)),
        ).reset_index()
        
        prefix = "Ben" if tag == 'ben' else "Ours"
        for cat, color in [('FF', '#1f77b4'), ('BB', '#ff7f0e'), ('Mic', '#2ca02c')]:
            ax.plot(ann['year'], ann[cat], ls=ls, color=color, lw=2*lw_mult,
                    label=f'{prefix} {cat}' if idx == 0 else '')
            if tag == 'ben':
                ax.fill_between(ann['year'], ann[f'{cat}_lo'], ann[f'{cat}_hi'],
                              alpha=0.1, color=color)
    
    ax.set_title(name, fontsize=12, fontweight='bold')
    ax.set_ylabel('Emissions (Tg/yr)')
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='k', ls=':', alpha=0.3)

axes[0, 0].legend(fontsize=7, ncol=2, loc='upper right')
for ax in axes[-1]:
    ax.set_xlabel('Year')
plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'all_versions_source_partitioning.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: all_versions_source_partitioning.png", flush=True)

# ---- PLOT 2: Δ(Ours - Ben) impact ----
fig, axes = plt.subplots(1, 3, figsize=(18, 6))
fig.suptitle('Impact of Improved δD: Δ(Our − Ben) Emissions by Source',
             fontsize=14, fontweight='bold')

for ci, cat in enumerate(['FF', 'BB', 'Mic']):
    ax = axes[ci]
    for idx, name in enumerate(model_names):
        key_b = f'{name}_ben'
        key_o = f'{name}_ours'
        db = results[key_b]
        do_ = results[key_o]
        if len(db) == 0 or len(do_) == 0: continue
        
        ab = db[db['physical']].groupby('year')[cat].mean()
        ao = do_[do_['physical']].groupby('year')[cat].mean()
        common = sorted(set(ab.index) & set(ao.index))
        if not common: continue
        
        delta = [ao[y] - ab[y] for y in common]
        ax.plot(common, delta, '-o', ms=3, lw=1.5, color=colors[idx], label=name)
    
    ax.set_title(f'{cat} Emissions', fontweight='bold')
    ax.set_ylabel('Δ Emissions (Tg/yr)')
    ax.set_xlabel('Year')
    ax.axhline(0, color='k', ls='--', alpha=0.5)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'delta_impact.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: delta_impact.png", flush=True)

# ---- PLOT 3: Physical solution rate ----
fig, ax = plt.subplots(figsize=(12, 6))
fig.suptitle('Physical Solution Rate: Ben\'s δD vs Our δD', fontsize=14)

x = np.arange(len(model_names))
w = 0.35
rates_b = []
rates_o = []
for name in model_names:
    rb = results.get(f'{name}_ben', pd.DataFrame())
    ro = results.get(f'{name}_ours', pd.DataFrame())
    rates_b.append(rb['physical'].mean()*100 if len(rb) > 0 else 0)
    rates_o.append(ro['physical'].mean()*100 if len(ro) > 0 else 0)

b1 = ax.bar(x - w/2, rates_b, w, label="Ben's δD", color='#1f77b4', alpha=0.8)
b2 = ax.bar(x + w/2, rates_o, w, label="Our δD", color='#ff7f0e', alpha=0.8)
ax.set_xticks(x)
ax.set_xticklabels(model_names, rotation=30, ha='right')
ax.set_ylabel('Physical Solutions (%)')
ax.legend()
ax.grid(True, alpha=0.3, axis='y')
for bars in [b1, b2]:
    for bar in bars:
        h = bar.get_height()
        if h > 0:
            ax.text(bar.get_x()+bar.get_width()/2, h+0.5, f'{h:.1f}%',
                    ha='center', va='bottom', fontsize=8)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'physical_rates.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: physical_rates.png", flush=True)

# ---- PLOT 4: Input δD comparison ----
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle('δD-CH₄ Input Data Comparison', fontsize=14, fontweight='bold')

ax = axes[0]
yrs = sorted(set(ben_dD_atm.keys()) & set(our_dD_atm.keys()))
ax.plot(yrs, [ben_dD_atm[y] for y in yrs], 'ko-', lw=2, ms=5, label="Ben's global")
ax.plot(yrs, [our_dD_atm[y] for y in yrs], 'r^-', lw=2, ms=5, label="Our global")
ax.set_ylabel('δD-CH₄ (‰)')
ax.set_title('Global Mean')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[1]
yrs = sorted(set(ben_dD_NH.keys()) & set(our_dD_NH.keys()))
ax.plot(yrs, [ben_dD_NH[y] for y in yrs], 'b-', lw=2, label="Ben NH")
ax.plot(yrs, [our_dD_NH[y] for y in yrs], 'b--', lw=2, label="Our NH")
ax.plot(yrs, [ben_dD_SH[y] for y in yrs], 'r-', lw=2, label="Ben SH")
ax.plot(yrs, [our_dD_SH[y] for y in yrs], 'r--', lw=2, label="Our SH")
ax.set_ylabel('δD-CH₄ (‰)')
ax.set_title('Hemispheric')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'dD_input_comparison.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: dD_input_comparison.png", flush=True)

# ---- PLOT 5: Best models detail (v3.2 + v4.0) ----
fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.suptitle('Best Models: v3.2 (BB fixed) vs v4.0 (Mic/NonMic)\nBen\'s δD (top) vs Our δD (bottom)',
             fontsize=14, fontweight='bold')

for col, name in enumerate(['v3.2 (BB fix)', 'v4.0 (Mic/NM)']):
    for row, tag in enumerate(['ben', 'ours']):
        ax = axes[row, col]
        key = f'{name}_{tag}'
        df = results[key]
        phys = df[df['physical']]
        if len(phys) == 0:
            ax.text(0.5, 0.5, 'No physical solutions', transform=ax.transAxes, ha='center')
            continue
        
        ann = phys.groupby('year').agg(
            FF=('FF', 'mean'), FF5=('FF', lambda x: np.percentile(x, 5)),
            FF95=('FF', lambda x: np.percentile(x, 95)),
            BB=('BB', 'mean'), BB5=('BB', lambda x: np.percentile(x, 5)),
            BB95=('BB', lambda x: np.percentile(x, 95)),
            Mic=('Mic', 'mean'), Mic5=('Mic', lambda x: np.percentile(x, 5)),
            Mic95=('Mic', lambda x: np.percentile(x, 95)),
        ).reset_index()
        
        for cat, color in [('FF', '#1f77b4'), ('BB', '#ff7f0e'), ('Mic', '#2ca02c')]:
            ax.plot(ann['year'], ann[cat], '-', color=color, lw=2, label=cat)
            ax.fill_between(ann['year'], ann[f'{cat}5'], ann[f'{cat}95'],
                          alpha=0.15, color=color)
        
        src_label = "Ben's δD" if tag == 'ben' else "Our δD"
        ax.set_title(f'{name} — {src_label}')
        ax.set_ylabel('Emissions (Tg/yr)')
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-50, 550)
        ax.axhline(0, color='k', ls=':', alpha=0.3)
        if row == 0 and col == 0:
            ax.legend()

for ax in axes[-1]:
    ax.set_xlabel('Year')
plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'best_models_detail.png', dpi=150, bbox_inches='tight')
plt.close()
print("  Saved: best_models_detail.png", flush=True)

print(f"\nAll outputs in: {OUTPUT_DIR}/")
print("Done!", flush=True)
