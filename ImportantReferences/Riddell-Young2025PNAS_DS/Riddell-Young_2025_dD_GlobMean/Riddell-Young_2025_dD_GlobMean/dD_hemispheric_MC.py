#!/usr/bin/env python3
"""
Construct hemispheric (NH/SH) δD-CH4 annual MC iterations from
Riddell-Young 2025's per-station smoothed MC output files.

Reads: output/*_smoothedMC.txt (from MBL_calc_Unc.py)
       data/siteinfo_all_ch4h2.txt (site coords)
       output/GlobMean_dD_iterations_DasguptaCal_noBUDS.xlsx (year vector)

Saves to:
  output/NHMean_dD_iterations_DasguptaCal_noBUDS.xlsx
  output/SHMean_dD_iterations_DasguptaCal_noBUDS.xlsx
  output/HemMean_dD_annual_DasguptaCal_noBUDS.csv
"""
import os, glob
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
os.chdir(BASE)

# ── Parse siteinfo (pipe-delimited) ──
site_lats = {}
with open('data/siteinfo_all_ch4h2.txt') as f:
    for line in f:
        parts = [p.strip() for p in line.split('|')]
        if len(parts) >= 4:
            site_code = parts[0]
            try:
                lat = float(parts[3])
                site_lats[site_code] = lat
            except ValueError:
                pass

print(f"Parsed {len(site_lats)} sites from siteinfo")

# ── Load station MC files ──
mc_files = sorted(glob.glob('output/*_01D0_smoothedMC.txt'))
print(f"Found {len(mc_files)} smoothedMC files")

stations = []
for fpath in mc_files:
    fname = os.path.basename(fpath)
    # Extract site code: everything before _01D0_smoothedMC.txt
    site_code = fname.replace('_01D0_smoothedMC.txt', '')
    # Try to match to siteinfo
    lat = site_lats.get(site_code)
    if lat is None:
        # Try removing network suffix (e.g., "brwIMAU" -> "brw" won't work, 
        # but "brwIMAU" might be in siteinfo directly)
        # Just skip if not found
        print(f"  WARNING: no lat for {site_code}, skipping")
        continue
    arr = np.loadtxt(fpath)
    stations.append({'code': site_code, 'lat': lat, 'data': arr})

print(f"Loaded {len(stations)} stations with coordinates")

# ── Get year vector ──
df_glob = pd.read_excel('output/GlobMean_dD_iterations_DasguptaCal_noBUDS.xlsx', header=None)
year = df_glob.iloc[:, 0].values.astype(int)
n_years = len(year)
iterations = df_glob.shape[1] - 1  # first col is year
print(f"Years: {year[0]}–{year[-1]} ({n_years} years), {iterations} MC iterations")

# Date vector — use the longest station file for maximum coverage
max_len = max(s['data'].shape[0] for s in stations)
ref_station = [s for s in stations if s['data'].shape[0] == max_len][0]
date_full = ref_station['data'][:, 0]
n_weeks = len(date_full)

# Split NH / SH
NH_stations = [s for s in stations if s['lat'] >= 0]
SH_stations = [s for s in stations if s['lat'] < 0]
print(f"NH: {len(NH_stations)} stations, SH: {len(SH_stations)} stations")
print(f"Weekly date vector length: {n_weeks}")

# ── Compute annual NH/SH means per MC iteration ──
smoothedNH = np.full((n_years, iterations), np.nan)
smoothedSH = np.full((n_years, iterations), np.nan)

for k in range(iterations):
    if k % 200 == 0:
        print(f"  MC iteration {k}/{iterations}")
    
    col = k + 1  # column 0 is date, 1..1000 are MC iterations
    
    # NH weekly mean — align all stations to the common date grid
    nh_vals = np.full((len(NH_stations), n_weeks), np.nan)
    for si, s in enumerate(NH_stations):
        if col < s['data'].shape[1]:
            slen = s['data'].shape[0]
            # Find where this station's dates fall in the reference grid
            s_dates = s['data'][:, 0]
            # Quick approach: find offset by matching first date
            offset = np.searchsorted(date_full, s_dates[0])
            end = min(offset + slen, n_weeks)
            actual_len = end - offset
            nh_vals[si, offset:end] = s['data'][:actual_len, col]
    nh_mean = np.nanmean(nh_vals, axis=0)
    
    # SH weekly mean
    sh_vals = np.full((len(SH_stations), n_weeks), np.nan)
    for si, s in enumerate(SH_stations):
        if col < s['data'].shape[1]:
            slen = s['data'].shape[0]
            s_dates = s['data'][:, 0]
            offset = np.searchsorted(date_full, s_dates[0])
            end = min(offset + slen, n_weeks)
            actual_len = end - offset
            sh_vals[si, offset:end] = s['data'][:actual_len, col]
    sh_mean = np.nanmean(sh_vals, axis=0)
    
    # Annual averages
    for yi, y in enumerate(year):
        mask = (date_full >= y) & (date_full < y + 1)
        if mask.sum() > 0:
            smoothedNH[yi, k] = np.nanmean(nh_mean[mask])
            smoothedSH[yi, k] = np.nanmean(sh_mean[mask])

# ── Save ──
# NH MC iterations
combined_NH = np.column_stack((year, smoothedNH))
pd.DataFrame(combined_NH).to_excel(
    'output/NHMean_dD_iterations_DasguptaCal_noBUDS.xlsx',
    index=False, header=False)

# SH MC iterations
combined_SH = np.column_stack((year, smoothedSH))
pd.DataFrame(combined_SH).to_excel(
    'output/SHMean_dD_iterations_DasguptaCal_noBUDS.xlsx',
    index=False, header=False)

# Annual summary
NH_mean_ann = np.nanmean(smoothedNH, axis=1)
SH_mean_ann = np.nanmean(smoothedSH, axis=1)
df_annual = pd.DataFrame({
    'Year': year,
    'NH_mean': NH_mean_ann,
    'NH_std': np.nanstd(smoothedNH, axis=1),
    'SH_mean': SH_mean_ann,
    'SH_std': np.nanstd(smoothedSH, axis=1),
    'NH_SH_diff': NH_mean_ann - SH_mean_ann,
})
df_annual.to_csv('output/HemMean_dD_annual_DasguptaCal_noBUDS.csv', index=False)

print("\n✓ Saved:")
print(f"  NHMean_dD_iterations_DasguptaCal_noBUDS.xlsx  ({n_years}×{iterations})")
print(f"  SHMean_dD_iterations_DasguptaCal_noBUDS.xlsx  ({n_years}×{iterations})")
print(f"  HemMean_dD_annual_DasguptaCal_noBUDS.csv")
print(f"\nNH–SH gradient summary:")
print(df_annual[['Year', 'NH_mean', 'SH_mean', 'NH_SH_diff']].to_string(index=False))
