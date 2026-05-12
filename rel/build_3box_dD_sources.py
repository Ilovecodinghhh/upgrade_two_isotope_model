#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_3box_dD_sources.py
========================
Construct 3-box δD source signatures for Mic, BB, and FF sectors.

3-box structure:
  Box 1: NH extratropical  (90°N – 30°N)  → rows 0:60 on 1°×1° grid
  Box 2: Tropical           (30°N – 30°S)  → rows 60:120
  Box 3: SH extratropical  (30°S – 90°S)  → rows 120:180

This matches standard atmospheric transport box models (e.g., Bousquet 2006,
Thanwerdas 2024) and aligns with the Riddell-Young semi-hemispheric bands:
  PN (~90°N–30°N) ≈ Box 1
  TN+TS (~30°N–30°S) ≈ Box 2  
  PS (~30°S–90°S) ≈ Box 3

Also produces 3-box atmospheric δD from the station-level MC data
(using the semi-hemispheric results from dD_globmean.py).

Outputs (saved to rel/output/ and rel/data/):
  - {Mic,BB,FF}_dD_{NHext,Trop,SHext}_MC.csv  (24yr × 1001 cols)
  - ThreeBox_dD_sources_summary.csv
  - ThreeBox_atm_dD_annual.csv  (atmospheric δD for each box)

Methodology: Same as build_hemispheric_dD_sources.py but with 3 latitude
bands instead of 2. See that script for full documentation.
"""

import numpy as np
import pandas as pd
import xarray as xr
import netCDF4 as nc
import tifffile
from pathlib import Path
import warnings
import sys

warnings.filterwarnings('ignore')
np.random.seed(42)

N_MC = 1000
YEARS_CTCH4 = np.arange(1998, 2022)  # 24 years
N_YEARS = len(YEARS_CTCH4)

# 3-box latitude boundaries on a 180-row grid (row 0 = 89.5°N)
# Box 1: NH extratropical = 90°N to 30°N → rows 0:60
# Box 2: Tropical = 30°N to 30°S → rows 60:120
# Box 3: SH extratropical = 30°S to 90°S → rows 120:180
BOX_SLICES = {
    'NHext': slice(0, 60),    # 90°N – 30°N
    'Trop':  slice(60, 120),  # 30°N – 30°S
    'SHext': slice(120, 180), # 30°S – 90°S
}
BOX_NAMES = ['NHext', 'Trop', 'SHext']

# Paths
SCRIPT_DIR = Path(__file__).resolve().parent
BASE = SCRIPT_DIR.parent / "ImportantReferences" / "Riddell-Young2025PNAS_DS" / \
       "Riddell-Young_2025_MassBalancePackage" / "Riddell-Young_2025_MassBalancePackage"
DATA = BASE / "data"
OUT = SCRIPT_DIR / "output"
OUT.mkdir(exist_ok=True)


# ============================================================================
# HELPER: Downsample TIFF to 1°×1° and fill ocean/missing with lat averages
# ============================================================================

def load_MAT_grid():
    """Load and process mean annual temperature TIFF to 1°×1° grid."""
    with tifffile.TiffFile(str(DATA / "d2h_MA.tif")) as tif:
        image = tif.pages[0].asarray()
    
    SPECIAL = np.float32(-3.4e38)
    new_rows = np.full((77, image.shape[1]), SPECIAL, dtype=np.float32)
    image = np.vstack((new_rows, image))
    
    block_size = 12
    n_rows = image.shape[0] // block_size
    n_cols = image.shape[1] // block_size
    reshaped = image.reshape(n_rows, block_size, n_cols, block_size)
    
    averaged = np.zeros((n_rows, n_cols))
    for i in range(n_rows):
        for j in range(n_cols):
            block = reshaped[i, :, j, :].flatten()
            valid = block[block >= -1000]
            if len(valid) > 0:
                averaged[i, j] = np.mean(valid)
            else:
                averaged[i, j] = SPECIAL
    
    return averaged, SPECIAL


def fill_ocean_with_lat_avg(grid):
    """Fill ocean/missing cells with latitudinal averages (interpolated)."""
    result = grid.copy()
    mask = result < -1000
    masked_data = np.where(mask, np.nan, result)
    row_means = np.nanmean(masked_data, axis=1)
    
    for i in range(len(row_means)):
        if np.isnan(row_means[i]):
            above = i - 1
            while above >= 0 and np.isnan(row_means[above]):
                above -= 1
            below = i + 1
            while below < len(row_means) and np.isnan(row_means[below]):
                below += 1
            if above >= 0 and below < len(row_means):
                row_means[i] = np.mean([row_means[above], row_means[below]])
            elif above >= 0:
                row_means[i] = row_means[above]
            elif below < len(row_means):
                row_means[i] = row_means[below]
    
    result[mask] = np.tile(row_means, (360, 1)).T[mask]
    return result


# ============================================================================
# HELPER: 3-box emission-weighted mean δD
# ============================================================================

def threebox_emission_weighted_mean(dD_grid, emission_weight_grid):
    """
    Compute 3-box emission-weighted mean δD.
    Returns: dict with keys 'NHext', 'Trop', 'SHext'
    """
    results = {}
    for box_name, box_slice in BOX_SLICES.items():
        box_dD = dD_grid[box_slice, :]
        box_w = emission_weight_grid[box_slice, :]
        total = box_w.sum()
        if total > 0:
            results[box_name] = (box_dD * box_w).sum() / total
        else:
            results[box_name] = np.nan
    return results


# ============================================================================
# LOAD CTCH4 EMISSION GRIDS
# ============================================================================

def load_CTCH4_fluxes():
    """Load CTCH4 fluxes at 1°×1°, row 0 = 90°N. Shape: (24, 180, 360)."""
    f = nc.Dataset(str(DATA / "CTCH4_2023_flux3x2.nc"), 'r')
    mic_flux = f.variables['microbial_flux'][:]
    fos_flux = f.variables['fossil_flux'][:]
    pyr_flux = f.variables['pyrogenic_flux'][:]
    f.close()
    
    def expand_to_1deg(flux_3x2):
        expanded = np.repeat(flux_3x2, 2, axis=1)
        expanded = np.repeat(expanded, 3, axis=2)
        expanded /= 6.0
        return expanded
    
    def annual_sum(monthly):
        return monthly.reshape(N_YEARS, 12, 180, 360).sum(axis=1)
    
    mic_ann = annual_sum(expand_to_1deg(mic_flux))[:, ::-1, :]
    fos_ann = annual_sum(expand_to_1deg(fos_flux))[:, ::-1, :]
    pyr_ann = annual_sum(expand_to_1deg(pyr_flux))[:, ::-1, :]
    
    return mic_ann, fos_ann, pyr_ann


# ============================================================================
# 1. MICROBIAL δD — 3-BOX MC
# ============================================================================

def compute_mic_dD_3box(MAT_1deg, mic_ann):
    """Douglas et al. 2021: δD = 0.6088*MAT - 285.7"""
    print("Computing Mic δD 3-box MC...")
    SLOPE, SLOPE_U = 0.6088, 0.072
    INTERCEPT, INTERCEPT_U = 285.7, 6.9
    
    results = {box: np.zeros((N_YEARS, N_MC)) for box in BOX_NAMES}
    
    for k in range(N_MC):
        s = SLOPE + np.random.normal() * SLOPE_U
        i = INTERCEPT + np.random.normal() * INTERCEPT_U
        dD_mc = fill_ocean_with_lat_avg(s * MAT_1deg - i)
        
        for yr_idx in range(N_YEARS):
            em = mic_ann[yr_idx]
            em_total = em.sum()
            if em_total <= 0:
                for box in BOX_NAMES:
                    results[box][yr_idx, k] = np.nan
                continue
            em_w = em / em_total
            box_vals = threebox_emission_weighted_mean(dD_mc, em_w)
            for box in BOX_NAMES:
                results[box][yr_idx, k] = box_vals[box]
        
        if (k + 1) % 200 == 0:
            print(f"  Mic MC: {k+1}/{N_MC}")
    
    return results


# ============================================================================
# 2. BIOMASS BURNING δD — 3-BOX MC
# ============================================================================

def compute_bb_dD_3box(MAT_1deg, pyr_ann):
    """Umezawa et al. 2011: δD = 1.16*MAT - 177"""
    print("Computing BB δD 3-box MC...")
    SLOPE, SLOPE_U = 1.16, 0.09
    INTERCEPT, INTERCEPT_U = 177.0, 6.5
    
    results = {box: np.zeros((N_YEARS, N_MC)) for box in BOX_NAMES}
    
    for k in range(N_MC):
        s = SLOPE + np.random.normal() * SLOPE_U
        i = INTERCEPT + np.random.normal() * INTERCEPT_U
        dD_mc = fill_ocean_with_lat_avg(s * MAT_1deg - i)
        
        for yr_idx in range(N_YEARS):
            em = pyr_ann[yr_idx]
            em_total = em.sum()
            if em_total <= 0:
                for box in BOX_NAMES:
                    results[box][yr_idx, k] = np.nan
                continue
            em_w = em / em_total
            box_vals = threebox_emission_weighted_mean(dD_mc, em_w)
            for box in BOX_NAMES:
                results[box][yr_idx, k] = box_vals[box]
        
        if (k + 1) % 200 == 0:
            print(f"  BB MC: {k+1}/{N_MC}")
    
    return results


# ============================================================================
# 3. FOSSIL FUEL δD — 3-BOX MC
# ============================================================================

def compute_ff_dD_3box(fos_ann):
    """Country-level ONG+coal δD, EDGAR-weighted, assigned to 3 boxes."""
    print("Computing FF δD 3-box MC...")
    
    coal_dD = pd.read_csv(DATA / "coal_dD.csv")
    ONG_dD = pd.read_csv(DATA / "ONG_dD.csv")
    Country_ONG_emis = pd.read_csv(DATA / "EDGAR8_ONG.csv")
    Country_Coal_emis = pd.read_csv(DATA / "EDGAR8_Coal.csv")
    
    NAME_MAP = {
        'russian federation': 'russia',
        'iran, islamic republic of': 'iran',
        'korea, republic of': 'south korea',
        'korea, dem. people\'s rep.': 'north korea',
        'taiwan, province of china': 'taiwan',
        'trinidad and tobago': 'trinidad',
        'the netherlands': 'netherlands',
        'viet nam': 'vietnam',
        'venezuela, bolivarian republic of': 'venezuela',
        'bolivia, plurinational state of': 'bolivia',
        'congo, democratic republic of the': 'democratic republic of the congo',
        'côte d\'ivoire': 'ivory coast',
        'czech republic': 'czechia',
        'lao people\'s democratic republic': 'laos',
        'united republic of tanzania': 'tanzania',
        'brunei darussalam': 'brunei',
        'syrian arab republic': 'syria',
        'myanmar/burma': 'myanmar',
    }
    
    def normalize_name(name):
        n = str(name).strip().lower()
        return NAME_MAP.get(n, n)
    
    # Build δD lookups
    ong_dD_lookup, ong_dD_unc_lookup = {}, {}
    coal_dD_lookup, coal_dD_unc_lookup = {}, {}
    
    for _, row in ONG_dD.iterrows():
        name = normalize_name(row.iloc[0])
        mean_val = row['mean']
        if pd.isna(mean_val):
            continue
        ong_dD_lookup[name] = float(mean_val)
        ong_dD_unc_lookup[name] = float(row['std']) if not pd.isna(row['std']) else 10.0
    
    for _, row in coal_dD.iterrows():
        name = normalize_name(row.iloc[0])
        mean_val = row['mean']
        if pd.isna(mean_val):
            continue
        coal_dD_lookup[name] = float(mean_val)
        coal_dD_unc_lookup[name] = float(row['std']) if not pd.isna(row['std']) else 10.0
    
    glob_ong_dD = np.nanmean(list(ong_dD_lookup.values()))
    glob_coal_dD = np.nanmean(list(coal_dD_lookup.values()))
    glob_ong_unc = np.nanmean(list(ong_dD_unc_lookup.values()))
    glob_coal_unc = np.nanmean(list(coal_dD_unc_lookup.values()))
    
    print(f"  ONG δD: {len(ong_dD_lookup)} countries, Coal δD: {len(coal_dD_lookup)} countries")
    print(f"  Global ONG δD: {glob_ong_dD:.1f}‰, Coal δD: {glob_coal_dD:.1f}‰")
    
    # Country → 3-box assignment based on centroid latitude
    try:
        import geopandas as gpd
        world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    except Exception:
        world = None
    
    # Known latitude bands for major emitters (fallback)
    COUNTRY_BOX_OVERRIDE = {
        # SH extratropical (< -30°)
        'australia': 'SHext', 'new zealand': 'SHext', 'argentina': 'SHext',
        'chile': 'SHext', 'south africa': 'SHext', 'uruguay': 'SHext',
        # Tropical (-30° to 30°)
        'brazil': 'Trop', 'indonesia': 'Trop', 'nigeria': 'Trop',
        'venezuela': 'Trop', 'angola': 'Trop', 'colombia': 'Trop',
        'india': 'Trop', 'thailand': 'Trop', 'egypt': 'Trop',
        'trinidad': 'Trop', 'mozambique': 'Trop', 'namibia': 'Trop',
        'ecuador': 'Trop', 'peru': 'Trop', 'bolivia': 'Trop',
        'malaysia': 'Trop', 'myanmar': 'Trop', 'vietnam': 'Trop',
        'qatar': 'Trop', 'united arab emirates': 'Trop',
        'oman': 'Trop', 'saudi arabia': 'Trop', 'yemen': 'Trop',
        'iraq': 'Trop',
        'mexico': 'Trop',
        # NH extratropical (> 30°)
        'russia': 'NHext', 'canada': 'NHext', 'united states': 'NHext',
        'china': 'NHext', 'germany': 'NHext', 'united kingdom': 'NHext',
        'norway': 'NHext', 'poland': 'NHext', 'japan': 'NHext',
        'iran': 'NHext', 'turkey': 'NHext', 'ukraine': 'NHext',
        'romania': 'NHext', 'italy': 'NHext', 'france': 'NHext',
        'austria': 'NHext', 'denmark': 'NHext', 'lithuania': 'NHext',
        'netherlands': 'NHext', 'taiwan': 'NHext', 'south korea': 'NHext',
        'turkmenistan': 'NHext', 'uzbekistan': 'NHext', 'kazakhstan': 'NHext',
        'algeria': 'NHext', 'libya': 'NHext',
    }
    
    def assign_box(country_name):
        cn = country_name.lower().strip()
        if cn in COUNTRY_BOX_OVERRIDE:
            return COUNTRY_BOX_OVERRIDE[cn]
        if world is not None:
            for _, row in world.iterrows():
                if row['name'].lower() == cn or normalize_name(row['name']) == cn:
                    lat = row.geometry.centroid.y
                    if lat > 30:
                        return 'NHext'
                    elif lat < -30:
                        return 'SHext'
                    else:
                        return 'Trop'
        return 'NHext'  # Default
    
    edgar_country_col_ong = Country_ONG_emis.columns[0]
    edgar_country_col_coal = Country_Coal_emis.columns[0]
    ong_year_cols = [c for c in Country_ONG_emis.columns[1:] if str(c).replace('.', '').isdigit()]
    coal_year_cols = [c for c in Country_Coal_emis.columns[1:] if str(c).replace('.', '').isdigit()]
    
    def get_edgar_year_col(year, year_cols):
        year_vals = [int(float(c)) for c in year_cols]
        closest = min(year_vals, key=lambda y: abs(y - year))
        return str(closest)
    
    results = {box: np.zeros((N_YEARS, N_MC)) for box in BOX_NAMES}
    
    for yr_idx, year in enumerate(YEARS_CTCH4):
        ong_yr_col = get_edgar_year_col(year, ong_year_cols)
        coal_yr_col = get_edgar_year_col(year, coal_year_cols)
        
        # Collect per-box data
        box_data = {box: {'emis': [], 'dD': [], 'unc': []} for box in BOX_NAMES}
        
        for _, row in Country_ONG_emis.iterrows():
            country = normalize_name(row[edgar_country_col_ong])
            emis = row[ong_yr_col] if ong_yr_col in row.index else 0
            if pd.isna(emis) or emis <= 0:
                continue
            box = assign_box(country)
            dD_val = ong_dD_lookup.get(country, glob_ong_dD)
            dD_unc = ong_dD_unc_lookup.get(country, glob_ong_unc)
            box_data[box]['emis'].append(emis)
            box_data[box]['dD'].append(dD_val)
            box_data[box]['unc'].append(dD_unc)
        
        for _, row in Country_Coal_emis.iterrows():
            country = normalize_name(row[edgar_country_col_coal])
            emis = row[coal_yr_col] if coal_yr_col in row.index else 0
            if pd.isna(emis) or emis <= 0:
                continue
            box = assign_box(country)
            dD_val = coal_dD_lookup.get(country, glob_coal_dD)
            dD_unc = coal_dD_unc_lookup.get(country, glob_coal_unc)
            box_data[box]['emis'].append(emis)
            box_data[box]['dD'].append(dD_val)
            box_data[box]['unc'].append(dD_unc)
        
        for k in range(N_MC):
            for box in BOX_NAMES:
                if not box_data[box]['emis']:
                    results[box][yr_idx, k] = np.nan
                    continue
                emis_arr = np.array(box_data[box]['emis'])
                dD_arr = np.array(box_data[box]['dD'])
                unc_arr = np.array(box_data[box]['unc'])
                dD_perturbed = dD_arr + np.random.normal(size=len(dD_arr)) * unc_arr
                results[box][yr_idx, k] = np.average(dD_perturbed, weights=emis_arr)
        
        if (yr_idx + 1) % 6 == 0:
            vals = {b: results[b][yr_idx].mean() for b in BOX_NAMES}
            print(f"  FF year {year}: NHext={vals['NHext']:.1f}, Trop={vals['Trop']:.1f}, SHext={vals['SHext']:.1f}‰")
    
    return results


# ============================================================================
# 3-BOX ATMOSPHERIC δD from semi-hemispheric data
# ============================================================================

def build_3box_atmospheric_dD():
    """
    Construct 3-box atmospheric δD from the semi-hemispheric results.
    PN ≈ NHext, TN+TS ≈ Trop, PS ≈ SHext.
    
    For Trop: weighted average of TN and TS (equal weight since they
    span similar latitude ranges: 30°N-0° and 0°-30°S).
    """
    print("\nBuilding 3-box atmospheric δD from semi-hemispheric data...")
    
    sem_file = SCRIPT_DIR / "data" / "SemiHemMean_dD_dei_DasguptaCal_noBUDS.csv"
    if not sem_file.exists():
        print(f"  WARNING: {sem_file} not found, skipping atmospheric δD")
        return
    
    df = pd.read_csv(sem_file)
    
    # Compute annual means
    df['year_int'] = df['Year'].astype(int)
    annual = df.groupby('year_int').agg(
        NHext=('PN_smooth_mean', 'mean'),
        Trop_TN=('TN_smooth_mean', 'mean'),
        Trop_TS=('TS_smooth_mean', 'mean'),
        SHext=('PS_smooth_mean', 'mean'),
    ).reset_index()
    
    # Tropical = average of TN and TS
    annual['Trop'] = (annual['Trop_TN'] + annual['Trop_TS']) / 2.0
    
    out = annual[['year_int', 'NHext', 'Trop', 'SHext']].rename(
        columns={'year_int': 'year'})
    
    out_path = OUT / "ThreeBox_atm_dD_annual.csv"
    out.to_csv(out_path, index=False, float_format='%.2f')
    
    # Also copy to data/
    out.to_csv(SCRIPT_DIR / "data" / "ThreeBox_atm_dD_annual.csv",
               index=False, float_format='%.2f')
    
    print(f"  Saved ThreeBox_atm_dD_annual.csv ({len(out)} years)")
    print(f"  NHext range: {out['NHext'].min():.1f} to {out['NHext'].max():.1f}‰")
    print(f"  Trop range:  {out['Trop'].min():.1f} to {out['Trop'].max():.1f}‰")
    print(f"  SHext range: {out['SHext'].min():.1f} to {out['SHext'].max():.1f}‰")
    
    # Also build MC iterations for atmospheric δD per box
    # We need the per-station MC data — check if NH/SH MC iterations exist
    # and split them further, or use the PN/TN/TS/PS weekly data as-is
    # For now, provide annual means (MC atmospheric iterations would require
    # re-running dD_globmean.py with 3-box output, which is a larger task)
    print("  NOTE: 3-box atmospheric δD MC iterations require re-running "
          "dD_globmean.py with per-box station groupings. For now, annual "
          "means from semi-hemispheric smoothed curves are provided.")


# ============================================================================
# SAVE
# ============================================================================

def save_mc_csv(filename, years, mc_data):
    out_data = np.column_stack([years, mc_data])
    header = "year," + ",".join([f"mc_{i}" for i in range(N_MC)])
    np.savetxt(OUT / filename, out_data, delimiter=',', header=header,
               comments='', fmt=['%d'] + ['%.3f'] * N_MC)
    # Also copy to data/
    np.savetxt(SCRIPT_DIR / "data" / filename, out_data, delimiter=',',
               header=header, comments='', fmt=['%d'] + ['%.3f'] * N_MC)
    print(f"  Saved {filename}: {out_data.shape}")


def save_summary(mic_res, bb_res, ff_res):
    rows = []
    for yr_idx, year in enumerate(YEARS_CTCH4):
        row = {'year': int(year)}
        for sector, res in [('Mic', mic_res), ('BB', bb_res), ('FF', ff_res)]:
            for box in BOX_NAMES:
                row[f'{sector}_dD_{box}_mean'] = np.nanmean(res[box][yr_idx])
                row[f'{sector}_dD_{box}_std'] = np.nanstd(res[box][yr_idx])
        rows.append(row)
    
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "ThreeBox_dD_sources_summary.csv", index=False, float_format='%.2f')
    df.to_csv(SCRIPT_DIR / "data" / "ThreeBox_dD_sources_summary.csv",
              index=False, float_format='%.2f')
    
    print(f"\n3-Box Source Signature Summary (1998-2021 means):")
    for sector, res in [('Mic', mic_res), ('BB', bb_res), ('FF', ff_res)]:
        means = {box: np.nanmean([np.nanmean(res[box][i]) for i in range(N_YEARS)]) for box in BOX_NAMES}
        stds = {box: np.nanmean([np.nanstd(res[box][i]) for i in range(N_YEARS)]) for box in BOX_NAMES}
        print(f"  {sector} δD:  NHext = {means['NHext']:.1f} ± {stds['NHext']:.1f}  "
              f"Trop = {means['Trop']:.1f} ± {stds['Trop']:.1f}  "
              f"SHext = {means['SHext']:.1f} ± {stds['SHext']:.1f}‰")


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("Building 3-box δD source signatures")
    print("  NHext (90°N–30°N) | Trop (30°N–30°S) | SHext (30°S–90°S)")
    print("=" * 70)
    
    print("\n[1/6] Loading MAT TIFF...")
    MAT_1deg, _ = load_MAT_grid()
    print(f"  MAT grid: {MAT_1deg.shape}, {(MAT_1deg > -1000).sum()} valid cells")
    
    print("\n[2/6] Loading CTCH4 fluxes...")
    mic_ann, fos_ann, pyr_ann = load_CTCH4_fluxes()
    
    # Show 3-box emission fractions
    for name, flux in [('Mic', mic_ann), ('Fossil', fos_ann), ('Pyrogenic', pyr_ann)]:
        total = flux.sum()
        for box, sl in BOX_SLICES.items():
            frac = flux[:, sl, :].sum() / total
            print(f"  {name} {box}: {frac:.1%}")
    
    print("\n[3/6] Microbial δD (3-box)...")
    mic_res = compute_mic_dD_3box(MAT_1deg, mic_ann)
    
    print("\n[4/6] Biomass Burning δD (3-box)...")
    bb_res = compute_bb_dD_3box(MAT_1deg, pyr_ann)
    
    print("\n[5/6] Fossil Fuel δD (3-box)...")
    ff_res = compute_ff_dD_3box(fos_ann)
    
    print("\n[6/6] Atmospheric δD (3-box)...")
    build_3box_atmospheric_dD()
    
    print("\nSaving source signature results...")
    for sector, res in [('Mic', mic_res), ('BB', bb_res), ('FF', ff_res)]:
        for box in BOX_NAMES:
            save_mc_csv(f"{sector}_dD_{box}_MC.csv", YEARS_CTCH4, res[box])
    save_summary(mic_res, bb_res, ff_res)
    
    print("\n" + "=" * 70)
    print("Done! 3-box δD source signatures saved to rel/output/ and rel/data/")
    print("=" * 70)


if __name__ == "__main__":
    main()
