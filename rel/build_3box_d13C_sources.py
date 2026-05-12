#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_3box_d13C_sources.py
==========================
Construct 3-box δ¹³C source signatures AND atmospheric δ¹³C for Mic, BB,
and FF sectors.

3-box structure (matching build_3box_dD_sources.py):
  Box 1: NH extratropical  (90°N – 30°N)  → rows 0:60 on 1°×1° grid
  Box 2: Tropical           (30°N – 30°S)  → rows 60:120
  Box 3: SH extratropical  (30°S – 90°S)  → rows 120:180

Methodology: Same as build_hemispheric_d13C_sources.py but with 3 latitude
bands instead of 2.

Outputs (saved to rel/data/):
  Source signatures:
    - {FF,Mic,BB}_d13C_{NHext,Trop,SHext}_MC.csv  (24yr × 1001 cols)
    - ThreeBox_d13C_sources_summary.csv
  Atmospheric δ¹³C:
    - ThreeBox_atm_d13C_annual.csv
    - ThreeBox_atm_d13C_MC.npz  (MC iterations per box)
"""

import numpy as np
import pandas as pd
import netCDF4 as nc
import xarray as xr
from pathlib import Path
import warnings
import sys

warnings.filterwarnings('ignore')
np.random.seed(42)

N_MC = 1000
YEARS = np.arange(1998, 2022)  # 24 years
N_YEARS = len(YEARS)

# 3-box latitude boundaries on a 180-row grid (row 0 = 89.5°N)
BOX_SLICES = {
    'NHext': slice(0, 60),    # 90°N – 30°N
    'Trop':  slice(60, 120),  # 30°N – 30°S
    'SHext': slice(120, 180), # 30°S – 90°S
}
BOX_NAMES = ['NHext', 'Trop', 'SHext']

SCRIPT_DIR = Path(__file__).resolve().parent
BASE = SCRIPT_DIR.parent / "ImportantReferences" / "Riddell-Young2025PNAS_DS" / \
       "Riddell-Young_2025_MassBalancePackage" / "Riddell-Young_2025_MassBalancePackage"
DATA = BASE / "data"
OUT = SCRIPT_DIR / "data"
OUT.mkdir(exist_ok=True)

# Suess effect trend
SUESS_TREND = -0.024  # ‰/yr
SUESS_TREND_UNC = 0.005


# ============================================================================
# DATA LOADING (shared with hemispheric builder)
# ============================================================================

def load_CTCH4():
    """Load CTCH4 3×2° → expand to 1°×1°, flip so row0=90°N."""
    f = nc.Dataset(str(DATA / "CTCH4_2023_flux3x2.nc"), 'r')
    mic = f.variables['microbial_flux'][:]
    fos = f.variables['fossil_flux'][:]
    pyr = f.variables['pyrogenic_flux'][:]
    f.close()

    def expand(flux_3x2):
        e = np.repeat(flux_3x2, 2, axis=1)
        e = np.repeat(e, 3, axis=2) / 6.0
        return e

    mic_1 = expand(mic); fos_1 = expand(fos); pyr_1 = expand(pyr)

    def annual(m):
        return m.reshape(N_YEARS, 12, 180, 360).sum(axis=1)

    mic_ann = annual(mic_1)[:, ::-1, :]
    fos_ann = annual(fos_1)[:, ::-1, :]
    pyr_ann = annual(pyr_1)[:, ::-1, :]
    return mic_ann, fos_ann, pyr_ann


def load_C3C4():
    """Load C4 fraction maps. Tries Luo 2023 first; falls back to Still 2003."""
    luo_path = DATA / "C4_distribution_NUS_v2.2.nc"
    # Auto-reassemble from split parts if needed
    if not luo_path.exists():
        reassemble = DATA / "reassemble_luo_c4.sh"
        if reassemble.exists():
            import subprocess
            subprocess.run(["bash", str(reassemble)], check=True)
    if luo_path.exists():
        C3C4_data = nc.Dataset(str(luo_path), 'r')
        C4_maps = C3C4_data.variables['C4_area'][:]
        C4_maps = np.nan_to_num(C4_maps)
        C4_c = C4_maps.reshape(19, 360, 2, 180, 2).mean(axis=(2, 4))
        C4_c = np.transpose(C4_c, (0, 2, 1))
        C4exp = np.concatenate([
            np.repeat(C4_c[0:1, :180, :360], 3, axis=0),
            C4_c[:, :180, :360],
            np.repeat(C4_c[-1:, :180, :360], 2, axis=0)
        ], axis=0)
        C4exp = np.flip(C4exp, axis=1)
        print("  Using Luo 2023 C4 map (time-varying)")
    else:
        C4_still = pd.read_excel(DATA / "Still2003_C4.xlsx", header=None).values
        C4_still = np.flipud(C4_still)
        C4_still = np.where(C4_still < 0, 0, C4_still)
        C4exp = np.repeat(C4_still[np.newaxis, :180, :360], N_YEARS, axis=0)
        print("  Using Still 2003 C4 map (static fallback)")
    return C4exp


# ============================================================================
# WETLAND δ¹³C LOADING (isotem per-year files or Oh 2022 fallback)
# ============================================================================

def load_wetland_d13C():
    """Load isotem wetland δ¹³C spatial maps (per-year NC4 files).
    Falls back to Oh 2022 global time series if unavailable.
    Returns (grid, unc_grid) each (24, 180, 360)."""
    isotem_dir = SCRIPT_DIR / "data" / "isotem_wetland_d13C-CH4"
    if isotem_dir.is_dir():
        result = np.full((N_YEARS, 180, 360), np.nan)
        unc_result = np.full((N_YEARS, 180, 360), np.nan)
        last_valid = None
        
        for yr_idx, year in enumerate(YEARS):
            fpath = isotem_dir / f"isotem_wetland_d13C-CH4_{year}.nc4"
            if fpath.exists():
                wf = nc.Dataset(str(fpath), 'r')
                d13c = wf.variables['wetland_d13C-CH4'][:]  # (12, 720, 360) = (month, lon, lat)
                wf.close()
                d13c = np.ma.filled(d13c, fill_value=np.nan)
                ann = np.nanmean(d13c, axis=0)  # (720, 360)
                ann_std = np.nanstd(d13c, axis=0)
                ann_ll = ann.T     # (360, 720), lat 90N→89.5S
                std_ll = ann_std.T
                ann_1deg = np.nanmean(ann_ll.reshape(180, 2, 360, 2), axis=(1, 3))
                std_1deg = np.nanmean(std_ll.reshape(180, 2, 360, 2), axis=(1, 3))
                # Keep NaN for cells with no wetlands — emission-weighting handles this
                result[yr_idx] = ann_1deg
                unc_result[yr_idx] = std_1deg
                last_valid = yr_idx
            elif last_valid is not None:
                result[yr_idx] = result[last_valid]
                unc_result[yr_idx] = unc_result[last_valid]
        
        first_valid = next((i for i in range(N_YEARS) if not np.all(np.isnan(result[i]))), 0)
        for i in range(first_valid):
            result[i] = result[first_valid]
            unc_result[i] = unc_result[first_valid]
        
        nh_mean = np.nanmean(result[:, :90, :])
        sh_mean = np.nanmean(result[:, 90:, :])
        print(f"  Using isotem wetland δ¹³C spatial maps")
        print(f"    NH mean: {nh_mean:.2f}‰, SH mean: {sh_mean:.2f}‰, gap: {nh_mean - sh_mean:.2f}‰")
        return result, unc_result
    
    oh_path = DATA / "Oh_2022_Wetlands.xlsx"
    if oh_path.exists():
        oh = pd.read_excel(oh_path).values
        wetland_d13C_ts = oh[:, 1]
        print(f"  Using Oh 2022 wetland δ¹³C global means: {np.mean(wetland_d13C_ts):.1f}‰")
        result = np.zeros((N_YEARS, 180, 360))
        unc = np.full((N_YEARS, 180, 360), 0.7)
        for yr in range(N_YEARS):
            result[yr, :, :] = wetland_d13C_ts[yr]
        return result, unc
    print("  WARNING: No wetland δ¹³C data, using -61.0‰")
    return np.full((N_YEARS, 180, 360), -61.0), np.full((N_YEARS, 180, 360), 1.0)


# ============================================================================
# HELPER: 3-box emission-weighted mean
# ============================================================================

def threebox_emission_weighted_mean(value_grid, weight_grid):
    """Compute 3-box emission-weighted mean. Returns dict {box: scalar}."""
    results = {}
    for box_name, box_slice in BOX_SLICES.items():
        box_v = value_grid[box_slice, :]
        box_w = weight_grid[box_slice, :]
        total = box_w.sum()
        if total > 0:
            results[box_name] = (box_v * box_w).sum() / total
        else:
            results[box_name] = np.nan
    return results


# ============================================================================
# 1. BB δ¹³C — 3-BOX MC
# ============================================================================

def compute_bb_d13C_3box(C4exp, pyr_ann):
    """BB δ¹³C = C3_frac × C3_δ¹³C + C4_frac × C4_δ¹³C, emission-weighted per box."""
    print("Computing BB δ¹³C 3-box MC...")

    C3_D13C = -26.8; C3_D13C_STD = 2.9
    C4_D13C = -12.7; C4_D13C_STD = 4.6

    results = {box: np.zeros((N_YEARS, N_MC)) for box in BOX_NAMES}

    for k in range(N_MC):
        c3_val = C3_D13C + np.random.normal() * C3_D13C_STD
        c4_val = C4_D13C + np.random.normal() * C4_D13C_STD

        for yr in range(N_YEARS):
            C4_frac = C4exp[yr] / 100.0
            d13C_map = C4_frac * c4_val + (1 - C4_frac) * c3_val
            em = pyr_ann[yr]
            em_total = em.sum()
            if em_total <= 0:
                for box in BOX_NAMES:
                    results[box][yr, k] = c3_val * 0.7 + c4_val * 0.3
                continue
            em_w = em / em_total
            box_vals = threebox_emission_weighted_mean(d13C_map, em_w)
            for box in BOX_NAMES:
                results[box][yr, k] = box_vals[box]

        if (k + 1) % 200 == 0:
            print(f"  BB MC: {k+1}/{N_MC}")

    return results


# ============================================================================
# 2. FF δ¹³C — 3-BOX MC
# ============================================================================

def compute_ff_d13C_3box(fos_ann):
    """FF δ¹³C: country-level ONG + coal, emission-weighted, 3-box split."""
    print("Computing FF δ¹³C 3-box MC...")

    coal_d13C = pd.read_csv(DATA / "coal_d13C.csv")
    ONG_d13C = pd.read_csv(DATA / "ONG_d13C.csv")
    Country_ONG_emis = pd.read_csv(DATA / "EDGAR8_ONG.csv")
    Country_Coal_emis = pd.read_csv(DATA / "EDGAR8_Coal.csv")

    # Temporal ONG trends
    us_ong = pd.read_csv(DATA / "US_ONG_trends.csv").iloc[:53, 2].values
    cc_ong = pd.read_csv(DATA / "China_Canada_ONG_Trends.csv")
    china_ong = cc_ong.iloc[:, 1].values
    canada_ong = cc_ong.iloc[:, 3].values

    NAME_MAP = {
        'russian federation': 'russia', 'iran, islamic republic of': 'iran',
        'korea, republic of': 'south korea', "korea, dem. people's rep.": 'north korea',
        'viet nam': 'vietnam', 'venezuela, bolivarian republic of': 'venezuela',
        'bolivia, plurinational state of': 'bolivia',
        "côte d'ivoire": 'ivory coast', 'czech republic': 'czechia',
        'brunei darussalam': 'brunei', 'syrian arab republic': 'syria',
        'the netherlands': 'netherlands',
    }

    def norm(name):
        n = str(name).strip().lower()
        return NAME_MAP.get(n, n)

    ong_lookup = {}; ong_unc = {}
    for _, row in ONG_d13C.iterrows():
        name = norm(row.iloc[0])
        if pd.notna(row['mean']):
            ong_lookup[name] = float(row['mean'])
            ong_unc[name] = float(row['std']) if pd.notna(row['std']) else 5.0

    coal_lookup = {}; coal_unc = {}
    for _, row in coal_d13C.iterrows():
        name = norm(row.iloc[0])
        if pd.notna(row['mean']):
            coal_lookup[name] = float(row['mean'])
            coal_unc[name] = float(row['std']) if pd.notna(row['std']) else 5.0

    glob_ong_d13C = np.nanmean(list(ong_lookup.values()))
    glob_coal_d13C = np.nanmean(list(coal_lookup.values()))
    glob_ong_u = np.nanmean(list(ong_unc.values()))
    glob_coal_u = np.nanmean(list(coal_unc.values()))

    print(f"  ONG δ¹³C: {len(ong_lookup)} countries, Coal δ¹³C: {len(coal_lookup)} countries")

    # Country → 3-box assignment by centroid latitude
    COUNTRY_BOX = {
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
        'iraq': 'Trop', 'mexico': 'Trop',
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

    try:
        import geopandas as gpd
        world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    except Exception:
        world = None

    def assign_box(country_name):
        cn = country_name.lower().strip()
        if cn in COUNTRY_BOX:
            return COUNTRY_BOX[cn]
        if world is not None:
            for _, row in world.iterrows():
                if row['name'].lower() == cn or norm(row['name']) == cn:
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

    def closest_year_col(year, cols):
        yrs = [int(float(c)) for c in cols]
        best = min(yrs, key=lambda y: abs(y - year))
        return str(best)

    def ong_temporal_idx(year):
        return min(max(year - 1970, 0), 52)

    results = {box: np.zeros((N_YEARS, N_MC)) for box in BOX_NAMES}

    for yr_idx, year in enumerate(YEARS):
        tidx = ong_temporal_idx(year)
        ong_yr = closest_year_col(year, ong_year_cols)
        coal_yr = closest_year_col(year, coal_year_cols)

        box_data = {box: {'emis': [], 'vals': [], 'uncs': []} for box in BOX_NAMES}

        for _, row in Country_ONG_emis.iterrows():
            country = norm(row[edgar_country_col_ong])
            emis = row[ong_yr] if ong_yr in row.index else 0
            if pd.isna(emis) or emis <= 0:
                continue
            box = assign_box(country)

            if country in ('united states', 'us', 'usa'):
                base_val = us_ong[tidx]
            elif country == 'china':
                base_val = china_ong[tidx]
            elif country == 'canada':
                base_val = canada_ong[tidx]
            else:
                base_val = ong_lookup.get(country, glob_ong_d13C)

            u = ong_unc.get(country, glob_ong_u)
            box_data[box]['emis'].append(emis)
            box_data[box]['vals'].append(base_val)
            box_data[box]['uncs'].append(u)

        for _, row in Country_Coal_emis.iterrows():
            country = norm(row[edgar_country_col_coal])
            emis = row[coal_yr] if coal_yr in row.index else 0
            if pd.isna(emis) or emis <= 0:
                continue
            box = assign_box(country)
            val = coal_lookup.get(country, glob_coal_d13C)
            u = coal_unc.get(country, glob_coal_u)
            box_data[box]['emis'].append(emis)
            box_data[box]['vals'].append(val)
            box_data[box]['uncs'].append(u)

        for box in BOX_NAMES:
            box_data[box]['emis'] = np.array(box_data[box]['emis'])
            box_data[box]['vals'] = np.array(box_data[box]['vals'])
            box_data[box]['uncs'] = np.array(box_data[box]['uncs'])

        for k in range(N_MC):
            for box in BOX_NAMES:
                d = box_data[box]
                if len(d['emis']) == 0:
                    results[box][yr_idx, k] = glob_ong_d13C
                    continue
                perturbed = d['vals'] + np.random.normal(size=len(d['vals'])) * d['uncs']
                results[box][yr_idx, k] = np.average(perturbed, weights=d['emis'])

        if (yr_idx + 1) % 6 == 0:
            vals = {b: results[b][yr_idx].mean() for b in BOX_NAMES}
            print(f"  FF year {year}: NHext={vals['NHext']:.1f}, "
                  f"Trop={vals['Trop']:.1f}, SHext={vals['SHext']:.1f}")

    return results


# ============================================================================
# 3. Mic δ¹³C — 3-BOX MC
# ============================================================================

def compute_mic_d13C_3box(C4exp, mic_ann, wetland_d13C_grid=None, wetland_d13C_unc_grid=None):
    """Microbial δ¹³C: subcategory mass balance per box.
    Uses isotem spatial wetland δ¹³C if provided."""
    use_isotem = wetland_d13C_grid is not None
    if use_isotem:
        print("Computing Mic δ¹³C 3-box MC (with isotem spatial wetland δ¹³C)...")
    else:
        print("Computing Mic δ¹³C 3-box MC (simplified)...")

    C3_RUM_D13C = -66.64; C3_RUM_STD = 3.39
    C4_RUM_D13C = -54.96; C4_RUM_STD = 3.43
    WASTE_D13C = -54.8; WASTE_STD = 4.4
    RICE_D13C = -59.9; RICE_STD = 4.5
    TERMITE_D13C = -65.2; TERMITE_STD = 7.6

    # Wetland δ¹³C — spatial (isotem) or global (Oh 2022)
    if use_isotem:
        wetland_d13C_ts = None
        wetland_d13C_unc = None
    else:
        oh_path = DATA / "Oh_2022_Wetlands.xlsx"
        if oh_path.exists():
            oh = pd.read_excel(oh_path).values
            wetland_d13C_ts = oh[:, 1]
            wetland_d13C_unc = 0.7
        else:
            wetland_d13C_ts = np.full(N_YEARS, -61.5)
            wetland_d13C_unc = 1.0

    # Ruminant δ¹³C from Chang 2019
    chang_path = DATA / "Chang_2019_ruminants.xlsx"
    if chang_path.exists():
        rum = pd.read_excel(chang_path).values
        rum_d13C_ts = rum[37:61, 1]
        rum_d13C_unc = 1.45
    else:
        rum_d13C_ts = np.full(N_YEARS, -64.8)
        rum_d13C_unc = 1.45

    FRAC_WETLAND = 0.35
    FRAC_RUMINANT = 0.23
    FRAC_RICE = 0.08
    FRAC_WASTE = 0.12
    FRAC_TERMITE = 0.03
    FRAC_WILD_ANIMAL = 0.02
    FRAC_OTHER_WETLAND = 1.0 - (FRAC_WETLAND + FRAC_RUMINANT + FRAC_RICE +
                                  FRAC_WASTE + FRAC_TERMITE + FRAC_WILD_ANIMAL)
    FRAC_WETLAND_TOTAL = FRAC_WETLAND + FRAC_OTHER_WETLAND

    results = {box: np.zeros((N_YEARS, N_MC)) for box in BOX_NAMES}

    for k in range(N_MC):
        c3_rum = C3_RUM_D13C + np.random.normal() * C3_RUM_STD
        c4_rum = C4_RUM_D13C + np.random.normal() * C4_RUM_STD
        waste_d13C = WASTE_D13C + np.random.normal() * WASTE_STD
        rice_d13C = RICE_D13C + np.random.normal() * RICE_STD
        termite_d13C = TERMITE_D13C + np.random.normal() * TERMITE_STD
        if not use_isotem:
            wetland_pert = np.random.normal() * wetland_d13C_unc
        rum_pert = np.random.normal() * rum_d13C_unc
        prop_pert = np.random.normal(loc=1, scale=0.1, size=5)
        suess = SUESS_TREND + np.random.normal() * SUESS_TREND_UNC

        for yr in range(N_YEARS):
            yr_offset = yr - 12
            C4_frac = C4exp[yr] / 100.0
            em = mic_ann[yr]

            for box_name, box_slice in BOX_SLICES.items():
                em_h = em[box_slice, :]
                em_total = em_h.sum()

                if em_total > 0:
                    c4_h = (C4_frac[box_slice, :] * em_h).sum() / em_total
                else:
                    # Fallback C4 fractions by box
                    c4_defaults = {'NHext': 0.15, 'Trop': 0.40, 'SHext': 0.20}
                    c4_h = c4_defaults[box_name]

                rum_d13C_h = (1 - c4_h) * c3_rum + c4_h * c4_rum

                # Wetland δ¹³C — spatial (isotem) or global (Oh 2022)
                if use_isotem:
                    wet_grid = wetland_d13C_grid[yr, box_slice, :]
                    unc_grid = wetland_d13C_unc_grid[yr, box_slice, :]
                    valid = ~np.isnan(wet_grid)
                    if em_total > 0 and valid.any():
                        w = em_h * valid
                        w_sum = w.sum()
                        if w_sum > 0:
                            wet_d13C_h = (wet_grid[valid] * em_h[valid]).sum() / w_sum
                            wet_unc_h = np.sqrt((unc_grid[valid]**2 * em_h[valid]).sum()) / w_sum
                        else:
                            wet_d13C_h = np.nanmean(wet_grid)
                            wet_unc_h = np.nanmean(unc_grid)
                    else:
                        wet_d13C_h = np.nanmean(wet_grid) if valid.any() else -61.0
                        wet_unc_h = np.nanmean(unc_grid) if valid.any() else 1.0
                    wet_d13C_h += np.random.normal() * wet_unc_h
                else:
                    wet_d13C_h = wetland_d13C_ts[yr] + wetland_pert

                waste_s = waste_d13C + yr_offset * suess
                rice_s = rice_d13C + yr_offset * suess
                termite_s = termite_d13C + yr_offset * suess

                fracs = np.array([FRAC_WETLAND_TOTAL, FRAC_RUMINANT + FRAC_WILD_ANIMAL,
                                  FRAC_RICE, FRAC_WASTE, FRAC_TERMITE])
                fracs = fracs * prop_pert
                fracs = np.maximum(fracs, 0)
                fracs /= fracs.sum()

                vals = np.array([wet_d13C_h, rum_d13C_h, rice_s, waste_s, termite_s])
                results[box_name][yr, k] = np.sum(fracs * vals)

        if (k + 1) % 200 == 0:
            print(f"  Mic MC: {k+1}/{N_MC}")

    return results


# ============================================================================
# 4. ATMOSPHERIC δ¹³C — 3-BOX
# ============================================================================

def build_3box_atmospheric_d13C():
    """
    Construct 3-box atmospheric δ¹³C from NH/SH observations.
    
    ch4c13_nh_sh_mean.xlsx has: year_frac, global, NH, SH (fortnightly).
    
    Approximations:
      NHext ≈ NH (extratropics dominate the hemisphere mean in δ¹³C)
      SHext ≈ SH
      Trop  ≈ (NH + SH) / 2  (mixing zone; δ¹³C IH gradient ~0.2‰)
    
    Also produces MC iterations by bootstrap-resampling fortnightly data.
    """
    print("\nBuilding 3-box atmospheric δ¹³C...")

    d13C_file = SCRIPT_DIR / "data" / "ch4c13_nh_sh_mean.xlsx"
    if not d13C_file.exists():
        print(f"  WARNING: {d13C_file} not found")
        return None, None

    df = pd.read_excel(d13C_file, header=None)
    df.columns = ['year_frac', 'global', 'NH', 'SH']
    df['year'] = df['year_frac'].astype(int)

    # Annual means
    annual = df.groupby('year').agg(
        NH=('NH', 'mean'),
        SH=('SH', 'mean'),
        NH_std=('NH', 'std'),
        SH_std=('SH', 'std'),
        n_obs=('NH', 'count'),
    ).reset_index()

    # Filter to our year range
    annual = annual[annual['year'].isin(YEARS)].sort_values('year').reset_index(drop=True)

    # Build annual means
    nhext = annual['NH'].values
    shext = annual['SH'].values
    trop = (annual['NH'].values + annual['SH'].values) / 2.0

    out = pd.DataFrame({
        'year': annual['year'].values.astype(int),
        'NHext': nhext,
        'Trop': trop,
        'SHext': shext,
    })

    out_path = OUT / "ThreeBox_atm_d13C_annual.csv"
    out.to_csv(out_path, index=False, float_format='%.4f')
    print(f"  Saved ThreeBox_atm_d13C_annual.csv ({len(out)} years)")
    print(f"  NHext range: {nhext.min():.3f} to {nhext.max():.3f}‰")
    print(f"  Trop range:  {trop.min():.3f} to {trop.max():.3f}‰")
    print(f"  SHext range: {shext.min():.3f} to {shext.max():.3f}‰")

    # MC iterations via bootstrap of fortnightly observations
    mc_nhext = np.zeros((len(YEARS), N_MC))
    mc_trop = np.zeros((len(YEARS), N_MC))
    mc_shext = np.zeros((len(YEARS), N_MC))

    for yr_idx, year in enumerate(YEARS):
        yr_data = df[df['year'] == year]
        if len(yr_data) == 0:
            mc_nhext[yr_idx, :] = np.nan
            mc_trop[yr_idx, :] = np.nan
            mc_shext[yr_idx, :] = np.nan
            continue

        nh_vals = yr_data['NH'].values
        sh_vals = yr_data['SH'].values

        for k in range(N_MC):
            idx = np.random.choice(len(nh_vals), size=len(nh_vals), replace=True)
            mc_nhext[yr_idx, k] = np.mean(nh_vals[idx])
            mc_shext[yr_idx, k] = np.mean(sh_vals[idx])
            mc_trop[yr_idx, k] = (mc_nhext[yr_idx, k] + mc_shext[yr_idx, k]) / 2.0

    np.savez(OUT / "ThreeBox_atm_d13C_MC.npz",
             years=YEARS, NHext=mc_nhext, Trop=mc_trop, SHext=mc_shext)
    print(f"  Saved ThreeBox_atm_d13C_MC.npz ({mc_nhext.shape})")

    # Also save as CSV (easier to inspect)
    for box_name, mc_data in [('NHext', mc_nhext), ('Trop', mc_trop), ('SHext', mc_shext)]:
        out_mc = np.column_stack([YEARS, mc_data])
        header = "year," + ",".join([f"mc_{i}" for i in range(N_MC)])
        fname = f"ThreeBox_atm_d13C_{box_name}_MC.csv"
        np.savetxt(OUT / fname, out_mc, delimiter=',', header=header,
                   comments='', fmt=['%d'] + ['%.4f'] * N_MC)

    return out, (mc_nhext, mc_trop, mc_shext)


# ============================================================================
# SAVE
# ============================================================================

def save_mc_csv(filename, mc_data):
    """Save as CSV: year + 1000 MC columns."""
    out = np.column_stack([YEARS, mc_data])
    header = "year," + ",".join([f"mc_{i}" for i in range(N_MC)])
    np.savetxt(OUT / filename, out, delimiter=',', header=header,
               comments='', fmt=['%d'] + ['%.4f'] * N_MC)
    print(f"  Saved {filename}: {out.shape}")


def save_summary(mic_res, bb_res, ff_res):
    rows = []
    for yr in range(N_YEARS):
        row = {'year': int(YEARS[yr])}
        for sector, res in [('Mic', mic_res), ('BB', bb_res), ('FF', ff_res)]:
            for box in BOX_NAMES:
                row[f'{sector}_d13C_{box}_mean'] = np.nanmean(res[box][yr])
                row[f'{sector}_d13C_{box}_std'] = np.nanstd(res[box][yr])
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "ThreeBox_d13C_sources_summary.csv", index=False, float_format='%.3f')

    print(f"\n3-Box δ¹³C Source Signature Summary (1998-2021 means):")
    for sector, res in [('Mic', mic_res), ('BB', bb_res), ('FF', ff_res)]:
        means = {box: np.nanmean([np.nanmean(res[box][i]) for i in range(N_YEARS)]) for box in BOX_NAMES}
        print(f"  {sector} δ¹³C:  NHext = {means['NHext']:.2f}  "
              f"Trop = {means['Trop']:.2f}  SHext = {means['SHext']:.2f}‰")


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("Building 3-box δ¹³C source signatures + atmospheric δ¹³C")
    print("  NHext (90°N–30°N) | Trop (30°N–30°S) | SHext (30°S–90°S)")
    print("=" * 70)

    print("\n[1/7] Loading CTCH4 fluxes...")
    mic_ann, fos_ann, pyr_ann = load_CTCH4()
    print(f"  Shapes: mic={mic_ann.shape}, fos={fos_ann.shape}, pyr={pyr_ann.shape}")

    # Show 3-box emission fractions
    for name, flux in [('Mic', mic_ann), ('Fossil', fos_ann), ('Pyrogenic', pyr_ann)]:
        total = flux.sum()
        for box, sl in BOX_SLICES.items():
            frac = flux[:, sl, :].sum() / total
            print(f"  {name} {box}: {frac:.1%}")

    print("\n[2/7] Loading C3/C4 vegetation...")
    C4exp = load_C3C4()
    print(f"  C4exp shape: {C4exp.shape}")

    print("\n[3/7] BB δ¹³C (3-box)...")
    bb_res = compute_bb_d13C_3box(C4exp, pyr_ann)

    print("\n[4/7] FF δ¹³C (3-box)...")
    ff_res = compute_ff_d13C_3box(fos_ann)

    print("\n[5/7] Mic δ¹³C (3-box)...")
    wetland_grid, wetland_unc_grid = load_wetland_d13C()
    mic_res = compute_mic_d13C_3box(C4exp, mic_ann, wetland_grid, wetland_unc_grid)

    print("\n[6/7] Atmospheric δ¹³C (3-box)...")
    atm_result, atm_mc = build_3box_atmospheric_d13C()

    print("\n[7/7] Saving source signature results...")
    for sector, res in [('Mic', mic_res), ('BB', bb_res), ('FF', ff_res)]:
        for box in BOX_NAMES:
            save_mc_csv(f"{sector}_d13C_{box}_MC.csv", res[box])
    save_summary(mic_res, bb_res, ff_res)

    print("\n" + "=" * 70)
    print("Done! 3-box δ¹³C source signatures + atmospheric δ¹³C saved to rel/data/")
    print("=" * 70)


if __name__ == "__main__":
    main()
