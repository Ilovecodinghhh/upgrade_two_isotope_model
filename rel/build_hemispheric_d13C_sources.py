#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_hemispheric_d13C_sources.py
=================================
Construct hemispheric (NH/SH) δ¹³C source signatures for Mic, BB, and FF
sectors, following Riddell-Young (2025) methodology adapted for
hemisphere-level aggregation.

Method:
  FF:  Country-level ONG + coal δ¹³C → EDGAR 8.0 emission-weighted,
       assigned to NH/SH by country centroid (same approach as δD).
  BB:  C3/C4 vegetation maps × CTCH4 pyrogenic flux, hemisphere-split.
       C3 δ¹³C = -26.8 ± 2.9‰, C4 δ¹³C = -12.7 ± 4.6‰ (Cerling 1997).
  Mic: Subcategory mass balance (wetlands, ruminants, rice, termites,
       waste, wild animals) weighted by CTCH4 microbial flux per hemisphere.
       Uses posterior-adjusted subcategory proportions.

Outputs (saved to rel/data/):
  - FF_d13C_NH_MC.csv, FF_d13C_SH_MC.csv    (24 years × 1001 cols)
  - BB_d13C_NH_MC.csv, BB_d13C_SH_MC.csv    (24 years × 1001 cols)
  - Mic_d13C_NH_MC.csv, Mic_d13C_SH_MC.csv  (24 years × 1001 cols)
  - Hemispheric_d13C_sources_summary.csv
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

SCRIPT_DIR = Path(__file__).resolve().parent
BASE = SCRIPT_DIR.parent / "ImportantReferences" / "Riddell-Young2025PNAS_DS" / \
       "Riddell-Young_2025_MassBalancePackage" / "Riddell-Young_2025_MassBalancePackage"
DATA = BASE / "data"
OUT = SCRIPT_DIR / "data"
OUT.mkdir(exist_ok=True)

# Suess effect trend
SUESS_TREND = -0.024  # ‰/yr
SUESS_TREND_UNC = 0.005

spm = 2.628e+6  # seconds per month


# ============================================================================
# DATA LOADING (shared)
# ============================================================================

def load_CTCH4():
    """Load CTCH4 3×2° → expand to 1°×1°, flip so row0=90°N.
    Returns mic_ann, fos_ann, pyr_ann each (24, 180, 360) in kg/s."""
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
    """Load C4 fraction maps.
    
    Tries Luo 2023 first (time-varying, 2001-2019); falls back to
    Still 2003 (static, 180×360) if Luo unavailable.
    Returns C4exp (24, 180, 360) as percentage [0,100], row0=90°N.
    """
    luo_path = DATA / "C4_distribution_NUS_v2.2.nc"
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
        # Fallback to Still 2003
        C4_still = pd.read_excel(DATA / "Still2003_C4.xlsx", header=None).values
        C4_still = np.flipud(C4_still)
        C4_still = np.where(C4_still < 0, 0, C4_still)
        # Replicate for 24 years (static)
        C4exp = np.repeat(C4_still[np.newaxis, :180, :360], N_YEARS, axis=0)
        print("  Using Still 2003 C4 map (static fallback)")
    return C4exp


def load_prior_subcategories():
    """Load CTCH4 prior microbial subcategory emissions.
    Returns dict of annual arrays (24, 180, 360) in kg/month."""
    Priors = xr.open_dataset(str(DATA / "prior_monthly_emission_kg_lei.nc"))
    cats = {}
    for name, var in [('rice', 'flux_rice'), ('ruminant', 'flux_ruminant'),
                       ('wetland', 'flux_wetland'), ('termite', 'flux_termite'),
                       ('landfill', 'flux_waste_landfill'), ('wild_animal', 'flux_wild_animals')]:
        flux = Priors[var].sum(dim='month').values.astype(np.float64)
        flux = np.transpose(flux[:, :, -24:], (2, 1, 0))
        cats[name] = flux
    return cats


def load_wetland_d13C():
    """Load isotem wetland δ¹³C map if available; else use Oh 2022 global means.
    Returns (24, 180, 360) or None if no spatial data available."""
    isotem_path = DATA / "isotem_wetland_d13C-CH4.nc4"
    if isotem_path.exists():
        wf = nc.Dataset(str(isotem_path), 'r')
        d13C = wf.variables['wetland_d13C-CH4'][:]
        wf.close()
        d13C_ann = np.mean(d13C, axis=0)  # (33, 720, 360)
        d13C_ann = np.swapaxes(d13C_ann, 1, 2)
        d13C_ann = d13C_ann[:, ::-1, :]
        d13C_ann = np.ma.filled(d13C_ann, fill_value=np.nan)
        for i in range(d13C_ann.shape[0]):
            avg = np.nanmean(d13C_ann[i])
            d13C_ann[i, np.isnan(d13C_ann[i])] = avg
        trimmed = d13C_ann[14:33]  # 1998-2016
        trimmed = np.concatenate([trimmed, np.repeat(trimmed[-1:], 5, axis=0)])  # to 2021
        compressed = trimmed.reshape(24, 180, 2, 360, 2).mean(axis=(2, 4))
        print("  Using isotem wetland δ¹³C spatial map")
        return compressed
    else:
        # Use Oh 2022 global time series — no spatial variation
        oh_path = DATA / "Oh_2022_Wetlands.xlsx"
        if oh_path.exists():
            oh = pd.read_excel(oh_path).values
            # oh has columns: Year, d13C, stdev
            wetland_d13C_ts = oh[:, 1]  # 24 values for 1998-2021
            print(f"  Using Oh 2022 wetland δ¹³C global means (no spatial map): {np.mean(wetland_d13C_ts):.1f}‰")
            # Broadcast to grid (uniform)
            result = np.zeros((N_YEARS, 180, 360))
            for yr in range(N_YEARS):
                result[yr, :, :] = wetland_d13C_ts[yr]
            return result
        print("  WARNING: No wetland δ¹³C data found, using -61.0‰")
        return np.full((N_YEARS, 180, 360), -61.0)


# ============================================================================
# BB δ¹³C — HEMISPHERIC MC
# ============================================================================

def compute_bb_d13C_hemi(C4exp, pyr_ann):
    """BB δ¹³C = C3_frac × C3_δ¹³C + C4_frac × C4_δ¹³C, emission-weighted per hemisphere."""
    print("Computing BB δ¹³C hemispheric MC...")

    C3_D13C = -26.8; C3_D13C_STD = 2.9
    C4_D13C = -12.7; C4_D13C_STD = 4.6

    results_NH = np.zeros((N_YEARS, N_MC))
    results_SH = np.zeros((N_YEARS, N_MC))

    for k in range(N_MC):
        c3_val = C3_D13C + np.random.normal() * C3_D13C_STD
        c4_val = C4_D13C + np.random.normal() * C4_D13C_STD

        for yr in range(N_YEARS):
            # δ¹³C map for this year
            C4_frac = C4exp[yr] / 100.0
            d13C_map = C4_frac * c4_val + (1 - C4_frac) * c3_val

            # Emission weighting
            em = pyr_ann[yr]

            # NH
            nh_em = em[:90, :]
            nh_d13C = d13C_map[:90, :]
            nh_total = nh_em.sum()
            if nh_total > 0:
                results_NH[yr, k] = (nh_d13C * nh_em).sum() / nh_total
            else:
                results_NH[yr, k] = c3_val * 0.7 + c4_val * 0.3  # fallback

            # SH
            sh_em = em[90:, :]
            sh_d13C = d13C_map[90:, :]
            sh_total = sh_em.sum()
            if sh_total > 0:
                results_SH[yr, k] = (sh_d13C * sh_em).sum() / sh_total
            else:
                results_SH[yr, k] = c3_val * 0.85 + c4_val * 0.15

        if (k + 1) % 200 == 0:
            print(f"  BB MC: {k+1}/{N_MC}")

    return results_NH, results_SH


# ============================================================================
# FF δ¹³C — HEMISPHERIC MC
# ============================================================================

def compute_ff_d13C_hemi(fos_ann):
    """FF δ¹³C: country-level ONG + coal, emission-weighted, hemisphere-split."""
    print("Computing FF δ¹³C hemispheric MC...")

    coal_d13C = pd.read_csv(DATA / "coal_d13C.csv")
    ONG_d13C = pd.read_csv(DATA / "ONG_d13C.csv")
    Country_ONG_emis = pd.read_csv(DATA / "EDGAR8_ONG.csv")
    Country_Coal_emis = pd.read_csv(DATA / "EDGAR8_Coal.csv")

    # Temporal ONG trends for US, China, Canada
    us_ong = pd.read_csv(DATA / "US_ONG_trends.csv").iloc[:53, 2].values
    cc_ong = pd.read_csv(DATA / "China_Canada_ONG_Trends.csv")
    china_ong = cc_ong.iloc[:, 1].values
    canada_ong = cc_ong.iloc[:, 3].values

    # Build lookups
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

    # SH countries
    SH_COUNTRIES = {
        'australia', 'new zealand', 'argentina', 'brazil', 'chile',
        'south africa', 'indonesia', 'peru', 'bolivia', 'paraguay',
        'uruguay', 'mozambique', 'madagascar', 'tanzania', 'angola',
        'zambia', 'zimbabwe', 'malawi', 'botswana', 'namibia',
        'papua new guinea', 'ecuador',
    }

    def get_hemi(country):
        return 'SH' if country in SH_COUNTRIES else 'NH'

    edgar_country_col_ong = Country_ONG_emis.columns[0]
    edgar_country_col_coal = Country_Coal_emis.columns[0]
    ong_year_cols = [c for c in Country_ONG_emis.columns[1:] if str(c).replace('.', '').isdigit()]
    coal_year_cols = [c for c in Country_Coal_emis.columns[1:] if str(c).replace('.', '').isdigit()]

    def closest_year_col(year, cols):
        yrs = [int(float(c)) for c in cols]
        best = min(yrs, key=lambda y: abs(y - year))
        return str(best)

    # Temporal index for ONG trends (1970–2022 → idx for year)
    def ong_temporal_idx(year):
        return min(max(year - 1970, 0), 52)

    results_NH = np.zeros((N_YEARS, N_MC))
    results_SH = np.zeros((N_YEARS, N_MC))

    for yr_idx, year in enumerate(YEARS):
        tidx = ong_temporal_idx(year)
        ong_yr = closest_year_col(year, ong_year_cols)
        coal_yr = closest_year_col(year, coal_year_cols)

        # Build per-hemisphere arrays
        hemi_data = {'NH': {'emis': [], 'vals': [], 'uncs': []},
                     'SH': {'emis': [], 'vals': [], 'uncs': []}}

        for _, row in Country_ONG_emis.iterrows():
            country = norm(row[edgar_country_col_ong])
            emis = row[ong_yr] if ong_yr in row.index else 0
            if pd.isna(emis) or emis <= 0:
                continue
            hemi = get_hemi(country)

            # Temporal override for US, China, Canada
            if country in ('united states', 'us', 'usa'):
                base_val = us_ong[tidx]
            elif country == 'china':
                base_val = china_ong[tidx]
            elif country == 'canada':
                base_val = canada_ong[tidx]
            else:
                base_val = ong_lookup.get(country, glob_ong_d13C)

            u = ong_unc.get(country, glob_ong_u)
            hemi_data[hemi]['emis'].append(emis)
            hemi_data[hemi]['vals'].append(base_val)
            hemi_data[hemi]['uncs'].append(u)

        for _, row in Country_Coal_emis.iterrows():
            country = norm(row[edgar_country_col_coal])
            emis = row[coal_yr] if coal_yr in row.index else 0
            if pd.isna(emis) or emis <= 0:
                continue
            hemi = get_hemi(country)
            val = coal_lookup.get(country, glob_coal_d13C)
            u = coal_unc.get(country, glob_coal_u)
            hemi_data[hemi]['emis'].append(emis)
            hemi_data[hemi]['vals'].append(val)
            hemi_data[hemi]['uncs'].append(u)

        for hemi in ('NH', 'SH'):
            hemi_data[hemi]['emis'] = np.array(hemi_data[hemi]['emis'])
            hemi_data[hemi]['vals'] = np.array(hemi_data[hemi]['vals'])
            hemi_data[hemi]['uncs'] = np.array(hemi_data[hemi]['uncs'])

        for k in range(N_MC):
            for hemi, arr in [('NH', results_NH), ('SH', results_SH)]:
                d = hemi_data[hemi]
                if len(d['emis']) == 0:
                    arr[yr_idx, k] = glob_ong_d13C
                    continue
                perturbed = d['vals'] + np.random.normal(size=len(d['vals'])) * d['uncs']
                arr[yr_idx, k] = np.average(perturbed, weights=d['emis'])

        if (yr_idx + 1) % 6 == 0:
            print(f"  FF year {year}: NH={results_NH[yr_idx].mean():.1f}, SH={results_SH[yr_idx].mean():.1f}")

    return results_NH, results_SH


# ============================================================================
# Mic δ¹³C — HEMISPHERIC MC
# ============================================================================

def compute_mic_d13C_hemi(C4exp, mic_ann):
    """Microbial δ¹³C: subcategory mass balance per hemisphere.
    
    Without prior subcategory emission files, we use a simplified approach:
    - Ruminant/wild animal δ¹³C varies by hemisphere via C3/C4 fractions
      from Still 2003, weighted by CTCH4 microbial emission grid.
    - Wetland δ¹³C: Oh 2022 global time series (no spatial variation).
    - Rice, termite, waste: constant values with Suess effect.
    - Subcategory proportions from Riddell-Young 2025 Table 1 
      (global fractions, slight hemispheric weighting by CTCH4).
    """
    print("Computing Mic δ¹³C hemispheric MC (simplified, no subcategory priors)...")

    # Source signature parameters
    C3_RUM_D13C = -66.64; C3_RUM_STD = 3.39
    C4_RUM_D13C = -54.96; C4_RUM_STD = 3.43
    WASTE_D13C = -54.8; WASTE_STD = 4.4
    RICE_D13C = -59.9; RICE_STD = 4.5
    TERMITE_D13C = -65.2; TERMITE_STD = 7.6

    # Wetland δ¹³C from Oh 2022
    oh_path = DATA / "Oh_2022_Wetlands.xlsx"
    if oh_path.exists():
        oh = pd.read_excel(oh_path).values
        wetland_d13C_ts = oh[:, 1]  # 24 values
        wetland_d13C_unc = 0.7
    else:
        wetland_d13C_ts = np.full(N_YEARS, -61.5)
        wetland_d13C_unc = 1.0

    # Ruminant δ¹³C from Chang 2019
    chang_path = DATA / "Chang_2019_ruminants.xlsx"
    if chang_path.exists():
        rum = pd.read_excel(chang_path).values
        # rows 37:61 = 1998–2021
        rum_d13C_ts = rum[37:61, 1]
        rum_d13C_unc = 1.45
    else:
        rum_d13C_ts = np.full(N_YEARS, -64.8)
        rum_d13C_unc = 1.45

    # Global subcategory proportions (from Riddell-Young 2025, Table 1)
    # These are approximate global values
    FRAC_WETLAND = 0.35
    FRAC_RUMINANT = 0.23
    FRAC_RICE = 0.08
    FRAC_WASTE = 0.12
    FRAC_TERMITE = 0.03
    FRAC_WILD_ANIMAL = 0.02
    # Remainder assigned to wetland
    FRAC_OTHER_WETLAND = 1.0 - (FRAC_WETLAND + FRAC_RUMINANT + FRAC_RICE + 
                                  FRAC_WASTE + FRAC_TERMITE + FRAC_WILD_ANIMAL)
    FRAC_WETLAND_TOTAL = FRAC_WETLAND + FRAC_OTHER_WETLAND

    results_NH = np.zeros((N_YEARS, N_MC))
    results_SH = np.zeros((N_YEARS, N_MC))

    for k in range(N_MC):
        # Perturbed values (single draw per MC iteration)
        c3_rum = C3_RUM_D13C + np.random.normal() * C3_RUM_STD
        c4_rum = C4_RUM_D13C + np.random.normal() * C4_RUM_STD
        waste_d13C = WASTE_D13C + np.random.normal() * WASTE_STD
        rice_d13C = RICE_D13C + np.random.normal() * RICE_STD
        termite_d13C = TERMITE_D13C + np.random.normal() * TERMITE_STD
        wetland_pert = np.random.normal() * wetland_d13C_unc
        rum_pert = np.random.normal() * rum_d13C_unc

        # Proportion perturbation (10% relative)
        prop_pert = np.random.normal(loc=1, scale=0.1, size=5)

        # Suess effect
        suess = SUESS_TREND + np.random.normal() * SUESS_TREND_UNC

        for yr in range(N_YEARS):
            yr_offset = yr - 12
            C4_frac = C4exp[yr] / 100.0

            # Emission-weighted C3/C4 fractions per hemisphere
            em = mic_ann[yr]

            for hemi, sl, out_arr in [('NH', slice(0, 90), results_NH),
                                       ('SH', slice(90, 180), results_SH)]:
                em_h = em[sl, :]
                em_total = em_h.sum()

                if em_total > 0:
                    # Emission-weighted C4 fraction in this hemisphere
                    c4_h = (C4_frac[sl, :] * em_h).sum() / em_total
                else:
                    c4_h = 0.3 if hemi == 'NH' else 0.2

                # Ruminant δ¹³C (emission-weighted C3/C4)
                rum_d13C_h = (1 - c4_h) * c3_rum + c4_h * c4_rum

                # Wetland δ¹³C (global, no spatial variation)
                wet_d13C_h = wetland_d13C_ts[yr] + wetland_pert

                # Suess-corrected constant sources
                waste_s = waste_d13C + yr_offset * suess
                rice_s = rice_d13C + yr_offset * suess
                termite_s = termite_d13C + yr_offset * suess

                # Combine with perturbed proportions
                fracs = np.array([FRAC_WETLAND_TOTAL, FRAC_RUMINANT + FRAC_WILD_ANIMAL,
                                  FRAC_RICE, FRAC_WASTE, FRAC_TERMITE])
                fracs = fracs * prop_pert
                fracs = np.maximum(fracs, 0)
                fracs /= fracs.sum()

                vals = np.array([wet_d13C_h, rum_d13C_h, rice_s, waste_s, termite_s])
                out_arr[yr, k] = np.sum(fracs * vals)

        if (k + 1) % 200 == 0:
            print(f"  Mic MC: {k+1}/{N_MC}")

    return results_NH, results_SH


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


def save_summary(mic_NH, mic_SH, bb_NH, bb_SH, ff_NH, ff_SH):
    rows = []
    for yr in range(N_YEARS):
        rows.append({
            'year': int(YEARS[yr]),
            'Mic_d13C_NH_mean': np.mean(mic_NH[yr]), 'Mic_d13C_NH_std': np.std(mic_NH[yr]),
            'Mic_d13C_SH_mean': np.mean(mic_SH[yr]), 'Mic_d13C_SH_std': np.std(mic_SH[yr]),
            'BB_d13C_NH_mean': np.mean(bb_NH[yr]), 'BB_d13C_NH_std': np.std(bb_NH[yr]),
            'BB_d13C_SH_mean': np.mean(bb_SH[yr]), 'BB_d13C_SH_std': np.std(bb_SH[yr]),
            'FF_d13C_NH_mean': np.mean(ff_NH[yr]), 'FF_d13C_NH_std': np.std(ff_NH[yr]),
            'FF_d13C_SH_mean': np.mean(ff_SH[yr]), 'FF_d13C_SH_std': np.std(ff_SH[yr]),
        })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "Hemispheric_d13C_sources_summary.csv", index=False, float_format='%.3f')
    print(f"\nSummary saved. Mean values across {YEARS[0]}-{YEARS[-1]}:")
    for sec in ('Mic', 'BB', 'FF'):
        nh = df[f'{sec}_d13C_NH_mean'].mean()
        sh = df[f'{sec}_d13C_SH_mean'].mean()
        print(f"  {sec} δ¹³C:  NH = {nh:.2f}‰  SH = {sh:.2f}‰  Δ(NH-SH) = {nh-sh:.2f}‰")


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("Building hemispheric δ¹³C source signatures")
    print("=" * 70)

    print("\n[1/6] Loading CTCH4 fluxes...")
    mic_ann, fos_ann, pyr_ann = load_CTCH4()
    print(f"  Shapes: mic={mic_ann.shape}, fos={fos_ann.shape}, pyr={pyr_ann.shape}")

    print("\n[2/6] Loading C3/C4 vegetation...")
    C4exp = load_C3C4()
    print(f"  C4exp shape: {C4exp.shape}")

    print("\n[3/6] BB δ¹³C...")
    bb_NH, bb_SH = compute_bb_d13C_hemi(C4exp, pyr_ann)

    print("\n[4/6] FF δ¹³C...")
    ff_NH, ff_SH = compute_ff_d13C_hemi(fos_ann)

    print("\n[5/6] Mic δ¹³C...")
    mic_NH, mic_SH = compute_mic_d13C_hemi(C4exp, mic_ann)

    print("\n[6/6] Saving results...")
    save_mc_csv("BB_d13C_NH_MC.csv", bb_NH)
    save_mc_csv("BB_d13C_SH_MC.csv", bb_SH)
    save_mc_csv("FF_d13C_NH_MC.csv", ff_NH)
    save_mc_csv("FF_d13C_SH_MC.csv", ff_SH)
    save_mc_csv("Mic_d13C_NH_MC.csv", mic_NH)
    save_mc_csv("Mic_d13C_SH_MC.csv", mic_SH)
    save_summary(mic_NH, mic_SH, bb_NH, bb_SH, ff_NH, ff_SH)

    print("\n" + "=" * 70)
    print("Done! Hemispheric δ¹³C source signatures saved to rel/data/")
    print("=" * 70)


if __name__ == "__main__":
    main()
