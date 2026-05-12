#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_hemispheric_dD_sources.py
===============================
Construct hemispheric (NH/SH) δD source signatures for Mic, BB, and FF
sectors, following Riddell-Young (2025) methodology.

Outputs (saved to rel/output/):
  - Mic_dD_NH_MC.csv   (24 years × 1001 cols: year + 1000 MC)
  - Mic_dD_SH_MC.csv
  - BB_dD_NH_MC.csv    (24 years × 1001 cols: year + 1000 MC)
  - BB_dD_SH_MC.csv
  - FF_dD_NH_MC.csv    (24 years × 1001 cols: year + 1000 MC)  [from EDGAR-based approach]
  - FF_dD_SH_MC.csv
  - Hemispheric_dD_sources_summary.csv  (human-readable summary)

Method:
  Mic: Douglas et al. 2021 regression (δD = 0.6088*MAT - 285.7) on
       d2h_MA.tif, emission-weighted by CTCH4 microbial_flux.
  BB:  Umezawa et al. 2011 regression (δD = 1.16*MAT - 177) on d2h_MA.tif,
       emission-weighted by CTCH4 pyrogenic_flux + GFED5.
  FF:  Country-level ONG + coal δD → EDGAR 8.0 emission-weighted,
       assigned to lat/lon via shapefile → hemisphere split.
       Uses the same approach as FF_dD_map_EDGAR.py but focused on
       hemisphere-level means rather than 1°×1° maps.

All grids are 1°×1° (180×360), row 0 = 90°N, row 179 = 89°S.
NH = rows 0:90 (90°N to 1°N), SH = rows 90:180 (0° to 89°S).
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
    
    # Add missing Arctic rows (77 rows of 5' in the Arctic)
    new_rows = np.full((77, image.shape[1]), SPECIAL, dtype=np.float32)
    image = np.vstack((new_rows, image))
    
    # Downsample from 5' (2160×4320) to 1° (180×360) by block averaging
    block_size = 12
    n_rows = image.shape[0] // block_size  # 180
    n_cols = image.shape[1] // block_size  # 360
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


def fill_ocean_with_lat_avg(grid, special=-3.4e38):
    """Fill ocean/missing cells with latitudinal averages (interpolated)."""
    result = grid.copy()
    mask = result < -1000
    
    masked_data = np.where(mask, np.nan, result)
    row_means = np.nanmean(masked_data, axis=1)
    
    # Interpolate NaN row means from nearest valid rows
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
# HELPER: Cosine-latitude weighting
# ============================================================================

def cos_weights(n_lat=180):
    """Cosine weighting for 1° grid, row 0 = 89.5°N, row 179 = 89.5°S."""
    lats = np.linspace(89.5, -89.5, n_lat)
    w = np.cos(np.radians(lats))
    return w


# ============================================================================
# HELPER: Emission-weighted hemispheric mean δD
# ============================================================================

def hemispheric_emission_weighted_mean(dD_grid, emission_weight_grid):
    """
    Compute NH and SH emission-weighted mean δD.
    
    dD_grid: (180, 360) — δD values per cell
    emission_weight_grid: (180, 360) — fractional emission weight per cell (sums to 1 globally)
    
    Returns: (dD_NH, dD_SH)
    """
    # NH = rows 0:90 (90°N to 1°N), SH = rows 90:180 (0° to 89°S)
    nh_dD = dD_grid[:90, :]
    nh_w = emission_weight_grid[:90, :]
    sh_dD = dD_grid[90:, :]
    sh_w = emission_weight_grid[90:, :]
    
    nh_total = nh_w.sum()
    sh_total = sh_w.sum()
    
    if nh_total > 0:
        dD_NH = (nh_dD * nh_w).sum() / nh_total
    else:
        dD_NH = np.nan
    
    if sh_total > 0:
        dD_SH = (sh_dD * sh_w).sum() / sh_total
    else:
        dD_SH = np.nan
    
    return dD_NH, dD_SH


# ============================================================================
# LOAD CTCH4 EMISSION GRIDS (3°×2° → 1°×1°)
# ============================================================================

def load_CTCH4_fluxes():
    """
    Load CTCH4 microbial, fossil, and pyrogenic monthly fluxes.
    Returns annual emission weights per sector at 1°×1° (180×360).
    Shape: (24, 180, 360) for each sector.
    Row 0 = 90°N (flipped from CTCH4's default lat ordering).
    """
    f = nc.Dataset(str(DATA / "CTCH4_2023_flux3x2.nc"), 'r')
    mic_flux = f.variables['microbial_flux'][:]   # (288, 90, 120) kg/s
    fos_flux = f.variables['fossil_flux'][:]
    pyr_flux = f.variables['pyrogenic_flux'][:]
    f.close()
    
    # Expand from 3°×2° to 1°×1° by repeating
    def expand_to_1deg(flux_3x2):
        """(288, 90, 120) → (288, 180, 360) by repeating cells."""
        expanded = np.repeat(flux_3x2, 2, axis=1)    # lat: 90 → 180
        expanded = np.repeat(expanded, 3, axis=2)     # lon: 120 → 360
        expanded /= 6.0  # conserve total flux
        return expanded
    
    mic_1deg = expand_to_1deg(mic_flux)
    fos_1deg = expand_to_1deg(fos_flux)
    pyr_1deg = expand_to_1deg(pyr_flux)
    
    # Sum to annual (288 months = 24 years × 12 months)
    def annual_sum(monthly):
        return monthly.reshape(N_YEARS, 12, 180, 360).sum(axis=1)
    
    mic_ann = annual_sum(mic_1deg)
    fos_ann = annual_sum(fos_1deg)
    pyr_ann = annual_sum(pyr_1deg)
    
    # Flip latitude so row 0 = 90°N (CTCH4 has row 0 = 90°S)
    mic_ann = mic_ann[:, ::-1, :]
    fos_ann = fos_ann[:, ::-1, :]
    pyr_ann = pyr_ann[:, ::-1, :]
    
    return mic_ann, fos_ann, pyr_ann


# ============================================================================
# LOAD GFED5 BB EMISSIONS (for BB weighting comparison)
# ============================================================================

def load_GFED5_fluxes():
    """
    Load GFED5 Beta monthly CH4 emissions, compress to 1°×1°.
    Returns annual emissions (n_years, 180, 360), covering 2002-2020.
    """
    gfed_dir = DATA / "GFED5_Beta"
    years_gfed = list(range(2002, 2021))
    data_list = []
    for year in years_gfed:
        fn = gfed_dir / f"GFED5_Beta_monthly_{year}.nc"
        if fn.exists():
            with xr.open_dataset(str(fn)) as ds:
                data_list.append(ds['CH4'].values)
    
    if not data_list:
        return None, None
    
    compiled = np.concatenate(data_list, axis=0)  # (N_months, 720, 1440)
    # Compress to 1°×1° by summing 4×4 blocks
    n_months = compiled.shape[0]
    compressed = compiled.reshape(n_months, 180, 4, 360, 4).sum(axis=(2, 4))
    n_years_gfed = n_months // 12
    annual = compressed.reshape(n_years_gfed, 12, 180, 360).sum(axis=1)
    # Flip lat (GFED5 has lat from 90°N to 90°S already? Check convention)
    annual = np.flip(annual, axis=1)  # Ensure row 0 = 90°N
    
    return annual, np.array(years_gfed[:n_years_gfed])


# ============================================================================
# 1. MICROBIAL δD — HEMISPHERIC MC
# ============================================================================

def compute_mic_dD_hemispheric(MAT_1deg, special, mic_ann):
    """
    Compute hemispheric microbial δD MC iterations.
    
    Douglas et al. 2021 regression: δD = slope * MAT - intercept
    Slope = 0.6088 ± 0.072, Intercept = 285.7 ± 6.9
    """
    print("Computing Mic δD hemispheric MC...")
    
    SLOPE = 0.6088
    SLOPE_U = 0.072
    INTERCEPT = 285.7
    INTERCEPT_U = 6.9
    
    # Base map
    dD_mic_base = SLOPE * MAT_1deg - INTERCEPT
    dD_mic_filled = fill_ocean_with_lat_avg(dD_mic_base)
    
    # Annual emission weights
    results_NH = np.zeros((N_YEARS, N_MC))
    results_SH = np.zeros((N_YEARS, N_MC))
    
    for k in range(N_MC):
        # Perturbed regression
        s = SLOPE + np.random.normal() * SLOPE_U
        i = INTERCEPT + np.random.normal() * INTERCEPT_U
        dD_mc = s * MAT_1deg - i
        dD_mc_filled = fill_ocean_with_lat_avg(dD_mc)
        
        for yr_idx in range(N_YEARS):
            # Emission weight for this year
            em = mic_ann[yr_idx]
            em_total = em.sum()
            if em_total <= 0:
                results_NH[yr_idx, k] = np.nan
                results_SH[yr_idx, k] = np.nan
                continue
            em_w = em / em_total
            
            nh, sh = hemispheric_emission_weighted_mean(dD_mc_filled, em_w)
            results_NH[yr_idx, k] = nh
            results_SH[yr_idx, k] = sh
        
        if (k + 1) % 100 == 0:
            print(f"  Mic MC: {k+1}/{N_MC}")
    
    return results_NH, results_SH


# ============================================================================
# 2. BIOMASS BURNING δD — HEMISPHERIC MC
# ============================================================================

def compute_bb_dD_hemispheric(MAT_1deg, special, pyr_ann):
    """
    Compute hemispheric BB δD MC iterations.
    
    Umezawa et al. 2011 regression: δD = 1.16 * MAT - 177
    Slope uncertainty: 0.09, Intercept uncertainty: 6.5
    """
    print("Computing BB δD hemispheric MC...")
    
    SLOPE = 1.16
    SLOPE_U = 0.09
    INTERCEPT = 177.0
    INTERCEPT_U = 6.5
    
    results_NH = np.zeros((N_YEARS, N_MC))
    results_SH = np.zeros((N_YEARS, N_MC))
    
    for k in range(N_MC):
        s = SLOPE + np.random.normal() * SLOPE_U
        i = INTERCEPT + np.random.normal() * INTERCEPT_U
        dD_mc = s * MAT_1deg - i
        dD_mc_filled = fill_ocean_with_lat_avg(dD_mc)
        
        for yr_idx in range(N_YEARS):
            em = pyr_ann[yr_idx]
            em_total = em.sum()
            if em_total <= 0:
                results_NH[yr_idx, k] = np.nan
                results_SH[yr_idx, k] = np.nan
                continue
            em_w = em / em_total
            
            nh, sh = hemispheric_emission_weighted_mean(dD_mc_filled, em_w)
            results_NH[yr_idx, k] = nh
            results_SH[yr_idx, k] = sh
        
        if (k + 1) % 100 == 0:
            print(f"  BB MC: {k+1}/{N_MC}")
    
    return results_NH, results_SH


# ============================================================================
# 3. FOSSIL FUEL δD — HEMISPHERIC MC
# ============================================================================

def compute_ff_dD_hemispheric(fos_ann):
    """
    Compute hemispheric FF δD MC iterations.
    
    Strategy: Use per-country ONG + coal δD values from Riddell-Young's
    compilation, weighted by EDGAR 8.0 emissions. Assign countries to
    NH/SH based on centroid latitude, then compute emission-weighted means.
    
    For MC: perturb each country's δD by its reported uncertainty.
    """
    print("Computing FF δD hemispheric MC...")
    
    # Load country-level data
    coal_dD = pd.read_csv(DATA / "coal_dD.csv")
    ONG_dD = pd.read_csv(DATA / "ONG_dD.csv")
    
    # Load EDGAR emissions
    Country_ONG_emis = pd.read_csv(DATA / "EDGAR8_ONG.csv")
    Country_Coal_emis = pd.read_csv(DATA / "EDGAR8_Coal.csv")
    
    print(f"  Coal δD: {len(coal_dD)} countries, ONG δD: {len(ONG_dD)} countries")
    print(f"  EDGAR ONG: {len(Country_ONG_emis)} countries, Coal: {len(Country_Coal_emis)} countries")
    
    # Inspect column names
    print(f"  coal_dD columns: {list(coal_dD.columns[:5])}")
    print(f"  ONG_dD columns: {list(ONG_dD.columns[:5])}")
    print(f"  EDGAR ONG columns: {list(Country_ONG_emis.columns[:5])}")
    print(f"  EDGAR Coal columns: {list(Country_Coal_emis.columns[:5])}")
    
    # We need to build a simple hemispheric approach:
    # 1. Match countries between δD and EDGAR
    # 2. Use CTCH4 fossil flux grid for NH/SH weighting
    # 3. Apply country δD values weighted by EDGAR emissions to get hemisphere means
    
    # Simpler approach: use CTCH4 fossil_flux for spatial weighting,
    # and compute NH/SH FF δD from the global MC time series with
    # a hemispheric correction factor.
    
    # The key insight: FF sources are overwhelmingly NH-dominated.
    # ONG δD varies by region: US/Canada ≈ -170 to -190, 
    # Middle East ≈ -180, Russia ≈ -200 to -220.
    # Coal δD is generally more depleted: -150 to -200.
    
    # We'll use the gridded approach: compute FF δD for NH and SH separately
    # using the EDGAR emission data and country δD values.
    
    # Build a mapping: for each country, get total ONG + coal emissions
    # and assign to hemisphere based on simple lat classification
    
    # Since FF_dD_map_EDGAR.py uses shapefiles (which we have via geopandas),
    # let's replicate that logic but just for NH/SH aggregation
    
    try:
        import geopandas as gpd
        import shapely.geometry
        
        # Load world boundaries
        world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    except Exception:
        # Fallback: use simple approach based on known country hemispheres
        world = None
    
    # Approach: For each EDGAR year, compute NH and SH emission-weighted FF δD
    
    # First, create country → hemisphere mapping
    country_hemisphere = {}
    if world is not None:
        for _, row in world.iterrows():
            centroid = row.geometry.centroid
            country_hemisphere[row['name']] = 'NH' if centroid.y >= 0 else 'SH'
    
    # Merge ONG and coal data with EDGAR emissions
    # Get country names from EDGAR (first column)
    edgar_country_col_ong = Country_ONG_emis.columns[0]
    edgar_country_col_coal = Country_Coal_emis.columns[0]
    
    # Get year columns (all numeric columns after the country name)
    ong_year_cols = [c for c in Country_ONG_emis.columns[1:] if str(c).replace('.', '').isdigit()]
    coal_year_cols = [c for c in Country_Coal_emis.columns[1:] if str(c).replace('.', '').isdigit()]
    
    print(f"  EDGAR ONG year columns: {ong_year_cols[:3]}...{ong_year_cols[-3:]}")
    print(f"  EDGAR Coal year columns: {coal_year_cols[:3]}...{coal_year_cols[-3:]}")
    
    # Name normalization: EDGAR uses formal names, δD CSVs use common names
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
        """Normalize country name for matching."""
        n = str(name).strip().lower()
        return NAME_MAP.get(n, n)
    
    # Build country δD lookup (using normalized names)
    ong_dD_lookup = {}
    coal_dD_lookup = {}
    ong_dD_unc_lookup = {}
    coal_dD_unc_lookup = {}
    
    for _, row in ONG_dD.iterrows():
        name = normalize_name(row.iloc[0])
        mean_val = row['mean'] if 'mean' in row.index else row.iloc[1]
        unc_val = row['std'] if 'std' in row.index else (row.iloc[2] if len(row) > 2 else 10.0)
        if pd.isna(mean_val):
            continue  # Skip countries without measurements
        ong_dD_lookup[name] = float(mean_val)
        ong_dD_unc_lookup[name] = float(unc_val) if not pd.isna(unc_val) else 10.0
    
    for _, row in coal_dD.iterrows():
        name = normalize_name(row.iloc[0])
        mean_val = row['mean'] if 'mean' in row.index else row.iloc[1]
        unc_val = row['std'] if 'std' in row.index else (row.iloc[2] if len(row) > 2 else 10.0)
        if pd.isna(mean_val):
            continue  # Skip countries without measurements
        coal_dD_lookup[name] = float(mean_val)
        coal_dD_unc_lookup[name] = float(unc_val) if not pd.isna(unc_val) else 10.0
    
    # Global defaults for countries without measurements
    glob_ong_dD = np.nanmean(list(ong_dD_lookup.values()))
    glob_coal_dD = np.nanmean(list(coal_dD_lookup.values()))
    glob_ong_unc = np.nanmean(list(ong_dD_unc_lookup.values()))
    glob_coal_unc = np.nanmean(list(coal_dD_unc_lookup.values()))
    
    print(f"  ONG δD lookup: {len(ong_dD_lookup)} countries with data")
    print(f"  Coal δD lookup: {len(coal_dD_lookup)} countries with data")
    
    print(f"  Global mean ONG δD: {glob_ong_dD:.1f} ± {glob_ong_unc:.1f}‰")
    print(f"  Global mean coal δD: {glob_coal_dD:.1f} ± {glob_coal_unc:.1f}‰")
    
    # For each MC iteration and year: 
    # 1. Perturb each country's δD
    # 2. Compute emission-weighted NH and SH means
    
    # Use EDGAR year closest to our CTCH4 years (1998-2021)
    # EDGAR8 columns are years as strings
    
    results_NH = np.zeros((N_YEARS, N_MC))
    results_SH = np.zeros((N_YEARS, N_MC))
    
    # Pre-assign hemispheres to EDGAR countries
    def assign_hemisphere(country_name):
        """Assign country to NH or SH based on centroid or known geography."""
        cn = country_name.lower().strip()
        
        # Try geopandas world dataset first
        if world is not None:
            for _, row in world.iterrows():
                if row['name'].lower() == cn or normalize_name(row['name']) == cn:
                    return 'NH' if row.geometry.centroid.y >= 0 else 'SH'
        
        # Known SH countries
        sh_countries = {
            'australia', 'new zealand', 'argentina', 'brazil', 'chile',
            'south africa', 'indonesia', 'peru', 'bolivia', 'paraguay',
            'uruguay', 'mozambique', 'madagascar', 'tanzania', 'angola',
            'zambia', 'zimbabwe', 'malawi', 'botswana', 'namibia',
            'congo', 'democratic republic of the congo', 'papua new guinea',
            'east timor', 'fiji', 'solomon islands', 'vanuatu', 'samoa',
            'tonga', 'comoros', 'mauritius', 'seychelles', 'reunion',
            'lesotho', 'eswatini', 'swaziland', 'burundi', 'rwanda',
            'equatorial guinea', 'gabon', 'sao tome and principe',
            'ecuador'  # straddles equator but mostly SH
        }
        if country_name.lower() in sh_countries:
            return 'SH'
        return 'NH'  # Default to NH (most major emitters are NH)
    
    # Pre-compute: for each CTCH4 year, find the closest EDGAR year column
    def get_edgar_year_col(year, year_cols):
        """Get closest EDGAR year column to target year."""
        year_vals = [int(float(c)) for c in year_cols]
        closest = min(year_vals, key=lambda y: abs(y - year))
        return str(closest)
    
    # Compute for each year
    for yr_idx, year in enumerate(YEARS_CTCH4):
        # Get EDGAR emissions for this year
        ong_yr_col = get_edgar_year_col(year, ong_year_cols)
        coal_yr_col = get_edgar_year_col(year, coal_year_cols)
        
        # Build per-hemisphere emission and δD arrays
        nh_ong_emis = []
        nh_ong_dD_vals = []
        nh_ong_dD_uncs = []
        sh_ong_emis = []
        sh_ong_dD_vals = []
        sh_ong_dD_uncs = []
        
        nh_coal_emis = []
        nh_coal_dD_vals = []
        nh_coal_dD_uncs = []
        sh_coal_emis = []
        sh_coal_dD_vals = []
        sh_coal_dD_uncs = []
        
        for _, row in Country_ONG_emis.iterrows():
            country = normalize_name(row[edgar_country_col_ong])
            emis = row[ong_yr_col] if ong_yr_col in row.index else 0
            if pd.isna(emis) or emis <= 0:
                continue
            
            hemi = assign_hemisphere(country)
            dD_val = ong_dD_lookup.get(country, glob_ong_dD)
            dD_unc = ong_dD_unc_lookup.get(country, glob_ong_unc)
            
            if hemi == 'NH':
                nh_ong_emis.append(emis)
                nh_ong_dD_vals.append(dD_val)
                nh_ong_dD_uncs.append(dD_unc)
            else:
                sh_ong_emis.append(emis)
                sh_ong_dD_vals.append(dD_val)
                sh_ong_dD_uncs.append(dD_unc)
        
        for _, row in Country_Coal_emis.iterrows():
            country = normalize_name(row[edgar_country_col_coal])
            emis = row[coal_yr_col] if coal_yr_col in row.index else 0
            if pd.isna(emis) or emis <= 0:
                continue
            
            hemi = assign_hemisphere(country)
            dD_val = coal_dD_lookup.get(country, glob_coal_dD)
            dD_unc = coal_dD_unc_lookup.get(country, glob_coal_unc)
            
            if hemi == 'NH':
                nh_coal_emis.append(emis)
                nh_coal_dD_vals.append(dD_val)
                nh_coal_dD_uncs.append(dD_unc)
            else:
                sh_coal_emis.append(emis)
                sh_coal_dD_vals.append(dD_val)
                sh_coal_dD_uncs.append(dD_unc)
        
        # Convert to arrays
        nh_ong_emis = np.array(nh_ong_emis)
        nh_ong_dD_vals = np.array(nh_ong_dD_vals)
        nh_ong_dD_uncs = np.array(nh_ong_dD_uncs)
        sh_ong_emis = np.array(sh_ong_emis)
        sh_ong_dD_vals = np.array(sh_ong_dD_vals)
        sh_ong_dD_uncs = np.array(sh_ong_dD_uncs)
        
        nh_coal_emis = np.array(nh_coal_emis)
        nh_coal_dD_vals = np.array(nh_coal_dD_vals)
        nh_coal_dD_uncs = np.array(nh_coal_dD_uncs)
        sh_coal_emis = np.array(sh_coal_emis)
        sh_coal_dD_vals = np.array(sh_coal_dD_vals)
        sh_coal_dD_uncs = np.array(sh_coal_dD_uncs)
        
        # MC iterations
        for k in range(N_MC):
            # Perturb δD for each country
            nh_ong_dD_perturbed = nh_ong_dD_vals + np.random.normal(size=len(nh_ong_dD_vals)) * nh_ong_dD_uncs
            sh_ong_dD_perturbed = sh_ong_dD_vals + np.random.normal(size=len(sh_ong_dD_vals)) * sh_ong_dD_uncs
            nh_coal_dD_perturbed = nh_coal_dD_vals + np.random.normal(size=len(nh_coal_dD_vals)) * nh_coal_dD_uncs
            sh_coal_dD_perturbed = sh_coal_dD_vals + np.random.normal(size=len(sh_coal_dD_vals)) * sh_coal_dD_uncs
            
            # Emission-weighted mean per hemisphere
            nh_all_emis = np.concatenate([nh_ong_emis, nh_coal_emis])
            nh_all_dD = np.concatenate([nh_ong_dD_perturbed, nh_coal_dD_perturbed])
            sh_all_emis = np.concatenate([sh_ong_emis, sh_coal_emis])
            sh_all_dD = np.concatenate([sh_ong_dD_perturbed, sh_coal_dD_perturbed])
            
            if nh_all_emis.sum() > 0:
                results_NH[yr_idx, k] = np.average(nh_all_dD, weights=nh_all_emis)
            else:
                results_NH[yr_idx, k] = np.nan
            
            if sh_all_emis.sum() > 0:
                results_SH[yr_idx, k] = np.average(sh_all_dD, weights=sh_all_emis)
            else:
                results_SH[yr_idx, k] = np.nan
        
        if (yr_idx + 1) % 6 == 0:
            print(f"  FF year {year}: NH={results_NH[yr_idx].mean():.1f}‰, SH={results_SH[yr_idx].mean():.1f}‰")
    
    return results_NH, results_SH


# ============================================================================
# SAVE RESULTS
# ============================================================================

def save_mc_csv(filename, years, mc_data):
    """Save MC results as CSV: year + 1000 MC columns."""
    out_data = np.column_stack([years, mc_data])
    header = "year," + ",".join([f"mc_{i}" for i in range(N_MC)])
    np.savetxt(OUT / filename, out_data, delimiter=',', header=header, 
               comments='', fmt=['%d'] + ['%.3f'] * N_MC)
    print(f"  Saved {filename}: {out_data.shape}")


def save_summary(mic_NH, mic_SH, bb_NH, bb_SH, ff_NH, ff_SH):
    """Save human-readable summary CSV."""
    rows = []
    for yr_idx, year in enumerate(YEARS_CTCH4):
        rows.append({
            'year': int(year),
            'Mic_dD_NH_mean': np.nanmean(mic_NH[yr_idx]),
            'Mic_dD_NH_std': np.nanstd(mic_NH[yr_idx]),
            'Mic_dD_SH_mean': np.nanmean(mic_SH[yr_idx]),
            'Mic_dD_SH_std': np.nanstd(mic_SH[yr_idx]),
            'BB_dD_NH_mean': np.nanmean(bb_NH[yr_idx]),
            'BB_dD_NH_std': np.nanstd(bb_NH[yr_idx]),
            'BB_dD_SH_mean': np.nanmean(bb_SH[yr_idx]),
            'BB_dD_SH_std': np.nanstd(bb_SH[yr_idx]),
            'FF_dD_NH_mean': np.nanmean(ff_NH[yr_idx]),
            'FF_dD_NH_std': np.nanstd(ff_NH[yr_idx]),
            'FF_dD_SH_mean': np.nanmean(ff_SH[yr_idx]),
            'FF_dD_SH_std': np.nanstd(ff_SH[yr_idx]),
        })
    
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "Hemispheric_dD_sources_summary.csv", index=False, float_format='%.2f')
    print(f"\nSummary saved. Mean values across 1998-2021:")
    print(f"  Mic δD:  NH = {df['Mic_dD_NH_mean'].mean():.1f} ± {df['Mic_dD_NH_std'].mean():.1f}‰  "
          f"SH = {df['Mic_dD_SH_mean'].mean():.1f} ± {df['Mic_dD_SH_std'].mean():.1f}‰  "
          f"Δ(NH-SH) = {(df['Mic_dD_NH_mean'] - df['Mic_dD_SH_mean']).mean():.1f}‰")
    print(f"  BB δD:   NH = {df['BB_dD_NH_mean'].mean():.1f} ± {df['BB_dD_NH_std'].mean():.1f}‰  "
          f"SH = {df['BB_dD_SH_mean'].mean():.1f} ± {df['BB_dD_SH_std'].mean():.1f}‰  "
          f"Δ(NH-SH) = {(df['BB_dD_NH_mean'] - df['BB_dD_SH_mean']).mean():.1f}‰")
    print(f"  FF δD:   NH = {df['FF_dD_NH_mean'].mean():.1f} ± {df['FF_dD_NH_std'].mean():.1f}‰  "
          f"SH = {df['FF_dD_SH_mean'].mean():.1f} ± {df['FF_dD_SH_std'].mean():.1f}‰  "
          f"Δ(NH-SH) = {(df['FF_dD_NH_mean'] - df['FF_dD_SH_mean']).mean():.1f}‰")


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("Building hemispheric δD source signatures")
    print("=" * 70)
    
    # Step 1: Load MAT grid
    print("\n[1/5] Loading MAT TIFF and downsampling to 1°×1°...")
    MAT_1deg, SPECIAL = load_MAT_grid()
    print(f"  MAT grid shape: {MAT_1deg.shape}")
    print(f"  Valid cells: {(MAT_1deg > -1000).sum()} / {MAT_1deg.size}")
    
    # Step 2: Load CTCH4 fluxes
    print("\n[2/5] Loading CTCH4 emission grids...")
    mic_ann, fos_ann, pyr_ann = load_CTCH4_fluxes()
    print(f"  Microbial annual flux shape: {mic_ann.shape}")
    print(f"  Fossil annual flux shape: {fos_ann.shape}")
    print(f"  Pyrogenic annual flux shape: {pyr_ann.shape}")
    
    # Check NH/SH emission fractions
    for name, flux in [('Mic', mic_ann), ('Fossil', fos_ann), ('Pyrogenic', pyr_ann)]:
        nh_frac = flux[:, :90, :].sum() / flux.sum()
        print(f"  {name} NH fraction: {nh_frac:.1%}")
    
    # Step 3: Mic δD
    print("\n[3/5] Microbial δD...")
    mic_NH, mic_SH = compute_mic_dD_hemispheric(MAT_1deg, SPECIAL, mic_ann)
    
    # Step 4: BB δD
    print("\n[4/5] Biomass Burning δD...")
    bb_NH, bb_SH = compute_bb_dD_hemispheric(MAT_1deg, SPECIAL, pyr_ann)
    
    # Step 5: FF δD
    print("\n[5/5] Fossil Fuel δD...")
    ff_NH, ff_SH = compute_ff_dD_hemispheric(fos_ann)
    
    # Save all results
    print("\nSaving results...")
    save_mc_csv("Mic_dD_NH_MC.csv", YEARS_CTCH4, mic_NH)
    save_mc_csv("Mic_dD_SH_MC.csv", YEARS_CTCH4, mic_SH)
    save_mc_csv("BB_dD_NH_MC.csv", YEARS_CTCH4, bb_NH)
    save_mc_csv("BB_dD_SH_MC.csv", YEARS_CTCH4, bb_SH)
    save_mc_csv("FF_dD_NH_MC.csv", YEARS_CTCH4, ff_NH)
    save_mc_csv("FF_dD_SH_MC.csv", YEARS_CTCH4, ff_SH)
    save_summary(mic_NH, mic_SH, bb_NH, bb_SH, ff_NH, ff_SH)
    
    print("\n" + "=" * 70)
    print("Done! All hemispheric δD source signatures saved to rel/output/")
    print("=" * 70)


if __name__ == "__main__":
    main()
