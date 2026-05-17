# 3-Box δ¹³C Data Preparation

## Overview

The 3-box model divides the atmosphere into three latitude bands:
- **NHext** (NH extratropical): 90°N – 30°N (grid rows 0–59)
- **Trop** (Tropical): 30°N – 30°S (grid rows 60–119)
- **SHext** (SH extratropical): 30°S – 90°S (grid rows 120–179)

Two build scripts produce all required per-box δ¹³C data:
- `rel/build_3box_d13C_sources.py` — source signatures + atmospheric δ¹³C
- `rel/build_hemispheric_d13C_sources.py` — 2-box (NH/SH) source signatures (analogous)

## Source Signature Construction

### 1. BB (Biomass Burning) δ¹³C

**Method:** C3/C4 vegetation mixing, weighted by CTCH4 pyrogenic emissions per box.

```
BB_δ¹³C(box) = Σ[C4_frac × C4_δ¹³C + (1-C4_frac) × C3_δ¹³C] × emissions(cell) / Σ emissions(box)
```

**Inputs:**
- **C4 vegetation:** Luo et al. (2024) *Nature Communications* 15:1219 — time-varying annual C4 distribution (2001–2019), 0.5° resolution, from `C4_distribution_NUS_v2.2.nc` (249 MB, split into 5×50 MB parts for GitHub)
  - Fallback: Still (2003) static C4 map
  - Luo → BB Trop shifts by −0.9‰ (lower tropical C4 fraction ≈ more C3)
- **Pyrogenic emissions:** CTCH4 2023 (`CTCH4_2023_flux3x2.nc`), expanded from 3°×2° to 1°×1°
- **Isotopic endmembers:** C3 = −26.8 ± 2.9‰, C4 = −12.7 ± 4.6‰ (MC-perturbed)

**Result:** Tropics distinctly less negative (~−25‰) due to C4 savanna burning; NHext/SHext ~−27‰.

### 2. FF (Fossil Fuel) δ¹³C

**Method:** Country-level gas (ONG) vs coal mixing, weighted by EDGAR8 national emissions.

```
FF_δ¹³C(box) = Σ[ONG_frac(country) × ONG_δ¹³C(country) + Coal_frac(country) × Coal_δ¹³C(country)]
             × emissions(country,cell) / Σ emissions(box)
```

**Inputs:**
- **EDGAR8 gridded emissions:** EDGARv8_CH4_1970-2022 at 0.1° (monthly by sector)
  - Sectors 1B1/1B2 (ONG), 1A1/1A2 (coal combustion)
  - Aggregated to countries via country mask
- **Country-level δ¹³C:** `FF_country_d13C.xlsx` (per-country ONG/coal values + uncertainties)
  - Gas-dominated: Russia (−42.3‰), US (−41.4‰), Middle East
  - Coal-dominated: China (−35.7‰), Australia, South Africa
- **Fossil total emissions:** CTCH4 fossil flux for per-cell weighting

**Result:** Strongest 3-box gradient (~6‰); NHext ≈ −43‰ (gas-heavy), SHext ≈ −49‰ (coal-heavy).

### 3. Mic (Microbial) δ¹³C

**Method:** Subcategory mass balance — wetland + ruminant + rice + waste + termite, each with hemisphere/box-specific values.

```
Mic_δ¹³C(box) = Σ frac(subcategory) × δ¹³C(subcategory, box)
```

**Subcategory δ¹³C sources:**
| Subcategory | Fraction | δ¹³C Source | Spatial variation |
|-------------|----------|-------------|-------------------|
| Wetland | 52% | Oh (2022) global ts or isotem spatial map | Global ts → ~0 gap; isotem → ~7‰ gap |
| Ruminant | 25% | Chang (2019) + C3/C4-weighted | Via C4 map → ~2‰ gap |
| Rice | 8% | −59.9 ± 4.5‰ constant | None (Suess only) |
| Waste | 12% | −54.8 ± 4.4‰ constant | None (Suess only) |
| Termite | 3% | −65.2 ± 7.6‰ constant | None (Suess only) |
| Wild animal | 2% | Same as ruminant C3/C4 | Via C4 map |

**Key issue:** With Oh (2022) global wetland δ¹³C, the Mic NH–SH gap is ~0‰ because:
- Wetlands (52%) use a single global value → zero contribution to gradient
- Only ruminants (25%) have C3/C4-weighted spatial variation (~2‰ gap)
- Diluted by 75% spatially-uniform values → net gap ≈ 0.04‰

**Isotem integration** (see separate update): replaces global wetland δ¹³C with spatially-resolved values, introducing a ~4‰ hemispheric gradient in the dominant subcategory.

### MC Uncertainty Propagation

All source signatures use 1000 Monte Carlo iterations:
1. **Isotopic endmembers** perturbed by their uncertainties (normal distribution)
2. **Subcategory fractions** perturbed ±10% relative
3. **Suess effect** perturbed: −0.024 ± 0.005 ‰/yr
4. **Inter-lab scale** uncertainties for δD stations

Output: `{source}_d13C_{box}_MC.csv` — 24 years × 1001 columns (year + 1000 iterations).

## Atmospheric δ¹³C (3-Box)

**Source:** INSTAAR `ch4c13_nh_sh_mean.xlsx` (fortnightly NH/SH means from flask network).

**Approximations:**
- NHext ≈ NH mean (extratropics dominate hemispheric average)
- SHext ≈ SH mean
- Trop ≈ (NH + SH) / 2 (mixing zone, IH gradient only ~0.2‰)

**MC iterations:** Bootstrap resampling of fortnightly observations within each year (1000 iterations).

**Output:**
- `ThreeBox_atm_d13C_annual.csv` (24 years × 3 boxes)
- `ThreeBox_atm_d13C_MC.npz` (24 years × 1000 iterations × 3 boxes)

## Data Flow

```
CTCH4_2023_flux3x2.nc  ─┐
                         ├─→ build_3box_d13C_sources.py ─→ 9 source MC CSVs
Luo C4 map (.nc parts)  ─┤                               + 3 atm d13C MC files
EDGARv8 emissions        ─┤                               + summary CSV
FF_country_d13C.xlsx     ─┤
Oh_2022_Wetlands.xlsx    ─┤
Chang_2019_ruminants.xlsx─┤
ch4c13_nh_sh_mean.xlsx   ─┘

common.py:load_data(three_box=True)  ←─── reads above outputs
  → LoadedData with per-box d13C fields
  → Used by 2x2_three.py, 3x3_three.py
```

## Emission Distribution by Box

| Source | NHext | Trop | SHext |
|--------|-------|------|-------|
| Fossil | 56.4% | 41.8% | 1.7% |
| Microbial | 28.9% | 69.4% | 1.6% |
| Pyrogenic | 29.2% | 69.8% | 1.0% |

## Summary Values (Luo 2024 C4 map)

| Source | NHext (‰) | Trop (‰) | SHext (‰) |
|--------|----------|---------|----------|
| FF | −43.0 | −45.2 | −48.9 |
| BB | −26.6 | −24.9 | −26.6 |
| Mic | −62.1 | −61.9 | −62.1 |

## Relevant Commits

- `1abf033`: Initial 3-box build script (Still 2003 C4)
- `6944502`: Integrated per-box δ¹³C into common.py + model scripts
- `36b61ac`: Added Luo 2024 C4 map, rebuilt all source signatures
