# Primary Annual Database

**Description of data sources, processing procedures, and citations for the methane dual-isotope box model.**

Model time span: **1999–2021** (24 years for concentrations, 23 for source partitioning).
Grid resolution: **1° × 1°** (180 × 360), row 0 = 90°N.

---

## 1. Atmospheric CH₄ Concentration

| Item | Value |
|------|-------|
| **Variable** | Annual-mean CH₄ mixing ratio (ppb) |
| **File** | `rel/data/GML_CH4_AnnualMean.xlsx` |
| **Source** | NOAA Global Monitoring Laboratory (GML) |
| **Scope** | Global (1984–2022); model uses rows 15–38 → 1999–2022 |
| **Spatial resolution** | Global mean; hemispheric NH/SH derived via inter-hemispheric gradient model |
| **3-box extension** | CH₄_NHext = global + 30 ppb, CH₄_Trop = global + 10 ppb, CH₄_SHext = global − 25 ppb (from NOAA latitudinal gradient) |

**Processing:**
- Direct read of annual means; no interpolation needed.
- NH/SH split computed via `compute_IH_gradient()` in `common.py` (latitude-dependent mixing).

**Citation:**
> Lan, X., K.W. Thoning, and E.J. Dlugokencky (2024), Trends in globally-averaged CH₄, N₂O, and SF₆ determined from NOAA Global Monitoring Laboratory measurements, Version 2024-09. https://doi.org/10.15138/P8XG-AA10

---

## 2. Atmospheric δ¹³C-CH₄

### 2.1 Global & Hemispheric (2-box)

| Item | Value |
|------|-------|
| **Variable** | Annual-mean δ¹³C-CH₄ (‰ VPDB) |
| **Files** | `rel/data/ch4c13_nh_sh_mean.xlsx` (means), `rel/data/d13C_dei_compiled.txt` (MC) |
| **Source** | NOAA/INSTAAR Cooperative Global Air Sampling Network |
| **Scope** | Global, NH, SH; fortnightly flask observations (1998–2022) |
| **MC iterations** | 1000 bootstrap resamples of fortnightly observations per year |

**Processing:**
- Sub-annual observations averaged to annual means (≥6 data points required per year).
- MC uncertainty: for each of 1000 iterations, resample fortnightly observations with replacement within each year, compute annual mean.
- Columns in `d13C_dei_compiled.txt`: col 0 = year, cols 1–1000 = MC iterations.

### 2.2 Three-Box (NHext / Trop / SHext)

| Item | Value |
|------|-------|
| **Variable** | Annual-mean δ¹³C-CH₄ per box (‰ VPDB) |
| **Files** | `rel/data/ThreeBox_atm_d13C_annual.csv`, `ThreeBox_atm_d13C_{NHext,Trop,SHext}_MC.csv` |
| **Source** | Derived from INSTAAR station data (>80 stations), latitude-binned |
| **Build script** | `rel/build_3box_d13C_sources.py` → Section 6 |

**Processing:**
- Stations assigned to NHext (>30°N), Trop (30°S–30°N), SHext (<30°S) by latitude.
- Annual means from fortnightly flask data; Trop = mean(NH, SH) as proxy where direct binning unavailable.
- MC: bootstrap of per-station fortnightly observations.

**Citation:**
> White, J.W.C., B.H. Vaughn, and S.E. Michel (2023), University of Colorado, Institute of Arctic and Alpine Research (INSTAAR), Stable Isotopic Composition of Atmospheric Methane (¹³C) from the NOAA GML Carbon Cycle Cooperative Global Air Sampling Network, 1998–2022. https://doi.org/10.15138/G3PM-4F05

---

## 3. Atmospheric δD-CH₄

### 3.1 Global

| Item | Value |
|------|-------|
| **Variable** | Annual-mean δD-CH₄ (‰ VSMOW) |
| **File** | `rel/data/GlobMean_dD_iterations_DasguptaCal_noBUDS.xlsx` |
| **Source** | Multi-network compilation: NOAA/INSTAAR, MPI, IMAU, NIPR |
| **Scope** | Global annual mean, 2005–2023 (19 years); pre-2005 front-padded with first value |
| **MC iterations** | ~998 iterations (network + atmospheric + measurement bias uncertainties) |

### 3.2 Hemispheric (2-box)

| Item | Value |
|------|-------|
| **Variable** | Annual-mean δD-CH₄ NH/SH (‰ VSMOW) |
| **Files** | `rel/data/NHMean_dD_iterations_DasguptaCal_noBUDS.xlsx`, `SHMean_dD_iterations_DasguptaCal_noBUDS.xlsx` |
| **Build script** | `dD_hemispheric_MC.py` (in Riddell-Young 2025 dD_GlobMean package) |

### 3.3 Three-Box (NHext / Trop / SHext)

| Item | Value |
|------|-------|
| **Variable** | Annual-mean δD-CH₄ per box (‰ VSMOW) |
| **File** | `rel/data/ThreeBox_atm_dD_annual.csv` (2005–2024, 20 years) |
| **Source** | Derived from semi-hemispheric (PN/TN/TS/PS) station-level MC data |
| **Build script** | `rel/build_3box_dD_sources.py` |

**Processing (all δD):**
1. Station data compiled from NOAA/INSTAAR, MPI, IMAU, NIPR networks; calibrated to MPI scale via Umezawa et al. (2018) inter-lab comparisons.
2. Dasgupta calibration applied; BUDS stations excluded (`noBUDS`).
3. Smoothed Marine Boundary Layer (MBL) fit per station; annual means computed.
4. MC uncertainties propagated: (a) network dropout (2 random stations removed per iteration), (b) atmospheric noise added to monthly means, (c) measurement bias from inter-lab comparison uncertainty.
5. Stations assigned to hemispheres/boxes by latitude; weighted mean of per-station annual values.
6. **Gap fix (commit `cd64c7a`):** Unified date grid spanning all stations (1988.6–2024.9) instead of single reference station. Prevents clipping of post-2020 data from newer MPI/IMAU/NIPR stations.

**Citations:**
> Riddell-Young, E., et al. (2025), Observational constraints on source contributions using a dual methane isotope mass balance, *PNAS*.
>
> Umezawa, T., et al. (2018), Inter-laboratory compatibility of δD measurements of atmospheric CH₄, *Atmospheric Measurement Techniques*, 11(2), 1059–1078. https://doi.org/10.5194/amt-11-1059-2018
>
> Dasgupta, P.K., et al. (2024), [Calibration methodology — see Riddell-Young 2025 SI for details].

---

## 4. CH₄ Emissions (Total & Sectoral)

### 4.1 CarbonTracker-CH₄ (CTCH4)

| Item | Value |
|------|-------|
| **Variable** | Posterior CH₄ fluxes: microbial, fossil, pyrogenic (kg/s, monthly, 3°×2°) |
| **File** | `ImportantReferences/.../data/CTCH4_2023_flux3x2.nc` |
| **Source** | CarbonTracker-CH₄ 2023 release (NOAA/GML) |
| **Scope** | 1997–2021 (288 months), 90×180 (3°lat × 2°lon) |
| **Used for** | Spatial emission weighting of source signatures; total BB emissions |

**Processing:**
- Expanded from 3°×2° to 1°×1° by repeating pixels (`np.repeat(flux, 2, axis=lat)`, `np.repeat(flux, 3, axis=lon) / 6`).
- Flipped so row 0 = 90°N.
- Monthly → annual sum.
- Sectoral annual totals: Microbial ~376 Tg/yr, Fossil ~89 Tg/yr, Pyrogenic ~29 Tg/yr (means).

### 4.2 BB Annual Emissions

| Item | Value |
|------|-------|
| **Variable** | Annual biomass burning (pyrogenic) CH₄ emissions (Tg/yr) |
| **File** | `rel/data/CarbonTracker_CH4.xlsx`, column 9 (= "Pyro Prior") |
| **Scope** | 1997–2021; mean ≈ 29.0 Tg/yr |
| **Used for** | BB-fixed models: BB is prescribed, FF + Mic solved from mass balance |

**Hemispheric/box split:**
- 2-box: NH 55%, SH 45% (GFED4 latitude fractions)
- 3-box: NHext 30%, Trop 55%, SHext 15% (GFED4 latitude fractions)

**Citation:**
> Bruhwiler, L., et al. (2023), CarbonTracker-CH₄, NOAA Global Monitoring Laboratory. https://gml.noaa.gov/ccgg/carbontracker-ch4/

---

## 5. Source Signature End-Members — δ¹³C

### 5.1 Fossil Fuel (FF) δ¹³C

| Item | Value |
|------|-------|
| **Variable** | Emission-weighted annual FF δ¹³C-CH₄ (‰ VPDB) |
| **Files** | Global: `rel/output/FF_d13C_GlobUnc.csv`, `FF_d13C_GlobMC_EDGAR.csv`; Hemispheric/3-box: `rel/data/FF_d13C_{NH,SH,NHext,Trop,SHext}_MC.csv` |
| **Build scripts** | `rel/build_hemispheric_d13C_sources.py`, `rel/build_3box_d13C_sources.py` |

**Method:**
- Country-level ONG (oil & natural gas) and coal δ¹³C end-members from published compilations.
- Emission-weighted by EDGAR 8.0 country-level CH₄ emissions.
- Countries assigned to NH/SH or 3-box by centroid latitude.
- MC: sampling country-level δ¹³C uncertainties.

**Typical values (1998–2021 means):**
- Global: −44.1‰; NH: −43.4‰; SH: −48.0‰
- NHext: −43.0‰; Trop: −45.2‰; SHext: −48.9‰

**Citations:**
> Sherwood, O.A., et al. (2017), Global inventory of gas geochemistry data from fossil fuel, microbial and burning sources, version 2017, *Earth System Science Data*, 9, 639–656. https://doi.org/10.5194/essd-9-639-2017
>
> Crippa, M., et al. (2024), EDGAR v8.0 Global Greenhouse Gas Emissions. https://edgar.jrc.ec.europa.eu/

### 5.2 Biomass Burning (BB) δ¹³C

| Item | Value |
|------|-------|
| **Variable** | Emission-weighted annual BB δ¹³C-CH₄ (‰ VPDB) |
| **Files** | Global: `rel/output/BB_d13C_annual.csv`; Hemispheric/3-box: `rel/data/BB_d13C_{NH,SH,NHext,Trop,SHext}_MC.csv` |
| **Build scripts** | `rel/build_hemispheric_d13C_sources.py`, `rel/build_3box_d13C_sources.py` |

**Method:**
1. Load C3/C4 vegetation distribution maps.
2. Weight C3 (δ¹³C = −26.8 ± 2.9‰) and C4 (δ¹³C = −12.7 ± 4.6‰) end-members by vegetation fraction per grid cell.
3. Emission-weight by CTCH4 pyrogenic flux per hemisphere/box.
4. MC: sample C3/C4 δ¹³C uncertainties + vegetation fraction.

**C3/C4 vegetation map:**
- **Primary:** Luo et al. (2024) time-varying C4 distribution, 0.5° resolution, 2001–2019.
  - File: `C4_distribution_NUS_v2.2.nc` (249 MB, split into 5×50 MB parts for GitHub).
  - Regridded 0.5° → 1° by averaging 2×2 blocks; years outside 2001–2019 padded with nearest year.
- **Fallback:** Still & Berry (2003) static C4 map.

**Typical values (1998–2021 means):**
- NH: −26.0‰; SH: −24.2‰ (SH has more tropical C4 savanna fires)
- NHext: −26.6‰; Trop: −24.9‰; SHext: −26.6‰

**Citations:**
> Cerling, T.E., et al. (1997), Global vegetation change through the Miocene/Pliocene boundary, *Nature*, 389, 153–158.
>
> Luo, X., et al. (2024), A global gridded C3/C4 vegetation distribution dataset at 0.5-degree resolution, *Nature Communications*, 15, 1219. https://doi.org/10.1038/s41467-024-45264-9 (Data: Zenodo, https://doi.org/10.5281/zenodo.10516423)
>
> Still, C.J., J.A. Berry, G.J. Collatz, and R.S. DeFries (2003), Global distribution of C3 and C4 vegetation: Carbon cycle implications, *Global Biogeochemical Cycles*, 17(1), 6-1–6-14.

### 5.3 Microbial (Mic) δ¹³C

| Item | Value |
|------|-------|
| **Variable** | Emission-weighted annual microbial δ¹³C-CH₄ (‰ VPDB) |
| **Files** | Global: `rel/output/Mic_d13C_annual.csv`, `Mic_d13C_MC.csv`; Hemispheric/3-box: `rel/data/Mic_d13C_{NH,SH,NHext,Trop,SHext}_MC.csv` |
| **Build scripts** | `rel/build_hemispheric_d13C_sources.py`, `rel/build_3box_d13C_sources.py` |

**Method — Subcategory Mass Balance:**

Microbial CH₄ = Wetlands + Ruminants + Rice + Termites + Waste + Wild Animals

For each hemisphere/box, the composite δ¹³C is:

δ¹³C_mic = Σᵢ (fᵢ × δ¹³Cᵢ)

where fᵢ are emission-weighted subcategory fractions.

| Subcategory | δ¹³C (‰) | Uncertainty (‰) | Source |
|-------------|-----------|-----------------|--------|
| Wetlands | Spatially varying (isotem) | Per-cell monthly σ | Parker et al. (2022) |
| Ruminants | C3: −66.8 ± 2.8; C4: −51.4 ± 3.4 | Per hemisphere via C4 map | Chang et al. (2019) |
| Rice | −63.0 | ±5.0 | Riddell-Young (2025) |
| Termites | −57.0 | ±10.0 | Riddell-Young (2025) |
| Waste/Landfill | −55.0 | ±5.0 | Riddell-Young (2025) |
| Wild Animals | −66.0 | ±5.0 | Riddell-Young (2025) |

**Wetland δ¹³C spatial source (isotem):**
- **Primary:** Per-year isotem wetland δ¹³C-CH₄ maps, 1984–2016.
  - Files: `rel/data/isotem_wetland_d13C-CH4/isotem_wetland_d13C-CH4_{YEAR}.nc4` (33 files)
  - Variable: `wetland_d13C-CH4`, dims (12, 720, 360) = (month, longitude, latitude), 0.5° resolution.
  - ~75.9% NaN (ocean/desert); emission-weighting uses only cells with valid data.
  - Regridded 0.5° → 1° via nanmean of 2×2 blocks; years 2017–2021 padded with 2016.
- **Fallback:** Oh et al. (2022) global time series (~−61‰ ± 0.7‰).

**Emission-weighted wetland δ¹³C (isotem):**
- Raw grid: NH = −60.4‰, SH = −55.7‰ (gap 4.7‰)
- Emission-weighted: NH = −57.4‰, SH = −57.2‰ (gap collapses — tropical wetlands dominate both hemispheres)
- 3-box: NHext = −59.3‰, Trop = −56.2‰, SHext = −55.4‰

**Suess effect correction:** −0.024 ± 0.005 ‰/yr applied to all subcategories relative to reference year 2010.

**Typical Mic δ¹³C values (1998–2021 means):**
- NH: −59.9‰; SH: −59.7‰
- NHext: −61.1‰; Trop: −59.1‰; SHext: −59.0‰

**Citations:**
> Parker, R.J., et al. (2022), Isotopically-resolved methane emissions from global wetland and non-wetland sources (isotem), Zenodo. (Files: `isotem_wetland_d13C-CH4_{YEAR}.nc4`)
>
> Oh, Y., et al. (2022), Improved global wetland carbon isotopic signatures of methane emissions, *Global Biogeochemical Cycles*, 36, e2021GB007049.
>
> Chang, J., et al. (2019), Revisiting enteric methane emissions from domestic ruminants and their δ¹³C-CH₄ source signature, *Nature Communications*, 10, 3420.
>
> Riddell-Young, E., et al. (2025), Observational constraints on source contributions using a dual methane isotope mass balance, *PNAS*.

---

## 6. Source Signature End-Members — δD

### 6.1 Fossil Fuel (FF) δD

| Item | Value |
|------|-------|
| **Variable** | Emission-weighted annual FF δD-CH₄ (‰ VSMOW) |
| **Files** | Global: `rel/output/FF_dD_GlobUnc.csv`, `FF_dD_GlobMC_EDGAR.csv`; Hemispheric/3-box: `rel/data/FF_dD_{NH,SH,NHext,Trop,SHext}_MC.csv` |
| **Build scripts** | `rel/build_hemispheric_dD_sources.py`, `rel/build_3box_dD_sources.py` |

**Method:**
- Country-level ONG and coal δD from published compilations.
- EDGAR 8.0 emission-weighted, same country-to-hemisphere/box assignment as δ¹³C.
- MC: 1000 iterations sampling country-level δD uncertainties.

### 6.2 Biomass Burning (BB) δD

| Item | Value |
|------|-------|
| **Variable** | Emission-weighted annual BB δD-CH₄ (‰ VSMOW) |
| **Files** | Global: `rel/output/BB_dD_annual.csv`; Hemispheric/3-box: `rel/data/BB_dD_{NH,SH,NHext,Trop,SHext}_MC.csv` |
| **Build scripts** | `rel/build_hemispheric_dD_sources.py`, `rel/build_3box_dD_sources.py` |

**Method:**
- Umezawa et al. (2011) regression: δD = 1.16 × MAT − 177 (‰), applied to mean annual temperature grid.
- Temperature from `d2h_MA.tif` (gridded δ²H of precipitation proxy for MAT).
- Emission-weighted by CTCH4 pyrogenic flux per hemisphere/box.

**Citation:**
> Umezawa, T., et al. (2011), Seasonally resolved source contributions to atmospheric methane using δ¹³C and δD isotope ratios, *Journal of Geophysical Research*, 116, D02308.

### 6.3 Microbial (Mic) δD

| Item | Value |
|------|-------|
| **Variable** | Emission-weighted annual microbial δD-CH₄ (‰ VSMOW) |
| **Files** | Global: `rel/output/Mic_dD_AnnGlob.csv`, `Mic_dD_MC.csv`; Hemispheric/3-box: `rel/data/Mic_dD_{NH,SH,NHext,Trop,SHext}_MC.csv` |
| **Build scripts** | `rel/build_hemispheric_dD_sources.py`, `rel/build_3box_dD_sources.py` |

**Method:**
- Douglas et al. (2021) regression: δD = 0.6088 × MAT − 285.7 (‰), applied to `d2h_MA.tif` temperature grid.
- Emission-weighted by CTCH4 microbial flux per hemisphere/box.
- MC: sampling regression uncertainty + temperature grid uncertainty.

**Citation:**
> Douglas, P.M.J., et al. (2021), Global estimates of methane δD, *Global Biogeochemical Cycles*, 35, e2020GB006858.

---

## 7. Auxiliary Datasets

### 7.1 Kinetic Isotope Effects (KIE)

| Sink | KIE ¹³C | KIE D | Distribution |
|------|---------|-------|-------------|
| OH | 1.0039–1.0054 | 1.294–1.327 | Uniform |
| Cl | 1.066 ± 0.002 | 1.52 ± 0.02 | Normal |
| Stratosphere | 1.003 ± 0.001 | 1.179 ± 0.01 | Normal |
| Soil | 1.0201 ± 0.003 | 1.083 ± 0.01 | Normal |

Bulk KIE computed as emission-weighted mean of individual sinks.

**Sink fractions (proportion of total loss):**

| | OH | Cl | Strat | Soil |
|---|---|---|---|---|
| Global | 0.835 | 0.035 | 0.070 | 0.060 |
| NH | 0.825 | 0.040 | 0.070 | 0.065 |
| SH | 0.850 | 0.028 | 0.070 | 0.052 |
| NHext | 0.810 | 0.040 | 0.080 | 0.070 |
| Trop | 0.860 | 0.035 | 0.055 | 0.050 |
| SHext | 0.840 | 0.025 | 0.080 | 0.055 |

### 7.2 CH₄ Lifetime

- **Varying mode** (default): τ(t) = 9.0 − 0.017 × (t − 2010) years, following He et al. (2026).
- **Fixed mode** (optional): τ = 9.0 years.
- 3-box lifetime ratios: NHext × 1.05, Trop × 0.90, SHext × 1.08 (relative to global τ).

### 7.3 Inter-Hemispheric Exchange

- **2-box:** Single exchange time (implicit in mass balance).
- **3-box:**
  - NHext ↔ Trop: τ_ex = 0.8 ± 0.1 yr (faster, no ITCZ barrier).
  - Trop ↔ SHext: τ_ex = 1.2 ± 0.1 yr (ITCZ barrier slows exchange).

---

## 8. Output File Inventory

All MC output files follow the format: col 0 = year (1998–2021), cols 1–1000 = MC iterations.

### δ¹³C Source Signatures
| File | Description |
|------|-------------|
| `rel/data/{Mic,BB,FF}_d13C_{NH,SH}_MC.csv` | 2-box hemispheric (6 files) |
| `rel/data/{Mic,BB,FF}_d13C_{NHext,Trop,SHext}_MC.csv` | 3-box (9 files) |
| `rel/data/Hemispheric_d13C_sources_summary.csv` | Summary means ± σ |
| `rel/data/ThreeBox_d13C_sources_summary.csv` | Summary means ± σ |

### δD Source Signatures
| File | Description |
|------|-------------|
| `rel/data/{Mic,BB,FF}_dD_{NH,SH}_MC.csv` | 2-box hemispheric (6 files) |
| `rel/data/{Mic,BB,FF}_dD_{NHext,Trop,SHext}_MC.csv` | 3-box (9 files) |
| `rel/data/Hemispheric_dD_sources_summary.csv` | Summary means ± σ |
| `rel/data/ThreeBox_dD_sources_summary.csv` | Summary means ± σ |

### Atmospheric δ¹³C
| File | Description |
|------|-------------|
| `rel/data/ch4c13_nh_sh_mean.xlsx` | INSTAAR NH/SH annual means |
| `rel/data/d13C_dei_compiled.txt` | Global MC (24yr × 1000 iter) |
| `rel/data/ThreeBox_atm_d13C_annual.csv` | 3-box annual means |
| `rel/data/ThreeBox_atm_d13C_{NHext,Trop,SHext}_MC.csv` | 3-box MC |

### Atmospheric δD
| File | Description |
|------|-------------|
| `rel/data/GlobMean_dD_iterations_DasguptaCal_noBUDS.xlsx` | Global MC |
| `rel/data/NHMean_dD_iterations_DasguptaCal_noBUDS.xlsx` | NH MC |
| `rel/data/SHMean_dD_iterations_DasguptaCal_noBUDS.xlsx` | SH MC |
| `rel/data/ThreeBox_atm_dD_annual.csv` | 3-box annual means |

---

## 9. Build Pipeline

```
                    ┌─────────────────────┐
                    │  Raw Station Data    │
                    │  (INSTAAR, MPI,      │
                    │   IMAU, NIPR)        │
                    └──────────┬──────────┘
                               │
                    ┌──────────▼──────────┐
                    │  dD_globmean.py      │  → Global + semi-hemispheric δD MC
                    │  dD_hemispheric_MC.py│  → NH/SH δD MC
                    └──────────┬──────────┘
                               │
         ┌─────────────────────┼─────────────────────┐
         │                     │                      │
┌────────▼─────────┐ ┌────────▼─────────┐ ┌─────────▼────────┐
│ build_hemispheric │ │ build_3box       │ │ build_3box       │
│ _d13C_sources.py  │ │ _d13C_sources.py │ │ _dD_sources.py   │
│ build_hemispheric │ │                  │ │                  │
│ _dD_sources.py    │ │                  │ │                  │
└────────┬─────────┘ └────────┬─────────┘ └─────────┬────────┘
         │                     │                      │
         │    Inputs:          │    Inputs:            │    Inputs:
         │    • CTCH4 fluxes   │    • CTCH4 fluxes     │    • CTCH4 fluxes
         │    • Luo C4 map     │    • Luo C4 map       │    • d2h_MA.tif
         │    • isotem wetland │    • isotem wetland    │    • EDGAR 8.0
         │    • EDGAR 8.0      │    • EDGAR 8.0        │    • SemiHem dD data
         │    • Chang ruminant │    • Chang ruminant    │
         │                     │                       │
         ▼                     ▼                       ▼
┌──────────────────────────────────────────────────────────────┐
│                      rel/data/                               │
│  {Mic,BB,FF}_{d13C,dD}_{NH,SH,NHext,Trop,SHext}_MC.csv     │
│  ThreeBox_atm_{d13C,dD}_annual.csv                           │
│  *_sources_summary.csv                                       │
└──────────────────────────────────────────────────────────────┘
         │
         ▼
┌──────────────────────────────────────────────────────────────┐
│  common.py :: load_data(base_dir, two_box, three_box)        │
│  → LoadedData dataclass with all fields populated            │
└──────────────────────────────────────────────────────────────┘
         │
         ▼
┌──────────────────────┐
│  Box Model Scripts   │
│  2x2_one.py          │  ← 1-box, dual isotope
│  2x2_two.py          │  ← 2-box (NH/SH), dual isotope
│  2x2_three.py        │  ← 3-box (NHext/Trop/SHext), dual isotope, separate δ¹³C + δD
│  3x3_one.py          │  ← 1-box, 3×3 matrix solve
│  3x3_two.py          │  ← 2-box, 3×3 matrix
│  3x3_three.py        │  ← 3-box, 3×3 matrix (full)
└──────────────────────┘
```

---

## 10. Full Citation List

1. **Bruhwiler, L., et al. (2023).** CarbonTracker-CH₄, NOAA GML. https://gml.noaa.gov/ccgg/carbontracker-ch4/
2. **Cerling, T.E., et al. (1997).** Global vegetation change through the Miocene/Pliocene boundary. *Nature*, 389, 153–158.
3. **Chang, J., et al. (2019).** Revisiting enteric methane emissions from domestic ruminants and their δ¹³C-CH₄ source signature. *Nature Communications*, 10, 3420.
4. **Crippa, M., et al. (2024).** EDGAR v8.0 Global Greenhouse Gas Emissions. https://edgar.jrc.ec.europa.eu/
5. **Douglas, P.M.J., et al. (2021).** Global estimates of methane δD. *Global Biogeochemical Cycles*, 35, e2020GB006858.
6. **He, J., et al. (2026).** [Methane lifetime trend] — referenced in model config.
7. **Lan, X., K.W. Thoning, and E.J. Dlugokencky (2024).** Trends in globally-averaged CH₄. NOAA GML. https://doi.org/10.15138/P8XG-AA10
8. **Luo, X., et al. (2024).** A global gridded C3/C4 vegetation distribution dataset. *Nature Communications*, 15, 1219.
9. **Oh, Y., et al. (2022).** Improved global wetland carbon isotopic signatures of methane emissions. *Global Biogeochemical Cycles*, 36, e2021GB007049.
10. **Parker, R.J., et al. (2022).** Isotopically-resolved methane emissions from global wetland and non-wetland sources (isotem). Zenodo.
11. **Riddell-Young, E., et al. (2025).** Observational constraints on source contributions using a dual methane isotope mass balance. *PNAS*.
12. **Sherwood, O.A., et al. (2017).** Global inventory of gas geochemistry data from fossil fuel, microbial and burning sources. *ESSD*, 9, 639–656.
13. **Still, C.J., et al. (2003).** Global distribution of C3 and C4 vegetation. *Global Biogeochemical Cycles*, 17(1), 6-1–6-14.
14. **Umezawa, T., et al. (2011).** Seasonally resolved source contributions to atmospheric methane. *JGR*, 116, D02308.
15. **Umezawa, T., et al. (2018).** Inter-laboratory compatibility of δD measurements of atmospheric CH₄. *AMT*, 11(2), 1059–1078.
16. **White, J.W.C., B.H. Vaughn, and S.E. Michel (2023).** INSTAAR δ¹³C-CH₄ from NOAA GML network. https://doi.org/10.15138/G3PM-4F05
