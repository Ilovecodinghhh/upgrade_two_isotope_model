# Sites Data Collection Report

**Date:** 2026-05-17  
**Task:** Collect site-level data for δD-CH₄, δ¹³C-CH₄, and CH₄ (ppb) from ImportantReferences and their supplement datasets.

---

## Summary

| Data Type | Sites Found | Source |
|-----------|-------------|--------|
| **CH₄ (ppb)** | 235 sites | Basu2022 ObsPack + NOAA GML 2025 |
| **δ¹³C-CH₄** | 64 sites | Basu2022 ObsPack + NOAA/INSTAAR 2023 |
| **δD-CH₄** | 18 sites | Literature (Riddell-Young2025, Dasgupta2025) |

---

## Data Successfully Collected & Downloaded

### 1. CH₄ Mixing Ratio (ppb)
- **✅ Basu et al. (2022) ObsPack** — 293 NC files, 235 unique stations, data through ~2020
  - Multi-lab compilation: NOAA/GML, MPI-BGC, NIWA, LSCE, EC, NILU, EMPA, KMA, IPEN, TU, RSE, etc.
  - Already in repo: `ImportantReferences/Basu2022ACP_DS/nc/`
  
- **✅ NOAA GML Flask Network (Aug 2025)** — 91 stations, latest release through 2024-2025
  - Downloaded from: https://gml.noaa.gov/aftp/data/trace_gases/ch4/flask/surface/
  - Most recent version available (2025-08-15)
  - Saved to: `/tmp/noaa_data/ch4_flask/`

### 2. δ¹³C-CH₄
- **✅ Basu et al. (2022) ObsPack** — 69 NC files, 63 unique stations, data through ~2020
  - Labs: NOAA/INSTAAR, MPI-BGC, NIWA, NIPR, TU
  - Already in repo: `ImportantReferences/Basu2022ACP_DS/nc/` (ch4c13_* files)

- **✅ NOAA/INSTAAR Flask Network (Sep 2023)** — 25 stations, data 1998–2022
  - Downloaded from: https://gml.noaa.gov/aftp/data/trace_gases/ch4c13/flask/surface/
  - Latest available version (2023-09-21)
  - Saved to: `/tmp/noaa_data/ch4c13_surface-flask_sil_text/`

- **✅ Chandra et al. (2024) Zenodo** — Downloaded but contains figure-level data only, not raw station data
  - Source: https://doi.org/10.5281/zenodo.10531749
  - Contains aggregated regional CH₄ and δ¹³C time series (NHL, Tropics, SHL)

### 3. δD-CH₄
- **✅ Site inventory compiled** from Riddell-Young et al. (2025, PNAS) SI Table S1 and Dasgupta et al. (2025, EGUsphere)
  - 18 unique station locations across 4 labs

- **✅ Fujita et al. (2025) Dataset S1** — Contains simulated global δD-CH₄ time series (1750–2015)
  - Already in repo: `ImportantReferences/Fujita2025JGR_DS/`
  - Note: These are model output, not station observations

### 4. Source Signature / Gridded Data (not station-level, but useful)
- **✅ Thanwerdas et al. (2024)** — Gridded δ¹³C and δD source signature maps (LMDz9696)
  - Already in repo: `ImportantReferences/Thanwerdas2024ACP_DS/`
  
- **✅ Schwietzke et al. (2016)** — Global fossil fuel δ¹³C source signatures
  - Already in repo: `ImportantReferences/Schwietzke2016Nature_DS/`

- **✅ He et al. (2026, JGR)** — Tagged tracer experiments with δ¹³C (1980–2017)
  - Already in repo: `ImportantReferences/He2026JGR_DS/`

---

## ⚠️ Data Referenced But NOT Publicly Downloadable

### δD-CH₄ Raw Station Data
**This is the most significant gap.** The raw station-level δD-CH₄ time series are NOT available through any standard public data portal. They were compiled directly from individual laboratory archives:

| Lab | Data Status | Contact |
|-----|-------------|---------|
| **INSTAAR** | Not publicly posted; compiled by Riddell-Young et al. for their 2025 PNAS paper | Sylvia Englund Michel (sylvia.michel@colorado.edu) |
| **MPI-BGC** | Available by request from MPI-BGC | Heiko Moossen (hmoossen@bgc-jena.mpg.de) |
| **IMAU (Utrecht)** | Available by request from Thomas Röckmann group | Thomas Röckmann (t.roeckmann@uu.nl) |
| **TU/NIPR** | Available by request; some data on WDCGG (https://gaw.kishou.go.jp/) | Shinji Morimoto; Ryo Fujita (ryo.fujita@mri-jma.go.jp) |

**Recommendation:** Contact the data providers directly to obtain raw station-level δD time series.

### Additional Sources Not Downloaded (Not Site-Level Data)
| Reference | Data Type | Reason Not Downloaded |
|-----------|-----------|----------------------|
| **He et al. (2026, Science)** | TROPOMI satellite CH₄ inversions | Satellite/gridded data, not station-level; partial download in `He2026Science_DS/` |
| **Maasakkers et al. (2019)** | GOSAT satellite inversions | Satellite/gridded data, not station-level |
| **Zhang et al. (2021)** | GOSAT inversion results | Available at https://doi.org/10.5281/zenodo.4052518 — gridded not site-level |
| **Worden et al. (2017)** | Satellite-based methane isotope estimates | Model output, not station-level |
| **Naus et al. (2019)** | OH anomalies from MCF | Not site-level CH₄/isotope data |
| **Rice et al. (2016)** | Isotopic fractionation study | Lab/model study, not site-level atmospheric data |
| **Nguyen et al. (2020)** | OH variability from inverse modeling | Not site-level data |
| **Yu et al. (2026)** | High-latitude isotope trends | Uses same NOAA/INSTAAR + MPI + IMAU + NIPR data already catalogued |
| **Zhao et al. (2023)** | 3D OH field analysis | Not site-level CH₄/isotope data; OH fields from GEOSCCM models |

---

## Data Update Paths

For the **most up-to-date and comprehensive** station data:

### CH₄ (ppb)
1. **NOAA ObsPack** (multi-lab, most comprehensive): https://gml.noaa.gov/ccgg/obspack/
   - Latest release: `obspack_ch4_1_GLOBALVIEWplus_v8.0_2024-10-31`
   - Includes NOAA + partner lab data in unified format

2. **NOAA GML Flask**: https://gml.noaa.gov/aftp/data/trace_gases/ch4/flask/surface/
   - Updated Aug 2025 ✅ (already downloaded)

3. **WDCGG (WMO)**: https://gaw.kishou.go.jp/
   - Additional non-NOAA stations (JMA, TU, NIPR)

### δ¹³C-CH₄
1. **NOAA/INSTAAR**: https://gml.noaa.gov/aftp/data/trace_gases/ch4c13/flask/surface/
   - Latest: Sep 2023 ✅ (already downloaded)
   - DOI: https://doi.org/10.15138/9p89-1x02

2. **ObsPack** (includes MPI, NIWA, etc.): Already captured in Basu2022

### δD-CH₄
1. **No public download portal exists** — contact individual labs
2. The harmonized global mean is available in `rel/data/GlobMean_dD_iterations_DasguptaCal_noBUDS.xlsx`
3. Dasgupta et al. (2025) preprint describes the harmonization: https://doi.org/10.5194/egusphere-2025-2439

---

## Directory Structure

```
sitesdata/
├── DATA_COLLECTION_REPORT.md     ← This file
├── methane_ppb/
│   ├── README.md                 ← Summary of CH₄ sites, abundance, timespan
│   ├── ch4_sites_inventory.csv   ← Full inventory (235 sites)
│   └── ch4_sites_map.png         ← Map of all CH₄ sites
├── isotope_d13C/
│   ├── README.md                 ← Summary of δ¹³C sites, abundance, timespan
│   ├── d13C_sites_inventory.csv  ← Full inventory (64 sites)
│   └── d13C_sites_map.png        ← Map of all δ¹³C sites
└── isotope_dD/
    ├── README.md                 ← Summary of δD sites, abundance, timespan
    ├── dD_sites_inventory.csv    ← Full inventory (18 sites)
    └── dD_sites_map.png          ← Map of all δD sites
```
