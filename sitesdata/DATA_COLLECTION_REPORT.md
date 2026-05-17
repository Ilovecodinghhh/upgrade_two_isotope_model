# Sites Data Collection Report

**Date:** 2026-05-17  
**Task:** Collect site-level data for δD-CH₄, δ¹³C-CH₄, and CH₄ (ppb) from ImportantReferences and their supplement datasets.

---

## Summary

| Data Type | Sites (copied) | Sites (catalogued) | Obs (copied) | Obs (catalogued) | Timespan |
|-----------|---------------|-------------------|-------------|-----------------|----------|
| **CH₄ (ppb)** | 91 (NOAA GML) | 235 (Basu2022) | 222,779 | 7,442,908 | 1983–2025 |
| **δ¹³C-CH₄** | 25 (NOAA/INSTAAR) | 63 (Basu2022) | 46,251 | 56,376 | 1998–2022 |
| **δD-CH₄** | 35 (Riddell-Young) | — | 7,691 | — | 1988–2024 |

---

## Data Successfully Collected & Copied

### 1. δD-CH₄ (Deuterium)

**Source:** Riddell-Young et al. (2025) dD_GlobMean package  
**Location found via:** `Old_files_before_organize/improved_dD_pipeline.py` → `BEN_DIR`

| What | Where |
|------|-------|
| 35 raw observation files | `sitesdata/isotope_dD/raw_observations/` |
| 35 best-fit curve files | `sitesdata/isotope_dD/smoothed_curves/` |
| Station metadata | `sitesdata/isotope_dD/siteinfo_all_ch4h2.txt` |
| Global/hemispheric means | `sitesdata/isotope_dD/GlobMean_dD_*.csv/xlsx`, `HemMean_dD_*.csv` |
| Site inventory | `sitesdata/isotope_dD/dD_sites_inventory.csv` |
| Site map | `sitesdata/isotope_dD/dD_sites_map.png` |

- **4 laboratories**: INSTAAR (16 stations), MPI-BGC (11), IMAU (6), TU/NIPR (2)
- **23 MBL sites** used for global mean construction
- All data harmonized to MPI scale via Dasgupta et al. (2025) calibration
- Full MC curves (148 MB) not copied — available at `ImportantReferences/Riddell-Young2025PNAS_DS/`

### 2. CH₄ Mixing Ratio (ppb)

**Sources:**

| Source | Status | Details |
|--------|--------|---------|
| NOAA GML Flask 2025 | ✅ **Copied** | 91 event + 96 monthly files |
| Basu2022 ObsPack | 📋 **Catalogued** | 292 NC files, ~500 MB total |

| What | Where |
|------|-------|
| 91 event data files | `sitesdata/methane_ppb/noaa_gml_2025_event/` |
| 96 monthly mean files | `sitesdata/methane_ppb/noaa_gml_2025_monthly/` |
| Basu2022 catalogue | `sitesdata/methane_ppb/basu2022_obspack_ch4_catalogue.csv` |
| Site inventory | `sitesdata/methane_ppb/ch4_sites_inventory.csv` |
| Site map | `sitesdata/methane_ppb/ch4_sites_map.png` |

- NOAA GML 2025: Downloaded 2026-05-17 from https://gml.noaa.gov/aftp/data/trace_gases/ch4/flask/surface/
- Basu2022 NC files available at: `ImportantReferences/Basu2022ACP_DS/nc/` (ch4_*.nc)

### 3. δ¹³C-CH₄

**Sources:**

| Source | Status | Details |
|--------|--------|---------|
| NOAA/INSTAAR 2023 | ✅ **Copied** | 25 event + 25 monthly files |
| Basu2022 ObsPack | 📋 **Catalogued** | 70 NC files |

| What | Where |
|------|-------|
| 25 event data files | `sitesdata/isotope_d13C/noaa_instaar_2023_event/` |
| 25 monthly mean files | `sitesdata/isotope_d13C/noaa_instaar_2023_monthly/` |
| Basu2022 catalogue | `sitesdata/isotope_d13C/basu2022_obspack_d13c_catalogue.csv` |
| Site inventory | `sitesdata/isotope_d13C/d13C_sites_inventory.csv` |
| Site map | `sitesdata/isotope_d13C/d13C_sites_map.png` |

- NOAA/INSTAAR 2023: Downloaded 2026-05-17 from https://gml.noaa.gov/aftp/data/trace_gases/ch4c13/flask/surface/
- Basu2022 NC files available at: `ImportantReferences/Basu2022ACP_DS/nc/` (ch4c13_*.nc)

### 4. Additional Gridded/Model Data (not station-level)
- **✅ Thanwerdas et al. (2024)** — Gridded δ¹³C and δD source signature maps
- **✅ Schwietzke et al. (2016)** — Global fossil fuel δ¹³C source signatures
- **✅ He et al. (2026, JGR)** — Tagged tracer experiments with δ¹³C (1980–2017)
- **✅ Fujita et al. (2025)** — Simulated global δD-CH₄ time series (1750–2015)

---

## Directory Structure

```
sitesdata/
├── DATA_COLLECTION_REPORT.md             ← This file
├── methane_ppb/
│   ├── README.md                         ← Documentation & citations
│   ├── ch4_sites_inventory.csv           ← 91 stations inventory
│   ├── ch4_sites_map.png                 ← Map of NOAA GML sites
│   ├── basu2022_obspack_ch4_catalogue.csv ← 292-file ObsPack catalogue
│   ├── noaa_gml_2025_event/              ← 91 discrete flask files
│   │   ├── ch4_alt_surface-flask_1_ccgg_event.txt
│   │   └── ... (91 files)
│   └── noaa_gml_2025_monthly/            ← 96 monthly mean files
│       ├── ch4_alt_surface-flask_1_ccgg_month.txt
│       └── ... (96 files)
├── isotope_d13C/
│   ├── README.md                         ← Documentation & citations
│   ├── d13C_sites_inventory.csv          ← 25 stations inventory
│   ├── d13C_sites_map.png                ← Map of NOAA/INSTAAR sites
│   ├── basu2022_obspack_d13c_catalogue.csv ← 70-file ObsPack catalogue
│   ├── noaa_instaar_2023_event/          ← 25 discrete flask files
│   │   ├── ch4c13_alt_surface-flask_7_sil_event.txt
│   │   └── ... (25 files)
│   └── noaa_instaar_2023_monthly/        ← 25 monthly mean files
│       ├── ch4c13_alt_surface-flask_7_sil_month.txt
│       └── ... (25 files)
└── isotope_dD/
    ├── README.md                         ← Documentation & citations
    ├── dD_sites_inventory.csv            ← 35 stations inventory with δD stats
    ├── dD_sites_map.png                  ← Map (color-coded by lab)
    ├── siteinfo_all_ch4h2.txt            ← Station metadata
    ├── glob_ann_dD.xlsx                  ← Published global annual mean
    ├── GlobMean_dD_*.csv/xlsx            ← Global mean time series + MC
    ├── HemMean_dD_*.csv                  ← NH/SH means
    ├── raw_observations/                 ← 35 files: {site}_01D0_dat.txt
    │   ├── alt_01D0_dat.txt
    │   └── ... (35 files total)
    └── smoothed_curves/                  ← 35 files: best-fit curves
        ├── alt_01D0_curves_rsd.txt
        └── ... (35 files; full MC in ImportantReferences/)
```

---

## References

1. **Basu, S., et al. (2022).** Estimating emissions of methane consistent with atmospheric measurements of methane and δ¹³C of methane. *Atmos. Chem. Phys.*, 22, 15351–15377.

2. **Lan, X., et al. (2025).** Atmospheric Methane from the NOAA GML CCGG Network, 1983–Present. Version 2025-08-15. https://doi.org/10.15138/VNCZ-M766

3. **Michel, S.E., et al. (2023).** INSTAAR δ¹³C-CH₄, 1998–2022. https://doi.org/10.15138/9p89-1x02

4. **Riddell-Young, B., et al. (2025).** Microbial driver of 2006–2023 CH₄ growth indicated by trends in atmospheric δD-CH₄ and δ¹³C-CH₄. *PNAS*.

5. **Dasgupta, B., et al. (2025).** Harmonisation of methane isotope ratio measurements from different laboratories. *EGUsphere*. https://doi.org/10.5194/egusphere-2025-2439

6. **Umezawa, T., et al. (2018).** Interlaboratory comparison of δ¹³C and δD measurements of atmospheric CH₄. *AMT*, 11, 1207–1231.
