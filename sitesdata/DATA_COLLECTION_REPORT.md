# Sites Data Collection Report

**Date:** 2026-05-17  
**Task:** Collect site-level data for δD-CH₄, δ¹³C-CH₄, and CH₄ (ppb) from ImportantReferences and their supplement datasets.

---

## Summary

| Data Type | Sites | Raw Obs | Timespan | Source |
|-----------|-------|---------|----------|--------|
| **CH₄ (ppb)** | 235 | ~500k+ | 1983–2025 | Basu2022 ObsPack + NOAA GML 2025 |
| **δ¹³C-CH₄** | 64 | ~50k+ | 1998–2022 | Basu2022 ObsPack + NOAA/INSTAAR 2023 |
| **δD-CH₄** | 35 (27 locations) | 7,691 | 1988–2024 | Riddell-Young 2025 dD_GlobMean package |

---

## Data Successfully Collected

### 1. CH₄ Mixing Ratio (ppb)
- **✅ Basu et al. (2022) ObsPack** — 293 NC files, 235 unique stations, multi-lab, through ~2020
  - Already in repo: `ImportantReferences/Basu2022ACP_DS/nc/`
- **✅ NOAA GML Flask Network (Aug 2025)** — 91 stations, latest release through 2024–2025
  - Downloaded from: https://gml.noaa.gov/aftp/data/trace_gases/ch4/flask/surface/

### 2. δ¹³C-CH₄
- **✅ Basu et al. (2022) ObsPack** — 69 NC files, 63 unique stations, through ~2020
  - Already in repo: `ImportantReferences/Basu2022ACP_DS/nc/` (ch4c13_* files)
- **✅ NOAA/INSTAAR Flask Network (Sep 2023)** — 25 stations, 1998–2022
  - Downloaded from: https://gml.noaa.gov/aftp/data/trace_gases/ch4c13/flask/surface/

### 3. δD-CH₄
- **✅ Riddell-Young et al. (2025) dD_GlobMean package** — **35 station files, 7,691 raw observations**
  - Source: `ImportantReferences/Riddell-Young2025PNAS_DS/Riddell-Young_2025_dD_GlobMean/`
  - Path found via: `Old_files_before_organize/improved_dD_pipeline.py` → `BEN_DIR` config
  - Includes raw observations, smoothed MC curves, and global/hemispheric means
  - **4 laboratories**: INSTAAR (16 stations), MPI-BGC (11), IMAU (6), TU/NIPR (2)
  - **23 MBL sites** used for global mean construction
  - **Timespan**: 1988 (IMAU Neumayer) to 2024.9 (MPI-BGC Cape Verde)
  - All data harmonized to MPI scale via Dasgupta et al. (2025) calibration

### 4. Additional Gridded/Model Data (not station-level)
- **✅ Thanwerdas et al. (2024)** — Gridded δ¹³C and δD source signature maps
- **✅ Schwietzke et al. (2016)** — Global fossil fuel δ¹³C source signatures
- **✅ He et al. (2026, JGR)** — Tagged tracer experiments with δ¹³C (1980–2017)
- **✅ Fujita et al. (2025)** — Simulated global δD-CH₄ time series (1750–2015)

---

## Directory Structure

```
sitesdata/
├── DATA_COLLECTION_REPORT.md          ← This file
├── methane_ppb/
│   ├── README.md                      ← Summary: 235 sites, 1983–2025
│   ├── ch4_sites_inventory.csv        ← Full inventory
│   └── ch4_sites_map.png              ← Map
├── isotope_d13C/
│   ├── README.md                      ← Summary: 64 sites, 1998–2022
│   ├── d13C_sites_inventory.csv       ← Full inventory
│   └── d13C_sites_map.png             ← Map
└── isotope_dD/
    ├── README.md                      ← Summary: 35 stations, 1988–2024
    ├── dD_sites_inventory.csv         ← Full inventory with δD stats
    ├── dD_sites_map.png               ← Map (color-coded by lab)
    ├── siteinfo_all_ch4h2.txt         ← Station metadata
    ├── glob_ann_dD.xlsx               ← Published global annual mean
    ├── GlobMean_dD_*.csv/xlsx         ← Global mean time series + MC
    ├── HemMean_dD_*.csv               ← NH/SH means
    ├── raw_observations/              ← 35 files: {site}_01D0_dat.txt
    │   ├── alt_01D0_dat.txt           ← INSTAAR Alert (2005–2009)
    │   ├── altMPI_01D0_dat.txt        ← MPI-BGC Alert (2011–2024)
    │   ├── gvnIMAU_01D0_dat.txt       ← IMAU Neumayer (1988–2020)
    │   ├── nyaNIPR_01D0_dat.txt       ← TU/NIPR Ny-Ålesund (2005–2023)
    │   └── ... (35 files total)
    └── smoothed_curves/               ← 35 files: best-fit curves
        ├── alt_01D0_curves_rsd.txt    ← Best-fit curve + trend + residual
        └── ... (35 files; full MC in ImportantReferences/Riddell-Young2025PNAS_DS/...)
```

---

## References

1. **Basu, S., et al. (2022).** Estimating emissions of methane consistent with atmospheric measurements of methane and δ¹³C of methane. *Atmos. Chem. Phys.*, 22, 15351–15377.

2. **Lan, X., et al. (2025).** Atmospheric Methane from the NOAA GML CCGG Network, 1983–Present. Version 2025-08-15. https://doi.org/10.15138/VNCZ-M766

3. **Michel, S.E., et al. (2023).** INSTAAR δ¹³C-CH₄, 1998–2022. https://doi.org/10.15138/9p89-1x02

4. **Riddell-Young, B., et al. (2025).** Microbial driver of 2006–2023 CH₄ growth indicated by trends in atmospheric δD-CH₄ and δ¹³C-CH₄. *PNAS*.

5. **Dasgupta, B., et al. (2025).** Harmonisation of methane isotope ratio measurements from different laboratories. *EGUsphere*. https://doi.org/10.5194/egusphere-2025-2439

6. **Umezawa, T., et al. (2018).** Interlaboratory comparison of δ¹³C and δD measurements of atmospheric CH₄. *AMT*, 11, 1207–1231.
