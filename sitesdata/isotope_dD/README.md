# δD-CH₄ (Deuterium) — Site Data Inventory

## Overview

| Metric | Value |
|--------|-------|
| **Total station entries** | 35 |
| **Unique physical locations** | 27 |
| **MBL (Marine Boundary Layer) sites** | 23 |
| **Northern Hemisphere** | 24 |
| **Southern Hemisphere** | 11 |
| **Earliest data** | 1988.6 (IMAU Neumayer, 1988) |
| **Latest data** | 2024.9 |
| **Total raw observations** | 7,691 |
| **Labs contributing** | 4 (INSTAAR, MPI-BGC, IMAU, TU/NIPR) |

## Data Source

All δD-CH₄ station data come from the **Riddell-Young et al. (2025) `dD_GlobMean` package**, which is the supplementary data for:

> **Riddell-Young, B., et al. (2025).** Microbial driver of 2006–2023 CH₄ growth indicated by trends in atmospheric δD-CH₄ and δ¹³C-CH₄. *Proceedings of the National Academy of Sciences (PNAS)*.

This package harmonizes δD-CH₄ measurements from four independent laboratories worldwide:

| Lab | Full Name | # Stations | Period | Reference |
|-----|-----------|------------|--------|-----------|
| **INSTAAR** | Institute of Arctic and Alpine Research, CU-Boulder / NOAA GML | 16 | 2004–2010 | Michel et al. (2023); White et al. (2023) |
| **MPI-BGC** | Max Planck Institute for Biogeochemistry, Jena | 11 | 2010–2024 | Moossen et al. (unpublished) |
| **IMAU** | Utrecht University, Institute for Marine and Atmospheric Research | 6 | 1988–2024 | Röckmann et al.; longest δD record |
| **TU/NIPR** | Tohoku University / National Institute of Polar Research | 2 | 1995–2023 | Morimoto & Fujita |

### Inter-Laboratory Calibration

All data are harmonized to the **MPI scale** using Dasgupta et al. (2025) calibration offsets (applied in the pipeline):

| Lab | Scale Adjustment (‰) | Uncertainty (‰) |
|-----|----------------------|-----------------|
| INSTAAR | +1.8 (Dasgupta) | ±1.6 |
| MPI-BGC | 0.0 (reference) | 0.0 |
| IMAU | +0.5 (Dasgupta) | ±2.2 |
| TU/NIPR | 0.0 | ±1.6 |

> **Dasgupta, B., et al. (2025).** Harmonisation of methane isotope ratio measurements from different laboratories using atmospheric samples. *EGUsphere* [preprint]. https://doi.org/10.5194/egusphere-2025-2439

## Station Details

### MBL (Marine Boundary Layer) Sites — Used for Global Mean Construction

| Site ID | Location | Country | Lat | Lon | Lab | Period | N obs |
|---------|----------|---------|-----|-----|-----|--------|-------|
| alt | Alert, Nunavut | Canada | 82.5 | -62.5 | INSTAAR | 2005–2010 | 237 |
| altMPI | Alert | Canada | 82.5 | -62.5 | MPI-BGC | 2011–2025 | 337 |
| vrsMPI | Villum Research Station, Stati | Greenland | 81.6 | -16.6 | MPI-BGC | 2020–2024 | 142 |
| nyaNIPR | nyvolesund station | Norway | 78.9 | 11.9 | TU/NIPR | 2006–2024 | 800 |
| zep | Ny-Alesund, Svalbard | Norway | 78.9 | 11.9 | INSTAAR | 2008–2010 | 35 |
| zepIMAU | Ny-Alesund, Svalbard | Norway | 78.9 | 11.9 | IMAU | 2013–2024 | 600 |
| brw | Barrow Atmospheric Baseline Ob | United States | 71.3 | -156.6 | INSTAAR | 2005–2010 | 106 |
| brwIMAU | Barrow Atmospheric Baseline Ob | United States | 71.3 | -156.6 | IMAU | 2022–2024 | 127 |
| sisMPI | Shetland Islands | United Kingdom | 59.9 | -1.3 | MPI-BGC | 2011–2025 | 646 |
| cba | Cold Bay, Alaska | United States | 55.2 | -162.7 | INSTAAR | 2005–2010 | 155 |
| mhd | Mace Head, County Galway | Ireland | 53.3 | -9.9 | INSTAAR | 2005–2010 | 83 |
| kum | Cape Kumukahi, Hawaii | United States | 19.6 | -154.9 | INSTAAR | 2005–2010 | 120 |
| cvoMPI | Cape Verde Atmospheric Observa | Cape Verde | 16.9 | -24.9 | MPI-BGC | 2011–2025 | 575 |
| asc | Ascension Island | United Kingdom | -8.0 | -14.4 | INSTAAR | 2005–2010 | 192 |
| smo | Tutuila | American Samoa | -14.2 | -170.6 | INSTAAR | 2005–2010 | 115 |
| smoIMAU | Tutuila | American Samoa | -14.2 | -170.6 | IMAU | 2022–2024 | 104 |
| namMPI | Gobabeb station | Namibia | -23.6 | 15.1 | MPI-BGC | 2013–2022 | 307 |
| cgo | Cape Grim, Tasmania | Australia | -40.7 | 144.7 | INSTAAR | 2005–2009 | 69 |
| cgoIMAU | Cape Grim, Tasmania | Australia | -40.7 | 144.7 | IMAU | 2023–2024 | 59 |
| syoNIPR | Syowa station | Antarctica | -69.0 | 39.6 | TU/NIPR | 1995–2023 | 162 |
| gvnMPI | Neumayer Station | Antarctica | -70.7 | -8.3 | MPI-BGC | 2018–2024 | 245 |
| gvnIMAU | Neumayer Station | Antarctica | -70.7 | -8.3 | IMAU | 1989–2020 | 163 |
| spo | South Pole | Antarctica | -90.0 | -24.8 | INSTAAR | 2005–2010 | 144 |

### Non-MBL Sites (continental/elevated, not used in global mean)

| Site ID | Location | Country | Lat | Lon | Lab | Period | N obs |
|---------|----------|---------|-----|-----|-----|--------|-------|
| kjnMPI | kjolnes | Norway | 70.8 | 29.2 | MPI-BGC | 2015–2021 | 176 |
| zotMPI | Zotino Tall Tower Observatory, | Russia | 60.8 | 89.3 | MPI-BGC | 2010–2022 | 633 |
| bal | Baltic Sea | Poland | 55.4 | 17.2 | INSTAAR | 2005–2010 | 221 |
| bikMPI | Bialystok | Poland | 53.2 | 23.0 | MPI-BGC | 2013–2021 | 151 |
| oxkMPI | Ochsenkopf | Germany | 50.0 | 11.8 | MPI-BGC | 2013–2019 | 243 |
| jfjMPI | Jungfraujoch | Switzerland | 46.5 | 8.0 | MPI-BGC | 2013–2025 | 253 |
| lef | Park Falls, Wisconsin | United States | 46.0 | -90.3 | INSTAAR | 2005–2008 | 75 |
| bsc | Black Sea, Constanta | Romania | 44.2 | 28.7 | INSTAAR | 2005–2008 | 64 |
| azr | Terceira Island, Azores | Portugal | 38.8 | -27.4 | INSTAAR | 2005–2010 | 71 |
| mlo | Mauna Loa, Hawaii | United States | 19.5 | -155.6 | INSTAAR | 2005–2010 | 202 |
| mloIMAU | Mauna Loa, Hawaii | United States | 19.5 | -155.6 | IMAU | 2023–2024 | 58 |
| ato | Amazonas | Brazil | -2.1 | -59.0 | INSTAAR | 2022–2023 | 21 |

## Data Files in This Directory

### Raw Observations (`raw_observations/`)
- **35 files**, one per station, named `{site_id}_01D0_dat.txt`
- Format: two columns — `decimal_year` and `δD_permil` (‰ VSMOW)
- These are the individual flask measurements before any smoothing

### Smoothed Curves (`smoothed_curves/`)
- **35 `_curves_rsd.txt` files** — Best-fit smoothed curve + trend + residual
  - Columns: `decimal_date`, `smoothed_curve`, `trend`, `residual_std_dev`
- **Note**: The full Monte Carlo smoothed curves (1000 iterations per station, `*_smoothedMC.txt`) are in the original package at:
  `ImportantReferences/Riddell-Young2025PNAS_DS/Riddell-Young_2025_dD_GlobMean/Riddell-Young_2025_dD_GlobMean/output/`

### Global/Hemispheric Means
- `glob_ann_dD.xlsx` — Published annual global mean δD-CH₄
- `GlobMean_dD_dei_DasguptaCal_noBUDS.csv` — Global mean time series (Dasgupta calibration)
- `GlobMean_dD_iterations_DasguptaCal_noBUDS.xlsx` — MC iterations for global mean
- `HemMean_dD_annual_DasguptaCal_noBUDS.csv` — NH/SH annual means
- `HemMean_dD_dei_DasguptaCal_noBUDS.csv` — NH/SH weekly means

### Metadata
- `siteinfo_all_ch4h2.txt` — Station metadata (site|name|country|lat|lon|elev|intake_ht|utc2lst|n/a|mbl_flag)
- `dD_sites_inventory.csv` — Full inventory CSV with all fields

## Spatial Coverage

The δD-CH₄ network covers:
- **High Arctic (70–83°N)**: ALT, BRW, NYA/ZEP, KJN, VRS, SIS (6 locations)
- **Northern midlatitudes (30–70°N)**: MHD, BAL, BIK, OXK, JFJ, LEF, CBA, AZR, BSC, ZOT (10 locations)
- **Subtropics/Tropics (30°S–30°N)**: KUM, MLO, CVO, ASC, ATO, SMO (6 locations)
- **Southern midlatitudes (30–70°S)**: CGO, NAM (2 locations)
- **Antarctic (>60°S)**: GVN, SYO, SPO (3 locations)

**⚠️ Key gap**: Limited tropical coverage, especially over continental tropics (only ATO/Amazonas, recently started 2021).

## Full Citations

1. **Riddell-Young, B., Englund Michel, S., Lan, X., Tans, P., Röckmann, T., Dasgupta, B., Oh, Y., Bruhwiler, L., Fujita, R., Umezawa, T., Morimoto, S., and Miller, J. (2025).** Microbial driver of 2006–2023 CH₄ growth indicated by trends in atmospheric δD-CH₄ and δ¹³C-CH₄. *PNAS*.

2. **Dasgupta, B., et al. (2025).** Harmonisation of methane isotope ratio measurements from different laboratories using atmospheric samples. *EGUsphere* [preprint]. https://doi.org/10.5194/egusphere-2025-2439

3. **Umezawa, T., et al. (2018).** Interlaboratory comparison of δ¹³C and δD measurements of atmospheric CH₄ for combined use of data sets from different laboratories. *Atmospheric Measurement Techniques*, 11, 1207–1231.

4. **Michel, S.E., Clark, J.R., Vaughn, B.H., et al. (2023).** INSTAAR Stable Isotopic Composition of Atmospheric Methane. NOAA GML Carbon Cycle Cooperative Global Air Sampling Network. https://doi.org/10.15138/G3PM-4F05

## Map

See `dD_sites_map.png` for the spatial distribution of all sites, color-coded by laboratory.
