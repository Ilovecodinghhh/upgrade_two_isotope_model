# Methane Mixing Ratio (CH₄ ppb) — Site Data Inventory

## Overview

| Metric | Value |
|--------|-------|
| **NOAA GML 2025 sites (copied)** | 91 |
| **Basu2022 ObsPack sites (catalogued)** | 235 |
| **Combined unique stations (est.)** | ~250 |
| **Earliest data** | 1983 |
| **Latest data** | 2024–2025 |
| **NOAA GML observations** | 222,779 |
| **Basu2022 ObsPack observations** | 7,442,908 |

## Data Sources & Files in This Directory

### 1. NOAA GML Flask Network (2025) — **Copied**

The latest release (2025-08-15) of discrete surface flask CH₄ from the NOAA Global Greenhouse Gas Reference Network.

| Folder | Contents | Files |
|--------|----------|-------|
| `noaa_gml_2025_event/` | Individual flask pair measurements | 91 files |
| `noaa_gml_2025_monthly/` | Monthly mean values | 96 files |

**File format:** Plain text with ~168 header lines (metadata, columns), then space-delimited data:
```
# site_code year month day hour minute second datetime value ...
ALT 1985 6 18 14 0 0 1985-06-18T14:00:00Z 1721.84 ...
```

- *Citation:* Lan, X., J.W. Mund, A.M. Crotwell, M.J. Crotwell, E. Moglia, M. Madronich, D. Neff, and K.W. Thoning (2025). Atmospheric Methane Dry Air Mole Fractions from the NOAA GML Carbon Cycle Cooperative Global Air Sampling Network, 1983–Present, Version 2025-08-15. https://doi.org/10.15138/VNCZ-M766
- *Downloaded:* 2026-05-17 from https://gml.noaa.gov/aftp/data/trace_gases/ch4/flask/surface/

### 2. Basu et al. (2022) ObsPack — **Catalogued Only (NC files too large)**

Multi-lab compilation harmonizing NOAA/GML, MPI-BGC, NIWA, NIPR, LSCE, EC, KMA, NILU, EMPA, IPEN, RSE, TU, and others. 293 CH₄ NetCDF files covering 235 unique stations with data through ~2020.

| File | Description |
|------|-------------|
| `basu2022_obspack_ch4_catalogue.csv` | Full catalogue: site, lab, lat/lon, n_obs, timespan, filename |

The original NC files (~500 MB) are in the repo at: `ImportantReferences/Basu2022ACP_DS/nc/` (files matching `ch4_*.nc`)

- *Citation:* Basu, S., et al. (2022). Estimating emissions of methane consistent with atmospheric measurements of methane and δ¹³C of methane. *Atmos. Chem. Phys.*, 22, 15351–15377.

## Station Inventory

The full station inventory for the copied NOAA GML 2025 data is in **`ch4_sites_inventory.csv`** with columns:
`site_code, site_name, country, latitude, longitude, elevation_m, n_observations, year_start, year_end, data_source, file`

## Spatial Coverage

The NOAA GML network (91 sites) provides excellent global coverage:
- **Arctic/boreal**: ALT (82°N), BRW (71°N), ZEP (79°N), SUM (73°N), NYA (79°N)
- **Northern midlatitudes**: Dense coverage across North America, Europe, and East Asia
- **Tropics**: ASC, BKT, CHR, GMI, KUM, MLO, SMO, RPB
- **Southern midlatitudes**: CGO, BHD, AMS, CRZ, MQA
- **Antarctic**: SPO (90°S), PSA (65°S), SYO (69°S), HBA (76°S)

The Basu2022 ObsPack adds ~150 more stations from non-NOAA labs worldwide.

## Timespan Distribution (NOAA GML 2025)

| Period | # of sites with data |
|--------|---------------------|
| Pre-1990 | ~20 |
| 1990–2000 | ~45 |
| 2000–2010 | ~70 |
| 2010–2020 | ~80 |
| 2020–2025 | ~60 |

## Labs Contributing Data

| Lab | Source | Sites |
|-----|--------|-------|
| NOAA/GML CCGG | NOAA GML 2025 (copied) | 91 surface sites |
| NOAA/GML + others | Basu2022 ObsPack (catalogued) | 235 unique stations |
| MPI-BGC | Basu2022 | European and Southern Hemisphere |
| NIWA | Basu2022 | New Zealand, Pacific |
| NIPR/TU | Basu2022 | Japanese Arctic/Antarctic |
| LSCE | Basu2022 | French sites and islands |
| EC (Canada) | Basu2022 | Canadian sites |
| Others | Basu2022 | KMA, NILU, EMPA, IPEN, RSE |

## Notes

- CH₄ values are in **nmol/mol (ppb)** on the **WMO X2004A** scale
- NOAA GML data: discrete flask pair measurements (event files) and monthly means
- Basu2022 ObsPack NC files include both flask and continuous in-situ measurements
- For the most up-to-date data, access NOAA GML: https://gml.noaa.gov/dv/data/
- The Basu2022 ObsPack catalogue CSV references the NC files in `ImportantReferences/Basu2022ACP_DS/nc/`

## Map

See `ch4_sites_map.png` for the spatial distribution of NOAA GML 2025 sites.
