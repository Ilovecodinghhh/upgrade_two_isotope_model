# δ¹³C-CH₄ — Site Data Inventory

## Overview

| Metric | Value |
|--------|-------|
| **NOAA/INSTAAR 2023 sites (copied)** | 25 |
| **Basu2022 ObsPack sites (catalogued)** | 63 |
| **Combined unique stations (est.)** | ~65 |
| **Earliest data** | 1998 |
| **Latest data** | 2022 |
| **NOAA/INSTAAR observations** | 46,251 |
| **Basu2022 ObsPack observations** | 56,376 |

## Data Sources & Files in This Directory

### 1. NOAA/INSTAAR Flask Network (2023) — **Copied**

The latest release (2023-09-21) of δ¹³C-CH₄ from the NOAA GML / INSTAAR Carbon Cycle Cooperative Global Air Sampling Network.

| Folder | Contents | Files |
|--------|----------|-------|
| `noaa_instaar_2023_event/` | Individual flask measurements | 25 files |
| `noaa_instaar_2023_monthly/` | Monthly mean values | 25 files |

**File format:** Plain text with ~170 header lines (metadata, columns), then space-delimited data:
```
# site_code year month day hour minute second datetime value ...
ALT 2000 1 12 12 35 0 2000-01-12T12:35:00Z -47.42 ...
```

- *Citation:* Michel, S.E., D. Butowicz, and B.H. Vaughn (2023). INSTAAR Stable Isotopic Composition of Atmospheric Methane (¹³C) from the NOAA GML Carbon Cycle Cooperative Global Air Sampling Network, 1998–2022, Version 2023-09-21. https://doi.org/10.15138/9p89-1x02
- *Downloaded:* 2026-05-17 from https://gml.noaa.gov/aftp/data/trace_gases/ch4c13/flask/surface/

### 2. Basu et al. (2022) ObsPack — **Catalogued Only (NC files too large)**

Multi-lab compilation (NOAA/INSTAAR, MPI-BGC, NIWA, NIPR, TU) with 69 δ¹³C-CH₄ NetCDF files covering 63 unique stations, data through ~2020.

| File | Description |
|------|-------------|
| `basu2022_obspack_d13c_catalogue.csv` | Full catalogue: site, lab, lat/lon, n_obs, timespan, filename |

The original NC files are in the repo at: `ImportantReferences/Basu2022ACP_DS/nc/` (files matching `ch4c13_*.nc`)

- *Citation:* Basu, S., et al. (2022). Estimating emissions of methane consistent with atmospheric measurements of methane and δ¹³C of methane. *Atmos. Chem. Phys.*, 22, 15351–15377.

### 3. Other References Using δ¹³C-CH₄ Station Data

- **Chandra et al. (2024)** — Used δ¹³C-CH₄ data from NOAA/INSTAAR and TU/NIPR networks for atmospheric inversion. Supplementary data at https://doi.org/10.5281/zenodo.1053174994.
- **Thanwerdas et al. (2024)** — 3-D variational inverse modeling using δ¹³C-CH₄ observations from the NOAA/INSTAAR network.

## Station Inventory

The full station inventory for the copied NOAA/INSTAAR 2023 data is in **`d13C_sites_inventory.csv`** with columns:
`site_code, site_name, country, latitude, longitude, elevation_m, n_observations, year_start, year_end, data_source, file`

## Key Long-Record Sites (>15 years of data)

| Site | Location | Lat | Lon | Period | N obs |
|------|----------|-----|-----|--------|-------|
| ALT | Alert, Nunavut | 82.5 | -62.5 | 2000–2022 | 2,999 |
| BRW | Barrow, Alaska | 71.3 | -156.6 | 1998–2022 | 4,174 |
| MLO | Mauna Loa, Hawaii | 19.5 | -155.6 | 1998–2022 | 4,058 |
| NWR | Niwot Ridge, Colorado | 40.1 | -105.6 | 1998–2022 | 2,055 |
| SMO | Tutuila, American Samoa | -14.2 | -170.6 | 1998–2022 | 3,568 |
| CGO | Cape Grim, Tasmania | -40.7 | 144.7 | 1998–2022 | 1,639 |
| SPO | South Pole, Antarctica | -90.0 | -24.8 | 1998–2022 | 3,215 |
| KUM | Cape Kumukahi, Hawaii | 19.7 | -155.0 | 1999–2022 | 2,117 |
| MHD | Mace Head, Ireland | 53.3 | -9.9 | 1999–2022 | 1,532 |
| ASC | Ascension Island | -8.0 | -14.4 | 2001–2022 | 3,647 |

## Spatial Coverage

The δ¹³C-CH₄ network is sparser than CH₄ but covers key latitudinal bands:
- **High Arctic**: ALT (82°N), BRW (71°N), ZEP (79°N), SUM (73°N)
- **Northern midlatitudes**: MHD, NWR, CBA, TAP, WLG, AZR, AMY
- **Tropics**: ASC, SMO, KUM, MLO
- **Southern midlatitudes**: CGO, BHD, EIC
- **Antarctic**: SPO (90°S), PSA (65°S)

The Basu2022 ObsPack adds ~40 more stations from non-NOAA labs (MPI-BGC, NIWA, NIPR, TU).

## Labs Contributing Data

| Lab | Source | # Sites | Notes |
|-----|--------|---------|-------|
| NOAA/INSTAAR | NOAA 2023 (copied) | 25 | Primary δ¹³C network, longest records |
| NOAA/INSTAAR | Basu2022 (catalogued) | ~40 | Earlier data compilation |
| MPI-BGC | Basu2022 (catalogued) | ~8 | European & SH high-latitude |
| NIWA | Basu2022 (catalogued) | ~3 | New Zealand, Pacific |
| NIPR/TU | Basu2022 (catalogued) | ~5 | Japanese program, Arctic/Pacific |

## Notes

- δ¹³C values are reported in **‰ VPDB** (Vienna Pee Dee Belemnite)
- All data harmonized to the NOAA/INSTAAR scale
- Primary application: distinguishing microbial vs. fossil fuel vs. pyrogenic CH₄ sources
- The Basu2022 ObsPack catalogue CSV references the NC files in `ImportantReferences/Basu2022ACP_DS/nc/`

## Map

See `d13C_sites_map.png` for the spatial distribution of NOAA/INSTAAR 2023 sites.
