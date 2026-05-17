# Methane Mixing Ratio (CH₄ ppb) — Site Data Inventory

## Overview

| Metric | Value |
|--------|-------|
| **Total sites** | 233 |
| **Northern Hemisphere (>0°N)** | 191 |
| **Southern Hemisphere (≤0°S)** | 42 |
| **Tropical (30°S–30°N)** | 52 |
| **Long records (>10 yr)** | 124 |
| **Earliest data** | 1983 |
| **Latest data** | 2024 |
| **Total observations** | 7,651,297 |

## Data Sources

1. **Basu et al. (2022) ObsPack** — Multi-lab compilation harmonizing NOAA/GML, MPI-BGC, NIWA, NIPR, LSCE, EC, KMA, NILU, EMPA, IPEN, RSE, TU, and others. 293 CH₄ files covering 235 unique stations with data through ~2020.
   - *Citation:* Basu, S., et al. (2022). Estimating emissions of methane consistent with atmospheric measurements of methane and δ13C of methane. *Atmos. Chem. Phys.*, 22, 15351–15377.
   - *Data:* ObsPack multi-species (obspack_multi-species_1_ch4_d13ch4_2022-10-17)

2. **NOAA GML Flask Network (2025)** — Latest release (2025-08-15) of discrete surface flask CH₄ from the NOAA Global Greenhouse Gas Reference Network. 91 surface stations, data through 2024/2025.
   - *Citation:* Lan, X., et al. (2025). Atmospheric Methane Dry Air Mole Fractions from the NOAA GML CCGG Network, 1983–Present. Version 2025-08-15. https://doi.org/10.15138/VNCZ-M766

## Spatial Coverage

The network provides excellent global coverage:
- **Arctic/boreal**: ALT (82°N), BRW (71°N), ZEP (79°N), SUM (73°N), NYA (79°N)
- **Northern midlatitudes**: Dense coverage across North America, Europe, and East Asia
- **Tropics**: ASC, BKT, CHR, GMI, KUM, MLO, SMO, RPB
- **Southern midlatitudes**: CGO, BHD, AMS, CRZ, MQA
- **Antarctic**: SPO (90°S), PSA (65°S), SYO (69°S), HBA (76°S)

## Timespan Distribution

| Period | # of sites with data |
|--------|---------------------|
| Pre-2000 | 73 |
| 2000–2010 | 169 |
| 2010–2020 | 193 |
| Post-2020 | 23 |

## Labs Contributing Data

| Lab | Sites |
|-----|-------|
| NOAA/GML | Primary network, ~90 surface sites |
| MPI-BGC | European and Southern Hemisphere stations |
| NIWA | New Zealand, Pacific sites |
| NIPR/TU | Japanese Arctic/Antarctic sites |
| LSCE | French sites and islands |
| EC (Environment Canada) | Canadian sites |
| Others | KMA, NILU, EMPA, IPEN, RSE, etc. |

## Notes

- CH₄ values are in **nmol/mol (ppb)** on the **WMO X2004A** scale
- Data includes both event (discrete flask) and hourly (continuous in-situ) measurements
- For the most up-to-date data, access NOAA GML: https://gml.noaa.gov/dv/data/
- The full inventory CSV is in `ch4_sites_inventory.csv`

## Map

See `ch4_sites_map.png` for the spatial distribution of all sites.
