# δ¹³C-CH₄ — Site Data Inventory

## Overview

| Metric | Value |
|--------|-------|
| **Total sites** | 63 |
| **Northern Hemisphere (>0°N)** | 52 |
| **Southern Hemisphere (≤0°S)** | 11 |
| **Long records (>10 yr)** | 21 |
| **Earliest data** | 1998 |
| **Latest data** | 2022 |
| **Total observations** | 102,161 |

## Data Sources

1. **Basu et al. (2022) ObsPack** — 69 δ¹³C-CH₄ files covering 63 unique stations. Multi-lab compilation (NOAA/INSTAAR, MPI-BGC, NIWA, NIPR, TU) with data through ~2020.
   - *Citation:* Basu, S., et al. (2022). *Atmos. Chem. Phys.*, 22, 15351–15377.

2. **NOAA/INSTAAR Flask Network (2023)** — Latest release (2023-09-21) of δ¹³C-CH₄ from the NOAA GML / INSTAAR Carbon Cycle Cooperative Global Air Sampling Network. 25 surface stations, data 1998–2022.
   - *Citation:* Michel, S.E., et al. (2023). INSTAAR Stable Isotopic Composition of Atmospheric Methane (¹³C), 1998–2022. https://doi.org/10.15138/9p89-1x02

3. **Chandra et al. (2024)** — Used δ¹³C-CH₄ data from NOAA/INSTAAR and TU/NIPR networks for atmospheric inversion. Supplementary data at https://doi.org/10.5281/zenodo.1053174994.

4. **Thanwerdas et al. (2024)** — 3-D variational inverse modeling using δ¹³C-CH₄ observations from the NOAA/INSTAAR network.

## Key Long-Record Sites (>15 years of data)

| Site | Location | Lat | Lon | Period | N obs |
|------|----------|-----|-----|--------|-------|
| comment | Ny-Alesund, Svalbard | 78.9 | 11.9 | 1998–2022 | 46,251 |
| BRW | Barrow Atmospheric Baseline Ob | 71.3 | -156.6 | 1998–2020 | 4,174 |
| MLO | Mauna Loa, Hawaii | 19.5 | -155.6 | 1998–2020 | 4,058 |
| NWR | Niwot Ridge, Colorado | 40.1 | -105.6 | 1998–2020 | 2,055 |
| SMO | Tutuila | -14.2 | -170.6 | 1998–2020 | 3,568 |
| CGO | Cape Grim, Tasmania | -40.7 | 144.7 | 1998–2020 | 1,639 |
| SPO | South Pole, Antarctica | -90.0 | -24.8 | 1998–2020 | 3,215 |
| KUM | Cape Kumukahi, Hawaii | 19.7 | -155.0 | 1999–2020 | 2,117 |
| MHD | Mace Head, County Galway | 53.3 | -9.9 | 1999–2020 | 1,532 |
| NYA | Ny-Alesund, Svalbard | 78.9 | 11.9 | 2000–2020 | 947 |
| BHD | Baring Head Station | -41.4 | 174.9 | 2000–2020 | 561 |
| ALT | Alert, Nunavut | 82.5 | -62.5 | 2001–2020 | 3,016 |
| CBA | Cold Bay, Alaska | 55.2 | -162.7 | 2001–2020 | 1,897 |
| ASC | Ascension Island | -8.0 | -14.4 | 2001–2020 | 3,174 |
| AZR | Terceira Island, Azores | 38.8 | -27.4 | 2001–2020 | 939 |
| TAP | Tae-ahn Peninsula | 36.7 | 126.1 | 2001–2020 | 1,406 |
| ZEP | Ny-Alesund, Svalbard | 78.9 | 11.9 | 2002–2020 | 1,168 |
| WLG | Mt. Waliguan | 36.3 | 100.9 | 2002–2020 | 1,948 |
| ARH | Arrival Heights, Antarctica | -77.8 | 166.7 | 2000–2016 | 308 |

## Spatial Coverage

The δ¹³C-CH₄ network is sparser than CH₄ but covers key latitudinal bands:
- **High Arctic**: ALT (82°N), BRW (71°N), ZEP (79°N), NYA (79°N), SUM (73°N)
- **Northern midlatitudes**: MHD, NWR, MLO, CBA, TAP, WLG, KUM
- **Tropics**: ASC, SMO, KUM, MLO
- **Southern midlatitudes**: CGO, BHD, EIC
- **Antarctic**: SPO (90°S)

## Labs Contributing Data

| Lab | # Sites | Notes |
|-----|---------|-------|
| NOAA/INSTAAR | ~40+ | Primary δ¹³C network, longest records |
| MPI-BGC | ~8 | European & SH high-latitude |
| NIWA | ~3 | New Zealand, Pacific |
| NIPR/TU | ~5 | Japanese program, Arctic/Pacific |

## Notes

- δ¹³C values are reported in **‰ VPDB** (Vienna Pee Dee Belemnite)
- All data harmonized to the NOAA/INSTAAR scale
- Primary application: distinguishing microbial vs. fossil fuel vs. pyrogenic CH₄ sources
- The full inventory CSV is in `d13C_sites_inventory.csv`

## Map

See `d13C_sites_map.png` for the spatial distribution of all sites.
