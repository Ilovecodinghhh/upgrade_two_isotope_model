# δD-CH₄ (Deuterium) — Site Data Inventory

## Overview

| Metric | Value |
|--------|-------|
| **Total sites** | 18 |
| **Northern Hemisphere** | 11 |
| **Southern Hemisphere** | 7 |
| **Earliest data** | 1988 |
| **Latest data** | 2024 |
| **Total approximate observations** | ~5,010 |

## ⚠️ Important Note on δD Data Availability

**δD-CH₄ measurements are far sparser than CH₄ and δ¹³C measurements.** Only four laboratories worldwide have produced marine boundary layer (MBL) δD-CH₄ time series:

1. **INSTAAR** (CU-Boulder / NOAA) — 10 MBL sites, 2005–2009
2. **MPI-BGC** (Max Planck, Jena) — 6 MBL sites, 2010–2024
3. **IMAU** (Utrecht University) — 5 MBL sites, 1988–2024 (longest record)
4. **TU/NIPR** (Tohoku U. / Nat. Inst. Polar Research) — 2 MBL sites, 1995–2023

These data have been recently harmonized by **Riddell-Young et al. (2025, PNAS)** and **Dasgupta et al. (2025, EGUsphere)** to produce the first comprehensive global δD-CH₄ record.

## Data Sources & References

1. **Riddell-Young et al. (2025)** — Primary reference for the harmonized global δD-CH₄ record.
   - *Citation:* Riddell-Young, B., et al. (2025). Microbial driver of 2006–2023 CH₄ growth indicated by trends in atmospheric δD-CH₄ and δ¹³C-CH₄. *PNAS*.
   - *Data package:* `Riddell-Young_2025_dD_GlobMean` (includes station-level MC data)
   - ⚠️ **Raw station-level δD data is NOT publicly available in a standard download portal** — it was compiled by the authors from individual lab datasets.

2. **Dasgupta et al. (2025)** — Harmonisation of methane isotope ratio measurements from different laboratories using atmospheric samples.
   - *Citation:* Dasgupta, B., et al. (2025). Harmonisation of methane isotope ratio measurements from different laboratories. *EGUsphere* [preprint]. https://doi.org/10.5194/egusphere-2025-2439
   - This paper provides the inter-laboratory calibration offsets used to harmonize δD data.

3. **Fujita et al. (2025)** — Used δD-CH₄ data in global CH₄ budget modeling.
   - *Data:* Supplementary Dataset S1 includes simulated global δD-CH₄ from posterior scenarios.

4. **Umezawa et al. (2018)** — Inter-laboratory comparison of δD measurements.
   - *Citation:* Umezawa, T., et al. (2018). *Atmos. Meas. Tech.*, 11, 1059–1078.

## Site Details

| Site | Location | Lat | Lon | Labs | Period | ~N obs |
|------|----------|-----|-----|------|--------|--------|
| ALT | Alert, Nunavut | 82.5 | -62.5 | INSTAAR; MPI-BGC | 2005–2024 | ~520 |
| NYA | Ny-Ålesund, Svalbard | 78.9 | 11.9 | IMAU; TU/NIPR | 1988–2024 | ~800 |
| ZEP | Ny-Ålesund, Svalbard (Zeppelin | 78.9 | 11.9 | MPI-BGC; IMAU | 2006–2024 | ~600 |
| BRW | Barrow (Utqiaġvik), Alaska | 71.3 | -156.6 | INSTAAR | 2005–2009 | ~120 |
| CBA | Cold Bay, Alaska | 55.2 | -162.7 | INSTAAR | 2005–2009 | ~80 |
| MHD | Mace Head, Ireland | 53.3 | -9.9 | INSTAAR; IMAU | 2005–2024 | ~300 |
| JFJ | Jungfraujoch, Switzerland | 46.5 | 8.0 | MPI-BGC | 2012–2017 | ~150 |
| NWR | Niwot Ridge, Colorado | 40.0 | -105.6 | INSTAAR | 2005–2009 | ~120 |
| KUM | Cape Kumukahi, Hawaii | 19.7 | -155.0 | INSTAAR | 2005–2009 | ~100 |
| MLO | Mauna Loa, Hawaii | 19.5 | -155.6 | INSTAAR | 2005–2009 | ~120 |
| CVO | Cape Verde Observatory | 16.9 | -24.9 | MPI-BGC | 2011–2017 | ~300 |
| SMO | Tutuila, American Samoa | -14.2 | -170.6 | INSTAAR | 2005–2009 | ~100 |
| NMB | Gobabeb, Namibia | -23.6 | 15.0 | MPI-BGC | 2013–2017 | ~200 |
| AMS | Amsterdam Island | -37.8 | 77.5 | IMAU | 2006–2024 | ~200 |
| CGO | Cape Grim, Tasmania | -40.7 | 144.7 | INSTAAR | 2005–2009 | ~100 |
| SYO | Syowa Station, Antarctica | -69.0 | 39.6 | TU/NIPR | 1995–2023 | ~500 |
| GVN | Neumayer Station, Antarctica | -70.7 | -8.2 | MPI-BGC; IMAU | 2006–2024 | ~600 |
| SPO | South Pole, Antarctica | -90.0 | -24.8 | INSTAAR | 2005–2009 | ~100 |

## Inter-Laboratory Calibration

All δD data are harmonized to the **MPI scale** using offsets from Umezawa et al. (2018) and Dasgupta et al. (2025):

| Lab | MPI Offset (Umezawa) | MPI Offset (Dasgupta) |
|-----|---------------------|----------------------|
| INSTAAR | −11.5 ± 1.5‰ | −10.9 ± 2.2‰ |
| MPI-BGC | = (reference) | = (reference) |
| IMAU | +4.2 ± 1.2‰ | +2.4 ± 1.6‰ |
| TU/NIPR | −8.9 ± 1.3‰ | −10.8 ± 1.6‰ |

## Spatial Coverage

The δD network is heavily weighted toward **high latitudes**:
- **Arctic**: ALT, BRW, NYA/ZEP (all >70°N)
- **Northern midlatitudes**: MHD, NWR, CBA, JFJ
- **Subtropics**: MLO, KUM, CVO
- **Southern midlatitudes**: CGO, AMS
- **Antarctic**: SPO, GVN, SYO

**⚠️ Gap: Very few tropical δD-CH₄ sites exist**, which is a key limitation for understanding tropical methane sources.

## Notes

- δD values are reported in **‰ VSMOW** (Vienna Standard Mean Ocean Water)
- δD is sensitive to both source water isotopic composition and atmospheric oxidation (KIE)
- Complementary to δ¹³C for distinguishing microbial vs. thermogenic sources
- **Raw station data availability**: Contact individual labs (MPI, IMAU, NIPR, INSTAAR) for access
- The harmonized global mean time series is available in the repository under `rel/data/`

## Map

See `dD_sites_map.png` for the spatial distribution of all sites.
