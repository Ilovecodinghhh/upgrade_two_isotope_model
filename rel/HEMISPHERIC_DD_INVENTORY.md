# Hemispheric δD Data Inventory

**Date:** 2026-05-12
**Purpose:** Document all available δD data (atmospheric + source signatures) for hemispheric (NH/SH) modeling.

---

## 1. Atmospheric δD-CH₄ — Hemispheric MC Iterations

### What we built (new, 2026-05-12)

| File | Location | Description |
|------|----------|-------------|
| `NHMean_dD_iterations_DasguptaCal_noBUDS.xlsx` | `rel/data/` | NH annual δD, 1000 MC iterations, 2005–2023 |
| `SHMean_dD_iterations_DasguptaCal_noBUDS.xlsx` | `rel/data/` | SH annual δD, 1000 MC iterations, 2005–2023 |
| `HemMean_dD_annual_DasguptaCal_noBUDS.csv` | `rel/data/` | NH/SH annual means + std + NH–SH gradient |
| `HemMean_dD_dei_DasguptaCal_noBUDS.csv` | `rel/data/` | NH/SH/Global weekly means (from Riddell-Young pipeline) |
| `GlobMean_dD_iterations_DasguptaCal_noBUDS.xlsx` | `rel/data/` | Global annual δD, 1000 MC iterations (recalculated) |
| `GlobMean_dD_dei_DasguptaCal_noBUDS.csv` | `rel/data/` | Global annual mean (recalculated) |
| `SemiHemMean_dD_dei_DasguptaCal_noBUDS.csv` | `rel/data/` | PN/TN/TS/PS weekly means (4 latitude bands) |

**Pipeline:** Riddell-Young 2025 `MBL_calc_Unc.py` (per-station smoothed MC) → `dD_hemispheric_MC.py` (our script, computes NH/SH annual MC iterations from station MC output).

**Station coverage:**
- 16 NH stations (8 NOAA, 3 MPI, 2 Tohoku, 3 IMAU) covering 83°N to 1°N
- 7 SH stations covering 14°S to 90°S

**Observed NH–SH gradient:** −12 to −16‰ (NH more depleted), consistent with NH fossil-fuel source dominance.

**Year range:** 2005–2019 have full coverage. 2020–2023 have partial NaN (station gaps from COVID + reporting lag). 2005 is the first year with adequate SH station coverage.

**Calibration:** Dasgupta (2025) calibration scale, no BUDS measurements.

### What existed before (the hack in `common.py`)

```python
DD_IH_OFFSET = 6.0  # ‰
dD_NH = dD_global − 6
dD_SH = dD_global + 6
```

This constant offset applied to a single global mean time series. No year-to-year hemispheric variability. Now replaced by actual observations.

---

## 2. Source Signature δD — Sector-by-Sector Inventory

### Fossil Fuel (FF)

| Product | Location | What it is | Hemispheric? |
|---------|----------|------------|:------------:|
| `FF_dD_GlobUnc.csv` | `rel/output/` | Global mean FF δD time series (annual, 24 years, with ONG/Coal/Shale breakdown) | ❌ Global only |
| `FF_dD_GlobMC_EDGAR.csv` | `rel/output/` | Global FF δD, 1000 MC iterations (EDGAR-weighted) | ❌ Global only |
| `FF_dD_GlobMC_CTCH4.csv` | `rel/output/` | Global FF δD, 1000 MC iterations (CarbonTracker-weighted) | ❌ Global only |
| `FF_dD_1x1.txt` | Riddell-Young `Output/` | **1°×1° gridded FF δD map** (time-averaged) | ✅ Splittable into NH/SH |
| `FF_dD_1x1_{year}.txt` | Riddell-Young `Output/` | **1°×1° annual FF δD maps** (per year, 2000–2023) | ✅ Splittable into NH/SH |
| `ONG_dD.csv`, `coal_dD.csv` | Riddell-Young `data/` | Raw per-country ONG and coal δD values | ✅ Country-level → can aggregate by hemisphere |

**Methodology:** Country-level ONG and coal δD (from Sherwood et al. 2017 + new compilations) → EDGAR 8.0 emission-weighted gridding → 1°×1° maps. Source: `FF_dD_map_EDGAR.py`.

**Key values:** Global mean FF δD ≈ −175 to −185‰ (varies by year as shale gas fraction changes).

### Microbial (Mic)

| Product | Location | What it is | Hemispheric? |
|---------|----------|------------|:------------:|
| `Mic_dD_AnnGlob.csv` | `rel/output/` | Global mean Mic δD (annual, 24 years, wetland/ruminant/rice/landfill/termite breakdown) | ❌ Global only |
| `Mic_dD_MC.csv` | `rel/output/` | Global Mic δD, 1000 MC iterations | ❌ Global only |
| `dD_90N` grid | Riddell-Young `Output/` | **180×360 gridded Mic δD map** (emission-weighted, from Douglas 2021 MAT regression) | ✅ Splittable into NH/SH |
| `std_dev_matrix_mic_90N` | Riddell-Young `Output/` | Mic δD uncertainty grid (1°×1°) | ✅ Splittable into NH/SH |

**Methodology:** δD_Mic = 0.6088 × MAT − 285.7 (Douglas et al. 2021 regression of methane δD on mean annual temperature). Applied per grid cell, then weighted by subcategory emissions (wetland, rice, ruminant, landfill, termite, wild animal) from EDGAR + CarbonTracker posteriors. Source: `Mic-BB_dD_Map_new.py`, `Mic_dD.py`.

**Key values:** Global mean Mic δD ≈ −310 ± 30‰.

### Biomass Burning (BB)

| Product | Location | What it is | Hemispheric? |
|---------|----------|------------|:------------:|
| `BB_dD_annual.csv` | `rel/output/` | Global mean BB δD (annual, 24 years) | ❌ Global only |
| `dDBB_90N` grid | Riddell-Young `Output/` | **180×360 gridded BB δD map** (from Umezawa 2011 MAT regression) | ✅ Splittable into NH/SH |
| `std_dev_matrix_90N` (BB) | Riddell-Young `Output/` | BB δD uncertainty grid (1°×1°) | ✅ Splittable into NH/SH |

**Methodology:** δD_BB = 1.16 × MAT − 177 (Umezawa et al. 2011 regression). Applied per grid cell, weighted by GFEDv4s fire emissions. Source: `Mic-BB_dD_Map_new.py`.

**Key values:** Global mean BB δD ≈ −210 to −230‰.

---

## 3. Summary: Can We Run a Proper 2-Box Agreement Filter?

| Requirement | Status | Source |
|-------------|--------|--------|
| NH/SH atmospheric δD (MC iterations) | ✅ **Done** | `rel/data/NHMean_dD_iterations_*` (built today) |
| NH/SH atmospheric δ¹³C (MC iterations) | ✅ **Already exists** | `rel/data/ch4c13_nh_sh_mean.xlsx` (loaded in `common.py`) |
| NH/SH FF δD source signature | ✅ **Available** | `FF_dD_1x1_{year}.txt` grids → split by latitude |
| NH/SH Mic δD source signature | ✅ **Available** | `dD_90N` grid → split by latitude |
| NH/SH BB δD source signature | ✅ **Available** | `dDBB_90N` grid → split by latitude |
| NH/SH FF δ¹³C source signature | ✅ **Already exists** | in `common.py` hemispheric config |
| NH/SH Mic δ¹³C source signature | ✅ **Already exists** | in `common.py` hemispheric config |
| NH/SH BB δ¹³C source signature | ✅ **Already exists** | in `common.py` hemispheric config |

**Answer: YES — all data exist to run a proper 2-box (NH/SH) dual-isotope agreement filter with real hemispheric observations and gridded source signatures, replacing the ±6‰ offset hack.**

The only processing needed is to latitude-weight the 1°×1° source signature grids into NH (0°–90°N) and SH (90°S–0°) means, accounting for emission weighting. This is straightforward from the existing gridded outputs.

---

## 4. Data Provenance

| Dataset | Reference |
|---------|-----------|
| Station δD-CH₄ measurements | NOAA GML, MPI, Tohoku/NIPR, IMAU networks |
| Global mean harmonization | Riddell-Young et al. (2025) PNAS; Umezawa (2012) calibration |
| FF δD (ONG, coal) | Sherwood et al. (2017); new ONG compilation in Riddell-Young (2025) |
| Mic δD (MAT regression) | Douglas et al. (2021) |
| BB δD (MAT regression) | Umezawa et al. (2011) |
| Emission grids | EDGAR 8.0 (ONG, coal); CarbonTracker CH4 v2025; GFEDv4s (fires) |
| Hemispheric MC pipeline | `dD_globmean.py` + `dD_hemispheric_MC.py` (this work) |
