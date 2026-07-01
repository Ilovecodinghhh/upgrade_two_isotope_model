# Session Summary — 2026-05-20

## Overview

Two tasks completed today:
1. **Processed the Li et al. (2026) ESSD wetland emission dataset** — exploratory data analysis of the newly uploaded `Li2026ESSD_DS.nc` NetCDF file to characterise wetland CH₄ seasonality at latitude bands relevant to our 12 KIE sites.
2. **Read and analysed Liu et al. (2025) Nature** — extracted methodological insights on how a state-of-the-art study handled the source–sink phase problem we face.

---

## 1. Li2026 ESSD Dataset — Exploratory Processing

### Dataset structure
- **File**: `ImportantReferences/Li2026ESSD_DS.nc` (20 MB)
- **Dimensions**: time = 312 months (Jan 2000 – Dec 2025), lat = 180 (−89.5 to 89.5°), lon = 360 (−179.5 to 179.5°)
- **Variable**: `wetch4` — monthly wetland CH₄ emissions (kg CH₄ per 1°×1° cell per month)
- **Source**: 35 Global Methane Budget model estimates × 10 ERA5 ensemble members, emulated with XGBoost

### 2005–2010 Climatology (matching INSTAAR measurement period)

| Latitude Band | Annual (Tg/yr) | Peak Month | Fractional Amplitude |
|---------------|---------------|------------|---------------------|
| NH high (60–90°N) | 10.9 | July | 3.02 (extreme) |
| NH mid (30–60°N) | 29.7 | July | 2.00 (strong) |
| Tropics (30°S–30°N) | 114.0 | August | 0.21 (nearly flat) |
| SH extra (90°S–30°S) | 2.9 | January | 0.63 (moderate) |
| **Global** | **157.5** | — | — |

Monthly values (Tg/month) for NH high-lat band:
```
Jan=0.08  Feb=0.07  Mar=0.09  Apr=0.19  May=0.93  Jun=1.98
Jul=2.81  Aug=2.42  Sep=1.47  Oct=0.58  Nov=0.18  Dec=0.10
```

### Point-extraction tests
- **ALT, SPO**: grid cells are NaN (no local wetlands) — confirms need for latitude-band aggregation
- **BRW**: strong signal peaking Jul–Aug, consistent with North Slope tundra wetlands
- **CGO**: modest year-round emissions peaking January (SH summer)

### Source region assignment per site

| Site | Lat | Assigned Source Region | Rationale |
|------|-----|-----------------------|-----------|
| ALT | +82.5° | NH high-lat (60–90°N) | Arctic; no local wetlands |
| ZEP | +78.9° | NH high-lat (60–90°N) | Svalbard; Eurasian boreal sources |
| BRW | +71.3° | NH high-lat (60–90°N) | North Slope; representative |
| CBA | +55.2° | NH mid-lat (30–60°N) | Downstream of boreal wetlands |
| MHD | +53.3° | NH mid-lat (30–60°N) | Atlantic; continental sources |
| AZR | +38.8° | NH mid-lat (30–60°N) | (Excluded site) |
| MLO | +19.5° | Tropics (30°S–30°N) | (Excluded site) |
| KUM | +19.6° | Tropics (30°S–30°N) | Tropical sources |
| ASC | −8.0° | Tropics (30°S–30°N) | (Excluded site) |
| SMO | −14.2° | Tropics (30°S–30°N) | (Excluded site) |
| CGO | −40.7° | Global weighted | SH sees interhemispheric transport |
| SPO | −90.0° | Global weighted | ~1 yr lag; global mean |

---

## 2. Liu et al. (2025) Nature — Phase Problem Analysis

### What they did
- Analysed CH₄ seasonal cycle amplitude (SCA) trends at 27 NOAA sites over 1984–2020
- Used GEOS-Chem v14.1.0 (full 3D CTM) with MERRA-2 meteorology
- Factorial experiments: T1 (all vary), T2 (emissions+sinks fixed at 1984), T3 (emissions vary, sinks fixed)
- SSA decomposition: SCA = SCA_E − SCA_S + SCA_T (Eqs 2–8)

### Key results
- BRW: SCA decreasing at −0.35 ppb/yr (wetland emission increase → smooths seasonal cycle)
- MLO, SPO: SCA increasing at +0.19, +0.18 ppb/yr (OH increase → steeper seasonal destruction)
- 10 ± 1% increase in tropospheric OH since 1984; 21 ± 1% increase in CH₄ sink
- **Transport effect (IAT_E)**: 35% of SPO SCA increase comes from transported emissions, not local chemistry

### The phase problem: did they face it?
**Yes, implicitly — but they avoided it by using a full 3D CTM.** The CTM inherently handles phase offsets between OH (peaks ~July) and wetland emissions (peaks ~Aug–Sep) because transport, chemistry, and emissions evolve simultaneously within the model at each time step. Their SSA decomposition then evaluates all components at the same peak/trough months.

### Implications for our work
1. **Validates our phasor approach**: Our vector (phasor) subtraction in complex B+iC space is the analytical equivalent of what the CTM does implicitly — necessary for us because no isotope-enabled CTM is available for this analysis.
2. **Validates CGO/SPO source region**: Liu's finding that 35% of SPO SCA change is from transported emissions confirms our decision to use a global-weighted (not local) source region for CGO and SPO.
3. **Phase offset is real**: Liu's wetland emissions (ORCHIDEE) peak 1–2 months after OH peak — simple amplitude subtraction would be wrong; vector decomposition is necessary.

---

## 3. Prior Work (completed in earlier sessions, referenced today)

### Site-specific δD source signatures (completed ~May 18–19)
- Built `data/dD_source_database.json` using OIPC precipitation δ²H + Douglas (2021) regressions
- Two methods compared: OIPC regression (site-specific) vs Douglas zonal means (emission-weighted)
- Recommended values per site (biggest shifts from old global −310‰: NH high-lat move −64 to −74‰)
- Script: `analysis/build_dD_source_db.py`; Figure: `fig6_dD_source_vs_latitude.png`

### Phases 1–5 (completed in earlier sessions)
| Phase | Script | Key Output |
|-------|--------|------------|
| Phase 1 | `phase1_data.py` | 12 co-located sites, monthly CSV files |
| Phase 2 | `phase2_harmonics.py` | Harmonic fits → B, C, amplitude, phase per site per isotope |
| Phase 3 | `phase3_synthesis.py` | Site classification (8 clean, 4 excluded); latitude gradient in R |
| Phase 4 | `phase4_deconv.py` | Source deconvolution; sink-only ratio is ~constant |
| Phase 5 | `phase5_kie.py` | **α¹³C_OH = 1.0034 [0.998, 1.009]** (SH direct); 1.0162 (source-corrected, biased high) |

### Current KIE results (from Phase 5)
| Approach | α¹³C_OH | 95% CI | Notes |
|----------|---------|--------|-------|
| SH direct (CGO+SPO) | 1.0034 | [0.998, 1.009] | Most robust |
| Source-corrected (all clean) | 1.0162 | [1.003, 1.031] | Biased high by wetland contamination |
| Saueressig (2001) lab | 1.0039 | — | |
| Cantrell (1990) lab | 1.0054 | — | |

---

## 4. Discussion.md (written ~May 19, pushed today)

Documented the full strategy for wetland source correction:
- **Problem**: Asymmetric contamination (source/sink ratio 2.5× for δ¹³C vs 0.8× for δD) inflates R
- **Solution**: Phasor decomposition: Z_sink = Z_obs − Z_source in complex B+iC space
- **Consistency check**: After subtraction, sink phases for both isotopes should agree (~July)
- **Data needs**: wetland flux seasonality (Li2026), δD source signatures (done), δ¹³C source (−62‰ adequate)

---

## Files Created/Modified Today

| File | Action | Status |
|------|--------|--------|
| `experiments/KIE_sites/Discussion.md` | Pushed to GitHub | ✅ |
| Li2026ESSD_DS.nc explored | Read & characterised | ✅ |
| Liu2025Nature & SI | Read & analysed | ✅ |
| Plan file (binary-conjuring-flurry.md) | Updated with full implementation plan | ✅ |
| No new analysis scripts today | — | Planned for next session |
