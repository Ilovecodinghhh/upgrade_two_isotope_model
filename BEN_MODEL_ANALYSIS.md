# Ben-BoxModel (Riddell-Young et al. 2025) — Analysis & Relevance

## What This Is

**Publication**: Ben Riddell-Young et al. (2025, GRL) — "Microbial driver of 2006–2023 CH₄ growth indicated by trends"

**Code package**: `Riddell-Young_2025_MassBalancePackage` — the complete one-box dual-isotope (δ¹³C + δD) Monte Carlo box model, including all upstream preprocessing scripts for source endmember time series.

This is essentially the **published version** of what your `TwoIsotopeBoxModel` is based on, but with several important methodological improvements and a complete, reproducible data pipeline.

---

## Architecture Comparison: Ben vs. Our Models

| Feature | Ben (Riddell-Young 2025) | Your Original Model | Our v3.0 Two-Hemi |
|---------|--------------------------|--------------------|--------------------|
| Spatial | Global one-box | Global one-box | NH/SH two-box |
| Isotopes | δ¹³C and δD (separate inversions) | δ¹³C + δD (combined) | δ¹³C + δD (combined, per hemi) |
| Source categories | FF, Mic, BB | FF, Mic, BB | FF, Mic, BB per hemisphere |
| MC iterations | 1000 | 1000 | 1000 |
| KIE treatment | **Fixed** per scenario | Fixed (v1), **Sampled** (v2+) | **Sampled** (from v2) |
| Source sigs | **Time-varying + MC** (key advance) | Time-varying + MC | Time-varying + MC |
| BB constraint | **Fixed from CT mean** (important!) | Fixed from CT mean | Freed (solved from isotopes) |
| Solver | **2 equations (mass + one isotope)** | 3×3 (mass + δ¹³C + δD) | 3×3 bounded lsq |
| Post-processing | **5-year smoothing** | None | None |
| Time extent | 1999–2022 (d13C); 2005–2022 (dD) | 1999–2022 | 1999–2021 |

### CRITICAL DIFFERENCE: Ben does NOT combine δ¹³C and δD

Ben's model runs **two separate inversions**:
1. `d13C_MassBalance_MC.py`: Mass + δ¹³C → solve for FF & Mic (BB fixed)
2. `dD_MassBalance_MC.py`: Mass + δD → solve for FF & Mic (BB fixed)

Each is a **2-equation system** (not 3×3), which is **much better conditioned** than the 3×3 system. The tradeoff is that BB emissions are not independently constrained by isotopes — they're prescribed from CarbonTracker prior.

Our approach (combining both isotopes in a 3×3 system) is more ambitious but suffers from the ill-conditioning documented in BOX_MODEL_ASSESSMENT.md.

---

## Key Things We Can Learn From Ben's Model

### 1. **BB Should Be Fixed, Not Solved** (Most Important)

Ben **fixes BB emissions from CarbonTracker** (BB ≈ 30 Tg/yr constant):
```python
BB = np.mean(bbCT)  # ~30 Tg/yr from GFED4/CT-CH4 prior
```

Then solves for FF and Mic using only 2 equations:
```python
FFS = (SumSource * d13C_source - mic_d13C * (SumSource - BB) - bb_d13C * BB) / (ff_d13C - mic_d13C)
MicS = SumSource - BB - FFS
```

**Why this matters**: BB is the smallest source (~30 Tg/yr) with the most uncertain isotopic signature. Trying to solve for it from isotopes alone (our 3×3 approach) makes the system ill-conditioned. Ben's approach of fixing BB and solving only for FF/Mic gives a robust 2×2 system.

**Recommendation for v3.0**: Consider switching to Ben's approach — fix BB from GFED/CT-CH4 per hemisphere, solve 2×2 for FF and Mic per hemisphere. This would:
- Eliminate the ill-conditioning problem entirely
- Remove the need for bounded least squares
- Give much tighter uncertainty bounds

### 2. **Source Signature Preprocessing Pipeline** (Directly Reusable!)

Ben's `Output/` directory contains the exact same files your model reads! The data pipeline is:

| Script | What it generates | Key inputs |
|--------|------------------|------------|
| `Mic_d13C.py` | `Mic_d13C_annual.csv`, `Mic_d13C_MC.csv` | Oh 2022 wetlands + Chang 2019 ruminants + CT-CH4 posterior fractions |
| `Mic_dD.py` | `Mic_dD_AnnGlob.csv`, `Mic_dD_MC.csv` | d2H water map + CT-CH4 posterior fractions |
| `BB_d13C.py` | `BB_d13C_annual.csv` | C3/C4 distribution + GFED5/CT-CH4 fire maps |
| `FF_d13C_GlobMean.py` | `FF_d13C_GlobUnc.csv`, `FF_d13C_GlobMC_EDGAR.csv` | Coal/ONG isotope databases + EDGAR8 |
| `FF_dD_GlobMean.py` | `FF_dD_GlobUnc.csv`, `FF_dD_GlobMC_EDGAR.csv` | Coal/ONG δD databases + EDGAR8 |

**These outputs ARE the data your model (and our v3.0) already uses.** Ben's package is the upstream generator.

### 3. **Sensitivity Scenarios** (Framework for Robustness Testing)

Ben tests several KIE/sink scenarios that we should also test:
- Reduced Cl sink (gradual reduction to 0.011)
- Increased OH sink (+0.3%/yr following Olaf et al.)
- Reduced BB emissions
- Stable source signatures (for comparison)
- CT-CH4 weighted vs EDGAR weighted FF signatures

We currently only have the "base case." Adding these scenarios to our two-hemisphere model would strengthen the paper.

### 4. **5-Year Smoothing** (Post-Processing Best Practice)

Ben applies 5-year moving average smoothing to all MC iterations **before** computing mean/std:
```python
# 5 year smooth of each MC iteration, then statistics
```

This is important because:
- Year-to-year variability in the inversion is partly noise from uncertain inputs
- The signal (decadal trends) emerges from smoothing
- Without smoothing, you get the "noisy" impression that individual years are meaningful

**Recommendation**: Add 5-year smoothing option to our v3.0 output.

### 5. **Updated KIE Values** (Important for Our Model)

Ben uses the **Cantrell (1990) OH KIE = 1.0054** as their primary value (with Saueressig 1.0039 as a sensitivity test). Our v2.0/v3.0 samples uniformly between these. Ben's code documents:

```python
OH_KIE = 1.0054    # Cantrell — "Best estimate in Ben Li's opinion"
OH_KIE_D = 1.294   # Saueressig (δD); also uses Whitehill-Joelson Avg = 1.327 as sensitivity
Cl_KIE = 1.066 (13C), 1.52 (D)
Strat_KIE = 1.003 (13C), 1.179 (D)  # Dyonisius 2020; Beck 2018
Soil_KIE = 1.0201 (13C), 1.083 (D)
```

Our sampling ranges are consistent (OH_13C: [1.0039, 1.0054], OH_D: [1.294, 1.327]).

### 6. **CarbonTracker v2025** (Newer Data Available!)

Ben's package includes `carbontracker_CH4_v2025.xlsx` and `ch4_life_time_CT2025.txt` — this is the latest CT-CH4 version. Our model only uses the older CarbonTracker xlsx. The v2025 data may have:
- Updated posterior fluxes
- Revised lifetimes per year (could replace our linear τ(t) approximation!)

### 7. **Global δD Processing Sub-Package** (`Riddell-Young_2025_dD_GlobMean`)

This is the most complete available derivation of the global mean δD-CH₄ time series from station data. It:
- Harmonizes multiple labs (MPI/IMAU/NIPR/INSTAAR) 
- Uses MC station-dropout for network sampling uncertainty
- Produces `GlobMean_dD_iterations_UmezawaCal_noBUDS.xlsx` which is exactly what our model reads

**For the two-hemisphere model**: This sub-package has station-level data with latitude info. We could extend it to produce **NH/SH δD means** directly (rather than our ±1.5‰ approximation). The stations span Alert (82°N) to Scott Base (-77°S).

---

## What to Extract/Use Immediately

### Data files we should use from Ben's package (newer/more complete than what's in TwoIsotopeBoxModel):

1. **`data/ch4_life_time_CT2025.txt`** — Replace our linear τ(t) approximation with actual CT-CH4-derived lifetimes
2. **`data/carbontracker_CH4_v2025.xlsx`** — Updated emissions by category
3. **`data/EDGAR8_Coal.csv` + `data/EDGAR8_ONG.csv`** — For hemisphere-specific FF signatures (split by country lat)
4. **`Output/` CSVs** — These may be slightly updated compared to what's in `TwoIsotopeBoxModel/rel/output`
5. **Station-level δD data** in `Riddell-Young_2025_dD_GlobMean/data/` — For deriving NH/SH δD

### Methodological improvements to adopt:

1. **Fix BB from GFED/CT-CH4** instead of solving freely → eliminates ill-conditioning
2. **Add 5-year smoothing** as post-processing option
3. **Add sensitivity scenarios** (Cl reduction, OH increase, stable sigs)
4. **Use CT-2025 lifetimes** instead of linear fit

---

## Comparison of Results

Ben's main result (from the paper title): **Microbial emissions are the primary driver of 2006–2023 CH₄ growth.**

Our v3.0 global results:
- FF: 189.9 Tg/yr, Mic: 371.1 Tg/yr, BB: 24.5 Tg/yr

Ben's typical results (from code examination):
- FF: ~100-150 Tg/yr (smoothed), Mic: ~380-430 Tg/yr (smoothed)

Our FF is higher and Mic lower than Ben's. This is likely because:
1. We solve for BB freely (getting only 24.5 Tg/yr vs Ben's fixed 30 Tg/yr — not huge)
2. Our bounded least squares biases toward some configurations
3. Different KIE treatment (fixed vs sampled)
4. Two-hemisphere structure changes the result slightly

---

## Recommended Next Steps for v3.1

**Priority 1**: Adopt Ben's BB-fixed approach → Switch to 2-equation system per hemisphere:
```
FF_NH = (S_NH × f13_source_NH - mic_sig × (S_NH - BB_NH) - bb_sig × BB_NH) / (ff_sig - mic_sig)
Mic_NH = S_NH - BB_NH - FF_NH
```
This eliminates all ill-conditioning issues and the need for bounded least squares.

**Priority 2**: Use Ben's δD station data to derive actual NH/SH δD means (replace ±1.5‰ approximation).

**Priority 3**: Add CT-2025 lifetime data instead of linear τ(t).

**Priority 4**: Add 5-year smoothing + sensitivity scenarios.

---

## Key Numerical Values from the Paper (PNAS 2025, Vol. 122, No. 50)

### Source Signatures (Table 1 — for MC uncertainty propagation)

| Parameter | Mean | Mean 1σ | Trend (1999→2022) | Trend 1σ |
|-----------|------|---------|-------------------|----------|
| **δ¹³C atmosphere** | −47.3‰ | ±0.01‰ | −0.66‰ | ±0.02‰ |
| **δ¹³C Fossil** | −44.0‰ | ±0.2‰ | +0.60‰ | ±0.14‰ |
| **δ¹³C Microbial** | −62.0‰ | ±1.3‰ | −0.37‰ | ±0.09‰ |
| **δ¹³C Pyrogenic** | −23.1‰ | ±2.5‰ | 0‰ | 0‰ |
| **δ¹³C Net sink KIE** | 1.0082 | 0 | 0 | 0 |
| **δD atmosphere** | −78.4‰ | ±0.5‰ | −6.0‰ | ±0.77‰ |
| **δD Fossil** | −193.4‰ | ±0.6‰ | +1.3‰ | ±0.7‰ |
| **δD Microbial** | −304.7‰ | ±8.2‰ | +1.8‰ | ±0.5‰ |
| **δD Pyrogenic** | −228.0‰ | ±7.6‰ | 0‰ | 0‰ |
| **δD Net sink KIE** | 1.281 | 0 | 0 | 0 |

### Sink Fractionation (SI Table S3)

| Sink | Proportion | ¹³C KIE | D KIE |
|------|-----------|---------|-------|
| OH | 0.835 | 1.0054 (Cantrell) | 1.294 (Saueressig) |
| Cl | 0.035 | 1.066 | 1.520 |
| Stratosphere | 0.06 | 1.003 | 1.179 (Dyonisius 2020) |
| Soil | 0.07 | 1.020 | 1.083 |
| **Net** | 1.00 | 1.0082 | 1.281 |

### Microbial Subcategory Signatures (SI Table S4, year 2022)

| Subcategory | δ¹³C | 1σ | δD | 1σ |
|-------------|------|----|----|-----|
| Wetlands | −62.1‰ | ±0.7‰ | −313.2‰ | ±7.4‰ |
| Ruminants | −65.0‰ | ±1.5‰ | −309.4‰ | ±7.2‰ |
| Waste | −55.1‰ | ±4.4‰ | −279.4‰ | ±27.4‰ |
| Rice | −60.2‰ | ±4.5‰ | −311.8‰ | ±7.3‰ |
| Wild Animals | −65.5‰ | ±3.1‰ | −301.4‰ | ±6.8‰ |
| Termites | −64.8‰ | ±7.6‰ | −323.1‰ | ±8.0‰ |
| **Net Microbial** | **−62.1‰** | **±1.3‰** | **−304.8‰** | **±8.0‰** |

### Key Results

- **Post-2006 CH₄ rise (140 ppb)** overwhelmingly driven by microbial emissions
- δ¹³C mass balance: Mic increase = **73 ± 5 Tg/yr**; FF stable
- δD mass balance: Mic increase = **73 ± 12 Tg/yr**; FF stable
- FF long-term means: **160 ± 29 Tg/yr** (δ¹³C), **133 ± 33 Tg/yr** (δD)
- 2013→2022: FF **declined** by 14±4 (δ¹³C) and 15±11 (δD) Tg/yr
- NH δD–CH₄ is **~12‰ lower** than SH δD–CH₄ (our ±1.5‰ approx was too conservative!)
- Global δD decreased **6.0 ± 0.8‰** between 2005 and 2023
- Global δ¹³C decreased **0.7 ± 0.02‰** between 2005 and 2023
- Geologic emissions: <5 Tg/yr, signature assumed = fossil fuel

### Critical Insight for Our Two-Hemisphere Model

The paper reports: "NH δD–CH₄ is ~12‰ lower than SH δD–CH₄" with seasonal range 8–15‰.

Our v3.0 uses ±1.5‰ hemispheric offset for δD — **this is far too small!** Should be ~6‰ per hemisphere from the global mean (i.e., NH ≈ global − 6‰, SH ≈ global + 6‰). This is a significant correction to implement.

---

## File Locations

```
Ben-BoxModel/
├── datasource.md                          ← Excellent Chinese-language documentation of full data pipeline
├── Ben Riddell-Young ... .pdf             ← The paper (needs MinerU for extraction)
├── Supporting Information.pdf             ← SI (needs MinerU for extraction)
├── Riddell-Young_2025_MassBalancePackage/
│   └── Riddell-Young_2025_MassBalancePackage/
│       ├── d13C_MassBalance_MC.py        ← Main d13C inversion (most relevant to us)
│       ├── dD_MassBalance_MC.py          ← Main dD inversion
│       ├── Mic_d13C.py, Mic_dD.py        ← Microbial signature preprocessing
│       ├── BB_d13C.py                     ← BB signature from C3/C4 maps
│       ├── FF_d13C_GlobMean.py, FF_dD_GlobMean.py  ← FF signature from EDGAR + national DBs
│       ├── data/                          ← All input data including EDGAR8, CT-CH4, station data
│       └── Output/                        ← Pre-computed source signatures + inversion results
└── Riddell-Young_2025_dD_GlobMean/
    └── Riddell-Young_2025_dD_GlobMean/
        ├── dD_globmean.py                ← Station→global δD derivation
        ├── data/                          ← Station-level δD observations (with latitudes!)
        └── output/                        ← Global mean δD + MC iterations
```
