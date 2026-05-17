# 2026-05-06 — Analysis of Ben's δD-CH₄ Global Mean Database

## Task
Investigate how Ben (Riddell-Young 2025) constructed the hydrogen isotope (δD-CH₄) atmospheric database, identify sources of uncertainty, and determine what drives the 2008 anomaly in our 3×3 box model.

## Key Findings: How Ben Built the δD Database

### Architecture: `Riddell-Young_2025_dD_GlobMean` Package

**Two-stage process:**
1. **`MBL_calc_Unc.py`** — Per-station curve fitting with MC uncertainty
   - Loads raw station δD data (33 station files from 4 labs: INSTAAR, MPI, IMAU, NIPR)
   - Applies NOAA `ccg_filter` smoothing (shortterm=150 days, longterm=667 days, 4 harmonics, 3 poly terms)
   - Computes residual standard deviation (rsd) for each station
   - Runs 1000 MC iterations: adds random monthly noise (drawn from N(0, rsd)) to data, re-fits curve → produces 1001-column matrix (date + 1000 smoothed curves)

2. **`dD_globmean.py`** — Global mean construction with MC propagation
   - Divides stations into 4 latitudinal bands: PN (30-90°N), TN (0-30°N), TS (0-30°S), PS (30-90°S)
   - For each of 1000 MC iterations:
     - **Network uncertainty**: randomly drops 2 stations per iteration
     - **Scale uncertainty**: adds N(0, σ) noise for inter-lab offsets (IMAU: σ=2.2‰, NIPR: σ=1.6‰, INSTAAR: σ=1.6‰)
     - **Scale correction**: Applies offsets to harmonize to MPI scale (INSTAAR -1.8‰, IMAU -0.5‰, NIPR 0‰ in Dasgupta calibration)
     - Averages stations in each band per week (52 weeks/yr)
     - Fills spatial gaps using mean weekly difference between bands
     - Calculates global mean = (PN + TN + TS + PS) / 4 (equal area weighting)
   - Applies "anomalous jump" filter: if point differs from both neighbors by >0.3‰, replace with mean of neighbors

### Uncertainty Sources (ranked by impact)

| Source | Magnitude | Notes |
|--------|-----------|-------|
| **Inter-lab scale offsets** | 1.6-2.2‰ (1σ) | Biggest systematic uncertainty. MPI vs INSTAAR vs IMAU vs NIPR |
| **Network sparsity (early years)** | Large | Only 1 TN station (kum/eom), 2 TS stations in 2005-2009 |
| **Station dropout MC** | Moderate | Removing 2/12 stations in 2007 = 17% of network |
| **Atmospheric/analytical noise** | ~1-3‰ per station | Captured by rsd from ccg_filter |
| **Gap-filling** | Moderate | TS band is entirely gap-filled 2010-2012 (0 stations!) |
| **Spatial weighting** | Systematic | Equal 4-band weighting regardless of # stations per band |

### Critical Station Coverage Around 2007-2009

**MBL stations only (used in published global mean):**
- **PN (30-90°N)**: 5-6 stations (alt, brw, cba, mhd, nyaNIPR, zep) — WELL SAMPLED
- **TN (0-30°N)**: 1 station (kum only!) — VERY SPARSE
- **TS (0-30°S)**: 2 stations (asc, smo) — SPARSE
- **PS (30-90°S)**: 4 stations (cgo, gvnIMAU, spo, syoNIPR) — OK

**Key problem**: The tropical bands (TN, TS) are dramatically under-sampled in 2007-2009. A single station (Kumukahi, Hawaii) represents all of 0-30°N. The tropics are where most methane sources are! This means the early δD record is essentially:
- Well-constrained at high latitudes
- Poorly constrained in the tropics where source mixing matters most

### NH-SH Difference (from hemispheric output)
- 2005-2009: NH-SH ≈ -10 to -11‰ (NH more depleted)
- 2010-2024: NH-SH ≈ -11 to -14‰ (grows over time)
- **2007-2008 specifically: NH-SH shifts from -9.88 to -10.08‰** — relatively small change

### Why 2008 Breaks Our 3×3 Model

The derived source δD at 2008 is anomalously HEAVY (-276.7‰ vs ~-282‰ typical), while δ¹³C goes LIGHT (-54.74‰). This contradiction is fatal for the 3×3.

Root causes in the δD data:
1. **Edge-of-record effects**: δD record starts 2005, so 2008 is only year 3. The ccg_filter curve fitting has edge instabilities (polynomial + harmonics extrapolate poorly at boundaries)
2. **The 2008→2009 "non-event"**: atmospheric δD barely changes (-0.05‰) despite 5.6 ppb/yr growth. This means the derived source δD must be isotopically heavy (~-277‰) — but there's no physical source mix that's that heavy while simultaneously being -54.7‰ in δ¹³C
3. **INSTAAR stations dominate 2005-2009** (then end ~2010). The scale correction (-1.8‰) applied uniformly doesn't account for potential time-dependent calibration drift
4. **Sparse tropical coverage**: only kum (19.5°N) represents all tropics → the global mean is essentially "high-latitude average minus gap-fill correction"

### What We Can Do to Improve

**Short-term fixes for our model:**
1. **Start δD constraint from 2010, not 2006** — the pre-2010 record has too few stations and edge effects. Use δ¹³C-only for 2006-2009.
2. **Apply 3-year running mean to δD** before inversion — smooths out the 2008-2009 noise
3. **Inflate δD uncertainty for 2006-2009** — reflect the actual network sparsity (12-13 stations vs 6-8 later)
4. **Use the hemispheric δD output directly** from Ben's package (NH/SH means are available!) instead of our ±6‰ approximation

**Medium-term improvements:**
5. **Non-MBL stations**: Include stations Ben excluded (ato, bal, bsc, lef, mlo, etc.) — adds ~8 more sites but introduces continental influence
6. **Weight by inverse-uncertainty** instead of equal band averaging
7. **Time-dependent inter-lab correction**: INSTAAR→MPI offset may drift over time
8. **Use ccg_filter trend (not smoothed)** for the mass balance — removes seasonal aliasing

**Structural insight:**
The fundamental limitation is that **δD-CH₄ atmospheric monitoring started ~15 years later than δ¹³C** (2005 vs 1990s). The network is still small and relies heavily on a few labs. Any combined δ¹³C+δD inversion inherits this weakness. Ben's choice to run them separately and compare trends (rather than combining in one inversion) is methodologically safer.

## File Locations
- Ben's dD package: `Ben-BoxModel/Riddell-Young_2025_dD_GlobMean/Riddell-Young_2025_dD_GlobMean/`
- Main script: `dD_globmean.py`
- Upstream MC: `MBL_calc_Unc.py`
- Station data: `data/*.txt`
- Station metadata: `data/siteinfo_all_ch4h2.txt`
- Hemispheric output: `output/HemMean_dD_dei_UmezawaCal_noBUDS.csv`
- Global mean output: `output/GlobMean_dD_iterations_UmezawaCal_noBUDS.xlsx`
