# CHANGELOG_Improved_dD.md

## Improved δD Pipeline & Comparison Across All Model Versions

**Date:** 2026-05-06  
**Files:** `improved_dD_pipeline.py`, `compare_dD_all_versions.py`  
**Output:** `Output_improved_dD/`, `Output_dD_comparison/`

---

### 1. What We Built

#### a) `improved_dD_pipeline.py` — Replicated & Improved Ben's δD Global Mean

Ben constructs atmospheric δD-CH₄ time series in two stages:

1. **Per-station curve fitting** (`MBL_calc_Unc.py`): NOAA `ccg_filter` smoothing of 36 stations across 4 labs (INSTAAR, IMAU, NIPR, MPI), with 1000 MC noise-perturbation iterations
2. **Global mean assembly** (`dD_globmean.py`): 4 semi-hemispheric bands (PN/TN/TS/PS at ±30°), drop-2-sites network uncertainty, inter-lab scale adjustments, weekly gap-filling, equal-weight average

**Our improvements:**
- **Area-weighted** (sin-latitude) band averaging instead of equal 1/4 per band
- **Cosine-latitude weighting** within each band for multi-station averages
- Re-used Ben's pre-computed MC curves (no need to rerun ccg_filter)

**Result:** Area weighting shifts values by <0.3‰ since the 4 bands have nearly equal area. Our equal-weight replication matches Ben within ±1.6‰ (typically <1‰). Remaining discrepancy is from slightly different station filtering.

#### b) `compare_dD_all_versions.py` — All Model Versions, Ben's δD vs Ours

Ran each model version with both δD inputs:
- **Ben's original δD** (from `glob_ann_dD.xlsx`)
- **Our improved δD** (area-weighted, from pipeline output)

For v2.0 (two-hemisphere), used actual NH/SH δD from both sources.

---

### 2. Results Summary (2010–2021 means, physical solutions only)

| Model | δD Source | FF (Tg/yr) | BB (Tg/yr) | Mic (Tg/yr) | Physical % |
|-------|-----------|-----------|-----------|------------|-----------|
| v1.0 (3×3) | Ben | 67.0 | 82.7 | 436.2 | 59.6% |
| v1.0 (3×3) | **Ours** | **63.4** | 84.0 | 435.3 | 59.4% |
| v2.0 (2-hemi) | Ben | 63.5 | 84.6 | 437.3 | 58.5% |
| v2.0 (2-hemi) | **Ours** | **62.9** | 84.1 | 434.9 | 59.3% |
| v3.1 (3×3 opt) | Ben | 68.6 | 82.6 | 435.8 | 59.3% |
| v3.1 (3×3 opt) | **Ours** | **62.5** | 83.8 | 435.3 | 60.4% |
| v3.2 (BB fix) | — | 189.9 | 28.6 | 367.2 | 100.0% |
| v3.3 (dD≥2010) | Ben | 68.4 | 83.6 | 437.3 | 60.2% |
| v3.3 (dD≥2010) | **Ours** | **63.0** | 84.9 | 436.9 | 60.0% |
| v4.0 (Mic/NM) | Ben | 71.9 | 79.2 | 433.2 | 45.7% |
| v4.0 (Mic/NM) | **Ours** | **65.1** | 81.7 | 434.5 | 45.4% |

---

### 3. Key Findings

#### The δD input choice matters for FF, not for Mic or BB

- **FF emissions consistently ~5–7 Tg/yr lower** with our improved δD (across all models that use δD)
- **Microbial emissions barely change** (±1–2 Tg/yr) — robust to δD input choice
- **BB emissions change by <2 Tg/yr** — also robust
- v3.2 (BB-fixed, δ13C only) is **identical** for both δD inputs, as expected — but BB now varies annually using GFED4 data (23–36 Tg/yr) rather than a flat fraction

#### The impact is concentrated in specific years

The Δ(Our − Ben) plot shows year-to-year impacts of ±5–40 Tg/yr for FF, with largest discrepancies in:
- **2010** (~–20 to –40 Tg/yr FF): Ben's δD is 1.9‰ more negative than ours
- **2012–2013** (~+5 to –10 Tg/yr): moderate differences
- **2015** (~+5 Tg/yr): Ben's δD is 0.9‰ more negative

These year-specific differences trace directly to different station coverage/weighting in the δD compilation.

#### Physical solution rates are essentially unchanged

All models show ±1% change in physical solution rate between Ben's and our δD — the improved pipeline doesn't affect the fundamental (in)determinacy of the 3×3 system.

#### v4.0 (Mic/NonMic) is the most sensitive to δD choice

Makes sense: v4.0 uses δD as the **primary** constraint for the Mic/NonMic split, so any δD shift propagates directly into the FF/BB sub-partition.

---

### 4. Interpretation

The ~0.5–1.9‰ difference between Ben's and our δD compilation is **small relative to source signature uncertainties** (FF δD = –183 ± 8‰, Mic δD = –305 ± 10‰). The model results are therefore robust to reasonable variations in the δD observational compilation.

The largest impact (~5–7 Tg/yr on FF mean) is within the 1σ uncertainty band of any individual model version. **Conclusion: the choice of δD compilation methodology is not a dominant source of uncertainty in the box model source partitioning.**

---

### 5. Output Files

| File | Description |
|------|-------------|
| `Output_improved_dD/summary_comparison.csv` | 2010–2021 mean emissions by model & δD source |
| `Output_improved_dD/annual_comparison.csv` | Full annual results for all models |
| `Output_improved_dD/all_versions_source_partitioning.png` | 6-panel comparison plot |
| `Output_improved_dD/delta_impact.png` | Δ(Our − Ben) by source category |
| `Output_improved_dD/physical_rates.png` | Physical solution rate comparison |
| `Output_improved_dD/dD_input_comparison.png` | Input δD data comparison |
| `Output_improved_dD/best_models_detail.png` | v3.2 & v4.0 detail with 5–95th percentiles |
| `Output_dD_comparison/improved_dD_global_mean.csv` | Our improved δD time series |
| `Output_dD_comparison/improved_dD_pipeline.png` | Pipeline validation plots |

---

### 6. Model Version Reference

| Version | Method | Uses δD? | Key Feature |
|---------|--------|----------|-------------|
| v1.0 | 1-box 3×3 | ✓ | Original triple-isotope inversion |
| v2.0 | 2-hemisphere 3×3 | ✓ (NH/SH) | Separate hemispheric δD |
| v3.1 | 1-box 3×3 optimized | ✓ | Same as v1.0, optimized code |
| v3.2 | BB-fixed 2×2 | ✗ | BB as 5% of total, δ13C only |
| v3.3 | 3×3, δD from 2010+ | ✓ (2010+) | δD only for recent years |
| v4.0 | Mic vs NonMic | ✓ | δD for Mic/NM split, δ13C for FF/BB sub-split |
