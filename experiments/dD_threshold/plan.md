# plan.md — Future Work for `experiments/dD_threshold`

## Current State (2026-05-12, v4)

Fully hemispheric model with both δ¹³C and δD source signatures.

Key results:
- **Baseline improvement:** 45.1% (CI = 57.6 vs 105.1 Tg/yr)
- **Threshold:** σ(Mic δD) ≈ 35‰ (multiplier ~4×)
- **NH drives the improvement** (82% NH, 7% SH from Phase 6 deep dive)
- **Robust** across all KIE (×3) and lifetime (×3) configurations

### Hemispheric source signature gaps

| Isotope | FF gap | BB gap | Mic gap |
|---------|--------|--------|---------|
| δD | −7‰ | −24‰ | −13‰ |
| δ¹³C | +4.5‰ | −2.4‰ | ~0‰ |

δD has 5–10× larger gaps → explains why δD adds hemispheric info that δ¹³C cannot.

---

## Completed Steps

### Data generation
- [x] Hemispheric δD source signatures (from `rel/build_hemispheric_dD_sources.py`)
- [x] Hemispheric δ¹³C source signatures (from `rel/build_hemispheric_d13C_sources.py`)
  - FF: country-level ONG + coal δ¹³C × EDGAR emissions, assigned to NH/SH
  - BB: C3/C4 vegetation maps × CTCH4 pyrogenic flux, hemisphere-split
  - Mic: subcategory mass balance (wetlands, ruminants, rice, etc.) with Oh 2022 + Still 2003

### Loading & integration
- [x] `common.py`: loads all 12 hemispheric MC CSV files (6 δD + 6 δ¹³C)
- [x] `common.py`: `sample_source_signatures_hemi()` returns 18 keys (6 global + 6 dD hemi + 6 d13C hemi)
- [x] All analysis scripts: NH A-matrix uses NH δ¹³C + NH δD; SH uses SH versions
- [x] δ¹³C-only mode also uses hemispheric δ¹³C for consistent comparison

### Analysis phases
- [x] Phase 1 baseline: CI = 55.6 Tg/yr, +45.2%
- [x] Phase 2 DFS: ΔDFS = +1.39 (unchanged by source sig changes)
- [x] Phase 3 threshold: crossover between 3× (+26.1%) and 5× (−37.2%)
- [x] Phase 3b Thanwerdas: +45.1% ours vs −138.5% theirs
- [x] Phase 5 sensitivity: threshold at 5× for all 6 configs
- [x] Phase 6 deep dive: exact crossover, hemispheric breakdown, bootstrap, year-range

---

## Version History

| Version | δD atm | δD src | δ¹³C src | Dual CI | Improvement | Threshold |
|---------|--------|--------|----------|---------|-------------|-----------|
| v1 | Umezawa, ±6‰ | Global | Global | 46.6 | +52% | ~25‰ |
| v2 | Dasgupta, real hemi | Global | Global | 37.8 | +60.8% | ~41‰ |
| v3 | Dasgupta, real hemi | Hemi | Global | 43.5 | +57.0% | ~41‰ |
| **v4** | Dasgupta, real hemi | **Hemi** | **Hemi** | **57.6** | **+45.1%** | **~35‰** |

---

## Open Questions / Remaining Work

### High Priority (for publication)

- [x] **Re-run Phase 6 deep dive with v4 data**: Exact crossover at 3.82× (σ=31.5‰), 10% threshold at 3.53× (σ=29.1‰). Bootstrap: 44.7±1.6%, P(>0)=100%, P(>30%)=100%
- [x] **Update figures**: `fig_comprehensive.py` re-run with Phase 6 fine-grid data, v4 version comparison, exact crossover annotation
- [x] **Paper Table 1**: Updated in RESULT.md with all v4 Phase 6 numbers

### Model extensions

- [ ] **3-box model**: ThreeBox data exists (NHext/Trop/SHext). Could test finer spatial resolution
- [ ] **Semi-hemispheric (4-box)**: SemiHemMean δD data available but has NaN gaps
- [ ] **Microbial d13C spatial variation**: Currently Mic δ¹³C shows ~0 hemispheric gap. If wetland d13C spatial map (isotem) becomes available, could improve this

### Sensitivity extensions

- [ ] **Per-source threshold**: Inflate Mic δD alone vs FF+BB alone to isolate which source's uncertainty matters
- [ ] **Time-varying threshold**: Does crossover shift over the study period?
- [ ] **W matrix optimization**: Current weights are fixed; could test sensitivity or optimize

### Data improvements

- [ ] **2020–2023 atmospheric δD gap**: Currently gap-filled; replace when station data arrives
- [ ] **Luo 2023 C4 map**: Currently using Still 2003 (static); Luo is time-varying and higher-res
- [ ] **Prior emission subcategory files**: Would improve Mic δ¹³C hemispheric calculation

---

## How to Re-run

```bash
cd upgrade_two_isotope_model

# Build hemispheric source signatures (if data changes):
python3 rel/build_hemispheric_dD_sources.py    # ~5 min
python3 rel/build_hemispheric_d13C_sources.py  # ~15 min

# Core phases:
python3 experiments/dD_threshold/analysis/phase1_baseline.py     # ~3 min
python3 experiments/dD_threshold/analysis/phase2_dfs.py          # ~10 sec
python3 experiments/dD_threshold/analysis/phase3_threshold.py    # ~15 min
python3 experiments/dD_threshold/analysis/phase3b_thanwerdas.py  # ~3 min
python3 experiments/dD_threshold/analysis/phase5_sensitivity.py  # ~20 min
python3 experiments/dD_threshold/analysis/phase6_deep_dive.py    # ~30 min

# Figures:
python3 experiments/dD_threshold/figures/fig_comprehensive.py    # ~5 sec
```
