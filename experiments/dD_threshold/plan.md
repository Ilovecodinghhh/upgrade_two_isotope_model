# plan.md — Future Work for `experiments/dD_threshold`

## Current State (2026-05-12)

All 6 phases complete with three data upgrades + deep-dive analysis.

Key results:
- **Exact crossover:** σ(Mic δD) = **32.0‰** (improvement = 0%)
- **10% criterion:** σ(Mic δD) = **29.4‰**
- **Baseline improvement:** 56.6% (CI = 44.0 vs 101.3 Tg/yr)
- **NH drives the threshold:** 82% improvement in NH, only 7% in SH at baseline
- **Bootstrap-confirmed:** 55 ± 2%, P(>30%) = 100%

---

## Completed Steps

### Data & Loading
- [x] `common.py`: Dasgupta calibration, hemispheric atm δD loading, NaN gap-fill (2020–2023)
- [x] `common.py`: `sample_atm_dD_hemi()` for real NH/SH atmospheric δD
- [x] `common.py`: `sample_source_signatures_hemi()` for hemisphere-specific δD source signatures
- [x] `common.py`: Load `{FF,Mic,BB}_dD_{NH,SH}_MC.csv` hemispheric source signature files

### Analysis Phases
- [x] Phase 1 (baseline): CI = 42.3 Tg/yr, +56.2%
- [x] Phase 2 (DFS): 2-box ΔDFS = +1.39
- [x] Phase 3 (threshold sweep): coarse grid, threshold at 5× (~41‰)
- [x] Phase 3b (Thanwerdas): confirmed uncertainty spec kills δD
- [x] Phase 5 (sensitivity): robust across all KIE + lifetime variants
- [x] **Phase 6A (fine grid)**: exact crossover at 3.88× (σ = 32.0‰); 10% at 3.57× (σ = 29.4‰)
- [x] **Phase 6B (hemispheric breakdown)**: NH = 82% improvement, SH = 7%; threshold driven by NH
- [x] **Phase 6C (bootstrap)**: 55.1 ± 1.7%, P(>0) = 100%, P(>30%) = 100%
- [x] **Phase 6D (year range)**: pre-2005 padding has <2% effect
- [x] **Phase 6E (bound diagnostics)**: bounds active ~100% (expected for constrained LSQ)

### Figures & Documentation
- [x] Updated `fig_comprehensive.py` → 6-panel + 4-panel figures
- [x] Updated `RESULT.md` with all Phase 6 numbers
- [x] Added hemispheric source signature gap panel (Panel E)
- [x] Added data version comparison panel (Panel F)

---

## Open Questions / Remaining Work

### Data & Model Extensions

- [ ] **3-box model**: `ThreeBox_atm_dD_annual.csv` and `ThreeBox_dD_sources_summary.csv` exist with NHext/Trop/SHext breakdown. Could test whether finer spatial resolution further improves constraint or just adds noise
- [ ] **Semi-hemispheric data**: `SemiHemMean_dD_dei_DasguptaCal_noBUDS.csv` (PN/TN/TS/PS) could enable 4-box, but sub-annual resolution and many NaN gaps make this challenging
- [ ] **δ¹³C hemispheric source signatures**: currently δ¹³C uses global sigs. If hemispheric δ¹³C data becomes available, add to `sample_source_signatures_hemi()`
- [ ] **2020–2023 gap**: hemispheric atm δD filled with global + offset; replace when station data arrives

### Analysis Refinements

- [ ] **Per-source threshold**: current analysis inflates all 3 source δD uncertainties simultaneously. Could inflate Mic alone (dominant) vs FF+BB alone to test which source's uncertainty matters most
- [ ] **Time-varying threshold**: does the crossover shift over the study period? (e.g., is δD more useful post-2010 when emissions changed faster?)
- [ ] **Optimal weighting (W matrix)**: currently W_NH = diag(100,1,0.5) and W_SH = diag(200,1,0.5). Could optimize these weights or test sensitivity
- [ ] **Posterior residual analysis**: quantify how much residual scatter exists after LSQ inversion; flag systematic biases

### For Publication

- [ ] **Table 1**: Fine-grid threshold numbers (Phase 6A) are now the definitive values — update manuscript
- [ ] **Key narrative elements**:
  - The NH drives the improvement (82% vs 7%) — this is new and publishable
  - The exact crossover is σ = 32‰, not the coarse ~41‰
  - Bootstrap gives ± 2% uncertainty on the improvement itself
- [ ] **Supplementary figure**: year-by-year CI breakdown for NH vs SH
- [ ] **Robustness table**: summarize all 6 sensitivity dimensions (KIE ×3, lifetime ×3, year range ×3, bootstrap) in one table

---

## How to Re-run

```bash
cd upgrade_two_isotope_model

# Core phases (~25 min total):
python3 experiments/dD_threshold/analysis/phase1_baseline.py     # ~3 min
python3 experiments/dD_threshold/analysis/phase2_dfs.py          # ~10 sec
python3 experiments/dD_threshold/analysis/phase3_threshold.py    # ~15 min
python3 experiments/dD_threshold/analysis/phase3b_thanwerdas.py  # ~3 min
python3 experiments/dD_threshold/analysis/phase5_sensitivity.py  # ~20 min

# Deep dive (~30 min):
python3 experiments/dD_threshold/analysis/phase6_deep_dive.py    # ~30 min

# Figures (fast, reads saved JSON):
python3 experiments/dD_threshold/figures/fig_comprehensive.py    # ~5 sec
```

---

## Key Files

| File | Purpose |
|------|---------|
| `common.py` | Shared loading, sampling, KIE, lifetime functions |
| `analysis/phase1_baseline.py` | 1-box/2-box baseline comparison |
| `analysis/phase2_dfs.py` | DFS information content |
| `analysis/phase3_threshold.py` | Coarse threshold sweep (8 multipliers) |
| `analysis/phase3b_thanwerdas.py` | Thanwerdas uncertainty replication |
| `analysis/phase5_sensitivity.py` | KIE + lifetime robustness |
| `analysis/phase6_deep_dive.py` | Fine grid + hemi breakdown + bootstrap + year-range |
| `figures/fig_comprehensive.py` | 4-panel and 6-panel publication figures |
| `RESULT.md` | Complete results summary |
| `plan.md` | This file |
