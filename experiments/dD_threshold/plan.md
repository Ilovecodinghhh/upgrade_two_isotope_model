# plan.md — Future Work for `experiments/dD_threshold`

## Current State (2026-05-12)

All 5 phases re-run with three data upgrades:
1. **Dasgupta calibration** (replacing Umezawa) for atmospheric δD
2. **Real hemispheric atmospheric δD** MC iterations (replacing ±6‰ offset)
3. **Hemispheric δD source signatures** (FF, Mic, BB for NH/SH separately)

Key result: **δD threshold at σ(Mic δD) ≈ 35–40‰**, improvement of ~57% at baseline.

---

## Completed Steps

- [x] `common.py`: Dasgupta calibration, hemispheric atm δD loading, NaN gap-fill (2020–2023)
- [x] `common.py`: `sample_atm_dD_hemi()` for real NH/SH atmospheric δD
- [x] `common.py`: `sample_source_signatures_hemi()` for hemisphere-specific δD source signatures
- [x] `common.py`: Load `{FF,Mic,BB}_dD_{NH,SH}_MC.csv` hemispheric source signature files
- [x] Phase 1 (baseline): uses hemi atm + hemi src sigs → CI = 42.3 Tg/yr, +56.2%
- [x] Phase 2 (DFS): unchanged (depends on source spread, not observations)
- [x] Phase 3 (threshold sweep): threshold at 5× (~41‰), crossover between 3× and 5×
- [x] Phase 3b (Thanwerdas): confirmed uncertainty spec kills δD (+57% ours vs −150% Thanwerdas)
- [x] Phase 5 (sensitivity): threshold robust across all KIE and lifetime variants
- [x] RESULT.md updated with all new numbers

---

## Open Questions / Next Steps

### Data & Model

- [ ] **Semi-hemispheric data**: `ThreeBox_dD_sources_summary.csv` exists but is unused — could enable a 3-box model
- [ ] **δ¹³C hemispheric source signatures**: currently only δD source signatures are hemisphere-specific; δ¹³C uses global. If hemispheric δ¹³C source data becomes available, update `sample_source_signatures_hemi()` to include `ff_d13C_NH/SH` etc.
- [ ] **2020–2023 gap**: hemispheric atmospheric δD filled with global + offset for those years; when station data becomes available, replace the fill
- [ ] **Year range alignment**: δD covers 2005–2023, δ¹³C/CH₄ covers 1999–2022; front-padding δD with 2005 values for 1999–2004 is crude — consider whether pre-2005 years should be excluded from CI calculations

### Analysis Improvements

- [ ] **Interpolate exact crossover**: phase 3 currently tests discrete multipliers (0.5, 1, 2, 3, 5, 8, 12, 16×). Add finer grid around 3–5× (e.g., 3.5, 4.0, 4.5) to pinpoint exact threshold
- [ ] **Posterior diagnostics**: add per-iteration residual checks to flag when the constrained LSQ hits bounds (currently silently clamps)
- [ ] **Bootstrap CI on CI**: the 90% CI width is itself estimated from 500–1000 MC draws. Could bootstrap to get uncertainty on the CI width itself
- [ ] **Hemispheric breakdown in Phase 3 output**: currently only reports global FF CI; could separately report NH vs SH CI to see which hemisphere drives the threshold

### Figures

- [ ] **Update `fig_threshold.py`** and `fig_comprehensive.py` with new numbers (currently generated from old results)
- [ ] **Add NH-SH source signature comparison panel** showing the hemispheric δD gaps
- [ ] **Add data provenance panel** showing Umezawa vs Dasgupta calibration difference

### Publication

- [ ] **Table 1**: Update all numbers in the manuscript draft
- [ ] **Discussion**: emphasize that hemispheric source-signature heterogeneity (especially BB: 24‰ gap) is a key information source for the two-box model
- [ ] **Robustness narrative**: threshold stable across all 3 data versions (~25‰ → ~41‰ → ~35–40‰); if anything, better data makes δD *more* useful

---

## How to Re-run

```bash
cd upgrade_two_isotope_model

# All phases (order matters for phase 2 which is fast):
python3 experiments/dD_threshold/analysis/phase1_baseline.py   # ~3 min
python3 experiments/dD_threshold/analysis/phase2_dfs.py        # ~10 sec
python3 experiments/dD_threshold/analysis/phase3_threshold.py  # ~15 min
python3 experiments/dD_threshold/analysis/phase3b_thanwerdas.py # ~3 min
python3 experiments/dD_threshold/analysis/phase5_sensitivity.py # ~20 min
```

All scripts auto-save to `experiments/dD_threshold/results/`.

---

## Key Files Modified

| File | What changed |
|------|-------------|
| `common.py` | Dasgupta cal, hemi atm δD, hemi src sigs, `sample_atm_dD_hemi()`, `sample_source_signatures_hemi()` |
| `phase1_baseline.py` | Uses `sample_source_signatures_hemi()`, separate NH/SH A-matrices |
| `phase3_threshold.py` | Same + `inflate_dD_uncertainty()` inflates hemi keys |
| `phase3b_thanwerdas.py` | Same + `apply_thanwerdas_uncertainties()` handles hemi keys |
| `phase5_sensitivity.py` | Same + inflation handles hemi keys |
| `RESULT.md` | Updated with all v3 numbers |
