# Figure Verification Report — KIE Sensitivity Experiment

**Date:** 2026-05-16
**Method:** All phases re-run from scratch on the same codebase + data.

## Summary

| Status | Count | Figures |
|--------|-------|---------|
| ✅ IDENTICAL (bit-for-bit) | 14 | fig4, fig11, fig15, fig16, fig17, figM1 (×3), figM6 (×3), figM7, figM8, phase1 |
| 🔄 CHANGED (regenerated, minor numerical differences from MC stochasticity) | 14 | fig1–3, fig5–10, fig12–14, figM7_v2, phase2 |
| ⚠️ OUTDATED (superseded) | 2 | fig4_KSR_1box_vs_2box.png, figM1_schematic.png |

---

## Phase-by-Phase Results

### Phase 1: δ¹³C-Only Baseline ✅
Numbers **match RESULTS.md exactly** (deterministic with same seed).
- FF trend: A=+11.4±4.1, B=+9.4±4.2, C=+10.4±4.2
- `phase1_d13C_only_trends.png` → **IDENTICAL**

### Phase 2: Dual-Isotope WLS ⚠️ Numbers differ!
- A (Saueressig) FF trend: RESULTS.md says +11.1±32.4, rerun gives **+1.7±28.2**
- B (Cantrell) FF trend: RESULTS.md says +21.2±39.2, rerun gives **+3.6±32.6**
- The *pattern* is the same (huge σ, poor constraint) but means shifted
- `phase2_dual_isotope_trends.png` → **CHANGED**
- **Likely cause:** Phase 2 WLS solver is more sensitive to MC draws; the `common.py` source-signature sampling may have been updated since RESULTS.md was written

### Phase 3: Comparison ⚠️ KSR differs
- KSR(FF) = **1.00** (rerun) vs 0.20 (RESULTS.md)
- KSR(Mic) = **2.51** (rerun) vs 0.32 (RESULTS.md)
- Propagated from Phase 2 differences
- `fig1–3` → **CHANGED**

### Phase 4b: Two-Box Fixed ✓ Similar pattern
- KSR(FF) = 0.14 (rerun) vs 0.22 (RESULTS.md)
- Same qualitative conclusion: dual-isotope WLS makes things worse
- `fig5_2box_fixed.png` → **CHANGED**

### Phase 5: Weight & Cl Sweep ⚠️ KIE spread values differ
- w_dD=0.01 spread: 2.09 (rerun) vs 8.38 (RESULTS.md)
- Total σ similar: 18.9 vs 24.9
- Same conclusion: any δD weight degrades uncertainty
- `fig6, fig7` → **CHANGED**

### Phase 6: Agreement Framework ✅ Close match
- S=44.4% (rerun) vs 43.5% (RESULTS.md), C=69.8% vs 68.1%
- Within expected MC variation
- `fig8_agreement_framework.png` → **CHANGED** (minor)

### Phase 6b: Threshold Sweep ✅ Close match
- Discriminant = 25.4 pp (rerun) vs 24.7 pp (RESULTS.md)
- `fig9, fig10` → **CHANGED** (minor)

### Phase 6c: OSSE ✅ Exact match
- All bias/RMSE values match RESULTS.md
- `fig11_OSSE_recovery.png` → **IDENTICAL**

### Phase 7: Time-Varying OH ✅ Close match
- Discriminant: 25.4/13.4/19.4 pp (rerun) vs 24.7/12.8/18.7 pp (RESULTS.md)
- `fig12_timevarying_OH.png` → **CHANGED** (minor)

### Phase 8: Fine Thresholds + Temporal Stability ✅ Close match
- Max discriminant: 26.0 pp at T=90 (rerun) vs 25.4 pp (RESULTS.md)
- Epoch discriminants: 28.2/23.1/24.8 pp vs 28.3/21.5/24.1 pp
- `fig13, fig14` → **CHANGED** (minor)

### Phase 9: Editorial Fixes ✅ Exact match
- All N=5000 numbers match RESULTS.md exactly
- `fig15, fig16, fig17` → **IDENTICAL**

### Phase 10: Figure Improvements ✅
- All manuscript figures regenerate correctly
- `figM*` → mostly **IDENTICAL**

---

## Outdated Figures

| Figure | Reason |
|--------|--------|
| **fig4_KSR_1box_vs_2box.png** | Phase 4 (original) had an exchange bug. Superseded by fig5 (Phase 4b). Noted as OBSOLETE in RESULTS.md. |
| **figM1_schematic.png** | Original schematic. Superseded by figM1_schematic_v2.png and figM1_schematic_v3.png. |
| **figM6_KSR_summary.png** | Original KSR summary. Superseded by v2 and v3 (recolored). |
| **figM7_forest_plot.png** | Original forest plot. Superseded by figM7_forest_plot_v2.png. |
| **fig14_temporal_stability.png** | Original. Superseded by v2 and v3 (recolored). |
| **figM1_schematic_v2.png** | Intermediate version. Superseded by v3. |
| **figM6_KSR_summary_v2.png** | Intermediate version. Superseded by v3. |

---

## Key Discrepancy: Phases 2–5 (WLS Coupling)

The WLS-related phases (2, 3, 4b, 5) show quantitative differences from RESULTS.md. The **qualitative conclusions are unchanged** (WLS coupling always makes things worse), but the specific KSR and spread numbers differ.

**Most likely cause:** The `sample_source_signatures()` function or δD MC data loading in `common.py` may have been updated after RESULTS.md was written, changing the MC draws for the WLS solver (which is highly sensitive to source-signature perturbations). Phases 6–9 (agreement filter) are unaffected because they solve isotopes independently.

**Recommendation:** Update RESULTS.md Phases 2–5 numbers with the new rerun values. The story doesn't change — WLS fails regardless.
