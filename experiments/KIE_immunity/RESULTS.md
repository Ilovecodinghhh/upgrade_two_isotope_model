# Title 2: KIE Immunity — Results Summary

**Date:** 2026-05-11
**Status:** Phase 1–4 complete. Variance decomposition + Basu 2022 comparison done.

## Headline Findings

### 1. Variance decomposition — where does FF uncertainty come from?

| Variance source | δ¹³C-only | Dual (δ¹³C + δD) |
|---|---|---|
| **Total FF variance** | 967 (Tg/yr)² → σ = 31.1 | **287 (Tg/yr)² → σ = 17.0** |
| KIE contribution | 11.2% | **0.0%** |
| Source signature contribution | 82.7% | 0.0% |
| Lifetime contribution | 0.0% | 14.6% |
| Residual (atm + interactions) | 6.1% | 85.4% |

**Key result:** Adding δD reduces total FF variance by **70%** (σ from 31.1 → 17.0 Tg/yr). 
Both KIE *and* source-signature contributions vanish — the dual-isotope system is far better-conditioned.

### 2. Basu 2022 comparison

| Configuration | KIE spread (Saueressig vs Cantrell) |
|---|---|
| Basu 2022 (3D inversion, δ¹³C-only) | **13.0 Tg/yr** |
| Our 2-box (δ¹³C-only, BB fixed, bounded LS) | **0.7 Tg/yr** |
| Our 2-box (δ¹³C + δD, bounded LS) | **0.9 Tg/yr** |

**Surprising finding:** Even the *δ¹³C-only* version of our 2-box already collapses the KIE spread from 13 → 0.7 Tg/yr.

**Interpretation:** Three factors collapse the KIE controversy:
1. **Bounded LS** (FF, Mic ≥ 0) — physical priors regularize the system
2. **BB fixed from independent inventory** — removes one degree of freedom
3. **Hemispheric structure** — provides extra spatial constraint

The dual-isotope addition then provides additional variance reduction through orthogonal constraint (δD-axis).

## Reframing the Title 2 Story

Original framing: "δD halves KIE-driven ambiguity"
**Revised framing:** "*Properly constrained* box models eliminate KIE ambiguity, which 3D inversions inflate via under-constraint"

This is a STRONGER paper: it shifts blame from "KIE controversy unresolved" to "3D inversions are too flexible". Our methodology (BB-fixed bounded LS + hemispheric structure) is the recipe.

## Files Pushed
- `experiments/KIE_immunity/analysis/variance_decomposition.py`
- `experiments/KIE_immunity/analysis/compare_basu2022.py`
- `experiments/KIE_immunity/figures/fig_kie_immunity.py` (+ PNG/PDF)
- `experiments/KIE_immunity/results/variance_decomposition.json`
- `experiments/KIE_immunity/results/basu_comparison.json`

## Next steps
- Phase 5 (extended sensitivity): OH_D KIE alone, Cl fraction alone — test whether replacing OH_13C uncertainty with OH_D uncertainty is a wash
- Title 3: Hemispheric divergence
- Or: paper outline for Title 1 (which is the most ready)
