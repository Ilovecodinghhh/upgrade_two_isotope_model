# Technical Review — Hemispheric Divergence Experiment

**Reviewer:** Independent Technical Review  
**Date:** 2026-05-14  
**Scope:** `experiments/Hemispheric_Divergence/` — code, methodology, results, and publication readiness  
**Target journals considered:** Nature Communications / ACP

---

## 1. Executive Summary

This experiment develops a two-box (NH/SH) dual-isotope (δ¹³C + δD) mass-balance model to partition global methane emissions into fossil fuel (FF), microbial (Mic), and biomass burning (BB) components per hemisphere. The central hypothesis — that "spatial aliasing" explains the long-standing disagreement between 3D inversions (FF increasing) and 1-box models (FF stable/declining) — is creative and potentially important. However, the experiment in its current state has several **fundamental methodological issues** that undermine the quantitative results and must be resolved before publication. The qualitative framework and scientific question are sound.

**Overall assessment: Major revisions required.**

---

## 2. Strengths

### 2.1 Scientific Framing
The reconciliation hypothesis (Section in SUMMARY.md) is the strongest asset. The idea that 1-box models and 3D inversions may both be "correct" because they are measuring different spatial aggregates is genuinely insightful and addresses the most contentious debate in contemporary methane science (Basu 2022 vs. Riddell-Young 2025 vs. He 2026).

### 2.2 Model Evolution (v1 → v2 → v3)
The progressive improvement from fraction-space solvers (v1/v2, condition number ~170,000) to the delta-space solver (v3, condition number ~27 after scaling) demonstrates good scientific self-correction. The `CRITICAL_ASSESSMENT.md` is unusually candid for an in-progress experiment and correctly identifies most of the issues raised below.

### 2.3 Comprehensive Sensitivity Testing
The experiment includes:
- 8 robustness configurations (KIE modes, lifetime modes, signature sets)
- 1000 Monte Carlo iterations per configuration
- Exchange rate sensitivity analysis
- 5 δD improvement approaches (A–E)
- Fisher information analysis for information content

This is more thorough than most comparable box-model studies.

### 2.4 Code Quality
The codebase is well-structured: shared infrastructure in `common.py` / `models/`, clean separation of analysis scripts, results serialized as both JSON and CSV. The `inputs.py` catalog with literature-referenced presets is excellent practice for reproducibility.

---

## 3. Major Issues (Must Fix)

### 3.1 Prescribed Interhemispheric CH₄ Gradient — Showstopper

**The Problem:**  
The function `realistic_IH_gradient()` in both `improved_model_v3.py` and `validation_phaseBC.py` returns a **hard-coded piecewise-linear ramp** (108 ppb in 2000 → 145 ppb in 2022) based on four anchor points. This is not data — it is an assumption. Since the NH/SH split of total source S depends directly on this gradient, the hemispheric partitioning of emissions is partially predetermined by the analyst's gradient choice.

**Impact:**  
Any trend in NH vs. SH emissions is confounded with the prescribed gradient trend. A reviewer at Nature Communications will immediately flag that the "hemispheric divergence" might simply reflect the assumed gradient divergence.

**Recommendation:**  
Replace with observed NOAA GML zonal-mean CH₄ data (published annually, freely available). The experiment already uses NOAA data for δ¹³C — the CH₄ concentration gradient should come from the same observational network. Propagate observational uncertainty in the gradient through the Monte Carlo.

### 3.2 Absolute FF Emissions Are Unrealistically Low

**The Problem:**  
From `trend_analysis.csv` and the v3 model output, global FF emissions are ~50–60 Tg/yr at the median. The EDGAR v7 inventory places FF at ~110 Tg/yr. This is a factor-of-two discrepancy that the authors acknowledge (CRITICAL_ASSESSMENT.md, Issue #3).

**Impact:**  
While the paper could pivot to a "relative trends only" framing, the spatial aliasing argument requires comparing *absolute* hemispheric trends to 3D inversion posteriors. If the absolute scale is off by 2×, the trend comparison loses meaning — a +1 Tg/yr² trend on a 50 Tg/yr base is very different from +1 Tg/yr² on a 110 Tg/yr base.

**Recommendation:**  
At minimum, add an informative prior on total FF (e.g., Gaussian centered on EDGAR with σ ≈ 20 Tg/yr) and demonstrate that the hemispheric divergence signal survives. Alternatively, present results as fractional changes (% per year) rather than absolute Tg/yr² trends, and discuss the absolute-level bias explicitly.

### 3.3 Weighting Matrix Is Ad Hoc and Result-Driving

**The Problem:**  
In the v1/v2 models, `W_NH = diag(100, 1, 0.5)` and `W_SH = diag(200, 1, 0.5)`. The mass-balance row receives 100–200× more weight than the isotope rows, effectively reducing the system to a mass-balance-only solver with isotopic "nudging." The δD constraint (weight 0.5 vs. 100) is nearly ignored.

The v3 delta-space solver uses `scale = [1.0, 1/50, 1/250]` which is better motivated but still somewhat arbitrary.

**Impact:**  
The claim that this is a "dual-isotope" model is undermined if δD contributes negligibly to the solution. The entire value proposition of δ¹³C + δD over δ¹³C alone must be demonstrated quantitatively.

**Recommendation:**  
Derive weights from measurement/representation uncertainties:
- σ(S) ≈ 20 Tg/yr → w₁ = 1/σ²  
- σ(δ¹³C) ≈ 0.05‰ → w₂ = 1/σ²  
- σ(δD) ≈ 3‰ → w₃ = 1/σ²  

Show results for both the current weighting and the uncertainty-based weighting. Report the δD contribution to the cost function — if it is <5%, acknowledge this.

### 3.4 The Spatial Aliasing Hypothesis Is Not Supported by the Results

**The Problem:**  
The `aliasing_test.json` output reads:
```json
{
  "aliasing_detected": false,
  "hypothesis_supported": false
}
```

The experiment's own diagnostic says the central hypothesis failed. Specifically:
- NH FF trend: +1.04 Tg/yr² (75% positive — not robustly significant)
- Global FF trend: +2.10 Tg/yr² (89% positive — *increasing*, not stable/declining)
- 1-box FF trend: −0.82 Tg/yr² (75% negative — opposite sign)

The 2-box *global aggregate* shows FF increasing, which is the *opposite* of what 1-box models find. This is not "spatial aliasing" but rather a fundamental disagreement between the 1-box and 2-box frameworks.

**Impact:**  
The paper cannot claim "both the 1-box (FF declining) and 3D inversions (FF increasing) are correct" if the 2-box itself shows FF increasing globally. The hypothesis as stated is falsified by the results.

**Recommendation:**  
Either: (a) investigate *why* the 1-box and 2-box disagree even at the global level and make that the paper's contribution, or (b) acknowledge the hypothesis is not supported and pivot to a more defensible claim about what hemispheric resolution adds to box-model source attribution.

---

## 4. Significant Issues (Should Fix)

### 4.1 Post-2019 Edge Effects

Year 2021 shows a +112 Tg/yr NH FF anomaly (noted in CRITICAL_ASSESSMENT.md). The mass-balance method requires M(t+1), making the last year unreliable. The `compute_trends()` function trims the last year (`end_trim=1`), which is good, but this should be stated explicitly in the results narrative and figures should shade the final 1–2 years.

### 4.2 BB ≈ 0 Throughout — Loss of a Degree of Freedom

From `trend_analysis.csv`, BB trends are near zero with symmetric confidence intervals spanning zero in both hemispheres. The SH BB median slope is literally 0 (to machine precision: −3.45e−16). This means the 3-source inversion effectively collapses to a 2-source (FF + Mic) inversion for most Monte Carlo iterations, with the non-negativity bound pushing BB to its lower limit.

**Recommendation:**  
Discuss this explicitly. The FF-BB degeneracy (similar δ¹³C signatures, ~20‰ separation in δD) is a known limitation. Consider whether a 2-source (FF vs. Mic) formulation with BB fixed from GFED4s might give more robust FF/Mic partitioning.

### 4.3 No Posterior Predictive Check in Published Results

The `validation_phaseBC.py` script implements forward-model validation, but I cannot find the output in the `results/` directory (no validation JSON/CSV). The code checks source fractions against EDGAR and computes a δD gradient discrepancy, but these results are not surfaced in RESULTS.md or RESULTS_v3.md.

**Recommendation:**  
Run the validation, save results, and include a dedicated validation section showing: (1) reconstructed vs. observed δ¹³C and δD time series, (2) residual RMSE, (3) χ² goodness-of-fit.

### 4.4 Robustness Test Is Self-Referential

All 8 robustness configurations share the same:
- Mass-balance equations
- Linear least-squares solver
- Prescribed IH CH₄ gradient
- Source signature databases

Varying KIE sampling mode or lifetime treatment tests *parametric* sensitivity, not *structural* robustness.

**Recommendation:**  
Add at least one structurally different comparison: a Bayesian MCMC inversion (the README mentions a PyMC branch exists), or compare directly to published hemispheric estimates from 3D inversions (e.g., Basu 2022 posteriors by hemisphere).

---

## 5. Minor Issues

### 5.1 Inconsistent Model Versions Across Scripts
The analysis directory contains `run_models.py`, `improved_model_v2.py`, `improved_model_v3.py`, and multiple figure-generation scripts (v1, v2, v3). It is unclear which version produced the final reported results. The `results/` directory has subdirectories for `onebox_reference`, `twobox_hemi`, `twobox_global_sigs`, `v2_improved`, and `v3_delta_space` — five different result sets.

**Recommendation:**  
Clearly designate the "canonical" result set and version. Archive or move deprecated versions to a `legacy/` folder.

### 5.2 Missing `requirements.txt` / Environment Specification
The experiment depends on `numpy`, `scipy`, `pandas`, `matplotlib`, `openpyxl`, and custom modules from `common.py` and `models/`. No environment file is provided.

### 5.3 Figure Quality
Figures are saved as both PNG and PDF (good). However, from the file listing, there is no figure showing the posterior predictive check, no figure for δD gradient validation, and no figure comparing the model to EDGAR hemispheric estimates — all of which are critical for the narrative.

### 5.4 Seed Hardcoding
`SEED = 42` is used throughout. Good for reproducibility, but the paper should verify that key conclusions are not seed-dependent by reporting results for at least 3 different seeds.

---

## 6. Numerical Concerns

### 6.1 Condition Number Improvement (v3)
The v3 delta-space solver claims condition number ~27 (after row scaling) vs. ~170,000 in v1/v2. This is a genuine improvement, but the scaling factors `[1.0, 1/50, 1/250]` are themselves somewhat arbitrary. The effective conditioning depends on the *weighted* system, including observation uncertainties.

**Suggestion:** Report both the unscaled and optimally-scaled (Jacobi preconditioned) condition numbers.

### 6.2 `lsq_linear` with Bounds
The solver uses `lsq_linear(A, b, bounds=(0, 1))` for fractions. When the unconstrained solution has a negative component (often BB), the bounded solver projects onto the feasible set, which introduces a systematic bias toward the interior of the simplex. This is the likely cause of the BB ≈ 0 result and the low FF absolute levels.

**Suggestion:** Compare with an unconstrained solve + post-hoc clipping, and report how often each source hits its bound.

---

## 7. Specific Questions for Authors

1. Has the IH CH₄ gradient been validated against NOAA surface flask data? What is the discrepancy between the prescribed ramp and actual observations?

2. The reconciliation.json shows `nh_fraction_of_ff_trend: 49.3%` — almost exactly 50/50 NH/SH split of the FF trend. EDGAR shows ~72% of FF in the NH. Why does the 2-box model produce a symmetric FF partition despite asymmetric observations?

3. The 1-box and 2-box disagree on the *sign* of the FF trend. This is not explained by spatial aliasing. What structural difference in the model causes this sign flip?

4. Why is BB pinned to zero in the hemispheric model but ~1 Tg/yr² positive in the 1-box? Is this a consequence of the non-negativity bound interacting differently with the hemispheric constraint?

5. The v3 model was motivated by numerical conditioning issues. Do the v3 results *qualitatively agree* with v2? The `fig_v1_v2_v3_comparison.png` exists but is not discussed in any results document.

---

## 8. Recommendation

**Verdict: Major revision — not ready for submission.**

The scientific question is excellent and timely. The codebase is well-engineered. But three fundamental issues must be resolved:

1. **Replace the prescribed IH gradient with observations** — this alone may change the hemispheric source partition significantly.
2. **Address the FF absolute-level bias** — either via priors or by reframing the paper around relative trends.
3. **Confront the failed aliasing hypothesis honestly** — the results contradict the central claim. Reframe or revise.

The self-critical `CRITICAL_ASSESSMENT.md` already identifies most of these issues, which suggests the authors are aware. The improvement plan outlined there (Phases A–D) would address most of the concerns raised in this review. I recommend executing that plan before manuscript preparation.

---

## Appendix: Files Reviewed

| File | Purpose | Reviewed |
|------|---------|----------|
| `PLAN.md` | Experiment design | ✓ |
| `SUMMARY.md` | Scientific framing & hypothesis | ✓ |
| `RESULTS.md` | v1 results | ✓ |
| `RESULTS_v3.md` | v3 delta-space results | ✓ |
| `CRITICAL_ASSESSMENT.md` | Self-assessment & improvement plan | ✓ |
| `analysis/improved_model_v3.py` | Core v3 solver | ✓ (detailed) |
| `analysis/hemispheric_trends.py` | Trend analysis | ✓ (detailed) |
| `analysis/validation_phaseBC.py` | Validation (Phases B+C) | ✓ (detailed) |
| `analysis/run_models.py` | v1 model runner | ✓ (scanned) |
| `results/trend_analysis.csv` | Trend slopes & significance | ✓ (detailed) |
| `results/aliasing_test.json` | Aliasing hypothesis test | ✓ (detailed) |
| `results/reconciliation.json` | Literature comparison | ✓ (detailed) |
| `README.md` (repo root) | Model taxonomy | ✓ |
