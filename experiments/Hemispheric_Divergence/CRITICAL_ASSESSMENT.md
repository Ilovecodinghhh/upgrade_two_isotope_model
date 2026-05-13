# Hemispheric Divergence — Critical Assessment & Improvement Plan

**Date:** 2026-05-13
**Objective:** Raise this experiment to publication quality (Nature Comm / ACP level)

---

## Honest Assessment: What's Wrong

### 🔴 Showstoppers (must fix)

#### 1. The "degeneracy breaking" claim is mathematically unsupported

The core claim — that the 2-box breaks a degeneracy the 1-box suffers from — is the paper's thesis. But the **condition numbers are nearly identical** (170k NH vs 180k global, ratio 1.06×). The A matrices in both 1-box and 2-box are equally ill-conditioned because:

- The 3 source signatures differ primarily in δ¹³C, not δD
- δD rows are O(10⁻⁴), δ¹³C rows are O(10⁻²), mass row is O(1) → 5 orders of magnitude span
- The weighting matrix W=diag(100, 1, 0.5) makes this WORSE (cond → 33M)

**What actually differs** between 1-box and 2-box is not conditioning but the **interplay between hemispheric total source S and isotopic mass balance**. The 2-box solves two coupled 3×3 systems with different S_NH and S_SH (from the hemispheric CH₄ gradient), which constrains the solutions differently.

**Fix:** Replace the "degeneracy breaking" narrative with a rigorous analysis of *why* the 2-box partitions differently. Compute the effective degrees of freedom, the sensitivity of each source to each isotope, and show that hemispheric coupling (exchange term) provides genuine additional constraint.

#### 2. CH₄ hemispheric gradient is PRESCRIBED, not observed

`compute_IH_gradient()` returns a linear ramp from 80 to 100 ppb — hard-coded, not from data. This means:
- The NH/SH split of total source S is assumption-driven
- Any trend in S_NH vs S_SH is partially an artifact of the prescribed gradient trend
- Reviewers will immediately flag this

**Fix:** Use actual NOAA GML zonal mean data. The Global Monitoring Laboratory publishes annual mean CH₄ by latitude (30-60°N, 0-30°N, etc.) — this is standard data. Alternatively, use the NH/SH flask means from AGAGE or NOAA surface network.

#### 3. FF is far too low (50 Tg/yr vs literature 100-120)

The model produces FF = 50-60 Tg/yr, roughly half the accepted value. This undermines credibility:
- EDGAR v7 FF emissions: ~110 Tg/yr
- Our model: 50-60 Tg/yr
- This suggests the isotopic inversion is not properly constrained

**Root cause:** The `lsq_linear` with non-negativity bounds and the ill-conditioned A matrix pushes solutions toward Mic-dominated answers because the microbial isotopic signature is furthest from the atmospheric mean.

**Fix:** Either (a) add an EDGAR-based prior that pulls FF toward ~100 Tg/yr, or (b) switch to Bayesian inversion with informative priors (Keeling-plot style), or (c) acknowledge this is a relative trend analysis and focus exclusively on trends, not absolute values.

#### 4. The W matrix is arbitrary and drives results

W_NH = diag(100, 1, 0.5) and W_SH = diag(200, 1, 0.5). The 100× and 200× weights on the mass row mean the solver prioritizes matching total S over matching isotopic compositions. This:
- Guarantees FF+Mic+BB ≈ S (good)
- But gives almost no weight to δD (0.5 vs 100) → the δD constraint is nearly ignored
- The entire dual-isotope advantage may be illusory with these weights

**Fix:** Derive weights from measurement uncertainties: σ(S) ≈ 20 Tg/yr → w₁ = 1/20² = 0.0025; σ(δ¹³C) ≈ 0.05‰ → w₂ = 1/0.05² = 400; σ(δD) ≈ 3‰ → w₃ = 1/3² = 0.11. This would properly balance the constraints.

### 🟡 Significant Issues (should fix)

#### 5. Post-2019 edge effects corrupt the last 2 years

Year 2021 shows NH_FF = +112 Tg/yr anomaly — clearly an artifact. The mass-balance method becomes unstable at the end of the time series because it needs M(t+1). The last year's result is always unreliable.

**Fix:** Trim the last year from all trend analyses. Flag this explicitly.

#### 6. δD gradient self-inconsistency is a real problem

The model predicts a −28‰ NH-SH δD gradient but uses −14.6‰ as input. If the source decomposition is right, the input data is wrong. If the input data is right, the source decomposition is wrong. Either way, something doesn't add up.

**Fix:** Check whether the observed hemispheric δD gradient from Dasgupta calibration files actually implies −14.6‰. If so, the source δD values may need revision. This is a consistency check that belongs in the paper.

#### 7. No model validation against independent data

We compare trends against literature qualitatively, but never validate the model's basic predictions:
- Does the model reproduce the observed global δ¹³C trend?
- Does it reproduce the observed global δD trend?
- Does the 2-box reproduce the observed IH δ¹³C gradient of −0.24‰?
- What are the residuals?

**Fix:** Add posterior predictive checks: feed the solved sources back through the mass balance and compare modeled vs observed isotope time series.

#### 8. Only 3 source categories is too coarse

FF+Mic+BB misses important source-types with distinct isotopic signatures:
- Coal vs gas vs oil (different δ¹³C within "FF")
- Wetlands vs rice vs ruminants (different δD within "Mic")
- The 3-category approach is standard but limits what the model can constrain

**Fix (for Discussion):** Acknowledge this limitation explicitly. Note that 3D inversions can use spatial distribution as a proxy for sub-source attribution.

### 🟢 Minor Issues (nice to fix)

#### 9. The robustness test is self-referential

Testing 8 variants of the same model doesn't prove robustness to model structure. All 8 share:
- Same mass-balance equations
- Same linear solve approach
- Same prescribed IH CH₄ gradient
- Same source signature datasets

**Fix:** Add at least one structurally different comparison: e.g., a Bayesian MCMC inversion, or comparison to existing published hemispheric estimates.

#### 10. No proper uncertainty propagation

The 90% CI from MC iterations captures parametric uncertainty but not:
- Structural uncertainty (model choice)
- Data uncertainty in IH gradient (it's prescribed!)
- Systematic errors in source signatures

**Fix:** At minimum, document what the MC captures and what it doesn't.

---

## Improvement Plan (Priority Order)

### Phase A: Fix the data foundation (2 hours)

1. **Replace prescribed IH CH₄ gradient with observed data**
   - Download NOAA GML zonal means (already in the dataset!)
   - Or use the NH/SH means from the data files we already have
   - Rerun all models with real hemispheric CH₄

2. **Fix the W matrix to be uncertainty-based**
   - w_mass = 1/σ²(S), w_d13C = 1/σ²(δ¹³C_src), w_dD = 1/σ²(δD_src)
   - Rerun with proper weighting

3. **Trim last year from all analyses**

### Phase B: Validate the model (1 hour)

4. **Posterior predictive check**
   - For each MC iteration, reconstruct atmospheric δ¹³C and δD from solved sources
   - Compare to observations
   - Report χ² or RMSE

5. **δD consistency check**
   - Compare predicted IH δD gradient to observed
   - Diagnose the −28‰ vs −14.6‰ discrepancy

### Phase C: Strengthen the narrative (2 hours)

6. **Replace "degeneracy breaking" with proper information-theoretic analysis**
   - Compute effective information gain from adding hemispheric resolution
   - Show that the 2-box has MORE independent constraints (6 observations vs 3)
   - The system goes from underdetermined (3 unknowns, ~2.5 effective constraints) to better-determined (6 unknowns, ~5 effective constraints)

7. **Add EDGAR cross-check**
   - Compare NH/SH FF partition to EDGAR 
   - If our model gets 50% NH (vs EDGAR 72%), discuss why

8. **Add GFED BB validation**
   - Compare our BB trend to GFED4s fire emissions
   - This strengthens the "BB as discriminator" argument

### Phase D: Rerun and finalize (2 hours)

9. **Rerun all models with fixes from Phase A**
10. **Regenerate all figures**
11. **Update RESULTS.md with validated findings**
12. **Rewrite the narrative around the corrected results**

---

## What This Paper Can Actually Claim (After Fixes)

The defensible contribution is:

> "A 2-box dual-isotope model provides 6 independent constraints on methane sources (3 per hemisphere × 2 hemispheres) compared to 3 in a 1-box model. With hemispheric source signatures and atmospheric observations, the 2-box consistently attributes more methane growth to fossil fuels and less to biomass burning than the equivalent 1-box, in qualitative agreement with spatially resolved 3D inversions. The key difference is the FF/BB partition, which is poorly constrained in 1-box models due to similar isotopic signatures of these two source categories."

This is a valid, publishable contribution. But it needs the data and methodology fixes above to be defensible.
