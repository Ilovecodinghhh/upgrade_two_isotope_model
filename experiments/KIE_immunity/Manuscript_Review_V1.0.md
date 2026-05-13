# Manuscript Review V1.0 — KIE Immunity Experiment

**Reviewer:** Senior Editor / Lead Reviewer  
**Manuscript:** "When Does δD-CH₄ Improve Methane Source Attribution? Thresholds, KIE Sensitivity, and the Limits of a Dual-Isotope Two-Box Framework"  
**Target Journal:** Atmospheric Chemistry and Physics  
**Repository:** `Ilovecodinghhh/upgrade_two_isotope_model`, branch `three-box`, path `experiments/KIE_immunity/`  
**Review Date:** 2026-05-13  
**Review Plan Reference:** `Review_Plan.md`

---

## Overall Assessment

**Recommendation: Major Revisions**

This manuscript investigates two timely and important questions: (1) under what conditions does δD-CH₄ improve fossil-fuel emission attribution, and (2) how much does the OH-¹³C KIE controversy still matter in a dual-isotope framework. The work is extensive — 13 analysis phases covering variance decomposition, multi-parameter sensitivity, bootstrap confidence intervals, and inventory validation. The codebase is clean, modular, and reproducible.

The δD threshold result (reconciling the contradiction between Thanwerdas et al. 2024 and Riddell-Young et al. 2025) is a genuine novel contribution and, if properly supported, is publishable in ACP.

However, three categories of issues must be addressed: **(A)** critical internal inconsistencies between the results database and the manuscript text that undermine confidence in every reported number; **(B)** an unjustified solver weighting matrix that may control the headline results; and **(C)** methodological gaps that reviewers at Nature-family or top-tier AGU/EGU journals will raise. Each is detailed below.

---

## A. Critical Issues (Must Fix Before Resubmission)

### A1. Pervasive Number Discrepancies Between RESULTS.md and Manuscript

**Severity: 🔴 Fatal — undermines every claim**

The manuscript and RESULTS.md report substantially different numbers for what should be identical model outputs. This is the single most damaging finding in this review, because a reviewer who cross-checks the code will lose all confidence in the reported results.

| Metric | RESULTS.md (v3) | Manuscript | Discrepancy |
|--------|:---:|:---:|:---:|
| σ(FF) dual real-hemi | **17.8** [16.7, 18.9] | **19.2** [18.0, 20.3] | 1.4 Tg/yr (8%) |
| ΔFF trend dual real-hemi | **+4.8** [+3.6, +5.8] | **−1.0** [−2.8, +0.4] | **Sign flip** |
| KIE% bootstrap | **24.8** [11.0, 34.9] | **24.9** [12.2, 33.8] | CI endpoints differ |
| Sig% bootstrap | **40.6** [28.5, 51.1] | **47.6** [37.7, 56.2] | 7 percentage points |
| τ% | **36.7%** (variance decomp) / **0.0%** (bootstrap) | **0.8%** [0.1, 1.5] | Internally inconsistent even within RESULTS.md |
| Residual% | **12.6%** (variance decomp) / **34.9%** (bootstrap) | **27.4%** | Three different values |
| KIE spread (Basu comparison) | **6.9** Tg/yr | **8.6** Tg/yr | 25% difference |
| Saueressig ΔFF | **+6.3** | **+2.9** | Different by 3.4 Tg/yr |
| Cantrell ΔFF | **−0.6** | **−5.6** | Different by 5.0 Tg/yr |
| EDGAR comparison ΔFF | **+4.8** (RESULTS) / **+8.8** (RESULTS v3 header) | **+1.3** | Three different values |
| Robustness matrix summary | **9/18 negative, 1/18 robust** | **15/18 negative, 6/18 robust** | Completely different |
| Cl = 6.5% ΔFF | **−4.7** [−16.1, +5.5] | **−8.5** [−19.4, −0.5] | Different sign on upper CI bound |
| Phase 11 τ_ex = 0.5 yr | σ = **3.9** | σ = **10.4** | Factor 2.7× |

**Diagnosis:** These are not rounding differences — they appear to be from entirely different model runs or code versions. The most likely explanation is that the manuscript was drafted using a different parameterization (perhaps with different W weights, different BB fractions, or a different v3/v4 of the source signatures) than what produced the JSON result files. The `table1.csv` matches RESULTS.md, not the manuscript.

**Required action:**
1. Identify which set of numbers is correct (rerun from clean state if necessary).
2. Regenerate ALL tables, figures, and text from a single definitive model run.
3. Add a `version.json` or commit hash to all result files so provenance is traceable.
4. The manuscript must not be submitted until every number in the text matches the underlying data files.

### A2. The W Matrix Is Unjustified and Potentially Controls All Results

**Severity: 🔴 Critical — may invalidate headline claims**

The 3×3 least-squares system is solved using `lsq_linear(W @ A, W @ B, bounds=(0, S*1.5))` with:

```python
W = np.diag([100.0, 1.0, 0.5])
```

This means:
- The **mass balance** constraint is weighted **100×** more than the **δ¹³C** constraint
- The **δ¹³C** constraint is weighted **2×** more than the **δD** constraint
- Effectively, δD receives 0.5/100 = **0.5% of the weight** assigned to mass balance

This weighting is never discussed, never justified, and never tested for sensitivity in any of the 13 analysis phases. Yet it directly controls:
- How much δD influences the solution (the 45% improvement claim)
- The relative importance of KIE vs. source-signature uncertainty
- Whether the FF trend is positive or negative (since changing W changes how much the δD constraint can "pull" the solution away from the δ¹³C-only answer)

**Specific concerns:**

1. **Physical justification is absent.** In a properly formulated inverse problem, the weight matrix should reflect the inverse of the observation error covariance. The mass balance is exact (by construction from atmospheric observations), so a large weight is defensible — but 100× is arbitrary. The δ¹³C and δD equations have different uncertainty levels (δ¹³C is measured to ±0.05‰, δD to ±1‰ at atmospheric levels, but source signatures differ by ±2‰ and ±15‰ respectively), so equal or inverse-variance weighting would be more principled.

2. **The 0.5 weight on δD is suspicious.** If δD gets half the weight of δ¹³C, the solver will preferentially fit δ¹³C and treat δD as a soft secondary constraint. This artificially limits the δD improvement — meaning the 45% number could be an underestimate if δD deserves equal weight, or an overestimate if the current weight happens to suppress noise that would otherwise degrade the solution.

3. **W likely explains the RESULTS.md / manuscript discrepancy.** If different W values were used in different runs, this alone could produce the observed sign flips and magnitude changes.

**Required action:**
1. Run a W sensitivity analysis: test W = I (identity), W = diag(100, 1, 1), W = diag(100, 2, 1), W = diag(S, 1/σ_δ¹³C, 1/σ_δD) (inverse-variance), and the current diag(100, 1, 0.5).
2. Report the full results table (σ(FF), KIE%, ΔFF, KIE spread) for each W.
3. Either (a) demonstrate that the results are W-insensitive (unlikely), (b) adopt a principled W with physical justification, or (c) present the W-dependent range as a structural uncertainty.
4. Discuss W explicitly in §2.2 of the manuscript.

### A3. Variance Decomposition Methodology: Selective Freezing Has Known Biases

**Severity: 🔴 Critical — affects interpretation of KIE%, Sig%, τ%**

The variance decomposition uses "selective freezing": fix one parameter group at its central value, rerun, and attribute the variance reduction to that parameter. This approach has two well-known problems:

**Problem 1: Non-additivity.** The decomposition computes:
```python
kie_contrib = max(0, var_total - var_no_kie)
sigs_contrib = max(0, var_total - var_no_sigs)
tau_contrib = max(0, var_total - var_no_tau)
residual = max(0, var_total - kie_contrib - sigs_contrib - tau_contrib)
```

If parameters interact nonlinearly (which they do — KIE and source signatures interact through the 3×3 solver), then:
- `kie_contrib + sigs_contrib + tau_contrib` can exceed `var_total` (the `max(0, ...)` on residual handles this by clamping)
- Or they can sum to less than `var_total`, leaving a large "residual" that is actually just interaction effects

The residual % varies wildly across configurations (6% for δ¹³C-only, 85% for dual-offset, 13–35% for dual real-hemi), suggesting strong parameter interactions that the decomposition cannot resolve.

**Problem 2: Midpoint sensitivity.** Fixing KIE at `(1.0039 + 1.0054)/2 = 1.00465` may not give the same variance reduction as fixing at 1.0039 or 1.0054. The "KIE contribution" depends on which value you fix to. In `variance_decomposition.py`, `fix_kie=True` uses the midpoint, but in `core.py`, `fix_kie=True` uses `KIE_FIXED` values — these should be verified to be identical.

Additionally, in `variance_decomposition.py`:
```python
if fix_kie:
    kie_fixed = {'OH_13C': 0.5*(1.0039+1.0054), ..., 'Soil_D': 1.103}
```
while in `core.py`:
```python
if fix_kie:
    kie_base = dict(KIE_FIXED)
```
where `KIE_FIXED` has `'Strat_D': 1.179` but variance_decomposition.py uses `'Strat_D': 1.050`. **These are different values** — the fix_kie conditions in the two scripts freeze to different KIE sets, meaning the variance decompositions from `core.py` (used by phases 5–13) and `variance_decomposition.py` (used for the primary decomposition) are not directly comparable.

**Specifically:**
- `variance_decomposition.py` fix_kie: Strat_D = **1.050**, Soil_D = **1.103**
- `core.py` / `KIE_FIXED`: Strat_D = **1.179**, Soil_D = **1.083**

These are substantially different values (1.050 vs 1.179 for stratospheric D KIE!). This is a bug.

**Required action:**
1. Fix the KIE_FIXED discrepancy between the two scripts — use a single source of truth.
2. Consider replacing selective freezing with Sobol sensitivity indices (first-order + total-order), which correctly handle interactions and don't require midpoint choices. At minimum, add a footnote acknowledging the limitation.
3. Report the "overshoot" (sum of components vs. total) explicitly; don't hide it in the residual.

---

## B. Major Issues (High Priority)

### B1. δ¹³C Atmospheric MC Sampling Logic

**Concern:** The atmospheric δ¹³C MC sampling applies a **global** offset to fixed hemispheric means:

```python
d13C_glob = sample_atm_d13C(data, k, n)           # MC global draw
d13C_off = d13C_glob[:nc] - c13_glob[:nc]          # offset = MC - mean
d13C_NH_MC = c13_NH[:nc] + d13C_off                # NH = fixed NH + global offset
d13C_SH_MC = c13_SH[:nc] + d13C_off                # SH = fixed SH + same offset
```

This means:
- The **NH–SH δ¹³C gradient** is fixed across all MC iterations (it equals `c13_NH - c13_SH` regardless of k)
- Only the **global level** varies between iterations
- In reality, NH and SH δ¹³C have partially independent uncertainties (different station coverage, different calibration histories)
- This likely **underestimates** the uncertainty in the hemispheric gradient, which is what the 2-box model actually exploits

**Impact:** If the gradient uncertainty is underestimated, then σ(FF) is biased low and the "45% improvement" may be too optimistic.

**Required action:** Either (a) generate independent NH and SH δ¹³C MC ensembles (preferred), or (b) add a sensitivity test where the NH–SH gradient is perturbed by ±0.1‰ independently, and discuss the limitation.

### B2. BB Emissions Are Prescribed Without Sensitivity Test

BB from GFEDv4s is treated as known:
```python
BB_NH = data.BB_global_mean * BB_NH_FRACTION   # BB_NH_FRACTION = 0.55
BB_SH = data.BB_global_mean * BB_SH_FRACTION   # BB_SH_FRACTION = 0.45
```

This uses the **mean** across all years (not annual values), and applies fixed 55/45 NH/SH fractions. In reality:
- GFEDv4s annual BB varies by ±30% (e.g., El Niño years)
- GFAS and GFEDv4s can disagree by 20–40% in specific years
- The NH/SH split varies annually (boreal fire years shift it)

Since BB is subtracted from total source before solving for FF and Mic, BB errors propagate directly into FF estimates. A +10 Tg/yr BB overestimate maps to a −10 Tg/yr FF bias.

**Required action:** Add a sensitivity test: ±20% BB perturbation and ±10% NH/SH split perturbation. Report the impact on ΔFF and σ(FF).

### B3. Time-Invariant Source Signatures Within Each MC Iteration

Source signatures are drawn once per iteration and held constant across all 24 years:

```python
sigs = sample_source_signatures_hemi(rng, data, k, n)
# sigs['ff_d13C_NH'] is shape (n,) but drawn from column k of the MC matrix
# The MC matrix has year-varying values, but they come from a single MC draw
```

Looking at the MC matrices: they are shape (24, 1000), where each column is a complete time series. So source signatures **do** vary by year within each iteration — this is actually correct. However, the manuscript text (§2.5) says "hemispheric MC draws from empirical distributions" without clarifying whether these are time-varying. The reader will assume static values unless told otherwise.

**Required action:** Clarify in the manuscript that source signatures vary annually within each MC iteration (if this is indeed the case). If they are static, add a time-varying sensitivity test.

### B4. Trend Metric: Difference-of-Means vs. Linear Regression

The post-2007 FF trend is defined as:
```python
def compute_trend(FF, years):
    """Post-2007 ΔFF: mean(2010–2018) − mean(2000–2006)"""
    FF_s = smooth_5yr(FF)
    pre = np.where((yrs >= 2000) & (yrs <= 2006))[0]
    post = np.where((yrs >= 2010) & (yrs <= 2018))[0]
    return np.nanmean(FF_s[post], axis=0) - np.nanmean(FF_s[pre], axis=0)
```

This is a step-change metric (difference of two period means), not a trend (slope). The distinction matters:
- A step-change captures both trends and one-time shifts
- It's sensitive to the choice of pre/post periods (why 2000–2006 and 2010–2018? why skip 2007–2009?)
- It has no associated p-value or formal significance test
- The 5-year smoothing before computing the step-change reduces effective degrees of freedom

**Required action:**
1. Add a linear regression slope (with standard error and p-value) as the primary trend metric, or justify the step-change approach.
2. Explain the gap years (2007–2009) — presumably to avoid the transition, but this should be stated.
3. Test sensitivity to period boundaries (e.g., 2000–2007 vs. 2008–2020).

### B5. MC Sample Size and Convergence

400 iterations is used throughout. For bootstrap resampling (Phase 9), 200 bootstrap resamples of 400 iterations are used. This is at the low end for reliable variance estimation.

**Required action:**
1. Produce a convergence plot: σ(FF), KIE%, ΔFF as functions of N_iter (100, 200, 400, 800, 1600). If the curves plateau by 400, the current sample size is adequate.
2. Increase to 1000 iterations for the final results if computationally feasible.

### B6. Strat_D KIE Discrepancy (Bug)

As identified in §A3, `variance_decomposition.py` freezes `Strat_D = 1.050` when `fix_kie=True`, while `KIE_FIXED` in `common.py` specifies `Strat_D = 1.179`. Rice et al. (2003) report the stratospheric D KIE as ~1.16–1.19, so 1.179 is the correct value and 1.050 appears to be a typo (perhaps confused with the soil sink).

Similarly, `Soil_D` is 1.103 in variance_decomposition.py but 1.083 in KIE_FIXED. The KIE_DISTRIBUTIONS have `Soil_D: normal(1.083, 0.01)`, confirming 1.083 is intended.

**Impact:** The primary variance decomposition (which produces the headline KIE% and Sig% numbers) uses wrong D KIE values when freezing KIE. This doesn't affect the "full MC" runs (which sample from the correct distributions), but it biases the "fix_kie" variance toward the wrong baseline, corrupting the decomposition.

**Required action:** Fix the hardcoded KIE values in `variance_decomposition.py` to match `KIE_FIXED` from `common.py`, then rerun the decomposition.

### B7. Bounds in the Least-Squares Solver

The solver uses `bounds=(0, S*1.5)`:
- Lower bound 0 prevents negative emissions (physical)
- Upper bound S*1.5 allows any single source to be 150% of total emissions

The upper bound is very loose. In realistic scenarios, FF emissions should not exceed ~200 Tg/yr (vs. total of ~560 Tg/yr). The 1.5× bound means FF could reach ~420 Tg/yr per hemisphere. This is unrealistic and could allow the solver to find unphysical solutions that happen to minimize the weighted residual.

**Required action:** Either tighten the bounds to physically plausible ranges (e.g., FF ∈ [0, 250], Mic ∈ [0, 400], BB already prescribed) or demonstrate that the bounds are rarely active. The `basu_comparison_v2.json` already tracks bound hits — report these statistics in the manuscript.

---

## C. Moderate Issues

### C1. The 2-Box Transport Limitation

The manuscript acknowledges (§4.4) that Naus et al. (2019) identified systematic biases in 2-box models. The τ_ex sensitivity analysis (Phase 11) confirms this: σ(FF) ranges from 3.9 to 23.7 Tg/yr depending on exchange time. This is a factor-of-6 uncertainty from transport alone — larger than the KIE contribution.

However, the manuscript frames the 2-box as a structural analysis tool rather than a competing inversion. This framing is appropriate but needs strengthening:

**Required action:** Add a sentence in §4.4 explicitly stating that the δD threshold (~35‰) is expected to be a lower bound in 3-D models (where transport adds noise), while the KIE contribution (~25%) may be similar because it is a parameter-space property. This guides readers on how to use the results.

### C2. No Formal Comparison with Dasgupta et al. (2025)

Dasgupta et al. (2025) use a Bayesian 2-box framework with δD — the closest methodological comparison to this study. The manuscript cites them but never compares quantitatively. Key questions:
- Do they use similar source signatures?
- What is their implicit W matrix?
- How does their FF trend compare?
- Do they find a similar KIE sensitivity?

**Required action:** Add a paragraph in §4 directly comparing methodology and results with Dasgupta et al. (2025).

### C3. The v2 → v3 Narrative

RESULTS.md documents that the "KIE immunity" finding from v2 (18/18 robustly negative) collapsed to 1/18 in v3 (or 6/18 in the manuscript — another discrepancy). This is scientifically honest but creates a narrative problem: the experiment's own title ("KIE Immunity") was invalidated by the authors' data upgrade.

The manuscript handles this by shifting the emphasis to "KIE still matters" and "source signatures are the binding constraint." This is the right approach but needs to be more explicit:

**Required action:**
1. Consider renaming the experiment/section from "KIE Immunity" to something more accurate (e.g., "KIE Sensitivity in the Dual-Isotope Framework").
2. Add a sentence in §3.2.1 explicitly noting that earlier versions with homogeneous source signatures showed KIE immunity, which was an artifact of the simplified data.
3. Frame the v2→v3 transition as a result itself: "The sensitivity of these results to source-signature assumptions is itself a key finding."

### C4. Silent Exception Handling

In the solver loop:
```python
try:
    res = lsq_linear(W@A, W@B, bounds=(0, S*1.5))
    FF_G[j, k] += res.x[1]
except:
    pass
```

Failed solver calls silently produce zero FF for that hemisphere-year-iteration. If failures are common, this biases FF downward. If rare, it's harmless — but we don't know which.

**Required action:** Log and report the failure rate. If >1%, discuss the bias.

### C5. The He et al. (2026) Lifetime Parameterization

The lifetime model is:
```python
return 9.0 - 0.017 * (years - 2010)
```

This gives τ = 9.19 yr in 1999 and τ = 8.81 yr in 2021 — a monotonic decline. The manuscript attributes this to He et al. (2026a), who use TROPOMI observations. However:
- The linear parameterization is a gross simplification of what is likely a nonlinear process
- It implies OH has been increasing monotonically, which is debated (Zhao et al. 2023)
- The sensitivity tests (Phase 5) only test fixed lifetimes, not alternative functional forms

**Required action:** Test an alternative lifetime parameterization (e.g., constant at 9.0 yr until 2007, then declining) to verify that the linear form doesn't control the trend calculation.

### C6. Residual Analysis in Basu Comparison

The Basu comparison (`compare_basu2022.py`) includes a residual analysis to determine which KIE "fits better." This is a powerful idea but is underexploited:

```python
if mean_res_C < mean_res_S:
    preferred = "Cantrell"
else:
    preferred = "Saueressig"
```

**Required action:** Report the residual analysis results in the manuscript. If the data genuinely prefer one KIE over the other, this is a significant finding. If residuals are indistinguishable, that itself is informative.

---

## D. Minor Issues

### D1. Notation

- §2.3 uses both `α` and `K` for the KIE. Standardize to one symbol.
- The code uses `a13_NH = 1.0/K13_NH` — i.e., `a` is the fractionation factor (= 1/KIE). The manuscript uses `α` for the KIE itself. Ensure consistency.
- "σ(Mic δD)" and "Mic δD uncertainty" are used interchangeably — pick one.

### D2. Reference List

- Gola et al. (2005) is cited for Cl-D KIE = 1.52, but the Cl-D KIE in `KIE_DISTRIBUTIONS` is `normal(1.52, 0.02)`. The manuscript says "1.52 (fixed; Gola et al., 2005)" — but the code samples it. Clarify whether it's fixed or sampled.
- Sherwood et al. (2017) is cited for hemispheric source-signature disaggregation, but the actual disaggregation procedure (using country-level gas/coal/oil fractions) is novel to this work. This should be described more explicitly.
- DOIs are missing for most references.

### D3. Figure Quality

Without rendering the figures, I note that `fig_variance_v2.py` and `fig_kie_immunity.py` exist with PNG/PDF outputs. Verify:
- 300 dpi minimum for PNG
- Colorblind-safe palette (avoid red-green only)
- Axis labels include units
- Figure numbers match manuscript text

### D4. Code Style

- `except: pass` should be `except Exception: pass` at minimum, or better, `except np.linalg.LinAlgError`
- Several scripts duplicate the entire 2-box solver (variance_decomposition.py, compare_basu2022.py) instead of using `core.py`. This creates divergence risk (as seen with the Strat_D bug). Refactor to use `core.py` everywhere.

### D5. Supplementary Material

The manuscript references Tables S1–S5 and Figures S1–S6, but these don't exist as files. They need to be generated before submission.

---

## E. Summary of Required Actions

### Priority 1 — Blocking (Week 1)

| # | Action | Manuscript Section |
|---|--------|-------------------|
| A1 | Resolve all number discrepancies; regenerate from single model run | All |
| A2 | W matrix sensitivity analysis + justification | §2.2 (new subsection) |
| A3 | Fix Strat_D/Soil_D bug in variance_decomposition.py; consider Sobol indices | §2.7, Table 3 |
| B6 | Fix KIE_FIXED discrepancy (same as A3) | Code |

### Priority 2 — Major Revisions (Week 2)

| # | Action | Manuscript Section |
|---|--------|-------------------|
| B1 | Address δ¹³C atmospheric MC gradient assumption | §2.6, Supplementary |
| B2 | BB sensitivity test (±20%) | §3.3 (new subsection) |
| B4 | Add linear regression trend + p-value | §3.2.3, all trend results |
| B5 | MC convergence analysis | Supplementary |
| B7 | Report solver bound-hit statistics | §2.2 or Supplementary |

### Priority 3 — Strengthening (Week 3)

| # | Action | Manuscript Section |
|---|--------|-------------------|
| C1 | Strengthen 2-box limitation framing | §4.4 |
| C2 | Quantitative comparison with Dasgupta et al. (2025) | §4 (new paragraph) |
| C3 | Rename "KIE Immunity"; frame v2→v3 as a finding | Title, §3.2 |
| C4 | Report solver failure rate | §2.2 or Supplementary |
| C5 | Alternative lifetime parameterization test | §3.3 or Supplementary |
| C6 | Report Basu residual analysis | §3.2.2 |

### Priority 4 — Minor (Week 4)

| # | Action |
|---|--------|
| D1 | Standardize notation |
| D2 | Fix reference details |
| D3 | Verify figure quality |
| D4 | Refactor duplicated solver code |
| D5 | Generate supplementary tables/figures |

---

## F. Specific Comments on Manuscript Sections

### Abstract

The abstract is well-structured but contains numbers that don't match the body text. Key claims:
- "45% [42%, 48%]" — RESULTS.md says 45.1% [41.7, 47.6] (from the dD_threshold experiment, not KIE_immunity). The KIE_immunity σ(FF) reduction is (31.0 − 17.8)/31.0 = 43%, not 45%. Clarify which experiment produces the 45% number.
- "25% [12%, 34%]" — RESULTS.md says 24.8 [11.0, 34.9]; manuscript Table 3 says 24.9 [12.2, 33.8]. These should be identical.
- "8.6 Tg yr⁻¹" KIE spread — RESULTS.md says 6.9. Different runs.
- "Saueressig: +2.9; Cantrell: −5.6" — RESULTS.md says +6.3 and −0.6. Different runs.

**The abstract cannot be evaluated until the number discrepancies are resolved (A1).**

### §1 Introduction

Well-written and comprehensive. The literature review correctly identifies the key tension between Thanwerdas et al. (2024) and Riddell-Young et al. (2025). The scope statement is clear. Two suggestions:

1. Add Chandra et al. (2024) as another recent dual-constraint study (already in reference list but not discussed).
2. The introduction promises to "reconcile" the contradiction — this is a strong claim. The manuscript delivers on it for the δD threshold, but the KIE findings are more nuanced. Consider softening to "investigate" or "provide a framework for understanding."

### §2 Methods

Generally clear. Specific comments:

- **§2.1:** The mass balance equation is standard. However, the phrase "total source for each hemisphere-year is computed from the observed CH₄ growth rate" obscures an important detail: the total source S depends on the assumed lifetime τ, so there is a circularity when τ is also treated as uncertain. This should be noted.

- **§2.2:** The over-determined system description is accurate but incomplete. State explicitly that the system is solved via bounded weighted least squares, name the algorithm (`scipy.optimize.lsq_linear`), and report W. Currently W is absent from the methods section entirely.

- **§2.3:** The source-weighted isotopic ratio derivation uses `a = 1/K` where K is the bulk KIE. Verify that this is the correct convention — some references define α as the fractionation factor (= 1/K) and others define α as the KIE itself (= K). The code computes `a13_NH = 1.0/K13_NH` and uses it as the fractionation factor in the numerator of the isotope budget, which is correct for the convention α = K.

- **§2.4:** Table of MC parameter ranges is useful. However, `Strat_D = 1.179 ± 0.01` disagrees with the fixed value used in variance_decomposition.py (1.050). See B6.

- **§2.5, Table 1:** The hemispheric source signatures are the key data input. The manuscript should provide more detail on how the FF δ¹³C gap (NH −44‰ vs. SH −48.5‰) was derived — this 4.5‰ gap drives much of the model behavior. Reference the build scripts or provide the gas/coal/oil fractions per hemisphere.

- **§2.7:** The variance decomposition method should cite prior use (e.g., Saltelli et al. 2004 for factor fixing) and explicitly acknowledge the interaction/non-additivity issue.

### §3 Results

- **§3.1:** The δD threshold analysis is the strongest section. The crossover at ~35‰ is clearly presented and the reconciliation of Thanwerdas/Riddell-Young is compelling. However, verify that the 45% improvement number comes from the dD_threshold experiment (not KIE_immunity), and clarify this in the text.

- **§3.2:** The variance decomposition and Basu comparison results depend entirely on which numbers are correct (A1). Cannot evaluate further until resolved.

- **§3.3:** The sensitivity analyses are comprehensive but all suffer from the step-change trend metric (B4). The Cl fraction sensitivity (Table 6) is a particularly strong result — the sign flip at Cl ≈ 5% is a clear and important finding.

- **§3.4:** The EDGAR comparison has three different ΔFF values across the repository (see A1). Additionally, the comparison is only qualitative — a more rigorous comparison would compute the correlation or RMSE between the model and EDGAR annual time series, not just the trends.

### §4 Discussion

- **§4.1:** The "diagnostic, not deterministic" framing is appropriate.
- **§4.2:** The claim that "a new laboratory measurement of the OH-¹³C KIE remains the single highest-priority experiment" is powerful but requires that the KIE% finding survives the W sensitivity test (A2) and the Strat_D bug fix (B6).
- **§4.3:** Good identification of source signatures as the next frontier.
- **§4.4:** The Naus et al. (2019) comparison needs the quantitative strengthening described in C1.
- **§4.5:** The discussion of OH trends is adequate.
- **§4.6:** The policy implications paragraph is measured and appropriate — neither over-claiming nor dismissive.

### §5 Conclusions

The five conclusions are clear and well-ordered. However:
- Conclusion 1 (45% improvement) comes from the dD_threshold experiment, not KIE_immunity
- Conclusion 2 (25% KIE) depends on the unresolved number discrepancies
- Conclusion 4 (Cl ≥ 6.5%) is robust and well-supported
- Conclusion 5 (new KIE measurement) is contingent on conclusion 2

---

## G. Verdict

This manuscript addresses a genuinely important question with a technically competent and extensive analysis. The δD threshold result is novel, clearly presented, and reconciles an active debate in the literature. The KIE sensitivity analysis, while compromised by data inconsistencies and the W matrix issue, would be a significant contribution if cleaned up.

**If the authors resolve the A-level issues (number consistency, W justification, KIE bug), this paper is suitable for ACP.** The B-level issues are standard major-revision requests. The C-level issues would strengthen the paper for a higher-impact venue.

The single most concerning finding in this review is **A1** — the pervasive number discrepancies suggest that the manuscript and the codebase have diverged, and neither set of numbers can currently be trusted. This must be the first thing fixed.

**Estimated revision effort:** 2–3 weeks for A+B items; 4 weeks for A+B+C.

---

*Review completed 2026-05-13. Based on Review_Plan.md §§2–6.*
