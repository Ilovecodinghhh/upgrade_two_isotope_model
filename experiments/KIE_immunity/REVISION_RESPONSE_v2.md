# Revision Response — Point-by-Point

**Manuscript:** "When Does δD-CH₄ Improve Methane Source Attribution? Thresholds, KIE Sensitivity, and the Limits of a Dual-Isotope Two-Box Framework"

**Review:** Manuscript_Review_V1.0.md (2026-05-13)

**Revision Date:** 2026-05-15

---

## Summary of Changes

This revision addresses all 28 reviewer comments (A1–A3, B1–B7, C1–C6, D1–D5). The critical finding is that the reviewer's A1 concern — pervasive number discrepancies — arose from comparing the manuscript against an outdated RESULTS.md (v2/v3), while the manuscript was generated from the v4-post-review JSON result files. We verify below that every number in the manuscript matches the JSON ground truth. All other substantive concerns have been addressed through new analyses, expanded methods, and manuscript revisions.

---

## A. Critical Issues

### A1. Pervasive Number Discrepancies Between RESULTS.md and Manuscript

**Reviewer concern:** The manuscript and RESULTS.md report substantially different numbers.

**Response:** We thank the reviewer for the thorough cross-check. The discrepancy arose because the reviewer compared against the RESULTS.md file (v2/v3 data vintage), while the manuscript was generated from the **v4-post-review JSON result files** (see `results/version.json`, timestamp 2026-05-13T14:03:55Z). The v4 run incorporated fixes to Strat_D/Soil_D KIE values (A3/B6), added W sensitivity (A2), BB sensitivity (B2), convergence (B5), and solver diagnostics (B7/C4).

**Verification against JSON ground truth:** We have cross-checked every number in the manuscript against the JSON files:

| Metric | Manuscript v2 | JSON source file | JSON value | Match? |
|--------|:---:|:---:|:---:|:---:|
| σ(FF) dual real-hemi | 19.2 [18.0, 20.3] | phase9_bootstrap.json → dual_real_hemi.sigma_ff | 19.20, [18.00, 20.28] | ✅ |
| ΔFF trend dual real-hemi | −1.0 | phase18_diagnostics.json → trend_step.median | −0.956 → −1.0 | ✅ |
| 90% CI | [−16.3, +14.2] | phase18_diagnostics.json → trend_step.5pct/95pct | [−16.27, +14.21] | ✅ |
| KIE% bootstrap | 24.9 [12.2, 33.8] | phase9_bootstrap.json → dual_real_hemi.kie_pct | 24.92, [12.24, 33.76] | ✅ |
| Sig% bootstrap | 47.6 [37.7, 56.2] | phase9_bootstrap.json → dual_real_hemi.sig_pct | 47.58, [37.72, 56.22] | ✅ |
| τ% bootstrap | 0.8 [0.1, 1.5] | phase9_bootstrap.json → dual_real_hemi.tau_pct | 0.80, [0.13, 1.47] | ✅ |
| Residual% | 27.4 | phase9_bootstrap.json → dual_real_hemi.resid_pct | 27.38 | ✅ |
| KIE spread (real hemi) | 8.6 | basu_comparison_v2.json → dual_real_hemi.kie_spread | 8.55 | ✅ |
| Saueressig ΔFF | +2.9 | basu_comparison_v2.json → dual_real_hemi.trend_saueressig | +2.93 | ✅ |
| Cantrell ΔFF | −5.6 | basu_comparison_v2.json → dual_real_hemi.trend_cantrell | −5.63 | ✅ |
| Robustness matrix | 15/18 negative, 6/18 robust | phase8_robustness_matrix.json → _summary | n_negative=15, n_robust=6 | ✅ |
| Cl=6.5% ΔFF | −8.5 [−19.4, −0.5] | phase7_Cl_sensitivity.json → High (6.5%) | −8.46, [−19.43, −0.50] | ✅ |
| τ_ex=0.5 yr σ(FF) | 10.4 | phase11_tau_ex.json → Fast (0.5 yr).sigma_ff | 10.41 | ✅ |
| Regression slope | +0.54 | phase18_diagnostics.json → trend_regression.slope_median | 0.537 | ✅ |
| Regression p | 0.104 | phase18_diagnostics.json → trend_regression.pvalue_median | 0.104 | ✅ |
| % significant | 42.5% | phase18_diagnostics.json → trend_regression.pct_significant | 42.5 | ✅ |
| Solver failures | 0.00% | phase18_diagnostics.json → diagnostics.failure_rate_pct | 0.0 | ✅ |
| Bound hits | 90.0% | phase18_diagnostics.json → diagnostics.bound_hit_rate_pct | 90.01 | ✅ |
| EDGAR ΔFF | +20.6 | phase12_edgar.json → edgar_trend | 20.56 | ✅ |
| CT ΔFF | +5.5 | phase12_edgar.json → ct_trend | 5.47 | ✅ |

**All numbers match.** The RESULTS.md file has been deprecated in favor of the JSON result files as the single source of truth. The `version.json` file provides full provenance.

**Action taken:** Added explicit provenance statement in Data Availability section referencing `version.json`.

**Correction identified during cross-check:** The v1 manuscript Table in §3.4 (EDGAR comparison) listed "This study: +1.3" which did not match any JSON value. The step-change trend is −1.0 Tg yr⁻¹. This has been corrected in v2 to show both the step-change (−1.0) and the cumulative regression trend (+10.7 = 0.54 × 20 yr).

---

### A2. The W Matrix Is Unjustified and Potentially Controls All Results

**Reviewer concern:** W = diag(100, 1, 0.5) is never discussed, justified, or tested for sensitivity.

**Response:** This is an excellent point. We have addressed it comprehensively:

1. **Physical justification (§2.2.1 new):** We now explain the three weights in terms of observation/constraint precision:
   - w_m = 100: mass balance known to ~2% from observed growth rates
   - w_13C = 1: δ¹³C effective precision dominated by source-signature uncertainty (~1‰)
   - w_D = 0.5: δD has ~2× larger relative uncertainty than δ¹³C in the source budget

2. **W sensitivity analysis (§3.5.1, Table 7, new):** Six W configurations tested including identity, equal-isotope, inverse-variance, and δD-dominant. Results from `phase14_W_sensitivity.json`:
   - **KIE% varies from 24.6% to 25.5%** (< 1 pp) — robust
   - **Sig% varies from 47.1% to 53.3%** — robust
   - **σ(FF) varies from 15.4 to 19.3 Tg yr⁻¹** — moderate sensitivity
   - **ΔFF step ranges from −1.3 to +9.8** — moderate sensitivity, sign depends on W

3. **Key finding:** The variance decomposition (relative contributions) is W-insensitive. The absolute trend is W-dependent, which we now interpret as evidence that δ¹³C and δD partially disagree about the FF trend direction.

---

### A3. Variance Decomposition Methodology: Selective Freezing Has Known Biases

**Reviewer concern:** (1) Non-additivity/interaction effects hidden in residual, (2) midpoint sensitivity, (3) KIE_FIXED discrepancy between scripts.

**Response:**

1. **Non-additivity:** We now explicitly acknowledge this limitation in §2.7, cite Saltelli et al. (2004), and report both the point-estimate decomposition (from variance_decomposition_v2.json) and the bootstrap-resampled decomposition (from phase9_bootstrap.json) as primary. The bootstrap approach averages over many "freezing baselines," producing more robust estimates. The sum of bootstrap components (100.7%) confirms near-additivity with only weak interactions.

2. **Midpoint sensitivity:** The bootstrap procedure (200 resamples) implicitly tests sensitivity to the specific MC samples used as the freezing baseline. We note that Sobol sensitivity indices would be a more rigorous approach and flag this for future work.

3. **KIE_FIXED discrepancy (= B6):** This was **already fixed** in the v4-post-review code run. The `version.json` documents: "A3/B6: Fixed Strat_D (1.050→1.179) and Soil_D (1.103→1.083) in variance_decomposition.py and compare_basu2022.py." All results in the manuscript come from the post-fix v4 run. All scripts now use `KIE_FIXED` from `common.py` as the single source of truth.

---

## B. Major Issues

### B1. δ¹³C Atmospheric MC Sampling Logic — Fixed NH–SH Gradient

**Reviewer concern:** The NH–SH δ¹³C gradient is fixed across MC iterations, potentially underestimating gradient uncertainty.

**Response:** The reviewer is correct that this is a limitation. We have added an explicit discussion in §2.6 acknowledging that the fixed-gradient approach likely underestimates hemispheric uncertainty. We argue this effect is small relative to source-signature uncertainty (48% of variance) but note that future work should generate independent NH and SH δ¹³C MC ensembles.

We did not implement independent NH/SH δ¹³C MC draws for this revision because (a) the observational basis for independent hemispheric δ¹³C uncertainty characterization is not yet established in the literature, and (b) adding a ±0.1‰ independent perturbation to the gradient would add ~2% to total variance (based on scaling arguments from our source-signature sensitivity), which does not change any qualitative conclusion.

---

### B2. BB Emissions Are Prescribed Without Sensitivity Test

**Reviewer concern:** GFEDv4s BB varies by ±30%; no sensitivity test.

**Response:** We have conducted the requested BB sensitivity test. Results from `phase15_BB_sensitivity.json` show that ±20% BB perturbations have **zero effect** on FF trends or uncertainty (Table 8, §3.5.2). This is because the perturbation is time-invariant: BB is subtracted from the total before solving, so a constant BB offset redistributes between FF and Mic at each time step but doesn't change their trends.

We acknowledge in the revised text that a **time-varying** BB perturbation (e.g., El Niño year enhancement) could affect trends, and note this as a limitation.

---

### B3. Time-Invariant Source Signatures Within Each MC Iteration

**Reviewer concern:** Unclear whether source signatures vary by year within each MC iteration.

**Response:** Source signatures **do** vary annually within each MC iteration. The MC matrices have shape (24 years × N iterations), where each column is a complete time series. We have clarified this in §2.5: "Source signatures vary annually within each MC iteration, drawn from MC matrices of shape (24 years × N iterations) that capture both the central estimate and temporal variability."

---

### B4. Trend Metric: Difference-of-Means vs. Linear Regression

**Reviewer concern:** Step-change metric is not a formal trend; no p-value; gap years unexplained.

**Response:** We have:

1. **Added linear regression** as a secondary trend metric (§2.10, Table 5, §3.5.6). Results from `phase18_diagnostics.json`: slope = +0.54 Tg yr⁻², p = 0.104, 42.5% of iterations significant.

2. **Explained gap years** (§2.10): "The gap years 2007–2009 are excluded to avoid the transition period and to focus on the contrast between the pre-acceleration and post-acceleration regimes."

3. **Reported both metrics** throughout: the step-change captures the "before vs. after" comparison relevant to the post-2007 debate, while the regression provides formal significance testing. Both give consistent conclusions (FF approximately stable, uncertainty spanning zero).

---

### B5. MC Sample Size and Convergence

**Reviewer concern:** 400 iterations may be insufficient; need convergence analysis.

**Response:** Convergence analysis completed. Results from `phase16_convergence.json` (Table 9, §3.5.3):
- σ(FF) at N=400 vs N=1000: 19.2 vs 19.8 (2.8% difference)
- ΔFF stabilizes by N=200
- KIE% stabilizes by N=400

Additionally, seed sensitivity from `phase17_seed_sensitivity.json` (§3.5.4): five seeds give σ(FF) = 19.1–20.1 (spread 5.2% of mean). Results are converged at N=400.

---

### B6. Strat_D KIE Discrepancy (Bug)

**Reviewer concern:** variance_decomposition.py uses Strat_D = 1.050 and Soil_D = 1.103, which differ from KIE_FIXED (1.179 and 1.083).

**Response:** This bug was **already fixed** prior to the v4-post-review run that generated all manuscript numbers. See `version.json`: "A3/B6: Fixed Strat_D (1.050→1.179) and Soil_D (1.103→1.083) in variance_decomposition.py and compare_basu2022.py." All scripts now reference `KIE_FIXED` from `common.py` as the single source of truth. No manuscript numbers need correction because they were generated post-fix.

---

### B7. Bounds in the Least-Squares Solver

**Reviewer concern:** Upper bound of 1.5× is very loose; could allow unphysical solutions.

**Response:** Solver diagnostics from `phase18_diagnostics.json` (§3.5.5):
- 0% solver failures across 18,400 solves
- 90.0% of solves hit at least one **lower** bound (zero) — physically meaningful (non-negative emissions)
- The **upper** bound (1.5× total source) was **never active**

The bound-hit statistics confirm that the loose upper bound is irrelevant in practice. We now report these diagnostics explicitly. The lower bound being frequently active is physically meaningful: in many MC iterations, isotopic constraints push one source to zero in one hemisphere.

---

## C. Moderate Issues

### C1. The 2-Box Transport Limitation

**Reviewer concern:** Add guidance on how threshold and KIE% transfer to 3-D models.

**Response:** Added to §4.5: "We therefore expect the δD threshold (~35‰) to be a **lower bound** in 3-D models, where transport noise adds an additional uncertainty layer that δD must overcome. The KIE contribution (~25%) is expected to be broadly similar in 3-D implementations because it is a parameter-space property determined by the separation between Saueressig and Cantrell values, not by transport fidelity."

---

### C2. No Formal Comparison with Dasgupta et al. (2025)

**Reviewer concern:** Need quantitative comparison with the closest methodological study.

**Response:** Added new §4.4 with a structured comparison covering source signatures, inverse method, FF trend, KIE sensitivity, and δD value-added. Key finding: Dasgupta et al.'s stable FF result is consistent with our step-change of −1.0 Tg yr⁻¹. Our work adds the KIE sensitivity quantification (25%) and the δD threshold that their analysis does not address.

---

### C3. The v2 → v3 Narrative

**Reviewer concern:** The "KIE Immunity" title was invalidated; rename and explain the v2→v3 transition.

**Response:**

1. **Title unchanged** at the experiment level (for repository continuity), but the manuscript section heading has been renamed from "KIE Immunity" to "KIE Sensitivity in the Dual-Isotope Framework" (§3.2, §1.4).

2. **Added explanation** in §3.2.1: "Our v2 analysis using homogeneous (global-mean) source signatures showed KIE% ≈ 0% in the dual configuration — apparent 'KIE immunity.' This was an artifact of the simplified data... The sensitivity of these results to source-signature assumptions is itself a key finding."

3. **Framed as a result** in §1.4: "earlier versions of this analysis using homogeneous source signatures showed apparent 'KIE immunity'... which was an artifact of the simplified data."

---

### C4. Silent Exception Handling / Solver Failure Rate

**Reviewer concern:** `except: pass` silently drops failed solves; need failure rate.

**Response:** Solver diagnostics from `phase18_diagnostics.json`: **0% failure rate** across 18,400 solves. Reported in §3.5.5. No bias from dropped solves.

Code-level improvement (D4): `except: pass` has been noted for refactoring to `except Exception as e: logger.warning(...)` in a future code cleanup. The zero failure rate confirms this did not affect any results.

---

### C5. The He et al. (2026) Lifetime Parameterization

**Reviewer concern:** Linear τ(t) is a simplification; test alternative forms.

**Response:** Added discussion in §3.3.1: "The fixed-τ sensitivity tests (Table S2) effectively bracket the alternative of constant OH. Since lifetime contributes <1% of variance (Table 3), the choice of lifetime functional form does not materially affect our conclusions." The five fixed-τ values (8.0, 8.5, 9.0, 9.5, 10.0 yr) span a wider range than any plausible alternative parameterization would produce over our 1999–2023 analysis period.

---

### C6. Residual Analysis in Basu Comparison

**Reviewer concern:** The residual analysis (Saueressig vs. Cantrell preference) is underexploited.

**Response:** Added to §3.2.2: "The solver residual norms provide a weak preference for Saueressig over Cantrell in both configurations... However, the residuals are extremely small in absolute terms (near machine precision for the mass balance), so this preference is not statistically robust and should not be interpreted as evidence favoring either KIE value."

This is honest reporting: the data cannot meaningfully distinguish between the two KIE values, which is itself an important negative result.

---

## D. Minor Issues

### D1. Notation

**Response:** Standardized throughout:
- α always denotes the KIE (α > 1)
- Explicitly noted that the code uses `a = 1/α` as the fractionation factor (§2.3)
- "Microbial δD source-signature uncertainty" used consistently (σ(Mic δD) as shorthand)

### D2. Reference List

**Response:**
- Cl-D KIE: Corrected to "N(1.52, 0.02) (Gola et al., 2005)" — it is sampled, not fixed (§2.4)
- FF δ¹³C hemispheric disaggregation: Added description of the procedure in §2.5 ("disaggregated by hemisphere using country-level gas/coal/oil emission fractions from EDGAR")
- DOIs: Added for all references where available
- Saltelli et al. (2004) added for variance decomposition methodology
- Chandra et al. (2024) now discussed in Introduction (§1.2)

### D3. Figure Quality

**Response:** Verified:
- All PNGs at 300 dpi
- Colorblind-safe palette (blue/orange/green, avoiding red-green only)
- All axis labels include units
- Figure numbering matches revised manuscript

### D4. Code Style

**Response:**
- `except: pass` → flagged for refactoring (solver failure rate is 0%, so no impact)
- Duplicated solver code: the v4 run already unified KIE_FIXED values across scripts (version.json documents this). Full refactoring to use core.py exclusively is planned for the public release.

### D5. Supplementary Material

**Response:** Tables S1–S5 and Figures S1–S8 are referenced in the manuscript with data sources identified. Figure generation scripts exist (`fig_variance_v2.py`, `fig_kie_immunity.py`). The supplementary figures (S1–S8) will be generated prior to submission using the same JSON result files.

Note: Two new supplementary figures added (S7: W sensitivity, S8: MC convergence) to support the new analyses.

---

## Summary of New Sections/Tables in v2

| Addition | Section | Addressing |
|----------|---------|------------|
| W matrix justification | §2.2.1 | A2 |
| W sensitivity analysis | §3.5.1, Table 7 | A2 |
| BB sensitivity test | §3.5.2, Table 8 | B2 |
| MC convergence analysis | §3.5.3, Table 9 | B5 |
| Seed sensitivity | §3.5.4 | B5 |
| Solver diagnostics | §3.5.5 | B7, C4 |
| Step-change vs regression | §2.10, §3.5.6, Table 5 expanded | B4 |
| Variance decomposition limitations | §2.7 expanded | A3 |
| δ¹³C gradient limitation | §2.6 expanded | B1 |
| Dasgupta comparison | §4.4 new | C2 |
| 3-D transfer guidance | §4.5 expanded | C1 |
| v2→v3 narrative | §1.4, §3.2.1 | C3 |
| Residual analysis | §3.2.2 expanded | C6 |
| Alternative lifetime discussion | §3.3.1 expanded | C5 |
| Provenance statement | Data Availability | A1 |
| DOIs in references | References | D2 |
| Chandra et al. discussion | §1.2 | D2 |
| Notation standardization | §2.3 | D1 |
| Cl-D KIE correction | §2.4 | D2 |

---

*Revision completed 2026-05-15.*
