# Review Plan — KIE Immunity Experiment

**Reviewer Role:** Senior Editor / Lead Reviewer (Nature Communications / Global Biogeochemical Cycles level)  
**Date:** 2026-05-13  
**Repository:** `Ilovecodinghhh/upgrade_two_isotope_model`, branch `three-box`, path `experiments/KIE_immunity/`  
**Manuscript:** `MANUSCRIPT_DUAL_ISOTOPE.md` (draft targeting ACP)

---

## 1. Executive Summary of Claims

The experiment makes five headline claims:

| # | Claim | Key Numbers |
|---|-------|-------------|
| C1 | δD improves FF attribution by ~45% in a 2-box framework, but only below a critical microbial δD uncertainty threshold (~35‰) | CI reduction 105→58 Tg/yr; crossover at 3.82× baseline σ |
| C2 | The OH-¹³C KIE controversy still accounts for ~25% of FF variance even with dual isotopes, and determines the **sign** of the post-2007 FF trend | KIE% = 24.8 [11.0, 34.9]; Saueressig +2.9 vs Cantrell −5.6 Tg/yr |
| C3 | Source signatures dominate the remaining uncertainty (~48%) | Sig% = 47.6 [37.7, 56.2] |
| C4 | The Cl fraction controls robustness — at Cl ≥ 6.5%, negative FF trend becomes robust | Sign flip at Cl ≈ 5% |
| C5 | A new OH-¹³C KIE lab measurement is the single highest-priority experiment | KIE spread = 8.6 Tg/yr (34% reduction from Basu's 13.0, but not elimination) |

---

## 2. Review Dimensions

### 2.1 Scientific Soundness

| Item | Question | Priority | Files to Examine |
|------|----------|----------|-----------------|
| **2.1.1** Model formulation | Is the 2-box mass balance correctly implemented? Are the isotopic mass-balance equations (§2.1–2.3) consistent with Naus et al. (2019) and standard formulations? | 🔴 Critical | `core.py`, `common.py` |
| **2.1.2** Variance decomposition methodology | Is "selective freezing" a valid approach? Does fixing one parameter at its midpoint correctly isolate its contribution, or do nonlinear interactions bias the decomposition? (Components don't have to sum to 100% if interactions exist — is the residual correctly capturing this?) | 🔴 Critical | `variance_decomposition.py`, `core.py` (fix_kie, fix_sigs) |
| **2.1.3** Solver choice & weighting | The 3×3 system uses `lsq_linear` with W = diag(100, 1, 0.5). Why these specific weights? How sensitive are results to W? The mass constraint is weighted 100× more than isotope constraints — is this justified or does it suppress δD's information content? | 🔴 Critical | `core.py` line `W = np.diag([100.0, 1.0, 0.5])` |
| **2.1.4** δ¹³C atmospheric MC sampling | The code applies a global MC offset to hemispheric means (`d13C_off = d13C_glob[:nc] - c13_glob[:nc]`). Does this correctly propagate hemispheric atmospheric δ¹³C uncertainty, or does it impose a correlated shift that underestimates the NH–SH gradient uncertainty? | 🟡 High | `core.py` lines 100–107 |
| **2.1.5** Over-determined vs. under-determined | The manuscript says "3 equations, 2 unknowns" (over-determined) for dual-isotope. But with BB prescribed, it's actually 3 eqs for 2 unknowns per hemisphere — the system is well-posed. Is this correctly stated and does the bounded least-squares handle edge cases (negative solutions, hitting bounds)? | 🟡 High | `core.py`, `MANUSCRIPT_DUAL_ISOTOPE.md` §2.2 |
| **2.1.6** Lifetime parameterization | `compute_lifetime()` implements He 2026's time-varying τ. Is this correctly translated from the original paper? Are the fixed-τ sensitivity tests meaningful given that τ interacts with KIE through the sink term? | 🟡 High | `common.py` (not in this dir) |
| **2.1.7** Interhemispheric exchange formulation | The exchange flux `(M_SH - M_NH)/tau_ex` should conserve mass — verify that the NH and SH terms are equal and opposite, and that the isotopic exchange fluxes are consistent. | 🟡 High | `core.py` lines 85–92 |

### 2.2 Statistical Rigor

| Item | Question | Priority | Files to Examine |
|------|----------|----------|-----------------|
| **2.2.1** MC sample size | 400 iterations — is this sufficient for convergence of the variance decomposition and trend statistics? Bootstrap uses 200 resamples of 400 — what are the bootstrap CIs on the CIs? | 🟡 High | `core.py`, `phase9_bootstrap_variance.py` |
| **2.2.2** Trend definition | "Post-2007 ΔFF" = mean(2010–2018) − mean(2000–2006). This is a difference-of-means, not a regression slope. Is this the right metric? A linear trend fit would be more standard and provide a p-value. | 🟡 High | `core.py` `compute_trend()` |
| **2.2.3** 5-year smoothing | `smooth_5yr` applied before trend calculation — does this reduce effective degrees of freedom? Is the smoothing window appropriate? Does it alias the 2007 breakpoint? | 🟢 Medium | `core.py`, `common.py` |
| **2.2.4** Seed dependence | All runs use `seed=42`. Has robustness to seed choice been tested? | 🟢 Medium | All analysis scripts |
| **2.2.5** Bootstrap validity | Bootstrap resampling of MC iterations assumes iterations are exchangeable (iid). If there are correlations between iterations (e.g., through shared atmospheric observations), the bootstrap CIs may be too narrow. | 🟡 High | `phase9_bootstrap_variance.py` |

### 2.3 Data Integrity

| Item | Question | Priority | Files to Examine |
|------|----------|----------|-----------------|
| **2.3.1** Hemispheric δ¹³C source signatures | These are the v3 addition that "changed everything." Verify the build pipeline: how are NH/SH FF, BB, Mic δ¹³C signatures constructed? Are the underlying datasets (Sherwood 2017, Luo 2024, isotem) correctly processed? | 🔴 Critical | `common.py`, `rel/build_*` scripts, MC CSV files |
| **2.3.2** Hemispheric δD atmospheric data | Station-level MC ensembles — how many stations per hemisphere? What is the spatial coverage? Are SH stations (sparse) introducing bias? | 🟡 High | `common.py` (`sample_atm_dD_hemi`) |
| **2.3.3** BB emissions from GFEDv4s | BB is prescribed, not solved. How sensitive are results to the BB time series? A sensitivity test with ±20% BB perturbation is missing. | 🟡 High | Data loading in `common.py` |
| **2.3.4** EDGAR comparison | The EDGAR FF trend (+20.6 Tg/yr) is compared to the model's +4.8 (RESULTS.md) or +1.3 (manuscript) — these numbers differ. Which is correct? | 🟡 High | `phase12_edgar.json`, `MANUSCRIPT_DUAL_ISOTOPE.md` §3.4 |

### 2.4 Internal Consistency

| Item | Question | Priority | Files to Examine |
|------|----------|----------|-----------------|
| **2.4.1** RESULTS.md vs. MANUSCRIPT numbers | Several numbers disagree between RESULTS.md and the manuscript draft. E.g.: | 🔴 Critical | Both files |
| | — RESULTS §Phase 9: KIE% = 24.8 [11.0, 34.9] vs. Manuscript Table 3: 24.9 [12.2, 33.8] | | |
| | — RESULTS: σ(FF) = 17.8 vs. Manuscript Table 3: 19.2 | | |
| | — RESULTS: ΔFF = +4.8 vs. Manuscript Table 5: −1.0 | | |
| | — RESULTS: KIE spread = 6.9 vs. Manuscript: 8.6 | | |
| | These suggest RESULTS.md and the manuscript may reflect different model runs or versions. **This must be resolved before any external review.** | | |
| **2.4.2** v2 → v3 labeling | The manuscript text doesn't mention version numbers (v1/v2/v3), but RESULTS.md documents major changes between versions. Are all manuscript numbers from v3, or is there a mix? | 🟡 High | All |
| **2.4.3** Figure–text agreement | Do the figures (`fig_variance_v2.png`, `fig_edgar_validation.png`, `fig_kie_immunity.png`) show the same numbers as the manuscript tables? | 🟡 High | `figures/` |

### 2.5 Methodological Concerns for Peer Review

| Item | Concern | Priority |
|------|---------|----------|
| **2.5.1** **W matrix arbitrariness** | The diagonal weight matrix W = diag(100, 1, 0.5) is not justified anywhere. W controls the relative importance of mass vs. δ¹³C vs. δD constraints. With mass weighted 200× more than δD, the solver will strongly enforce total emissions and treat isotopes as soft constraints. A reviewer will ask: (a) what happens with W = I (identity)? (b) what if δD is weighted more heavily? (c) is there a principled way to set W (e.g., inverse-variance weighting)? **This is potentially the most important unaddressed issue.** | 🔴 Critical |
| **2.5.2** **BB not solved** | Prescribing BB from GFEDv4s removes one degree of freedom but also removes a source of uncertainty. If BB is uncertain (and it is — GFEDv4s vs. GFAS can differ by 30%), how much does this affect FF? | 🟡 High |
| **2.5.3** **No formal model selection** | The 2-box model is compared to a 1-box (implicitly, via the "offset" config) and the manuscript claims the 2-box is better. But there's no BIC/AIC/DIC comparison, no formal information criterion. | 🟢 Medium |
| **2.5.4** **Time-invariant source signatures** | Source signatures are drawn per-MC-iteration but are constant across years within each iteration. In reality, the FF δ¹³C signature has likely shifted as the coal/gas mix changed. How sensitive is the trend to allowing time-varying source sigs? | 🟡 High |
| **2.5.5** **No OH trend / CH₄-OH feedback** | The manuscript acknowledges this (§4.5) but doesn't test it. Given He 2026b's finding of 25% bias from ignoring the feedback, this is a significant limitation. | 🟡 High |
| **2.5.6** **2-box transport limitations** | Naus et al. (2019) showed that 2-box models introduce systematic biases. The τ_ex sensitivity (§3.3.5) confirms this — σ(FF) ranges from 3.9 to 23.7 Tg/yr depending on τ_ex. A reviewer may argue the 2-box is too simple for the claims being made. | 🟡 High |

### 2.6 Presentation & Clarity

| Item | Question | Priority |
|------|----------|----------|
| **2.6.1** Abstract length | ~250 words — appropriate for ACP. But the abstract quotes specific numbers that differ from the manuscript body (see 2.4.1). | 🟡 High |
| **2.6.2** Figure quality | Are figures publication-ready (300 dpi, colorblind-safe, Nature-family formatting)? | 🟢 Medium |
| **2.6.3** Notation consistency | Are δ¹³C, δD, α, τ, σ used consistently throughout? | 🟢 Medium |
| **2.6.4** Reference completeness | 24 references cited — reasonable for ACP. Check: are all cited works in the reference list? Are any key omissions? (e.g., Sherwood et al. 2017 is cited for source sigs but the hemispheric disaggregation is novel — needs more explanation) | 🟢 Medium |

---

## 3. Prioritized Review Checklist

### Phase A — Blocking Issues (must fix before submission)

- [ ] **A1.** Resolve all number discrepancies between RESULTS.md and MANUSCRIPT (§2.4.1). Determine which numbers are correct and ensure a single consistent dataset underpins all documents.
- [ ] **A2.** Justify the W matrix (§2.5.1). Run sensitivity tests with W = I, W = diag(1, 1, 1), and inverse-variance-weighted W. Report impact on all headline numbers. If W matters (likely), discuss in the manuscript.
- [ ] **A3.** Verify the variance decomposition methodology (§2.1.2). Confirm that "selective freezing" correctly attributes variance when parameters interact nonlinearly. Consider an alternative (e.g., Sobol indices or ANOVA decomposition) as a robustness check.
- [ ] **A4.** Verify the δ¹³C atmospheric MC sampling logic (§2.1.4). The global offset approach may underestimate gradient uncertainty.
- [ ] **A5.** Verify hemispheric δ¹³C source signature build pipeline end-to-end (§2.3.1). This is the v3 change that flipped all results — it must be bulletproof.

### Phase B — High-Priority Improvements

- [ ] **B1.** Add BB sensitivity test (§2.5.2): ±20% perturbation of prescribed BB emissions.
- [ ] **B2.** Test time-varying source signatures (§2.5.4): allow FF δ¹³C to drift with changing coal/gas mix.
- [ ] **B3.** Increase MC iterations to 1000+ and verify convergence (§2.2.1). Plot σ(FF) vs. N_iter to show convergence.
- [ ] **B4.** Replace difference-of-means trend with linear regression + p-value (§2.2.2), or justify the current choice.
- [ ] **B5.** Test seed sensitivity (§2.2.4): run with 5 different seeds, report spread.
- [ ] **B6.** Clarify the EDGAR comparison discrepancy (§2.3.4): +4.8 vs +1.3 vs +8.8 Tg/yr appear in different documents for the same model output.

### Phase C — Strengthening for High-Impact Journal

- [ ] **C1.** Add a formal comparison with Dasgupta et al. (2025) and Fujita et al. (2025) — both recent dual/multi-isotope studies. The manuscript cites them but doesn't compare quantitatively.
- [ ] **C2.** Add a "what would it take?" section: what source-signature precision would eliminate the remaining ambiguity? Back-calculate the required σ(Mic δ¹³C), σ(FF δ¹³C), σ(Mic δD) for a robust trend.
- [ ] **C3.** Consider adding a 3-box (NHext/Trop/SHext) version as supplementary — the `three-box` branch name suggests this was planned. Even preliminary 3-box results would strengthen the paper.
- [ ] **C4.** Supplementary table of all station-level δD data used, with coverage years and hemispheric assignment.
- [ ] **C5.** Code review of `common.py` — the shared engine imported by all scripts. This is the single point of failure for the entire experiment.

### Phase D — Minor / Cosmetic

- [ ] **D1.** Unify notation: manuscript uses both "σ(Mic δD)" and "Mic δD uncertainty" — pick one.
- [ ] **D2.** Add DOIs to all references.
- [ ] **D3.** Check figure axes labels, units, and legends for consistency.
- [ ] **D4.** Manuscript §2.4 lists MC parameters — ensure all match what `core.py` actually samples.

---

## 4. Key Risks for Rejection

| Risk | Severity | Mitigation |
|------|----------|------------|
| **W matrix sensitivity** — if results change substantially with different W, the headline numbers are arbitrary | 🔴 Fatal | Run W sensitivity now; if results are W-dependent, either find a principled W or acknowledge and present the range |
| **Internal inconsistency** — reviewers who spot different numbers in abstract vs. tables vs. figures will lose trust | 🔴 Fatal | Single-source-of-truth: regenerate all tables/figures from one definitive model run |
| **"2-box is too simple" objection** — given that 3-D inversions exist (Basu, Thanwerdas), why trust a 2-box? | 🟡 Major | Emphasize the 2-box as a structural/information-theoretic analysis, not a competitor to 3-D inversions. The threshold and variance decomposition are about information content, not about producing the best emission estimate. |
| **v2→v3 narrative** — the paper's strongest finding (KIE immunity) was invalidated by the authors' own upgrade. Reviewers may ask: what happens in v4? | 🟡 Major | Frame honestly: v3 shows the result is sensitive to source-sig assumptions, which is itself a finding. Don't oversell robustness. |
| **No CH₄-OH feedback** — He 2026b showed 25% bias | 🟡 Major | Acknowledge explicitly; argue that the variance decomposition (relative contributions) is less sensitive than absolute magnitudes |

---

## 5. Suggested Review Timeline

| Week | Task | Output |
|------|------|--------|
| 1 | Resolve Phase A blocking issues (A1–A5) | Consistent dataset, W sensitivity results |
| 2 | Phase B improvements (B1–B6) | Additional sensitivity tests, convergence check |
| 3 | Phase C strengthening (C1–C5) | Enhanced discussion, supplementary material |
| 4 | Final manuscript revision, figure polish (D1–D4) | Submission-ready draft |

---

## 6. Bottom Line Assessment

**The experiment asks an important and timely question** — the OH-¹³C KIE controversy is the single largest identified source of disagreement in the methane isotope community, and quantifying whether δD can resolve it is directly relevant to multiple recent high-profile papers (Basu 2022, Riddell-Young 2025, Thanwerdas 2024).

**The analysis is extensive and technically competent** — 13 phases covering variance decomposition, multi-parameter sensitivity, bootstrap CIs, and validation against inventories. The code is clean, modular, and reproducible.

**However, three issues must be addressed before submission:**

1. **Internal number inconsistencies** between RESULTS.md and the manuscript suggest multiple model versions may be conflated. This is the most urgent fix.
2. **The W matrix** is unjustified and potentially controls the headline result (the balance between δ¹³C and δD information). Without a sensitivity test, the 45% improvement claim and the 25% KIE contribution are parameter choices, not findings.
3. **The v2→v3 narrative** — where the authors' own improvement invalidated their prior conclusion — is scientifically honest but needs careful framing. The paper should lead with "source signatures matter more than model complexity" rather than "KIE immunity," since the immunity claim was weakened.

If these are addressed, this is a strong candidate for ACP or GBC. The δD threshold result (reconciling Thanwerdas and Riddell-Young) is the clearest novel contribution and should be the paper's anchor.
