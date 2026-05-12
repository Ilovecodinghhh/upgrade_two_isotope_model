# PLAN.md — Next Steps & Prompts

**Last updated:** 2026-05-12  
**Status:** Three experiments complete (dD_threshold, KIE_sensitivity, KIE_immunity).  
Below: concrete next actions with copy-paste prompts.

---

## Priority 1: Strengthen the FF Trend Reversal Finding

The most provocative result is that real hemispheric δD reverses the post-2007 FF trend from +13 → −7 Tg/yr. This needs robustness testing before publication.

### 1a. Lifetime sensitivity on FF trend

> **Prompt:** In `experiments/KIE_immunity`, add a Phase 5 script that tests FF trend sensitivity to lifetime assumptions. Run the dual real-hemi-dD model at τ = 8.0, 8.5, 9.0, 9.5, 10.0 yr (fixed), and also with the time-varying τ from He 2026. For each, report the post-2007 ΔFF trend (median ± CI). The key question: does the negative FF trend survive across all reasonable lifetime assumptions, or does it flip sign at some τ? Save results to `results/phase5_tau_sensitivity.json` and update `RESULTS.md`.

### 1b. OH_D KIE sensitivity on FF trend

> **Prompt:** In `experiments/KIE_immunity`, add a Phase 6 script testing whether the FF trend reversal depends on the OH-D KIE. Run dual real-hemi-dD with OH_D = 1.294, 1.310 (midpoint), 1.327, and 1.35 (He 2026 upper). Report post-2007 ΔFF for each. If the trend is robust to OH_D, this strengthens the finding; if it flips, we need to report the OH_D threshold analogous to the σ(Mic δD) threshold from Exp 1.

### 1c. Cl fraction sensitivity

> **Prompt:** Test FF trend under Cl fractions of 0.6% (Thanwerdas), 3.5% (default), 6.5% (high), and 10% (Allan 2007 upper). Use dual real-hemi-dD. The Cl sink has the largest δD KIE (α=1.52), so its fraction could dominate the δD constraint. Report ΔFF trend for each Cl fraction.

---

## Priority 2: 3-Box Model Implementation

We have 3-box source signatures and atmospheric δD ready. Building the 3-box model tests whether finer spatial resolution changes conclusions.

### 2a. Write the 3-box model

> **Prompt:** Write `3x3_three.py` — a 3-box (NHext / Trop / SHext) version of the 3×3 simultaneous model. Use the 3-box δD source signatures (`{Mic,BB,FF}_dD_{NHext,Trop,SHext}_MC.csv`) and 3-box atmospheric δD (`ThreeBox_atm_dD_annual.csv`). Inter-box exchange: NHext↔Trop (τ_ex ~ 1.0 yr), Trop↔SHext (τ_ex ~ 1.0 yr). No direct NHext↔SHext exchange. Use bounded LS like `3x3_two.py`. Run 1000 MC iterations, save to `results/`. Compare FF total and trend with the 2-box version.

### 2b. Build 3-box atmospheric δD MC iterations

> **Prompt:** The 3-box atmospheric δD currently uses annual means only (no MC iterations). We need MC iterations. Modify the Riddell-Young `dD_globmean.py` approach: group stations by latitude band (NHext: >30°N, Trop: 30°S–30°N, SHext: <30°S), compute MBL surface + uncertainty + bootstrap resampling per band. Save as `NHext_dD_MC.xlsx`, `Trop_dD_MC.xlsx`, `SHext_dD_MC.xlsx` in `rel/data/`. Then update `common.py` to load these in 3-box mode.

---

## Priority 3: Publication Figures

### 3a. Master comparison figure

> **Prompt:** Create a publication-quality 3-panel figure for the KIE_immunity experiment: Panel A: variance decomposition bar chart (δ¹³C-only / dual offset / dual real hemi), stacked by KIE/sig/τ/residual. Panel B: FF time series with uncertainty bands for all 3 configs (smoothed, 2σ). Panel C: Basu 2022 comparison — our KIE spread vs theirs. Use matplotlib, export PNG+PDF to `experiments/KIE_immunity/figures/`. Style: Nature-family format, single column width (~88mm), 8pt font.

### 3b. Cross-experiment synthesis figure

> **Prompt:** Create a single figure that synthesizes all 3 experiments: Panel A: the δD threshold curve from Exp 1 (σ_Mic_dD vs improvement). Panel B: the agreement-rate discriminant from Exp 2 (threshold vs Cantrell−Saueressig pp). Panel C: the variance decomposition from Exp 3 (3 bars). Panel D: FF trend comparison (δ¹³C-only vs dual-offset vs dual-real-hemi). This is the "elevator pitch" figure for a potential cover letter or conference talk.

---

## Priority 4: Paper Draft Skeleton

### 4a. Structure the narrative

> **Prompt:** Write a paper outline (section headers + 2–3 sentence summaries) for a paper titled something like: "Hemispheric δD resolves the methane source attribution paradox: KIE immunity, trend reversal, and a 25‰ threshold". Structure: Abstract → Introduction (16-paper lit review, identify the 3 contradictions from THEMATIC_GAP_ANALYSIS.md) → Methods (2-box model, MC framework, hemispheric δD pipeline) → Results (Exp 1 → Exp 2 → Exp 3 in order of increasing sophistication) → Discussion (resolve the 3 contradictions, implications for FF trend) → Conclusion. Target: Nature Communications or ACP.

### 4b. Write the abstract

> **Prompt:** Write a 250-word abstract for the paper. Key numbers to include: σ(FF) reduction from 31.1 → 14.3 Tg/yr (54%), KIE contribution drops from 11% → 0%, the 25‰ threshold, the 25.4 pp agreement-rate discriminant, and the FF trend reversal (+13 → −7 Tg/yr). Tone: assertive but not overclaiming. Acknowledge the trend reversal needs further validation.

---

## Priority 5: Model Robustness / Extensions

### 5a. Bootstrap uncertainty on variance decomposition

> **Prompt:** The current variance decomposition uses n_iter=400. Add bootstrap confidence intervals: resample the 400 MC iterations 1000 times, compute σ(FF) for each bootstrap sample, report 95% CI on each variance component percentage. This tells us whether the "42.6% source signatures" is significantly different from "85% residual" in the offset version.

### 5b. Independent validation with EDGAR/GAINS inventories

> **Prompt:** Compare our model's FF time series (median + CI from dual real-hemi-dD) against EDGAR 8.0 and GAINS bottom-up inventories. Plot our top-down FF estimate alongside the bottom-up data. If our negative trend is consistent with EDGAR's reported leveling-off of fossil CH₄ after ~2010, that's strong validation. Load EDGAR data from `rel/data/CarbonTracker_CH4.xlsx` or download from EDGAR website.

### 5c. Sensitivity to interhemispheric exchange time

> **Prompt:** τ_ex is currently N(1.0, 0.1) yr. Test extreme values: 0.5 yr (fast mixing) and 2.0 yr (slow mixing). If τ_ex matters, the hemispheric constraint is genuinely providing new information. If results are insensitive, the 2-box advantage is from source-sig separation, not transport.

---

## Priority 6: Data Gaps to Fill

### 6a. Post-2019 hemispheric δD

The NH/SH atmospheric δD MC data has all-NaN for 2020–2023 (forward-filled from 2019). This is the biggest data gap. Options:
- Contact Riddell-Young group for updated iterations
- Use the semi-hemispheric annual means (which exist through 2023) as point estimates with inflated uncertainty
- Wait for next NOAA data release

### 6b. Tropospheric Cl fraction constraint

> **Prompt:** The Cl sink fraction (default 3.5%) is poorly constrained but has the largest δD KIE. Search the recent literature (2023–2026) for updated constraints on tropospheric Cl as a CH₄ sink. Key papers to check: Allan et al. 2007, Hossaini et al. 2016, Wang et al. 2021, and any newer estimates. Summarize the current range and whether our default 3.5% is still reasonable.

---

## Completed ✅

- [x] Experiment 1: δD threshold (Phases 1–5, 4-panel figure)
- [x] Experiment 2: KIE sensitivity (Phases 1–8, 14 figures)
- [x] Experiment 3: KIE immunity (v1 + v2 with real hemi δD)
- [x] Hemispheric δD source signature pipeline (2-box + 3-box)
- [x] `common.py` upgrade with real hemispheric δD loading
- [x] Thematic gap analysis of 16 reference papers
- [x] Full sensitivity test matrix documentation
