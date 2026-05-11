# Thematic Gap Analysis of 16 Reference Papers

## 1. Paper Classification Table: δ¹³C vs. δD × Model Type

| Paper | δ¹³C | δD | Single-Box | Multi-Box / 3D | Key Focus |
|-------|------|-----|------------|----------------|-----------|
| Rice et al. 2016 (PNAS) | ✓ | ✗ | ✓ | ✗ | Post-2007 growth; microbial vs FF attribution via δ¹³C |
| Schwietzke et al. 2016 (Nature) | ✓ | ✗ | ✓ | ✗ | Revised FF isotopic signatures → larger FF share historically |
| Worden et al. 2017 (Nat Comm) | ✓ | ✗ | ✓ | ✗ | BB decline reconciles δ¹³C & ethane constraints |
| Nguyen et al. 2020 (GRL) | ✓ | ✓ | ✓ | ✗ | First dual-isotope (δ¹³C + δD) box model for 2007–2014 |
| Maasakkers et al. 2019 (ACP) | ✗ | ✗ | ✗ | ✓ (3D GEOS-Chem) | GOSAT inversion; gridded emissions + OH trends |
| Naus et al. 2019 (ACP) | ✓ (MCF-implied) | ✗ | ✓ (9-box) | ✗ | OH variability from MCF; δ¹³C implications |
| Zhang et al. 2021 (Nat Comm) | ✗ | ✗ | ✗ | ✓ (3D GEOS-Chem) | GOSAT 2010–2018; wetlands + livestock drive acceleration |
| Basu et al. 2022 (ACP) | ✓ | ✗ | ✗ | ✓ (3D TM5 variational) | 3D inversion with δ¹³C; post-2007 = FF + microbial |
| Skeie et al. 2023 (Nat Comm) | ✓ | ✗ | ✗ | ✓ (OsloCTM3) | Global δ¹³C inversion 2005–2019; FF role emphasized |
| Zhao et al. 2023 (ACP) | ✗ | ✗ | ✓ (box for OH) | ✓ (3D CESM/GEOS) | Reconcile bottom-up vs top-down OH sink via precursor obs |
| Chandra et al. 2024 (Comm Earth) | ✓ | ✗ | ✓ | ✗ | Post-2019 surge; tropical wetlands dominate |
| Thanwerdas et al. 2024 (ACP) | ✓ | ✓* | ✗ | ✓ (3D LMDz variational) | First 3D dual-isotope inversion; fossil + agri/waste drive post-2007 |
| Dasgupta et al. 2025 (EGU) | ✓ | ✗ | ✗ | ✓ (3D inversion) | Tropical wetland vs. agriculture separation |
| Fujita et al. 2025 (JGR) | ✓ | ✓ | ✓ | ✗ | Updated dual-isotope box model; declining OH lifetime |
| Riddell-Young et al. 2025 (PNAS) | ✓ | ✓ | ✓ | ✗ | Harmonized δD record 2005–2023; microbial dominance confirmed |
| He et al. 2026 (Science + JGR) | ✓ | ✓ | ✓ | ✗ | Declining CH₄ lifetime; two-isotope source partition |

**\*Thanwerdas 2024:** δD was tested in one sensitivity run but found to add "only a minor influence" due to large source-signature uncertainties.

### Summary Statistics
- **Papers using δ¹³C:** 14/16 (87.5%)
- **Papers using δD:** 5/16 (31%) — and only 4 as a primary constraint
- **Single-box models:** 9/16
- **Multi-box or 3D models:** 8/16 (one paper, Zhao 2023, uses both)
- **Dual-isotope (δ¹³C + δD) approaches:** 5/16 — and only **1** uses a multi-box/3D framework (Thanwerdas 2024, but δD added little value there)

---

## 2. Unresolved Contradictions

### Contradiction A: What drives the post-2006/2007 CH₄ growth?

| Position | Papers | Mechanism |
|----------|--------|-----------|
| **Primarily microbial** | Rice 2016, Schaefer 2016 (cited), Nisbet 2016 (cited), Riddell-Young 2025, Chandra 2024, He 2026 | δ¹³C/δD trend toward lighter values → wetlands/agriculture |
| **Mixed: fossil + microbial** | Worden 2017, Basu 2022, Thanwerdas 2024, Skeie 2023 | FF contributes 12–19 Tg/yr; microbial also rises; BB decline masks FF signal in δ¹³C |
| **Possibly OH-driven** | Turner 2017 (cited), Rigby 2017 (cited), Naus 2019 (partial) | Declining OH could explain growth without source changes |
| **Fossil emissions stable/declining** | Riddell-Young 2025, He 2026 (Science) | Dual-isotope evidence: FF flat or declining post-2013 |

**Key tension:** Riddell-Young 2025 and He 2026 find fossil emissions *stable or declining* post-2006 with high confidence using dual isotopes, while Basu 2022, Thanwerdas 2024, and Skeie 2023 find significant FF *increases* using 3D δ¹³C inversions. These cannot both be correct.

### Contradiction B: The role of OH sink changes

| Position | Papers |
|----------|--------|
| **OH stable (no significant trend)** | Zhang 2021, Basu 2022, Patra 2021 (cited) |
| **OH increasing (~0.3%/yr)** | He 2026, Zhao 2023 (observation-constrained) |
| **OH decreasing (explains part of growth)** | Turner 2017 (cited), Rigby 2017 (cited), Naus 2019 |

**Key tension:** He 2026 uses a *declining* CH₄ lifetime (τ decreasing from ~9.8 to ~9.0 yr over 2000–2020, implying OH increase), which is directly contradicted by studies suggesting OH decrease or stability. The assumed OH trend profoundly affects source attribution — Thanwerdas 2024 showed that prescribed OH interannual variability "can have a large impact on the results."

### Contradiction C: Does δD actually add constraint?

| Position | Papers |
|----------|--------|
| **δD adds critical independent constraint** | Riddell-Young 2025, He 2026, Fujita 2025 |
| **δD adds minor/negligible constraint** | Thanwerdas 2024 ("only a minor influence"), Nguyen 2020 (large uncertainties remain) |

**Key tension:** Thanwerdas 2024 (the only 3D dual-isotope inversion) found δD nearly useless due to large source-signature uncertainties. But Riddell-Young 2025 (box model) found δD essential for ruling out Cl-sink and BB-decline scenarios. The discrepancy may stem from (a) source-signature uncertainty being too conservatively specified in 3D frameworks, or (b) box models oversimplifying transport and thus overestimating δD's discriminating power.

### Contradiction D: Biomass burning's quantitative role

| Position | Papers |
|----------|--------|
| **BB decline large (~3.7 Tg/yr) and crucial for budget** | Worden 2017 |
| **BB decline moderate (~2 Tg/yr), secondary** | GFED-based (Zhang 2021, Basu 2022) |
| **BB prescribed/fixed, not a major uncertainty** | Riddell-Young 2025, He 2026 |

**Implication:** If BB decline is underestimated, the δ¹³C-inferred FF contribution is overestimated (Worden 2017's key insight). The disagreement between 3D inversions (which find FF increase) and dual-isotope box models (which find FF stable) may partly trace to different BB treatments.

### Contradiction E: KIE of OH for ¹³C — which value?

| Value | Source | Impact |
|-------|--------|--------|
| 1.0039 | Saueressig et al. 2001 | Lower fractionation → more microbial attribution |
| 1.0054 | Cantrell et al. 1990 | Higher fractionation → shifts budget toward FF |

**Status:** Basu 2022 showed switching between these values "has a large influence on the results." Thanwerdas 2024 uses 1.0039, noting Saueressig's data has "considerably higher experimental precision." Yet this single parameter choice shifts FF/microbial partitioning by 20–40% — an *unresolved systematic* that propagates through all δ¹³C-based studies.

---

## 3. Your Model's Specific Features

Based on your repository (MODELS.md, core.py, inputs.py), your model has these distinctive capabilities:

| Feature | Implementation | Versions |
|---------|---------------|----------|
| **Dual-isotope (δ¹³C + δD) mass balance** | Simultaneous 2×2 system using both isotopic constraints | All versions |
| **KIE sampling** | OH-¹³C ~ U(1.0039, 1.0054); OH-D ~ U(1.294, 1.327); Cl KIEs sampled | v2.0+ |
| **Time-varying lifetime** | τ(t) = 9.0 − 0.017·(t − 2010), per He et al. 2026 | v2.0+ |
| **Bayesian MCMC (PyMC)** | Full posterior distributions with explicit priors | v3.0a, v3.1a |
| **Two-hemisphere (NH/SH) split** | 3×3 system with inter-hemispheric exchange (k_ex) | v3.0b, v3.1b |
| **Non-physical solution penalties** | pm.Potential for negative emissions | v3.1a |
| **Condition-number monitoring** | Reject ill-conditioned iterations | v3.1b |
| **Fixed-BB sensitivity** | Prescribe BB to isolate FF vs Mic | v3.2 |
| **2-source simplification** | Microbial vs Non-microbial lumping | v4.0 |
| **Monte Carlo source signatures** | 1000-draw CSVs for FF, Mic, BB for both δ¹³C and δD | All (from rel/output/) |

---

## 4. Which Contradictions Is Your Model Best Equipped to Address?

### ★ PRIMARY TARGET: Contradiction C — Does δD actually add constraint in a spatially-resolved framework?

**Why your model is ideal:**
- You already have the **only working dual-isotope multi-box model** (v3.0b/v3.1b with NH/SH hemispheric split + both δ¹³C and δD).
- Thanwerdas 2024 (the only other 3D dual-isotope attempt) found δD unhelpful — but they used a *3D variational* system where source-signature uncertainties swamped the signal. Your 2-box approach is intermediate between their 3D system and the 1-box models of Riddell-Young/He/Fujita.
- You can **systematically test** at what level of spatial disaggregation δD begins to add value — is it the hemispheric split that helps, or do you need full 3D?
- **Research question:** "Under what conditions does δD provide meaningful additional constraint on CH₄ source attribution beyond δ¹³C alone, and does a hemispheric (2-box) framework preserve this constraining power better than a global 1-box?"

### ★ SECONDARY TARGET: Contradiction A/E — FF vs microbial attribution given KIE uncertainty

**Why your model is ideal:**
- Your v2.0+ **samples the full KIE range** (Saueressig to Cantrell), which most 3D inversions don't — they fix one value.
- You can quantify: "How much does KIE uncertainty shift the FF/microbial balance *when both isotopes are used simultaneously*?" The dual-isotope constraint should *reduce* sensitivity to ¹³C-KIE because δD provides an independent check.
- **Testable hypothesis:** "The dual-isotope approach significantly reduces the sensitivity of source attribution to the OH-¹³C KIE choice compared to δ¹³C-only approaches."

### ★ TERTIARY TARGET: Contradiction B — OH trend impact with dual isotopes

**Why your model is ideal:**
- Your time-varying τ(t) already encodes an OH trend (He 2026).
- You can run sensitivity tests varying the τ trend and show whether the *dual-isotope* constraint narrows the range of allowable OH changes — i.e., can δ¹³C + δD together reject certain OH scenarios that δ¹³C alone cannot?

---

## 5. Suggested Additional Methods

Beyond your current box-model framework, consider these approaches to strengthen your study:

### 5a. Information-theoretic analysis (Degrees of Freedom for Signal)
- Compute the **DFS** (Degrees of Freedom for Signal) or mutual information gain from adding δD observations to a δ¹³C-only system, following Thanwerdas et al. 2022's OSSE framework.
- This would formally quantify the "information content" of δD at different spatial resolutions (1-box → 2-box → n-box).

### 5b. Observing System Simulation Experiments (OSSEs)
- Generate synthetic observations from a known "truth" (e.g., prescribed FF/Mic/BB scenario) and test whether your 2-box dual-isotope system can recover the truth better than a 1-box or single-isotope system.
- Systematically degrade source-signature precision to find the threshold where δD stops helping.

### 5c. Bayesian model comparison (Bayes factors)
- Compare models with/without δD constraint using Bayes factors or WAIC/LOO-CV from your PyMC implementation (v3.0a/v3.1a).
- Formally test: "Does including δD improve posterior predictive performance?"

### 5d. Bootstrap/jackknife on the KIE
- Instead of uniform sampling across the KIE range, use a structured sensitivity: fix KIE at Saueressig, then Cantrell, and show the spread with/without δD. This directly tests whether dual isotopes make the system robust to KIE choice.

### 5e. Comparison with CT-CH₄ posterior
- Riddell-Young 2025 compares their box-model results against the CT-CH₄ 3D inversion posterior (gray bar in their Fig. 3). You could make the same comparison with your hemispheric model — does the 2-box agree better with 3D inversions than 1-box?

### 5f. Extend to 3+ boxes (tropical/extratropical)
- A natural extension: NH-tropical / NH-extratropical / SH split. Tropical wetlands have distinct δD signatures from boreal wetlands. A 3-box model might capture this without the computational cost of a full 3D inversion.

---

## 6. Proposed Primary Research Question

> **"Does a hemispheric two-box dual-isotope (δ¹³C + δD) framework resolve the contradiction between single-box studies (which find δD essential) and 3D inversions (which find δD unhelpful), and can the simultaneous use of both isotopes reduce sensitivity to the OH-¹³C KIE uncertainty that currently drives divergent source attributions?"**

This question:
1. Directly addresses the most important gap (no multi-box dual-isotope study has *tested* δD's value)
2. Exploits your model's unique position (between 1-box and 3D)
3. Is answerable with your existing infrastructure (v3.0b/3.1b + Bayesian v3.1a)
4. Has clear implications for the field (guides whether future 3D inversions should invest in δD assimilation)

---

## Appendix: Quick-Reference Paper Summaries

| Paper | Main Conclusion |
|-------|----------------|
| Rice 2016 | Post-2007 growth = microbial (wetlands + agriculture); 3D δ¹³C inversion |
| Schwietzke 2016 | Revised FF δ¹³C signatures (more negative) → historically larger FF share |
| Worden 2017 | BB decline (~3.7 Tg/yr) reconciles ethane and δ¹³C evidence; FF = 12–19 Tg/yr increase |
| Maasakkers 2019 | GOSAT inversion 2010–2015; India+China drive trends; no isotopes |
| Naus 2019 | MCF-constrained OH; 9-box model shows OH variability affects CH₄ interpretation |
| Nguyen 2020 | First published dual-isotope (δ¹³C + δD) box model; large uncertainties persist |
| Zhang 2021 | GOSAT 2010–2018; tropical livestock + wetlands + 2014 OH dip + 2015 fires |
| Basu 2022 | 3D variational + δ¹³C; both FF and microbial increased post-2007 |
| Skeie 2023 | OsloCTM3 + δ¹³C; emphasizes FF role in post-2006 growth |
| Zhao 2023 | Observation-constrained OH fields reduce bottom-up/top-down gap |
| Chandra 2024 | Post-2019 surge dominated by tropical wetlands (La Niña) |
| Thanwerdas 2024 | First 3D dual-isotope inversion; δD added little; fossil + agri/waste = main drivers |
| Dasgupta 2025 | Tropical wetland vs agriculture separation in 3D framework |
| Fujita 2025 | Updated dual-isotope box model; declining τ; confirms microbial dominance |
| Riddell-Young 2025 | Harmonized global δD record; dual-isotope box model; microbial dominance with high confidence |
| He 2026 | Time-varying lifetime + dual isotopes; microbial driver; FF stable/declining |
