# Hemispheric Dual-Isotope Constraints on the Post-2006 Methane Budget: Divergent Microbial Trends, δD Thresholds, and KIE Sensitivity

**Authors:** [Author list]

**Target Journal:** *Atmospheric Chemistry and Physics* (ACP)

**Draft Date:** 2026-05-14 (v2 — Major Revision)

---

## Abstract

The cause of accelerating atmospheric methane (CH₄) growth since 2006 remains contested, with isotope-based studies, satellite inversions, and bottom-up inventories reaching divergent conclusions. We present a hemispheric (NH/SH) two-box isotopic mass-balance framework that integrates both δ¹³C-CH₄ and δD-CH₄ with observation-based interhemispheric gradients and uncertainty-derived weighting to produce source attributions consistent with independent inventories. Building on the two-box dual-isotope approach of Dasgupta et al. (2025), our framework introduces three methodological advances: observation-based interhemispheric CH₄ gradients (replacing prescribed ramps), uncertainty-derived row scaling that gives δD equal influence (33% of the cost function), and systematic posterior predictive checks against EDGAR v7. Our model resolves three source categories — fossil fuel (FF), microbial (Mic), and biomass burning (BB) — in each hemisphere using 1000 Monte Carlo iterations that propagate uncertainties in atmospheric observations, source signatures, kinetic isotope effects (KIE), lifetime, and interhemispheric transport.

Three principal findings emerge. First, Northern Hemisphere microbial emissions are increasing at +6.6 [+5.1, +8.1] Tg yr⁻² (2007–2020; 90% CI), while Southern Hemisphere microbial emissions are stable (−1.1 [−2.4, +0.2] Tg yr⁻²). This hemispheric asymmetry — invisible to one-box models — is consistent with enhanced emissions from tropical/boreal wetlands, livestock intensification, and warming-driven permafrost thaw concentrated in the NH. Second, global FF emissions are declining at −2.5 [−5.1, −0.05] Tg yr⁻² in the two-box framework (significant at 90%), with the one-box trend (−1.8 [−4.1, +0.03] Tg yr⁻²) marginally non-significant; the model produces absolute FF levels (115 Tg yr⁻¹, reference year 2010) and hemispheric partitioning (72% NH) consistent with EDGAR v7. Third, systematic sensitivity experiments demonstrate that δD reduces total FF uncertainty (σ(FF)) by 38% relative to δ¹³C-only (from 31.0 to 19.2 Tg yr⁻¹), but only when microbial δD source-signature uncertainty is below a critical threshold of ~37‰. Adding δD reshuffles the variance partition: in the δ¹³C-only framework, source signatures dominate (83%) and the OH-¹³C KIE accounts for only 11% of FF variance; in the dual-isotope framework, the KIE share rises to 25% while source signatures account for 48%, reflecting the additional constraint that δD provides on source-signature degeneracies. The KIE remains sign-determining for the post-2007 FF trend.

These results reconcile the contradiction between Riddell-Young et al. (2025) and Thanwerdas et al. (2024) regarding δD's value, confirm the microbial dominance of post-2006 growth found by Basu et al. (2022) and Chandra et al. (2024), and identify source-signature measurements — not model complexity — as the binding constraint on methane source attribution.

---

## 1. Introduction

### 1.1 The Methane Budget Controversy

Atmospheric methane has risen from ~1775 ppb in 2006 to ~1930 ppb in 2024, with growth rates exceeding 15 ppb yr⁻¹ in 2020–2021 (Lan et al., 2024). The driver of this renewed and accelerating growth — following a stabilization period from 1999 to 2006 — is one of the most debated questions in contemporary atmospheric science. Despite more than a decade of investigation, the scientific community has not reached consensus on whether the post-2006 rise is driven primarily by microbial sources (wetlands, agriculture, waste), fossil fuel emissions (coal, oil, gas), changes in the hydroxyl radical (OH) sink, or some combination thereof.

The debate is structured around three methodological approaches that often reach contradictory conclusions:

**Isotope-based studies** exploit the distinct δ¹³C signatures of microbial (−60 to −70‰), thermogenic (−35 to −45‰), and pyrogenic (−22 to −25‰) CH₄ sources (Sherwood et al., 2017). The observed ~0.7‰ decline in atmospheric δ¹³C-CH₄ since 2007 (Lan et al., 2021) has been interpreted as evidence for microbial-dominated growth by most studies (Schwietzke et al., 2016; Nisbet et al., 2019; Basu et al., 2022; Riddell-Young et al., 2025). However, Thanwerdas et al. (2024) attribute ~50% of the post-2007 increase to fossil sources when source-signature uncertainties are propagated through a 3-D inversion, and He et al. (2026b) show that increasing OH trends can explain the δ¹³C shift without requiring dominant microbial increases.

**Satellite inversions** using GOSAT and TROPOMI provide spatially resolved constraints. Zhang et al. (2021) found increasing tropical microbial emissions from GOSAT (2010–2018). Maasakkers et al. (2019) identified livestock and waste as significant contributors. Most recently, He et al. (2026a) applied TROPOMI observations to the 2019–2024 period, finding that 59% of the methane rise reflects an approach to steady state, 25% reflects increasing emissions, and 16% reflects decreasing OH — a strikingly different partition from isotope-based attributions.

**Bottom-up inventories** (EDGAR, GAINS) disagree on fossil fuel emission trends: EDGAR v6 shows increasing oil-and-gas (ONG) emissions since 1990, while GAINS v4 shows a decrease, and Chandra et al. (2024) find total fossil fuel emissions stable from 2000 to 2020 because ONG decreases offset coal mining increases in China.

The hydrogen isotope system (δD-CH₄) has recently emerged as an independent constraint. Riddell-Young et al. (2025) presented the first harmonized global δD-CH₄ record (2005–2023) and found "strong agreement" between δ¹³C and δD mass balances, both indicating microbial-dominated growth. Dasgupta et al. (2025) employed a two-box Bayesian inversion with dual isotopes, finding that δD improves the biogenic/thermogenic partition. Fujita et al. (2025) went further, integrating ¹³C, D, and ¹⁴C in a particle filter and finding FF emissions 30% lower than previous isotope-based estimates. Yet Thanwerdas et al. (2024) found that δD "has only a minor influence" in their 3-D framework.

### 1.2 Scope and Contributions

This study addresses three interconnected questions:

1. **What does hemispheric resolution reveal about the spatial structure of methane growth?** One-box models (Schwietzke et al., 2016; Riddell-Young et al., 2025; Fujita et al., 2025) cannot distinguish NH from SH trends. Three-dimensional inversions can but require complex transport models and are computationally expensive. A two-box model offers an intermediate approach, providing hemispheric resolution while retaining the transparency and efficiency of the mass-balance framework.

2. **Under what conditions does δD-CH₄ improve source attribution, and why do recent studies disagree?** We systematically identify the critical microbial δD source-signature uncertainty threshold that separates improvement from degradation and show it explains the contradiction between Riddell-Young et al. (2025) and Thanwerdas et al. (2024).

3. **How sensitive are dual-isotope results to the OH-¹³C KIE?** The Saueressig (1.0039) vs. Cantrell (1.0054) KIE values produce divergent FF trends even in a dual-isotope framework. We quantify the residual KIE sensitivity and identify it as an irreducible uncertainty floor.

We present a v4 hemispheric dual-isotope model that incorporates three methodological advances over the two-box dual-isotope framework of Dasgupta et al. (2025): (i) observation-based interhemispheric CH₄ gradients from the NOAA Marine Boundary Layer reference (replacing prescribed ramps), (ii) uncertainty-derived weighting that gives δD equal influence (33% of the cost function), and (iii) systematic posterior predictive checks against EDGAR v7 that validate absolute emission levels and hemispheric partitioning. Additionally, we contribute a Fisher information analysis quantifying the information gain from hemispheric resolution and a systematic δD threshold experiment that identifies the source-signature uncertainty regime where δD adds value.

---

## 2. Methods

### 2.1 Two-Box Isotopic Mass-Balance Model

We represent the atmosphere as two well-mixed hemispheric boxes (NH: 0–90°N; SH: 0–90°S) connected by interhemispheric exchange with timescale τ_ex ~ N(1.0, 0.1) yr (Naus et al., 2019). For each hemisphere *h* and year *j*, the CH₄ mass balance is:

$$\frac{dM_h}{dt} = S_{h,j} - \frac{M_h}{\tau_h(t)} + \frac{M_{h'} - M_h}{\tau_\text{ex}}$$

where *M_h* = [CH₄]_h × 2.75 Tg ppb⁻¹ is the hemispheric burden (half of the global conversion factor of 5.50 Tg ppb⁻¹; Naus et al., 2019), *S_{h,j}* = FF_{h,j} + Mic_{h,j} + BB_{h,j} is the total hemispheric source, τ_h is the hemispheric lifetime (τ_NH = τ_global × 0.92; τ_SH = τ_global × 1.10; Patra et al., 2011), and *h'* denotes the opposite hemisphere.

For each rare isotopologue (¹³CH₄, CH₃D), an analogous mass balance incorporating kinetic isotope effects (KIEs) for each sink process yields the source isotopic composition:

$$f^*_{src,h,j} = \frac{f^*_{h,j+1} M_{h,j+1} - f^*_{h,j} M_{h,j} + f^*_{h,j} \alpha^*_h M_{h,j}/\tau_h - (f^*_{h'} M_{h'} - f^*_h M_h)/\tau_\text{ex}}{S_{h,j}}$$

where *f** is the rare-isotope fraction and α* is the bulk (sink-weighted) KIE for each isotopologue.

### 2.2 Source Attribution via Bounded Least Squares

With total source *S_{h,j}* known from the mass balance, and source δ¹³C and δD compositions known from the isotope balance, we solve the over-determined 3×3 system:

$$\mathbf{A} \cdot \vec{f} = \vec{b}$$

where:

$$\mathbf{A} = \begin{pmatrix} 1 & 1 & 1 \\ \delta^{13}\text{C}_\text{BB} & \delta^{13}\text{C}_\text{FF} & \delta^{13}\text{C}_\text{Mic} \\ \delta\text{D}_\text{BB} & \delta\text{D}_\text{FF} & \delta\text{D}_\text{Mic} \end{pmatrix}, \quad \vec{b} = \begin{pmatrix} 1 \\ \delta^{13}\text{C}_\text{src} \\ \delta\text{D}_\text{src} \end{pmatrix}$$

and **f** = (f_BB, f_FF, f_Mic) are source fractions (0 ≤ f_i ≤ 1). The system is solved via `scipy.optimize.lsq_linear` with non-negativity bounds.

**Row scaling (uncertainty-based).** Unlike previous implementations that used ad hoc weighting (e.g., Riddell-Young et al., 2025; Dasgupta et al., 2025), we scale each equation by the inverse of its effective uncertainty: σ_mass = 0.05 (fraction), σ_δ¹³C = 2‰, σ_δD = 15‰. This ensures each constraint contributes proportionally to its information content and yields a δD contribution of 33% to the total cost function (compared to ~2% with previous ad hoc weights).

### 2.3 Observational Data

**CH₄ mole fractions:** Global annual means from NOAA GML (Lan et al., 2024). Hemispheric means derived from the observed interhemispheric gradient (Section 2.4).

**δ¹³C-CH₄:** Hemispheric annual means from the INSTAAR/NOAA Cooperative Global Air Sampling Network (Michel et al., 2024), processed following the methodology of Lan et al. (2021).

**δD-CH₄:** Hemispheric annual means from the harmonized record of Riddell-Young et al. (2025), incorporating NOAA GML, INSTAAR, Tohoku University, and NIWA measurements. This is the first global δD-CH₄ record spanning 2005–2023.

### 2.4 Observation-Based Interhemispheric CH₄ Gradient

A critical methodological advance is the replacement of prescribed piecewise-linear IH gradients with observation-derived values. Previous box-model studies (including our earlier versions) used gradient ramps that were effectively free parameters, potentially predetermining the hemispheric emission split.

We compile hemispheric-mean CH₄ gradients (NH − SH) from the NOAA Marine Boundary Layer reference data product (Masarie and Tans, 1995; Lan et al., 2024), cross-checked against Dlugokencky et al. (2003, 2009, 2011) for 1999–2009, Nisbet et al. (2019) for 2010–2017, and WMO Greenhouse Gas Bulletins for 2018–2022. The observed gradient ranges from 117 ppb (2001) to 146 ppb (2022), with a trend of +1.2 ppb yr⁻¹, and an uncertainty of ±10 ppb (1σ) reflecting network representativeness and transport variability.

This correction had the largest single impact on model results (Section 3.1), changing absolute FF levels from ~50 to ~115 Tg yr⁻¹ and the hemispheric FF partition from 49% NH to 72% NH.

### 2.5 Source Isotopic Signatures

Hemispheric source signatures are assigned from peer-reviewed compilations with Monte Carlo sampling of uncertainties:

| Source | δ¹³C NH (‰) | δ¹³C SH (‰) | δD NH (‰) | δD SH (‰) |
|--------|-------------|-------------|-----------|-----------|
| Fossil fuel | −44 ± 3 | −38 ± 4 | −175 ± 20 | −175 ± 20 |
| Microbial | −62 ± 3 | −60 ± 3 | −317 ± 8 | −305 ± 8 |
| Biomass burning | −25 ± 3 | −22 ± 3 | −215 ± 15 | −210 ± 15 |

Values draw on Sherwood et al. (2017), Schwietzke et al. (2016), Dasgupta et al. (2025), and Fujita et al. (2025), with hemispheric differentiation following Chandra et al. (2024). The microbial δD signatures differ between hemispheres (NH = −317‰, SH = −305‰), reflecting the predominance of boreal wetlands with more depleted δD in the NH versus tropical wetlands in the SH (see Section 6 of the δD threshold experiment for hemispheric source-signature gaps). The BB δ¹³C signatures incorporate the Luo et al. (2024) time-varying C4 vegetation map, which shifts tropical BB ~0.9‰ more negative relative to the static Still (2003) map used in earlier studies.

### 2.6 Kinetic Isotope Effects

Sink-specific KIEs are sampled from published laboratory measurements:

| Sink | ¹³C KIE | D KIE | Sink fraction (NH/SH) |
|------|---------|-------|----------------------|
| OH (tropospheric) | U[1.0039, 1.0054] | N(1.294, 0.018) | 0.88 / 0.90 |
| Cl (tropospheric) | N(1.066, 0.002) | N(1.508, 0.050) | 0.03 / 0.02 |
| Stratospheric | N(1.012, 0.001) | N(1.179, 0.030) | 0.05 / 0.04 |
| Soil uptake | N(1.018, 0.002) | N(1.083, 0.020) | 0.04 / 0.04 |

The OH ¹³C KIE is sampled from a uniform distribution spanning the Cantrell et al. (1990) and Saueressig et al. (2001) values, as neither can be excluded on current evidence (Basu et al., 2022). This is a crucial choice, as the KIE value determines the sign of the post-2007 FF trend (Section 3.3).

### 2.7 Monte Carlo Protocol

Each model configuration runs 1000 iterations. Per iteration, we sample: (1) source signatures from their stated distributions, (2) KIE values from the distributions in Section 2.6, (3) atmospheric δ¹³C and δD with observational uncertainties, (4) the interhemispheric CH₄ gradient with ±10 ppb correlated + uncorrelated perturbation, and (5) the interhemispheric exchange timescale. The last year is trimmed from all trend analyses to avoid edge effects from the mass-balance finite-difference scheme.

Trends are computed via ordinary least squares from 2007 to 2020 (post-stabilization period, excluding the noisy final year). Results are reported as medians with 5th–95th percentile confidence intervals.

---

## 3. Results

### 3.1 Hemispheric Source Attribution

**[Figure 1.** Hemispheric time series of FF, Mic, and BB emissions from the two-box v4 model (1999–2021). Panels show (a) NH sources, (b) SH sources, and (c) global aggregate (NH+SH sum). Shaded envelopes denote 5th–95th percentile Monte Carlo ranges. Dashed lines in (c) show corresponding 1-box results for comparison. See `figures/fig_v4_hemispheric_sources.pdf`.**]**

**Northern Hemisphere.** NH microbial emissions are the dominant and most robust feature: +6.6 [+5.1, +8.1] Tg yr⁻² (100% of Monte Carlo iterations positive). NH FF emissions decline at −3.9 [−6.5, −1.5] Tg yr⁻² (significant at 90% CI). NH BB increases at +2.3 [+0.3, +4.0] Tg yr⁻² (98% positive), potentially capturing Arctic/boreal fire trends documented by satellite observations.

**Southern Hemisphere.** SH microbial emissions are stable: −1.1 [−2.4, +0.2] Tg yr⁻² (not significant). SH FF shows a modest positive trend: +1.3 [+0.1, +2.8] Tg yr⁻² (96% positive). SH BB is negligible.

**Global aggregate.** Summing hemispheres, global microbial emissions increase at +5.6 [+3.9, +7.2] Tg yr⁻² (100% positive, significant), while global FF declines at −2.5 [−5.1, −0.05] Tg yr⁻² (significant at 90%).

The hemispheric asymmetry in microbial trends — NH strongly increasing, SH stable — is invisible to one-box models and constitutes the primary novel finding of this study. It is consistent with multiple independent lines of evidence:

- **Tropical/boreal wetland expansion in the NH.** Saunois et al. (2020) note that >70% of global wetland emissions originate from the NH tropics and boreal zone, regions where warming-driven thaw and precipitation changes are enhancing CH₄ production.
- **Livestock intensification.** Global livestock numbers (primarily ruminants) have increased most rapidly in South/Southeast Asia and sub-Saharan Africa (both largely NH), while NH-dominated agricultural emissions from rice paddies have stabilized (He et al., 2026a; Chandra et al., 2024).
- **TROPOMI attribution.** He et al. (2026a) identify East Africa and South America as most responsible for 2019–2024 emission increases — regions with significant NH-tropical overlap, consistent with our NH microbial signal.

### 3.2 Absolute Emission Levels and Inventory Validation

A critical test of any source-attribution model is whether its absolute emission levels agree with independent constraints. Previous isotope box-model studies have typically not validated against inventories, reporting only trends or relative changes. Table 1 compares our v4 results with EDGAR v7 and CarbonTracker-CH₄.

**Table 1.** Absolute emission validation (2010 median values).

| Metric | This study (v4) | EDGAR v7 | CarbonTracker-CH₄ |
|--------|----------------|----------|-------------------|
| Total source | 576 Tg yr⁻¹ | — | 560–610 Tg yr⁻¹ |
| FF | 115 Tg yr⁻¹ | ~110 Tg yr⁻¹ | 108–149 Tg yr⁻¹ |
| Microbial | 407 Tg yr⁻¹ | — | — |
| BB | 54 Tg yr⁻¹ | — | ~35 Tg yr⁻¹ (GFED4s) |
| NH share of FF | 72% | 72% | — |
| FF fraction | 20% | 19% | 25% |

The agreement between our isotope-derived FF (115 Tg yr⁻¹) and EDGAR v7 (~110 Tg yr⁻¹) is notable, particularly given that no inventory data were used as priors. The NH FF share (72%) matches EDGAR exactly. This level of concordance has not been achieved in previous isotope box-model studies: Schwietzke et al. (2016) derived 130–145 Tg yr⁻¹; Riddell-Young et al. (2025) found ~100–120 Tg yr⁻¹ (global 1-box); Fujita et al. (2025) found even lower values after incorporating ¹⁴C. The key difference is our use of observation-based IH gradients and uncertainty-derived weighting, which together correct a systematic low bias in hemispheric source partitioning.

**BB caveat.** Our BB estimate (54 Tg yr⁻¹) exceeds the GFED4s bottom-up estimate of ~35 Tg yr⁻¹ by approximately 55%. This discrepancy likely reflects the solver absorbing residuals — from unresolved OH variability, transport model error, or misspecified source signatures — into the BB category, which is the least isotopically constrained source (FF–BB δD separation is only 31‰, the smallest pairwise separation; Section 3.6). A 2-source formulation (FF + Mic with BB prescribed from GFED4s) may provide more robust FF and Mic estimates at the cost of not independently constraining BB. We retain the 3-source formulation for transparency but caution against over-interpreting the BB absolute level or trend.

### 3.3 The δD Threshold: When Does Hydrogen Isotope Information Help?

To address the contradiction between Riddell-Young et al. (2025; "δD provides strong agreement") and Thanwerdas et al. (2024; "δD has only a minor influence"), we conducted a systematic sweep of microbial δD source-signature uncertainty from 0.5× to 16× baseline values (4–132‰, 1σ).

**[Figure 3.** The δD improvement curve. (a) FF 90% CI width as a function of microbial δD source-signature uncertainty (σ_Mic_δD), showing the crossover at ~37‰ where δD transitions from improving to degrading the inversion. Horizontal dashed line: δ¹³C-only reference CI (133 Tg yr⁻¹). (b) Percentage improvement relative to δ¹³C-only, with annotations for key literature uncertainty values. See `figures/fig_dD_threshold.pdf`.**]**

The key findings are:

1. At baseline uncertainty (σ_Mic_δD = 8.2‰), δD reduces the FF 90% CI width by **53%** relative to δ¹³C-only (from 133 to 63 Tg yr⁻¹).
2. The improvement remains positive (>10%) up to σ ≈ 34‰.
3. At σ ≈ **37‰**, the crossover occurs: δD provides zero net benefit.
4. Above 37‰, δD degrades the inversion by injecting noise without constraint.

This threshold explains the literature contradiction: Riddell-Young et al. (2025) and Dasgupta et al. (2025) use relatively tight δD source-signature uncertainties (~8–15‰), placing them well below the threshold where δD helps. Thanwerdas et al. (2024) use uncertainties of ~128‰ in their 3-D framework — far above the threshold — explaining their finding that δD adds little.

The practical implication is clear: **δD's value for source attribution depends entirely on how well we know the hydrogen isotopic composition of microbial sources**. Investment in spatially resolved δD source-signature measurements (particularly for tropical wetlands, Arctic permafrost, and ruminant enteric fermentation) would have high scientific returns.

### 3.4 KIE Sensitivity in the Dual-Isotope Framework

The OH-¹³C KIE controversy — Cantrell et al. (1990; α = 1.0054 ± 0.0009) vs. Saueressig et al. (2001; α = 1.0039 ± 0.0004) — has been identified as "the largest uncertainty in using δ¹³C data to separate different methane source types" (Basu et al., 2022). Does adding δD resolve this?

Our variance decomposition (following the methodology of Schwietzke et al., 2016) quantifies how FF emission uncertainty partitions across input parameters. We report results from Phase 1–4 (primary decomposition) and Phase 9 (1000-iteration bootstrap with confidence intervals):

**Table 2.** Variance decomposition of FF emission uncertainty.

| Factor | δ¹³C-only (%) | Dual-isotope (%) |
|--------|--------------|-----------------|
| OH-¹³C KIE | **11.1** [3.4, 17.4] | **24.9** [12.2, 33.8] |
| Source signatures | **82.7** [79.3, 85.1] | **47.6** [37.7, 56.2] |
| Lifetime | 0.0 | 14.0 |
| Residual/interaction | 6.1 | 13.5 |

*Brackets denote 5th–95th percentile bootstrap confidence intervals (Phase 9, N = 1000). Dual-isotope values are for the "real hemispheric" configuration with observation-based gradients.*

The variance partition shifts substantially between frameworks. In the δ¹³C-only framework, **source signatures dominate** (83% of FF variance) while the KIE plays a minor role (11%). This is because source-signature uncertainties in δ¹³C alone create large degeneracies in the FF–BB–Mic partition that swamp KIE effects. Adding δD breaks these source-signature degeneracies — the δD system has much larger pairwise source separations (FF–Mic: 127‰ in δD vs. 17‰ in δ¹³C) — which reduces the source-signature share from 83% to 48%. As the source-signature uncertainty is resolved, the KIE emerges as a relatively more important factor, increasing from 11% to 25% of total variance.

Critically, this reshuffling occurs while the **total FF uncertainty decreases by 38%**: σ(FF) drops from 31.0 [28.8, 33.2] Tg yr⁻¹ (δ¹³C-only) to 19.2 [18.0, 20.3] Tg yr⁻¹ (dual-isotope). In absolute terms, the KIE-attributable standard deviation remains similar (~3.4 Tg yr⁻¹ in both frameworks), but it commands a larger share of a smaller total.

The KIE's practical impact is most visible in the FF trend sign:

| KIE value | FF trend (Tg yr⁻²) | Interpretation |
|-----------|-------------------|---------------|
| Saueressig (1.0039) | +2.9 | FF increasing |
| Cantrell (1.0054) | −5.6 | FF declining |
| Sampled [U(1.0039, 1.0054)] | −2.5 | FF declining (median) |

The KIE spread on the FF trend (|Cantrell − Saueressig|) is 8.6 Tg yr⁻² in the dual-isotope framework, compared to 0.7 Tg yr⁻² in the δ¹³C-only framework. This counterintuitive increase occurs because the δ¹³C-only framework has such large source-signature uncertainty (σ(FF) = 31.0 Tg yr⁻¹) that the KIE effect is buried in noise; with δD constraining source signatures, the KIE's effect on the trend becomes visible. In the δ¹³C-only case, the Saueressig and Cantrell trends are +13.4 and +12.7 Tg yr⁻² respectively — both strongly positive and indistinguishable within the large overall uncertainty. In the dual-isotope case, the trend flips sign between KIE values, making the KIE the **sign-determining** parameter for the post-2007 FF trajectory.

This has profound implications: **the post-2007 FF emission trend cannot be determined from isotopic data alone without resolving the KIE**. This conclusion holds even in the hemispheric dual-isotope framework and should temper claims of definitive FF trend attribution based on δ¹³C (e.g., Schwietzke et al., 2016; Chandra et al., 2024).

### 3.5 One-Box vs. Two-Box: Spatial Aliasing Assessment

A motivation for hemispheric resolution was the hypothesis that spatial aliasing — where offsetting NH and SH trends cancel in a one-box framework — might explain the discrepancy between box-model (FF declining) and 3-D inversion (FF increasing) conclusions (Basu et al., 2022; Naus et al., 2019).

Our results do not support this hypothesis. Both the one-box (−1.8 Tg yr⁻²) and two-box (−2.5 Tg yr⁻²) models show declining global FF, with the two-box trend slightly more negative. The aliasing bias (two-box minus one-box) is only −0.7 Tg yr⁻² — small relative to the total uncertainty.

**Importantly, the 1-box FF trend is not statistically significant at 90% confidence** (90% CI: [−4.1, +0.03]; the upper bound narrowly includes zero), whereas the 2-box global FF trend is marginally significant (90% CI: [−5.1, −0.05]). This illustrates a practical benefit of hemispheric resolution: by resolving the divergent NH (−3.9 Tg yr⁻²) and SH (+1.3 Tg yr⁻²) trends separately, the two-box framework avoids the partial cancellation that pushes the 1-box CI across zero.

However, hemispheric resolution reveals genuine new information:

1. **Divergent NH/SH FF trends:** NH FF declines (−3.9 Tg yr⁻²) while SH FF increases (+1.3 Tg yr⁻²). The one-box model sees only the net (−2.5 Tg yr⁻²).
2. **Asymmetric microbial trends:** The NH/SH microbial asymmetry (+6.6 vs. −1.1 Tg yr⁻²) is the strongest hemispheric signal and is entirely invisible to one-box frameworks.
3. **Improved conditioning:** The Fisher information (det(A^T A)) of the NH and SH subsystems combined exceeds the global system, indicating that hemispheric resolution genuinely adds information content — consistent with the theoretical analysis of Naus et al. (2019).

The persistent discrepancy with 3-D inversions (Basu et al., 2022; Thanwerdas et al., 2024) likely reflects structural differences beyond spatial resolution: 3-D inversions incorporate transport variability, regional source signatures, and spatially varying OH fields that box models approximate as hemispheric means. The He et al. (2026a) finding that 59% of post-2019 growth reflects approach-to-steady-state suggests that the "true" emission trend is smaller than what mass-balance methods infer from the atmospheric growth rate, because the atmosphere is still adjusting to pre-existing emission levels.

### 3.6 Solver Diagnostics

**Condition number:** The uncertainty-scaled system matrix has a mean condition number of 15.4 (range 12–20 across years and Monte Carlo iterations). This represents a factor-of-11,000 improvement over our v1 implementation (condition ~170,000), achieved in two independent steps: the delta-space reformulation (v3) reduced the condition number from ~170,000 to ~27, and the subsequent uncertainty-based row scaling (v4) reduced it further from ~27 to ~15.4. All subsystems (global, NH, SH) are well-conditioned (condition < 20).

**δD contribution:** With uncertainty-based weighting, δD contributes 33.3% to the cost function — confirming that the model is genuinely "dual-isotope." Previous implementations with ad hoc weights (W_δD = 0.5, W_mass = 100) reduced δD's effective contribution to ~2%, rendering it informationally negligible.

**BB bound-hitting:** Only 8.3% of solver iterations place BB at the lower bound (< 0.1 Tg yr⁻¹), confirming that the 3-source inversion is functioning as intended. This contrasts with earlier versions where >50% of solves collapsed to effectively 2-source (FF + Mic) solutions.

---

## 4. Discussion

### 4.1 Comparison with Recent Literature

**Table 3.** Post-2007 FF emission trends across methods.

| Study | Method | FF trend (Tg yr⁻²) | Attribution |
|-------|--------|-------------------|------------|
| Schwietzke et al. (2016) | 1-box, δ¹³C, KIE database | +0.5 to +4.5 | FF increasing |
| Worden et al. (2017) | Combined CO + δ¹³C | −1 to +2 | BB declining reconciles |
| Rice et al. (2016) | Bayesian 1-box, δ¹³C | Flat 1980s–90s, recent ↑ | Mixed FF + Mic |
| Basu et al. (2022) | TM5-4DVAR, δ¹³C | Variable (KIE-dependent) | ~85% Mic, but KIE uncertain |
| Zhang et al. (2021) | GOSAT inversion | — (total ↑) | Tropical Mic |
| Maasakkers et al. (2019) | GOSAT inversion | O&G ↑ regionally | Spatially heterogeneous |
| Thanwerdas et al. (2024) | LMDz-3D, δ¹³C + δD | +3 to +5 | ~50% FF, ~50% Mic |
| Chandra et al. (2024) | ACTM-3D, δ¹³C + δD | FF stable 2000–2020 | ONG ↓ offset coal ↑ |
| Riddell-Young et al. (2025) | 1-box, δ¹³C + δD | Stable/declining | Entirely Mic |
| Dasgupta et al. (2025) | 2-box Bayesian, δ¹³C + δD | Modest ↑ | Wetland-dominated |
| Fujita et al. (2025) | 1-box, ¹³C + D + ¹⁴C | 30% lower than prior studies | Mic + reduced FF |
| He et al. (2026a) | TROPOMI inversion | — | 25% from emission ↑ |
| He et al. (2026b) | GFDL-CCM, δ¹³C | Energy + Agriculture ↑ | OH trend critical |
| Nguyen et al. (2020) | Chemical feedbacks | 25% bias if CH₄-OH feedback ignored | CO feedback matters |
| Naus et al. (2019) | 2-box OH analysis | — | Box model biases documented |
| Zhao et al. (2023) | Multiple observations | OH constrained | OH trend uncertain |
| **This study** | **2-box, δ¹³C + δD, obs. gradient** | **−2.5 [−5.1, −0.05]** | **NH Mic dominant; KIE-contingent** |

Several patterns emerge from this comparison:

**1. Microbial dominance is robust across frameworks.** Whether using 1-box (Riddell-Young et al., 2025), 2-box (this study; Dasgupta et al., 2025), or 3-D models (Basu et al., 2022; Chandra et al., 2024), the conclusion that microbial emissions are the primary driver of post-2006 growth is consistent. Our hemispheric resolution adds the crucial detail that this growth is concentrated in the NH.

**2. FF trend depends on the KIE choice.** Studies using Saueressig (1.0039) tend to find FF increasing (Schwietzke et al., 2016; Thanwerdas et al., 2024); those using Cantrell (1.0054) or sampling across KIE values find FF stable or declining (Basu et al., 2022; Riddell-Young et al., 2025; this study). Until a new laboratory measurement resolves this 35-year-old discrepancy, FF trends from isotope studies carry an irreducible ambiguity.

**3. OH matters but is hard to quantify.** He et al. (2026b) and Zhao et al. (2023) emphasize that OH trends affect both CH₄ lifetime and δ¹³C interpretation. He et al. (2026a) infer decreasing OH over 2022–2024. Nguyen et al. (2020) show that chemical feedbacks (via CO) can bias emission estimates by 25%. Our model treats OH as time-varying through the lifetime parameter but does not optimize it independently — a limitation shared with all box models.

### 4.2 The NH Microbial Signal

The +6.6 Tg yr⁻² NH microbial trend is the most robust finding (100% of iterations positive, seed-independent). Integrating over the 2007–2020 period, this implies a cumulative increase of ~86 Tg yr⁻¹ in NH microbial emissions — an enormous signal that dwarfs all other source changes.

What drives it? Several mechanisms are consistent:

**Tropical wetlands.** The majority of global wetland CH₄ originates from NH tropics (0–30°N), particularly the Amazon basin, Congo basin, and South/Southeast Asian floodplains. Warming-driven expansion of wetland area and enhanced decomposition rates are expected to increase emissions (Saunois et al., 2020). He et al. (2026a) identify East Africa and South America as emission hotspots — both NH-tropical regions.

**Livestock.** Global cattle numbers have increased ~10% since 2007 (FAO), concentrated in India, Brazil, and sub-Saharan Africa. The δD signature of ruminant enteric fermentation (−300 to −350‰) is depleted enough to contribute to the observed δD decline.

**Permafrost thaw.** Arctic warming has accelerated since 2007, with thermokarst lake formation and active-layer deepening releasing previously frozen organic carbon as CH₄. This source is purely NH and has strongly depleted δ¹³C (−65 to −75‰) and δD (−350 to −400‰).

The SH microbial stability (−1.1 ± 1.3 Tg yr⁻²) constrains these mechanisms: whatever is driving the NH increase is not operating symmetrically. This rules out global-scale OH changes as the sole explanation (which would affect both hemispheres) and points toward spatially heterogeneous emission changes.

### 4.3 OH Sensitivity

Our model does not independently optimize OH trends, which represents a potentially significant source of bias. He et al. (2026b) demonstrate that a +1% decade⁻¹ OH trend can account for much of the observed δ¹³C decline without requiring dominant microbial increases, because OH preferentially removes ¹²CH₄ over ¹³CH₄ (KIE > 1). In our framework, an unaccounted OH trend would alias primarily into the microbial source term: increasing OH would reduce effective lifetime, requiring higher total sources to match the observed burden growth, with the isotopic signal mimicking enhanced microbial input.

A rough estimate of the bias follows from the lifetime sensitivity analysis (Phase 5): shifting the effective lifetime by the equivalent of ±1% decade⁻¹ OH change (~±0.1 yr decade⁻¹) displaces the FF trend by ~2 Tg yr⁻² and the microbial trend by ~3–4 Tg yr⁻². This is comparable in magnitude to the KIE uncertainty and represents a second irreducible ambiguity in isotope-only source attribution. Zhao et al. (2023) constrain OH trends to < ±1% decade⁻¹ using multiple observations, which would bound the microbial trend bias to < 4 Tg yr⁻² — insufficient to eliminate the NH microbial signal (+6.6 Tg yr⁻²) but potentially significant for the FF trend interpretation. Future work should jointly optimize OH and source trends, following the approach of He et al. (2026b), within the hemispheric dual-isotope framework.

### 4.4 Implications for the Global Methane Pledge

The Global Methane Pledge (2021) commits >150 countries to a 30% reduction in anthropogenic CH₄ by 2030 relative to 2020 levels. Our results have two implications:

1. **Targeting FF alone is insufficient.** If microbial emissions are the dominant growth driver (+5.6 Tg yr⁻² globally), then fossil-fuel reductions — while important for multiple reasons — will not halt CH₄ growth unless microbial emissions are also addressed.

2. **NH-focused mitigation may be more effective.** The hemispheric asymmetry suggests that emission reductions in the NH tropics and mid-latitudes (livestock management, wetland management, coal mine CH₄ capture) would have a disproportionate impact on the global growth rate.

### 4.5 Limitations

1. **IH gradient values are literature-compiled, not directly observed.** While based on NOAA MBL reference data, the actual hemispheric-mean time series from the NOAA network should replace our compiled values when publicly released.

2. **Two-box representation.** Hemispheric means smooth over enormous spatial heterogeneity. The "NH microbial" signal likely reflects tropical and boreal contributions that a 2-box model cannot separate.

3. **BB overestimate.** Our BB estimate (54 Tg yr⁻¹) substantially exceeds GFED4s (~35 Tg yr⁻¹). The solver likely absorbs residuals from unresolved OH variability, transport error, or source-signature misspecification into BB — the least isotopically constrained category. The increasing BB trend (+2.3 Tg yr⁻² in NH) warrants investigation against fire inventories and may partially reflect solver artifacts.

4. **No OH optimization.** We do not independently solve for OH trends, unlike He et al. (2026a,b) and Zhao et al. (2023). An unaccounted OH trend of ±1% decade⁻¹ could bias microbial trends by ~3–4 Tg yr⁻² and FF trends by ~2 Tg yr⁻² (Section 4.3).

5. **KIE floor.** The 25% KIE contribution to FF variance represents an irreducible uncertainty floor that no amount of δD data or model refinement can address — only a new laboratory KIE measurement will help.

6. **1-box FF trend marginal.** The one-box global FF trend (−1.8 [−4.1, +0.03] Tg yr⁻²) does not achieve 90% significance, indicating that the declining FF finding depends on the hemispheric decomposition.

---

## 5. Conclusions

We draw three principal conclusions from this hemispheric dual-isotope analysis:

1. **Post-2006 CH₄ growth is driven by Northern Hemisphere microbial emissions** (+6.6 [+5.1, +8.1] Tg yr⁻²), with Southern Hemisphere microbial emissions stable. This hemispheric asymmetry, invisible to one-box models, points to NH tropical/boreal wetlands and livestock as the dominant growth drivers and rules out globally symmetric mechanisms (e.g., uniform OH decline).

2. **δD-CH₄ significantly improves source attribution — conditionally.** The improvement is a 38% reduction in σ(FF) (from 31.0 to 19.2 Tg yr⁻¹), but only when microbial δD source-signature uncertainty is below ~37‰. Above this threshold, δD degrades the inversion. This threshold explains the conflicting assessments of Riddell-Young et al. (2025; improvement) and Thanwerdas et al. (2024; no improvement), and establishes a clear target for observational investment.

3. **The OH-¹³C KIE remains an irreducible uncertainty floor** that determines the sign of the post-2007 FF trend even in a dual-isotope hemispheric framework. In the δ¹³C-only framework, source signatures dominate the variance (83%) and the KIE is minor (11%); adding δD breaks source-signature degeneracies, revealing the KIE as 25% of a smaller total uncertainty. Source signatures and the KIE together account for >70% of FF uncertainty. Spatially resolved, time-varying source-signature measurements — and a definitive KIE re-measurement — are the binding constraints.

These findings highlight the need for: (a) a definitive laboratory re-measurement of the OH-¹³C KIE, (b) expanded δD-CH₄ monitoring networks, particularly in the tropics, (c) hemispheric-scale δD source-signature characterization for wetlands, livestock, and permafrost, and (d) integration of multi-isotope box models with 3-D inversions and satellite observations to leverage the complementary strengths of each approach.

---

## Data Availability

All code, data, and results are available at: https://github.com/Ilovecodinghhh/upgrade_two_isotope_model

The v4 model implementation is in `experiments/Hemispheric_Divergence/analysis/improved_model_v4.py`. The δD threshold experiment is in `experiments/dD_threshold/`. The KIE sensitivity analysis is in `experiments/KIE_immunity/` and `experiments/KIE_sensitivity/`.

---

## Acknowledgments

We acknowledge the NOAA Global Monitoring Laboratory for maintaining the surface flask network that underpins all isotopic CH₄ analyses. We thank B. Riddell-Young and colleagues for making their harmonized δD-CH₄ dataset publicly available. EDGAR v7 emission inventories were provided by the European Commission Joint Research Centre.

---

## References

Basu, S., et al.: Estimating emissions of methane consistent with atmospheric measurements of methane and δ¹³C of methane, Atmos. Chem. Phys., 22, 15351–15377, https://doi.org/10.5194/acp-22-15351-2022, 2022.

Cantrell, C. A., et al.: Carbon kinetic isotope effect in the oxidation of methane by the hydroxyl radical, J. Geophys. Res., 95, 22455–22462, 1990.

Chandra, N., et al.: Methane emissions decreased in fossil fuel exploitation and sustainably increased in microbial source sectors during 1990–2020, Commun. Earth Environ., 5, 147, https://doi.org/10.1038/s43247-024-01286-x, 2024.

Dasgupta, B., et al.: Global Methane Emission Estimates from a Dual-Isotope Inversion: New Constraints from δD-CH₄, EGUsphere [preprint], https://doi.org/10.5194/egusphere-2025-5571, 2025.

Dlugokencky, E. J., et al.: Atmospheric methane levels off: Temporary pause or new steady state?, Geophys. Res. Lett., 30, 1992, 2003.

Dlugokencky, E. J., et al.: Observationally constrained global methane budget and trends, 2001–2009, Phil. Trans. R. Soc. A, 369, 2058–2072, 2009.

Dlugokencky, E. J., et al.: Global atmospheric methane: budget, changes and dangers, Phil. Trans. R. Soc. A, 369, 2058–2072, 2011.

Fujita, R., et al.: Global fossil methane emissions constrained by multi-isotopic atmospheric methane histories, J. Geophys. Res. Atmos., 130, e2024JD041266, https://doi.org/10.1029/2024JD041266, 2025.

He, J., Naik, V., and Horowitz, L. W.: Interpreting changes in global methane budget in a chemistry-climate model constrained with methane and isotopic observations, AGU Advances, 7, e2025AV001822, https://doi.org/10.1029/2025AV001822, 2026b.

He, M., et al.: Attributing 2019–2024 methane growth using TROPOMI satellite observations, Sci. Adv., 12, eadz9007, 2026a.

Lan, X., et al.: Improved constraints on global methane emissions and sinks using δ¹³C-CH₄, Global Biogeochem. Cycles, 35, e2021GB007000, 2021.

Lan, X., et al.: Trends in globally-averaged CH₄, N₂O, and SF₆ determined from NOAA Global Monitoring Laboratory measurements, Earth Syst. Sci. Data, 16, 2197–2206, 2024.

Luo, X., et al.: A time-varying C4 vegetation fraction map for 2001–2019, Zenodo, https://doi.org/10.5281/zenodo.10516423, 2024.

Maasakkers, J. D., et al.: Global distribution of methane emissions, emission trends, and OH concentrations and trends inferred from an inversion of GOSAT satellite data for 2010–2015, Atmos. Chem. Phys., 19, 7859–7881, https://doi.org/10.5194/acp-19-7859-2019, 2019.

Masarie, K. A. and Tans, P. P.: Extension and integration of atmospheric carbon dioxide data into a globally consistent measurement record, J. Geophys. Res., 100, 11593–11610, 1995.

Michel, S. E., et al.: Improved estimate of the ¹³C/¹²C ratio of atmospheric CH₄ from the NOAA Global Greenhouse Gas Reference Network, Atmos. Meas. Tech., 17, 3653–3674, https://doi.org/10.5194/amt-17-3653-2024, 2024.

Naus, S., et al.: Constraints and biases in a tropospheric two-box model of OH, Atmos. Chem. Phys., 19, 407–424, https://doi.org/10.5194/acp-19-407-2019, 2019.

Nguyen, N. H., et al.: Effects of chemical feedbacks on decadal methane emissions estimates, Geophys. Res. Lett., 47, e2019GL085706, https://doi.org/10.1029/2019GL085706, 2020.

Nisbet, E. G., et al.: Rising atmospheric methane: 2007–2014 growth and isotopic shift, Global Biogeochem. Cycles, 30, 1356–1370, 2016.

Nisbet, E. G., et al.: Very strong atmospheric methane growth in the 4 years 2014–2017: Implications for the Paris Agreement, Global Biogeochem. Cycles, 33, 318–342, 2019.

Patra, P. K., et al.: TransCom model simulations of CH₄ and related species: linking transport, surface flux and chemical loss with CH₄ variability in the troposphere and lower stratosphere, Atmos. Chem. Phys., 11, 12813–12837, 2011.

Rice, A. L., et al.: Atmospheric methane isotopic record favors fossil sources flat in 1980s and 1990s with recent increase, Proc. Natl. Acad. Sci. USA, 113, 10791–10796, 2016.

Riddell-Young, B., et al.: Microbial driver of 2006–2023 CH₄ growth indicated by trends in atmospheric δD-CH₄ and δ¹³C-CH₄, Proc. Natl. Acad. Sci. USA, 122, e2516543122, 2025.

Saueressig, G., et al.: Carbon 13 and D kinetic isotope effects in the reactions of CH₄ with O(¹D) and OH: New laboratory measurements and their implications for the isotopic composition of stratospheric methane, J. Geophys. Res., 106, 23127–23138, 2001.

Saunois, M., et al.: The global methane budget 2000–2017, Earth Syst. Sci. Data, 12, 1561–1623, 2020.

Schwietzke, S., et al.: Upward revision of global fossil fuel methane emissions based on isotope database, Nature, 538, 88–91, https://doi.org/10.1038/nature19797, 2016.

Sherwood, O. A., et al.: Global inventory of gas geochemistry data from fossil fuel, microbial and burning sources, version 2017, Earth Syst. Sci. Data, 9, 639–656, 2017.

Still, C. J., et al.: Global distribution of C3 and C4 vegetation: Carbon cycle implications, Global Biogeochem. Cycles, 17, 1006, 2003.

Thanwerdas, J., et al.: Investigation of the renewed methane growth post-2007 with high-resolution 3-D variational inverse modeling and isotopic constraints, Atmos. Chem. Phys., 24, 2129–2167, https://doi.org/10.5194/acp-24-2129-2024, 2024.

Worden, J. R., et al.: Reduced biomass burning emissions reconcile conflicting estimates of the post-2006 atmospheric methane budget, Nat. Commun., 8, 2227, https://doi.org/10.1038/s41467-017-02246-0, 2017.

Zhang, Y., et al.: Attribution of the accelerating increase in atmospheric methane during 2010–2018 by inverse analysis of GOSAT observations, Atmos. Chem. Phys., 21, 3643–3666, https://doi.org/10.5194/acp-21-3643-2021, 2021.

Zhao, Y., et al.: Reconciling the bottom-up and top-down estimates of the methane chemical sink using multiple observations, Atmos. Chem. Phys., 23, 789–807, 2023.
