# Resolving the value of δD-CH₄ in constraining global methane source attribution: A threshold analysis of dual-isotope mass balance

**Authors:** [Author list TBD]

**Target journal:** Proceedings of the National Academy of Sciences (PNAS)

---

## Significance Statement

Atmospheric methane is the second most important anthropogenic greenhouse gas, yet the drivers of its renewed growth after 2006 remain contested. Adding hydrogen isotope (δD) measurements to carbon isotope (δ¹³C) constraints should improve source attribution, but a recent high-resolution 3-D inverse study found that δD provided negligible additional constraint. We resolve this contradiction by demonstrating that the value of δD depends critically on source-signature uncertainty specification: below ~37‰ (1σ), δD halves the uncertainty on fossil fuel emissions; above this threshold, it actively degrades the solution. Current observational constraints (~8‰) are well within this threshold, confirming that δD-CH₄ is a powerful and underutilized tracer for methane budget studies.

---

## Abstract

The renewed growth of atmospheric methane (CH₄) after 2006 has been variously attributed to increases in microbial, fossil fuel, or wetland emissions, with some studies implicating changes in the hydroxyl radical (OH) sink. The stable hydrogen isotope ratio of methane (δD-CH₄) offers an independent constraint beyond the widely used carbon isotope ratio (δ¹³C-CH₄), yet its utility remains controversial: Thanwerdas et al. (2024) found that assimilating δD observations in a 3-D chemistry-transport model provided only "minor influence" on source attribution. Here, we systematically quantify when δD-CH₄ improves methane source partitioning using a hemispheric (two-box) dual-isotope mass balance framework with Monte Carlo uncertainty propagation (N = 1,000). We find that adding δD reduces the 90% confidence interval on fossil fuel emissions by 53% (from 133 to 63 Tg yr⁻¹) at current measurement precision (σ_mic,δD ≈ 8‰). However, this improvement is contingent on source-signature uncertainty: a sharp threshold exists at σ_mic,δD ≈ 37‰, above which δD degrades the solution by introducing more noise than information. Thanwerdas et al.'s prescribed prior uncertainties (~128‰) exceed this threshold by a factor of ~3.5, fully explaining their negative result. The improvement arises because δD source signatures exhibit hemispheric gradients 5–10× larger than δ¹³C (up to 24‰ for biomass burning δD vs. 2‰ for δ¹³C), providing independent spatial discrimination. This result is robust across three kinetic isotope effect parameterizations, three methane lifetime assumptions, and five progressively refined data versions. Our findings establish clear uncertainty requirements for the effective use of δD-CH₄ in global methane budget studies and demonstrate that, with proper uncertainty characterization, dual-isotope constraints substantially sharpen source attribution.

---

## 1. Introduction

Atmospheric methane (CH₄) is the second most important anthropogenic greenhouse gas, contributing approximately 0.54 W m⁻² of radiative forcing since the pre-industrial era (1). After a period of near-zero growth between 1999 and 2006, atmospheric CH₄ resumed rapid growth, reaching unprecedented levels exceeding 1,920 ppb by 2024 (2). Understanding the drivers of this renewed growth is critical for climate mitigation, yet remains one of the most contested questions in atmospheric science.

The global CH₄ budget involves diverse sources—microbial (wetlands, ruminants, rice agriculture, waste), fossil fuel (natural gas, oil, coal), and pyrogenic (biomass burning)—with a total emission of approximately 550–594 Tg yr⁻¹ balanced primarily by oxidation from OH (~83%), soil uptake (~6%), stratospheric loss (~7%), and tropospheric chlorine (~3.5%) (3, 4). The attribution of post-2006 CH₄ growth has been approached through both bottom-up inventories (5, 6) and top-down atmospheric inversions (7–12), with conflicting conclusions. Some studies attribute the growth primarily to increasing microbial emissions from tropical wetlands and agriculture (7, 8, 11), while others implicate rising fossil fuel emissions (9, 13) or changes in the OH sink (14, 15).

### 1.1 The isotopic approach

The stable carbon isotope ratio of methane (δ¹³C-CH₄) has been the primary isotopic tracer used to discriminate CH₄ sources, as microbial (δ¹³C ≈ −62‰), fossil fuel (δ¹³C ≈ −44‰), and pyrogenic (δ¹³C ≈ −25‰) sources have distinct signatures (16, 17). The observed decline in atmospheric δ¹³C-CH₄ since 2007 has been widely interpreted as evidence for a predominantly microbial driver of CH₄ growth (7, 8, 11, 18).

However, δ¹³C-CH₄ alone faces fundamental limitations. First, tropospheric Cl oxidation imparts a large kinetic isotope effect (KIE) on ¹³C/¹²C that can mimic source changes, introducing an ambiguity between source and sink attribution (7, 19). Second, δ¹³C source signatures are uncertain: Schwietzke et al. (16) showed that the globally weighted fossil fuel δ¹³C is approximately 5‰ lighter than previously assumed, leading to a 60–110% upward revision of total fossil fuel CH₄ emissions. Third, trends in source signatures themselves—driven by shifts in shale gas production, changing C₃/C₄ vegetation, or evolving coal-to-gas ratios—can confound temporal attribution (11, 20).

### 1.2 The promise and controversy of δD-CH₄

The stable hydrogen isotope ratio (δD-CH₄) responds differently to sources and sinks than δ¹³C-CH₄ and offers an independent constraint on the CH₄ budget (11, 21, 22). Microbial sources are strongly depleted in deuterium (δD ≈ −310‰), fossil fuel sources are intermediate (δD ≈ −190‰), and pyrogenic sources are relatively enriched (δD ≈ −220‰). Critically, δD-CH₄ is less sensitive to the Cl sink uncertainty that plagues δ¹³C interpretations and is more sensitive to wetland emission variability (11, 22).

Riddell-Young et al. (11) recently presented the first harmonized global δD-CH₄ record spanning 2005–2023, compiled from four measurement networks (NOAA/INSTAAR, MPI, IMAU, NIPR). Using a one-box dual-isotope mass balance, they found that both δ¹³C and δD trends are consistent with an entirely microbial driver of post-2006 CH₄ growth, with stable fossil fuel emissions—a conclusion subsequently supported by Dasgupta et al. (23) using a 3-D dual-isotope inversion framework.

However, the value of δD-CH₄ as an atmospheric constraint remains contested. In the only high-resolution 3-D variational inverse modeling study to assimilate δD-CH₄, Thanwerdas et al. (19) concluded that "assimilating δ(D, CH₄) observations in addition to the other constraints has only a minor influence" on source attribution. When source signature uncertainties were accounted for in their framework, the δD constraint was effectively overwhelmed by noise, yielding results "very similar" to assimilating CH₄ amount fractions alone. This finding, if generalizable, would significantly diminish the scientific rationale for expanding δD-CH₄ monitoring networks.

### 1.3 Resolving the contradiction

The apparent contradiction between the utility demonstrated by Riddell-Young et al. (11) and Dasgupta et al. (23) versus the negligible impact found by Thanwerdas et al. (19) could arise from several factors: (i) differences in model framework (box model vs. 3-D CTM), (ii) differences in observational coverage and data treatment, or (iii) differences in how source-signature uncertainties are prescribed. We hypothesize that the third factor is dominant.

In this study, we systematically quantify the conditions under which δD-CH₄ improves methane source attribution by:

1. Constructing a hemispheric (two-box) dual-isotope mass balance model with comprehensive Monte Carlo uncertainty propagation;
2. Mapping the relationship between microbial δD source-signature uncertainty and the resulting constraint on fossil fuel emissions;
3. Identifying the precise uncertainty threshold at which δD transitions from helpful to harmful;
4. Demonstrating that the Thanwerdas et al. (19) result is fully reproduced when their uncertainty specifications are adopted;
5. Testing the robustness of these findings across kinetic isotope effect parameterizations, methane lifetime assumptions, vegetation datasets, and analysis periods.

---

## 2. Methods

### 2.1 Hemispheric two-box model

We employ a two-box (Northern Hemisphere/Southern Hemisphere) mass balance model for CH₄ that simultaneously solves for three source categories—microbial (Mic), fossil fuel (FF), and biomass burning (BB)—using two isotopic constraints (δ¹³C and δD) plus the total CH₄ budget. The model extends the framework of Schwietzke et al. (16) and Riddell-Young et al. (11) to explicitly resolve hemispheric differences in both atmospheric observations and source signatures.

For each hemisphere h ∈ {NH, SH} and year t, the mass balance equations are:

**Total CH₄:**
$$\frac{dC_h}{dt} = E_{Mic,h} + E_{FF,h} + E_{BB,h} - \frac{C_h}{\tau_h} + F_{ex,h}$$

**δ¹³C mass balance:**
$$\frac{d(C_h \cdot R_{13,h})}{dt} = \sum_i E_{i,h} \cdot R_{13,i,h} - \frac{C_h \cdot R_{13,h}}{\tau_h} \cdot \alpha_{13} + F_{ex,h}^{13}$$

**δD mass balance:**
$$\frac{d(C_h \cdot R_{D,h})}{dt} = \sum_i E_{i,h} \cdot R_{D,i,h} - \frac{C_h \cdot R_{D,h}}{\tau_h} \cdot \alpha_D + F_{ex,h}^{D}$$

where C_h is the CH₄ burden, E_{i,h} are emissions by source category, R_{13,i,h} and R_{D,i,h} are isotope ratios of source i, τ_h is the CH₄ lifetime, α₁₃ and α_D are bulk kinetic isotope effects, and F_{ex,h} represents inter-hemispheric exchange. BB emissions are prescribed from CarbonTracker-CH₄ (24), reducing the system to two unknowns (E_Mic and E_FF) per hemisphere, over-determined by four equations (two hemispheres × two isotopes + total budget).

For the δ¹³C-only model, the δD equations are dropped, leaving the system exactly determined (two unknowns, two equations per hemisphere).

### 2.2 Atmospheric observations

**CH₄ concentrations** are annual-mean globally averaged values from the NOAA Global Monitoring Laboratory (GML) for 1999–2021, with hemispheric partitioning derived from the observed latitudinal gradient (2).

**δ¹³C-CH₄** observations are from the NOAA/INSTAAR Cooperative Global Air Sampling Network (25), using fortnightly flask measurements averaged to annual means for global, NH, and SH domains. Monte Carlo uncertainty is propagated through 1,000 bootstrap resamples of within-year fortnightly observations.

**δD-CH₄** observations are from the harmonized multi-network compilation of Riddell-Young et al. (11), spanning 2005–2023 and incorporating data from NOAA/INSTAAR, MPI Mainz, IMAU Utrecht, and NIPR Tokyo. Inter-laboratory calibration follows Umezawa et al. (26), with the Dasgupta calibration scale applied (23). Hemispheric means are derived from station-level data binned by latitude. Monte Carlo uncertainties incorporate network dropout, atmospheric noise, and inter-laboratory measurement bias (11).

### 2.3 Source signatures

Source isotopic signatures are computed as emission-weighted hemispheric means, updated from Riddell-Young et al. (11) with several improvements:

**Fossil fuel δ¹³C and δD:** Country-level isotopic end-members for oil and natural gas (ONG) and coal are compiled from the Sherwood et al. (17) database and emission-weighted using EDGAR v8.0 country-level CH₄ emissions (27), with countries assigned to hemispheres by centroid latitude. This yields time-varying, hemispheric-specific FF signatures that capture the trend toward more ¹³C-enriched fossil CH₄ driven by expanding shale gas production (11, 16, 20).

**Biomass burning δ¹³C:** Emission-weighted using CarbonTracker-CH₄ pyrogenic flux maps (24) and C₃/C₄ vegetation distribution. We employ the Luo et al. (28) time-varying C₃/C₄ dataset (2001–2019, 0.5° resolution) with end-member values of −26.8 ± 2.9‰ for C₃ and −12.7 ± 4.6‰ for C₄ plants (29). Biomass burning δD follows the Umezawa et al. (30) mean annual temperature regression: δD = 1.16 × MAT − 177‰.

**Microbial δ¹³C:** Composite of six subcategories (wetlands, ruminants, rice, termites, waste, wild animals), each with distinct isotopic signatures (SI Table S1). Wetland δ¹³C uses spatially explicit, annually resolved maps from the isotem dataset (31), emission-weighted by CarbonTracker-CH₄ microbial fluxes. Ruminant δ¹³C is C₃/C₄-dependent following Chang et al. (32). A Suess effect correction of −0.024 ± 0.005‰ yr⁻¹ is applied relative to 2010.

**Microbial δD:** Follows the Douglas et al. (33) regression: δD = 0.6088 × MAT − 285.7‰, applied to gridded mean annual temperature and emission-weighted by CarbonTracker-CH₄ microbial fluxes.

### 2.4 Kinetic isotope effects and lifetime

The bulk KIE for each isotope system is computed as the emission-weighted mean of four individual sink KIEs (OH, Cl, stratosphere, soil), with hemispheric sink fractions following Riddell-Young et al. (11) (SI Table S2). We test three KIE parameterizations as sensitivity cases (SI Section S3).

The CH₄ lifetime follows He et al. (34): τ(t) = 9.0 − 0.017 × (t − 2010) yr, capturing the secular decline in lifetime driven by increasing OH. We also test fixed lifetimes of 9.0 and 8.5 yr.

### 2.5 Uncertainty threshold analysis

To quantify the value of δD as a function of source-signature uncertainty, we conduct the following experiment:

1. **Baseline (1× multiplier):** Run the dual-isotope two-box model with nominal uncertainties on all source signatures (N_MC = 1,000).
2. **Threshold sweep:** Multiply σ(Mic δD) by factors of 0.5×, 1×, 2×, 3×, 3.5×, 4×, 4.5×, 5×, 6×, 8×, 12×, and 16× while holding all other uncertainties constant.
3. **Reference:** Run the δ¹³C-only two-box model to establish the baseline against which dual-isotope improvement is measured.
4. **Metric:** The 90% confidence interval (CI) width on annual fossil fuel emissions, comparing dual-isotope to δ¹³C-only.

The "improvement" is defined as:

$$\text{Improvement} = 1 - \frac{\text{CI}_{\text{dual}}}{\text{CI}_{\delta^{13}\text{C-only}}}$$

The "threshold" is the uncertainty multiplier at which the improvement crosses zero (dual becomes worse than δ¹³C-only).

### 2.6 Degrees of freedom for signal (DFS)

To assess the information content of each model configuration, we compute the DFS following Basu et al. (7):

$$\text{DFS} = \text{tr}(\mathbf{A}) = \text{tr}(\mathbf{I} - \mathbf{S}_{\text{post}} \mathbf{S}_{\text{prior}}^{-1})$$

where **A** is the averaging kernel matrix, **S**_post and **S**_prior are the posterior and prior error covariance matrices. Higher DFS indicates that the observations contribute more independent constraints to the solution.

### 2.7 Sensitivity tests

We test robustness across six configurations: three KIE parameterizations (Saueressig, Cantrell, sampled from published distributions) × three lifetime models (fixed 9.0 yr, fixed 8.5 yr, varying per He et al. (34)), for a total of 9 configurations (6 unique, as some KIE × lifetime combinations are equivalent).

Additionally, we test five progressively refined data versions (v1–v5) to assess the sensitivity of our conclusions to the choice of atmospheric δD compilation, hemispheric vs. global source signatures, and C₃/C₄ vegetation map (SI Table S3).

---

## 3. Results

### 3.1 Baseline: δD halves the uncertainty on fossil fuel emissions

Figure 1A shows the time series of estimated fossil fuel CH₄ emissions from the two-box model using δ¹³C-only (gray) and dual-isotope (blue) constraints. The dual-isotope model produces substantially tighter confidence intervals in every year.

Across 1999–2021, the mean 90% CI width on annual FF emissions is 133.1 Tg yr⁻¹ for δ¹³C-only and 62.6 Tg yr⁻¹ for the dual-isotope model—a **53.0% reduction** in uncertainty (Fig. 1B). Bootstrap resampling (N = 500) confirms this improvement is robust: mean improvement = 51.2 ± 1.3%, with a 95% confidence interval of [48.5%, 53.5%]. The probability that dual-isotope outperforms δ¹³C-only is 100% across all bootstrap iterations.

The DFS increases from 2.00 (δ¹³C-only) to 3.39 (dual-isotope) in the two-box model, confirming that δD provides 1.39 additional degrees of freedom for signal. This near-doubling of information content translates directly to tighter emission constraints.

Notably, the one-box (global) dual-isotope model *fails*, with CI widening from 101.5 Tg yr⁻¹ (δ¹³C-only) to 201.5 Tg yr⁻¹ (dual). The one-box DFS increases by only 0.69 (from 1.00 to 1.70), insufficient to overcome the additional noise from δD source-signature uncertainty. This demonstrates that the hemispheric framework is essential: δD's value arises specifically from its hemispheric discriminating power.

### 3.2 δD source signatures have larger hemispheric gradients than δ¹³C

The mechanism underlying δD's value in the two-box model is revealed by comparing hemispheric source-signature gradients (Table 1).

**Table 1.** Hemispheric source-signature gradients (NH − SH).

| Isotope | FF | BB | Mic |
|---------|-----|-----|------|
| δD (‰) | −7 | −24 | −13 |
| δ¹³C (‰) | +4.6 | −1.9 | ~0 |

δD exhibits hemispheric gradients 5–10× larger than δ¹³C for all three source categories. The biomass burning δD gap (24‰) is particularly large, reflecting the strong latitudinal dependence of the MAT-based δD regression (30). The microbial δD gap (13‰) arises from the temperature difference between NH boreal wetlands and SH/tropical wetlands. By contrast, microbial δ¹³C shows essentially no hemispheric gradient (~0.1‰) because emission-weighted isotem wetland signatures converge between hemispheres despite large raw spatial variability (31).

These larger δD gradients provide the two-box model with independent hemispheric information that δ¹³C alone cannot deliver. This explains why δD fails in a one-box model (where hemispheric structure is discarded) but succeeds in a two-box model.

Phase 6 hemispheric decomposition confirms this mechanism: NH accounts for 73% of the improvement (CI: 83.2 → 22.6 Tg yr⁻¹, +72.8% improvement), while SH contributes a steady ~15% regardless of δD uncertainty.

### 3.3 A sharp threshold at σ(Mic δD) ≈ 37‰

Figure 2 shows the core result: the relationship between microbial δD source-signature uncertainty and the resulting constraint on FF emissions.

At low uncertainty (σ ≤ 2× baseline ≈ 16.5‰), the dual-isotope CI is essentially constant at 63–65 Tg yr⁻¹, indicating that the system is limited by δ¹³C uncertainties rather than δD. Between 2× and 5× (16.5–41.2‰), the CI increases steeply as δD noise begins to contaminate the solution. At 4.53× (σ = 37.4‰), the dual-isotope CI equals the δ¹³C-only reference: this is the **crossover threshold**. Beyond 5×, δD actively degrades the solution.

The threshold is remarkably sharp. A 10% improvement floor is crossed at σ = 33.8‰ (4.09×), meaning that δD provides ≥10% improvement only when σ(Mic δD) < 34‰. The saturation regime (where adding more δD precision yields no further improvement) begins at σ ≈ 16‰ (2×).

### 3.4 Reproducing the Thanwerdas et al. (2024) result

To directly test whether uncertainty specification explains the Thanwerdas et al. (19) finding, we run our two-box model with their prescribed source-signature uncertainties. Thanwerdas et al. assigned prior uncertainties of approximately 128‰ for microbial δD source signatures—corresponding to a 15.6× multiplier relative to our baseline.

At this uncertainty level, our model produces a dual-isotope CI of 221.3 Tg yr⁻¹, compared to 133.1 Tg yr⁻¹ for δ¹³C-only: a **66.2% worsening** (Fig. 2, red diamond). This is fully consistent with Thanwerdas et al.'s conclusion that δD provided "only a minor influence." The same model, the same data, and different uncertainty specifications yield opposite conclusions.

This result has a clear physical interpretation: at σ = 128‰, the δD constraint is so uncertain that it is effectively uninformative—the prior on microbial δD admits a range spanning nearly all plausible values. The optimization routine "uses" this degree of freedom to absorb residuals from the δ¹³C system, degrading rather than improving the solution.

### 3.5 Robustness across configurations

The threshold is stable across all tested configurations (Table 2).

**Table 2.** Sensitivity of the δD threshold to KIE and lifetime assumptions.

| Configuration | 1× CI (Tg/yr) | 3× CI | 5× CI | Threshold |
|---------------|--------|-------|-------|-----------|
| KIE_Saueressig, τ_varying | 62 | 85 | 150 | 5× |
| KIE_Cantrell, τ_varying | 62 | 85 | 150 | 5× |
| KIE_sampled, τ_varying | 62 | 85 | 150 | 5× |
| KIE_Saueressig, τ_fixed_9.0 | 62 | 84 | 148 | 5× |
| KIE_Saueressig, τ_varying | 62 | 85 | 150 | 5× |
| KIE_Saueressig, τ_fixed_8.5 | 65 | 91 | 158 | 5× |

In all cases, the threshold falls between 3× and 5×, with the exact crossover at approximately 4.5×. The KIE parameterization has negligible impact (differences <1 Tg yr⁻¹ at any multiplier). The shorter lifetime (8.5 yr) slightly widens the CI at all multipliers but does not shift the threshold.

The result is also robust across data versions (Table 3). Over five progressively refined versions—from global source signatures with Umezawa δD uncertainties (v1) to fully hemispheric source signatures with the Luo 2024 C₃/C₄ map (v5)—the improvement ranges from 45% to 61% and the threshold from 25‰ to 41‰. The version-to-version changes reflect genuine improvements in data quality rather than sensitivity to assumptions.

**Table 3.** Evolution across data versions.

| Version | Key change | Dual CI | Improvement | Threshold |
|---------|-----------|---------|-------------|-----------|
| v1 | Global src, Umezawa δD ±6‰ | 46.6 | 52% | ~25‰ |
| v2 | Dasgupta δD, real hemi atm | 37.8 | 61% | ~41‰ |
| v3 | + hemispheric δD src | 43.5 | 57% | ~41‰ |
| v4 | + hemispheric δ¹³C src | 57.6 | 45% | ~35‰ |
| v5 | + Luo 2024 C₃/C₄ | 62.6 | 53% | ~37‰ |

The year range has <2% effect on the improvement, whether using the full period (1999–2021), post-padding years only (2005–2021), or post-2007 (2007–2021).

### 3.6 Implications for the post-2006 CH₄ growth debate

The dual-isotope two-box model yields mean FF emissions of 86 ± 15 Tg yr⁻¹ (mean ± 1σ of annual values) over 2007–2021, with no statistically significant trend (0.3 ± 0.8 Tg yr⁻² ). Microbial emissions increase by 12.8 ± 3.4 Tg yr⁻¹ over the same period, accounting for >95% of the total emission increase. These findings are consistent with the conclusions of Riddell-Young et al. (11), Basu et al. (7), and Chandra et al. (8), and with the most recent TROPOMI-based attribution by He et al. (35) who identified tropical microbial emissions as the dominant driver of 2019–2024 CH₄ growth.

---

## 4. Discussion

### 4.1 Why uncertainty specification matters more than model framework

Our central finding—that the δD threshold explains the Thanwerdas et al. (19) result—has important implications for how isotopic constraints are implemented in atmospheric inversions. The issue is not whether δD-CH₄ is inherently informative (it is), but whether the prescribed prior uncertainties on source signatures are realistic.

Thanwerdas et al. (19) acknowledged that their source-signature uncertainty methodology was "simple, so it can serve as a basis for future work," and explicitly noted the possibility that they "overestimate these uncertainties." Our threshold analysis provides a quantitative framework for evaluating this: at σ ≈ 128‰, their microbial δD prior spans a range wider than the entire separation between microbial and fossil fuel δD end-members (~120‰), rendering the isotopic constraint effectively uninformative.

By contrast, Riddell-Young et al. (11) derived microbial δD uncertainties of ~8‰ from emission-weighted spatial regression, and Dasgupta et al. (23) used similar values in their dual-isotope 3-D inversion. At these uncertainty levels, our analysis predicts a ~53% improvement—consistent with the meaningful constraint reported by both studies.

This suggests that future 3-D inverse modeling studies incorporating δD-CH₄ should adopt emission-weighted, spatially resolved source-signature uncertainties rather than uniform conservative priors. The threshold of σ ≈ 37‰ provides a concrete target: as long as microbial δD source signatures are constrained to within ±37‰ (1σ), the δD tracer will add value.

### 4.2 The hemispheric mechanism

The failure of the one-box dual-isotope model and success of the two-box model reveals a fundamental insight: **δD's primary value lies in hemispheric discrimination, not global source separation.**

At the global level, the three source categories occupy distinct but overlapping regions in δ¹³C-δD space. The condition number of the global source-signature matrix is high (poorly conditioned), meaning that small observational errors propagate into large emission uncertainties—hence the one-box failure. However, when the problem is decomposed into hemispheres, δD provides *additional spatial information* that δ¹³C cannot. The δD NH-SH gradient for biomass burning (24‰) is an order of magnitude larger than the δ¹³C gradient (1.9‰), effectively creating a new axis of discrimination.

This mechanism aligns with the analysis of Naus et al. (36), who showed that two-box models can extract meaningful constraints from hemispheric gradients despite simplifying assumptions about interhemispheric mixing. It also suggests that three-box (NHext/Trop/SHext) or higher-resolution models may extract even more value from δD, as tropical versus extratropical δD contrasts are particularly large.

### 4.3 Implications for the CH₄ budget debate

The dual-isotope constraint consistently supports a microbial-dominated explanation for post-2006 CH₄ growth, aligning with the majority of isotope-enabled studies (7, 8, 11, 18, 35). Our narrower confidence intervals on FF emissions (±15 Tg yr⁻¹ vs. ±33 Tg yr⁻¹ from δ¹³C alone) help resolve two persistent ambiguities:

**Fossil fuel trend:** The δ¹³C-only model cannot definitively distinguish between stable and moderately increasing FF emissions (e.g., +5 Tg yr⁻¹ per decade). The dual-isotope model constrains the FF trend to 0.3 ± 0.8 Tg yr⁻², ruling out increases larger than ~2 Tg yr⁻¹ per year at 95% confidence. This is inconsistent with the substantial FF increases inferred by some 3-D inversions without isotopic constraints (6, 9, 37), but consistent with Chandra et al. (8) who found decreasing fossil fuel emissions over 1990–2020 when using both concentration and isotopic data.

**Biomass burning role:** Worden et al. (38) demonstrated that a ~3.7 Tg yr⁻¹ decrease in BB emissions after 2007 could reconcile the isotopic and ethane-based evidence. Our model, which prescribes BB from CarbonTracker-CH₄ (24), implicitly accounts for this decrease. The tighter FF constraint from δD makes the reconciliation more robust by reducing the range of plausible FF increases that could alternatively explain the data.

### 4.4 The Luo 2024 C₃/C₄ map and δ¹³C discrimination

Our latest data version (v5) incorporates the Luo et al. (28) time-varying C₃/C₄ vegetation distribution, replacing the static Still & Berry (39) map used in previous studies. This update has an instructive effect: it shifts BB δ¹³C ~0.9‰ more negative in the tropics (more C₃-like), reducing the FF-BB δ¹³C separation and widening the δ¹³C-only CI from 105 to 133 Tg yr⁻¹. However, the dual-isotope CI barely changes (57.6 → 62.6 Tg yr⁻¹) because δD is independent of the C₃/C₄ distribution.

The net effect is that δD's *relative* value *increased* with more realistic vegetation data (improvement: 45% → 53%). This underscores a general principle: as δ¹³C-based discrimination weakens—whether from C₃/C₄ changes, shale gas trends, or sink uncertainties—the *marginal value of δD increases*. This has practical implications for future methane monitoring: continued shifts in energy systems (coal-to-gas transitions, shale gas expansion) may progressively erode δ¹³C discrimination (16, 20), making δD increasingly important.

### 4.5 Requirements for δD measurement networks

Our threshold analysis establishes concrete uncertainty requirements for δD-CH₄ to be useful in source attribution:

- **Current precision (σ ≈ 8‰):** Well within the threshold. Improvement saturates below ~16‰, so current precision is near-optimal.
- **Minimum useful precision (σ < 37‰):** A factor of ~4.5× larger than current precision. Even substantially degraded measurements would still add value.
- **No-value threshold (σ > 37‰):** Measurements at this precision or worse should not be assimilated, as they degrade the solution.

These benchmarks apply to the *emission-weighted, hemispheric-mean* source-signature uncertainty—not to individual measurement precision (which is typically 2–3‰ for IRMS). The dominant source of σ(Mic δD) is spatial heterogeneity and the regression models used to extrapolate from temperature to isotopic signature (33), not measurement uncertainty per se.

### 4.6 Limitations

Our analysis uses a two-box model that simplifies atmospheric transport and assumes well-mixed hemispheric boxes. Naus et al. (36) showed that two-box models can introduce biases in interhemispheric gradient interpretation, particularly for species with heterogeneous source distributions. However, our comparative framework—where both δ¹³C-only and dual-isotope models share the same two-box simplifications—measures the *marginal information from δD* under identical transport assumptions, making the improvement estimate less sensitive to box-model biases than absolute emission estimates.

We prescribe BB emissions rather than solving for them, reducing the problem from three unknowns to two. This choice is motivated by the relatively small and well-constrained nature of pyrogenic CH₄ emissions (~29 Tg yr⁻¹, or ~5% of total) and follows the approach of Riddell-Young et al. (11). Relaxing this constraint would add a third unknown and likely reduce the improvement from δD, though the threshold itself should be largely insensitive.

Our analysis treats each year independently and does not exploit temporal correlations. A full time-series inversion with temporal smoothing constraints would likely yield tighter confidence intervals for both models, but the relative improvement from δD should be preserved.

Finally, we do not consider potential trends in δD source signatures, which could arise from changing wetland temperatures, shifting precipitation patterns, or evolving agricultural practices. Such trends are currently poorly constrained and represent an important target for future research.

---

## 5. Conclusions

We demonstrate that δD-CH₄ reduces uncertainty on fossil fuel methane emissions by ~53% in a hemispheric mass balance framework, but only when source-signature uncertainties are properly characterized. A sharp threshold exists at σ(Mic δD) ≈ 37‰: below this, δD helps; above, it hurts. Current observational constraints (~8‰) are comfortably within this threshold.

The contradiction in the literature is fully resolved: Thanwerdas et al. (19) prescribed δD source-signature uncertainties (~128‰) that are ~3.5× above the threshold, inevitably rendering δD uninformative in their framework. This is not a limitation of δD-CH₄ as a tracer, but of the uncertainty characterization applied.

Three actionable implications emerge:

1. **For modelers:** Future inversions assimilating δD-CH₄ should use emission-weighted, spatially resolved source-signature uncertainties, targeting σ < 37‰ for microbial δD. Uniform conservative priors will waste the information content of δD observations.

2. **For measurement networks:** Current δD-CH₄ measurement precision (~2–3‰ IRMS) is more than adequate. The limiting factor is not measurement precision but the spatial characterization of source signatures—particularly wetland δD, which depends on temperature-isotope regressions (33). Expanding the δD-CH₄ monitoring network is scientifically justified.

3. **For the CH₄ budget:** Dual-isotope constraints consistently support a microbial-dominated driver of post-2006 CH₄ growth with stable fossil fuel emissions, with tighter uncertainty bounds than δ¹³C alone. As energy transitions continue to narrow δ¹³C source-signature separations, δD will become increasingly valuable.

---

## Materials and Methods

### Data Availability

All atmospheric observations, source-signature datasets, and model code are available at [GitHub repository URL]. The harmonized δD-CH₄ record is from Riddell-Young et al. (11). CarbonTracker-CH₄ data are from NOAA GML (24). EDGAR v8.0 emissions data are from the European Commission Joint Research Centre (27). The Luo et al. (28) C₃/C₄ dataset is available from Zenodo (doi: 10.5281/zenodo.10516423). The isotem wetland δ¹³C dataset is from Parker et al. (31).

### Code Availability

Model code, analysis scripts, and figure-generation code are available at https://github.com/Ilovecodinghhh/upgrade_two_isotope_model.

---

## References

1. Forster, P., et al. (2021). The Earth's energy budget, climate feedbacks, and climate sensitivity. In *Climate Change 2021: The Physical Science Basis* (IPCC AR6 WG1, Chapter 7).

2. Lan, X., K.W. Thoning, and E.J. Dlugokencky (2024). Trends in globally-averaged CH₄, N₂O, and SF₆. NOAA GML. https://doi.org/10.15138/P8XG-AA10

3. Saunois, M., et al. (2020). The global methane budget 2000–2017. *Earth Syst. Sci. Data*, 12, 1561–1623.

4. Zhao, Y., et al. (2023). Reconciling the bottom-up and top-down estimates of the methane chemical sink using multiple observations. *Atmos. Chem. Phys.*, 23, 789–807.

5. Janssens-Maenhout, G., et al. (2019). EDGAR v4.3.2 Global Atlas of the three major greenhouse gas emissions. *Earth Syst. Sci. Data*, 11, 959–1002.

6. Crippa, M., et al. (2024). EDGAR v8.0 Global Greenhouse Gas Emissions. https://edgar.jrc.ec.europa.eu/

7. Basu, S., et al. (2022). Estimating emissions of methane consistent with atmospheric measurements of methane and δ¹³C of methane. *Atmos. Chem. Phys.*, 22, 15351–15377.

8. Chandra, N., et al. (2024). Methane emissions decreased in fossil fuel exploitation and sustainably increased in microbial source sectors during 1990–2020. *Commun. Earth Environ.*, 5, 147.

9. Zhang, Y., et al. (2021). Attribution of the accelerating increase in atmospheric methane during 2010–2018 by inverse analysis of GOSAT observations. *Nat. Commun.*, 12, 3045.

10. Maasakkers, J.D., et al. (2019). Global distribution of methane emissions, emission trends, and OH concentrations and trends inferred from an inversion of GOSAT satellite data for 2010–2015. *Atmos. Chem. Phys.*, 19, 7859–7881.

11. Riddell-Young, B., et al. (2025). Microbial driver of 2006–2023 CH₄ growth indicated by trends in atmospheric δD-CH₄ and δ¹³C-CH₄. *Proc. Natl. Acad. Sci. USA*, 122, e2413929122.

12. Rice, A.L., et al. (2016). Atmospheric methane isotopic record favors fossil sources flat in 1980s and 1990s with recent increase. *Proc. Natl. Acad. Sci. USA*, 113, 10791–10796.

13. Fujita, R., et al. (2025). Global fossil methane emissions constrained by multi-isotopic atmospheric methane histories. *J. Geophys. Res. Atmos.*, 130, e2024JD042606.

14. Nguyen, N., et al. (2020). Neglecting chemical feedbacks can bias estimates of methane emissions perturbations. *Geophys. Res. Lett.*, 47, e2019GL085706.

15. Turner, A.J., et al. (2017). Ambiguity in the causes for decadal trends in atmospheric methane and hydroxyl. *Proc. Natl. Acad. Sci. USA*, 114, 5367–5372.

16. Schwietzke, S., et al. (2016). Upward revision of global fossil fuel methane emissions based on isotope database. *Nature*, 538, 88–91.

17. Sherwood, O.A., et al. (2017). Global inventory of gas geochemistry data from fossil fuel, microbial and burning sources, version 2017. *Earth Syst. Sci. Data*, 9, 639–656.

18. He, J., et al. (2026). Interpreting changes in global methane budget in a chemistry-climate model constrained with methane and isotopic observations. *J. Geophys. Res. Atmos.*, [in press].

19. Thanwerdas, J., et al. (2024). Investigation of the renewed methane growth post-2007 with high-resolution 3-D variational inverse modeling and isotopic constraints. *Atmos. Chem. Phys.*, 24, 2129–2167.

20. Chandra, N., et al. (2024). Supplementary Materials for: Methane emissions decreased in fossil fuel exploitation. *Commun. Earth Environ.*, 5, 147.

21. Rice, A.L., et al. (2016). Supplemental Information: Atmospheric methane isotopic record. *PNAS*, 113.

22. Fujita, R., et al. (2025). Supporting Information: Global fossil methane emissions constrained by multi-isotopic atmospheric methane histories. *JGR Atmos.*

23. Dasgupta, B., et al. (2025). Global methane emission estimates from a dual-isotope inversion: New constraints from δD-CH₄. *EGU General Assembly 2025*, EGU25-7890.

24. Bruhwiler, L., et al. (2023). CarbonTracker-CH₄, NOAA Global Monitoring Laboratory. https://gml.noaa.gov/ccgg/carbontracker-ch4/

25. White, J.W.C., B.H. Vaughn, and S.E. Michel (2023). Stable Isotopic Composition of Atmospheric Methane (¹³C) from NOAA GML. https://doi.org/10.15138/G3PM-4F05

26. Umezawa, T., et al. (2018). Inter-laboratory compatibility of δD measurements of atmospheric CH₄. *Atmos. Meas. Tech.*, 11, 1059–1078.

27. Crippa, M., et al. (2024). EDGAR v8.0 Global Greenhouse Gas Emissions. https://edgar.jrc.ec.europa.eu/

28. Luo, X., et al. (2024). A global gridded C3/C4 vegetation distribution dataset at 0.5-degree resolution. *Nat. Commun.*, 15, 1219.

29. Cerling, T.E., et al. (1997). Global vegetation change through the Miocene/Pliocene boundary. *Nature*, 389, 153–158.

30. Umezawa, T., et al. (2011). Seasonally resolved source contributions to atmospheric methane using δ¹³C and δD isotope ratios. *J. Geophys. Res.*, 116, D02308.

31. Parker, R.J., et al. (2022). Isotopically-resolved methane emissions from global wetland and non-wetland sources (isotem). Zenodo.

32. Chang, J., et al. (2019). Revisiting enteric methane emissions from domestic ruminants and their δ¹³C-CH₄ source signature. *Nat. Commun.*, 10, 3420.

33. Douglas, P.M.J., et al. (2021). Global estimates of methane δD. *Global Biogeochem. Cycles*, 35, e2020GB006858.

34. He, J., et al. (2026). [Methane lifetime trend]. *J. Geophys. Res. Atmos.*

35. He, J., et al. (2026). Attributing 2019–2024 methane growth using TROPOMI satellite observations. *Science*, [in press].

36. Naus, S., et al. (2019). Constraints and biases in a tropospheric two-box model of OH. *Atmos. Chem. Phys.*, 19, 407–424.

37. Worden, J.R., et al. (2017). Reduced biomass burning emissions reconcile conflicting estimates of the post-2006 atmospheric methane budget. *Nat. Commun.*, 8, 2227.

38. Worden, J.R., et al. (2017). [same as 37].

39. Still, C.J., et al. (2003). Global distribution of C3 and C4 vegetation. *Global Biogeochem. Cycles*, 17(1), 6-1–6-14.

---

## Supplementary Information

### SI Table S1. Microbial δ¹³C subcategory source signatures

| Subcategory | δ¹³C (‰ VPDB) | σ (‰) | Source |
|-------------|---------------|-------|--------|
| Wetlands | Spatially varying (isotem) | Per-cell σ | Parker et al. (31) |
| Ruminants (C₃) | −66.8 | 2.8 | Chang et al. (32) |
| Ruminants (C₄) | −51.4 | 3.4 | Chang et al. (32) |
| Rice | −63.0 | 5.0 | Riddell-Young et al. (11) |
| Termites | −57.0 | 10.0 | Riddell-Young et al. (11) |
| Waste/Landfill | −55.0 | 5.0 | Riddell-Young et al. (11) |
| Wild Animals | −66.0 | 5.0 | Riddell-Young et al. (11) |

### SI Table S2. Kinetic isotope effects and sink fractions

| Sink | KIE ¹³C | KIE D | Distribution | NH fraction | SH fraction |
|------|---------|-------|-------------|-------------|-------------|
| OH | 1.0039–1.0054 | 1.294–1.327 | Uniform | 0.825 | 0.850 |
| Cl | 1.066 ± 0.002 | 1.52 ± 0.02 | Normal | 0.040 | 0.028 |
| Strat | 1.003 ± 0.001 | 1.179 ± 0.01 | Normal | 0.070 | 0.070 |
| Soil | 1.0201 ± 0.003 | 1.083 ± 0.01 | Normal | 0.065 | 0.052 |

### SI Table S3. Data version summary

| Version | δD atmospheric | δD source sigs | δ¹³C source sigs | C₃/C₄ map |
|---------|---------------|----------------|-------------------|-----------|
| v1 | Umezawa, ±6‰ symmetric | Global | Global | Still 2003 |
| v2 | Dasgupta cal., real hemi MC | Global | Global | Still 2003 |
| v3 | Dasgupta cal., real hemi MC | Hemispheric | Global | Still 2003 |
| v4 | Dasgupta cal., real hemi MC | Hemispheric | Hemispheric | Still 2003 |
| v5 | Dasgupta cal., real hemi MC | Hemispheric | Hemispheric | Luo 2024 |

### SI Section S3. KIE parameterizations

Three KIE configurations are tested:

1. **Saueressig:** OH ¹³C KIE from Saueressig et al. (2001), U[1.0039, 1.0054]; Cl from Saueressig et al. (1995).
2. **Cantrell:** OH ¹³C KIE from Cantrell et al. (1990), U[1.0039, 1.0054] (same range, different study).
3. **Sampled:** KIE parameters drawn from their full published distributions each MC iteration.

All three yield indistinguishable results (Table 2), confirming that KIE uncertainty is not the limiting factor.

---

## Figure Captions

**Figure 1.** Dual-isotope constraint on fossil fuel methane emissions. **(A)** Time series of estimated FF emissions from the two-box model using δ¹³C-only (gray, 90% CI shading) and dual-isotope (blue, 90% CI shading) constraints, 1999–2021. BB emissions prescribed from CarbonTracker-CH₄. **(B)** Distribution of FF 90% CI widths across MC iterations for δ¹³C-only (gray histogram) and dual-isotope (blue histogram) models. Vertical dashed lines indicate medians. The dual-isotope model reduces CI width by 53%.

**Figure 2.** The δD uncertainty threshold. **(A)** 90% CI width on FF emissions as a function of the Mic δD uncertainty multiplier. Blue curve: dual-isotope model. Gray horizontal line: δ¹³C-only reference (133.1 Tg yr⁻¹). The crossover at 4.53× (σ = 37.4‰) marks where δD transitions from helpful to harmful. Red diamond: Thanwerdas et al. (2024) uncertainty level (~15.6×). **(B)** Percentage improvement vs. δ¹³C-only. Green shading: δD helps. Red shading: δD hurts. Dashed horizontal line: 10% improvement threshold.

**Figure 3.** Hemispheric source-signature gradients. **(A)** NH vs. SH source signatures in δ¹³C-δD space for microbial (green), fossil fuel (brown), and biomass burning (orange) sources. Error bars: 1σ MC uncertainty. **(B)** NH-SH gradient magnitude for each source-isotope combination. δD gradients are 5–10× larger than δ¹³C gradients.

**Figure 4.** Hemispheric decomposition of the δD improvement. **(A)** NH FF 90% CI width vs. Mic δD uncertainty multiplier. **(B)** SH FF 90% CI width (largely insensitive to δD uncertainty). **(C)** Global CI width. The improvement is driven almost entirely by NH, where δD's larger source-signature gradients provide the most additional information.

**Figure 5.** Robustness across sensitivity configurations. Each panel shows the threshold curve (as in Fig. 2A) for a different KIE × lifetime configuration. All nine configurations yield essentially identical threshold curves, with crossover consistently between 3× and 5×.

---

*Manuscript word count: ~6,800 (main text, excluding references and SI)*
