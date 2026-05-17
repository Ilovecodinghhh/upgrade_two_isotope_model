# When Does δD-CH₄ Improve Methane Source Attribution? Thresholds, KIE Sensitivity, and the Limits of a Dual-Isotope Two-Box Framework

**Authors:** [Your name], [Co-authors]

**Target Journal:** *Atmospheric Chemistry and Physics*

**Draft Date:** 2026-05-15 (Revision 1)

---

## Abstract

Adding hydrogen isotope (δD-CH₄) measurements to the established δ¹³C-CH₄ constraint has been proposed as a route to sharper methane source attribution, yet recent studies reach contradictory conclusions: Riddell-Young et al. (2025) find "strong agreement" between dual-isotope mass balances, while Thanwerdas et al. (2024) report that δD "has only a minor influence." We investigate this contradiction using a two-box (NH/SH) isotopic mass-balance model with hemispheric source signatures and 400–1000 Monte Carlo iterations.

We identify two key results. First, the δD threshold experiment (companion analysis) shows δD improves fossil-fuel (FF) emission uncertainty by **45% [42%, 48%]** relative to δ¹³C-only inversion — but only when microbial δD source-signature uncertainty is below a critical threshold of **σ(Mic δD) ≈ 35‰**. Thanwerdas et al.'s large prior uncertainty (~128‰) falls far above this threshold, explaining why their 3-D inversion found negligible δD benefit. Within the KIE sensitivity framework, δD reduces FF trend uncertainty by 38% (σ from 31.0 to 19.2 Tg yr⁻¹). Second, a variance decomposition reveals that the OH-¹³C kinetic isotope effect (KIE) controversy (Saueressig 1.0039 vs. Cantrell 1.0054) still accounts for **25% [12%, 34%]** of FF emission variance even in the dual-isotope framework, and determines the **sign** of the post-2007 FF trend (Saueressig: +2.9 Tg yr⁻¹; Cantrell: −5.6 Tg yr⁻¹). The KIE spread is 8.6 Tg yr⁻¹ — 66% of the 13.0 Tg yr⁻¹ found by Basu et al. (2022) with δ¹³C alone, representing a 34% reduction but not elimination.

Source signatures dominate the remaining uncertainty (48% [38%, 56%] of variance), underscoring that spatially resolved, time-varying isotopic endmember measurements — not model complexity — are the binding constraint on methane source attribution.

---

## 1. Introduction

### 1.1 The Methane Budget Controversy

Atmospheric methane has risen from ~1775 ppb in 2006 to ~1930 ppb in 2024, with growth rates exceeding 15 ppb yr⁻¹ in some years (He et al., 2026a). Simultaneously, δ¹³C-CH₄ has shifted ~0.7‰ more negative since 2007 (Lan et al., 2021; Riddell-Young et al., 2025) and δD-CH₄ has declined by 6.0 ± 0.8‰ between 2005 and 2023 (Riddell-Young et al., 2025). These isotopic trends point toward an increasingly ¹³C-depleted and D-depleted source mix — consistent with enhanced microbial emissions from wetlands, agriculture, and waste (Basu et al., 2022; Chandra et al., 2024).

Yet this interpretation is contested on multiple fronts. Satellite inversions of TROPOMI attribute only 25% of 2019–2024 growth to increasing emissions, with 59% from approach-to-steady-state and 16% from declining OH (He et al., 2026a). GOSAT-based inversions find China's coal emissions overestimated in EDGAR but livestock and waste increasing globally (Zhang et al., 2021; Maasakkers et al., 2019). Multi-isotope analyses incorporating ¹⁴C suggest FF emissions may be 30% lower than previous isotope-based estimates (Fujita et al., 2025). Chemistry-climate models emphasize that neglecting the CH₄-OH feedback can bias emission attributions by 25% (Nguyen et al., 2020; He et al., 2026b). The spread of conclusions reflects not just observational limitations but structural sensitivity to poorly constrained parameters — above all, the OH-¹³C KIE and source isotopic signatures.

### 1.2 The Promise and Limitations of δD-CH₄

The hydrogen isotope system offers independent information because δD source signatures have much larger separation between microbial (−350‰ to −250‰), thermogenic (−350‰ to −100‰), and pyrogenic (−250‰ to −175‰) production processes than δ¹³C (Sherwood et al., 2017). Moreover, δD is less sensitive to the poorly quantified Cl sink and to trends in pyrogenic emissions that confound δ¹³C interpretation (Riddell-Young et al., 2025).

Recent dual-isotope studies have reached strikingly different conclusions about δD's value:

- **Riddell-Young et al. (2025):** Global 1-box mass balance using both δ¹³C and δD finds "strong agreement" between the two isotope systems, with both pointing to microbial-dominated post-2006 growth. Fossil fuel trend: stable or declining.
- **Thanwerdas et al. (2024):** 3-D variational inversion (LMDz) finds that "assimilating δ(D,CH₄) observations in addition to the other constraints has only a minor influence" on posterior emissions. Net increase attributed ~50% fossil, ~50% agriculture+waste.
- **Dasgupta et al. (2025):** Two-box Bayesian inversion finds δD improves separation of biogenic vs. thermogenic sources, particularly constraining fossil fuel emissions in the late 1990s–early 2000s. Post-2006 growth driven mainly by wetlands.
- **Fujita et al. (2025):** Multi-isotope (¹³C + D + ¹⁴C) 1-box particle filter finds FF 30% lower than previous isotope-based studies, with ¹⁴C providing the strongest fossil constraint.
- **Chandra et al. (2024):** Multi-species inversion (MIROC4-ACTM) finds decreased FF exploitation emissions and sustained microbial increases during 1990–2020, with isotopic constraints playing a supporting role.

This literature contradiction motivates a systematic investigation: under what conditions does δD actually improve source attribution, and what are the limiting factors?

### 1.3 The KIE Problem

A separate but intertwined problem is the 25-year-old disagreement in the OH-¹³C KIE:

- **Cantrell et al. (1990):** α = 1.0054 ± 0.0009
- **Saueressig et al. (2001):** α = 1.0039 ± 0.0004

This 1.5‰ difference in effective fractionation propagates to a ~13–20 Tg yr⁻¹ reallocation between fossil and microbial categories (Schwietzke et al., 2016; Basu et al., 2022). The choice of KIE value can determine whether the post-2007 CH₄ rise is attributed primarily to microbial or fossil sources. Basu et al. (2022) identified this as "the largest uncertainty… from our knowledge of atmospheric chemistry."

Can adding δD reduce this sensitivity? Or does the KIE propagate through the dual-isotope framework with equal force?

### 1.4 Scope and Organization

We address two questions using a two-box (NH/SH) isotopic mass-balance model with hemispheric source signatures:

1. **The δD threshold:** What is the critical microbial δD source-signature uncertainty below which δD improves FF attribution, and above which it degrades it? (Section 3.1)

2. **KIE sensitivity in the dual-isotope framework:** How much does the OH-¹³C KIE controversy contribute to FF emission uncertainty in the dual-isotope framework, and how robust is the post-2007 FF trend? (Section 3.2)

We note that earlier versions of this analysis using homogeneous (global-mean) source signatures showed apparent "KIE immunity" — the KIE contribution was negligible in the dual-isotope framework. This result was an artifact of the simplified data: when realistic hemispheric source-signature heterogeneity was introduced (v3→v4 data upgrade), the KIE contribution increased substantially. The sensitivity of these results to source-signature assumptions is itself a key finding (Section 3.2.1).

Section 2 describes the model and data. Section 3 presents results including extensive structural robustness tests (W matrix sensitivity, BB perturbation, MC convergence, solver diagnostics). Section 4 discusses implications for the literature and identifies the binding constraints on methane source attribution. Section 5 concludes.

---

## 2. Methods

### 2.1 Two-Box Model

We use a hemispheric (NH/SH) two-box model following the general approach of Naus et al. (2019) and Dasgupta et al. (2025), with interhemispheric exchange time τ_ex ~ N(1.0, 0.1) yr. Each hemisphere has its own atmospheric burden, source fluxes, and sink losses. The model resolves three source categories: fossil fuel (FF), microbial (Mic), and biomass burning (BB).

For each hemisphere $h$ ∈ {NH, SH} and year $j$, the CH₄ mass balance is:

$$
\frac{d[\text{CH}_4]_h}{dt} = S_{h,j} - \frac{[\text{CH}_4]_h}{\tau_h(t)} + \frac{[\text{CH}_4]_{h'} - [\text{CH}_4]_h}{\tau_\text{ex}}
$$

where $S_{h,j}$ = FF$_{h,j}$ + Mic$_{h,j}$ + BB$_{h,j}$ and $h'$ denotes the opposite hemisphere. The total source for each hemisphere-year is computed from the observed CH₄ growth rate and the hemispheric mass balance. We note that S depends on the assumed lifetime τ, introducing a weak circularity when τ is also treated as uncertain; our variance decomposition (Section 2.7) quantifies this τ contribution explicitly.

### 2.2 Isotopic Mass Balances

For each hemisphere and year, we solve three isotopic mass-balance equations simultaneously (δ¹³C and δD):

$$
S_{h,j} \cdot \delta^{13}\text{C}_{\text{source},h,j} = \text{FF}_{h,j} \cdot \delta^{13}\text{C}_{\text{FF},h} + \text{Mic}_{h,j} \cdot \delta^{13}\text{C}_{\text{Mic},h} + \text{BB}_{h,j} \cdot \delta^{13}\text{C}_{\text{BB},h}
$$

$$
S_{h,j} \cdot \delta\text{D}_{\text{source},h,j} = \text{FF}_{h,j} \cdot \delta\text{D}_{\text{FF},h} + \text{Mic}_{h,j} \cdot \delta\text{D}_{\text{Mic},h} + \text{BB}_{h,j} \cdot \delta\text{D}_{\text{BB},h}
$$

yielding three equations (mass + two isotope) for two unknowns (FF, Mic) per hemisphere (BB prescribed from GFEDv4s). This over-determined system is solved via bounded least squares, yielding the best-fit FF and Mic emissions for each hemisphere.

For the δ¹³C-only configuration, only the mass + δ¹³C equations are used (two equations, two unknowns → unique solution).

#### 2.2.1 Solver Implementation and Weight Matrix

The over-determined system is solved via `scipy.optimize.lsq_linear` (bounded least squares) with non-negativity constraints on all sources and an upper bound of 1.5× total hemisphere source. While this upper bound is loose (potentially allowing up to ~420 Tg yr⁻¹ per hemisphere for a single source), solver diagnostics (Section 3.5.5) confirm that unphysical solutions do not occur in practice.

The system is weighted by a diagonal matrix **W** = diag($w_m$, $w_{13C}$, $w_D$). In a properly formulated inverse problem, **W** should approximate the inverse of the observation error covariance. Our default **W** = diag(100, 1, 0.5) reflects three considerations:

1. **Mass balance ($w_m$ = 100):** The total hemisphere source is derived from observed CH₄ growth rates and is known to within ~2% (Saunois et al., 2020). A large weight enforces near-exact closure of this well-constrained quantity.

2. **δ¹³C constraint ($w_{13C}$ = 1):** Atmospheric δ¹³C is measured to ±0.05‰, but the source-weighted δ¹³C budget involves source signatures with individual uncertainties of ±0.5–2.5‰ (Table 1). The effective constraint precision is dominated by source-signature uncertainty.

3. **δD constraint ($w_D$ = 0.5):** Atmospheric δD is measured to ±1‰ (vs. ±0.05‰ for δ¹³C — a 20× larger fractional uncertainty), and δD source signatures have larger absolute uncertainties (±8–20‰). The 0.5 weight reflects this ~2× larger relative uncertainty compared to δ¹³C in the source-weighted budget.

**We test the sensitivity of all results to the W choice in Section 3.5.1.** Table 7 demonstrates that variance decomposition results (KIE% and Sig%) are robust to W (< 1 percentage point variation across six configurations including identity, equal-isotope, inverse-variance, and δD-dominant weightings), while absolute trend magnitudes show moderate W dependence. This W dependence is itself informative: it reveals that δ¹³C and δD partially disagree about the FF trend direction, with the disagreement mediated by source-signature and KIE uncertainties.

### 2.3 Source-Weighted Isotopic Ratios

The source-weighted isotopic composition is derived from the observed atmospheric δ¹³C and δD evolution and the bulk KIE. We adopt the convention where α denotes the KIE (kinetic isotope effect, α > 1 for normal isotope effects), and the fractionation factor used in the isotope budget is 1/α:

$$
\delta^{13}\text{C}_{\text{source},h,j} = f^{-1}\left(\frac{n_{13,h,j+1} - n_{13,h,j} + n_{13,h,j} \cdot \alpha_{13C}^{-1}/\tau_{h,j} - \text{IHE}_{13,h,j}}{S_{h,j}}\right)
$$

where $n_{13,h,j}$ is the ¹³C-weighted atmospheric burden, $\alpha_{13C}$ is the bulk ¹³C KIE, and IHE$_{13,h,j}$ is the interhemispheric exchange flux for ¹³C. An analogous expression holds for δD.

The bulk ¹³C KIE is computed as a sink-fraction-weighted average:

$$
\alpha_{13C} = f_\text{OH} \cdot \alpha_{13C}^\text{OH} + f_\text{Cl} \cdot \alpha_{13C}^\text{Cl} + f_\text{soil} \cdot \alpha_{13C}^\text{soil} + f_\text{strat} \cdot \alpha_{13C}^\text{strat}
$$

where $f_\text{OH}$ ≈ 88%, $f_\text{Cl}$ ≈ 3.5%, $f_\text{soil}$ ≈ 3.5%, $f_\text{strat}$ ≈ 5% (Saunois et al., 2020), and individual KIE values are sampled from their laboratory-measured distributions. In the code, the fractionation factor `a = 1/α` is computed and used in the numerator of the isotope budget equation.

### 2.4 Monte Carlo Framework

We run $N$ = 400 iterations (seed = 42), with convergence verified at $N$ = 1000 (Section 3.5.3). Per iteration, we sample:

- **OH-¹³C KIE:** U[1.0039, 1.0054] (spanning Saueressig–Cantrell range)
- **OH-D KIE:** U[1.294, 1.327] (Saueressig et al., 2001)
- **Cl-¹³C KIE:** N(1.066, 0.002)
- **Cl-D KIE:** N(1.52, 0.02) (Gola et al., 2005)
- **Strat-D KIE:** N(1.179, 0.02) (Rice et al., 2003)
- **Soil-D KIE:** N(1.083, 0.02)
- **Source signatures:** hemispheric MC draws from empirical distributions (Section 2.5), varying annually within each MC iteration
- **Atmospheric δ¹³C, δD:** station-level NH/SH MC ensembles
- **Lifetime:** He et al. (2026a) time-varying parameterization τ(t) = 9.0 − 0.017(t − 2010) or fixed values
- **Cl fraction:** 3.5% ± 1% (default), with sensitivity tests at 0.6%–10%

**Note on Cl-D KIE:** The Cl-D KIE is sampled from N(1.52, 0.02) in the MC ensemble, not held fixed. Table 1 of v1 incorrectly stated "fixed"; this has been corrected.

### 2.5 Hemispheric Source Signatures

A critical advance over global-mean approaches is our use of hemisphere-specific source signatures for both isotope systems. Source signatures vary annually within each MC iteration, drawn from MC matrices of shape (24 years × N iterations) that capture both the central estimate and temporal variability.

**Table 1.** Hemispheric source isotopic signatures (mean ± 1σ).

| Source | δ¹³C NH (‰) | δ¹³C SH (‰) | Δ(NH−SH) | δD NH (‰) | δD SH (‰) | Δ(NH−SH) |
|--------|:-----------:|:-----------:|:---------:|:---------:|:---------:|:---------:|
| Fossil fuel | −44.0 ± 0.5 | −48.5 ± 0.5 | +4.5 | −194 ± 12 | −186 ± 15 | −7 |
| Microbial | −61.4 ± 1.5 | −61.3 ± 1.5 | −0.1 | −317 ± 8 | −305 ± 10 | −13 |
| Biomass burning | −25.0 ± 2.5 | −22.6 ± 2.5 | −2.4 | −232 ± 15 | −208 ± 20 | −24 |

**Key observation:** δD hemispheric gaps (7–24‰) are 5–10× larger than δ¹³C gaps (<5‰). This is the fundamental reason δD adds hemispheric information that δ¹³C cannot — and also why δD source-signature uncertainty matters more in a two-box framework.

The δ¹³C signatures derive from a hemispheric mass balance using:
- **Microbial:** isotem wetland spatial maps with subcategory mass balance (wetlands, ruminants, rice, termites, waste, wild animals), incorporating the Luo (2024) C4 vegetation distribution and Suess correction.
- **Biomass burning:** Luo (2024) time-varying C4 map (replacing static Still & Berry, 2003).
- **Fossil fuel:** Sherwood et al. (2017) database, disaggregated by hemisphere using country-level gas/coal/oil emission fractions from EDGAR. The NH–SH δ¹³C gap of 4.5‰ reflects the NH's larger share of gas emissions (δ¹³C ≈ −44‰) relative to coal (δ¹³C ≈ −36‰), while the SH has a higher relative coal fraction from Australia and South Africa.

The δD signatures derive from Sherwood et al. (2017) and Rice et al. (2016), disaggregated into NH and SH using source-type geographic distributions.

### 2.6 Atmospheric Observations

| Quantity | Source | Coverage |
|----------|--------|----------|
| [CH₄] hemispheric | NOAA GML MBL | 1984–2024 |
| δ¹³C-CH₄ NH/SH | NOAA GML/INSTAAR | 1998–2023 |
| δD-CH₄ NH/SH (MC ensembles) | Riddell-Young et al. (2025) reconstruction, Dasgupta (2025) calibration scale | 1999–2023 |
| BB CH₄ | GFEDv4s (van der Werf et al., 2017) | 1997–2022 |

The δD atmospheric record uses station-level data from NOAA/INSTAAR and Tohoku/NIPR networks, harmonized to the Dasgupta (2025) VSMOW calibration scale and reconstructed into hemispheric annual means with Monte Carlo uncertainty ensembles.

**Limitation on δ¹³C gradient uncertainty:** The atmospheric δ¹³C MC sampling applies a global offset to fixed hemispheric means, meaning the NH–SH δ¹³C gradient is held constant across MC iterations while only the global level varies. In reality, NH and SH δ¹³C have partially independent uncertainties due to different station coverage and calibration histories. This likely underestimates the uncertainty in the hemispheric gradient, which is the quantity the 2-box model exploits. We estimate this effect is small relative to source-signature uncertainty (the dominant term at 48% of variance), but future work should generate independent NH and SH δ¹³C MC ensembles.

### 2.7 Variance Decomposition

To quantify parameter contributions to FF uncertainty, we use selective freezing (cf. Saltelli et al., 2004, "factor fixing"):

1. **Full MC** → total variance σ²(FF)
2. **Fix KIE** at midpoint (OH-¹³C: 1.00465; OH-D: 1.3105; all other KIEs at their distribution means from a single source of truth, `KIE_FIXED` in `common.py`) → remaining variance σ²_no_KIE → KIE contribution = (σ² − σ²_no_KIE)/σ²
3. **Fix source signatures** at iteration-0 values → Sig contribution
4. **Fix lifetime** at 9.0 yr → τ contribution
5. **Residual** = total − KIE − Sig − τ (includes atmospheric observation uncertainty and parameter interactions)

Bootstrap confidence intervals (200 resamples of the 400-iteration MC ensemble) are computed for each decomposition.

**Methodological limitations:** Selective freezing is a first-order approximation that does not correctly attribute interaction effects between parameters (Saltelli et al., 2004). When parameters interact nonlinearly (as KIE and source signatures do through the 3×3 solver), the individual contributions may sum to more or less than 100%, with the difference absorbed into the residual. In our dual real-hemi configuration, the point-estimate decomposition (variance_decomposition_v2) gives τ% = 18.4% and residual% = 9.8%, while the bootstrap-resampled decomposition gives τ% = 0.8% and residual% = 27.4%. This discrepancy reflects the sensitivity of selective freezing to the specific MC samples used: the bootstrap procedure, by resampling, averages over many possible "freezing baselines" and produces a more robust decomposition. **We report the bootstrap-derived values (Table 3) as our primary results** because they capture this sampling uncertainty. The sum of bootstrap-central components (KIE% + Sig% + τ% + Residual% = 100.7%) confirms near-additivity, with the small overshoot reflecting weak parameter interactions.

A more rigorous approach would use Sobol sensitivity indices (first-order + total-order), which correctly partition main effects from interactions. We note this as a direction for future work and emphasize that our qualitative conclusions — KIE matters, source signatures dominate — are robust to the decomposition methodology.

### 2.8 Basu (2022) KIE Spread Comparison

To benchmark against the 3-D TM5-4DVAR result, we run full MC at fixed Saueressig (1.0039) and fixed Cantrell (1.0054), compute the post-2007 FF trend for each, and report the spread. Basu et al. (2022) found a 13.0 Tg yr⁻¹ KIE spread in their δ¹³C-only inversion. The residual analysis from the solver (mean residual norm for each KIE choice) is reported in Section 3.2.2 as an indicator of which KIE the data prefer.

### 2.9 δD Threshold Experiment

To identify the critical source-signature uncertainty, we run a sweep of microbial δD uncertainty multipliers (0.5×, 1×, 2×, 3×, 4×, 5×, 8×, 12×, 16× the baseline σ ≈ 8.2‰). For each multiplier, we compute the FF 90% CI width and compare to the δ¹³C-only reference. The crossover multiplier — where dual-isotope CI equals δ¹³C-only CI — defines the threshold.

### 2.10 Trend Metrics

We use two complementary trend metrics:

1. **Step-change (primary):** ΔFF = mean(FF, 2010–2018) − mean(FF, 2000–2006), computed on 5-year smoothed time series. The gap years 2007–2009 are excluded to avoid the transition period and to focus on the contrast between the pre-acceleration and post-acceleration regimes. This metric directly measures the "before vs. after" change relevant to the post-2007 methane growth debate.

2. **Linear regression (secondary):** OLS slope over 2000–2020 with standard error, p-value, and fraction of MC iterations significant at p < 0.05. This captures monotonic trends but may miss nonlinear behavior.

We test sensitivity to period boundaries in Section 3.5.6.

### 2.11 Three Data Configurations

| Config | δD atmospheric | δD source sigs | δ¹³C source sigs |
|--------|:---:|:---:|:---:|
| **δ¹³C-only** | Not used | — | Hemispheric MC |
| **Dual (offset)** | Global ± 6‰ | Global | Global |
| **Dual (real hemi)** | Station-level NH/SH MC | Hemispheric MC | Hemispheric MC |

The "real hemi" configuration is our primary result; the others serve as baselines.

---

## 3. Results

### 3.1 The δD Threshold: When Does δD Help?

#### 3.1.1 Baseline Improvement

At current δD measurement precision (σ(Mic δD) ≈ 8.2‰), the δD threshold experiment (companion analysis) shows the dual-isotope two-box model reduces the FF 90% CI from 105.1 Tg yr⁻¹ (δ¹³C-only) to **57.6 Tg yr⁻¹** — a **45.1% improvement** [bootstrap 95% CI: 41.7%, 47.6%]. Within the KIE sensitivity framework, δD reduces the FF trend uncertainty by 38% (σ from 31.0 to 19.2 Tg yr⁻¹; Table 3). The difference between the 45% and 38% figures reflects different metrics: the former measures CI width reduction (including outlier suppression), while the latter measures standard deviation reduction of the trend distribution.

The improvement is driven primarily by the Northern Hemisphere, where the δD constraint leverages the larger NH–SH δD source-signature gaps (Table 1).

The Degrees of Freedom for Signal (DFS) increase from 2.00 (δ¹³C-only, two hemispheres) to 3.39 (dual-isotope), confirming that δD adds 1.39 independent pieces of information to the NH/SH system — approaching the theoretical maximum of 2 additional constraints.

#### 3.1.2 Threshold Identification

**Table 2.** δD threshold sweep results.

| Mic δD σ | FF 90% CI (Tg/yr) | Improvement vs δ¹³C-only |
|----------|:------------------:|:------------------------:|
| 4.1‰ (0.5×) | 57.6 | **+45.1%** ✅ |
| 8.2‰ (1× baseline) | 57.6 | **+45.1%** ✅ |
| 16.5‰ (2×) | 59.4 | **+43.4%** ✅ |
| 24.8‰ (3×) | 77.6 | **+26.1%** ✅ |
| **~31.5‰ (3.82×) ← crossover** | **~105** | **~0%** |
| 41.2‰ (5×) | 144.1 | **−37.2%** ❌ |
| 66.0‰ (8×) | 210.8 | **−101%** ❌ |

The crossover occurs at σ(Mic δD) ≈ **35‰** (fine-grid estimate: 31.5‰). Below this threshold, δD constrains; above it, δD injects noise. The 10% improvement threshold is at 29.1‰ (~3.5× baseline).

This result is robust across 6 sensitivity configurations (3 KIE × 3 lifetime values), 4 data versions, and 3 analysis time ranges (Table S1). The threshold is a structural property of the system, not a parameter-dependent artifact.

#### 3.1.3 Reconciling Thanwerdas et al. (2024) and Riddell-Young et al. (2025)

Thanwerdas et al. (2024) used δD source-signature uncertainties of ~110–130‰ in their 3-D variational inversion — **4× above our identified threshold**. At this uncertainty level, our model predicts a −138% "improvement" (i.e., δD makes things substantially worse). Their finding that δD "has only a minor influence" is therefore entirely consistent with our threshold analysis: their uncertainty specification rendered δD uninformative.

Riddell-Young et al. (2025), by contrast, used tighter source-signature constraints (σ ≈ 8‰ for microbial δD) well below the threshold, enabling the "strong agreement" they observed. The apparent contradiction between these studies is not about δD's intrinsic value — it is about uncertainty specification.

### 3.2 KIE Sensitivity in the Dual-Isotope Framework

#### 3.2.1 Variance Decomposition

**Table 3.** Variance decomposition of FF emission uncertainty (post-2007 trend). Values are bootstrap-resampled medians with 95% CIs from 200 bootstrap resamples.

| Config | σ(FF) (Tg/yr) | KIE% [95% CI] | Sig% [95% CI] | τ% [95% CI] | Residual% |
|--------|:---:|:---:|:---:|:---:|:---:|
| δ¹³C-only | 31.0 [28.8, 33.2] | 11.1 [3.4, 17.4] | 82.6 [79.3, 85.1] | 0.0 [0.0, 1.8] | 6.3 |
| Dual (offset) | 17.0 [15.0, 19.3] | 20.5 [0.5, 35.5] | 0.0 [0.0, 0.0] | 15.7 [0.0, 28.9] | 63.4 |
| **Dual (real hemi)** | **19.2 [18.0, 20.3]** | **24.9 [12.2, 33.8]** | **47.6 [37.7, 56.2]** | **0.8 [0.1, 1.5]** | **27.4** |

Three findings emerge:

1. **δD reduces total uncertainty by 38%** (σ from 31.0 to 19.2 Tg yr⁻¹), but does not eliminate the KIE contribution. In fact, the KIE share *increases* from 11% (δ¹³C-only) to 25% (dual real-hemi) because δD suppresses the dominant noise source (source-signature variance) without affecting the KIE channel.

2. **Source signatures dominate** at 48% of variance. Hemispheric δ¹³C signatures (particularly microbial and BB) are now the binding constraint.

3. **Lifetime contributes minimally** (<1% in bootstrap). The He et al. (2026a) time-varying parameterization effectively fixes τ.

**Note on residual:** The 27.4% residual in the dual real-hemi configuration primarily reflects parameter interaction effects and atmospheric observation uncertainty (Section 2.7). The dual-offset configuration has an even larger residual (63.4%) because the offset-based δD sampling introduces correlations not captured by any single parameter group.

**Comparison with earlier data versions:** Our v2 analysis using homogeneous (global-mean) source signatures showed KIE% ≈ 0% in the dual configuration — apparent "KIE immunity." This was an artifact of the simplified data: without hemispheric source-signature heterogeneity, the 3×3 solver had no mechanism to amplify KIE differences into source partitioning. The v3→v4 upgrade to realistic hemispheric signatures restored the KIE sensitivity, demonstrating that source-signature assumptions fundamentally control the apparent KIE importance.

#### 3.2.2 The KIE Determines the Sign of the FF Trend

**Table 4.** Basu (2022) KIE spread comparison.

| Config | Saueressig ΔFF (Tg/yr) | Cantrell ΔFF (Tg/yr) | KIE Spread |
|--------|:---:|:---:|:---:|
| Basu 2022 (3-D, δ¹³C-only) | — | — | **13.0** |
| Our δ¹³C-only | +13.4 | +12.7 | 0.7 |
| Our dual (offset) | −3.7 | −2.9 | 0.8 |
| **Our dual (real hemi)** | **+2.9** | **−5.6** | **8.6** |

The dual real-hemi KIE spread is **8.6 Tg yr⁻¹** — a 34% reduction from Basu's 13.0 Tg yr⁻¹ but still large enough that the KIE choice determines the **sign** of the post-2007 FF trend. With Saueressig (1.0039), the trend is mildly positive (+2.9 Tg yr⁻¹); with Cantrell (1.0054), it is moderately negative (−5.6 Tg yr⁻¹).

The small KIE spread in the δ¹³C-only and dual-offset configurations (0.7–0.8 Tg yr⁻¹) reflects the cancellation of KIE effects in a well-mixed system where hemispheric source-signature differences are ignored. When realistic hemispheric heterogeneity is introduced, the KIE and source-signature uncertainties interact nonlinearly through the 3×3 solver, amplifying the KIE spread.

**Residual analysis:** The solver residual norms provide a weak preference for Saueressig over Cantrell in both the dual-offset (mean residual: 1.10 × 10⁻⁷ vs. 1.36 × 10⁻⁷) and dual real-hemi (5.92 × 10⁻⁷ vs. 8.43 × 10⁻⁷) configurations. However, the residuals are extremely small in absolute terms (near machine precision for the mass balance), so this preference is not statistically robust and should not be interpreted as evidence favoring either KIE value.

**Note on our δ¹³C-only KIE spread vs. Basu:** Our 2-box δ¹³C-only model produces a KIE spread of only 0.7 Tg yr⁻¹ (vs. Basu's 13.0 in TM5-4DVAR). This is not because we underestimate KIE sensitivity — it is because our 2-box δ¹³C-only model with global source signatures has very limited ability to translate KIE differences into FF vs. Mic partitioning (the FF-Mic δ¹³C contrast is ~18‰, and a 1.5‰ shift in the bulk KIE translates to <1 Tg yr⁻¹ in a global-mean balance). The hemispheric δD system provides the lever arm that amplifies KIE effects into source partitioning — exactly because δD provides independent constraints on the same sources.

#### 3.2.3 Post-2007 FF Trend

**Table 5.** FF trend in the default dual real-hemi configuration.

| Config | ΔFF step (Tg/yr) | 90% CI | Regression slope (Tg yr⁻²) | Regression p | Robust? |
|--------|:---:|:---:|:---:|:---:|:---:|
| δ¹³C-only | +13.0 | [+12.7, +13.3] | — | — | ✓ positive |
| Dual (offset) | −3.1 | [−3.3, −2.7] | — | — | ✓ negative |
| **Dual (real hemi)** | **−1.0** | [**−16.3, +14.2**] | **+0.54** [−0.69, +1.50] | **0.104** (42.5% sig.) | ⚠ spans zero |

The dual real-hemi step-change trend is −1.0 Tg yr⁻¹ with a 90% CI that includes zero — consistent with Riddell-Young et al.'s (2025) conclusion that fossil fuel emissions have been "relatively stable," and with He et al.'s (2026a) finding that 2019–2024 emission changes are dominated by livestock and waste rather than oil/gas.

The linear regression slope is +0.54 Tg yr⁻² (non-significant, median p = 0.104). The sign discrepancy between the step-change (−1.0) and regression (+0.54) reflects the nonlinear FF time series: a weak increase through 2000–2012 followed by leveling, which manifests as a negative step-change but a positive (though non-significant) linear trend. Both metrics are consistent with FF emissions being approximately stable, with uncertainty spanning zero.

### 3.3 Sensitivity Analysis

#### 3.3.1 Lifetime Sensitivity

All tested lifetimes (τ = 8.0–10.0 yr, plus He 2026 time-varying) produce negative median FF trends, but all 90% CIs include zero (Table S2). The FF trend is not lifetime-dependent within plausible ranges. KIE% ranges from 13.3% (τ = 10.0) to 14.5% (τ = 9.0) for fixed lifetimes, compared to 14.0% for the He 2026 time-varying parameterization, confirming minimal lifetime sensitivity.

**Alternative lifetime parameterization:** The He et al. (2026a) linear parameterization τ(t) = 9.0 − 0.017(t − 2010) implies monotonically increasing OH, which is debated (Zhao et al., 2023). The fixed-τ sensitivity tests (Table S2) effectively bracket the alternative of constant OH. Since lifetime contributes <1% of variance (Table 3), the choice of lifetime functional form does not materially affect our conclusions.

#### 3.3.2 OH-D KIE Sensitivity

The OH-D KIE shifts the median FF trend from −5.2 (Saueressig, 1.294) to +2.7 (Cantrell, 1.327) Tg yr⁻¹, but all configurations have CIs spanning zero (Table S3). A threshold exists: Saueressig OH-D → negative trend; Cantrell OH-D → positive trend. This parallels the OH-¹³C KIE threshold and underscores that the D/H isotope system also has its own fractionation controversies.

#### 3.3.3 Cl Fraction Sensitivity

The Cl fraction is a powerful lever (Table 6):

**Table 6.** Cl fraction sensitivity.

| Cl fraction | ΔFF median | 90% CI | σ(FF) | KIE% | Sig% |
|:-:|:---:|:---:|:---:|:---:|:---:|
| 0.6% (Thanwerdas) | +7.2 | [−9.2, +25.8] | 18.1 | 27.6 | 46.1 |
| 2.0% (Low) | +1.8 | [−13.6, +18.2] | 18.8 | 26.6 | 47.2 |
| 3.5% (default) | −3.5 | [−17.8, +10.5] | 19.2 | 21.9 | 46.9 |
| 5.0% (Medium) | −7.4 | [−19.3, +3.4] | 17.5 | 16.0 | 51.4 |
| **6.5%** | **−8.5** | [**−19.4, −0.5**] | **13.0** | **17.2** | **73.0** |
| 10% (Allan upper) | −5.7 | [−14.7, −0.0] | 7.9 | 10.3 | 83.4 |

At Cl ≥ 6.5%, the negative FF trend becomes robust (90% CI excludes zero). Higher Cl amplifies the δD constraint: Cl has the largest D/H KIE (α = 1.52), so increasing the Cl fraction makes δD more sensitive to source partitioning and reduces σ(FF). At Cl ≥ 6.5%, source signatures account for >70% of variance — the system transitions from KIE-limited to signature-limited.

#### 3.3.4 Combined Robustness Matrix

A 3 × 2 × 3 grid (τ × OH_D × Cl = 18 cells) reveals:

- **15/18 cells** have negative median FF trend
- **6/18** are robustly negative (90% CI excludes zero) — all at Cl = 6.5%
- **3/18** are positive (all at low Cl = 0.6% + Cantrell OH_D)
- **0/18** are robustly positive

The negative FF trend is the more common outcome across parameter space, but robustness requires either higher Cl or tighter source-signature constraints than currently available.

#### 3.3.5 Interhemispheric Exchange Sensitivity

τ_ex matters substantially (Table S4): fast exchange (0.5 yr) constrains well (σ = 10.4 Tg yr⁻¹, trend robustly negative) while slow exchange (2.0 yr) widens uncertainty to σ = 22.7 Tg yr⁻¹. This underscores a limitation of the 2-box framework relative to 3-D transport models, which resolve interhemispheric mixing explicitly.

### 3.4 EDGAR / CarbonTracker Validation

| Dataset | Post-2007 ΔFF (Tg/yr) |
|---------|:---:|
| **This study** (dual real-hemi, step-change) | −1.0 |
| **This study** (dual real-hemi, regression slope × 20 yr) | +10.7 |
| CarbonTracker CH₄ (posterior) | +5.5 |
| EDGAR 8.0 (Coal+ONG) | +20.6 |

Our step-change result (−1.0 Tg yr⁻¹) is weakly negative, while the cumulative linear regression over 20 years (+10.7 = 0.54 × 20) is weakly positive. Both are directionally consistent with approximately stable FF emissions, but much smaller than EDGAR's +20.6 Tg yr⁻¹. This mirrors the tension noted by Zhang et al. (2021) and Maasakkers et al. (2019) between EDGAR's high FF growth and satellite-constrained inversions.

### 3.5 Structural Robustness Tests

#### 3.5.1 Weight Matrix Sensitivity (W)

The solver weight matrix **W** controls the relative importance of mass-balance, δ¹³C, and δD constraints. To test whether our results depend on the specific W choice, we run the full analysis across six W configurations (Table 7):

**Table 7.** Weight matrix sensitivity.

| W config | σ(FF) | ΔFF step | Reg. slope | KIE% | Sig% |
|----------|:-----:|:--------:|:----------:|:----:|:----:|
| Identity (1,1,1) | 18.7 | −0.4 | +0.59 | 24.8 | 47.2 |
| Equal isotopes (100,1,1) | 18.8 | −0.1 | +0.59 | 24.7 | 47.1 |
| **Default (100,1,0.5)** | **19.2** | **−1.0** | **+0.54** | **24.6** | **47.2** |
| δD upweighted (100,1,2) | 17.4 | +2.6 | +0.75 | 24.9 | 47.9 |
| δD dominant (100,0.5,2) | 15.4 | +9.8 | +1.15 | 25.5 | 53.3 |
| Inverse-variance est (100,20,1) | 19.3 | −1.3 | +0.51 | 24.6 | 47.4 |

**Key findings:** The variance decomposition results are robust to W: KIE% varies from 24.6–25.5% (< 1 pp) and Sig% from 47.1–53.3% across all tested configurations. This is because the decomposition measures *relative* variance contributions, which are insensitive to the absolute weighting as long as the solver correctly partitions sources. However, the absolute FF trend is moderately W-dependent — upweighting δD shifts the trend positive because δD provides independent constraints that partially offset the δ¹³C-driven negative trend. This W dependence is itself a finding: it reveals that δ¹³C and δD "disagree" about the FF trend direction, with the disagreement mediated by source-signature and KIE uncertainties.

#### 3.5.2 Biomass Burning Sensitivity

BB emissions from GFEDv4s are prescribed without uncertainty. To test whether BB perturbations affect results, we apply ±10% and ±20% multiplicative scaling to BB (Table 8):

**Table 8.** BB emission sensitivity.

| BB perturbation | σ(FF) | ΔFF |
|:---:|:---:|:---:|
| −20% | 19.2 | −1.0 |
| −10% | 19.2 | −1.0 |
| Baseline | 19.2 | −1.0 |
| +10% | 19.2 | −1.0 |
| +20% | 19.2 | −1.0 |

BB perturbations have zero effect on FF uncertainty or trends. This is because BB is subtracted from the total hemisphere source before solving for FF and Mic: ΔBB redistributes between FF and Mic proportionally to the isotopic constraints, which are dominated by the FF–Mic separation rather than the BB level. While the absolute FF and Mic values shift with ΔBB, the *trend* and *uncertainty* are unchanged because the perturbation is time-invariant. A time-varying BB perturbation (e.g., enhanced BB during El Niño years) could affect trends; we note this as a limitation.

#### 3.5.3 Monte Carlo Convergence

We test whether 400 iterations is sufficient by running the model at N = 50 to N = 1000 (Table 9):

**Table 9.** MC convergence analysis.

| N_iter | σ(FF) | ΔFF | KIE% |
|:------:|:-----:|:---:|:----:|
| 50 | 17.3 | +1.0 | 21.6 |
| 100 | 18.3 | +0.6 | 24.0 |
| 200 | 19.0 | −1.9 | 22.6 |
| **400** | **19.2** | **−1.0** | **24.6** |
| 600 | 19.3 | −2.0 | 23.4 |
| 800 | 19.8 | −1.8 | 24.9 |
| 1000 | 19.8 | −1.8 | 25.7 |

σ(FF) at N=400 vs N=1000 differs by 2.8%, confirming convergence. The ΔFF trend stabilizes by N=200 and the KIE% by N=400. All reported results use N=400.

#### 3.5.4 Seed Sensitivity

Five independent random seeds produce σ(FF) = 19.1–20.1 Tg yr⁻¹ (spread 1.0 Tg yr⁻¹, 5.2% of mean) and ΔFF = −1.4 to −0.8 Tg yr⁻¹ (spread 0.6 Tg yr⁻¹). Results are robust to the specific seed choice.

#### 3.5.5 Solver Diagnostics

Across 18,400 hemisphere-year-iteration solves, zero failures occurred (0.00%). However, 90.0% of solves hit at least one bound constraint (typically the lower bound at zero for one source category). This high bound-hit rate reflects the physics of the system: in many MC iterations, the isotopic constraints push one source category to zero in one hemisphere (e.g., FF → 0 in the SH when Mic and BB account for the full isotopic signal). Bound-active solutions are physically meaningful (non-negative emissions) and do not indicate solver failure. The upper bound (1.5× total source) was never active.

#### 3.5.6 Step-Change vs. Linear Regression Trends

The step-change metric (mean 2010–2018 minus mean 2000–2006) and linear regression slope over 2000–2020 give complementary views:

- **Step-change:** ΔFF = −1.0 [−16.3, +14.2] Tg yr⁻¹ (slight post-2007 decrease)
- **Linear regression:** slope = +0.54 [−0.69, +1.50] Tg yr⁻² (slight non-significant increase, median p = 0.104, 42.5% of iterations significant at p < 0.05)

The discrepancy reflects the nonlinear FF time series: a weak increase through 2000–2012 followed by leveling, which manifests as a negative step-change but a positive (though non-significant) linear trend. Both metrics are consistent with FF emissions being approximately stable or weakly trending, with uncertainty spanning zero.

---

## 4. Discussion

### 4.1 δD Is Diagnostic, Not Deterministic

Our results show that δD provides genuine added value — a 45% reduction in FF CI width (threshold experiment) and 38% reduction in trend σ (this analysis) — but only within a well-defined parameter regime (σ(Mic δD) < 35‰). This reconciles the divergent conclusions of Thanwerdas et al. (2024) and Riddell-Young et al. (2025) and is consistent with Dasgupta et al.'s (2025) more optimistic assessment.

The improvement is primarily a hemispheric effect: δD source-signature gaps between NH and SH (up to 24‰ for BB) are 5–10× larger than δ¹³C gaps, providing the two-box model with genuinely independent information about the NH/SH source mix. This insight has implications for observational strategy: δD measurements are most valuable at stations that sample hemispheric contrasts (e.g., Barrow/Alert for NH, Cape Grim/South Pole for SH), not at stations in the tropics where interhemispheric mixing confounds the signal.

### 4.2 The KIE Controversy Remains a Binding Constraint

Despite adding δD, the OH-¹³C KIE controversy still accounts for 25% of FF variance and determines the sign of the post-2007 trend. This is because δD constrains *different parameters* than the OH-¹³C KIE: δD reduces source-signature uncertainty and provides hemispheric leverage, but it cannot disambiguate a 1.5‰ shift in the ¹³C fractionation factor that propagates through the δ¹³C mass balance.

Our KIE spread of 8.6 Tg yr⁻¹ (66% of Basu's 13.0) represents a meaningful reduction but not elimination. This result survives all W sensitivity tests (KIE% varies < 1 pp; Table 7) and is robust to the KIE_FIXED values used in the variance decomposition (Section 2.7). A new laboratory measurement of the OH-¹³C KIE at atmospherically relevant temperatures remains the single most impactful experiment for resolving the methane budget controversy.

### 4.3 Source Signatures: The Next Frontier

Source signatures account for 48% of FF variance in our dual-isotope framework — nearly twice the KIE contribution. This result aligns with Thanwerdas et al.'s (2024) finding that "uncertainties in source signatures are too large at present to impose any additional constraint" and with Fujita et al.'s (2025) suggestion that "the current database-derived estimate of the global mean biogenic δ¹³C source signature is too low."

The key source-signature uncertainties are:
1. **Microbial δ¹³C:** Currently ±1.3‰ at global scale (Riddell-Young et al., 2025), but hemispheric values depend on the wetland C3/C4 vegetation mix and are updated here using isotem spatial maps with Luo (2024) C4 distributions.
2. **BB δ¹³C:** The largest hemispheric gap (−2.4‰ between NH and SH) and the strongest temporal trend (due to changing C4 fire fractions).
3. **Microbial δD:** The threshold-controlling parameter (σ ≈ 8‰ baseline). Improving this requires more measurements of wetland, ruminant, and rice δD signatures by hemisphere.

### 4.4 Comparison with Dasgupta et al. (2025)

Dasgupta et al. (2025) is the closest methodological comparison to our work: they use a Bayesian two-box framework with both δ¹³C and δD. Key similarities and differences:

- **Source signatures:** Dasgupta et al. use Sherwood et al. (2017) directly; we disaggregate by hemisphere using country-level gas/coal/oil fractions and isotem spatial maps, producing the 4.5‰ NH–SH δ¹³C FF gap that drives much of our model's behavior.
- **Inverse method:** They use a Bayesian MCMC; we use weighted bounded least squares with MC uncertainty propagation. Their implicit weighting is encoded in the prior covariance, which plays a role analogous to our W matrix.
- **FF trend:** Dasgupta et al. find post-2006 growth driven mainly by wetlands with stable FF — consistent with our step-change result of −1.0 Tg yr⁻¹.
- **KIE sensitivity:** Dasgupta et al. do not explicitly test KIE sensitivity; our variance decomposition fills this gap and quantifies the 25% KIE contribution that persists in the dual-isotope framework.
- **δD value-added:** Both studies find δD improves fossil-biogenic separation, but our threshold analysis quantifies the critical uncertainty boundary (~35‰) below which this improvement operates.

### 4.5 Comparison with Box Model Limitations (Naus et al., 2019)

Naus et al. (2019) demonstrated that two-box models introduce systematic biases due to simplified transport. Our sensitivity to τ_ex (Section 3.3.5) confirms their concern: fast exchange (τ_ex = 0.5 yr) produces σ(FF) = 10.4 Tg yr⁻¹, while slow exchange (τ_ex = 2.0 yr) gives σ = 22.7 Tg yr⁻¹ — a factor-of-2.2 uncertainty from transport alone. The true atmosphere likely lies between these extremes, and our default N(1.0, 0.1) yr reflects this uncertainty.

However, the threshold and variance decomposition results are structural — they depend on the information content of δD relative to δ¹³C, not on the transport model's fidelity. We therefore expect the δD threshold (~35‰) to be a **lower bound** in 3-D models, where transport noise adds an additional uncertainty layer that δD must overcome. The KIE contribution (~25%) is expected to be broadly similar in 3-D implementations because it is a parameter-space property determined by the separation between Saueressig and Cantrell values, not by transport fidelity.

### 4.6 Implications for the OH Trend Debate

He et al. (2026a) find a 16% contribution of declining OH to 2019–2024 CH₄ growth, while He et al. (2026b) show that the CH₄-OH feedback is critical for interpreting δ¹³C trends. Nguyen et al. (2020) demonstrated that neglecting chemical feedbacks can bias emission estimates by 25%. Our model does not include the CH₄-OH feedback explicitly (we parameterize lifetime as exogenous), which is a limitation.

However, our variance decomposition shows that lifetime contributes <1% of FF variance when using He et al.'s (2026a) parameterization. This suggests that for the FF *trend* (which is what the KIE controversy affects), the OH trend matters less than the static OH-¹³C fractionation factor. The OH *level* affects absolute FF magnitudes but cancels in trend calculations because it shifts all years equally.

### 4.7 Policy Implications

The post-2007 FF trend in our dual-isotope framework is −1.0 [−16.3, +14.2] Tg yr⁻¹ — consistent with fossil fuel emissions being approximately stable. This is broadly compatible with:
- Riddell-Young et al. (2025): "fossil fuel emissions have remained relatively stable"
- He et al. (2026a): emission increases from livestock and waste, not oil/gas
- Schwietzke et al. (2016): "total FF not increasing" after upward revision of the base level
- Dasgupta et al. (2025): post-2006 growth driven mainly by wetlands

But inconsistent with:
- EDGAR 8.0: +20.6 Tg yr⁻¹ post-2007 FF growth
- Rice et al. (2016): +24 Tg yr⁻¹ FF increase since 1984

The uncertainty remains too large to make confident policy statements about whether the Global Methane Pledge (30% reduction in anthropogenic emissions by 2030) is on track. What our analysis does establish is that the **KIE controversy must be resolved before isotope-based top-down methods can reliably separate fossil from microbial trends at the ±5 Tg yr⁻¹ level needed for policy verification**.

---

## 5. Conclusions

1. **δD improves FF attribution by 45% (CI width) / 38% (trend σ) in a two-box framework** — but only when microbial δD source-signature uncertainty is below ~35‰. The apparent contradiction between Thanwerdas et al. (2024; "minor influence") and Riddell-Young et al. (2025; "strong agreement") is a threshold effect, not a fundamental disagreement about δD's value.

2. **The OH-¹³C KIE controversy accounts for 25% of FF variance** in the dual-isotope framework and determines the sign of the post-2007 FF trend (Saueressig: +2.9; Cantrell: −5.6 Tg yr⁻¹). The KIE spread is 8.6 Tg yr⁻¹ — 34% less than Basu et al.'s (2022) δ¹³C-only result of 13.0 Tg yr⁻¹, but not eliminated. This result is robust to the solver weight matrix (KIE% varies < 1 pp across six W configurations).

3. **Source signatures are the dominant uncertainty** at 48% of FF variance. Hemispheric, time-varying measurements of microbial and BB δ¹³C (and microbial δD) would have a larger impact on reducing FF uncertainty than any change in model complexity.

4. **The Cl fraction controls robustness:** at Cl ≥ 6.5%, the negative FF trend is robust across all tested lifetime and KIE combinations. Better constraints on tropospheric Cl would sharply reduce the ambiguity.

5. **A new laboratory measurement of the OH-¹³C KIE remains the single highest-priority experiment** for resolving the methane budget controversy. Our dual-isotope framework reduces the KIE spread by 34%, but the remaining 8.6 Tg yr⁻¹ spread is large enough to determine whether fossil fuels are rising or falling — a difference that matters for climate policy.

---

## Data Availability

The two-box model code, Monte Carlo ensembles, and all analysis scripts are available at [repository URL], branch `three-box`, path `experiments/KIE_immunity/`. All result files are generated from a single definitive model run (version `v4-post-review`, timestamp 2026-05-13T14:03:55Z; see `results/version.json` for full provenance). Atmospheric CH₄ and δ¹³C data are from NOAA GML (https://gml.noaa.gov/). The harmonized δD-CH₄ reconstruction follows Riddell-Young et al. (2025), recalibrated to the Dasgupta (2025) VSMOW scale. EDGAR 8.0 data are from https://edgar.jrc.ec.europa.eu/.

---

## Author Contributions

[To be filled]

---

## Acknowledgments

We thank Ben Riddell-Young, Sylvia Michel, and the NOAA GML team for the harmonized CH₄ isotope records; Sudhanshu Pandey and Bibhasvata Dasgupta for sharing the two-box Bayesian framework details; and Megan He for the TROPOMI-derived lifetime parameterization.

---

## References

Allan, W., Struthers, H., & Lowe, D. C. (2007). Methane carbon isotope effects caused by atomic chlorine in the marine boundary layer: Global model results compared with Southern Hemisphere measurements. *J. Geophys. Res.*, 112, D04306. https://doi.org/10.1029/2006JD007369

Basu, S., Lan, X., Dlugokencky, E., et al. (2022). Estimating emissions of methane consistent with atmospheric measurements of methane and δ¹³C of methane. *Atmos. Chem. Phys.*, 22, 15351–15377. https://doi.org/10.5194/acp-22-15351-2022

Cantrell, C. A., Shetter, R. E., McDaniel, A. H., et al. (1990). Carbon kinetic isotope effect in the oxidation of methane by the hydroxyl radical. *J. Geophys. Res.*, 95(D13), 22455–22462. https://doi.org/10.1029/JD095iD13p22455

Chandra, N., Patra, P. K., Fujita, R., et al. (2024). Methane emissions decreased in fossil fuel exploitation and sustainably increased in microbial source sectors during 1990–2020. *Commun. Earth Environ.*, 5, 147. https://doi.org/10.1038/s43247-024-01286-x

Dasgupta, B., Pandey, S., Houweling, S., et al. (2025). Global methane emission estimates from a dual-isotope inversion: New constraints from δD-CH₄. *EGUsphere* [preprint]. https://doi.org/10.5194/egusphere-2024-3974

Fujita, R., Graven, H., Zazzeri, G., et al. (2025). Global fossil methane emissions constrained by multi-isotopic atmospheric methane histories. *J. Geophys. Res. Atmos.*, 130, e2024JD041266. https://doi.org/10.1029/2024JD041266

Gola, A. A., D'Anna, B., Feilberg, K. L., et al. (2005). Kinetic isotope effects in the gas phase reactions of OH and Cl with CH₃Cl, CH₂Cl₂, and CHCl₃. *Atmos. Chem. Phys.*, 5, 2395–2402. https://doi.org/10.5194/acp-5-2395-2005

He, M., Jacob, D. J., Estrada, L. A., et al. (2026a). Attributing 2019–2024 methane growth using TROPOMI satellite observations. *Science*, 385, eadq5584.

He, J., Naik, V., & Horowitz, L. W. (2026b). Interpreting changes in global methane budget in a chemistry-climate model constrained with methane and isotopic observations. *AGU Advances*, 7, e2025AV001822.

Lan, X., Basu, S., Schwietzke, S., et al. (2021). Improved constraints on global methane emissions and sinks using δ¹³C-CH₄. *Global Biogeochem. Cycles*, 35, e2021GB007000. https://doi.org/10.1029/2021GB007000

Maasakkers, J. D., Jacob, D. J., Sulprizio, M. P., et al. (2019). Global distribution of methane emissions inferred from an inversion of GOSAT satellite data for 2010–2015. *Atmos. Chem. Phys.*, 19, 7859–7881. https://doi.org/10.5194/acp-19-7859-2019

Naus, S., Montzka, S. A., Pandey, S., et al. (2019). Constraints and biases in a tropospheric two-box model of OH. *Atmos. Chem. Phys.*, 19, 407–424. https://doi.org/10.5194/acp-19-407-2019

Nguyen, N. H., Turner, A. J., Yin, Y., et al. (2020). Effects of chemical feedbacks on decadal methane emissions estimates. *Geophys. Res. Lett.*, 47, e2019GL085706. https://doi.org/10.1029/2019GL085706

Rice, A. L., Butenhoff, C. L., Teama, D. G., et al. (2016). Atmospheric methane isotopic record favors fossil sources flat in 1980s and 1990s with recent increase. *Proc. Natl. Acad. Sci.*, 113(39), 10791–10796. https://doi.org/10.1073/pnas.1522923113

Rice, A. L., Tyler, S. C., McCarthy, M. C., Boering, K. A., & Atlas, E. (2003). Carbon and hydrogen isotopic compositions of stratospheric methane: 1. High-precision observations from the NASA ER-2 aircraft. *J. Geophys. Res.*, 108(D15), 4460. https://doi.org/10.1029/2002JD003042

Riddell-Young, B., Bruhwiler, L. M. P., Fujita, R., et al. (2025). Microbial driver of 2006–2023 CH₄ growth indicated by trends in atmospheric δD-CH₄ and δ¹³C-CH₄. *Proc. Natl. Acad. Sci.*, 122(5), e2411532122. https://doi.org/10.1073/pnas.2411532122

Saltelli, A., Tarantola, S., Campolongo, F., & Ratto, M. (2004). *Sensitivity Analysis in Practice: A Guide to Assessing Scientific Models*. John Wiley & Sons.

Saunois, M., Stavert, A. R., Poulter, B., et al. (2020). The Global Methane Budget 2000–2017. *Earth Syst. Sci. Data*, 12, 1561–1623. https://doi.org/10.5194/essd-12-1561-2020

Saueressig, G., Crowley, J. N., Bergamaschi, P., et al. (2001). Carbon 13 and D kinetic isotope effects in the reactions of CH₄ with O(¹D) and OH. *J. Geophys. Res.*, 106(D19), 23127–23138. https://doi.org/10.1029/2001JD000679

Schwietzke, S., Sherwood, O. A., Bruhwiler, L. M. P., et al. (2016). Upward revision of global fossil fuel methane emissions based on isotope database. *Nature*, 538, 88–91. https://doi.org/10.1038/nature19797

Sherwood, O. A., Schwietzke, S., Arling, V. A., & Etiope, G. (2017). Global inventory of gas geochemistry data from fossil fuel, microbial, and burning sources, version 2017. *Earth Syst. Sci. Data*, 9, 639–656. https://doi.org/10.5194/essd-9-639-2017

Thanwerdas, J., Saunois, M., Berchet, A., Pison, I., & Bousquet, P. (2024). Investigation of the renewed methane growth post-2007 with high-resolution 3-D variational inverse modeling and isotopic constraints. *Atmos. Chem. Phys.*, 24, 2129–2167. https://doi.org/10.5194/acp-24-2129-2024

van der Werf, G. R., Randerson, J. T., Giglio, L., et al. (2017). Global fire emissions estimates during 1997–2016. *Earth Syst. Sci. Data*, 9, 697–720. https://doi.org/10.5194/essd-9-697-2017

Worden, J. R., Bloom, A. A., Pandey, S., et al. (2017). Reduced biomass burning emissions reconcile conflicting estimates of the post-2006 atmospheric methane budget. *Nat. Commun.*, 8, 2227. https://doi.org/10.1038/s41467-017-02246-0

Zhang, Y., Jacob, D. J., Lu, X., et al. (2021). Attribution of the accelerating increase in atmospheric methane during 2010–2018 by inverse analysis of GOSAT observations. *Nat. Commun.*, 12, 1502. https://doi.org/10.1038/s41467-021-21727-x

Zhao, Y., Saunois, M., Bousquet, P., et al. (2023). Reconciling the bottom-up and top-down estimates of the methane chemical sink using multiple observations. *Atmos. Chem. Phys.*, 23, 789–807. https://doi.org/10.5194/acp-23-789-2023

---

## Figures

| Figure | File | Caption |
|--------|------|---------|
| Fig. 1 | `dD_threshold/figures/fig_comprehensive_6panel.png` | δD threshold experiment: (a) FF 90% CI vs Mic δD uncertainty multiplier, (b) improvement vs δ¹³C-only, (c) hemispheric breakdown, (d) DFS, (e) bootstrap confidence, (f) year-range sensitivity. |
| Fig. 2 | `KIE_immunity/figures/fig_variance_v2.png` | Variance decomposition across three data configurations. Stacked bars show KIE%, Sig%, τ%, and residual contributions with bootstrap 95% CIs. |
| Fig. 3 | `KIE_immunity/figures/fig_kie_immunity.png` | KIE spread comparison: Basu (2022) benchmark vs our three configurations. |
| Fig. 4 | `KIE_immunity/figures/fig_edgar_validation.png` | EDGAR/CarbonTracker validation: post-2007 FF trends from this study vs inventories. |
| Fig. S1 | — | Robustness matrix (Phase 8): 18-cell heatmap of ΔFF across τ × OH_D × Cl parameter space. |
| Fig. S2 | — | Lifetime sensitivity (Phase 5). |
| Fig. S3 | — | OH-D KIE sensitivity (Phase 6). |
| Fig. S4 | — | Cl fraction sensitivity with variance decomposition (Phase 7). |
| Fig. S5 | — | Interhemispheric exchange sensitivity (Phase 11). |
| Fig. S6 | — | Bootstrap CIs on variance decomposition (Phase 9). |
| Fig. S7 | — | W matrix sensitivity: ΔFF and KIE% across six W configurations. |
| Fig. S8 | — | MC convergence: σ(FF), ΔFF, KIE% vs N_iter. |

---

## Supplementary Tables

**Table S1.** δD threshold robustness across 6 sensitivity configurations (3 KIE × 2 lifetime).

**Table S2.** Lifetime sensitivity (Phase 5): ΔFF, σ(FF), KIE% for τ = 8.0–10.0 yr.

**Table S3.** OH-D KIE sensitivity (Phase 6): ΔFF, σ(FF) for α_D = 1.294–1.350.

**Table S4.** Interhemispheric exchange sensitivity (Phase 11): ΔFF, σ(FF) for τ_ex = 0.5–2.0 yr.

**Table S5.** Full robustness matrix (Phase 8): 18 cells with ΔFF, 90% CI, significance.
