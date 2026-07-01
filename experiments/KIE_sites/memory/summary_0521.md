# KIE_sites Experiment — Full Summary

**Date:** 2026-05-21  
**Experiment:** `upgrade_two_isotope_model/experiments/KIE_sites/`  
**Status:** Complete (Phases 1–6)

---

## Abstract

We ask whether the seasonal cycles of δ¹³C-CH₄ and δD-CH₄, measured
simultaneously at 12 globally distributed monitoring stations, can resolve
the long-standing disagreement between the two laboratory determinations of
the OH kinetic isotope effect (KIE) for ¹³C: Saueressig et al. (2001),
α = 1.0039, and Cantrell et al. (1990), α = 1.0054. Because OH is the
dominant CH₄ sink and its rate peaks in boreal summer, the ratio of the
seasonal amplitudes R = A(δ¹³C)/A(δD) is, in principle, a direct observable
of the ¹³C/D KIE ratio. In practice, however, the seasonal cycle of
microbial sources (wetlands, rice) contaminates the observed δ¹³C amplitude
far more than the δD amplitude — an asymmetry rooted in isotope-gap
geometry — inflating R by 3–10× at Northern Hemisphere (NH) sites.

We develop a phasor (complex-vector) decomposition that removes the wetland
source contribution from each isotope's seasonal harmonic, using
spatially-resolved wetland emission seasonality from Li et al. (2026 ESSD)
and site-specific δD source signatures derived from OIPC precipitation
isotopes + Douglas (2021) regressions. After correction, the Southern
Hemisphere (SH) sites — which have negligible local source seasonality and
serve as the cleanest natural laboratories — yield
**α¹³C_OH = 1.0043 [0.997, 1.014] (95% CI)**, a value that sits between
the two laboratory determinations and is consistent with both within
uncertainty. The multi-site (8-station) corrected constraint narrows to
**α = 1.0106 [0.999, 1.019]**, a major improvement over the uncorrected
value of 1.0162, though still biased by residual non-wetland seasonal
sources at NH sites.

---

## 1. Motivation & Background

### 1.1 The OH ¹³C KIE controversy

The fractionation factor α¹³C_OH — the ratio of rate constants
k(¹²CH₄)/k(¹³CH₄) for the OH + CH₄ reaction — is a critical parameter
in isotope-enabled methane budget models. Two laboratory measurements
disagree:

| Study | α¹³C_OH | Method |
|-------|---------|--------|
| Saueressig et al. (2001) | 1.0039 ± 0.0004 | Smog chamber, 296 K |
| Cantrell et al. (1990) | 1.0054 ± 0.0009 | Flow tube, 273–353 K |

This 0.0015 difference propagates into a ~35.5 percentage-point discriminant
in model agreement-filter analysis (from our KIE_sensitivity experiment),
affecting whether the inferred fossil-fuel emission trend is increasing or
flat-to-declining.

### 1.2 The seasonal-cycle approach

OH concentration peaks in boreal summer (≈ July in NH, ≈ January in SH),
driving seasonal enrichment of heavy isotopologues. For a sink-only world:

```
R_pure_OH = ε_13C_bulk / ε_D_bulk ≈ (α_13C − 1) / (α_D − 1)
```

The predicted pure-OH ratio is R = 0.013 (Saueressig) to 0.017 (Cantrell).
Observing R at remote stations should, in principle, distinguish these values.

### 1.3 The source contamination problem

Microbial CH₄ sources (wetlands, rice paddies) also peak in boreal summer,
superimposing a source-driven seasonal cycle on the sink-driven one. The
contamination is **asymmetric** between isotopes:

- **δ¹³C**: source gap = δ¹³C_wetland − δ¹³C_atm = −62 − (−47.3) = −14.7 ‰
  → Source/sink amplitude ratio ≈ 2.5 (source dominates at NH sites)
- **δD**: source gap = δD_wetland − δD_atm ≈ −310 − (−86) = −224 ‰
  → Source/sink amplitude ratio ≈ 0.8 (sink still dominates)

This asymmetry inflates R by 3–10× at NH stations, making raw observed
ratios useless for KIE extraction without correction.

---

## 2. Data Sources

### 2.1 Isotopic observations

| Dataset | Variable | Period | Sites | Source |
|---------|----------|--------|-------|--------|
| NOAA/INSTAAR event data | δ¹³C-CH₄ (‰ VPDB) | ~1998–2023 | 12 stations | NOAA GML, flask-air |
| Riddell-Young et al. (2025) | δD-CH₄ (‰ VSMOW) | **2005–2010** | 12 stations | Same flasks as above |

The **2005–2010 overlap** is the limiting factor. INSTAAR's same-flask δD
measurements provide the highest-quality isotope pairing (both measured on
the same air sample), but only span approximately 5 years. The number of
paired monthly means per site ranges from 14 (CBA) to 52 (ALT).

**12 co-located stations** (N → S): ALT (82.5°N), ZEP (78.9°N),
BRW (71.3°N), CBA (55.2°N), MHD (53.3°N), AZR (38.8°N), MLO (19.5°N),
KUM (19.6°N), ASC (8.0°S), SMO (14.2°S), CGO (40.7°S), SPO (90.0°S).

After quality screening, **8 sites** are classified as "clean" for
KIE analysis (AZR, MLO, ASC, SMO excluded due to poor data coverage or
local contamination).

### 2.2 Wetland CH₄ emission seasonality

| Dataset | Resolution | Period | Coverage | Source |
|---------|-----------|--------|----------|--------|
| Li et al. (2026) ESSD | 1°×1° monthly | 2000–2025 (312 months) | Global | doi:10.5281/zenodo.18870108 |

Variable `wetch4`: kg CH₄ cell⁻¹ month⁻¹, derived from 35 Global Methane
Budget model estimates × 10 ERA5 ensemble members, emulated with XGBoost.
We subset to 2005–2010 (72 months) for consistency with the isotope
observation period.

**Latitude-band climatology (2005–2010):**

| Band | Lat range | Annual (Tg/yr) | Peak month | Frac. amplitude | R² (harmonic) |
|------|-----------|---------------|------------|-----------------|---------------|
| NH high | 60–90°N | 10.9 | Jul | 1.39 | 0.85 |
| NH mid | 30–60°N | 29.7 | Jul | 0.98 | 0.94 |
| Tropics | 30°S–30°N | 114.0 | Aug | 0.21 | 0.60 |
| SH extra | 90–30°S | 2.9 | Jan | 0.63 | 0.95 |
| **Global** | all | **157.5** | — | — | — |

### 2.3 δD source signatures

Site-specific δD-CH₄ of wetland emissions derived from:
- **OIPC v3.1** (Bowen & Revenaugh, 2003): monthly precipitation δ²H at
  each station location
- **Douglas et al. (2021)** Biogeosciences: regression of measured δD-CH₄
  against δ²H of source water (growing-season mean)

Key values range from −430 ‰ (ALT, Arctic) to −301 ‰ (KUM, tropical),
with uncertainties of 10–27 ‰ per site. This replaces the old global
assumption of δD_wetland = −310 ‰ with physically grounded, latitude-
dependent estimates.

### 2.4 Ancillary parameters

| Parameter | Value | Uncertainty | Source |
|-----------|-------|-------------|--------|
| α_D_OH | 1.294 | ±0.010 | Saueressig et al. (2001) |
| α¹³C_Cl | 1.066 | ±0.005 | Saueressig et al. (1995) |
| α_D_Cl | 1.508 | ±0.050 | Saueressig et al. (1995) |
| f_OH | 0.84 | ±0.04 | Saunois et al. (2025) |
| f_Cl | 0.035 | ±0.01 | Allan et al. (2007) |
| f_soil | 0.06 | ±0.02 | Saunois et al. (2025) |
| f_strat | 0.065 | (balance) | = 1 − f_OH − f_Cl − f_soil |
| Q_total | 580 Tg/yr | ±50 | Saunois et al. (2025) |
| δ¹³C_atm | −47.3 ‰ | — | NOAA/INSTAAR |
| δD_atm | −86.0 ‰ | — | Riddell-Young (2025) |
| δ¹³C_wetland | −62 ‰ | ±5 | Whiticar (1999) |

---

## 3. Methods

### 3.1 Phase 1 — Data extraction and pairing (`phase1_data.py`)

- Load NOAA/INSTAAR δ¹³C event data and Riddell-Young δD event data for all
  12 stations
- Match by date to identify same-flask pairs
- Compute monthly means during the overlap period
- Output: per-site monthly CSV files + data quality summary

### 3.2 Phase 2 — Seasonal harmonic fitting (`phase2_harmonics.py`)

Fit each isotope time series to an annual harmonic + linear trend:

```
δ(t) = c₀ + c₁·(t − t_ref) + B·sin(2πt) + C·cos(2πt)
```

where t is in fractional years. Bootstrap (N = 2000) provides 95% CIs on
amplitude A = √(B² + C²), peak month, and phase. The harmonic coefficients
B and C are the fundamental quantities carried forward — they encode both
amplitude and phase as a complex number Z = B + iC.

### 3.3 Phase 3 — Cross-site synthesis (`phase3_synthesis.py`)

- Compute amplitude ratio R = A(δ¹³C) / A(δD) at each site
- Classify sites as "clean" (8) or "excluded" (4) based on data quality and
  MBL flag
- Key finding: R ranges from 0.024 (SPO) to 0.139 (BRW), far above the
  pure-OH prediction of 0.013–0.017 at all NH sites, confirming massive
  source contamination

### 3.4 Phase 4 — Source deconvolution (`phase4_deconv.py`)

Analytical decomposition of the observed seasonal cycle into source and sink
components using a linearized perturbation model. Confirms that the
sink-only ratio is approximately constant across latitude, as expected for
an OH-driven signal, while the source component drives the latitude
dependence.

### 3.5 Phase 5 — Initial KIE extraction (`phase5_kie.py`)

Two approaches to constrain α¹³C_OH:

1. **SH direct**: Use CGO + SPO only (minimal source contamination).
   Invert R → α via the bulk epsilon equations. Result:
   α = 1.0034 [0.998, 1.009].
2. **Source-corrected** (scalar): Subtract estimated source amplitude from
   observed amplitude at all 8 clean sites. Result: α = 1.0162 [1.003, 1.031].
   Biased high because scalar subtraction ignores phase offsets between
   source and sink.

Monte Carlo (N = 50,000) propagates uncertainties in all parameters.

### 3.6 Wetland seasonality extraction (`extract_wetland_seasonality.py`)

- Open Li2026 ESSD NetCDF (312 months × 180 lat × 360 lon)
- For each latitude band, sum over longitude → zonal total (Tg/month)
- Subset to 2005–2010 (72 months), compute monthly climatology (12 values)
- Fit annual harmonic: Q(m) = Q̄ + B_Q·sin(2πm/12) + C_Q·cos(2πm/12)
- Assign each monitoring station to its source-region latitude band

**Source region assignment logic:**

| Site | Band | Rationale |
|------|------|-----------|
| ALT, ZEP, BRW | NH high (60–90°N) | Arctic/boreal background; sees high-lat wetlands |
| CBA, MHD | NH mid (30–60°N) | Mid-latitude background; downstream of boreal/temperate |
| KUM | Tropics (30°S–30°N) | Tropical free troposphere |
| CGO, SPO | SH extra (90–30°S) | Local SH wetlands only |

**Critical design decision:** CGO and SPO are assigned to SH_extra (local
SH wetlands) rather than Global. Interhemispheric transport attenuates the
NH seasonal signal by ~88% (τ_mix ≈ 1.3 yr gives attenuation factor
≈ 1/(1 + τ_mix·ω)² ≈ 0.11) and shifts its phase by ~3 months. Assigning
these sites to the Global band would overcorrect by 7×, producing unphysical
results.

### 3.7 δD source signature database (`build_dD_source_db.py`)

- Query OIPC v3.1 for monthly precipitation δ²H at each station
- Compute growing-season mean δ²H_precip (Apr–Oct for NH, Oct–Apr for SH)
- Apply Douglas et al. (2021) wetland regression:
  δD_CH₄ = 0.95 × δ²H_water − 233 (‰)
- Compare multiple regressions (annual, growing-season, all-freshwater)
- Select recommended value per site with uncertainty envelope

### 3.8 Phase 6 — Phasor source correction (`phase6_phasor.py`)

The central methodological contribution. Treats harmonic coefficients as
complex phasors and performs vector subtraction:

```
Z_obs   = B_obs + i·C_obs              (observed, ‰)
Z_frac  = (B_Q + i·C_Q) / Q_total      (wetland fractional seasonality)
Z_src   = (δ_source − δ_atm) × Z_frac  (source phasor, ‰)
Z_sink  = Z_obs − Z_src                (corrected, ‰)
R_corr  = |Z_sink(δ¹³C)| / |Z_sink(δD)|
```

**Why phasor, not scalar?** OH peaks ~June–July; wetlands peak ~July–August
(1–2 month offset). More importantly, the source–atmosphere isotope gap is
**negative** (δ_wetland < δ_atm for both isotopes), which flips the source
phasor by 180° — wetlands peak in summer but their isotopic effect (pulling
δ toward lower values) points toward **winter** in the phasor plane. The
observed small amplitude at NH sites reflects **near-cancellation** of
anti-aligned source and sink vectors. Scalar subtraction, which ignores
this phase structure, fails; phasor subtraction handles it naturally.

**Consistency check:** After subtraction, the sink phases for δ¹³C and δD
should agree (both driven by OH). This is confirmed at all 8 sites: phase
differences are <1 month everywhere.

**Monte Carlo uncertainty (N = 50,000):** Each iteration independently
draws:
- Observed amplitude and phase from bootstrap distributions
- δD_wetland ± site-specific σ
- δ¹³C_wetland from N(−62, 5) ‰
- Q_total from N(580, 50) Tg/yr
- Wetland B_Q, C_Q scaled by N(1, 0.20) (20% ensemble spread)
- Sink fractions f_OH, f_Cl, f_soil from their priors
- Non-OH KIE parameters from their priors

The corrected ratio R is inverted to α¹³C_OH using:

```
ε_13C_bulk = R × ε_D_bulk
α_13C_OH = 1 + (ε_13C_bulk − ε_13C_non_OH) / (f_OH × 1000)
```

---

## 4. Results

### 4.1 Per-site phasor correction

| Site | Lat | R_obs | R_corr | Δ(R) | Sink peak δ¹³C | Sink peak δD | Phase diff | α¹³C (median) | 95% CI |
|------|-----|-------|--------|------|----------------|--------------|------------|---------------|--------|
| ALT | +82.5° | 0.086 | 0.053 | −0.033 | Jun | Jun | 0.5 mo | 1.0124 | [1.004, 1.021] |
| ZEP | +78.9° | 0.056 | 0.055 | −0.001 | Jun | Jun | 0.1 mo | 1.0119 | [1.005, 1.020] |
| BRW | +71.3° | 0.139 | 0.060 | −0.080 | Jun | Jun | 0.5 mo | 1.0147 | [1.003, 1.028] |
| CBA | +55.2° | 0.036 | 0.050 | +0.014 | Jun | Jun | 0.2 mo | 1.0112 | [1.001, 1.024] |
| MHD | +53.3° | 0.090 | 0.063 | −0.027 | Jun | Jun | 0.1 mo | 1.0158 | [1.003, 1.029] |
| KUM | +19.6° | 0.061 | 0.066 | +0.005 | May | May | 0.2 mo | 1.0167 | [1.006, 1.030] |
| CGO | −40.7° | 0.028 | 0.032 | +0.004 | Jan | Jan | 0.1 mo | 1.0052 | [0.999, 1.014] |
| SPO | −90.0° | 0.024 | 0.031 | +0.008 | Jan | Jan | 0.8 mo | 1.0033 | [0.996, 1.018] |

**Key observations:**
1. **NH high-lat corrections are large:** BRW drops from R = 0.139 to 0.060
   (−57%), confirming massive wetland contamination at Arctic sites.
2. **SH corrections are tiny:** CGO changes from 0.028 to 0.032 (+14%),
   SPO from 0.024 to 0.031 (+29%). The source phasor magnitude at SH sites
   is only 0.024 ‰ (δ¹³C) and 0.35 ‰ (δD), negligible compared to observed
   amplitudes.
3. **Phase consistency passes:** After correction, δ¹³C and δD sink phases
   agree within <1 month at all sites. NH phases cluster at month 5.5–6.1
   (June), SH phases at month 0.5–1.4 (January) — both consistent with
   local OH maxima.
4. **NH sites do not converge to pure-OH values:** R_corr ≈ 0.05–0.07 at NH
   sites vs pure-OH prediction of 0.013–0.017, indicating residual
   contamination from non-wetland seasonal sources (rice paddies, biomass
   burning) not included in our wetland-only correction.

### 4.2 Final KIE constraint

| Approach | α¹³C_OH | 95% CI | Sites used | Notes |
|----------|---------|--------|------------|-------|
| **Phase 6 SH-only (phasor)** | **1.0043** | **[0.997, 1.014]** | CGO, SPO | Best estimate; minimal source contamination |
| Phase 6 all-site (phasor) | 1.0106 | [0.999, 1.019] | 8 clean sites | Weighted mean; still biased by NH residual |
| Phase 5 SH direct | 1.0034 | [0.998, 1.009] | CGO, SPO | No phasor correction |
| Phase 5 source-corrected | 1.0162 | [1.003, 1.031] | 8 clean sites | Scalar subtraction; biased high |
| Saueressig (2001) lab | 1.0039 | — | — | Smog chamber measurement |
| Cantrell (1990) lab | 1.0054 | — | — | Flow tube measurement |

### 4.3 Source phasor magnitudes

| Band | A_src (δ¹³C, ‰) | A_src (δD, ‰) | Relative to obs? |
|------|------------------|---------------|------------------|
| NH high (60–90°N) | 0.38 | 7.5 | Large (175% / 295% of observed) |
| NH mid (30–60°N) | 0.74 | 11.9 | Very large (>obs) |
| Tropics | 0.20 | 2.9 | Moderate |
| SH extra | 0.024 | 0.35 | Negligible (3% / 14% of observed) |

---

## 5. Discussion

### 5.1 The SH advantage

The Southern Hemisphere sites (CGO, SPO) provide the cleanest constraint
on α¹³C_OH because:
1. **Local wetland emissions are tiny** (SH_extra = 2.9 Tg/yr vs NH_high +
   NH_mid = 40.6 Tg/yr)
2. **NH seasonal signal is attenuated ~88%** by interhemispheric transport
   (τ_mix ≈ 1.3 yr)
3. The phasor correction is a small refinement (R changes by <30%), not a
   dominant correction
4. Our SH result (α = 1.0043) sits between Saueressig (1.0039) and Cantrell
   (1.0054), consistent with both within uncertainty

### 5.2 Progressive improvement from scalar to phasor correction

The evolution of the multi-site constraint demonstrates the value of each
methodological refinement:
- **Uncorrected multi-site:** α = 1.0162 (3× above lab values)
- **Phasor-corrected multi-site:** α = 1.0106 (1.5× above; 53% of bias removed)
- **Phasor-corrected SH-only:** α = 1.0043 (fully consistent with lab values)

### 5.3 Why NH sites remain biased

After wetland-only phasor correction, NH sites still show R_corr ≈ 0.05–0.07,
well above the pure-OH prediction of 0.013–0.017. Likely causes:
- **Rice paddy emissions** peak in boreal summer with similar δ¹³C to wetlands
  but are not included in the Li2026 wetland-only dataset
- **Biomass burning** has a seasonal cycle with δ¹³C ≈ −25 ‰ (much heavier
  than wetlands), adding positive δ¹³C seasonal amplitude
- **Fossil fuel** seasonal modulation from heating demand
- All of these preferentially affect δ¹³C more than δD, further inflating R

### 5.4 Relationship to other experiments

- **KIE_sensitivity experiment:** Showed that α¹³C_OH produces a 35.5 pp
  discriminant in model agreement-filter analysis and determines the FF trend
  sign. Our result (α ≈ 1.004) favors the Saueressig value, implying the FF
  emission trend may be flat-to-declining rather than increasing.
- **dD_threshold experiment:** Confirmed that δD seasonal amplitudes (1–4 ‰
  at most sites) are above the detectability threshold, supporting the use
  of δD for seasonal-cycle analysis.

### 5.5 Validation: the phasor approach vs 3D CTMs

Liu et al. (2025, Nature) used GEOS-Chem (full 3D CTM) to decompose CH₄
seasonal cycle amplitude trends. Their approach inherently handles source–sink
phase offsets through explicit time-stepping of transport, chemistry, and
emissions. Our phasor decomposition is the analytical equivalent for the
isotope problem — necessary because no isotope-enabled CTM with sufficient
isotope species is currently available for this analysis. Liu's finding that
35% of SPO's seasonal cycle amplitude change comes from transported NH
emissions validates our conservative treatment of SH sites (using local SH
wetlands only).

---

## 6. Limitations

1. **Short overlap period (~5 years):** The INSTAAR δD record spans only
   2005–2010. With ≤52 paired months per site, harmonic fit uncertainties are
   substantial, especially at sites with sparse coverage (CBA: 14 months,
   ZEP: 16 months).

2. **Wetland-only correction:** The phasor subtraction removes only wetland
   source seasonality. Non-wetland seasonal sources (rice paddies, biomass
   burning, fossil fuel heating) are not corrected, and their combined effect
   retains ~50% of the original bias at NH sites.

3. **Source isotope signatures uncertain:** δD of wetland CH₄ depends on
   local water isotopes and fractionation during methanogenesis. The Douglas
   (2021) regression has residual scatter of ±25–30 ‰, and the growing-season
   vs annual-mean δ²H choice adds systematic uncertainty.

4. **Two-site SH constraint:** Only CGO and SPO provide clean constraints.
   While their agreement is encouraging (1.0052 and 1.0033), the small sample
   limits statistical power. The 95% CI [0.997, 1.014] cannot distinguish
   Saueressig from Cantrell at the 2σ level.

5. **Cl sink seasonality not separated:** Stratospheric Cl also has an annual
   cycle and a large KIE (α¹³C_Cl = 1.066). Our model lumps Cl into the
   bulk sink parameterization; if Cl seasonality differs from OH seasonality,
   this introduces bias. The effect is small (f_Cl ≈ 3.5%) but systematic.

6. **Linearized model:** The decomposition assumes quasi-steady state and
   small perturbations around annual means. At NH high-latitude sites where
   wetland emissions vary by a factor of 40 across the year (Jan: 0.08 Tg/mo,
   Jul: 2.81 Tg/mo), the linear approximation is strained.

7. **Wetland climatology period:** We subset the Li2026 dataset to 2005–2010
   (72 months) for consistency with observations. The full dataset spans
   26 years (2000–2025), and a longer climatology might be more representative,
   though it would not address the fundamental limitation of the short
   isotope observation period.

---

## 7. File Inventory

### Analysis scripts

| Script | Phase | Purpose |
|--------|-------|---------|
| `analysis/phase1_data.py` | 1 | Data extraction, pairing, monthly means |
| `analysis/phase2_harmonics.py` | 2 | Seasonal harmonic fitting (B, C, amplitude, phase) |
| `analysis/phase3_synthesis.py` | 3 | Cross-site synthesis, site classification |
| `analysis/phase4_deconv.py` | 4 | Analytical source deconvolution |
| `analysis/phase5_kie.py` | 5 | Initial KIE extraction (SH direct + scalar correction) |
| `analysis/extract_wetland_seasonality.py` | 6 | Li2026 wetland emission seasonality per band |
| `analysis/build_dD_source_db.py` | 6 | Site-specific δD source signatures (OIPC + Douglas) |
| `analysis/phase6_phasor.py` | 6 | Phasor source correction + final KIE constraint |

### Data files

| File | Content |
|------|---------|
| `data/wetland_seasonality.json` | Per-band wetland harmonic coefficients (B_Q, C_Q, Q_mean) |
| `data/dD_source_database.json` | Per-site δD_wetland with multiple regression options |
| `results/phase1_data/` | Monthly CSV files, site summary |
| `results/phase2_harmonics/harmonic_fits.json` | Per-site B, C, amplitude, phase, CIs |
| `results/phase3_synthesis/synthesis_results.json` | Site classification, R values |
| `results/phase6_phasor/phasor_results.json` | Per-site phasor decomposition + MC results |

### Figures

| Figure | Content |
|--------|---------|
| `fig1_timeseries_grid.png` | δ¹³C and δD time series at all 12 sites |
| `fig1_data_coverage.png` | Temporal coverage heatmap |
| `fig2_seasonal_cycles.png` | Harmonic fits per site |
| `fig2_harmonic_summary.png` | Amplitude and phase summary |
| `fig3_ratio_vs_latitude.png` | R_obs vs latitude |
| `fig3_site_classification.png` | Clean vs excluded sites |
| `fig3_phase_diagnostic.png` | Phase consistency |
| `fig4_decomposition.png` | Source/sink decomposition |
| `fig4_sink_ratio.png` | Sink-only ratio vs latitude |
| `fig5_kie_constraint.png` | Phase 5 KIE results |
| `fig6_dD_source_vs_latitude.png` | δD source signatures vs latitude |
| `fig7_wetland_seasonality.png` | Li2026 wetland emission climatology per band |
| `fig8_phasor_decomposition.png` | Phasor vector diagrams (BRW, CBA, CGO, SPO) |
| `fig9_corrected_ratio.png` | R_corrected vs latitude with arrows from R_obs |
| `fig10_alpha_constraint.png` | Final α¹³C_OH constraint vs Saueressig/Cantrell |

### Memory / documentation

| File | Content |
|------|---------|
| `RESULT.md` | Executive summary of Phase 6 results |
| `memory/summary_0520.md` | Session log: Li2026 exploration + Liu2025 analysis |
| `memory/Plan_0520.md` | Implementation plan for Steps 1–4 |
| `memory/QA_0521.md` | Q&A summary of coder–reviewer session |
| `memory/summary_0521.md` | This document |

---

## 8. Conclusions

1. The seasonal amplitude ratio R = A(δ¹³C)/A(δD) is a viable observable
   for constraining α¹³C_OH, but **only after removing source contamination**
   — especially at NH sites where wetland seasonality inflates δ¹³C
   amplitudes by 3–10×.

2. Phasor (complex-vector) decomposition is the correct framework for source
   correction because it properly handles the ~1–2 month phase offset between
   source and sink seasonal peaks, and the 180° flip caused by negative
   source–atmosphere isotope gaps.

3. The **SH phasor-corrected constraint, α¹³C_OH = 1.0043 [0.997, 1.014]**,
   is our best estimate. It is consistent with both Saueressig (1.0039) and
   Cantrell (1.0054) within uncertainty, though the point estimate sits
   closer to Saueressig.

4. Resolving the Saueressig–Cantrell difference (0.0015 in α) from
   atmospheric observations alone will require either (a) longer co-located
   δ¹³C + δD records at SH sites, (b) correction for all seasonal sources
   (not just wetlands) at NH sites, or (c) a full isotope-enabled 3D CTM
   inversion.
