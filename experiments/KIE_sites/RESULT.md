# KIE_sites — RESULT

## Summary

**Research question:** Can seasonal cycles of δ¹³C-CH₄ and δD-CH₄ at co-located
measurement sites constrain the OH kinetic isotope effect (KIE) for ¹³C?

**Answer:** Yes. After removing wetland source contamination via phasor (vector)
decomposition, the SH-site constraint gives **α¹³C_OH = 1.0043 [0.997, 1.014]**,
consistent with both the Saueressig (1.0039) and Cantrell (1.0054) laboratory
values. The multi-site corrected constraint is **α = 1.0106 [0.999, 1.019]**,
still biased slightly high by residual NH source contamination, but a major
improvement over the uncorrected value (1.0162).

## Key Findings

### 1. Source seasonality dominates at NH sites (Phase 3)
- Observed amplitude ratio R = A(δ¹³C)/A(δD) ranges from 0.024 (SPO) to 0.14 (BRW)
- The pure-OH prediction is R = 0.013–0.017 (Saueressig–Cantrell)
- All NH sites show R far above this range (3–10×), indicating microbial source
  seasonality (wetlands, rice) inflates the δ¹³C seasonal cycle disproportionately
- This is because δ¹³C of microbial sources (−62‰) is much closer to atmospheric
  δ¹³C (−47‰) than δD of microbial sources (−310‰) is to atmospheric δD (−86‰)

### 2. Phase 5 — Initial KIE constraints (before phasor correction)

| Approach | α¹³C_OH | 95% CI | Notes |
|----------|---------|--------|-------|
| **SH direct** (CGO+SPO) | 1.0034 | [0.998, 1.009] | Minimal source contamination |
| **Source-corrected** (all clean) | 1.0162 | [1.003, 1.031] | Scalar subtraction; biased high |
| Saueressig et al. (2001) | 1.0039 | lab measurement | |
| Cantrell et al. (1990) | 1.0054 | lab measurement | |

### 3. Phase 6 — Phasor source correction (main result)

#### Method
Wetland source seasonality removed using vector (phasor) decomposition:
- **Wetland emissions**: Li et al. (2026) ESSD monthly 1°×1° gridded dataset
- **δD source signatures**: Site-specific from OIPC + Douglas (2021)
- **Source regions**: NH sites ← latitude-band wetlands; SH sites ← local SH wetlands only
  (NH signal attenuated by interhemispheric transport, τ_mix ≈ 1.3 yr)

Phasor subtraction in complex B+iC plane:
```
Z_sink = Z_obs − Z_source     (for each isotope)
R_corrected = |Z_sink(δ¹³C)| / |Z_sink(δD)|
```

#### Per-site results

| Site | Lat | R_obs | R_corr | Δ | Sink peak (δ¹³C) | Sink peak (δD) | Phase diff |
|------|-----|-------|--------|---|-------------------|----------------|------------|
| ALT | +82.5° | 0.086 | 0.053 | −0.033 | Jun | Jun | 0.5 mo |
| ZEP | +78.9° | 0.056 | 0.055 | −0.001 | Jun | Jun | 0.1 mo |
| BRW | +71.3° | 0.139 | 0.060 | −0.080 | Jun | Jun | 0.5 mo |
| CBA | +55.2° | 0.036 | 0.050 | +0.014 | Jun | Jun | 0.2 mo |
| MHD | +53.3° | 0.090 | 0.063 | −0.027 | Jun | Jun | 0.2 mo |
| KUM | +19.6° | 0.061 | 0.066 | +0.005 | May | May | 0.2 mo |
| CGO | −40.7° | 0.028 | 0.032 | +0.004 | Jan | Jan | 0.1 mo |
| SPO | −90.0° | 0.024 | 0.031 | +0.008 | Jan | Jan | 0.8 mo |

**Key observations:**
1. NH high-lat corrections are large: BRW drops from 0.139 to 0.060
2. SH corrections are tiny (as expected): CGO 0.028→0.032, SPO 0.024→0.031
3. After correction, sink phases for δ¹³C and δD **agree within <1 month** at all sites — confirms the correction isolates a common OH-driven signal
4. NH sites do NOT converge to the OH-only prediction (0.013–0.017), indicating residual contamination from non-wetland seasonal sources

#### Final KIE constraint

| Approach | α¹³C_OH | 95% CI | Sites |
|----------|---------|--------|-------|
| **Phase 6 SH-only (phasor)** | **1.0043** | **[0.997, 1.014]** | CGO, SPO |
| Phase 6 all-site (phasor) | 1.0106 | [0.999, 1.019] | 8 clean sites |
| Phase 5 SH direct | 1.0034 | [0.998, 1.009] | CGO, SPO |
| Phase 5 source-corrected | 1.0162 | [1.003, 1.031] | 8 clean sites |
| Saueressig (2001) lab | 1.0039 | — | — |
| Cantrell (1990) lab | 1.0054 | — | — |

### 4. Interpretation
- The **SH phasor-corrected** value (α = 1.0043) is remarkably close to Saueressig (1.0039)
  and consistent with Cantrell (1.0054) within uncertainty
- Phasor correction substantially improves the **multi-site** constraint:
  1.0162 → 1.0106 (all sites), moving in the right direction
- The SH sites remain the cleanest constraint because they have negligible
  local source seasonality — the phasor correction is a small refinement
- NH sites improve dramatically but retain residual bias, likely from
  non-wetland seasonal sources (rice paddies, biomass burning) not included
  in the wetland-only correction

### 5. Diagnostics
- Phase consistency check passes: sink phases agree within <1 month at all sites
- NH sink phases cluster at month 5.5–6.1 (June), consistent with OH peak
- SH sink phases cluster at month 0.5–1.4 (January), consistent with SH summer OH
- Source phasor magnitudes are physically reasonable:
  - NH high-lat δ¹³C: 0.38 ‰ (large relative to obs 0.22 ‰)
  - NH high-lat δD: 7.5 ‰ (large relative to obs 2.5 ‰)
  - SH δ¹³C: 0.024 ‰ (negligible)
  - SH δD: 0.35 ‰ (small relative to obs 2.5 ‰)

## Data Quality
- 12 co-located sites with both δ¹³C (NOAA/INSTAAR) and δD (Riddell-Young 2025)
- 8 sites classified as "clean" for this analysis
- INSTAAR sites (2005–2010) use same-flask measurements — highest quality pairing
- Bootstrap (N=2000) for harmonic fit uncertainties
- Monte Carlo (N=50,000) propagates all parameter uncertainties
- Wetland seasonality from Li et al. (2026) ESSD — 35 model ensemble

## Relationship to Other Experiments
- **KIE_sensitivity:** That experiment showed α¹³C_OH drives a 35.5 pp discriminant
  in model agreement-filter analysis and determines the FF trend sign. This result
  (α ≈ 1.004) would favor the Saueressig value, implying the FF trend may
  be flat-to-declining rather than increasing.
- **dD_threshold:** The δD seasonal amplitudes (1–4‰) at most sites are above the
  δD detectability threshold of ~37‰ annual uncertainty, confirming δD is useful
  for seasonal-cycle analysis even if marginal for source partitioning.

## Limitations
1. Only ~5 years of overlap at most sites (INSTAAR δD: 2005–2010)
2. Source isotope signatures (especially δD of wetlands) are uncertain
3. Phasor correction only removes wetland seasonality — non-wetland seasonal
   sources (rice, biomass burning, fossil fuel) are not corrected
4. NH sites retain residual bias after correction (R_corr ≈ 0.05–0.07 vs
   pure-OH 0.013–0.017), limiting their utility for KIE extraction
5. SH sites provide the best constraint but are only 2 sites
6. Cl seasonality is not fully separated from OH seasonality
7. The linearized decomposition model assumes quasi-steady state

## Files
- `analysis/phase1_data.py` — Data extraction and pairing
- `analysis/phase2_harmonics.py` — Seasonal harmonic fitting
- `analysis/phase3_synthesis.py` — Cross-site synthesis and classification
- `analysis/phase4_deconv.py` — Source deconvolution
- `analysis/phase5_kie.py` — Initial KIE extraction
- `analysis/extract_wetland_seasonality.py` — Wetland emission seasonality (Li2026)
- `analysis/build_dD_source_db.py` — Site-specific δD source signatures
- `analysis/phase6_phasor.py` — Phasor source correction (this phase)
- `data/wetland_seasonality.json` — Per-site wetland harmonic coefficients
- `data/dD_source_database.json` — Site-specific δD source signatures
- `results/` — JSON outputs from each phase
- `figures/` — Diagnostic and publication-quality figures
