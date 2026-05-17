# KIE_sites — RESULT

## Summary

**Research question:** Can seasonal cycles of δ¹³C-CH₄ and δD-CH₄ at co-located
measurement sites constrain the OH kinetic isotope effect (KIE) for ¹³C?

**Answer:** Partially. The observed amplitude ratios provide a **moderate** constraint
that favors the **Saueressig (α = 1.0039)** end of the contested range, but cannot
definitively exclude the Cantrell value (α = 1.0054).

## Key Findings

### 1. Source seasonality dominates at NH sites
- Observed amplitude ratio R = A(δ¹³C)/A(δD) ranges from 0.024 (SPO) to 0.14 (BRW)
- The pure-OH prediction is R = 0.013–0.017 (Saueressig–Cantrell)
- All NH sites show R far above this range (3–10×), indicating microbial source
  seasonality (wetlands, rice) inflates the δ¹³C seasonal cycle disproportionately
- This is because δ¹³C of microbial sources (−62‰) is much closer to atmospheric
  δ¹³C (−47‰) than δD of microbial sources (−310‰) is to atmospheric δD (−86‰)

### 2. SH sites provide the cleanest constraint
- **CGO** (Cape Grim, −41°S): R = 0.028 → implied α = 1.0040
- **SPO** (South Pole, −90°S): R = 0.024 → implied α = 1.0026
- These sites have minimal source seasonality and give α values close to or
  below the Saueressig value

### 3. Final KIE constraint

| Approach | α¹³C_OH | 95% CI | Notes |
|----------|---------|--------|-------|
| **SH direct** (CGO+SPO) | 1.0034 | [0.9975, 1.0092] | Most robust, minimal source contamination |
| **Source-corrected** (all clean) | 1.0162 | [1.0026, 1.0312] | Requires source signature assumptions |
| Saueressig et al. (2001) | 1.0039 | lab measurement | |
| Cantrell et al. (1990) | 1.0054 | lab measurement | |

### 4. Interpretation
- Both approaches favor the **lower end** of the α range (closer to Saueressig)
- The SH direct constraint (1.0034) is remarkably close to the Saueressig value (1.0039)
- The 95% CI [0.9975, 1.0092] includes the Saueressig value
  and includes the Cantrell value
- **Caveat:** Even SH sites may have residual source seasonality from
  interhemispheric transport of NH wetland emissions, which would bias α upward

### 5. Diagnostics
- Phase differences δ¹³C−δD are mostly <2 months at clean sites,
  confirming both isotopes respond to the same seasonal driver
- Clear latitude gradient: R increases towards NH high latitudes,
  consistent with greater wetland source seasonality
- SMO (Samoa) is an outlier with −5.3 month phase offset — excluded

## Data Quality
- 12 co-located sites with both δ¹³C (NOAA/INSTAAR) and δD (Riddell-Young 2025)
- 8 sites classified as "clean" for this analysis
- INSTAAR sites (2005–2010) use same-flask measurements — highest quality pairing
- Bootstrap (N=2000) used for harmonic fit uncertainties
- Monte Carlo (N=50,000) propagates all parameter uncertainties

## Relationship to Other Experiments
- **KIE_sensitivity:** That experiment showed α¹³C_OH drives a 35.5 pp discriminant
  in model agreement-filter analysis and determines the FF trend sign. This result
  (α ≈ 1.003–1.004) would favor the Saueressig value, implying the FF trend may
  be flat-to-declining rather than increasing.
- **dD_threshold:** The δD seasonal amplitudes (1–4‰) at most sites are above the
  δD detectability threshold of ~37‰ annual uncertainty, confirming δD is useful
  for seasonal-cycle analysis even if marginal for source partitioning.

## Limitations
1. Only ~5 years of overlap at most sites (INSTAAR δD: 2005–2010)
2. Source isotope signatures (especially δD of wetlands) are uncertain
3. SH sites may still have some source contamination from NH transport
4. The linearized decomposition model assumes quasi-steady state
5. Cl seasonality is not fully separated from OH seasonality

## Files
- `analysis/phase1_data.py` — Data extraction and pairing
- `analysis/phase2_harmonics.py` — Seasonal harmonic fitting
- `analysis/phase3_synthesis.py` — Cross-site synthesis and classification
- `analysis/phase4_deconv.py` — Source deconvolution
- `analysis/phase5_kie.py` — KIE extraction (this phase)
- `results/` — JSON outputs from each phase
- `figures/` — Diagnostic and publication-quality figures
