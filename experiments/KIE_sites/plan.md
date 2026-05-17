# KIE_sites — Constraining the OH KIE Ratio from Seasonal Isotope Cycles at Co-located Sites

## Research Question

Can we use the seasonal cycles of δ¹³C-CH₄ and δD-CH₄ at co-located measurement sites to observationally constrain the ratio of OH kinetic isotope effects (KIE) for ¹³C and D?

## Motivation

The OH-¹³C KIE has been contested for 25+ years:
- **Saueressig et al. (2001):** α = 1.0039
- **Cantrell et al. (1990):** α = 1.0054

This 1.5‰ difference propagates to ~13–20 Tg/yr reallocation between fossil and microbial CH₄ sources (Schwietzke et al. 2016; Basu et al. 2022) and can flip the **sign** of the post-2007 FF trend (KIE_sensitivity experiment, this repo).

Lab measurements alone haven't resolved this. An **atmospheric observational constraint** — independent of lab conditions — would be valuable.

## Core Idea

In summer, OH concentrations are 2–3× higher than winter. OH preferentially destroys the lighter isotopologue, so the remaining CH₄ becomes enriched in both heavy isotopes:
- δ¹³C shifts **more positive** (heavier)
- δD shifts **more positive** (heavier)

If OH dominates the seasonal isotope cycle, the ratio of seasonal amplitudes directly reflects the KIE ratio:

```
Δδ¹³C_seasonal / ΔδD_seasonal ≈ (α_13C − 1) / (α_D − 1)
```

With α_D ≈ 1.31 (relatively well-constrained), extracting this ratio from observations would constrain α_13C.

Predicted ratio from pure OH:
- Saueressig: (1.0039 − 1) / (1.31 − 1) ≈ 0.013
- Cantrell:   (1.0054 − 1) / (1.31 − 1) ≈ 0.017

These differ by ~30%, which may be resolvable with enough sites.

## Data Availability

### Co-located sites with both δ¹³C and δD measurements

From `sitesdata/` inventories, 12 sites have co-located measurements:

| # | Site | Location | Lat | δ¹³C range | δD range (lab) | Overlap | δD obs/yr | Suitability |
|---|------|----------|-----|-----------|----------------|---------|-----------|-------------|
| 1 | **ALT** | Alert, Canada | 82°N | 2000–2022 | 2005–2010 (INSTAAR) + 2011–2025 (MPI) | ~17 yr | 47 + 24 | ⭐ Best |
| 2 | **ZEP** | Ny-Ålesund, Svalbard | 79°N | 2001–2022 | 2008–2010 (INSTAAR) + 2006–2024 (NIPR) + 2013–2024 (IMAU) | ~16 yr | 18/50/35 | ⭐ Best |
| 3 | **BRW** | Barrow, Alaska | 71°N | 1998–2022 | 2005–2010 (INSTAAR) + 2022–2024 (IMAU) | 5 + 2 yr | 21 + 42 | Good |
| 4 | **CBA** | Cold Bay, Alaska | 55°N | 2000–2022 | 2005–2010 (INSTAAR) | 5 yr | 31 | Good |
| 5 | **MHD** | Mace Head, Ireland | 53°N | 1999–2022 | 2005–2010 (INSTAAR) | 5 yr | 21 | Good |
| 6 | **AZR** | Azores, Portugal | 39°N | 2000–2022 | 2005–2010 (INSTAAR) | 5 yr | 14 | Marginal |
| 7 | **MLO** | Mauna Loa, Hawaii | 20°N | 1998–2022 | 2005–2010 (INSTAAR) + 2023–2024 (IMAU) | 5 + 2 yr | 40 | ⭐ Best |
| 8 | **KUM** | Cape Kumukahi, Hawaii | 20°N | 1999–2022 | 2005–2010 (INSTAAR) | 5 yr | 24 | Good |
| 9 | **ASC** | Ascension Island | 8°S | 2000–2022 | 2005–2010 (INSTAAR) | 5 yr | 38 | Good |
| 10 | **SMO** | Samoa | 14°S | 1998–2022 | 2005–2010 (INSTAAR) + 2022–2024 (IMAU) | 5 + 2 yr | 29 | Good |
| 11 | **CGO** | Cape Grim, Tasmania | 41°S | 1998–2022 | 2005–2009 (INSTAAR) + 2023–2024 (IMAU) | 4 + 1 yr | 17 | OK |
| 12 | **SPO** | South Pole | 90°S | 1998–2022 | 2005–2010 (INSTAAR) | 5 yr | 29 | Poor* |

*SPO: minimal seasonal OH cycle at the pole — weak signal, likely not useful for this experiment.

### Key data properties

- **INSTAAR sites (2005–2010):** δ¹³C and δD measured from the **same flask samples** — truly paired measurements with no inter-calibration issues. This is the cleanest subset.
- **Extended records (MPI, IMAU, NIPR):** Different labs, different calibration scales. Cross-lab offsets must be characterized before combining with INSTAAR data.
- **Latitude coverage:** 82°N to 41°S — spans a wide range of OH seasonality regimes.
- **Usable sites:** ~10–11 (excluding SPO and possibly AZR due to sparse δD sampling).

## Known Complications

### 1. Source seasonality (highest priority)

Microbial sources (wetlands, rice) peak in **the same season** as OH (summer). They push δ¹³C and δD **more negative** — the opposite direction from OH enrichment. The observed seasonal amplitude is the **net** of sink enrichment and source depletion:

```
Δδ_observed = Δδ_OH_sink + Δδ_sources
```

These partially cancel. Worse, the cancellation may not be proportional for ¹³C vs D because microbial source signatures have different separations from the atmospheric mean for the two isotopes. This means the observed amplitude ratio ≠ KIE ratio unless sources are negligible or accounted for.

**Mitigations:**
- Select remote marine boundary layer sites where source seasonality is minimal (MLO, KUM, ASC, SMO)
- Use CH₄ concentration seasonal cycle as a third constraint (same source/sink balance, no fractionation)
- Fit a simple seasonal box model with both source and sink terms
- Compare tropical vs high-latitude sites: if the ratio is consistent despite very different source environments, the source contamination is likely small

### 2. Non-OH sinks

Even in summer, soil uptake (~6%, also peaks in summer), Cl (~3.5%, marine boundary layer), and stratospheric loss (~7%) contribute. The observed seasonal fractionation reflects the **bulk** KIE:

```
ε_bulk_13C = f_OH × ε_OH_13C + f_Cl × ε_Cl_13C + f_soil × ε_soil_13C + f_strat × ε_strat_13C
```

The amplitude ratio constrains ε_bulk_13C / ε_bulk_D, not the OH-specific ratio. Since OH accounts for ~84% of the total sink, the correction is modest but not negligible.

**Mitigations:**
- First report the bulk ratio — this is itself a useful observational result
- Then correct for non-OH contributions using established sink fractions and KIEs
- Compare sites with different non-OH sink environments (continental with more soil vs marine with more Cl) to see if the ratio shifts — this could help disentangle contributions

### 3. Analytical precision

| Isotope | Measurement precision (1σ) | Typical seasonal amplitude | SNR |
|---------|---------------------------|---------------------------|-----|
| δ¹³C | ~0.05‰ | ~0.2–0.5‰ | 4–10 |
| δD | ~2–3‰ | ~5–15‰ | 2–5 |

δD has worse SNR despite the larger signal. With only ~20–40 δD observations per year, the harmonic fit uncertainty on the δD amplitude may be significant. Need careful error propagation (bootstrap or MCMC on the harmonic fit).

### 4. Phase alignment as a diagnostic

If OH dominates, **both isotopes should peak at the same time** (both enriched in summer). A phase offset between δ¹³C and δD seasonal peaks is a red flag that sources or other sinks are shifting one more than the other. Reporting the phase difference site-by-site is both a quality control and a publishable result in itself.

### 5. Cross-lab calibration

For sites with δD from multiple labs (e.g., ALT has INSTAAR + MPI, ZEP has INSTAAR + IMAU + NIPR), systematic offsets between labs exist. Options:
- Analyze each lab segment independently (most conservative)
- Estimate offsets from overlapping periods and then combine
- Start with the INSTAAR-only data (same flasks as δ¹³C, cleanest pairing)

### 6. Interannual variability

With only ~5 years of overlap at most sites, the seasonal cycle stability matters. El Niño/La Niña modulate both OH concentrations and wetland emissions. Could test:
- Fit individual-year seasonal cycles and check if the ratio is stable
- Flag ENSO years and test sensitivity

### 7. Latitude-dependent seasonality

"Summer" differs by latitude:
- **NH sites:** summer = JJA, OH peaks ~July
- **SH sites:** summer = DJF, OH peaks ~January
- **Tropical sites:** OH has a weaker or double-peaked seasonal cycle near equinoxes

Don't hardcode seasons — use harmonic fitting (annual sinusoid) to extract amplitude and phase at each site regardless of hemisphere.

## Proposed Analysis Phases

### Phase 1: Data extraction and pairing
- Load event-level δ¹³C and δD data from `sitesdata/` for all 12 co-located sites
- For INSTAAR sites, pair measurements from the same flasks (match by date)
- Compute monthly means for each isotope at each site
- Basic time series plots and data quality summary

### Phase 2: Seasonal harmonic fitting
- Fit annual harmonic (A·sin(2π·t + φ) + trend) to δ¹³C and δD monthly anomalies at each site
- Extract: amplitude (A_13C, A_D), phase (φ_13C, φ_D) with bootstrap uncertainties
- Also fit CH₄ concentration seasonal cycle for reference
- Report amplitude ratio R = A_13C / A_D and phase difference Δφ = φ_13C − φ_D per site

### Phase 3: Cross-site synthesis
- Plot amplitude ratio vs latitude — is it consistent across sites?
- Plot phase difference vs latitude — where does the OH-only assumption hold?
- Identify "clean" sites (small phase offset, remote MBL) vs "contaminated" sites
- Compute weighted mean ratio from clean sites with uncertainty

### Phase 4: Source deconvolution (if needed)
- Build a simple seasonal box model: dC/dt = S(t) − L(t), where S has source seasonality and L = C/τ with OH seasonality
- Fit simultaneously to CH₄, δ¹³C, and δD seasonal cycles
- Extract the sink-only KIE ratio after accounting for source terms
- Compare with the simpler Phase 3 result

### Phase 5: KIE extraction
- From the observed bulk ratio, correct for non-OH sinks to extract OH-specific ratio
- Propagate all uncertainties (analytical, harmonic fit, sink fraction, non-OH KIE)
- Compare with Saueressig (1.0039) vs Cantrell (1.0054) predictions
- Assess whether the data can discriminate between the two

## Expected Output

```
experiments/KIE_sites/
├── plan.md                  # This file
├── analysis/
│   ├── phase1_data.py       # Data extraction and pairing
│   ├── phase2_harmonics.py  # Seasonal harmonic fitting
│   ├── phase3_synthesis.py  # Cross-site comparison
│   ├── phase4_deconv.py     # Source deconvolution (if needed)
│   └── phase5_kie.py        # KIE extraction
├── figures/
│   ├── fig1_seasonal_cycles.png       # δ¹³C and δD time series per site
│   ├── fig2_amplitude_ratio_vs_lat.png # Key result
│   ├── fig3_phase_difference.png      # QC diagnostic
│   └── fig4_kie_constraint.png        # Final KIE comparison
├── results/                 # .npz/.json output from each phase
└── RESULT.md                # Summary of findings
```

## What Success Looks Like

- **Strong result:** Amplitude ratios are consistent across ≥6 clean sites, phase offsets are small (<1 month), and the derived KIE ratio clearly favors one end of the Saueressig–Cantrell range (or falls outside it).
- **Moderate result:** Ratios are consistent at tropical MBL sites but diverge at high latitudes (implying source contamination). Tropical-only estimate constrains the ratio but with large uncertainty.
- **Negative result:** Ratios scatter widely with no coherent pattern, or phase offsets are large everywhere, meaning source seasonality dominates and the OH signal can't be isolated from seasonal cycles alone. This is still publishable — it quantifies the source-contamination problem.

## Two-Agent Workflow

Each phase is executed through a **coder → reviewer** loop:

### Roles

**Coder agent:**
- Writes the analysis script for the current phase
- Runs it and verifies it produces output without errors
- Commits code + results to the branch

**Reviewer agent:**
- Reviews the code for **clarity** (readability, comments, naming, structure, scientific correctness of implementation)
- Gives a score from 1–10
- If score **≥ 9**: phase is complete, move to the next phase
- If score **< 9**: provides specific feedback; coder revises and resubmits

### Loop

```
For each phase:
  1. Coder writes + runs the script
  2. Reviewer scores clarity (1–10) with line-specific feedback
  3. If score < 9 → Coder revises → back to step 2
  4. If score ≥ 9 → Commit, push, advance to next phase
```

### Review Criteria (clarity, scored 1–10)

| Aspect | Weight | What to look for |
|--------|--------|-----------------|
| Code readability | 25% | Clear variable names, logical flow, no magic numbers |
| Comments & docstrings | 25% | Functions documented, scientific reasoning explained inline |
| Output quality | 20% | Figures labeled, axes titled, units shown, legends present |
| Scientific correctness | 20% | Equations match plan, data handling is sound |
| Structure | 10% | Imports organized, functions factored, no copy-paste blocks |

---

## Relationship to Other Experiments

- **KIE_sensitivity:** That experiment showed the OH-¹³C KIE drives a 35.5 pp discriminant in agreement-filter analysis and determines the FF trend sign. This experiment attempts to resolve the controversy from the **observational side** rather than the model side.
- **dD_threshold:** Showed δD improves source attribution by 53% when uncertainty is below ~37‰. The seasonal analysis here uses the same δD station data but asks a different question (KIE ratio, not source partitioning).
- **KIE_immunity manuscript:** If this experiment yields a clean KIE ratio estimate, it would strengthen the manuscript's discussion of the KIE controversy.
