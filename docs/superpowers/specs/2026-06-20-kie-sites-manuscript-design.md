# KIE Sites Manuscript Design

## Scope

This document defines the planned structure for the first English manuscript draft from the `KIE_sites` experiment.

The approved framing is a combined constraint-and-diagnostic paper:

```text
Seasonal methane isotope phasors provide an atmospheric constraint on the OH 13C kinetic isotope effect, but the constraint is only robust when wetland source seasonality and transport-sensitive source-region assumptions are treated explicitly.
```

The main manuscript should be written as a submission-ready draft, not as an internal project note. It should focus on the phasor method, wetland source correction, the Southern Hemisphere KIE constraint, and the broad limitation that unresolved source phasors affect Northern Hemisphere sites. Biomass burning should be mentioned only briefly in the main text and documented in detail in the Supplementary Information.

The main text should lead with the strongest, most defensible result: the Southern Hemisphere phasor-corrected constraint. Southern Hemisphere source-region sensitivity should be framed as a robustness check. The tropical-wetland contribution can be highlighted as reassuringly modest; delayed Northern Hemisphere wetland residuals should be treated as a boundary-case structural sensitivity and described mainly in the Supplementary Information.

## Core Message

The paper should make three claims, in this order:

1. Co-located atmospheric `delta13C-CH4` and `deltaD-CH4` seasonal cycles can in principle constrain the OH `13C` kinetic isotope effect because OH imprints a common seasonal sink phasor on both isotope systems.

2. Raw amplitude ratios are not clean OH KIE constraints because seasonal microbial sources, especially wetlands, contribute isotope source phasors that must be removed as vectors rather than scalar amplitudes.

3. After wetland phasor correction, Southern Hemisphere sites provide the cleanest current atmospheric constraint, while Northern Hemisphere sites remain valuable diagnostics of unresolved source-seasonality structure rather than direct evidence for a very large OH KIE.

## Citation And Provenance Rule

Each quantitative or factual claim must be labeled during drafting as:

```text
CITED    = supported by external literature or data-source documentation
DERIVED  = computed from local scripts/results
UNCERTAIN = plausible but requiring verification, stronger citation, or softer wording
```

Every numerical claim in the manuscript should trace to one of:

1. A cited external source.
2. A table or figure generated from the analysis and described in Methods or Supplementary Information.
3. A reproducible analysis output listed in the code/data availability section.

Do not cite local files in the manuscript body. Local scripts and result files are internal provenance and reproducibility artifacts. The manuscript should cite the original data products and peer-reviewed literature wherever possible. Local files may be named only in the code/data availability statement or SI reproducibility table.

Use a manuscript claim ledger while drafting:

| Claim | Value | Status | Source | Manuscript location |
|---|---:|---|---|---|
| SH phasor-corrected KIE | verify before drafting | DERIVED | Main analysis table/figure generated from described method | Results |
| Wetland emissions are gridded local emissions, not station footprints | qualitative | CITED | Li et al. 2026 ESSD article and Zenodo data set | Data/Discussion |
| Effective SH source-region sensitivity is a structural robustness test | qualitative | CITED + DERIVED | Interhemispheric transport literature; SI sensitivity analysis | SI, brief Discussion |
| BB correction is secondary to the main manuscript framing | qualitative | DERIVED | SI source-sector sensitivity analysis | SI |

## Primary Source Targets

The manuscript should prioritize these original data products and primary literature sources.

| Use | Preferred primary source target | Notes |
|---|---|---|
| Atmospheric `CH4` mole fractions | NOAA/GML Global Greenhouse Gas Reference Network methane data. | Cite the NOAA/GML data product/version actually used. |
| Atmospheric `delta13C-CH4` | INSTAAR / NOAA GML stable carbon isotope data product: Michel et al. 2023, `Stable Isotopic Composition of Atmospheric Methane (13C) from the NOAA GML Carbon Cycle Cooperative Global Air Sampling Network, 1998-2022`, DOI `10.15138/9p89-1x02`. | This is the preferred current primary data citation if it covers the sites and years used. |
| Atmospheric `deltaD-CH4` | White, Vaughn, and Michel 2016, INSTAAR / NOAA ESRL `Stable Isotopic Composition of Atmospheric Methane (D/H) from the NOAA ESRL Carbon Cycle Cooperative Global Air Sampling Network, 2005-2009`, version `2016-04-26`. | Needs final access/citation wording because the source appears as a data path rather than a DOI. |
| Wetland methane emissions | Li et al. 2026 ESSD article, `Machine-learning-based estimates of global natural vegetated wetland methane emissions (2000-2025)`, and Zenodo data set DOI `10.5281/zenodo.18870108`. | Cite the article in the main text; include the data set DOI in Data Availability. |
| Wetland `deltaD-CH4` source signatures | Stell, Douglas, Rigby, and Ganesan 2021, `The impact of spatially varying wetland source signatures on the atmospheric variability of deltaD-CH4`, Phil. Trans. R. Soc. A. | Verify author order and exact bibliographic details before final draft. |
| OIPC precipitation isotope estimates | OIPC citation guidance: Bowen and Wilkinson 2002; Bowen and Revenaugh 2003; GNIP/IAEA-WMO database as underlying data source. | Use only for the OIPC-derived source-signature sensitivity. |
| OH `13C` KIE lab references | Cantrell et al. 1990, JGR, DOI `10.1029/JD095iD13p22455`; Saueressig et al. 2001, JGR, DOI `10.1029/2000JD000120`. | Present both neutrally unless the final Results justify stronger language. |
| OH `D` KIE reference | Saueressig et al. 2001, plus any more recent laboratory comparison used by the final method. | Verify final value and citation before drafting. |

Known unresolved source item: the best citation for the simplified `delta13C_wetland = -62 permil` prior still needs final verification. Possible sources include spatially resolved wetland `delta13C-CH4` source-signature studies and global methane isotope-budget inventories.

## Working Title

Primary option:

```text
Seasonal methane isotope phasors constrain the OH 13C kinetic isotope effect and reveal source-seasonality limits at northern sites
```

Shorter option:

```text
Seasonal methane isotope phasors constrain the methane OH 13C kinetic isotope effect
```

The longer title is preferred for the first draft because it signals both the positive constraint and the diagnostic caution.

## Abstract Skeleton

1. Methane isotope budgets depend on the OH kinetic isotope effect, but laboratory `13C` KIE values differ enough to affect source attribution. Cite Cantrell et al. 1990 and Saueressig et al. 2001.

2. We use co-located atmospheric `delta13C-CH4` and `deltaD-CH4` seasonal cycles to derive an observational constraint on `alpha_13C_OH`.

3. Raw isotope amplitude ratios exceed pure-OH expectations, especially at Northern Hemisphere sites, indicating strong source-seasonality contamination.

4. We correct wetland source effects using complex phasor subtraction based on primary gridded wetland emission data and published wetland isotope source signatures.

5. Southern Hemisphere sites provide the cleanest current constraint, while Northern Hemisphere sites retain residual source structure after wetland correction.

6. Seasonal isotope phasors can constrain OH KIE, while source-region phasor structure remains the main limitation for extending the inference to all sites.

## Main Text Structure

### 1. Introduction

Paragraph 1: Methane isotope KIE motivation.

- `delta13C-CH4` and `deltaD-CH4` are used to infer methane sources and sinks.
- OH is the dominant atmospheric methane sink, so OH isotope fractionation affects methane isotope budgets.
- Laboratory OH `13C` KIE values differ by enough to matter for atmospheric interpretation.
- Citation targets: Cantrell et al. 1990; Saueressig et al. 2001; recent isotope-budget literature.
- Status: CITED.

Paragraph 2: Seasonal-cycle opportunity.

- OH has a strong seasonal cycle.
- If seasonal isotope cycles were sink-dominated, the amplitude ratio

```text
R = A(delta13C) / A(deltaD)
```

would constrain the relative OH KIEs.
- Status: DERIVED + method assumption. Needs concise derivation in Methods.

Paragraph 3: Source-seasonality problem.

- Wetlands are seasonal, isotopically depleted methane sources.
- Wetland seasonality adds isotope source phasors to observed atmospheric isotope cycles.
- Source contamination is stronger for `delta13C` relative to its OH sink signal than for `deltaD`.
- Citation targets: Li et al. 2026 for wetland seasonality; Douglas et al. 2021 for wetland `deltaD`; final wetland `delta13C` source-signature citation still needs verification.
- Status: CITED + DERIVED.

Paragraph 4: Paper contribution.

- We combine harmonic fitting, wetland source phasors, and Monte Carlo uncertainty propagation.
- We derive a conditional Southern Hemisphere atmospheric KIE constraint.
- We show why Northern Hemisphere sites should be treated as diagnostics of source-seasonality structure rather than clean KIE constraints.
- Status: DERIVED from local Phase 2-14 outputs.

### 2. Data And Provenance

Subsection 2.1: Atmospheric observations.

- Co-located `CH4`, `delta13C-CH4`, and `deltaD-CH4`.
- Monthly paired observations, with emphasis on the 2005-2010 overlap period.
- Site-selection workflow: 12 candidate sites, 8 clean sites.
- Manuscript citations should point to original NOAA/GML and INSTAAR data products, not local processed files.
- The exact `deltaD-CH4` data citation/access statement needs final verification.
- Local processed files can be listed only in SI reproducibility or code/data availability.

Subsection 2.2: Wetland emission seasonality.

- Li et al. 2026 provides monthly gridded natural wetland methane emissions for 2000-2025.
- We extract 2005-2010 climatological annual phasors by latitude/source band.
- Li wetland fields are local emissions, not transported station footprints.
- Main-text citation target: Li et al. 2026 ESSD article.
- Data availability target: Zenodo data set DOI `10.5281/zenodo.18870108`.
- Local processed files can be listed only in SI reproducibility or code/data availability.

Subsection 2.3: Wetland isotope source signatures.

- `deltaD_wetland` uses Douglas et al. 2021 and OIPC-derived estimates.
- `delta13C_wetland = -62 permil` is a simplifying prior.
- Citation targets: Stell/Douglas et al. 2021 for spatially varying wetland `deltaD-CH4`; OIPC/GNIP references for precipitation isotope estimates; final wetland/global microbial `delta13C` source-signature citation.
- Local processed files can be listed only in SI reproducibility or code/data availability.

Recommended main-text table:

```text
Table 1. Data and parameter provenance.
```

Columns:

```text
Variable | Original data source | Years used | Processing step | Manuscript use | Availability
```

### 3. Methods

Subsection 3.1: Seasonal harmonic fitting.

- Fit annual harmonics to detrended monthly anomalies.
- Store each isotope seasonal cycle as a phasor:

```text
Z = B + iC
```

- Local scripts and outputs are reproducibility artifacts and should be summarized in SI, not cited as evidence in the main text.

Subsection 3.2: Observed isotope amplitude ratio.

Define:

```text
R_obs = A(delta13C) / A(deltaD)
```

Compare `R_obs` with the pure-OH expectation implied by laboratory `alpha_13C_OH` and `alpha_D_OH`.

Subsection 3.3: Wetland source phasor construction.

Define:

```text
Z_source = (delta_source - delta_atm) * Z_wetland_flux_fraction
```

State explicitly that source bands represent approximate upstream source regions, not local station grid cells and not CTM footprints.

Subsection 3.4: Vector source correction.

Core equation:

```text
Z_sink = Z_obs - Z_source
R_corrected = |Z_sink_13C| / |Z_sink_D|
```

Emphasize that scalar amplitude subtraction is not valid when source and sink seasonal phases differ.

Subsection 3.5: Conversion to OH `13C` KIE.

- Convert `R_corrected` to `alpha_13C_OH`.
- Account for non-OH methane sinks using specified sink fractions and KIEs.
- All sink fractions and non-OH KIEs need citation or explicit model-source provenance before final manuscript submission.
- All numerical values must be reported through manuscript tables/figures, with original citations for externally specified parameters.

Subsection 3.6: Sensitivity experiments.

Main-text sensitivity:

- Southern Hemisphere source-region sensitivity should be presented as a short robustness check rather than a central result.
- Emphasize that adding tropical wetland influence has a modest effect on the Southern Hemisphere KIE constraint under the tested assumptions.
- Delayed Northern Hemisphere wetland residuals should be discussed mainly in SI as a structural upper-bound sensitivity.
- Treat as diagnostic, not as a replacement main model.
- Transport parameters use effective annual phasor attenuation and phase lag.
- Main-text citation/provenance should point to interhemispheric transport literature and the SI sensitivity analysis, not local result files.

Supplementary-only sensitivity:

- Biomass burning correction.
- OSSE validation.
- Block bootstrap.
- Harmonic model sensitivity.
- Uncertainty attribution.

### 4. Results

Result 1: Raw seasonal isotope ratios exceed pure-OH expectations.

- Show latitude dependence of `R_obs`.
- Northern sites sit far above pure-OH prediction.
- Southern sites are closer.
- Main figure candidate:

```text
Figure 1. Site coverage and/or observed seasonal isotope amplitude ratios.
```

The reported values should appear in a manuscript table or figure. The underlying atmospheric observations should be cited to the original NOAA/GML and INSTAAR data products, with analysis outputs listed only in SI reproducibility or data/code availability.

Result 2: Wetland phasor correction reduces source contamination.

- Wetland source phasors are large at NH wetland-influenced sites and small for SH sites under the base `SH_extra-only` assumption.
- Vector correction changes both amplitude and phase.
- Main figure candidate:

```text
Figure 2. Wetland phasor decomposition for representative sites.
```

Internal figure candidates:

```text
experiments/KIE_sites/figures/fig8_phasor_decomposition.png
experiments/KIE_sites/figures/fig9_corrected_ratio.png
```

Result 3: SH sites provide the cleanest current atmospheric KIE constraint.

- Report the SH-only phasor-corrected estimate from Phase 6 after rechecking the JSON.
- Compare with Saueressig and Cantrell lab values.
- Lead this as the primary result of the paper.
- Main figure candidate:

```text
Figure 3. Corrected site ratios and inferred alpha_13C_OH.
```

Numerical values should be taken from the final manuscript results table after verification, not cited to local JSON in the prose.

Result 4: All-site corrected estimate remains elevated.

- Interpret as residual source structure at NH sites, not as direct evidence for very large OH `13C` KIE.
- Keep biomass burning out of main text except one sentence: additional source-sector sensitivity is evaluated in SI.
- Do not let this result dominate the Results section; it should support the decision not to over-interpret Northern Hemisphere sites.

Result 5: SH constraint is robust to modest tropical wetland influence but remains structurally conditional.

- Main text should emphasize the reassuring result: tested tropical wetland influence does not strongly disrupt the Southern Hemisphere constraint.
- Delayed Northern Hemisphere wetland residuals can move `CGO/SPO`-based `R_corrected` and `alpha_13C_OH`, but this should be framed as a structural sensitivity and mostly placed in SI.
- Main figure only if needed:

```text
Figure 4 or SI Figure. Southern Hemisphere source-region sensitivity, with tropical and delayed-NH cases visually separated.
```

If the existing figure overemphasizes delayed NH influence, create a new figure that leads with `SH_extra-only` and tropical scenarios, with delayed NH cases in muted colors or moved to SI.

### 5. Discussion

Theme 1: Atmospheric seasonal isotope phasors are useful KIE constraints.

- The method recovers a physically interpretable SH constraint.
- The method is strongest when source phasors are small or independently constrained.

Theme 2: NH sites are diagnostics, not clean KIE constraints yet.

- Residual high corrected ratios likely reflect missing or mis-specified seasonal source phasors.
- Avoid wording that implies an extremely large OH KIE is required by the NH data.

Theme 3: Transport matters for background stations.

- Li wetland data are local emissions.
- Station influence requires transport, mixing, and footprint assumptions.
- Keep this concise in the main Discussion.
- Phase 14 is an effective phasor sensitivity anchored by interhemispheric transport literature, not a CTM result; detailed transport discussion belongs in SI.

Theme 4: What must improve.

- Longer co-located isotope records.
- Better wetland isotope priors.
- CTM or footprint constraints on source-region phasors.
- Additional source-sector tracers or inventories.

### 6. Conclusions

1. Paired `delta13C-CH4` and `deltaD-CH4` seasonal cycles provide an atmospheric route to constraining `alpha_13C_OH`.

2. Wetland source seasonality must be corrected in phasor space.

3. Southern Hemisphere background sites provide the cleanest current constraint, but the result remains conditional on source-region and transport assumptions.

4. Northern Hemisphere sites reveal unresolved source-seasonality structure and should not yet be pooled as clean OH KIE constraints.

5. Future constraints require better transport-aware source phasors and longer paired isotope records.

## Main Figure Plan

Recommended first-draft figure set:

| Figure | Purpose | Candidate local file/status |
|---|---|---|
| Figure 1 | Data coverage and/or raw seasonal isotope ratio vs latitude | `fig1_data_coverage.png`, `fig2_harmonic_summary.png`, or combined remake |
| Figure 2 | Phasor correction schematic/decomposition | `fig8_phasor_decomposition.png` or simplified remake |
| Figure 3 | Corrected ratios/alpha by site and SH constraint | `fig9_corrected_ratio.png`, `fig10_alpha_constraint.png`, or remake |
| Figure 4 | Final KIE comparison with lab values | `fig10_alpha_constraint.png` or new compact summary |
| Optional Figure 5 or SI Figure | SH source-region robustness, emphasizing modest tropical influence and moving delayed-NH sensitivity into SI if visually dominant | Existing `fig24` or a new, less NH-centered remake |

Biomass-burning figures should be SI-only unless the manuscript story changes:

```text
fig15_bb_seasonality_by_band.png
fig16_bb_source_phasor_comparison.png
fig17_bb_correction_comparison.png
```

## Supplementary Information Design

The SI should be a structured technical archive rather than a collection of leftover material.

Recommended SI sections:

1. Site list, site metadata, data coverage, and clean-site criteria.
2. Harmonic fitting details and convention checks.
3. Wetland seasonality extraction from Li et al. 2026.
4. Wetland `deltaD-CH4` source-signature database.
5. Monte Carlo parameter table and uncertainty propagation.
6. Biomass burning correction and interpretation.
7. OSSE validation.
8. Block bootstrap and individual-year stability.
9. Harmonic model sensitivity.
10. Uncertainty attribution.
11. Southern Hemisphere wetland source-region sensitivity.
12. Expanded discussion notes from existing experiments, edited into readable subsections.

The SI should preserve the valuable experiment discussions, but rewritten in manuscript style and connected to primary citations. It should not read like a local lab notebook.

Useful existing discussion sources to mine for SI:

```text
experiments/KIE_sites/Discussion.md
experiments/KIE_sites/RESULT.md
experiments/KIE_sites/NH_high_latitude_discussion.md
experiments/KIE_sites/phase14_transport_lag_literature_review.md
experiments/KIE_sites/memory/*.md
```

## Known Uncertainties To Resolve Before Drafting Final Text

1. Exact citation and data-access statement for atmospheric `deltaD-CH4`.

2. Final citation for the `delta13C_wetland = -62 permil` prior.

3. Exact Phase 6 and Phase 14 numerical values from current JSON outputs.

4. Whether Li et al. 2026 should be cited as the ESSD article, the Zenodo dataset, or both.

5. Literature provenance for sink fractions and non-OH sink KIE values used in `ratio_to_alpha_13c`.

6. Whether the final manuscript should report both Saueressig and Cantrell values neutrally or use one as the primary lab reference and the other as a comparison.

7. Whether the current Phase 14 figure overemphasizes delayed Northern Hemisphere wetland influence. If so, make a new figure that foregrounds the stable Southern Hemisphere constraint and the modest tropical sensitivity, and move the delayed-NH envelope to SI.

## Drafting Rule

Do not write final-sounding prose around an unverified number. Mark the sentence as a temporary verification note in the working draft until the claim ledger has a source, then replace it with the verified value and citation/provenance entry.
