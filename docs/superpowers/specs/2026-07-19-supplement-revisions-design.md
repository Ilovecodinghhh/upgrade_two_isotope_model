# KIE Sites Supplement Revisions Design

## Goal

Revise the KIE_sites Supplementary Information so that Figures S8 and S14, the sink-parameter comparison, the Southern Hemisphere transport sensitivity, and the biomass-burning discussion are internally consistent, scientifically justified, and reproducible from the analysis outputs.

## Scope

The work modifies the Supplementary Information and the manuscript-figure generator. It does not change the primary atmospheric observations, the Phase 6 preferred KIE constraint, or the Phase 14 mass-conserving sensitivity calculations.

## Figure S8

Panel (b) currently compares absolute deterministic corrected ratios against a laboratory reference range. CGO and SPO retain the uniform `delta13C_wetland` value in the banded sensitivity, so their grey uniform markers and blue banded markers are identical; the blue markers are drawn later and hide the grey markers. The deterministic values also differ from earlier figures that display Monte Carlo medians, especially at SPO.

Panel (b) will instead display the deterministic sensitivity shift

```text
Delta R = R_banded - R_uniform.
```

The uniform baseline will be shown at zero and the banded result at its site-specific shift. Slight vertical marker offsets will keep coincident zero-shift markers visible for CGO and SPO. The bulk-sink laboratory band will be removed because an absolute KIE reference range is not meaningful on a delta axis. The caption will explicitly state that CGO and SPO have zero shift because the banded sensitivity retains the uniform source signature for the Southern Hemisphere/ambiguous category.

## S6 Sink Parameters

Section S6 will distinguish the component-specific parameter set used by the analysis from the compact harmonized set in Riddell-Young et al. (2025). Differing values will be tied to their original sources:

- chlorine sink magnitude: Hossaini et al. (2016);
- soil sink magnitude: Dutaur and Verchot (2007) and the global methane budget;
- `alphaD_Cl = 1.508`: room-temperature measurement of Saueressig et al. (1996);
- soil carbon and hydrogen KIEs: King et al. (1989) and Snover and Quay (2000);
- apparent stratospheric carbon and hydrogen fractionation: Rice et al. (2003) and McCarthy et al. (2003).

The text will not claim that the base values supersede Riddell-Young. It will state that the parameterizations make different representative choices from experimentally or observationally supported ranges, while their sink-weighted net KIEs remain close.

## S7 Southern Hemisphere Sensitivity

The section will explain that nominal deterministic CGO/SPO values use fixed observed and source phasors, convert each site-level corrected ratio to alpha, and summarize the two site values without uncertainty sampling. The preferred main-text result instead samples observational, wetland, and sink uncertainties; combines the corrected-ratio samples using inverse-variance site weights; applies the nonlinear bulk-sink-to-OH conversion; and reports the resulting median and 95% interval. These different estimators need not share the same central value.

Only mass-conserving source-region mixtures will remain in the Supplementary Information. All discussion and tables for the older additive transport stress test will be removed. The three plotted panels remain mass-conserving: the full grid and the two one-dimensional slices all reassign weight among local Southern Hemisphere, tropical, and delayed Northern Hemisphere wetland phasors while maintaining weights that sum to one.

The `2.8 month` Northern Hemisphere phase shift will be justified as an annual-harmonic response rather than a literal parcel travel time. For a first-order mixing response,

```text
H(omega) = 1 / (1 + i omega tau),
```

interhemispheric exchange or transit times of approximately `1.3-1.5 yr` imply annual-cycle attenuation near `0.10-0.12` and a phase lag near `2.8 months`. Geller et al. (1997), Patra et al. (2011), and Holzer and Waugh (2015) will provide the transport evidence. The text will state that the tested weights are effective phasor response weights, not transported mass fractions or station footprints.

## S8 Biomass Burning

The rationale for excluding biomass burning from the primary correction will not rely on the previous diagnostic uncertainty ranking. Instead, it will explain that fire emissions are episodic in space and time, have region-dependent fire seasons, and exhibit substantial interannual variability. A short 2005-2010 latitude-band climatology is therefore not a robust universal source phasor for CGO and SPO. Biomass-burning methane is also a relatively small part of the global methane source compared with microbial emissions. These points will be supported by van der Werf et al. (2017) and Saunois et al. (2020).

The existing biomass-burning figures will be presented as structural diagnostics. Event-resolved emissions, isotope signatures, transport, and station footprints will be identified as future work.

## Figure S14

The current Figure S14 will be removed from the Supplementary Information because it converts deterministic corrected ratios with a pure-OH formula, omits non-OH sink contributions and propagated uncertainty, and duplicates conclusions already shown more completely elsewhere. The underlying Phase 10 artifact will remain in the repository as an internal diagnostic.

The existing uncertainty-attribution Figure S15 will be renumbered to Figure S14 in the body, caption, and figure list.

## Verification

- Add a regression test for the Figure S8 shift data, including zero shifts for CGO and SPO.
- Regenerate the manuscript figures and visually inspect Figure S8.
- Verify that the Supplementary Information contains no additive second-diagnostic discussion and no Figure S15 reference.
- Verify that the new references appear in the bibliography with correct DOIs.
- Run the complete test suite.

