# KIE Sites Manuscript Citation Audit

This note distills the Markdown literature in `ImportantReferences` and
`Other References` into manuscript-ready citation decisions. It is an
internal audit document: manuscript text should cite the primary papers,
dataset DOIs, and official data products listed here, not the local Markdown
files.

## Scope

- Manuscript language: English.
- Citation rule: every observational or gridded data product must point to its
  original dataset, DOI, archive page, or official data citation.
- Scientific emphasis: present the strongest robust result first; keep the
  Northern Hemisphere wetland-to-Southern Hemisphere transport sensitivity as
  a structural sensitivity/SI item rather than the main claim.
- Notation: this document uses `delta13C`, `deltaD`, and `permil` for
  manuscript drafting consistency.

## Literature Read

The local literature was reviewed in the following thematic groups.

| Theme | Key local Markdown sources read | Main use in manuscript |
|---|---|---|
| Wetland emissions and seasonality | Li2026ESSD, Saunois2020ESSD, Saunois2025ESSD, Liu2025Nature, Lin2024 via Liu references | Primary wetland emission fields, tropical/high-latitude context, and budget-scale uncertainty |
| Wetland source signatures | Ganesan2018GRL, Douglas2021BG/EGU, Riddell-Young2025PNAS/SI, Sherwood2017ESSD | `delta13C` and `deltaD` source priors and uncertainty ranges |
| Classic isotope mechanisms | Whiticar1986GCA, Whiticar1999CG, Waldron1999GCA, Chanton2006JGR, Quay et al. references | Historical foundation for methanogenesis, oxidation, and environmental-water controls |
| Atmospheric isotope constraints | Rice2016PNAS, Basu2022ACP, Chandra2024CommEarth, Yu2026NatureCom, Warwick2016ACP | How atmospheric CH4, `delta13C`, and `deltaD` constrain sector attribution |
| OH, sink, and KIEs | Cantrell1990, Saueressig1996/2001, Snover & Quay2000, Montzka2011, Turner2017, Naus2019, Zhao2023, Nguyen2020 | Sink fractionation, OH uncertainty, and chemical feedback sensitivity |
| Biomass burning and fire isotopes | Worden2017NatComm, Nguyen2020GRL, GFED references, Umezawa2011/2012 | SI treatment of biomass burning corrections and CO-OH feedback |
| Transport and lag | Geller1997, Levin & Hesshaimer1996, Patra2011, Holzer & Waugh2015, Yang2019, Schuck2024, phase14 review | Justification for effective phasor attenuation and lag as sensitivity terms |

## Resolved Uncertainty Map

| Uncertainty item | Decision for the draft | Primary citations and data products | Confidence and remaining checks |
|---|---|---|---|
| Natural wetland emission dataset | Use Li et al. (2026) as the primary gridded natural vegetated wetland CH4 emission seasonality product. Describe it as local emissions, not source-receptor transport. | Li et al., ESSD, DOI `10.5194/essd-18-3507-2026`; Zenodo data DOI `10.5281/zenodo.18870108`. Use Saunois et al. (2020, 2025) for budget context: `10.5194/essd-12-1561-2020`, `10.5194/essd-17-1873-2025`. | High. The current extraction reads `ImportantReferences/Li2026ESSD_DS.nc`, variable `wetch4`, monthly 1 degree x 1 degree emissions in kg CH4 cell-1 month-1, and uses the 2005-2010 climatology. Cite both article and data DOI. |
| Uniform `delta13C_wetland = -62 permil` | Defensible only as a global wetland/microbial prior. Do not present it as site-specific. | Ganesan et al. (2018), first spatially resolved wetland map, DOI `10.1002/2018GL077536`; global flux/area-weighted wetland mean about `-62 permil`, high-latitude mean about `-67.8 permil`, tropical mean about `-56.7 permil`. Riddell-Young et al. (2025), PNAS, DOI `10.1073/pnas.2516543122`, gives global microbial `delta13C = -62.0 +/- 1.3 permil`. | High for a base-case global prior; medium for regional interpretation. The SI should show sensitivity to the Ganesan spatial gradient or latitudinal bands. |
| Wetland `deltaD-CH4` source signature | Use environmental-water/precipitation controls as the base physical rationale. Avoid saying `deltaD` uniquely diagnoses methanogenic pathway. | Douglas et al. (2021), Biogeosciences, DOI `10.5194/bg-18-3505-2021`; data DOI `10.6084/m9.figshare.13194833.v2`. Waldron et al. (1999), GCA; Chanton et al. (2006), JGR; Whiticar et al. (1986), GCA; Whiticar (1999), Chem. Geol. Riddell-Young et al. (2025) gives a modern wetland `deltaD` regression and global microbial `deltaD` estimate. | High for geographic dependence; medium for low-latitude and wetland-type detail. If we use a precipitation isotope field, cite the precipitation isotope product explicitly. |
| Precipitation/environmental water isotope field | Use the OIPC product/version explicitly, not the local CSV alone. In the current source-signature database this is mainly a predictor/cross-check; many recommended site wetland `deltaD` values are Douglas et al. zonal means. | OIPC v3.1 monthly precipitation `delta2H` values retrieved 2026-05-20 from the University of Utah OIPC page; Bowen and Revenaugh (2003), WRR, DOI `10.1029/2003WR002086`; Bowen et al. (2005), WRR; IAEA/WMO GNIP for measurement basis where used. | High for product identity. In methods, state coordinates/elevations came from station metadata and avoid citing the local CSV as a source. |
| Atmospheric CH4 mole fractions | Cite the NOAA/GML station flask/surface product actually used by the pipeline, not the global trend product. | Lan et al. (2025), NOAA GML Global Greenhouse Gas Reference Network CH4 surface flask data, version 2025-08-15, DOI `10.15138/VNCZ-M766`; calibration scale `CH4_X2004A`. | High. Current files are under `sitesdata/methane_ppb/noaa_gml_2025_event/`, downloaded 2026-05-17 from the NOAA GML surface flask archive. The global trend DOI `10.15138/P8XG-AA10` should only be cited if the trend product itself is used. |
| Atmospheric `delta13C-CH4` observations | Cite the NOAA/INSTAAR data product, not local processed files. | Michel et al. (2023), NOAA GML stable isotopic composition of atmospheric methane (`13C`), version 2023-09-21, DOI `10.15138/9p89-1x02`; calibration scale VPDB. Rice et al. (2016) is useful for historical archive-air context. | High. Current phase1 reads NOAA/INSTAAR event files under `sitesdata/isotope_d13C/noaa_instaar_2023_event/` and applies the file QC flag. No non-NOAA `delta13C` harmonization is used in the present phase1 extraction. |
| Atmospheric `deltaD-CH4` observations | Treat as a sensitive but now pinned citation item. The current pipeline uses Riddell-Young et al. raw station observations rather than the older INSTAAR-only product. | Riddell-Young et al. (2025), PNAS, DOI `10.1073/pnas.2516543122`; NOAA GML archive DOI `10.15138/setb-jy31`; Riddell-Young `dD_GlobMean` package copied into `sitesdata/isotope_dD/`; Dasgupta et al. (2025), AMT, DOI `10.5194/amt-18-6591-2025`, for interlaboratory harmonization. | High. DOI `10.15138/setb-jy31` resolves to the NOAA GML archive page titled "Atmospheric data analyzed in 'Microbial driver of 2006-2023 CH4 growth indicated by trends in atmospheric deltaD-CH4 and delta13C-CH4', PNAS". The code prefers the bare site-code/INSTAAR lab file if present and otherwise chooses the lab with most observations. |
| Global source-signature database | Use NOAA/GML/Sherwood database for source-signature provenance, not a derived local table. | Sherwood et al. (2017), ESSD article DOI `10.5194/essd-9-639-2017`, data DOI `10.15138/G3201T`; updated NOAA source-signature database DOI `10.15138/92GG-EY58`; Schwietzke et al. (2016), Nature, for historical isotope-database use. | High. Cite the updated database if drawing modern source signatures; cite Sherwood/Schwietzke for historical context. |
| OH `delta13C` KIE | Keep as a sensitivity, not a single settled constant. | Cantrell et al. (1990), JGR, DOI `10.1029/JD095iD13p22455`; Saueressig et al. (2001), JGR, DOI `10.1029/2000JD000120`; Basu et al. (2022), ACP, DOI `10.5194/acp-22-15351-2022`, explicitly discusses the Cantrell/Saueressig choice. | High. Report both values or a range; this is a defensible uncertainty item. |
| OH `deltaD` KIE and other sink KIEs | Use Saueressig et al. for OH/O(1D) and Cl kinetic isotope effects; use soil-uptake KIE literature for soil sink. State the values actually used in phase5/phase6. | Saueressig et al. (1996), GRL, `CH4 + Cl` D/H KIE; Saueressig et al. (2001), JGR; Snover and Quay (2000), GBC, DOI `10.1029/1999GB900089`; King et al. (1989) where soil `delta13C` KIE is used; Riddell-Young SI Table S3 for a compact modern parameter set. | High. The phase5/phase6 values are now aligned against Riddell-Young SI Table S3 in the table below; differences should be shown as SI sensitivity unless the analysis is rerun. |
| OH magnitude, trend, and sink attribution | Treat OH as a major uncertainty in attribution, but not as a free explanation for every pattern. | Montzka et al. (2011), Science; Turner et al. (2017), PNAS; Naus et al. (2019), ACP, DOI `10.5194/acp-19-407-2019`; Zhao et al. (2023), ACP, DOI `10.5194/acp-23-789-2023`; Saunois et al. (2020, 2025). | High for uncertainty framing. Avoid claiming a unique OH trend unless tied to a cited inversion or model. |
| Chemical feedback and perturbation lifetime | Use as a caution for decadal interpretations and biomass-burning/CO events. | Prather (1994, 1996); Prather et al. (2012); Nguyen et al. (2020), GRL, DOI `10.1029/2019GL085706`; Holmes et al. (2013), ACP. | High for SI discussion. This supports not over-interpreting lagged seasonal signals as direct emissions alone. |
| Transport attenuation and phase lag | Keep as an effective annual phasor sensitivity, not literal transported mass fraction or physical travel time. | Geller et al. (1997), interhemispheric exchange from SF6/CFCs; Levin and Hesshaimer (1996), SF6 and CFC constraints; Patra et al. (2011), TransCom CH4 interhemispheric exchange; Holzer and Waugh (2015), transit-time distributions; Yang et al. (2019) and Schuck et al. (2024) for modern transport constraints. | Medium. Good for SI sensitivity; do not make it a central result unless transport model evidence is added. |
| Biomass burning correction | Put detailed discussion in SI. Use Worden et al. to show biomass burning can influence isotope budget; use GFED as the primary fire emission product if used numerically. | Worden et al. (2017), Nature Communications, DOI `10.1038/s41467-017-02246-0`; van der Werf/Randerson GFED4.1s references; van Marle et al. (2017) BB4CMIP; Nguyen et al. (2020) for CO-OH feedback; Umezawa et al. (2011/2012) for atmospheric isotope observations and pyrogenic signature context. | Medium-high. Keep main text light unless biomass burning materially changes headline estimates. |
| Tropical wetland influence on Southern Hemisphere stations | Present as modest and physically plausible; do not overstate as direct transport evidence. | Liu et al. (2025), Nature, DOI `10.1038/s41586-025-08900-8`; Lin et al. (2024), Nature Communications, DOI `10.1038/s41467-024-55266-y`; Saunois et al. (2025); He et al. (2026), Science Advances, DOI `10.1126/sciadv.adz9007`. | Medium-high. Good main-text framing: tropical emissions and OH shape broad seasonal amplitude; station response remains a convolution with transport and sink seasonality. |
| Isotopologue inversions for tropical/subtropical attribution | Use as supporting context, not as our data source unless we reproduce their setup. | Yu et al. (2026), Nature Communications, DOI `10.1038/s41467-026-72668-2`; Chandra et al. (2024), Communications Earth & Environment, DOI `10.1038/s43247-024-01286-x`; Basu et al. (2022), ACP. | Medium. Useful to show current literature disagreement and the value of isotope constraints. |

## Current Pipeline Provenance Confirmed

These are the data products actually traced in the current analysis code and
local data headers. Manuscript citations should use these original products,
not local processed filenames.

| Pipeline use | Current product/version | How it enters the analysis | Manuscript citation action |
|---|---|---|---|
| Site CH4 mole fraction | NOAA GML Global Greenhouse Gas Reference Network surface flask CH4 event files; Lan et al. (2025), version 2025-08-15, DOI `10.15138/VNCZ-M766`; scale `CH4_X2004A`. | `phase2_harmonics.py` reads `sitesdata/methane_ppb/noaa_gml_2025_event/ch4_*_event.txt`, filters NOAA QC flag `.`, then forms monthly means over the isotope-overlap interval. | Cite Lan et al. (2025) station flask product. Do not substitute the global trend DOI unless that product is separately used for a figure or context statement. |
| Site `delta13C-CH4` | NOAA/INSTAAR surface flask isotope event files; Michel et al. (2023), version 2023-09-21, DOI `10.15138/9p89-1x02`; scale VPDB. | `phase1_data.py` reads `sitesdata/isotope_d13C/noaa_instaar_2023_event/ch4c13_*_event.txt`, filters event QC flag `.`, and computes paired monthly means only where `delta13C` and `deltaD` are both available. | Cite Michel et al. (2023) data product. Note that current phase1 does not mix in non-NOAA `delta13C` labs. |
| Site `deltaD-CH4` | Riddell-Young et al. (2025) `dD_GlobMean` raw observation files, harmonized across INSTAAR, MPI-BGC, IMAU, and TU/NIPR with Dasgupta et al. offsets. | `phase1_data.py` reads `sitesdata/isotope_dD/raw_observations/{site}_01D0_dat.txt` or equivalent lab file, prefers the bare site-code/INSTAAR record if available, otherwise selects the lab with most observations, then pairs monthly means with `delta13C`. | Cite Riddell-Young et al. (2025) PNAS, NOAA archive DOI `10.15138/setb-jy31`, and Dasgupta et al. (2025) AMT. Also note the NOAA page's data-provider contact/fair-use language if these data are central to the paper. |
| Wetland emission seasonality | Li et al. (2026) ESSD gridded natural vegetated wetland CH4 emissions; article DOI `10.5194/essd-18-3507-2026`; data DOI `10.5281/zenodo.18870108`. | `extract_wetland_seasonality.py` reads `ImportantReferences/Li2026ESSD_DS.nc`, variable `wetch4`, monthly 1 degree x 1 degree kg CH4 cell-1 month-1, and extracts 2005-2010 annual phasors for NH high, NH mid, tropics, SH extratropics, and global bands. | Cite both article and data DOI for quantitative wetland phasors. Describe bands as emission-seasonality predictors, not transport footprints. |
| Precipitation `delta2H` predictor | OIPC v3.1 monthly precipitation `delta2H` values from the University of Utah OIPC page, retrieved 2026-05-20; based on Bowen and Revenaugh (2003) and Bowen et al. (2005). | `build_dD_source_db.py` uses these station monthly values as an environmental-water predictor/cross-check. The recommended wetland `deltaD-CH4` values often come from Douglas et al. zonal wetland means rather than direct OIPC regression output. | Cite OIPC v3.1 and Bowen papers for the predictor, plus Douglas et al. (2021) for source-signature regressions/zonal means. |

## Current Pipeline Constants To Report

The manuscript methods table should match these values unless we deliberately
rerun the analysis. A short SI sensitivity table can compare them with
Riddell-Young SI Table S3 or alternative literature values.

| Quantity | Value used | Code location | Citation/status note |
|---|---:|---|---|
| Atmospheric `delta13C-CH4` reference | `-47.3 permil` | `phase6_phasor.py` | Report as analysis baseline; cite the observational data product and explain it as a representative atmospheric value. |
| Atmospheric `deltaD-CH4` reference | `-86.0 permil` | `phase6_phasor.py` | Report as analysis baseline; cite Riddell-Young/NOAA `deltaD` data source. |
| Wetland `delta13C-CH4` base case | `-62.0 +/- 5.0 permil` | `phase6_phasor.py` | Cite Ganesan et al. (2018) and Riddell-Young et al. (2025); keep regional sensitivity in SI. |
| Total CH4 source | `580 +/- 50 Tg yr-1` | `phase6_phasor.py` | Cite Saunois et al. (2025) global methane budget. |
| Wetland phasor uncertainty | `20%` on wetland B/C coefficients | `phase6_phasor.py` | Treat as analysis uncertainty, tied to Li et al. ensemble spread or documented as an assumed sensitivity. |
| OH sink fraction | `0.84 +/- 0.04` | `phase5_kie.py`, `phase6_phasor.py` | Cite global methane budget/sink literature; compare with Riddell-Young SI if using their compact parameter set. |
| Cl sink fraction | `0.035 +/- 0.01` | `phase5_kie.py`, `phase6_phasor.py` | Cite sink literature and SI parameter table. |
| Soil sink fraction | `0.06 +/- 0.02` | `phase5_kie.py`, `phase6_phasor.py` | Cite soil-uptake literature and global budget context. |
| Stratospheric sink fraction | `0.065` | `phase5_kie.py`, `phase6_phasor.py` | Current code uses closure with OH, Cl, and soil; state explicitly if used. |
| OH `alphaD` | `1.294 +/- 0.01` | `phase5_kie.py`, `phase6_phasor.py` | Cite Saueressig et al. (2001). Note that the earlier phase3 pure-OH diagnostic used `1.31`, but production phase5/6 uses `1.294`. |
| OH `alpha13C` reference values | `1.0039` and `1.0054` | `phase5_kie.py`, `phase6_phasor.py` | Cite Saueressig et al. (2001) and Cantrell et al. (1990); present as sensitivity. |
| Cl `alpha13C`, `alphaD` | `1.066 +/- 0.005`, `1.508 +/- 0.05` | `phase5_kie.py`, `phase6_phasor.py` | Cite Saueressig Cl KIE literature; note exact value differs slightly from some compact SI parameter sets. |
| Soil `alpha13C`, `alphaD` | `1.022`, `1.066` | `phase5_kie.py`, `phase6_phasor.py` | Cite King et al. (1989) or soil KIE literature for `alpha13C`, and Snover and Quay (2000) for soil D/H fractionation. |
| Stratospheric `alpha13C`, `alphaD` | `1.013`, `1.16` | `phase5_kie.py`, `phase6_phasor.py` | Cite stratospheric isotope observations/model literature; compare with Riddell-Young SI before final table submission. |

## Final KIE Table Alignment

Decision: the manuscript methods table should report the values actually used
in phase5/phase6, with a side-by-side SI comparison to Riddell-Young et al.
(2025) Table S3. Do not silently replace code values with Table S3 values
unless the analysis is rerun.

| Sink term | Phase5/phase6 value | Riddell-Young SI Table S3 | Final manuscript handling |
|---|---:|---:|---|
| OH fraction | `0.84 +/- 0.04` | `0.835` | Aligned within rounding. Cite Saunois/Montzka/global-budget context plus Riddell-Young Table S3. |
| Cl fraction | `0.035 +/- 0.01` | `0.035` | Aligned. Cite Hossaini et al. (2016) for global chlorine chemistry where using Riddell-Young's sink-fraction basis. |
| Stratospheric fraction | `0.065` | `0.06` | Close; code uses closure after OH, Cl, and soil. State as closure in Methods. |
| Soil fraction | `0.06 +/- 0.02` | `0.07` | Close but not identical. Keep code value in main analysis; note Riddell-Young value in SI sensitivity table. |
| OH `alpha13C` | `1.0039` and `1.0054` | `1.0054` | Keep both as sensitivity. Cite Saueressig et al. (2001), DOI `10.1029/2000JD000120`, and Cantrell et al. (1990), DOI `10.1029/JD095iD13p22455`. |
| OH `alphaD` | `1.294 +/- 0.01` | `1.294` | Aligned. Cite Saueressig et al. (2001). |
| Cl `alpha13C` | `1.066 +/- 0.005` | `1.066` | Aligned. Cite Saueressig et al. (1995), DOI `10.1029/95GL00881`. |
| Cl `alphaD` | `1.508 +/- 0.05` | `1.520` | Very close; keep code value and uncertainty. Cite Saueressig et al. (1996), DOI `10.1029/96GL03292`. |
| Soil `alpha13C` | `1.022` | `1.020` | Very close. Cite King et al. (1989), DOI `10.1029/JD094iD15p18273`, or Snover and Quay (2000) if using the paired soil KIE source. |
| Soil `alphaD` | `1.066` | `1.083` | Difference is material for the soil term but small in bulk because soil is a minor sink. Keep code value; include Riddell-Young value as a sensitivity. Cite Snover and Quay (2000), DOI `10.1029/1999GB900089`. |
| Stratospheric `alpha13C` | `1.013` | `1.003` | Difference is material for the stratospheric term but small in bulk. Keep code value unless rerun; cite stratospheric isotope literature and compare with Dyonisius et al. (2020), DOI `10.1126/science.aax0504`. |
| Stratospheric `alphaD` | `1.16` | `1.179` | Close. Keep code value unless rerun; compare with Riddell-Young Table S3. |
| Net sink `alpha13C` | `1.0078` using Saueressig OH; `1.0090` using Cantrell OH | `1.0082` | The code's range brackets Riddell-Young's net value. Report this as a useful consistency check. |
| Net sink `alphaD` | `1.2791` | `1.281` | Aligned to within about `0.002`; no concern for the current draft. |

Recommended wording:

"The base analysis uses the sink fractions and kinetic isotope effects listed
in Table Sx. These values reproduce the net deuterium sink fractionation in
Riddell-Young et al. (2025) to within about 0.002 in `alphaD`; for
`alpha13C`, the reported range reflects the laboratory OH KIE difference
between Saueressig et al. (2001) and Cantrell et al. (1990)."

## Data Availability Text

Use the following draft in the manuscript Data Availability section, adjusted
to match the final set of files actually used:

"NOAA GML surface-flask CH4 dry-air mole fractions are available from the NOAA
Global Monitoring Laboratory archive (Lan et al., 2025, DOI
`10.15138/VNCZ-M766`). NOAA/INSTAAR surface-flask `delta13C-CH4` measurements
are available from Michel et al. (2023, version 2023-09-21, DOI
`10.15138/9p89-1x02`). Atmospheric `deltaD-CH4` observations were obtained
from the NOAA GML archive record for Riddell-Young et al. (2025), DOI
`10.15138/setb-jy31`, and were harmonized following the laboratory-offset
framework described by Riddell-Young et al. and Dasgupta et al. (2025). The
NOAA archive requests appropriate citation of the data record and contact with
the relevant data providers when these isotope measurements are central to a
publication."

Important submission note: because the NOAA archive page includes
provider-contact and fair-use language for `deltaD-CH4`, `delta13C-CH4`, and
CH4 data, the final Data Availability section should not simply cite the DOI.
It should also acknowledge the NOAA archive record and the relevant data
providers.

## SI Decision On Ganesan Spatial Wetland Signature

Decision: include a compact Ganesan et al. (2018) spatial wetland
`delta13C-CH4` sensitivity in the SI, not in the main-text core result.

Rationale:

- The main analysis can keep the uniform `delta13C_wetland = -62 permil` base
  case because it is a defensible global wetland/microbial prior and matches
  the modern global microbial estimate in Riddell-Young et al. (2025).
- Ganesan et al. (2018) shows a meaningful regional gradient, with a more
  depleted high-latitude wetland mean around `-67.8 permil` and a less depleted
  tropical mean around `-56.7 permil`. This can matter most at NH sites where
  wetland source phasors are large.
- Keeping the sensitivity in SI protects the draft from overclaiming regional
  specificity while showing that the final KIE inference is not an artifact of
  using one global `delta13C` wetland source value.

Recommended SI implementation:

1. Base case: uniform `delta13C_wetland = -62 permil`.
2. Banded Ganesan case: assign high-latitude/boreal wetland influence near
   `-67.8 permil`, tropical wetland influence near `-56.7 permil`, and use the
   global `-62 permil` value where source-region attribution is ambiguous.
3. Report a small table showing the change in corrected seasonal ratio and
   inferred `alpha13C_OH` for each clean site.
4. Keep the text focused on robustness. Do not make the NH-to-SH transported
   wetland signal a central claim.

## Recommended Manuscript Wording

The following wording is drafted to be defensible and well cited.

1. Wetland emission seasonality:

   "We used monthly gridded natural vegetated wetland CH4 emissions from Li
   et al. (2026) to define emission-seasonality phasors. These fields are
   interpreted as local source seasonality rather than source-receptor
   footprints."

2. Wetland `delta13C` base case:

   "The base case used a uniform wetland `delta13C-CH4` source signature of
   `-62 permil`, equal to the global area/flux-weighted wetland mean in the
   spatially resolved wetland source-signature map of Ganesan et al. (2018).
   We treated the high-latitude to tropical gradient in that map as a structural
   uncertainty."

3. Wetland `deltaD` source model:

   "Wetland `deltaD-CH4` signatures were linked to precipitation or
   environmental-water `deltaD`, following the observed geographic dependence
   compiled by Douglas et al. (2021) and the earlier freshwater relationships
   of Waldron et al. (1999) and Chanton et al. (2006). We therefore interpret
   `deltaD-CH4` as reflecting hydrologic and pathway controls, not pathway alone."

4. Sink fractionation:

   "Because laboratory estimates of the OH `delta13C` kinetic isotope effect
   differ between Cantrell et al. (1990) and Saueressig et al. (2001), we report
   sensitivity to the OH KIE choice rather than treating one value as exact."

5. Transport lag:

   "Delayed remote-source terms are treated as effective annual-phasor
   sensitivities. Their amplitude and phase parameters should not be read as
   transported mass fractions or literal transit times."

6. Southern Hemisphere framing:

   "The weak response of Southern Hemisphere sites to imposed Northern
   Hemisphere wetland phasors is consistent with interhemispheric mixing and
   sink-seasonality attenuation. We therefore use these tests as robustness
   checks rather than as the central evidence for the manuscript."

## Reference Backbone For The Draft

### Main text anchors

- Global methane budget and uncertainty: Kirschke et al. (2013); Saunois et al.
  (2020, 2025); Dlugokencky et al. (2003, 2009, 2011); Nisbet et al. (2016,
  2019, 2020).
- Atmospheric observations: Lan et al. (2025), NOAA/GML station flask CH4
  product, DOI `10.15138/VNCZ-M766`; Michel et al. (2023) for
  `delta13C-CH4`, DOI `10.15138/9p89-1x02`; Riddell-Young et al. (2025) for
  `deltaD-CH4`; Dasgupta et al. (2025), AMT, for interlaboratory isotope
  harmonization.
- Wetland emissions: Li et al. (2026); Saunois et al. (2020, 2025); Bloom et al.
  (2012); Spahni et al. (2011); Koffi et al. (2020); Liu et al. (2025); Lin et
  al. (2024).
- Wetland source signatures: Ganesan et al. (2018); Douglas et al. (2021);
  Sherwood et al. (2017); Riddell-Young et al. (2025).
- Isotope mechanism and interpretation: Whiticar et al. (1986); Whiticar
  (1999); Waldron et al. (1999); Chanton et al. (2006); Quay et al. (1999);
  Bowen and Revenaugh (2003) and Bowen et al. (2005) for OIPC precipitation
  isotope fields.
- Modern isotope inversions and attribution: Rice et al. (2016); Schaefer et al.
  (2016); Worden et al. (2017); Basu et al. (2022); Chandra et al. (2024); Yu et
  al. (2026).
- Sink and KIE uncertainty: Cantrell et al. (1990); Saueressig et al. (1996,
  2001); Snover and Quay (2000); Montzka et al. (2011); Turner et al. (2017);
  Naus et al. (2019); Zhao et al. (2023); Nguyen et al. (2020).

### SI anchors

- Biomass burning: Worden et al. (2017); van der Werf/Randerson GFED4.1s;
  van Marle et al. (2017); Nguyen et al. (2020); Umezawa et al. (2011/2012).
- Transport lag and attenuation: Geller et al. (1997); Levin and Hesshaimer
  (1996); Patra et al. (2011); Holzer and Waugh (2015); Yang et al. (2019);
  Schuck et al. (2024).
- Paleo/historical context if needed: Wahlen et al. (1989); Ferretti et al.
  (2005); Mischler et al. (2009); Mitchell et al. (2013); Rubino et al. (2019);
  Hmiel et al. (2020).
- Satellite inversion context: Maasakkers et al. (2019), DOI
  `10.5194/acp-19-7859-2019`; He et al. (2026), DOI `10.1126/sciadv.adz9007`;
  Yu et al. (2026), DOI `10.1038/s41467-026-72668-2`.

## Items Now Mostly Resolved

1. The `-62 permil` `delta13C_wetland` value is usable, but only as a global
   wetland/microbial base-case prior. Ganesan et al. (2018) is the best primary
   wetland citation; Riddell-Young et al. (2025) is a strong modern global
   microbial consistency check.
2. Wetland `deltaD` should be tied to precipitation/environmental water and
   cited through Douglas et al. (2021), Waldron et al. (1999), and Chanton et
   al. (2006). This resolves the source-model rationale better than a fixed
   global `deltaD` assumption.
3. OH `delta13C` KIE remains a legitimate sensitivity because Cantrell et al.
   (1990) and Saueressig et al. (2001) differ and Basu et al. (2022) explicitly
   treats the choice as consequential.
4. The NOAA CH4 product is now pinned to the station flask/surface event files:
   Lan et al. (2025), version 2025-08-15, DOI `10.15138/VNCZ-M766`.
5. The OIPC precipitation isotope field is now pinned to OIPC v3.1, retrieved
   2026-05-20, with Bowen and Revenaugh (2003) and Bowen et al. (2005) as the
   underlying product citations.
6. The Riddell-Young atmospheric `deltaD-CH4` data DOI is now pinned to the
   NOAA GML archive DOI `10.15138/setb-jy31`; the archive page title and
   PNAS citation should be copied carefully into the final Data Availability
   or Methods text.
7. The phase5/phase6 KIE and sink constants are now extracted into a table
   above and aligned against Riddell-Young et al. (2025) SI Table S3. The code
   values should be reported as the analysis values; Riddell-Young differences
   should be shown as SI sensitivity, not silently substituted.
8. The NOAA archive data-access/fair-use requirements have been converted into
   a Data Availability draft that cites Lan et al. (2025), Michel et al. (2023),
   Riddell-Young et al. (2025), the NOAA archive DOI `10.15138/setb-jy31`, and
   Dasgupta et al. (2025).
9. The Ganesan spatial wetland `delta13C` sensitivity should be included in the
   SI as a compact robustness test, while the main text keeps the uniform
   `-62 permil` wetland base case.
10. The NH-to-SH transport lag/attenuation experiment is literature-supported as
   an effective phasor sensitivity, but not as a direct physical travel-time
   estimate.
11. Biomass burning is worth a detailed SI section, but it does not need to be a
   main-text emphasis unless its correction changes our headline result.

## Remaining Checks Before Submission

1. If the manuscript changes any sink/KIE value to match Riddell-Young SI Table
   S3 exactly, rerun the relevant phase5/phase6 analysis and update the results
   before submission.
2. If the OIPC precipitation values are used quantitatively beyond a
   predictor/cross-check, record the exact query method and station metadata
   source in Methods or SI for reproducibility.
3. If citing papers stored locally as "Article in Press" (for example Yu et al.
   2026), use the final official publication details and DOI, not the local
   folder title.
