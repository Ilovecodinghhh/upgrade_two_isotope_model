# Draft Results: Seasonal Cycle Amplitude of NOAA/INSTAAR CH4 Carbon Isotopes

This standalone draft summarizes the current results of the
`SCA_carbon_isotopes` experiment after the Phase 1-8 update. It is intended as
source text for the manuscript Results section and Supplementary Information,
not as a final polished manuscript section.

## Data Provenance and Analysis Set

The analysis uses the NOAA/INSTAAR monthly `CH4C13 SIL` individual-flask site
files archived with NOAA GML trace-gas data documentation and access guidance
at the NOAA GML trace-gas archive.[^noaa-trace-gases] Station coordinates and
site descriptions were taken from the NOAA GML Measurement Sites table.[^noaa-sites]
The local Phase 1-3 pipeline reads 24 site files and 4607 monthly records
spanning 1998-2022. After applying the site-year coverage rule
(`n_months >= 8` and `n_quarters >= 3`), 381 of 423 site-years are retained as
usable.

Fixed-period station eligibility was evaluated using the fraction of expected
calendar years that pass the site-year coverage rule. The 2002-2022 period
retains 12 of 24 sites and 248 usable site-year rows. The 2016-2022 period
retains 20 of 24 sites and 136 usable site-year rows. The 2020-2022 period
retains 19 of 24 sites and 57 usable site-year rows, but is treated as
diagnostic only because the three-year window is too short to provide a primary
trend constraint.

Primary result tables:

- `../results/site_period_trends.csv`
- `../results/site_metadata.csv`
- `../results/site_period_trends_with_metadata.csv`
- `../results/period_metric_summary.csv`
- `../results/trend_robustness_summary.csv`
- `../results/legacy_sca_trend_comparison.csv`

Primary figures:

- `../figures/phase6_sca_harmonic_trends_by_latitude.png`
- `../figures/phase6_sca_harmonic_period_comparison.png`
- `../figures/phase6_sca_method_sensitivity.png`
- `../figures/phase6_period_coverage_heatmap.png`

## Long-Period Trends, 2002-2022

Using the harmonic definition of seasonal cycle amplitude (SCA), the
2002-2022 fixed-period analysis shows a weak positive tendency in annual
delta13C-CH4 SCA. Across the 12 eligible sites, the median harmonic SCA trend is
+0.034 per mil per decade, with an interquartile range from +0.010 to +0.042
per mil per decade. Eleven of the twelve site-level trends are positive, but
only one site-level trend is significant at p < 0.05 under the current
ordinary-least-squares diagnostic. The only negative harmonic trend is at SMO
(-0.005 per mil per decade; p = 0.869), while the largest positive long-period
trend is at MHD (+0.067 per mil per decade; p = 0.042).

The robust Theil-Sen diagnostic supports the same qualitative long-period
interpretation while giving a slightly smaller median slope. For 2002-2022
harmonic SCA, the median Theil-Sen trend is +0.028 per mil per decade, compared
with the OLS median of +0.034 per mil per decade. The OLS and Theil-Sen trends
have the same sign at all 12 eligible sites, and 11 of 12 Theil-Sen site-level
slopes remain positive. However, no 2002-2022 site has a Theil-Sen bootstrap
confidence interval that excludes zero, so the robust analysis favors a weak
network-level tendency rather than strong single-station claims.

The latitudinal pattern is modest. The median 2002-2022 harmonic trend is
+0.038 per mil per decade for northern mid/high-latitude sites, +0.033 per mil
per decade for southern mid/high-latitude sites, and +0.022 per mil per decade
for tropical sites. This pattern is consistent with a broadly weak increase in
SCA over the long period, rather than a sharply localized or hemisphere-limited
signal.

## Recent-Period Trends, 2016-2022

The recent 2016-2022 period shows a stronger positive tendency. Under the OLS
diagnostic, all 20 eligible sites have positive harmonic SCA trends. The median
OLS trend is +0.126 per mil per decade, with an interquartile range from +0.062
to +0.200 per mil per decade. The largest positive recent OLS trends occur at
BHD (+0.278), ALT (+0.249), ZEP (+0.219), WLG (+0.204), TAP (+0.200), and MHD
(+0.200 per mil per decade).

The Theil-Sen median for 2016-2022 is very similar (+0.128 per mil per decade),
but the site-level sign coherence is weaker than in the OLS diagnostic: 15 of
20 Theil-Sen slopes are positive. No 2016-2022 Theil-Sen bootstrap interval
excludes zero at individual sites. Thus, the robust estimator preserves the
positive network median while tempering the claim that every station shows a
robust positive recent trend.

The recent-period increase is visible across latitude bands. The median
harmonic trend is +0.171 per mil per decade for northern mid/high-latitude
sites, +0.147 per mil per decade for southern mid/high-latitude sites, and
+0.106 per mil per decade for tropical sites. Although no individual site-level
2016-2022 harmonic trend reaches p < 0.05, the OLS sign coherence and the
positive Theil-Sen network median together suggest that this interval contains
a stronger network-wide increase in SCA than the longer 2002-2022 period.

The contrast between 2002-2022 and 2016-2022 should be described cautiously.
The recent window has fewer years and therefore larger trend uncertainty at
individual stations. Nevertheless, the period-comparison figure indicates that
many sites have steeper recent-period trends than their long-period trends.
This pattern suggests either a recent acceleration in the processes controlling
the isotopic seasonal cycle or a sensitivity of short-window trend estimates to
interannual variability.

## Diagnostic 2020-2022 Window

The 2020-2022 diagnostic window behaves differently from the two primary
periods. The median harmonic trend is -0.164 per mil per decade; 14 of 19 sites
are negative, and the range is large (-1.441 to +1.129 per mil per decade).
This behavior is not interpreted as a robust trend because each station has at
most three annual SCA values. Instead, the 2020-2022 result is useful as a
stress test: it demonstrates that very short windows can produce large apparent
slopes and even reverse the sign of the network median. For the manuscript,
2020-2022 should remain in the Supplementary Information or be clearly labeled
as diagnostic-only.

## Sensitivity to SCA Definition

The three SCA definitions give similar OLS median trends for the two primary
periods. For 2002-2022, the median trends are +0.0349 per mil per decade for
the raw annual range, +0.0302 for the detrended annual range, and +0.0341 for
the harmonic SCA. For 2016-2022, the corresponding medians are +0.1286,
+0.1258, and +0.1259 per mil per decade. The Theil-Sen medians are slightly
more conservative for most metric-period combinations but preserve the same
main contrast between the weak long-period tendency and the stronger
recent-period tendency. This agreement supports using the harmonic SCA as the
primary metric while presenting the range-based metrics and robust trend
estimates as sensitivity tests.

The diagnostic 2020-2022 period is also internally consistent across SCA
definitions, but in the opposite direction: all three definitions have negative
median trends. This reinforces the conclusion that the short-window result is a
window-length effect or a short-term feature rather than a stable long-term
constraint.

## Comparison With the Legacy Exploratory Workflow

The updated fixed-period harmonic trends were compared with the legacy
exploratory trend table (`sca_trend_summary.csv`). This is a diagnostic
comparison rather than a strict like-for-like method swap: the legacy table uses
each site's available record length, whereas the updated values below use the
fixed 2002-2022 period and the Phase 3 eligibility screen. The old workflow also
used an SSA reconstruction step that effectively retained the original time
series before computing annual ranges, so it should be treated as an exploratory
range-based diagnostic rather than as the final seasonal component extraction.

For the 12 common 2002-2022 sites, the updated harmonic trends are usually
smaller than the legacy slopes. Nine of twelve sites have lower updated slopes,
with a median new-minus-legacy difference of about -0.017 per mil per decade.
The sign remains positive for most common sites, but SMO changes from a small
positive legacy slope to a slightly negative updated harmonic slope. This
comparison suggests that the updated workflow preserves the broad impression of
weak-to-positive long-period SCA changes while reducing sensitivity to the
legacy reconstruction/range method.

## Working Interpretation

The current results support a conservative manuscript narrative. The long
2002-2022 record indicates a weak, mostly positive change in delta13C-CH4 SCA
across the eligible network, with limited station-level significance. The
2016-2022 period shows a stronger and sign-coherent positive tendency, but it
should be framed as a recent-window feature with larger uncertainty rather than
as a standalone long-term trend estimate. The 2020-2022 diagnostic window
demonstrates why short periods should not be used as primary constraints.

In manuscript language, the strongest defensible statement is that the revised
coverage-controlled workflow identifies a weak long-period increase in isotopic
SCA and a stronger recent-period positive tendency. This result is robust at
the network-median level to the choice of SCA definition and to replacing OLS
with a Theil-Sen slope estimator, but individual site-level intervals are wide.
The result should not yet be framed as evidence for a single dominant
source-region mechanism. Instead, it provides an observational constraint on
changes in the seasonal isotopic cycle that can be compared against
source-region sensitivity experiments in the main analysis and Supplementary
Information.

## Items to Confirm Before Manuscript Integration

- Decide whether OLS remains the main reported estimator with Theil-Sen in the
  Supplementary Information, or whether both estimators should be reported side
  by side in a main-text table.
- Decide whether the 2016-2022 sign coherence should be emphasized in the main
  text or shown primarily as a recent-period diagnostic.
- Add exact figure numbering after the manuscript figure set is finalized.
- Add formal citations for the NOAA/INSTAAR isotope data product according to
  the README and citation text associated with the trace-gas archive.

[^noaa-trace-gases]: NOAA Global Monitoring Laboratory, trace-gas data archive:
    https://gml.noaa.gov/aftp/data/trace_gases/

[^noaa-sites]: NOAA Global Monitoring Laboratory, Measurement Sites:
    https://gml.noaa.gov/dv/site/
