# SCA Carbon Isotopes Phase 4-8 Outputs

## Scope

This document records the reproducible Phase 4-8 additions to the
`SCA_carbon_isotopes` experiment. These phases start from the Phase 3 fixed
period input tables and generate trend estimates, station metadata, diagnostic
figures, a legacy-method comparison, and a standalone results draft.

## Inputs

Required Phase 3 inputs:

```text
experiments/SCA_carbon_isotopes/results/period_inputs_2002_2022.csv
experiments/SCA_carbon_isotopes/results/period_inputs_2016_2022.csv
experiments/SCA_carbon_isotopes/results/period_inputs_2020_2022.csv
experiments/SCA_carbon_isotopes/results/site_period_eligibility.csv
```

Station metadata are fetched from the NOAA GML Measurement Sites table:

```text
https://gml.noaa.gov/dv/site/
```

The legacy diagnostic comparison uses:

```text
experiments/SCA_carbon_isotopes/SCA_noaa_ch4_isotopes/sca_trend_summary.csv
```

## Trend Method

For each site, period, and SCA metric, the Phase 4-8 module fits:

```text
SCA(year) = intercept + slope * year
```

The output reports `slope_per_year`, `slope_per_decade`, t-based 95% confidence
intervals, p values, R2, and the number of usable years. The primary metric is
`sca_harmonic`; `sca_detrended_range` and `sca_raw_range` are retained as method
sensitivity checks.

The module also reports a robust Theil-Sen slope estimate for each site-period
metric. The robust columns are:

```text
sen_slope_per_year
sen_slope_per_decade
sen_intercept
sen_ci95_low
sen_ci95_high
sen_ci95_low_per_decade
sen_ci95_high_per_decade
```

The Theil-Sen confidence intervals are estimated by bootstrap resampling of the
site-year SCA pairs. These robust estimates are intended as sensitivity checks
against short-window leverage and non-Gaussian residual behavior, not as a
replacement for reporting the primary OLS diagnostic table.

## Outputs

Tables:

```text
experiments/SCA_carbon_isotopes/results/site_period_trends.csv
experiments/SCA_carbon_isotopes/results/site_metadata.csv
experiments/SCA_carbon_isotopes/results/site_period_trends_with_metadata.csv
experiments/SCA_carbon_isotopes/results/period_metric_summary.csv
experiments/SCA_carbon_isotopes/results/trend_robustness_summary.csv
experiments/SCA_carbon_isotopes/results/legacy_sca_trend_comparison.csv
```

Figures:

```text
experiments/SCA_carbon_isotopes/figures/phase6_sca_harmonic_trends_by_latitude.png
experiments/SCA_carbon_isotopes/figures/phase6_sca_harmonic_period_comparison.png
experiments/SCA_carbon_isotopes/figures/phase6_sca_method_sensitivity.png
experiments/SCA_carbon_isotopes/figures/phase6_period_coverage_heatmap.png
```

Standalone result draft:

```text
experiments/SCA_carbon_isotopes/drafts/sca_experiment_results_draft.md
```

## Reproduction Command

Run from the repository root:

```bash
python experiments/SCA_carbon_isotopes/analysis/sca_pipeline.py
python experiments/SCA_carbon_isotopes/analysis/sca_phase4_to8.py
```

The first command regenerates Phase 1-3 tables. The second command regenerates
Phase 4-8 trend tables, metadata, comparison tables, and figures.
