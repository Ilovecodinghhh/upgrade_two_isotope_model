# SCA Carbon Isotopes Phase 1-3 Method

## Scope

This document records the reproducible Phase 1-3 pipeline for the
`SCA_carbon_isotopes` experiment.

Phase 1 makes the project runnable from the repository workspace, replacing
absolute desktop paths with project-relative paths and writing standard output
tables under `results/`.

Phase 2 separates three definitions of seasonal cycle amplitude (SCA) so later
trend analyses can compare old and revised methods.

Phase 3 adds formal coverage screening for fixed analysis periods. It separates
site-year usability from period-level station eligibility, so later trend fits
can use only years and stations that satisfy predefined coverage rules.

## Input Data

Raw input files are NOAA/INSTAAR monthly `CH4C13 SIL` site files stored in:

```text
experiments/SCA_carbon_isotopes/SCA_noaa_ch4_isotopes/
```

Each file is read with comment lines skipped and the following data fields:

```text
site year month value
```

The pipeline writes a normalized monthly table:

```text
experiments/SCA_carbon_isotopes/results/site_monthly_clean.csv
```

## Coverage Rule

Each site-year is marked usable when it has:

```text
n_months >= 8
n_quarters >= 3
```

The coverage table is:

```text
experiments/SCA_carbon_isotopes/results/site_year_coverage.csv
```

Years failing this rule are retained in the output table for diagnostics, but
later trend calculations should filter to `usable == True`.

## Fixed-Period Eligibility

The Phase 3 station-level screening table is:

```text
experiments/SCA_carbon_isotopes/results/site_period_eligibility.csv
```

The default fixed periods are:

| Period | Years | Role |
|---|---:|---|
| `2002_2022` | 2002-2022 | Long-period analysis |
| `2016_2022` | 2016-2022 | Recent-period analysis |
| `2020_2022` | 2020-2022 | Diagnostic only |

For each site and period, the pipeline reports:

```text
expected_years
observed_years
usable_years
usable_fraction
eligible
diagnostic_only
```

The eligibility criterion is:

```text
usable_fraction = usable_years / expected_years
eligible = usable_fraction >= 0.70
```

The denominator is the number of calendar years expected in the period, not the
number of years observed in the local file. This prevents sparsely sampled
stations from appearing artificially complete. Periods shorter than 5 years are
marked `diagnostic_only == True`; the default `2020_2022` period should be used
for recent-pattern checks, not as the primary trend constraint.

The pipeline also writes period-specific input tables:

```text
experiments/SCA_carbon_isotopes/results/period_inputs_2002_2022.csv
experiments/SCA_carbon_isotopes/results/period_inputs_2016_2022.csv
experiments/SCA_carbon_isotopes/results/period_inputs_2020_2022.csv
```

These tables contain only usable site-years from stations that pass the
period-level eligibility screen.

## SCA Definitions

The Phase 2 output table is:

```text
experiments/SCA_carbon_isotopes/results/site_yearly_sca.csv
```

It reports three SCA columns:

| Column | Definition | Role |
|---|---|---|
| `sca_raw_range` | Annual `max(value) - min(value)` using raw monthly values | Reproduces the exploratory range-style metric |
| `sca_detrended_range` | Annual range after removing each site's linear trend | Reduces long-term trend leakage into annual amplitude |
| `sca_harmonic` | `2 * sqrt(B^2 + C^2)` from an annual sine/cosine fit to detrended monthly values within each site-year | Preferred candidate for later trend analysis |

The annual harmonic model is:

```text
y(month) = c0 + B sin(2 pi (month - 0.5) / 12)
                + C cos(2 pi (month - 0.5) / 12)
```

## Reproduction Command

Run from the repository root:

```bash
python experiments/SCA_carbon_isotopes/analysis/sca_pipeline.py
```

This writes:

```text
results/site_monthly_clean.csv
results/site_year_coverage.csv
results/site_yearly_sca.csv
results/site_yearly_sca_usable.csv
results/site_period_eligibility.csv
results/period_inputs_2002_2022.csv
results/period_inputs_2016_2022.csv
results/period_inputs_2020_2022.csv
```

## Current Phase 1-3 Output Summary

The current run reads 24 sites and 4607 monthly records spanning 1998-2022.
It produces 423 site-year rows, of which 381 pass the coverage rule.

The fixed-period screens currently retain:

| Period | Eligible sites | Period input rows |
|---|---:|---:|
| `2002_2022` | 12 / 24 | 248 |
| `2016_2022` | 20 / 24 | 136 |
| `2020_2022` | 19 / 24 | 57 |
