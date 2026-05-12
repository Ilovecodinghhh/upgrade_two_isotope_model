# KIE Immunity — Plan & Prompts

**Last updated:** 2026-05-12  
**Depends on:** Phases 1–4 + v2 (real hemispheric δD) complete.

---

## Priority 1: FF Trend Robustness

The post-2007 FF trend reversal (+13 → −7 Tg/yr) is the strongest finding
but also the most provocative. It must survive every reasonable sensitivity
test before we can publish it.

### Phase 5 — Lifetime sensitivity

> In `experiments/KIE_immunity/analysis/`, write `phase5_tau_sensitivity.py`.
> Run the dual real-hemi-dD model at fixed τ = 8.0, 8.5, 9.0, 9.5, 10.0 yr,
> and also with He 2026's time-varying τ (declining from ~9.8 to ~9.0 over
> 2000–2020). For each, report: (a) post-2007 ΔFF trend (median ± 90% CI),
> (b) σ(FF), (c) KIE contribution %. The key question: does the negative FF
> trend survive across all reasonable τ values, or does it flip sign at some
> τ? If it flips, report the critical τ value. Save results to
> `results/phase5_tau_sensitivity.json`. Update RESULTS.md with a new
> "Phase 5" section.

### Phase 6 — OH-D KIE sensitivity

> Write `phase6_OHD_sensitivity.py`. The OH-D KIE (1.294–1.327) is the
> second-largest δD uncertainty after source signatures. Run dual real-hemi-dD
> with OH_D fixed at: 1.294 (Saueressig), 1.310 (midpoint), 1.327 (Cantrell),
> and 1.35 (He 2026 upper bound). Report post-2007 ΔFF and σ(FF) for each.
> If the FF trend reversal is robust to OH_D, that's a strong result. If it
> flips, report the OH_D threshold — analogous to the σ(Mic δD) = 25‰
> threshold from the dD_threshold experiment. Save to
> `results/phase6_OHD_sensitivity.json`.

### Phase 7 — Cl fraction sensitivity

> Write `phase7_Cl_sensitivity.py`. The Cl sink has the largest δD KIE
> (α = 1.52). Run dual real-hemi-dD with Cl fraction = 0.6% (Thanwerdas
> 2024), 2.0%, 3.5% (default), 5.0%, 6.5%, 10% (Allan 2007 upper). Report
> ΔFF trend, σ(FF), and variance decomposition for each. The Cl fraction
> amplifies the δD constraint weight — higher Cl means δD matters more. Save
> to `results/phase7_Cl_sensitivity.json`.

### Phase 8 — Combined robustness matrix

> Write `phase8_robustness_matrix.py`. Run a grid: τ ∈ {8.5, 9.0, 9.5} ×
> OH_D ∈ {1.294, 1.327} × Cl ∈ {0.6%, 3.5%, 6.5%} = 18 combinations.
> For each, report ΔFF trend. Present as a 3D table. The question: in how
> many of the 18 cells is ΔFF negative? If all 18 → the trend reversal is
> iron-clad. If only some → report the boundary conditions. This is the
> "how confident should we be?" table for the paper. Save to
> `results/phase8_robustness_matrix.json`.

---

## Priority 2: Bootstrap Confidence Intervals

The current variance decomposition has no uncertainty on the percentages.
"42.6% source signatures" could be 42.6 ± 2% or 42.6 ± 20%.

### Phase 9 — Bootstrap variance decomposition

> Write `phase9_bootstrap_variance.py`. Take the 400 MC iterations from
> the v2 variance decomposition. Bootstrap-resample 1000 times (sample 400
> with replacement). For each bootstrap sample, compute σ(FF) and the
> variance component percentages. Report 95% CI on each component. The key
> test: is the real-hemi "43% source sigs" significantly different from the
> offset "0% source sigs"? Save to `results/phase9_bootstrap.json`.

---

## Priority 3: Publication Figure

### Phase 10 — 3-panel variance decomposition figure

> Write `figures/fig_variance_v2.py`. Three-panel figure:
>
> **Panel A:** Stacked bar chart — σ²(FF) decomposed into KIE / Sig / τ /
> Residual for all 3 configs. Use distinct colors, add σ(FF) value labels
> on top of each bar.
>
> **Panel B:** FF time series (1999–2021, 5-yr smoothed) — median + 2σ
> shading for all 3 configs. Vertical dashed line at 2007. Annotate ΔFF
> trend for each.
>
> **Panel C:** Basu 2022 comparison — horizontal bar chart. Basu (13.0),
> our δ¹³C-only (0.7), dual offset (0.9), dual real hemi (0.8). Label
> "KIE spread (Tg/yr)".
>
> Style: Nature-family format, 180mm wide (2-column), 8pt font, colorblind-
> safe palette. Export PNG (300 dpi) + PDF. Save to
> `figures/fig_variance_v2.{png,pdf}`.

---

## Priority 4: Interhemispheric Exchange Sensitivity

### Phase 11 — τ_ex sensitivity

> Write `phase11_tau_ex.py`. The interhemispheric exchange time τ_ex is
> currently N(1.0, 0.1) yr. Test extreme values: τ_ex = 0.5 yr (fast
> mixing, hemispheres nearly coupled) and τ_ex = 2.0 yr (slow mixing,
> hemispheres nearly independent). If σ(FF) changes a lot → the hemispheric
> δD constraint genuinely exploits transport separation. If insensitive →
> the 2-box advantage comes purely from source-sig separation. Save to
> `results/phase11_tau_ex.json`.

---

## Priority 5: EDGAR Validation

### Phase 12 — Compare with bottom-up inventories

> Write `phase12_edgar_validation.py`. Plot our FF time series (median + CI
> from dual real-hemi-dD) against EDGAR 8.0 bottom-up. Load EDGAR data from
> `rel/data/CarbonTracker_CH4.xlsx` (column 9 = total, or find FF-specific
> columns). If our negative trend is consistent with EDGAR's leveling-off
> post-2010, that's independent validation. Also compare with IEA/GAINS if
> available. Save figure to `figures/fig_edgar_validation.{png,pdf}` and
> comparison stats to `results/phase12_edgar.json`.

---

## Priority 6: Paper-Ready Analysis

### Phase 13 — Table of all headline numbers with CIs

> Write `phase13_summary_table.py`. Compile every headline number from this
> experiment into one table with bootstrap 95% CIs:
> - σ(FF) for each config ± CI
> - KIE% ± CI
> - KIE spread ± CI
> - ΔFF trend ± CI
> - Agreement rate ± CI (if applicable)
>
> Export as LaTeX table (`results/table1.tex`) and CSV
> (`results/table1.csv`). This becomes Table 1 in the paper.

---

## Completed ✅

- [x] Phase 1–2: Variance decomposition (δ¹³C-only + dual offset) — v1
- [x] Phase 3–4: Basu 2022 comparison (δ¹³C-only + dual offset) — v1
- [x] v2 upgrade: real hemispheric δD (atm MC + source sig MC)
- [x] v2 variance decomposition (3-config comparison)
- [x] v2 Basu comparison (+ residual analysis)
- [x] common.py: `sample_atm_dD_hemi()`, `sample_source_sigs_hemi()`
