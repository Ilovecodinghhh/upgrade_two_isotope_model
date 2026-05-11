# Title 2: Implementation Plan & Coding Agent Prompts

## Phase 1: Structured KIE Sensitivity — δ¹³C-only vs. Dual-Isotope

### Task 1.1: Fixed-KIE endpoint runs (4 models × 2 KIE values × 2 isotope modes)

**Prompt for coding agent:**
```
In the repository at https://github.com/Ilovecodinghhh/upgrade_two_isotope_model (branch: master),
create `analysis/kie_sensitivity_matrix.py` that:

1. Runs all 4 model variants (2x2_one, 2x2_two, 3x3_one, 3x3_two) under 6 configurations:

   Config A: OH_13C = 1.0039 (Saueressig), dual isotopes
   Config B: OH_13C = 1.0054 (Cantrell), dual isotopes
   Config C: OH_13C = 1.0039, δ¹³C-only (disable δD equation)
   Config D: OH_13C = 1.0054, δ¹³C-only
   Config E: OH_13C sampled U(1.0039, 1.0054), dual isotopes (= default)
   Config F: OH_13C sampled U(1.0039, 1.0054), δ¹³C-only

   For δ¹³C-only mode in 3x3 models: drop the δD row, reduce to 2×2 system.
   For δ¹³C-only mode in 2x2 models: use only the δ¹³C-derived FF/Mic solution.

2. N = 1000 MC iterations per configuration. All other parameters at `default` preset.
   OH_D KIE remains sampled normally (it only matters when δD equation is active).

3. Output per configuration:
   - CSV: year, FF_median, FF_p5, FF_p95, Mic_median, Mic_p5, Mic_p95, BB_median (if solved)
   - Store in `Output_KIE_sensitivity/<model>/<config>/results.csv`

4. Summary statistics:
   - "KIE spread" = |FF_median(Config A) - FF_median(Config B)| for dual-isotope
   - "KIE spread" = |FF_median(Config C) - FF_median(Config D)| for δ¹³C-only
   - "Immunity ratio" = KIE_spread(δ¹³C-only) / KIE_spread(dual) per year
   - Save to `Output_KIE_sensitivity/immunity_summary.csv`

Import model functions from `models/core.py`. Use presets from `models/inputs.py`.
```

### Task 1.2: KIE sensitivity for the 2×2 cross-validation diagnostic

**Prompt for coding agent:**
```
Create `analysis/kie_cross_validation.py`:

The 2×2 models solve δ¹³C and δD INDEPENDENTLY for FF and Mic.
This means you get TWO estimates of FF: one from δ¹³C, one from δD.

1. For 2x2_one.py and 2x2_two.py, extract both:
   - FF_from_d13C (per MC iteration, per year)
   - FF_from_dD (per MC iteration, per year)

2. Compute for each KIE endpoint (Saueressig vs. Cantrell):
   - FF_from_d13C shifts substantially (because KIE directly affects δ¹³C equation)
   - FF_from_dD shifts less (δD KIE is independently sampled, OH_D range is smaller)

3. Quantify:
   - δ¹³C-derived FF sensitivity to OH_13C KIE (Tg/yr per 0.001 KIE unit)
   - δD-derived FF sensitivity to OH_13C KIE (should be near zero — δD doesn't use OH_13C)
   - δD-derived FF sensitivity to OH_D KIE

4. Key insight: δD provides an INDEPENDENT CHECK that doesn't share the OH_13C KIE vulnerability.
   When the two estimates agree, the answer is robust regardless of KIE choice.

5. Output:
   - Figure: "δ¹³C-derived FF vs. δD-derived FF" scatter plot, colored by OH_13C KIE value
   - Figure: time series showing δ¹³C-derived and δD-derived FF bands at each KIE endpoint
   - CSV: per-year correlation and bias between the two estimates
```

---

## Phase 2: Quantifying the "Immunity Factor"

### Task 2.1: Bootstrap confidence interval comparison

**Prompt for coding agent:**
```
Create `analysis/immunity_factor.py`:

1. Define "KIE-driven ambiguity" for a given model and isotope configuration:
   - Run 1000 MC with OH_13C fixed at 1.0039 → get FF distribution D_low
   - Run 1000 MC with OH_13C fixed at 1.0054 → get FF distribution D_high
   - Ambiguity = |median(D_high) - median(D_low)| (in Tg/yr)
   - Also compute: overlap coefficient between D_low and D_high (Bhattacharyya)

2. Compute ambiguity for:
   - Each model variant (2x2_one, 2x2_two, 3x3_one, 3x3_two)
   - Each isotope mode (δ¹³C-only, dual)
   - Each year (2000–2022)

3. "Immunity factor" = Ambiguity(δ¹³C-only) / Ambiguity(dual)
   - If immunity_factor > 2.0: "dual isotopes halve the ambiguity" (title claim)
   - Report year-by-year and mean over 2007–2022

4. Bootstrap: resample MC iterations (with replacement) 500 times to get 
   confidence interval on the immunity factor itself.

5. Output:
   - `results/immunity_factors.csv`: year, model, ambiguity_d13C_only, ambiguity_dual, 
     immunity_factor, immunity_CI_low, immunity_CI_high
   - Figure: time series of immunity factor (with CI bands) for each model
   - Summary: "The dual-isotope approach reduces KIE-driven FF ambiguity by X ± Y%"
```

### Task 2.2: FF signature inventory interaction

**Prompt for coding agent:**
```
Create `analysis/inventory_kie_interaction.py`:

Test whether the KIE immunity depends on which FF signature dataset is used:

1. Run the immunity analysis (Task 2.1) three times:
   - FF signatures from EDGAR (`FF_d13C_GlobMC_EDGAR.csv` + `FF_dD_GlobMC_EDGAR.csv`)
   - FF signatures from CarbonTracker-CH₄ (`FF_d13C_GlobMC_CTCH4.csv` + `FF_dD_GlobMC_CTCH4.csv`)
   - FF signatures from GlobUnc (`FF_d13C_GlobUnc.csv` + `FF_dD_GlobUnc.csv`)

2. Compare immunity factors across inventories.
   - If similar: result is robust (stronger paper)
   - If different: report which inventory is more/less sensitive to KIE

3. Also test interaction with lifetime:
   - fixed_lifetime (τ=9.0) vs. default (time-varying)
   - Does time-varying lifetime change the immunity factor?

4. Output: 2D heatmap figure (inventory × lifetime → immunity factor)
```

---

## Phase 3: Policy-Relevant Attribution Robustness

### Task 3.1: "Robust" post-2007 source attribution

**Prompt for coding agent:**
```
Create `analysis/robust_attribution.py`:

The paper's key deliverable: a source attribution that is ROBUST to KIE choice.

1. For the 3x3_two model (most complete: dual-isotope, hemispheric):
   - Run with full KIE sampling (default preset)
   - Report FF, Mic, BB trends (2000–2022) with:
     a. Total uncertainty (from all MC sources combined)
     b. Breakdown: KIE contribution vs. source-signature contribution vs. lifetime contribution
   
2. Variance decomposition:
   - Run with KIE fixed (at midpoint) but all else sampled → variance_minus_KIE
   - Run with all sampled (default) → variance_total
   - KIE contribution to variance = variance_total - variance_minus_KIE
   - Report: "KIE accounts for X% of total FF uncertainty in δ¹³C-only, but only Y% in dual-isotope"

3. Trend significance test:
   - For FF emissions 2007–2022: fit linear trend, report slope ± SE
   - Is the trend significantly different from zero?
   - Does the answer change between δ¹³C-only and dual-isotope?
   - Key claim: "With dual isotopes, FF trend is [stable/declining] regardless of KIE"

4. Output:
   - Figure: FF and Mic time series with "envelope of KIE uncertainty" shown explicitly
   - Figure: Variance decomposition pie charts (δ¹³C-only vs. dual)
   - Table: Trend slopes and p-values for each source × isotope mode × KIE setting
```

### Task 3.2: Comparison with Basu 2022 posteriors

**Prompt for coding agent:**
```
Create `analysis/compare_basu2022.py`:

Basu 2022 (ACP) demonstrated large KIE sensitivity in their 3D inversion.
They switched between Saueressig and Cantrell and found major shifts.

1. Extract Basu 2022's reported FF emissions for their two KIE experiments
   (from their paper/SI — digitize Figure 5 or use reported values).
   
   Basu 2022 key numbers (approximate from their paper):
   - With Saueressig KIE: FF post-2007 increase ≈ +12 Tg/yr
   - With Cantrell KIE: FF post-2007 increase ≈ +25 Tg/yr
   - Spread ≈ 13 Tg/yr

2. Compare to your model's KIE spread:
   - Your δ¹³C-only spread (should be similar magnitude to Basu — ~10-15 Tg/yr)
   - Your dual-isotope spread (should be much smaller — this is the paper's point)

3. Generate comparison figure:
   - Basu 2022's KIE spread (gray bar, from literature)
   - Your δ¹³C-only KIE spread (blue bar)
   - Your dual-isotope KIE spread (red bar, much narrower)
   - Annotate with "X% reduction"

4. Implication: "If Basu et al. (2022) had incorporated δD constraints, their FF 
   attribution would have been Y Tg/yr more certain, potentially resolving the 
   disagreement with Riddell-Young et al. (2025)."
```

---

## Phase 4: Main Figures

### Task 4.1: Figure 1 — "KIE immunity demonstration"

**Prompt for coding agent:**
```
Create `figures/fig_kie_immunity.py`:

Layout: 3-panel figure (GBC style, single column 85mm or double column 170mm)

Panel A: "FF emissions under Saueressig vs. Cantrell KIE — δ¹³C-only"
- Time series 2000–2022
- Two shaded bands (blue = Saueressig, orange = Cantrell) with medians
- Large gap between bands = high sensitivity
- Y-axis: Fossil fuel emissions (Tg CH₄ yr⁻¹)

Panel B: "FF emissions under Saueressig vs. Cantrell KIE — dual-isotope (δ¹³C + δD)"
- Same format as Panel A
- Bands overlap substantially = low sensitivity
- Same y-axis scale as Panel A for visual comparison

Panel C: "Immunity factor over time"
- X-axis: Year (2000–2022)
- Y-axis: Immunity factor (ratio of KIE spreads)
- Horizontal dashed line at 2.0 ("halved")
- Lines for each model variant (2x2_one, 3x3_two)
- Shaded CI band from bootstrap

Style: matplotlib + seaborn. Nature-style clean aesthetics.
Export PDF + PNG (300 dpi).
```

### Task 4.2: Figure 2 — "Variance decomposition"

**Prompt for coding agent:**
```
Create `figures/fig_variance_decomposition.py`:

Layout: 2×2 panel figure

Top row: δ¹³C-only models (2x2_one left, 3x3_two right)
Bottom row: Dual-isotope models (2x2_one left, 3x3_two right)

Each panel: Stacked bar chart showing variance contributions to FF uncertainty:
- KIE uncertainty (red)
- Source-signature uncertainty (blue) 
- Lifetime uncertainty (green)
- Observational uncertainty (gray)

Message: In δ¹³C-only (top), KIE dominates. In dual-isotope (bottom), KIE shrinks dramatically.

Average over 2007–2022.
```

---

## Phase 5: Extended Sensitivity Tests

### Task 5.1: OH_D KIE sensitivity (the "other" KIE)

**Prompt for coding agent:**
```
Create `analysis/oh_d_kie_sensitivity.py`:

The OH_D KIE (1.294–1.327) also has uncertainty. Test whether:

1. When using dual isotopes, does OH_D KIE uncertainty REPLACE OH_13C KIE as the 
   dominant systematic? (i.e., are we just trading one KIE problem for another?)

2. Compute:
   - FF spread from OH_D KIE alone (fix OH_13C at midpoint, vary OH_D)
   - FF spread from OH_13C KIE alone (fix OH_D at midpoint, vary OH_13C)
   - FF spread from both KIEs sampled simultaneously

3. If OH_D spread << OH_13C spread: excellent — dual isotopes genuinely help
   If OH_D spread ≈ OH_13C spread: caution — we've just moved the problem

4. Physical reasoning: OH_D should matter less because the δD signal is dominated 
   by source-signature differences (which span ~200‰) rather than sink fractionation.
   Verify this quantitatively.

5. Output: comparison table + figure showing "KIE attribution" to each reaction.
```

### Task 5.2: Cl sink fraction sensitivity

**Prompt for coding agent:**
```
Create `analysis/cl_fraction_sensitivity.py`:

Riddell-Young 2025 used δD to rule out "high Cl" scenarios.
Test whether your model confirms this:

1. Run dual-isotope models with Cl fraction varied:
   - Default: Cl = 0.035 (3.5% of total sink)
   - Thanwerdas: Cl = 0.006 (0.6%)
   - High-Cl: Cl = 0.05 (5%)
   - Very-high-Cl: Cl = 0.08 (8%, extreme scenario)

2. For each Cl fraction:
   - Does the dual-isotope system produce physical (positive) solutions?
   - What % of MC iterations fail quality checks?
   - How do FF/Mic emissions change?

3. Test: Can your model reject high-Cl scenarios using δD (as Riddell-Young claims)?
   - If quality-monitor rejects >50% of iterations at Cl=0.08: δD rules it out
   - Compare quality-pass rates in δ¹³C-only vs. dual-isotope modes

4. Output: Table of quality-pass rates × Cl fraction × isotope mode
   Figure: "Solution quality vs. Cl fraction" showing dual-isotope's ability to reject
```

---

## File Structure (Expected)

```
proposals/title2_KIE_immunity/
├── SUMMARY.md
├── PLAN.md
├── analysis/
│   ├── kie_sensitivity_matrix.py
│   ├── kie_cross_validation.py
│   ├── immunity_factor.py
│   ├── inventory_kie_interaction.py
│   ├── robust_attribution.py
│   ├── compare_basu2022.py
│   ├── oh_d_kie_sensitivity.py
│   └── cl_fraction_sensitivity.py
├── figures/
│   ├── fig_kie_immunity.py
│   └── fig_variance_decomposition.py
└── results/
    ├── immunity_factors.csv
    ├── variance_decomposition.json
    ├── robust_attribution.csv
    └── basu_comparison.json
```

## Dependencies

```
numpy>=1.24
pandas>=2.0
scipy>=1.10
matplotlib>=3.7
seaborn>=0.12
```

## Estimated Runtime

- Phase 1: ~30 min (4 models × 6 configs × 1000 MC)
- Phase 2: ~20 min (bootstrap resampling)
- Phase 3: ~15 min (variance decomposition runs)
- Phase 4: ~5 min (figures)
- Phase 5: ~30 min (Cl and OH_D sensitivity)

Total: ~2 hours of compute
