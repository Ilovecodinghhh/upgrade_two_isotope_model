# Title 3: Implementation Plan & Coding Agent Prompts

## Phase 1: Core Hemispheric Model Runs

### Task 1.1: Run the full 2-box dual-isotope model (v3.1b)

**Prompt for coding agent:**
```
In the repository at https://github.com/Ilovecodinghhh/upgrade_two_isotope_model (branch: master),
run `3x3_two.py` with the `default` preset (N=1000 MC iterations).

Then modify `3x3_two.py` to output ADDITIONAL diagnostics beyond the current results.csv:

Create `Output_3x3_two/hemispheric_detail.csv` with columns:
- year
- NH_FF_median, NH_FF_p5, NH_FF_p95, NH_FF_mean, NH_FF_std
- NH_Mic_median, NH_Mic_p5, NH_Mic_p95
- NH_BB_median, NH_BB_p5, NH_BB_p95
- SH_FF_median, SH_FF_p5, SH_FF_p95
- SH_Mic_median, SH_Mic_p5, SH_Mic_p95
- SH_BB_median, SH_BB_p5, SH_BB_p95
- Global_FF_median (= NH_FF + SH_FF), Global_FF_p5, Global_FF_p95
- Global_Mic_median, Global_Mic_p5, Global_Mic_p95
- Global_BB_median, Global_BB_p5, Global_BB_p95

Also save per-iteration raw data:
`Output_3x3_two/all_iterations.npz` with arrays:
- NH_FF[1000, n_years], NH_Mic[1000, n_years], NH_BB[1000, n_years]
- SH_FF[1000, n_years], SH_Mic[1000, n_years], SH_BB[1000, n_years]
- years[n_years]

This enables post-hoc analysis without re-running.
```

### Task 1.2: Run the equivalent 1-box model for comparison

**Prompt for coding agent:**
```
Run `3x3_one.py` with `default` preset (N=1000 MC).

Save similarly detailed output:
`Output_3x3_one/global_detail.csv` with columns:
- year, FF_median, FF_p5, FF_p95, Mic_median, Mic_p5, Mic_p95, BB_median, BB_p5, BB_p95

Also save: `Output_3x3_one/all_iterations.npz` with FF[1000, n_years], Mic[1000, n_years], BB[1000, n_years]

This is the "Riddell-Young-like" reference result.
```

### Task 1.3: Run 2x2_two (fixed-BB hemispheric) for BB-sensitivity

**Prompt for coding agent:**
```
Run `2x2_two.py` with `default` preset (N=1000 MC).

Save: `Output_2x2_two/hemispheric_detail.csv` with same format as Task 1.1 
(but BB is fixed, so BB columns are the prescribed values).

Also run with `BB_declining` preset to test Worden 2017's BB-decline scenario hemisphericaly.
Save to `Output_2x2_two_BBdecline/`.
```

---

## Phase 2: Trend Analysis — Does NH ≠ SH?

### Task 2.1: Hemispheric trend extraction and significance testing

**Prompt for coding agent:**
```
Create `analysis/hemispheric_trends.py`:

1. Load `Output_3x3_two/all_iterations.npz`

2. For each source (FF, Mic, BB) and hemisphere (NH, SH):
   - Fit linear trend over 2007–2022 for EACH of the 1000 MC iterations
   - Get distribution of slopes: slope[1000]
   - Report: median slope, 90% CI, p-value (fraction of iterations with slope > 0 or < 0)

3. Also compute global aggregate trends:
   - Global_FF[i, t] = NH_FF[i, t] + SH_FF[i, t]
   - Fit trend to Global_FF → get slope distribution

4. Key comparison:
   - NH_FF trend: positive or negative? Significant?
   - SH_FF trend: positive or negative? Significant?
   - Global_FF trend: positive or negative? Significant?
   - If NH_FF positive but Global_FF ~ zero: spatial aliasing confirmed

5. Same for 1-box model (3x3_one): does its global FF trend match the 2-box global aggregate?

6. Output:
   - `results/trend_analysis.csv`: source, hemisphere, slope_median, slope_p5, slope_p95, 
     p_value_positive, p_value_negative
   - Figure: "Hemispheric FF trends" — violin plots of slope distributions for NH, SH, Global
```

### Task 2.2: Is the NH/SH divergence robust to model assumptions?

**Prompt for coding agent:**
```
Create `analysis/divergence_robustness.py`:

Run trend analysis (Task 2.1) under multiple presets:

1. `default` — full stochastic (KIE sampled + time-varying τ)
2. `fixed_lifetime` — τ = 9.0 yr
3. `cantrell_only` — OH_13C KIE = 1.0054
4. `saueressig_only` — OH_13C KIE = 1.0039
5. `thanwerdas_sinks` — low Cl fraction
6. `CTCH4_FF` — CarbonTracker FF signatures
7. `tau_ex_fixed` — fixed interhemispheric exchange (no sampling)

For each preset:
- Run 3x3_two.py (1000 MC)
- Extract NH_FF and SH_Mic trends
- Report whether the "NH stable/increasing FF, SH increasing Mic" pattern persists

Output:
- Table: preset, NH_FF_slope, NH_FF_significant, SH_Mic_slope, SH_Mic_significant, pattern_holds
- Figure: Forest plot of NH_FF slopes across all presets (with CI bars)
- Conclusion: "The hemispheric divergence is robust to [X of 7] sensitivity assumptions"
```

---

## Phase 3: Reconciliation Test — Does 2-Box Bridge 1-Box and 3D?

### Task 3.1: Compare NH-FF to Basu 2022 regional estimates

**Prompt for coding agent:**
```
Create `analysis/reconcile_with_3D.py`:

1. From Basu 2022 (ACP), extract their posterior FF emission estimates for:
   - Global total
   - NH contribution (if reported; otherwise estimate as ~70% of global based on their Figure 7)
   - Post-2007 FF trend

   Basu 2022 key numbers (from their paper):
   - Total FF: ~115 Tg/yr average (2000-2017)
   - FF increase post-2007: significant positive trend
   - They report by latitude band in some figures

2. From Riddell-Young 2025 (PNAS):
   - Global FF: relatively stable post-2007
   - Their Fig. 3 shows FF time series with uncertainty

3. Compare:
   a. Your 2-box Global_FF vs. Riddell-Young Global_FF (should agree → both show FF stable)
   b. Your 2-box NH_FF vs. Basu's NH estimates (should agree → both show NH-FF increase)
   c. Your 1-box (3x3_one) Global_FF vs. your 2-box Global_FF (should agree → aggregation works)

4. Quantify the "aliasing effect":
   - aliasing_bias = mean(NH_FF_trend + SH_FF_trend) vs. trend(NH_FF + SH_FF)
   - These should be identical mathematically, but the INTERPRETATION changes:
     - "FF is stable globally" (1-box view) vs. "NH-FF rises, SH-FF declines" (2-box view)

5. Generate reconciliation diagram (conceptual figure for paper):
   - Show how 1-box, 2-box, and 3D all see different aspects of the same truth

Output:
- Comparison table with literature values and your estimates
- Figure: "Reconciliation" — your results bridging between Riddell-Young and Basu
```

### Task 3.2: Validate SH-Microbial against tropical studies

**Prompt for coding agent:**
```
Create `analysis/validate_sh_microbial.py`:

1. Your model's SH_Mic trend should be dominated by tropical wetlands + livestock.
   Compare against:

   a. Zhang 2021 (Nat Comm): GOSAT-derived tropical emission increases 2010-2018
      - Their wetland increase: ~5 Tg/yr over 2010-2018
      - Their livestock increase: ~5 Tg/yr over 2010-2018
      
   b. Chandra 2024 (Comm Earth): Post-2019 tropical wetland surge
      - La Niña years (2020-2022) drove record wetland emissions
      
   c. Dasgupta 2025 (EGU): Tropical wetland vs. agriculture separation
      - Their tropical total microbial increase

2. Your SH_Mic includes: wetlands + livestock + rice + waste (all microbial-signature sources)
   Sum of Zhang's tropical wetland + livestock ≈ your SH_Mic? (rough validation)

3. Also check: Does your SH_BB match GFED fire data for SH?
   - SH fires are dominated by South America, Africa, SE Asia
   - Should show interannual variability (El Niño/La Niña)

4. Output:
   - Comparison figure: Your SH_Mic time series overlaid with Zhang/Chandra estimates
   - Correlation analysis
   - Note: imperfect match expected (SH ≠ tropics exactly; box model is simpler)
```

---

## Phase 4: The Interhemispheric Exchange Test

### Task 4.1: Sensitivity to exchange rate

**Prompt for coding agent:**
```
Create `analysis/exchange_rate_sensitivity.py`:

The interhemispheric exchange rate (τ_ex) controls how much NH and SH "communicate."
If τ_ex is too fast, hemispheres equilibrate and you lose information.
If τ_ex is too slow, you overestimate hemispheric differences.

1. Run 3x3_two.py with τ_ex varied:
   - τ_ex = 0.8 yr (fast exchange — hemispheres nearly mixed)
   - τ_ex = 1.0 yr (default — Patra et al. 2011)
   - τ_ex = 1.2 yr (slow exchange — more hemispheric independence)
   - τ_ex = 1.5 yr (extreme slow — maximum hemispheric divergence)
   - τ_ex = 2.0 yr (unphysically slow — sanity check)

2. For each τ_ex:
   - Extract NH_FF and SH_Mic trends
   - Compute the NH-SH gradient in δ¹³C and δD (model-predicted vs. observed)
   - Does the model reproduce the observed interhemispheric δ¹³C gradient? (~0.2-0.3‰)

3. Use the δ¹³C gradient as a constraint:
   - Which τ_ex values produce a realistic NH-SH δ¹³C gradient?
   - Report: "τ_ex is constrained to [X, Y] yr by the interhemispheric gradient"
   - Within this range, how much does the hemispheric divergence conclusion change?

4. Output:
   - Figure: "NH_FF trend vs. τ_ex" and "Modeled vs. observed IH gradient"
   - Table: τ_ex, IH_gradient_d13C, IH_gradient_dD, NH_FF_trend, SH_Mic_trend
   - Conclusion on robustness
```

### Task 4.2: NH/SH δD gradient (novel prediction)

**Prompt for coding agent:**
```
Create `analysis/dD_gradient_prediction.py`:

Your 2-box model predicts a NH/SH δD-CH₄ gradient that has NEVER been published.
(δD observations are too sparse to directly observe a reliable hemispheric gradient.)

1. From your 3x3_two model results, compute:
   - Predicted δD_NH - δD_SH gradient for each year (2000-2022)
   - What sign and magnitude does the model predict?
   
2. The predicted gradient depends on:
   - Source distribution (more FF in NH → higher δD; more Mic in SH → lower δD)
   - Sink fractionation differences (NH Cl is higher → larger δD fractionation)
   - Exchange rate

3. Compare to the ±6‰ offset prescribed in your `inputs.py` (dD_offset parameter).
   - Is this consistent with what the model predicts should exist?
   - If your model predicts a different gradient, which is correct?

4. This becomes a TESTABLE PREDICTION for future observations:
   - "Our model predicts that the NH-SH δD-CH₄ gradient should be X ± Y ‰"
   - "This prediction can be tested as the δD monitoring network expands"

5. Output:
   - Time series of predicted NH-SH δD gradient
   - Discussion of observational requirements to test this
```

---

## Phase 5: Main Figures

### Task 5.1: Figure 1 — "Hemispheric source divergence"

**Prompt for coding agent:**
```
Create `figures/fig_hemispheric_divergence.py`:

Layout: 3-row × 2-column panel figure (Nature Communications format, 180mm wide)

Row 1: NH sources
- Panel A (left): NH Fossil Fuel emissions (Tg/yr) with uncertainty band, 2000-2022
- Panel B (right): NH Microbial emissions with uncertainty band

Row 2: SH sources
- Panel C (left): SH Fossil Fuel emissions
- Panel D (right): SH Microbial emissions

Row 3: Global aggregate + comparison
- Panel E (left): Global FF (= NH + SH) with uncertainty
  - Overlay: Riddell-Young 2025 estimate (dashed, from their Fig. 3)
  - Overlay: Basu 2022 estimate (dotted)
- Panel F (right): Global Mic with uncertainty
  - Overlay: same literature comparisons

Key visual: NH_FF shows slight increase (Panel A), SH_Mic shows strong increase (Panel D),
but Global_FF (Panel E) appears stable — demonstrating the aliasing effect.

Add trend lines (linear fits) with slope ± SE annotated.
Color scheme: FF = warm (orange/red), Mic = cool (blue/teal), BB = gray
Uncertainty: 90% CI as shaded band, median as solid line.

Style: Clean, minimal. No gridlines. 8pt labels.
Export PDF + PNG (300 dpi).
```

### Task 5.2: Figure 2 — "Reconciliation schematic"

**Prompt for coding agent:**
```
Create `figures/fig_reconciliation.py`:

A conceptual/data figure showing HOW the reconciliation works.

Layout: Single wide panel (or 1×3)

Panel A: "What 1-box models see"
- Global FF time series (flat) with arrow saying "FF stable"
- Citation: Riddell-Young 2025, He 2026

Panel B: "What our 2-box reveals"  
- NH-FF (slight increase, red) + SH-FF (slight decrease, blue) = Global (flat, black)
- Arrow: "Spatial aliasing hides NH trend"

Panel C: "What 3D inversions see"
- Regional FF showing NH hotspots (conceptual map or regional bars)
- Citation: Basu 2022, Thanwerdas 2024
- Arrow: "Regional increases in NH"

Bottom: Summary arrow connecting all three:
"1-box, 2-box, and 3D all correct — viewing different spatial scales of same phenomenon"

This is the paper's conceptual contribution — make it visually clear and compelling.
```

### Task 5.3: Figure 3 — "Exchange rate constraint and robustness"

**Prompt for coding agent:**
```
Create `figures/fig_robustness.py`:

Layout: 2-panel figure

Panel A: "Interhemispheric δ¹³C gradient constrains τ_ex"
- X-axis: τ_ex (yr)
- Y-axis: Modeled NH-SH δ¹³C gradient (‰)
- Horizontal band: Observed gradient (from data ± uncertainty)
- Vertical shading: Acceptable τ_ex range
- Message: Exchange rate is observationally constrained

Panel B: "Hemispheric divergence robust across assumptions"
- Forest plot: NH_FF trend slope ± CI for each sensitivity preset
  (default, fixed_lifetime, cantrell, saueressig, thanwerdas_sinks, CTCH4_FF, tau_ex_fixed)
- Vertical zero line
- Mark significant trends with stars
- Message: NH-FF positive trend persists across most assumptions
```

---

## Phase 6: Extended Analysis

### Task 6.1: Post-2019 acceleration — hemispheric decomposition

**Prompt for coding agent:**
```
Create `analysis/post2019_acceleration.py`:

Chandra 2024 showed a surge in CH₄ growth post-2019 (La Niña tropical wetlands).
Does your hemispheric model attribute this correctly?

1. Focus on 2019–2022 period:
   - Compute year-on-year change in NH_FF, NH_Mic, SH_FF, SH_Mic
   - Is the 2020-2022 acceleration primarily SH-Mic? (matching Chandra)

2. Compare magnitudes:
   - Chandra 2024 tropical wetland surge: ~10-15 Tg/yr above 2010-2018 average
   - Your SH_Mic 2020-2022 anomaly: should be similar order of magnitude

3. Test: Did NH sources also change post-2019?
   - COVID-19 (2020): potential FF dip due to reduced fossil fuel activity
   - Does NH_FF show a 2020 dip followed by 2021-2022 rebound?

4. Output:
   - Figure: Year-on-year source changes, split by hemisphere (bar chart)
   - Table: Anomalies relative to 2010-2018 mean for each source × hemisphere
```

### Task 6.2: Extend to 3-box (NH-extratropical / Tropical / SH-extratropical)

**Prompt for coding agent:**
```
Create `analysis/three_box_extension.py`:

As a forward-looking analysis (for Discussion section), sketch a 3-box extension:

1. Define three boxes:
   - NH extratropical (30°N–90°N): dominated by FF + boreal wetlands
   - Tropical (30°S–30°N): dominated by tropical wetlands + livestock + some FF
   - SH extratropical (90°S–30°S): minimal sources, clean background

2. Set up the mass-balance system:
   - 5×5 per box? Or solve each box's 3×3 with exchange between adjacent boxes?
   - Exchange rates: NH-Trop = faster (~0.8 yr), Trop-SH = slower (~1.2 yr)
   
3. This is EXPLORATORY — don't need full implementation, just:
   - Write the mathematical framework (equations in comments)
   - Identify what additional data would be needed (tropical δ¹³C/δD stations)
   - Estimate feasibility given current observation network

4. Key question: Would a 3-box model further resolve the tropical wetland signal?

5. Output:
   - Mathematical framework document (LaTeX-formatted equations in markdown)
   - Data requirements assessment
   - Feasibility conclusion (for paper's Discussion/Future Work section)
```

---

## File Structure (Expected)

```
proposals/title3_hemispheric_divergence/
├── SUMMARY.md
├── PLAN.md
├── analysis/
│   ├── hemispheric_trends.py
│   ├── divergence_robustness.py
│   ├── reconcile_with_3D.py
│   ├── validate_sh_microbial.py
│   ├── exchange_rate_sensitivity.py
│   ├── dD_gradient_prediction.py
│   ├── post2019_acceleration.py
│   └── three_box_extension.py
├── figures/
│   ├── fig_hemispheric_divergence.py
│   ├── fig_reconciliation.py
│   └── fig_robustness.py
└── results/
    ├── trend_analysis.csv
    ├── robustness_table.csv
    ├── reconciliation_comparison.csv
    └── exchange_rate_sensitivity.csv
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

- Phase 1: ~15 min (3 model runs × 1000 MC)
- Phase 2: ~20 min (trend fitting + robustness with 7 presets)
- Phase 3: ~10 min (comparison analysis)
- Phase 4: ~30 min (exchange rate: 5 values × 1000 MC)
- Phase 5: ~10 min (figures)
- Phase 6: ~20 min (post-2019 + 3-box sketch)

Total: ~2 hours of compute
