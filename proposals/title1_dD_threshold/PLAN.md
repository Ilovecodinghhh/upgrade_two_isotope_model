# Title 1: Implementation Plan & Coding Agent Prompts

## Phase 1: Baseline Runs (δ¹³C-only vs. dual-isotope)

### Task 1.1: Modify v2.0 (1-box) to run in δ¹³C-only mode

**Prompt for coding agent:**
```
In the repository at https://github.com/Ilovecodinghhh/upgrade_two_isotope_model (branch: master),
modify `2x2_one.py` to add a CLI flag `--isotopes d13C_only | dual`.

When `--isotopes d13C_only`:
- Solve only the δ¹³C mass-balance equation (1×1 system: solve for FF+Mic combined, then 
  attribute using a δ¹³C-only linear system with BB fixed from CarbonTracker).
- Store results in `Output_2x2_one_d13C_only/`

When `--isotopes dual` (default, current behavior):
- Keep existing behavior unchanged.

Run both modes with `default` preset for N=1000 MC iterations.
Save results as CSV with columns: year, FF_median, FF_p5, FF_p95, Mic_median, Mic_p5, Mic_p95.
```

### Task 1.2: Modify v3.1b (2-box) for δ¹³C-only mode

**Prompt for coding agent:**
```
In `3x3_two.py`, add `--isotopes d13C_only | dual` flag.

When `--isotopes d13C_only`:
- Drop the δD row from the 3×3 system, making it a 2×2 per hemisphere.
- Solve for FF and Mic per hemisphere with BB fixed (like 2x2_two.py but hemispheric).
- Use `lsq_linear` with bounds [0, inf) for non-negativity.

When `--isotopes dual` (default):
- Keep existing 3×3 behavior.

Run both modes with `default` preset, N=1000.
Output: `Output_3x3_two_d13C_only/results.csv` with columns:
year, NH_FF_median, NH_FF_p5, NH_FF_p95, NH_Mic_median, ..., SH_FF_median, ..., Global_FF_median, ...
```

---

## Phase 2: Information Gain Quantification

### Task 2.1: Implement DFS (Degrees of Freedom for Signal) calculator

**Prompt for coding agent:**
```
Create a new script `analysis/dfs_calculator.py` that:

1. For a given model variant (1-box or 2-box) and a given set of source-signature uncertainties:
   a. Constructs the Jacobian matrix H (∂observations/∂sources) numerically using finite differences.
   b. Constructs the prior covariance matrix B from the source-signature MC spread.
   c. Constructs the observation error covariance R from measurement uncertainties.
   d. Computes DFS = trace(HBH^T (HBH^T + R)^{-1})

2. Compute DFS for four configurations:
   - 1-box, δ¹³C-only
   - 1-box, δ¹³C + δD
   - 2-box, δ¹³C-only
   - 2-box, δ¹³C + δD

3. Report ΔDFS = DFS(dual) - DFS(δ¹³C-only) for each spatial resolution.

4. Output: JSON file with DFS values and a matplotlib figure showing DFS vs. spatial resolution.

Use numpy and scipy. Import source-signature data from `rel/output/` CSVs.
Atmospheric observation uncertainties: δ¹³C = 0.04‰ (annual mean), δD = 3‰ (annual mean).
```

### Task 2.2: Bayesian information gain via PyMC

**Prompt for coding agent:**
```
Create `analysis/bayesian_information_gain.py` using PyMC (v5+):

1. Implement the 1-box mass-balance as a PyMC model with:
   - Priors: FF ~ N(100, 30), Mic ~ N(230, 50), BB ~ N(30, 15) Tg/yr
   - Source signatures drawn from their MC CSVs
   - KIE sampled: OH_13C ~ U(1.0039, 1.0054), OH_D ~ U(1.294, 1.327)
   - Likelihood: observed δ¹³C (and optionally δD) with measurement uncertainties

2. Run NUTS sampler (2000 draws, 4 chains) for:
   - Model A: δ¹³C observation only
   - Model B: δ¹³C + δD observations

3. Compute:
   - Posterior entropy for FF and Mic in each model
   - Mutual information gain = H(prior) - H(posterior)
   - WAIC or LOO-CV for model comparison
   - 95% HDI width ratio: HDI_width(Model A) / HDI_width(Model B)

4. Repeat for 2-box (NH/SH split) using hemispheric observations.

5. Output: comparison table + trace plots + HDI comparison figure.

Dependencies: pymc>=5.0, arviz, numpy, pandas, matplotlib.
Use year 2015 as a representative single-year test case first, then loop over 2000–2022.
```

---

## Phase 3: Threshold Discovery (Core Experiment)

### Task 3.1: Systematic uncertainty inflation experiment

**Prompt for coding agent:**
```
Create `analysis/dD_threshold_experiment.py`:

1. Define a grid of microbial δD uncertainty multipliers: [0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
   where 1.0 = current uncertainty (~7‰ for Mic δD in our MC draws).
   
   At each multiplier level:
   - Scale the spread of `Mic_dD_MC.csv` by the multiplier (inflate around the mean)
   - Also scale `FF_dD_GlobMC_EDGAR.csv` and `BB_dD_annual.csv` by the same multiplier
   
2. For each uncertainty level, run:
   - 2x2_one.py (1-box) with dual isotopes → get FF/Mic uncertainty bands
   - 3x3_two.py (2-box) with dual isotopes → get FF/Mic uncertainty bands
   - Compare to δ¹³C-only baseline (from Phase 1)

3. Define "δD adds value" criterion:
   - Criterion A: 95% CI width for FF emissions is ≥10% narrower with δD than without
   - Criterion B: DFS gain > 0.1 (from Phase 2 DFS calculator)

4. Find the uncertainty multiplier where both criteria fail → this is the threshold.

5. Convert multiplier back to absolute δD uncertainty (‰) for microbial sources.

6. Output:
   - CSV: multiplier, mic_dD_uncertainty_permil, FF_CI_width_dual, FF_CI_width_d13C_only, 
     DFS_gain, criterion_A_met, criterion_B_met
   - Figure: "Information gain vs. source-signature uncertainty" (line plot with threshold marked)
   - Figure: "FF emission uncertainty bands" at 3 representative uncertainty levels

Run with N=1000 MC iterations per configuration. Use `default` preset as base.
```

### Task 3.2: Replicate Thanwerdas 2024 uncertainty specification

**Prompt for coding agent:**
```
Create `analysis/thanwerdas_replication.py`:

1. Read Thanwerdas 2024 Table 2 uncertainty values:
   - Microbial (WET+AGW) δD source signature uncertainty: ~40% of -320‰ = ±128‰ (WET), 
     ~30% of -310‰ = ±93‰ (AGW)
   - FF δD: ~20% of -183‰ = ±37‰
   - BB δD: ~35% of -200‰ = ±70‰

2. Convert these to equivalent MC spread and apply to your source-signature CSVs:
   - For Mic_dD_MC.csv: inflate to match Thanwerdas's ~100‰ 1σ for microbial
   - For FF_dD: inflate to ~37‰ 1σ
   - For BB_dD: inflate to ~70‰ 1σ

3. Run your 2-box model (3x3_two.py) with these inflated uncertainties.

4. Compare:
   - Does δD still add constraint? (use criteria from Task 3.1)
   - Compare your 2-box result under Thanwerdas uncertainties to your 2-box with original uncertainties
   - Compare to your 1-box under Thanwerdas uncertainties

5. Key test: Is it the SPATIAL FRAMEWORK or the UNCERTAINTY SPECIFICATION that kills δD?
   - If 2-box + Thanwerdas uncertainties → δD useless: it's the uncertainties
   - If 2-box + Thanwerdas uncertainties → δD still helps: it's the 3D framework that dilutes the signal

6. Output: comparison table + figure showing "δD information gain" across configurations.
```

---

## Phase 4: Validation and Figure Generation

### Task 4.1: Cross-validation against published results

**Prompt for coding agent:**
```
Create `analysis/validation_against_literature.py`:

1. Load results from Riddell-Young 2025 (PNAS) Fig. 3 data (digitize if needed):
   - Their FF, Mic, BB time series with uncertainty bands
   - Their "δD rules out" scenarios (high-Cl, BB-decline)

2. Load your 1-box dual-isotope results (2x2_one default) and compare:
   - Correlation, bias, RMSE for FF and Mic time series
   - Do your uncertainty bands encompass theirs?

3. Load your 2-box results (3x3_two default) and compare:
   - Does your global aggregate (NH+SH) match Riddell-Young's global estimate?
   - Does your NH-FF match Basu 2022's hemispheric posterior? (if extractable)

4. Generate validation figure (3 panels):
   - Panel A: Your 1-box vs. Riddell-Young 2025 (should agree closely)
   - Panel B: Your 2-box global aggregate vs. Riddell-Young 2025
   - Panel C: Your 2-box hemispheric split (new result — no direct comparison available)

5. Statistical tests: Kolmogorov-Smirnov test on posterior distributions for FF (your model vs. theirs).
```

### Task 4.2: Main paper figure — "The δD Threshold"

**Prompt for coding agent:**
```
Create `figures/fig_dD_threshold.py` to generate the key paper figure:

Layout: 2×2 panel figure (Nature Communications style, 180mm wide)

Panel A (top-left): "DFS gain from δD vs. spatial resolution"
- X-axis: Model complexity (1-box, 2-box)
- Y-axis: ΔDFS (DFS_dual - DFS_d13C_only)
- Color: source-signature uncertainty level (3 representative levels)
- Message: δD adds more DFS in simpler models

Panel B (top-right): "Information gain vs. source-signature uncertainty"
- X-axis: Microbial δD source-signature uncertainty (‰, 1σ)
- Y-axis: % reduction in FF uncertainty from adding δD
- Two lines: 1-box (blue), 2-box (red)
- Vertical dashed line: threshold where δD becomes uninformative
- Gray shading: "Thanwerdas 2024 uncertainty range"
- Message: threshold exists; Thanwerdas was above it

Panel C (bottom-left): "FF emissions — δ¹³C-only vs. dual-isotope"
- Time series 2000–2022
- Shaded bands: dual-isotope (narrow) vs. δ¹³C-only (wide)
- For 2-box model at standard uncertainties
- Message: δD narrows FF constraint substantially

Panel D (bottom-right): "Source-signature precision needed for δD to help in 3D"
- Bar chart: current measurement precision vs. required precision for each source category
- Arrows showing required improvement factor
- Message: practical roadmap for the community

Style: Use seaborn + matplotlib. Colors from ColorBrewer "Set2". 
Font: 8pt for labels, 9pt for axis titles. No gridlines. 
Export as both PDF (vector) and PNG (300 dpi).
```

---

## Phase 5: Sensitivity and Robustness

### Task 5.1: KIE sensitivity under threshold conditions

**Prompt for coding agent:**
```
Create `analysis/kie_sensitivity_at_threshold.py`:

Test whether the δD threshold depends on KIE assumptions:

1. Run the threshold experiment (Task 3.1) three times:
   - KIE = Saueressig (OH_13C = 1.0039)
   - KIE = Cantrell (OH_13C = 1.0054)
   - KIE = sampled (full range)

2. For each KIE setting, find the threshold where δD stops helping.

3. Question: Does the threshold shift significantly with KIE choice?
   - If yes: report that the threshold is KIE-dependent (additional uncertainty)
   - If no: report that the threshold is robust (stronger result)

4. Output: Table of thresholds per KIE setting + overlay figure.
```

### Task 5.2: Lifetime sensitivity

**Prompt for coding agent:**
```
Create `analysis/lifetime_sensitivity_threshold.py`:

Test whether fixed vs. time-varying lifetime affects the δD threshold:

1. Run threshold experiment under:
   - τ = 9.0 yr (fixed, IPCC)
   - τ(t) = 9.0 − 0.017·(t−2010) (He 2026 trend)
   - τ = 8.5 yr (short-lifetime scenario)

2. For each, find the δD threshold.

3. Output: sensitivity table + figure.

Use the `fixed_lifetime` preset and `default` preset from inputs.py.
```

---

## File Structure (Expected)

```
proposals/title1_dD_threshold/
├── SUMMARY.md          (this file's companion)
├── PLAN.md             (this file)
├── analysis/
│   ├── dfs_calculator.py
│   ├── bayesian_information_gain.py
│   ├── dD_threshold_experiment.py
│   ├── thanwerdas_replication.py
│   ├── validation_against_literature.py
│   ├── kie_sensitivity_at_threshold.py
│   └── lifetime_sensitivity_threshold.py
├── figures/
│   └── fig_dD_threshold.py
└── results/
    ├── dfs_comparison.json
    ├── threshold_results.csv
    ├── thanwerdas_replication.csv
    └── validation_statistics.json
```

---

## Dependencies

```
numpy>=1.24
pandas>=2.0
scipy>=1.10
matplotlib>=3.7
seaborn>=0.12
pymc>=5.0
arviz>=0.15
```

## Estimated Runtime

- Phase 1: ~10 min (4 model runs × 1000 MC × ~0.15s each)
- Phase 2: ~30 min (PyMC NUTS sampling)
- Phase 3: ~2 hours (8 uncertainty levels × 2 models × 1000 MC + Thanwerdas replication)
- Phase 4: ~5 min (plotting)
- Phase 5: ~3 hours (3 KIE settings × 8 uncertainty levels × 2 models × 1000 MC)

Total: ~6 hours of compute
