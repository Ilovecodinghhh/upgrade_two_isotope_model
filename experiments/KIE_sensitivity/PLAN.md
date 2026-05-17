# Coding Plan: KIE Sensitivity Experiment

## Overview

Four phases, each with a self-contained coding task suitable for a coding agent.
All scripts should import from `../../common.py` (the shared model infrastructure).
All output goes to `results/` and `figures/` subdirectories.

---

## Phase 1: Baseline δ¹³C-Only KIE Sensitivity

### Task 1.1 — Create `phase1_d13c_only.py`

**Prompt for coding agent:**

> Create a Python script `experiments/KIE_sensitivity/phase1_d13c_only.py` that measures how sensitive the δ¹³C-only source attribution is to the OH-¹³C KIE choice.
>
> **Context:**
> - Import from `../../common.py` which provides: `ModelConfig`, `LoadedData`, `load_data`, `sample_KIE`, `compute_bulk_KIE`, `compute_lifetime`, `delta_to_fraction_d13C`, `fraction_to_delta_d13C`, `sample_source_signatures`, `sample_atm_d13C`, `smooth_5yr`, `trend_change`, `pad_to_length`, `SINK_FRACTIONS_GLOBAL`, `PT`, `C13_STD`
> - Look at `../../2x2_one.py` for the existing single-isotope inversion logic (the δ¹³C branch specifically)
> - The model solves: given total source S, BB (fixed from CarbonTracker), and the isotopic source composition δ¹³C_src derived from atmospheric observations and KIE, partition into FF and Mic using:
>   ```
>   FF = (S × δ¹³C_src − δ¹³C_Mic × (S − BB) − δ¹³C_BB × BB) / (δ¹³C_FF − δ¹³C_Mic)
>   Mic = S − BB − FF
>   ```
>
> **What the script should do:**
> 1. Run the δ¹³C-only mass balance 3 times with N=1000 MC iterations each:
>    - Run A: OH_13C KIE fixed at 1.0039 (Saueressig), all other KIEs sampled
>    - Run B: OH_13C KIE fixed at 1.0054 (Cantrell), all other KIEs sampled
>    - Run C: OH_13C KIE sampled uniformly from [1.0039, 1.0054], all other KIEs sampled
> 2. For each run, use `lifetime_mode='varying'` (He 2026 time-varying τ)
> 3. Store 5-yr smoothed FF and Mic emission arrays [years × iterations] for each run
> 4. Compute trend metrics: Δ(mean 2020–2022 vs mean 2005–2007) for FF and Mic per iteration
> 5. Save results to `results/phase1_d13C_only/`:
>    - `run_A_saueressig.npz` — FF_d13C, Mic_d13C arrays + trend stats
>    - `run_B_cantrell.npz` — same
>    - `run_C_sampled.npz` — same
>    - `summary.json` — means, stds, percent positive for each run's trends
> 6. Generate a figure `figures/phase1_d13C_only_trends.png`:
>    - 2 panels (FF | Mic)
>    - Each panel: overlapping histograms of Δ(2020–2022 vs 2005–2007) for runs A, B, C
>    - Colors: A=blue, B=red, C=gray
>    - Vertical dashed line at Δ=0
>    - Legend with mean ± std for each
>
> **Key implementation notes:**
> - Use `sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))` to import from `common`
> - BASE_DIR should point to the repo root: `Path(__file__).resolve().parent.parent.parent`
> - To fix OH_13C while sampling everything else, modify the KIE dict after calling `sample_KIE(rng, 'sampled')` by overwriting `kies['OH_13C']` with the fixed value
> - Use `seed=42` for reproducibility
> - Print progress every 200 iterations
> - Use matplotlib with `Agg` backend

---

## Phase 2: Dual-Isotope Joint Inversion

### Task 2.1 — Create `phase2_dual_isotope.py`

**Prompt for coding agent:**

> Create a Python script `experiments/KIE_sensitivity/phase2_dual_isotope.py` that implements a **joint dual-isotope (δ¹³C + δD) weighted least squares inversion** and measures its sensitivity to the OH-¹³C KIE choice.
>
> **Context:**
> - Same imports from `../../common.py` as Phase 1, plus: `delta_to_fraction_dD`, `fraction_to_delta_dD`, `sample_atm_dD`, `D_STD`
> - Look at `../../2x2_one.py` for how δ¹³C and δD branches each separately solve for FF and Mic
>
> **The joint system:**
> Instead of solving δ¹³C and δD independently, combine them into a single **over-determined** system (3 equations, 2 unknowns: FF and Mic):
>
> ```
> Equation 1 (mass):     FF + Mic = S − BB
> Equation 2 (¹³C):      FF·δ¹³C_FF + Mic·δ¹³C_Mic = S·δ¹³C_src − BB·δ¹³C_BB
> Equation 3 (D):         FF·δD_FF + Mic·δD_Mic = S·δD_src − BB·δD_BB
> ```
>
> In matrix form: **A·x = b** where:
> ```
> A = [[1,              1             ],
>      [δ¹³C_FF,        δ¹³C_Mic      ],
>      [δD_FF,          δD_Mic        ]]
>
> x = [FF, Mic]^T
>
> b = [S − BB,
>      S·δ¹³C_src − BB·δ¹³C_BB,
>      S·δD_src − BB·δD_BB]
> ```
>
> Solve using `scipy.optimize.lsq_linear(A, b, bounds=(0, np.inf))` to enforce non-negativity, or `np.linalg.lstsq(A, b)` with post-hoc clamping.
>
> **Weighting:** Apply inverse-variance weights to the three equations:
> - Mass equation weight: 1.0 (exact)
> - δ¹³C equation weight: 1.0 / σ²(δ¹³C_src), where σ ≈ 1.0‰ (typical observational uncertainty)
> - δD equation weight: 1.0 / σ²(δD_src), where σ ≈ 5.0‰ (typical observational uncertainty)
>
> To apply weights, multiply each row of A and b by sqrt(weight).
>
> **What the script should do:**
> 1. Run the joint dual-isotope inversion 3 times with N=1000 MC iterations each:
>    - Run A: OH_13C = 1.0039 (Saueressig), OH_D sampled, all others sampled
>    - Run B: OH_13C = 1.0054 (Cantrell), OH_D sampled, all others sampled
>    - Run C: OH_13C sampled, OH_D sampled, all others sampled
> 2. For each year and iteration:
>    - Compute total source S from CH₄ mass balance
>    - Compute δ¹³C_src and δD_src from atmospheric observations + KIE
>    - Build the 3×2 weighted system and solve
>    - Record FF, Mic (clamped to ≥ 0)
> 3. Smooth (5-yr), compute trends, save in same format as Phase 1:
>    - `results/phase2_dual_isotope/run_A_saueressig.npz`
>    - `results/phase2_dual_isotope/run_B_cantrell.npz`
>    - `results/phase2_dual_isotope/run_C_sampled.npz`
>    - `results/phase2_dual_isotope/summary.json`
> 4. Generate `figures/phase2_dual_isotope_trends.png` (same layout as Phase 1)
>
> **Key implementation notes:**
> - Use `from scipy.optimize import lsq_linear` for bounded least squares
> - The δD_src computation follows the same logic as δ¹³C_src but using δD atmospheric data, αD = 1/KIE_D, and D fractions
> - When using `lsq_linear`, set `bounds=(np.zeros(2), np.full(2, np.inf))`
> - If `lsq_linear` returns a solution where FF+Mic > S (which shouldn't happen with mass equation), cap FF+Mic = S−BB
> - Track the residual norm from the least-squares solve as a diagnostic

---

## Phase 3: Comparison & KSR Calculation

### Task 3.1 — Create `phase3_comparison.py`

**Prompt for coding agent:**

> Create a Python script `experiments/KIE_sensitivity/phase3_comparison.py` that loads Phase 1 and Phase 2 results and produces the key comparison metrics and publication-quality figures.
>
> **What the script should do:**
>
> 1. **Load results** from `results/phase1_d13C_only/` and `results/phase2_dual_isotope/` (.npz files)
>
> 2. **Compute the KIE Sensitivity Ratio (KSR):**
>    ```
>    For FF emissions:
>      spread_single = |mean_trend(Run_B_d13C) − mean_trend(Run_A_d13C)|
>      spread_dual   = |mean_trend(Run_B_dual) − mean_trend(Run_A_dual)|
>      KSR_FF = spread_single / spread_dual
>
>    Same for Mic emissions → KSR_Mic
>    ```
>    Also compute bootstrap confidence intervals for KSR (1000 bootstrap samples):
>    - Resample the MC iteration indices (with replacement)
>    - Recompute KSR each time
>    - Report 95% CI
>
> 3. **Compute uncertainty reduction:**
>    ```
>    For each year:
>      σ_single = std across iterations (Run_C_d13C FF)
>      σ_dual   = std across iterations (Run_C_dual FF)
>      reduction_pct = (1 − σ_dual/σ_single) × 100
>    ```
>
> 4. **Statistical test:**
>    - Kolmogorov-Smirnov test comparing Run_A vs Run_B trend distributions for both methods
>    - Report D-statistic and p-value
>    - If p < 0.05 for single-isotope but p > 0.05 for dual-isotope, that's strong evidence that dual isotopes make the system insensitive to KIE choice
>
> 5. **Generate figures:**
>
>    **Figure 1 — `figures/fig1_KSR_summary.png`** (publication quality, Nature-style):
>    - 2×2 grid:
>      - (a) FF trend histograms: δ¹³C-only Saueressig (blue) vs Cantrell (red)
>      - (b) FF trend histograms: dual-isotope Saueressig (blue) vs Cantrell (red)
>      - (c) Same as (a) for Microbial
>      - (d) Same as (b) for Microbial
>    - Each panel: vertical line at 0, annotate KSR value, KS p-value
>    - Figure size: (10, 8), dpi=300, font size 10
>
>    **Figure 2 — `figures/fig2_uncertainty_timeseries.png`**:
>    - 2 panels (FF | Mic), time on x-axis
>    - Each panel: 2σ uncertainty band width for δ¹³C-only (red) vs dual-isotope (blue) over time
>    - Shaded area showing the uncertainty reduction
>
>    **Figure 3 — `figures/fig3_emission_timeseries.png`**:
>    - 2 panels (FF | Mic)
>    - Each panel: median + 2σ band for Run_C (sampled KIE) of both methods
>    - Shows how the central estimate and uncertainty change
>
> 6. **Save summary** to `results/phase3_comparison/summary.json`:
>    ```json
>    {
>      "KSR_FF": {"value": ..., "CI_95": [low, high]},
>      "KSR_Mic": {"value": ..., "CI_95": [low, high]},
>      "uncertainty_reduction_FF_pct": {"mean": ..., "by_year": [...]},
>      "uncertainty_reduction_Mic_pct": {"mean": ..., "by_year": [...]},
>      "KS_test_FF": {"single_D": ..., "single_p": ..., "dual_D": ..., "dual_p": ...},
>      "KS_test_Mic": {"single_D": ..., "single_p": ..., "dual_D": ..., "dual_p": ...}
>    }
>    ```
>
> **Key implementation notes:**
> - Use `scipy.stats.ks_2samp` for the KS test
> - For bootstrap: `rng.choice(n_iterations, size=n_iterations, replace=True)` to get index arrays
> - Matplotlib style: use `plt.style.use('seaborn-v0_8-paper')` if available, else default
> - Add panel labels (a), (b), (c), (d) in upper-left corner

---

## Phase 4: Extension to 2-Box (NH/SH)

### Task 4.1 — Create `phase4_two_box.py`

**Prompt for coding agent:**

> Create a Python script `experiments/KIE_sensitivity/phase4_two_box.py` that repeats the Phase 1–3 analysis using the hemispheric 2-box framework.
>
> **Context:**
> - The 2-box model uses a 3×3 system per isotope (see `../../3x3_two.py` for reference)
> - NH and SH are treated as separate boxes with inter-hemispheric exchange (τ_ex ~ 1.0 yr)
> - NH and SH have different sink fractions (see `SINK_FRACTIONS_NH`, `SINK_FRACTIONS_SH` in common.py)
> - CH₄ NH/SH split uses an interhemispheric gradient (80–100 ppb)
> - δ¹³C NH/SH from `data.c13_NH`, `data.c13_SH`
> - δD NH/SH: global ± DD_IH_OFFSET (6‰, Riddell-Young 2025)
> - Load data with `load_data(BASE_DIR, two_box=True)`
>
> **The 2-box joint dual-isotope system** becomes a 5×4 over-determined system:
> ```
> Unknowns: [FF_NH, Mic_NH, FF_SH, Mic_SH]
>
> Equations:
> 1. NH mass: FF_NH + Mic_NH = S_NH − BB_NH − exchange(NH→SH) + exchange(SH→NH)
> 2. SH mass: FF_SH + Mic_SH = S_SH − BB_SH − exchange(SH→NH) + exchange(NH→SH)
> 3. NH δ¹³C balance
> 4. SH δ¹³C balance (or replace 3+4 with combined hemispheric δ¹³C)
> 5. NH δD balance
> (optionally 6. SH δD balance → 6×4 system)
> ```
>
> Alternatively, solve per-hemisphere 3×2 systems (mass + ¹³C + D for each hemisphere separately), which is simpler and follows the same WLS approach as Phase 2 but applied to each hemisphere.
>
> **Recommended approach: Per-hemisphere WLS (simpler)**
> - For NH: 3×2 WLS system (mass + δ¹³C + δD) → [FF_NH, Mic_NH]
> - For SH: 3×2 WLS system (mass + δ¹³C + δD) → [FF_SH, Mic_SH]
> - Inter-hemispheric exchange is absorbed into the S_NH, S_SH source terms (as in existing 3x3_two.py)
> - Global totals: FF = FF_NH + FF_SH, Mic = Mic_NH + Mic_SH
>
> **What the script should do:**
> 1. Run per-hemisphere dual-isotope WLS for 3 KIE settings × N=1000
> 2. Also run per-hemisphere δ¹³C-only for the same 3 settings (for comparison)
> 3. Compute KSR for the hemispheric model
> 4. Compare KSR(1-box) from Phase 3 vs KSR(2-box) from this phase
> 5. Save to `results/phase4_two_box/`
> 6. Generate `figures/phase4_KSR_comparison.png`:
>    - Bar chart: KSR(1-box) vs KSR(2-box) for FF and Mic
>    - Error bars from bootstrap CI
>    - Dashed line at KSR=1 (no improvement)
>
> **Key implementation notes:**
> - Use `load_data(BASE_DIR, two_box=True)` to get NH/SH data
> - For exchange terms: `F_ex = CH4_box * PT_HEMI / tau_ex` where `tau_ex ~ N(1.0, 0.1)`
> - NH lifetime = global_tau × 0.95; SH lifetime = global_tau × 1.05
> - δD NH ≈ δD_global − 6‰; δD SH ≈ δD_global + 6‰ (from `DD_IH_OFFSET`)
> - For the δ¹³C-only 2-box comparison, use a 2×2 system per hemisphere (mass + δ¹³C)

---

## Execution Order

```
Phase 1 → Phase 2 → Phase 3 → Phase 4 → Phase 4b → Phase 5 → Phase 6
                                                                  │
                                          ┌───────────────────────┼───────────────────────┐
                                          ▼                       ▼                       ▼
                                   Phase 6b              Phase 6c                Phase 7 / 8
                                (threshold sweep)    (OSSE recovery)        (time-varying KIE,
                                                                            fine thresholds,
                                                                            temporal stability)
```

Phase 7 (time-varying OH-¹³C KIE) and Phase 8 (fine-resolution sweep +
temporal stability) extend Phase 6 / 6b and are independent of each other.

---

## Running the Experiment

```bash
# From repo root
cd experiments/KIE_sensitivity

# Phase 1: baseline
python phase1_d13c_only.py

# Phase 2: dual isotope
python phase2_dual_isotope.py

# Phase 3: comparison (requires Phase 1 + 2 outputs)
python phase3_comparison.py

# Phase 4 / 4b: two-box extension
python phase4_two_box.py
python phase4b_two_box_fixed.py

# Phase 5: weight + Cl sweep
python phase5_weight_Cl_sweep.py

# Phase 6 family: agreement framework
python phase6_agreement.py
python phase6b_threshold_sweep.py
python phase6c_OSSE.py

# Phase 7: time-varying OH-¹³C KIE
python phase7_timevarying_OH.py

# Phase 8: fine threshold sweep + temporal stability
python phase8_fine_thresholds.py
```

---

## Phase 7: Time-Varying OH-¹³C KIE

### Task 7.1 — `phase7_timevarying_OH.py`

**Goal:** Test whether the agreement-rate discriminant from Phase 6b
collapses if the *bulk* OH-¹³C KIE drifts year-by-year (e.g. because of
the [OH] / temperature changes implied by He 2026 Science).

**Scenarios (each compared at threshold = 100 Tg/yr):**
1. `const_saueressig` — fixed 1.0039 (baseline)
2. `const_cantrell` — fixed 1.0054 (baseline)
3. `drift_saueressig` — linear 1.0039 (1999) → 1.00465 (2022)
4. `drift_cantrell` — linear 1.0054 (1999) → 1.00465 (2022)
5. `convergent` — Saueressig → midpoint (most aggressive damping case)

**For each scenario** compute: overall agreement rate, per-year rate,
95% bootstrap CI, n_good iterations, KSR. Then run discriminant tests
on three pairs (constant baseline, symmetric drift, convergent vs Cantrell).

**Output:**
- `results/phase7_timevarying_OH/summary.json` — full scenario+discriminant data
- `figures/fig12_timevarying_OH.png` — 1×3: trajectories | per-year rates | discriminant bars

### Implementation notes
- Vary `oh13c_trajectory[j]` *inside* the year-loop; sample all other KIEs
  per iteration as in Phase 6b.
- Reuse the same `seed=42` so cross-scenario comparisons are paired.
- Bootstrap (n_boot=2000) over the flattened (years × iterations) agreement
  matrix, restricted to valid entries.

---

## Phase 8: Fine Threshold Sweep + Temporal Stability

### Task 8a — `phase8_fine_thresholds.py` (Part 1)

**Goal:** Resolve the threshold-vs-discriminant curve at 10 Tg/yr resolution
(Phase 6b only sampled 7 thresholds with non-uniform spacing). Bootstrap
the discriminant difference at every threshold; locate (a) the threshold
maximising KSR, (b) the threshold maximising the discriminant, and (c)
the contiguous range over which the difference is statistically significant.

### Task 8b — `phase8_fine_thresholds.py` (Part 2)

**Goal:** Test temporal stability of the agreement-rate discriminant.
Split 1999–2022 into three 8-year epochs and recompute the
Cantrell–Saueressig difference for each. If the discriminant survives in
all three (with non-overlapping bootstrap CIs), the signal is *not* an
artifact of a particular atmospheric regime.

**Epochs:**
- 1999–2006 — pre-renewed-growth plateau
- 2007–2014 — renewed growth phase (post-2007 inflection)
- 2015–2022 — post-2014 acceleration

**Output:**
- `results/phase8_fine_thresholds/summary.json`
- `figures/fig13_fine_threshold.png` — 1×3: rates | discriminant ± CI | KSR
- `figures/fig14_temporal_stability.png` — 1×2: rates by epoch | discriminant by epoch

### Implementation notes
- Threshold grid: `list(range(30, 221, 10))` (20 values, 30..220).
- For each threshold, bootstrap `n_boot=2000` to get a 95% CI on
  Δ = rate(Cantrell) − rate(Saueressig); the threshold is "significant"
  when that CI is fully above zero.
- Save MC arrays internally and slice them by epoch mask — no need to
  re-run the inversions per epoch.

---

## Phase 9: Editorial Assessment Fixes

### Task 9a — `phase9_editorial_fixes.py` (High-N + KSR bootstrap)

**Goal:** Increase MC iterations from N=1000 to N=5000. Compute bootstrap 95% CIs
on both the discriminant and the KSR at key thresholds (50, 90, 100, 150 Tg/yr).

### Task 9b — `phase9_editorial_fixes.py` (Cl sensitivity)

**Goal:** Test whether the agreement-filter discriminant (35.5 pp at T=90) survives
under Thanwerdas low-Cl (0.6%) and high-Cl (6.5%) scenarios. Phase 5 only tested
Cl sensitivity for the WLS approach.

### Task 9c — `phase9_editorial_fixes.py` (Year-agree sweep)

**Goal:** Sweep the 80% year-agreement parameter across {60%, 70%, 80%, 90%, 95%}
to verify discriminant insensitivity.

**Output:**
- `results/phase9_editorial_fixes/high_n_summary.json`
- `results/phase9_editorial_fixes/cl_sensitivity.json`
- `results/phase9_editorial_fixes/year_agreement_sweep.json`
- `figures/fig15_high_n.png`, `fig16_cl_sensitivity.png`, `fig17_year_agree_sweep.png`

### Implementation notes
- Single script (`phase9_editorial_fixes.py`) with three task functions.
- Reuses run_inversions() with parameterised `cl_fraction` and `n_iter`.
- For Cl override: adjusts OH fraction to keep sink fractions summing to 1.
- Bootstrap KSR uses 2000 bootstrap samples of MC iteration indices.

---

## Execution Order

```
Phase 1 → Phase 2 → Phase 3 → Phase 4 → Phase 4b → Phase 5 → Phase 6
                                                                  │
                                          ┌───────────────────────┼───────────────────────┐
                                          ▼                       ▼                       ▼
                                   Phase 6b              Phase 6c                Phase 7 / 8
                                (threshold sweep)    (OSSE recovery)        (time-varying KIE,
                                                                            fine thresholds,
                                                                            temporal stability)
                                          │
                                          ▼
                                    Phase 9
                              (editorial fixes:
                               N=5000, Cl sweep,
                               year-agree sweep)
```

```
numpy
scipy
pandas
matplotlib
openpyxl          # for .xlsx reading
```

All already available in the repo environment.

---

## Notes for the Coding Agent

- **DO NOT modify** `common.py`, `2x2_one.py`, or any existing model scripts
- **DO** import from `common.py` — it has everything you need
- The data lives in `../../rel/data/` and `../../rel/output/` (relative to repo root), or in a sibling `TwoIsotopeBoxModel/rel/` directory — `load_data()` handles this automatically
- Use `Path(__file__).resolve().parent.parent.parent` to get the repo root
- Create output directories with `mkdir(parents=True, exist_ok=True)`
- All figures should be saved at 300 dpi for publication quality
- Use `matplotlib.use('Agg')` at the top (no display needed)
- Target Python 3.10+ (f-strings, type hints, `|` union syntax are fine)
