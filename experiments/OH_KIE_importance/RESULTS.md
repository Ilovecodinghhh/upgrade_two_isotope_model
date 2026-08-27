# OH KIE Importance Experiment — Results

**Experiment:** `experiments/OH_KIE_importance/`
**Date:** 2026-05-17
**Data vintage:** Corrected (RNG-stream-isolated, v2)

---

## 1. Objective

Quantify how the two OH kinetic isotope effects — OH-¹³C (Saueressig 1.0039 vs Cantrell 1.0054) and OH-D (1.294–1.327) — propagate into fossil-fuel (FF) emission estimates across two model architectures:

| Architecture | Solver | Free sources | Isotope coupling |
|:---:|:---:|:---:|:---:|
| **3×3 one-box** | `np.linalg.solve` (3×3 matrix) | BB + FF + Mic | Coupled (simultaneous δ¹³C + δD) |
| **2×2 one-box** | 2×2 analytic per isotope | FF + Mic (BB fixed) | Decoupled (δ¹³C and δD solved independently) |

---

## 2. Method

### 2.1 Selective Freezing

Seven KIE configurations, each with 1000 MC iterations, identical random seeds:

| # | Config | OH-¹³C | OH-D | Other KIE | Purpose |
|:-:|--------|:------:|:----:|:---------:|---------|
| 1 | ALL_SAMPLED | sampled | sampled | sampled | Baseline σ(FF) |
| 2 | FIX_OH13C | **fixed** (1.00465) | sampled | sampled | Remove OH-¹³C uncertainty |
| 3 | FIX_OHD | sampled | **fixed** (1.3105) | sampled | Remove OH-D uncertainty |
| 4 | FIX_BOTH_OH | **fixed** | **fixed** | sampled | Remove all OH-KIE uncertainty |
| 5 | ALL_KIE_FIXED | **fixed** | **fixed** | **fixed** | Remove all KIE uncertainty |
| 6 | OH13C_SAUERESSIG | **fixed** (1.0039) | sampled | sampled | Saueressig endpoint |
| 7 | OH13C_CANTRELL | **fixed** (1.0054) | sampled | sampled | Cantrell endpoint |

Variance attribution: `% of σ²(FF) = 100 × (σ²_baseline − σ²_config) / σ²_baseline`

### 2.2 RNG Isolation (Bug Fix)

A critical bug was identified and fixed in this experiment. The original code used a single RNG stream for both KIE sampling and data/source-signature sampling. When a KIE parameter was frozen, the draw was skipped, shifting the RNG stream and causing all downstream data draws to differ between configs. This produced **spurious cross-isotope contamination** (e.g., OH-D falsely appearing to contribute 10% to δ¹³C variance).

**Fix:** Two independent RNG streams:
- `rng` (seed=42) — data, source signatures, atmospheric observations
- `rng_kie` (seed=1042) — KIE sampling only

The KIE sampler **always consumes 8 draws** regardless of freezing (draw then override), ensuring the data stream is identical across all configs.

---

## 3. Results — 3×3 Dual-Isotope

### 3.1 Summary Table

| Config | σ(FF) [Tg/yr] | Mean FF [Tg/yr] | ΔFF trend | Nonphysical % |
|--------|:--------------:|:----------------:|:---------:|:-------------:|
| ALL_SAMPLED | **62.7** | 46.5 | +1.9 ± 14.9 | 25.9% |
| FIX_OH13C | 62.6 | 46.5 | +1.9 ± 14.9 | 25.7% |
| FIX_OHD | **52.9** | 45.5 | +1.8 ± 14.7 | 23.1% |
| FIX_BOTH_OH | 52.6 | 45.5 | +1.8 ± 14.6 | 22.9% |
| ALL_KIE_FIXED | **52.1** | 45.4 | +1.8 ± 14.6 | 22.9% |
| OH13C_SAUERESSIG | 62.6 | **37.9** | +1.2 ± 14.9 | 28.9% |
| OH13C_CANTRELL | 62.7 | **55.0** | +2.5 ± 14.9 | 23.6% |

*Source: `results/summary.json`*

### 3.2 Variance Attribution

| Component | σ² reduction | % of σ²(FF) |
|-----------|:------------:|:-----------:|
| **OH-D KIE** | 1133 | **28.8%** |
| OH-¹³C KIE | 9 | 0.2% |
| Both OH combined | 1162 | 29.5% |
| All 8 KIE parameters | 1214 | **30.9%** |
| Non-OH KIE (Cl, Strat, Soil) | 52 | 1.3% |
| **Source sigs + data noise** | 2719 | **69.1%** |

**→ Figure 1:** σ(FF) bar chart by config
**→ Figure 3:** Pie chart + cumulative reduction

### 3.3 OH-¹³C Level Shift

| OH-¹³C value | Mean FF [Tg/yr] |
|:---:|:---:|
| Saueressig (1.0039) | 37.9 |
| Midpoint (1.00465) | 46.5 |
| Cantrell (1.0054) | 55.0 |

**ΔFF = 17.1 Tg/yr** — a systematic level shift with no change in spread (σ identical at ~62.6–62.7).

**→ Figure 2:** Time series + bar chart showing level shift

### 3.4 Interpretation

The 3×3 one-box system solves:

```
[  1        1        1     ] [BB ]   [ S_total      ]
[ f13_BB   f13_FF   f13_Mic] [FF ] = [ S·f13_source ]
[ fD_BB    fD_FF    fD_Mic ] [Mic]   [ S·fD_source  ]
```

The δD row has entries ~100× smaller than the mass balance row, making it the **binding constraint** in this ill-conditioned matrix (κ ~ 10⁴–10⁵). Perturbing the D-system KIE (OH-D) destabilizes this weakest row, inflating σ(FF). Perturbing OH-¹³C shifts the intermediate-scale ¹³C row, which translates the solution without inflating spread.

**→ Figure 4:** FF time series under 4 progressive freezing configs
**→ Figure 5:** Summary diagram showing the dual role

---

## 4. Results — 2×2 Separate-Isotope

### 4.1 Summary Table — δ¹³C-Derived FF

| Config | σ(FF) [Tg/yr] | Mean FF [Tg/yr] | ΔFF trend | Neg % |
|--------|:--------------:|:----------------:|:---------:|:-----:|
| ALL_SAMPLED | **30.3** | 177.4 | +10.4 ± 4.2 | 0.0% |
| FIX_OH13C | 28.4 | 177.5 | +10.4 ± 4.2 | 0.0% |
| FIX_OHD | **30.3** | 177.4 | +10.4 ± 4.2 | 0.0% |
| FIX_BOTH_OH | 28.4 | 177.5 | +10.4 ± 4.2 | 0.0% |
| ALL_KIE_FIXED | **27.7** | 177.8 | +10.4 ± 4.1 | 0.0% |
| OH13C_SAUERESSIG | 27.1 | **196.3** | +11.4 ± 4.1 | 0.0% |
| OH13C_CANTRELL | 29.7 | **158.7** | +9.4 ± 4.2 | 0.0% |

### 4.2 Summary Table — δD-Derived FF

| Config | σ(FF) [Tg/yr] | Mean FF [Tg/yr] | ΔFF trend | Neg % |
|--------|:--------------:|:----------------:|:---------:|:-----:|
| ALL_SAMPLED | **39.4** | 89.2 | +5.0 ± 9.5 | 4.7% |
| FIX_OH13C | **39.4** | 89.2 | +5.0 ± 9.5 | 4.7% |
| FIX_OHD | **33.6** | 87.9 | +4.8 ± 9.5 | 3.2% |
| FIX_BOTH_OH | 33.6 | 87.9 | +4.8 ± 9.5 | 3.2% |
| ALL_KIE_FIXED | **33.4** | 87.9 | +4.9 ± 9.5 | 3.1% |
| OH13C_SAUERESSIG | **39.4** | 89.2 | +5.0 ± 9.5 | 4.7% |
| OH13C_CANTRELL | **39.4** | 89.2 | +5.0 ± 9.5 | 4.7% |

*Source: `results/2x2_summary.json`*

### 4.3 Variance Attribution

**δ¹³C inversion:**

| Component | σ² reduction | % of σ²(FF) |
|-----------|:------------:|:-----------:|
| **OH-¹³C KIE** | 109 | **11.9%** |
| OH-D KIE | 0 | **0.0%** |
| Both OH | 109 | 11.9% |
| All KIE | 148 | 16.1% |
| Source sigs + noise | 767 | **83.9%** |

**δD inversion:**

| Component | σ² reduction | % of σ²(FF) |
|-----------|:------------:|:-----------:|
| OH-¹³C KIE | 0 | **0.0%** |
| **OH-D KIE** | 428 | **27.5%** |
| Both OH | 428 | 27.5% |
| All KIE | 440 | 28.3% |
| Source sigs + noise | 1114 | **71.7%** |

**→ Figure 6:** Side-by-side bar charts

### 4.4 Zero Cross-Isotope Contamination

The critical validation result:

| KIE parameter | → δ¹³C FF variance | → δD FF variance |
|:---:|:---:|:---:|
| **OH-¹³C** | **11.9%** | **0.0%** |
| **OH-D** | **0.0%** | **27.5%** |

Each isotope is sensitive **exclusively** to its own OH KIE. The off-diagonal entries are exactly zero (to numerical precision), confirming:
1. The RNG isolation fix works correctly
2. There is no physical pathway for cross-isotope KIE contamination in decoupled inversions
3. The previous 10.1% "OH-D → δ¹³C" contribution was entirely a bug artifact

**→ Figure 10:** Cross-contamination matrix heatmap

### 4.5 Saueressig vs Cantrell Level Shift

| Isotope | Saueressig FF | Cantrell FF | **ΔFF** |
|:---:|:---:|:---:|:---:|
| **δ¹³C** | 196.3 | 158.7 | **−37.6** |
| **δD** | 89.2 | 89.2 | **0.0** |

The δ¹³C-only inversion amplifies the Saueressig–Cantrell shift to **37.6 Tg/yr** (>2× the 3×3 shift of 17.1 Tg/yr) because there is no δD constraint to anchor the solution. The δD inversion is completely immune because OH-¹³C never enters the δD mass balance.

**→ Figure 7:** Time series panels for δ¹³C and δD
**→ Figure 9:** Grouped bar chart across architectures

---

## 5. Cross-Architecture Comparison

### 5.1 Master Summary Table

| Architecture | σ(FF) [Tg/yr] | OH-¹³C → σ² | OH-D → σ² | All KIE → σ² | Saueressig–Cantrell ΔFF |
|:---:|:---:|:---:|:---:|:---:|:---:|
| **2×2 δ¹³C-only** | 30.3 | **11.9%** | 0.0% | 16.1% | **37.6 Tg/yr** |
| **2×2 δD-only** | 39.4 | 0.0% | **27.5%** | 28.3% | **0.0 Tg/yr** |
| **3×3 dual-isotope** | 62.7 | 0.2% | **28.8%** | 30.9% | **17.1 Tg/yr** |

**→ Figure 8:** Stacked horizontal bars showing variance decomposition across all 3 architectures

### 5.2 Key Observations

1. **OH-D dominance in 3×3 is a matrix-conditioning effect.** The 3×3 system couples δ¹³C and δD through a shared linear solve. The ill-conditioned δD row is the binding constraint, so OH-D perturbations (which affect that row) dominate the spread. OH-¹³C perturbations affect the well-conditioned ¹³C row, so they shift the level but not the spread.

2. **Decoupling confirms each isotope "owns" its OH KIE.** In the 2×2, OH-¹³C drives 12% of δ¹³C FF variance and exactly 0% of δD FF variance. OH-D drives 28% of δD FF variance and exactly 0% of δ¹³C FF variance. The zero off-diagonal confirms there is no hidden cross-isotope pathway.

3. **The Saueressig–Cantrell controversy is a δ¹³C-specific problem.** The level shift is:
   - **37.6 Tg/yr** in δ¹³C-only (maximum amplification, no anchoring)
   - **17.1 Tg/yr** in 3×3 (diluted by δD constraint)
   - **0.0 Tg/yr** in δD-only (completely immune)

   This means δD provides value not just for uncertainty reduction but also for **bias mitigation**: it partially anchors the solution against the OH-¹³C ambiguity.

4. **The δ¹³C–δD FF discrepancy is large.** The 2×2 gives FF = 177 Tg/yr from δ¹³C vs 89 Tg/yr from δD — a factor of 2 difference. The 3×3 compromise is 47 Tg/yr, pulled toward the δD side because the δD row determines the solution in the ill-conditioned system. This discrepancy reflects fundamental tension between the two isotope constraints in a one-box framework.

5. **Source signatures dominate in all architectures** (69–84% of variance). KIE contributes 16–31%. Non-OH KIE (Cl, Strat, Soil) contributes only 1–4%. Lifetime was not varied in this experiment but contributes <1% based on the KIE_immunity experiment.

6. **Theoretical evidence favors Cantrell.** Melissas & Truhlar (1993) calculated k₁₂/k₁₃ = 1.005 at 273–353 K using ab initio variational transition state theory (IVTST/SCT), in excellent agreement with Cantrell (1.0054) and notably higher than Saueressig (1.0039). If the Cantrell value is correct, the 3×3 FF estimate shifts to ~55 Tg/yr (vs ~38 Tg/yr for Saueressig). Their physical decomposition also explains why OH-D has ~60× larger fractionation than OH-¹³C: for ¹³C, inverse vibrational/rotational contributions nearly cancel the normal translational/tunneling terms, producing a KIE near unity; for D, all contributions are cumulatively normal. See `ImportantReferences/Melissas1993ACP/`.

---

## 6. Figure Inventory

| # | Figure | Description | File |
|:-:|--------|-------------|------|
| 1 | 3×3 σ(FF) bar chart | 5 configs, variance reduction annotations | `fig1_variance_attribution.png` |
| 2 | OH-¹³C level shift | Saueressig vs Cantrell time series + bar chart | `fig2_oh13c_level_shift.png` |
| 3 | Variance decomposition | Pie chart + cumulative σ reduction | `fig3_variance_decomposition.png` |
| 4 | 3×3 FF time series | 4 progressive freezing configs with ±σ bands | `fig4_ff_timeseries_comparison.png` |
| 5 | Dual-role diagram | OH-D = spread, OH-¹³C = level | `fig5_dual_role_summary.png` |
| 6 | 2×2 σ(FF) bar chart | Side-by-side δ¹³C and δD panels | `fig6_2x2_variance_attribution.png` |
| 7 | 2×2 Saueressig vs Cantrell | δ¹³C panel (37 Tg shift) + δD panel (0 shift) | `fig7_2x2_saueressig_cantrell.png` |
| 8 | Grand comparison | Stacked horizontal bars, 3 architectures | `fig8_grand_comparison.png` |
| 9 | Level shift comparison | Grouped bars across architectures | `fig9_level_shift_comparison.png` |
| 10 | Cross-contamination proof | Heatmap matrix + primary vs cross pathway | `fig10_cross_contamination.png` |

All figures at 250 dpi, colorblind-safe palette.

---

## 7. Technical Details

### 7.1 KIE Distributions Sampled

| Parameter | Distribution | Range/Parameters |
|-----------|:---:|:---:|
| OH-¹³C | Uniform | [1.0039, 1.0054] |
| OH-D | Uniform | [1.294, 1.327] |
| Cl-¹³C | Normal | μ=1.066, σ=0.002 |
| Cl-D | Normal | μ=1.52, σ=0.02 |
| Strat-¹³C | Normal | μ=1.003, σ=0.001 |
| Strat-D | Normal | μ=1.179, σ=0.01 |
| Soil-¹³C | Normal | μ=1.0201, σ=0.003 |
| Soil-D | Normal | μ=1.083, σ=0.01 |

Sink fractions (global): OH=0.835, Cl=0.035, Strat=0.070, Soil=0.060

### 7.2 Frozen Midpoint Values

| Parameter | Value |
|-----------|:---:|
| OH-¹³C | 1.00465 |
| OH-D | 1.3105 |
| Cl-¹³C | 1.066 |
| Cl-D | 1.52 |
| Strat-¹³C | 1.003 |
| Strat-D | 1.179 |
| Soil-¹³C | 1.0201 |
| Soil-D | 1.083 |

### 7.3 Run Parameters

- MC iterations: 1000 per config
- Seed: 42 (data), 1042 (KIE)
- Lifetime: time-varying (He et al. 2026 parameterization)
- Smoothing: 5-year running mean
- Trend metric: mean(last 3 years) − mean(2005–2007)
- Years: 1999–2021 (23 annual steps)

---

## 8. Files

| File | Description |
|------|-------------|
| `run_experiment.py` | 3×3 one-box experiment driver |
| `run_2x2_experiment.py` | 2×2 one-box experiment driver |
| `make_all_figures.py` | Comprehensive figure generation (all 10 figures) |
| `make_figures.py` | Legacy: 3×3 figures only (superseded by make_all_figures.py) |
| `make_2x2_figures.py` | Legacy: 2×2 figures only (superseded by make_all_figures.py) |
| `results/summary.json` | 3×3 results (7 configs) |
| `results/2x2_summary.json` | 2×2 results (7 configs) |
| `results/*_timeseries.csv` | Per-year FF/Mic/BB time series for each config |
| `figures/fig{1-10}_*.png` | All 10 publication figures |

---

## 9. Conclusions

1. **OH-D KIE is the dominant source of KIE-driven FF uncertainty** across all architectures — 29% in the coupled 3×3, 28% in the δD-only 2×2. Reducing OH-D KIE uncertainty would reduce FF spread by ~10 Tg/yr in the 3×3 and ~6 Tg/yr in the δD system.

2. **OH-¹³C KIE is a pure systematic bias**, not a random uncertainty. It shifts FF by 17 Tg/yr (3×3) to 38 Tg/yr (δ¹³C-only) between the Saueressig and Cantrell values. Resolving this controversy would not reduce spread but would eliminate the largest systematic offset in isotope-based source attribution.

3. **The 2×2 decoupled results confirm zero cross-isotope KIE contamination**, validating both the experimental design and the RNG isolation fix. Each isotope responds exclusively to its own OH KIE.

4. **δD provides bias mitigation**: in the 3×3, δD halves the OH-¹³C level shift (38 → 17 Tg/yr) by providing an independent constraint that is immune to the Saueressig–Cantrell controversy.

5. **Source signatures remain the dominant uncertainty** (69–84% across architectures). Improving source-signature characterization would have 2–5× more impact on FF uncertainty than resolving any single KIE parameter.

6. **KIE importance is architecture-dependent**. Cross-study comparisons of "KIE sensitivity" must specify the inversion framework. The 3×3 amplifies OH-D; the 2×2 δ¹³C amplifies OH-¹³C.

---

*Generated: 2026-05-17. All numbers from `results/summary.json` and `results/2x2_summary.json`.*
