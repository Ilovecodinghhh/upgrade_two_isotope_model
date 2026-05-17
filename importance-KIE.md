# importance-KIE: Relative Importance of OH-¹³C vs OH-D KIE in One-Box Dual-Isotope Source Partitioning

## 1. Introduction

The hydroxyl radical (OH) is the dominant CH₄ sink (~84% of total loss), and the kinetic isotope effects (KIE) for the CH₄ + OH reaction — separately for ¹³C/¹²C and D/H — are critical but contested parameters. In the ¹³C system, two laboratory values bracket a range that has dominated the methane budget debate for two decades. In the D/H system, three measurements cluster within ~1.5% relative spread.

**This report quantifies how OH-¹³C and OH-D KIE uncertainties actually propagate into fossil-fuel (FF) emission estimates** using both a one-box dual-isotope (3×3) experiment and a one-box decoupled (2×2) experiment, each with 7000 Monte Carlo inversions.

Key finding: **In the coupled 3×3, OH-D KIE dominates the FF uncertainty spread (28% of variance), while OH-¹³C KIE controls the FF level (17 Tg/yr shift). In the decoupled 2×2, each isotope is sensitive to its own OH KIE: OH-¹³C drives 18% of δ¹³C-derived FF variance, OH-D drives 27% of δD-derived FF variance. The Saueressig–Cantrell level shift is amplified to 38 Tg/yr in the δ¹³C-only inversion but vanishes entirely in the δD-only inversion.**

---

## 2. The OH KIE Parameters

### 2.1 OH-¹³C KIE: The Saueressig–Cantrell Controversy

| Study | KIE^C_OH | ε_OH = (α−1)×1000 |
|-------|:--------:|:------------------:|
| Saueressig et al. (2001) | **1.0039** | 3.9‰ |
| Cantrell et al. (1990) | **1.0054** | 5.4‰ |

**Δε = 1.5‰** — a 38% relative difference in fractionation factor.

This repository samples OH-¹³C as U(1.0039, 1.0054) (`common.py` line 104).

### 2.2 OH-D KIE: Tighter Constraint, Larger Absolute Fractionation

| Study | KIE^D_OH |
|-------|:--------:|
| Gierczak et al. (1997) | **1.292** |
| Saueressig et al. (2001) | **1.294** |
| Joelsson et al. (2016) | **1.311** |

This repository samples OH-D as U(1.294, 1.327) (`common.py` line 105).

The absolute fractionation is ~60× larger for D (ε_D ≈ 300‰ vs ε_C ≈ 5‰), meaning OH-D has far more leverage on the δD mass balance.

---

## 3. Experiment Design

**Location:** `experiments/OH_KIE_importance/`

We ran the one-box 3×3 dual-isotope model (`3x3_one.py` architecture) under **7 KIE configurations**, each with 1000 MC iterations, time-varying lifetime, and identical random seeds:

| Config | OH-¹³C | OH-D | Other KIE | Purpose |
|--------|:------:|:----:|:---------:|---------|
| ALL_SAMPLED | sampled | sampled | sampled | Baseline |
| FIX_OH13C | **fixed** (1.00465) | sampled | sampled | Remove OH-¹³C uncertainty |
| FIX_OHD | sampled | **fixed** (1.3105) | sampled | Remove OH-D uncertainty |
| FIX_BOTH_OH | **fixed** | **fixed** | sampled | Remove all OH uncertainty |
| ALL_KIE_FIXED | **fixed** | **fixed** | **fixed** | Remove all KIE uncertainty |
| OH13C_SAUERESSIG | **fixed** (1.0039) | sampled | sampled | Saueressig endpoint |
| OH13C_CANTRELL | **fixed** (1.0054) | sampled | sampled | Cantrell endpoint |

All configs sample source signatures and atmospheric observations identically. Differences in σ(FF) are attributable solely to the KIE treatment.

---

## 4. Results

### 4.1 Summary Table

| Config | σ(FF) [Tg/yr] | Mean FF [Tg/yr] | ΔFF trend [Tg/yr] | Nonphysical % |
|--------|:--------------:|:----------------:|:------------------:|:-------------:|
| ALL_SAMPLED | **61.9** | 46.3 | +1.9 ± 14.8 | 25.6% |
| FIX_OH13C | 62.0 | 45.8 | +1.9 ± 14.9 | 25.5% |
| FIX_OHD | **52.4** | 44.9 | +1.8 ± 14.6 | 23.6% |
| FIX_BOTH_OH | 53.0 | 45.5 | +1.8 ± 14.6 | 22.9% |
| ALL_KIE_FIXED | **52.1** | 45.4 | +1.8 ± 14.6 | 22.9% |
| OH13C_SAUERESSIG | 62.0 | **37.1** | +1.2 ± 14.9 | 28.8% |
| OH13C_CANTRELL | 62.1 | **54.4** | +2.5 ± 14.9 | 23.4% |

*Source: `experiments/OH_KIE_importance/results/summary.json`*

### 4.2 Variance Attribution

Using the σ² reduction method (baseline variance minus config variance, divided by baseline):

| Component | σ² reduction | % of baseline σ² |
|-----------|:------------:|:-----------------:|
| **OH-D KIE** | **1087** | **28.4%** |
| OH-¹³C KIE | −16 | **< 1%** (noise) |
| Both OH combined | 1027 | 26.8% |
| All 8 KIE parameters | 1113 | **29.0%** |
| Non-OH KIE (Cl, Strat, Soil) | 86 | 2.2% |
| Source signatures + data noise | 2719 | **71.0%** |

**OH-D alone accounts for 28.4% of FF variance** — nearly all of the KIE-attributable uncertainty. OH-¹³C contributes effectively nothing to the spread.

![Variance attribution bar chart](figures/fig1_variance_attribution.png)
*Figure 1: σ(FF) under each KIE configuration. Fixing OH-D reduces σ by 9.5 Tg/yr (61.9 → 52.4). Fixing OH-¹³C changes nothing.*

### 4.3 OH-¹³C: Level Shift, Not Spread

While OH-¹³C doesn't affect the *uncertainty*, it controls the *level*:

| OH-¹³C value | Mean FF [Tg/yr] |
|:------------:|:----------------:|
| Saueressig (1.0039) | **37.1** |
| Midpoint (1.00465) | 45.8 |
| Cantrell (1.0054) | **54.4** |

**ΔFF = 17.3 Tg/yr** between the two endpoints — a systematic shift without changing the spread.

This occurs because in the 3×3 system, the ¹³C row determines the *position* of the FF–Mic partition (more fractionation → more apparent fossil contribution), while the D row provides the *constraint power* that resolves the three-source system.

![OH-13C level shift](figures/fig2_oh13c_level_shift.png)
*Figure 2: (a) FF time series under Saueressig vs Cantrell KIE — parallel tracks separated by ~17 Tg/yr with identical spread. (b) Mean FF ± σ for the three OH-¹³C values.*

### 4.4 Variance Decomposition Pie Chart

![Variance decomposition](figures/fig3_variance_decomposition.png)
*Figure 3: (a) Pie chart of σ²(FF) attribution. OH-D dominates the KIE contribution. Source signatures and data noise account for 71%. (b) Cumulative σ reduction by progressively fixing KIE parameters.*

### 4.5 FF Time Series Comparison

![FF time series comparison](figures/fig4_ff_timeseries_comparison.png)
*Figure 4: FF emissions (smoothed, ±1σ bands) under four configurations. The orange band (fix OH-D) is visibly narrower than the blue baseline. Fixing additional KIE parameters beyond OH-D gives minimal further reduction.*

---

## 5. 2×2 One-Box Results: Separate δ¹³C and δD Inversions

To disentangle how OH KIE propagates through each isotope system independently, we repeated the same 7-config experiment using the **2×2 one-box model** (`2x2_one.py`), where BB is fixed from CarbonTracker and FF+Mic are solved separately for δ¹³C and δD.

### 5.1 Summary Table — δ¹³C-Derived FF

| Config | σ(FF) [Tg/yr] | Mean FF [Tg/yr] | ΔFF trend [Tg/yr] | Neg % |
|--------|:--------------:|:----------------:|:------------------:|:-----:|
| ALL_SAMPLED | **31.0** | 177.7 | +10.4 ± 4.2 | 0.0% |
| FIX_OH13C | 28.1 | 178.1 | +10.4 ± 4.2 | 0.0% |
| FIX_OHD | 29.4 | 178.4 | +10.4 ± 4.2 | 0.0% |
| FIX_BOTH_OH | 28.1 | 177.9 | +10.4 ± 4.1 | 0.0% |
| ALL_KIE_FIXED | **27.7** | 177.8 | +10.4 ± 4.1 | 0.0% |
| OH13C_SAUERESSIG | 26.8 | **196.9** | +11.4 ± 4.1 | 0.0% |
| OH13C_CANTRELL | 29.4 | **159.4** | +9.4 ± 4.2 | 0.0% |

### 5.2 Summary Table — δD-Derived FF

| Config | σ(FF) [Tg/yr] | Mean FF [Tg/yr] | ΔFF trend [Tg/yr] | Neg % |
|--------|:--------------:|:----------------:|:------------------:|:-----:|
| ALL_SAMPLED | **39.1** | 88.8 | +5.0 ± 9.4 | 4.5% |
| FIX_OH13C | 38.8 | 88.9 | +4.9 ± 9.6 | 4.5% |
| FIX_OHD | **33.4** | 87.9 | +4.8 ± 9.5 | 3.1% |
| FIX_BOTH_OH | 33.6 | 88.0 | +4.8 ± 9.5 | 3.2% |
| ALL_KIE_FIXED | **33.4** | 87.9 | +4.9 ± 9.5 | 3.1% |
| OH13C_SAUERESSIG | 38.8 | 88.9 | +4.9 ± 9.6 | 4.5% |
| OH13C_CANTRELL | 38.8 | 88.9 | +4.9 ± 9.6 | 4.5% |

### 5.3 Variance Attribution — 2×2

**δ¹³C inversion:**

| Component | σ² reduction | % of baseline σ² |
|-----------|:------------:|:-----------------:|
| **OH-¹³C KIE** | 173 | **18.0%** |
| OH-D KIE | 97 | 10.1% |
| Both OH combined | 175 | 18.2% |
| All 8 KIE parameters | 194 | **20.2%** |
| Source signatures + data noise | 767 | **79.8%** |

**δD inversion:**

| Component | σ² reduction | % of baseline σ² |
|-----------|:------------:|:-----------------:|
| OH-¹³C KIE | 27 | **1.7%** |
| **OH-D KIE** | 421 | **27.4%** |
| Both OH combined | 400 | 26.1% |
| All 8 KIE parameters | 419 | **27.3%** |
| Source signatures + data noise | 1113 | **72.7%** |

![2×2 variance attribution](figures/fig6_2x2_variance_attribution.png)
*Figure 6: σ(FF) under each KIE config for the 2×2 inversion. Left: δ¹³C-derived FF shows OH-¹³C as primary KIE driver (18%). Right: δD-derived FF shows OH-D as dominant (27.4%).*

### 5.4 Saueressig vs Cantrell in 2×2

The OH-¹³C controversy produces starkly different impacts depending on which isotope constrains FF:

| Isotope | Saueressig FF [Tg/yr] | Cantrell FF [Tg/yr] | **ΔFF** |
|---------|:---------------------:|:-------------------:|:-------:|
| **δ¹³C** | 196.9 | 159.4 | **37.5** |
| **δD** | 88.9 | 88.9 | **0.0** |

The δ¹³C inversion shows a massive **37.5 Tg/yr level shift** — more than double the 17 Tg shift in the 3×3 system (which dilutes it across three equations). The δD inversion is completely immune because OH-¹³C never enters the δD mass balance.

![2×2 Saueressig vs Cantrell](figures/fig7_2x2_saueressig_cantrell.png)
*Figure 7: (a) δ¹³C-derived FF under Saueressig vs Cantrell — 38 Tg/yr separation. (b) δD-derived FF — perfectly overlapping, zero sensitivity to OH-¹³C.*

### 5.5 Grand Comparison Across Architectures

![Grand comparison](figures/fig8_grand_comparison.png)
*Figure 8: Variance attribution stacked bar across three architectures — 2×2 δ¹³C-only, 2×2 δD-only, and 3×3 dual-isotope. Each architecture has a different KIE vulnerability profile.*

![Level shift comparison](figures/fig9_level_shift_comparison.png)
*Figure 9: OH-¹³C level shift (Saueressig→Cantrell ΔFF) across architectures. Decoupled δ¹³C shows the largest shift (38 Tg); δD shows zero; 3×3 shows the intermediate coupled effect (17 Tg).*

### 5.6 Key 2×2 Insights

1. **Each isotope is sensitive to its own OH KIE** — OH-¹³C dominates δ¹³C variance (18%), OH-D dominates δD variance (27%). This is physically intuitive: in the decoupled 2×2 system, cross-isotope KIE contamination is minimal.

2. **The δ¹³C-only inversion has an unexpectedly large OH-D sensitivity (10.1%)** — this arises because OH-D affects the bulk KIE calculation through the sink-weighted average, which feeds back into the source δ¹³C via the total fractionation budget. This cross-talk disappears in the δD system (OH-¹³C → δD: only 1.7%).

3. **The Saueressig–Cantrell shift is amplified in δ¹³C-only** (38 Tg/yr vs 17 Tg in 3×3) because the 2×2 δ¹³C solve has no δD constraint to anchor the solution. The 3×3 dilutes the ¹³C bias across three equations.

4. **δD-derived FF is much lower than δ¹³C-derived FF** (89 vs 178 Tg/yr). This large discrepancy between isotope-specific estimates highlights a fundamental tension: the two isotope systems do not agree on FF magnitude in a one-box framework. The 3×3 compromise (46 Tg/yr) is even lower, suggesting the coupled system is pulled toward the δD constraint.

5. **BB prescription eliminates ~25% of 3×3 nonphysical solutions** — the 2×2 δ¹³C gives 0% negative FF (vs 25.6% in 3×3), confirming that the ill-conditioning of the 3×3 matrix (not just KIE uncertainty) drives the nonphysical fraction.

*Source: `experiments/OH_KIE_importance/results/2x2_summary.json`*

---

## 6. Why OH-D Dominates in the 3×3 (and Each Isotope "Owns" Its KIE in the 2×2)

### 6.1 The Condition Number Argument

The 3×3 system solves:

```
[ 1         1         1      ] [BB ]   [S_total          ]
[ f13_BB    f13_FF    f13_Mic] [FF ] = [S_total × f13_src]
[ fD_BB     fD_FF     fD_Mic ] [Mic]   [S_total × fD_src ]
```

The δD row entries (fD) are ~100× smaller in absolute scale than the mass balance row, making the matrix ill-conditioned (mean κ ≈ 10⁴–10⁵ from `3x3_one.py` quality monitor). Because the D row is the *weakest constraint*, any perturbation to the D-system KIE (which changes fD_src) directly destabilizes the solution. OH-D KIE uncertainty thus propagates as *spread*.

In contrast, OH-¹³C KIE perturbs the ¹³C row, which is intermediate in scale between mass balance and D. This row is better conditioned, so KIE perturbations shift the entire solution (level) without inflating the spread.

### 6.2 Fractionation Budget Comparison

| System | Source-to-atm shift | KIE^OH uncertainty (Δε) | Δε / shift |
|--------|:-------------------:|:-----------------------:|:----------:|
| ¹³C | ~6‰ | 1.5‰ | **25%** |
| D | ~230‰ | ~33‰ | **14%** |

Despite OH-¹³C having a larger *relative* KIE uncertainty (25% vs 14% of the fractionation budget), the one-box 3×3 system is more sensitive to OH-D because the δD row is the binding constraint in the ill-conditioned matrix.

### 6.3 Comparison with Two-Box Results

The `KIE_immunity` experiment (two-box, WLS formulation) found KIE accounts for ~25% of FF variance, with OH-¹³C driving the Saueressig/Cantrell sign-change controversy. The difference from our one-box result arises because:

1. **The two-box uses a WLS formulation** (not a direct 3×3 solve), where the W matrix explicitly weights δ¹³C vs δD constraints
2. **Hemispheric resolution** introduces additional δ¹³C sensitivity through the NH–SH gradient
3. **The two-box solves FF and Mic only** (BB is prescribed), changing the problem dimensionality

In the one-box 3×3, the δD row carries more relative weight because all three sources are free, and the δD row provides the decisive third equation.

### 6.4 Why the 2×2 Confirms the Physical Intuition

The 2×2 results provide crucial validation: when isotopes are solved independently, each system is sensitive primarily to **its own** OH KIE. OH-¹³C accounts for 18% of δ¹³C-derived FF variance but only 1.7% of δD-derived FF variance. OH-D accounts for 27.4% of δD-derived FF variance but only 10.1% of δ¹³C-derived FF variance.

The asymmetry (10.1% cross-talk from OH-D into δ¹³C, but only 1.7% from OH-¹³C into δD) arises because OH-D has a larger absolute effect on the bulk sink fractionation. When OH-D varies by Δε_D ≈ 33‰ (weighted by f_OH = 0.835), this perturbs the total removal rate, which has second-order effects on the ¹³C mass balance through the total source magnitude. OH-¹³C's Δε_C ≈ 1.5‰ has negligible second-order effect on the D balance.

---

## 7. Source Signature Context

From `rel/data/*_MC.csv` (1000 MC iterations, hemispheric):

| Source | δ¹³C MC std (‰) | δD MC std (‰) |
|--------|:----------------:|:--------------:|
| Fossil fuel (NH) | 2.4 | 5.6 |
| Fossil fuel (SH) | 2.7 | 8.1 |
| Microbial (NH) | 1.1 | 7.8 |
| Microbial (SH) | 1.1 | 7.3 |
| BB (NH) | 2.7 | 8.2 |
| BB (SH) | 2.4 | 7.1 |

Source signature uncertainties are 3–7× larger for δD than δ¹³C, contributing to the 71% "source sigs + data noise" share of FF variance.

---

## 8. Atmospheric Data Noise

From `rel/data/` MC ensembles:

| Observable | Annual MC std | Trend (2005–2023) | SNR/yr |
|-----------|:------------:|:-----------------:|:------:|
| δ¹³C (global) | 0.006–0.021‰ | −0.035‰/yr | ~2–6 |
| δD (global) | 0.41–0.92‰ | −0.30‰/yr | ~0.3–0.7 |

δD has ~50× larger measurement noise than δ¹³C, consistent with the finding that the D-system is the binding constraint in the one-box inversion.

---

## 9. Summary: Architecture-Dependent OH KIE Impacts

![Dual-role summary](figures/fig5_dual_role_summary.png)
*Figure 5: OH-D KIE controls FF spread (uncertainty); OH-¹³C KIE controls FF level (systematic bias).*

| Property | OH-D KIE | OH-¹³C KIE |
|----------|:--------:|:-----------:|
| **Experimental range** | U(1.294, 1.327) | U(1.0039, 1.0054) |
| **Fractionation (ε)** | ~300‰ | ~5‰ |
| **3×3 effect on σ(FF)** | **28% of variance** | **< 1%** |
| **3×3 effect on mean FF** | < 2 Tg/yr | **17 Tg/yr** |
| **2×2 δ¹³C effect on σ(FF)** | 10% (cross-talk) | **18%** |
| **2×2 δD effect on σ(FF)** | **27%** | < 2% |
| **2×2 δ¹³C Saueressig→Cantrell ΔFF** | — | **38 Tg/yr** |
| **2×2 δD Saueressig→Cantrell ΔFF** | — | **0 Tg/yr** |
| **Nature of impact** | Spread (stochastic) | Level shift (systematic) |

### Key Takeaways

1. **OH-D KIE is the dominant KIE uncertainty source** in the one-box 3×3 inversion (28% of variance) and in the δD-only 2×2 inversion (27%).

2. **OH-¹³C KIE controls the FF level** — shifting it by 17 Tg/yr in the 3×3, and by a massive 38 Tg/yr in the δ¹³C-only 2×2. The δD inversion is completely immune to OH-¹³C (ΔFF = 0).

3. **When isotopes are decoupled** (2×2), each is sensitive to its own OH KIE — confirming the physical intuition that cross-isotope KIE contamination requires coupled solving.

4. **Source signatures and data noise dominate** at 71–80% of variance across all architectures.

5. **Non-OH KIE parameters** (Cl, Strat, Soil) contribute only ~2% combined in all architectures.

6. **KIE importance is model-architecture-dependent**: OH-D dominates in the 3×3 (ill-conditioned δD row is binding constraint), while OH-¹³C dominates in the 2×2 δ¹³C-only (direct pathway, no δD dilution). This means cross-study comparisons of "KIE sensitivity" must account for the inversion architecture.

7. **The δ¹³C–δD FF discrepancy** (178 vs 89 Tg/yr in one-box) highlights the fundamental tension between isotope systems — the 3×3 result (46 Tg/yr) is a compromise weighted toward the more constraining isotope.

---

## 10. Figures

| Figure | Description | File |
|--------|-------------|------|
| Fig. 1 | σ(FF) bar chart by KIE config (3×3) | `figures/fig1_variance_attribution.png` |
| Fig. 2 | OH-¹³C level shift — Saueressig vs Cantrell (3×3) | `figures/fig2_oh13c_level_shift.png` |
| Fig. 3 | Variance decomposition pie + cumulative reduction (3×3) | `figures/fig3_variance_decomposition.png` |
| Fig. 4 | FF time series under 4 configs (3×3) | `figures/fig4_ff_timeseries_comparison.png` |
| Fig. 5 | Dual-role summary diagram (3×3) | `figures/fig5_dual_role_summary.png` |
| Fig. 6 | 2×2 σ(FF) bar chart — δ¹³C and δD panels | `figures/fig6_2x2_variance_attribution.png` |
| Fig. 7 | 2×2 Saueressig vs Cantrell — δ¹³C and δD panels | `figures/fig7_2x2_saueressig_cantrell.png` |
| Fig. 8 | Grand comparison: variance attribution across architectures | `figures/fig8_grand_comparison.png` |
| Fig. 9 | OH-¹³C level shift comparison across architectures | `figures/fig9_level_shift_comparison.png` |

---

## 11. References (Local Repository Files)

| Source | Path |
|--------|------|
| 3×3 experiment code | `experiments/OH_KIE_importance/run_experiment.py` |
| 2×2 experiment code | `experiments/OH_KIE_importance/run_2x2_experiment.py` |
| 3×3 figure generation | `experiments/OH_KIE_importance/make_figures.py` |
| 2×2 figure generation | `experiments/OH_KIE_importance/make_2x2_figures.py` |
| 3×3 results JSON | `experiments/OH_KIE_importance/results/summary.json` |
| 2×2 results JSON | `experiments/OH_KIE_importance/results/2x2_summary.json` |
| Time series CSVs | `experiments/OH_KIE_importance/results/*_timeseries.csv` |
| KIE distributions | `common.py` lines 103–126 |
| 3×3 one-box model | `3x3_one.py` |
| 2×2 one-box model | `2x2_one.py` |
| KIE_immunity experiment | `experiments/KIE_immunity/` |
| dD_threshold experiment | `experiments/dD_threshold/` |
| δ¹³C MC data | `rel/data/d13C_dei_compiled.txt` |
| δD global data | `rel/data/GlobMean_dD_dei_DasguptaCal_noBUDS.csv` |
| Source signature MC | `rel/data/{FF,Mic,BB}_{d13C,dD}_{NH,SH}_MC.csv` |
| Thanwerdas et al. (2024) | `ImportantReferences/Thanwerdas2024ACP/` |
| Fujita et al. (2025) SI | `ImportantReferences/Fujita2025JGR_SI/` |
| Riddell-Young et al. (2025) | `ImportantReferences/Riddell-Young2025PNAS/` |

