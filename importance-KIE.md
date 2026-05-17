# importance-KIE: Relative Importance of OH-¹³C vs OH-D KIE in One-Box Dual-Isotope Source Partitioning

## 1. Introduction

The hydroxyl radical (OH) is the dominant CH₄ sink (~84% of total loss), and the kinetic isotope effects (KIE) for the CH₄ + OH reaction — separately for ¹³C/¹²C and D/H — are critical but contested parameters. In the ¹³C system, two laboratory values bracket a range that has dominated the methane budget debate for two decades. In the D/H system, three measurements cluster within ~1.5% relative spread.

**This report quantifies how OH-¹³C and OH-D KIE uncertainties actually propagate into fossil-fuel (FF) emission estimates** using a newly conducted one-box dual-isotope (3×3) experiment with 7000 Monte Carlo inversions.

Key finding: **OH-D KIE dominates the FF uncertainty spread (28% of variance), while OH-¹³C KIE controls the FF level (17 Tg/yr shift) but not the spread (<1% of variance).**

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

## 5. Why OH-D Dominates in the One-Box 3×3

### 5.1 The Condition Number Argument

The 3×3 system solves:

```
[ 1         1         1      ] [BB ]   [S_total          ]
[ f13_BB    f13_FF    f13_Mic] [FF ] = [S_total × f13_src]
[ fD_BB     fD_FF     fD_Mic ] [Mic]   [S_total × fD_src ]
```

The δD row entries (fD) are ~100× smaller in absolute scale than the mass balance row, making the matrix ill-conditioned (mean κ ≈ 10⁴–10⁵ from `3x3_one.py` quality monitor). Because the D row is the *weakest constraint*, any perturbation to the D-system KIE (which changes fD_src) directly destabilizes the solution. OH-D KIE uncertainty thus propagates as *spread*.

In contrast, OH-¹³C KIE perturbs the ¹³C row, which is intermediate in scale between mass balance and D. This row is better conditioned, so KIE perturbations shift the entire solution (level) without inflating the spread.

### 5.2 Fractionation Budget Comparison

| System | Source-to-atm shift | KIE^OH uncertainty (Δε) | Δε / shift |
|--------|:-------------------:|:-----------------------:|:----------:|
| ¹³C | ~6‰ | 1.5‰ | **25%** |
| D | ~230‰ | ~33‰ | **14%** |

Despite OH-¹³C having a larger *relative* KIE uncertainty (25% vs 14% of the fractionation budget), the one-box 3×3 system is more sensitive to OH-D because the δD row is the binding constraint in the ill-conditioned matrix.

### 5.3 Comparison with Two-Box Results

The `KIE_immunity` experiment (two-box, WLS formulation) found KIE accounts for ~25% of FF variance, with OH-¹³C driving the Saueressig/Cantrell sign-change controversy. The difference from our one-box result arises because:

1. **The two-box uses a WLS formulation** (not a direct 3×3 solve), where the W matrix explicitly weights δ¹³C vs δD constraints
2. **Hemispheric resolution** introduces additional δ¹³C sensitivity through the NH–SH gradient
3. **The two-box solves FF and Mic only** (BB is prescribed), changing the problem dimensionality

In the one-box 3×3, the δD row carries more relative weight because all three sources are free, and the δD row provides the decisive third equation.

---

## 6. Source Signature Context

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

## 7. Atmospheric Data Noise

From `rel/data/` MC ensembles:

| Observable | Annual MC std | Trend (2005–2023) | SNR/yr |
|-----------|:------------:|:-----------------:|:------:|
| δ¹³C (global) | 0.006–0.021‰ | −0.035‰/yr | ~2–6 |
| δD (global) | 0.41–0.92‰ | −0.30‰/yr | ~0.3–0.7 |

δD has ~50× larger measurement noise than δ¹³C, consistent with the finding that the D-system is the binding constraint in the one-box inversion.

---

## 8. Summary: The Dual Role of OH KIE

![Dual-role summary](figures/fig5_dual_role_summary.png)
*Figure 5: OH-D KIE controls FF spread (uncertainty); OH-¹³C KIE controls FF level (systematic bias).*

| Property | OH-D KIE | OH-¹³C KIE |
|----------|:--------:|:-----------:|
| **Experimental range** | U(1.294, 1.327) | U(1.0039, 1.0054) |
| **Fractionation (ε)** | ~300‰ | ~5‰ |
| **Effect on σ(FF)** | **28% of variance** | **< 1%** |
| **Effect on mean FF** | < 2 Tg/yr | **17 Tg/yr** |
| **Nature of impact** | Spread (stochastic) | Level shift (systematic) |
| **Reducible by averaging?** | No (drives width) | No (drives bias) |

### Key Takeaways

1. **OH-D KIE is the dominant KIE uncertainty source** in the one-box 3×3 inversion, contributing 28% of FF variance — nearly all of the KIE-attributable uncertainty.

2. **OH-¹³C KIE shifts the FF level by 17 Tg/yr** (Saueressig: 37 vs Cantrell: 54 Tg/yr) but does not change the spread. This makes it a *systematic bias*, not a *random uncertainty*.

3. **Source signatures and data noise dominate** at 71% of variance, consistent with the two-box `KIE_immunity` finding of ~48% for source signatures alone.

4. **Non-OH KIE parameters** (Cl, Strat, Soil) contribute only ~2% combined.

5. **The one-box 3×3 amplifies OH-D importance** relative to the two-box WLS framework because the D row is the weakest constraint in the ill-conditioned 3×3 matrix. This geometry-dependence means KIE importance is model-architecture-dependent — an important caveat for cross-study comparisons.

---

## 9. Figures

| Figure | Description | File |
|--------|-------------|------|
| Fig. 1 | σ(FF) bar chart by KIE config | `figures/fig1_variance_attribution.png` |
| Fig. 2 | OH-¹³C level shift (Saueressig vs Cantrell) | `figures/fig2_oh13c_level_shift.png` |
| Fig. 3 | Variance decomposition pie + cumulative reduction | `figures/fig3_variance_decomposition.png` |
| Fig. 4 | FF time series under 4 configs | `figures/fig4_ff_timeseries_comparison.png` |
| Fig. 5 | Dual-role summary diagram | `figures/fig5_dual_role_summary.png` |

---

## 10. References (Local Repository Files)

| Source | Path |
|--------|------|
| Experiment code | `experiments/OH_KIE_importance/run_experiment.py` |
| Figure generation | `experiments/OH_KIE_importance/make_figures.py` |
| Results JSON | `experiments/OH_KIE_importance/results/summary.json` |
| Time series CSVs | `experiments/OH_KIE_importance/results/*_timeseries.csv` |
| KIE distributions | `common.py` lines 103–126 |
| 3×3 one-box model | `3x3_one.py` |
| KIE_immunity experiment | `experiments/KIE_immunity/` |
| dD_threshold experiment | `experiments/dD_threshold/` |
| δ¹³C MC data | `rel/data/d13C_dei_compiled.txt` |
| δD global data | `rel/data/GlobMean_dD_dei_DasguptaCal_noBUDS.csv` |
| Source signature MC | `rel/data/{FF,Mic,BB}_{d13C,dD}_{NH,SH}_MC.csv` |
| Thanwerdas et al. (2024) | `ImportantReferences/Thanwerdas2024ACP/` |
| Fujita et al. (2025) SI | `ImportantReferences/Fujita2025JGR_SI/` |
| Riddell-Young et al. (2025) | `ImportantReferences/Riddell-Young2025PNAS/` |

