# RESULTS.md — Dual-Isotope Methane Source Attribution

**Repository:** `Ilovecodinghhh/upgrade_two_isotope_model`  
**Last updated:** 2026-05-12

---

## Project Overview

This project investigates whether adding δD (hydrogen isotope) to δ¹³C-based methane source attribution improves or degrades the result, and under what conditions. We use a 2-box (NH/SH) mass-balance framework with Monte Carlo uncertainty propagation, benchmarked against Basu et al. (2022 ACP) and other published inversions.

The work is organized into three experiments, each addressing a distinct question:

---

## Experiment 1: δD Threshold (`experiments/dD_threshold/`)

### Question
When does adding δD actually help, and when does it make things worse?

### Core Result
**σ(Mic δD) ≈ 25‰ is the critical threshold.**

| Mic δD uncertainty | FF 90% CI width | vs δ¹³C-only |
|---|---|---|
| 8.3‰ (baseline) | 54.4 Tg/yr | **+46% better** |
| 16.5‰ (2×) | 73.8 Tg/yr | **+27% better** |
| **24.8‰ (3× — threshold)** | **115.3 Tg/yr** | **−14% worse** |
| 41.3‰ (5×) | 206.6 Tg/yr | −104% worse |

**δ¹³C-only reference:** 101.3 Tg/yr CI width.

### Key Findings
1. Below ~25‰ microbial δD uncertainty, δD tightens FF constraints by 27–46%
2. Above ~25‰, δD is *actively counterproductive* — injecting noise faster than signal
3. This explains the Thanwerdas 2024 paradox: their 3D inversion used σ(Mic δD) ≈ 110‰ → δD added "only minor influence"
4. The threshold is robust to KIE choice and lifetime assumptions
5. DFS analysis: δD nearly doubles information content in 2-box (DFS: 2.0 → 3.39)

### Narrative
The literature contradiction on δD utility is NOT about model complexity — it's about a single number: σ(Mic δD). Both sides are correct within their uncertainty assumptions. **Measuring and reducing microbial δD uncertainty below 25‰ is the prerequisite for δD to be useful.**

→ See: `experiments/dD_threshold/RESULT.md` for full phase-by-phase details.

---

## Experiment 2: KIE Sensitivity (`experiments/KIE_sensitivity/`)

### Question
Can dual-isotope mass balance reduce sensitivity to the OH-¹³C KIE controversy (Saueressig 1.0039 vs Cantrell 1.0054)?

### Core Result — Two-Part Answer

**❌ WLS coupling fails:** Treating δD as a hard algebraic constraint in a coupled least-squares system makes KIE sensitivity **5× worse** (KSR = 0.20).

**✅ Agreement filtering succeeds:** Using δ¹³C–δD consistency as a quality filter yields **KSR up to 3.21** and a statistically significant **25.4 pp discriminant** between Cantrell and Saueressig KIE values.

### Agreement Filter Results

| Threshold | Rate (Saueressig) | Rate (Cantrell) | Discriminant |
|---|---|---|---|
| 50 Tg/yr | 14.6% | 33.4% | 18.8 pp (KSR=3.21) |
| **90 Tg/yr (optimal)** | 40.5% | 65.9% | **25.4 pp** |
| 100 Tg/yr | 43.5% | 68.1% | 24.7 pp |

### Robustness
- **Time-varying KIE:** Discriminant survives even under symmetric drift (12.8 pp, still significant) — Phase 7
- **Temporal stability:** Significant across all 3 epochs (1999–2006: 28.3 pp, 2007–2014: 21.5 pp, 2015–2022: 24.1 pp) — Phase 8
- **WLS failure root cause:** δ¹³C and δD respond asymmetrically to OH-¹³C KIE change; WLS resolver amplifies this mismatch 5× instead of damping it

### Publication-Ready Finding
**The δ¹³C–δD agreement rate is a novel observational discriminant for the OH-¹³C KIE.** No previous study has used dual-isotope consistency to discriminate between KIE values. The real atmosphere is more consistent with Cantrell's KIE (1.0054).

→ See: `experiments/KIE_sensitivity/RESULTS.md` for 8-phase detail.

---

## Experiment 3: KIE Immunity (`experiments/KIE_immunity/`)

### Question
How much does the OH-¹³C KIE controversy actually contribute to FF uncertainty, and does real hemispheric δD eliminate it?

### Core Result — Variance Decomposition (Three Configurations)

| Config | σ(FF) Tg/yr | KIE | Source sigs | Lifetime | Residual |
|---|:---:|:---:|:---:|:---:|:---:|
| **δ¹³C-only** | 31.1 | 11.2% | 82.7% | 0.0% | 6.1% |
| **Dual (±6‰ offset hack)** | 17.0 | 0.0% | 0.0% | 14.6% | 85.4% |
| **Dual (real hemi δD)** | **14.3** | **0.0%** | **42.6%** | **27.1%** | **30.3%** |

### Key Findings

1. **KIE contribution = 0% in all dual-isotope configs.** Adding δD completely eliminates KIE as a variance source. The "KIE controversy" is an artifact of under-constrained inversions.

2. **Real hemispheric δD reduces σ(FF) by 54%** (31.1 → 14.3 Tg/yr), vs 45% with the offset hack. The true NH–SH δD gradient (~15‰) is 2.5× larger than the ±6‰ offset assumed previously.

3. **Variance budget becomes interpretable.** With offset δD, 85% was unattributed "residual" (a black box). With real hemispheric data, source signatures (43%) and lifetime (27%) are the dominant identifiable sources — a methodological advance.

4. **FF trend reversal.** δ¹³C-only: ΔFF = +13 Tg/yr post-2007. Dual (real hemi): ΔFF = **−6 to −7 Tg/yr**. Real hemispheric δD *reverses* the inferred FF trend direction, challenging the narrative that fossil fuels drive recent CH₄ growth.

### Basu 2022 Comparison — KIE Spread

| Config | KIE spread (Saueressig vs Cantrell) |
|---|---|
| Basu 2022 (3D, δ¹³C-only) | **13.0 Tg/yr** |
| Our δ¹³C-only | 0.7 Tg/yr |
| Our dual (offset) | 0.9 Tg/yr |
| Our dual (real hemi) | 0.8 Tg/yr |

KIE spread < 1 Tg/yr in ALL our configurations. The entire "KIE controversy" contributing 13 Tg/yr in Basu 2022 is an artifact of their under-constrained system, not a fundamental observational ambiguity.

→ See: `experiments/KIE_immunity/RESULTS.md` for full detail.

---

## Cross-Experiment Synthesis

### The Three Papers' Thesis

These experiments tell a coherent story in 3 parts:

1. **δD helps — but conditionally** (Exp 1): Below σ(Mic δD) ≈ 25‰, δD tightens FF constraints by ~46%. Modern process-based estimates (Douglas 2021, Riddell-Young 2025) achieve σ ≈ 8‰, placing us firmly in the "δD helps" regime.

2. **δD's value is as a diagnostic, not a hard constraint** (Exp 2): Naively coupling δD into WLS makes things worse. The correct use is as an independent consistency filter — which yields a novel KIE discriminant (25.4 pp, significant across all epochs).

3. **Real hemispheric δD eliminates KIE uncertainty entirely** (Exp 3): With station-level hemispheric δD and gridded source signatures, KIE contributes 0% of FF variance. The remaining uncertainty is in source signatures (43%) and lifetime (27%) — both reducible with better measurements.

### Unifying Insight

The literature contradictions on δD's utility, the KIE controversy, and the post-2007 FF trend all resolve when you:
- Use the right uncertainty (σ < 25‰ for microbial δD)
- Use δD as a filter/diagnostic (not WLS coupling)
- Use real hemispheric data (not crude offsets)

---

## Infrastructure: Real Hemispheric δD Data Pipeline

Built 2026-05-12. Key data products:

### Hemispheric atmospheric δD (from Riddell-Young station-level pipeline)
- `rel/data/NHMean_dD_iterations_DasguptaCal_noBUDS.xlsx` — NH MC (19yr × 1000)
- `rel/data/SHMean_dD_iterations_DasguptaCal_noBUDS.xlsx` — SH MC
- True NH–SH gradient: ~15‰ (NH more depleted)
- Coverage: 2005–2019 (real), 2020+ forward-filled

### Hemispheric δD source signatures (1000 MC × 24 years)

| Sector | NH (‰) | SH (‰) | Δ(NH−SH) | Method |
|---|---|---|---|---|
| Mic | −316.9 ± 7.8 | −304.9 ± 7.3 | −11.9 | Douglas 2021 MAT regression × CTCH4 flux |
| BB | −236.7 ± 8.2 | −210.3 ± 7.1 | −26.4 | Umezawa 2011 MAT regression × CTCH4 flux |
| FF | −193.1 ± 5.6 | −189.6 ± 8.1 | −3.5 | Country ONG+coal δD × EDGAR 8.0 |

Files: `rel/data/{Mic,BB,FF}_dD_{NH,SH}_MC.csv`

### 3-box source signatures (also built, not yet used in experiments)

| Sector | NHext (‰) | Trop (‰) | SHext (‰) |
|---|---|---|---|
| Mic | −328.1 | −304.4 | −305.0 |
| BB | −269.6 | −210.5 | −215.6 |
| FF | −194.8 | −187.0 | −202.1 |

Files: `rel/data/{Mic,BB,FF}_dD_{NHext,Trop,SHext}_MC.csv`

---

## Model Code

| Script | Description |
|---|---|
| `common.py` | Shared data loading, MC sampling, KIE, lifetime, hemispheric δD |
| `2x2_one.py` | BB-fixed, separate δ¹³C & δD, global 1-box |
| `2x2_two.py` | BB-fixed, separate δ¹³C & δD, NH/SH 2-box |
| `3x3_one.py` | 3×3 simultaneous (δ¹³C + δD), global 1-box |
| `3x3_two.py` | 3×3 simultaneous, NH/SH 2-box, bounded LS |

See `MODELS.md` for full model descriptions, `sensitivity_test.md` for tunable parameters.
