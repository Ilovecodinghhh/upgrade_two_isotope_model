# CHANGELOG — v3.1 and v3.2 (2026-05-05)

## Overview

Two new model variants, both incorporating insights from Riddell-Young et al. (2025, PNAS):

| Version | Approach | File | Output Dir |
|---------|----------|------|------------|
| **v3.1** | Optimized 3×3 (our method + Ben improvements) | `v3.1_optimized_3x3.py` | `Output_v3.1_3x3/` |
| **v3.2** | 2×2 BB-fixed (Ben's method + our specialties) | `v3.2_bb_fixed_2x2.py` | `Output_v3.2_2x2/` |
| v3.0 (prev) | Original 3×3 two-hemisphere | `two_hemisphere_box_model.py` | `Output_2Hemi/` |

---

## v3.1 — Optimized 3×3 (Our Method + Ben Improvements)

### Changes from v3.0

1. **δD hemispheric offset corrected: ±6‰** (was ±1.5‰)
   - Riddell-Young (2025): "NH δD–CH₄ is ~12‰ lower than SH δD–CH₄"
   - Annual mean: ~12‰ total offset → ±6‰ from global mean
   - This is the single biggest change and dramatically shifts results

2. **τ_ex sampled as Normal(1.0, 0.1)** per MC iteration
   - Previously fixed at 1.0 yr
   - Source strengths now vary per iteration (both mass balance and isotopes)

3. **All 4 sinks now sampled** (was only OH + Cl)
   - Strat: N(1.003, 0.001) for ¹³C; N(1.179, 0.01) for D
   - Soil: N(1.0201, 0.003) for ¹³C; N(1.083, 0.01) for D
   - References: Dyonisius 2020; Beck 2018; Snover & Quay

4. **5-year moving average smoothing** on all MC iterations before stats
   - Reduces interannual noise, extracts robust trends
   - Edge-padded (3-point at ends, 4-point near-ends)

5. **Trend analysis** (Ben's method): Δ(2020–2022 avg vs 2005–2007 avg)
   - Reports fraction of MC runs with positive/negative trends

6. **Microbial δD uncertainty updated**: σ = 8.2‰ (from Ben's Table 1; was 7.0‰)

### v3.1 Results

```
GLOBAL (Smoothed):
  BB:  88.2 ± 16.2 Tg/yr
  FF:  52.9 ± 14.4 Tg/yr
  Mic: 444.2 ± 20.5 Tg/yr

TRENDS Δ(2020–2022 vs 2005–2007):
  FF:  +4.4 ± 9.3 Tg/yr (positive in 61% of MC runs) → STABLE
  Mic: +65.2 ± 5.2 Tg/yr (positive in 100% of MC runs) → STRONG INCREASE
  BB:  +4.5 ± 5.2 Tg/yr (positive in 83% of MC runs) → SLIGHT INCREASE

Quality: 0% non-physical (bounded LS), mean cond = 178,412
```

### Key Observation

The corrected δD offset (±6‰) causes the 3×3 system to allocate almost all FF to SH and almost all BB to NH. This is physically questionable — BB should be 55/45 NH/SH, and FF should be 85/15 NH/SH. The large δD offset combined with the ill-conditioned 3×3 system causes the bounded least squares to find unusual hemisphere-specific solutions (even though global totals are reasonable).

**Implication**: The 3×3 approach struggles with the correct δD hemispheric gradient. The ±1.5‰ in v3.0 was "accidentally" making it work. This further motivates the 2×2 approach in v3.2.

---

## v3.2 — 2×2 BB-Fixed (Ben's Method + Our Specialties)

### Approach

Following Riddell-Young et al. (2025, PNAS):
1. Fix BB from CarbonTracker/GFED4 prior (~29 Tg/yr)
2. Solve **separately** for FF and Mic using δ¹³C alone, then δD alone
3. Compare the two independent estimates (cross-validation)

### Our Specialties Retained (beyond Ben's one-box)

1. Two-hemisphere structure (NH/SH + interhemispheric exchange)
2. KIE sampling in MC loop (all 4 sinks)
3. Time-varying τ: 9.0 - 0.017*(t - 2010)
4. τ_ex sampled: Normal(1.0, 0.1)
5. Hemisphere-specific sink fractions and lifetimes
6. δD hemispheric offset = ±6‰
7. NH/SH BB split from GFED4 (55/45%)

### v3.2 Results

```
δ¹³C-derived (Global, Smoothed):
  FF:  180.4 ± 29.4 Tg/yr
  Mic: 376.3 ± 29.4 Tg/yr
  BB:  29.0 Tg/yr (fixed)

δD-derived (Global, Smoothed):
  FF:  93.4 ± 29.4 Tg/yr
  Mic: 463.2 ± 29.4 Tg/yr
  BB:  29.0 Tg/yr (fixed)

δ¹³C vs δD agreement:
  FF difference = 87.0 Tg/yr (Ben's paper: ~27 Tg/yr)
  → Our δD offset may be too large, or δD data processing needs refinement

TRENDS Δ(2020–2022 vs 2005–2007):
  FF (δ¹³C): +10.6 ± 4.0 (positive in 100%) → SLIGHT INCREASE
  Mic (δ¹³C): +63.5 ± 4.0 (positive in 100%) → STRONG INCREASE
  FF (δD): +13.4 ± 12.2 (positive in 84%) → AMBIGUOUS
  Mic (δD): +60.7 ± 12.2 (positive in 100%) → STRONG INCREASE

Negative solutions:
  δ¹³C: 6.3% (low — good!)
  δD: 33.7% (high — δD still problematic)
```

### Key Observations

1. **δ¹³C results (FF=180, Mic=376) agree well with Ben's paper** (FF=160±29, Mic=~400)
2. **δD results diverge more** (FF=93 vs Ben's FF=133±33)
3. **The δ¹³C/δD divergence (87 Tg/yr) is larger than Ben's (~27 Tg/yr)**
   - Likely because our hemispheric δD offset (±6‰) is too crude
   - Ben used actual station-level δD data with proper spatial averaging
   - The ±6‰ may over-correct because Ben's 12‰ includes seasonal extremes
4. **Trend results are robust**: Both isotopes agree Mic increased ~60-65 Tg/yr
5. **33.7% negative solutions from δD** confirms δD is inherently less constraining

---

## Comparison: v3.0 vs v3.1 vs v3.2

| Metric | v3.0 | v3.1 | v3.2 (δ¹³C) | v3.2 (δD) | Ben (2025) |
|--------|------|------|-------------|-----------|------------|
| FF Global (Tg/yr) | 189.9 | 52.9 | **180.4** | 93.4 | 160±29 |
| Mic Global (Tg/yr) | 371.1 | 444.2 | **376.3** | 463.2 | ~400 |
| BB Global (Tg/yr) | 24.5 | 88.2 | **29.0** (fixed) | 29.0 (fixed) | ~29 |
| δD hemispheric offset | ±1.5‰ | ±6.0‰ | ±6.0‰ | ±6.0‰ | ~12‰ total |
| Non-physical % | 0% | 0% | 6.3% | 33.7% | N/A |
| Mic trend (Tg/yr) | N/A | +65.2 | +63.5 | +60.7 | +73±5 |
| FF trend (Tg/yr) | N/A | +4.4 | +10.6 | +13.4 | stable |

### Conclusions

1. **v3.2 δ¹³C is our best result** — matches Ben's published values closely
2. **v3.1 is unreliable** — the ±6‰ δD offset breaks the 3×3 system's ability to partition hemisphericaly
3. **δD is always weaker** — even with the 2×2 approach, 33.7% negative solutions
4. **The key scientific finding is robust across all approaches**: Microbial emissions increased ~60-73 Tg/yr since 2005-2007; FF emissions are stable or slightly increasing

---

## Recommended Path Forward

### For publication-quality results:
1. Use **v3.2 (δ¹³C-derived)** as the primary result
2. Report v3.2 (δD-derived) as independent confirmation of trends
3. Show δ¹³C/δD divergence explicitly (as Ben does) to discuss remaining uncertainties

### To reduce δ¹³C/δD divergence:
1. Use Ben's station-level δD data to derive actual NH/SH means (not ±6‰ approximation)
2. The ±6‰ may be too large for annual means (seasonal max is 15‰, min is 8‰; annual mean difference may be ~10‰ → ±5‰)
3. Try ±4‰ or ±5‰ and compare against Ben's results

### To improve hemispheric partitioning in v3.1:
- The 3×3 with correct δD offset doesn't work well hemisphericaly
- Would need hemisphere-specific source signatures (not just global)
- Or switch to informative Bayesian priors on hemispheric allocations

---

## Files Created

```
upgrade_two_isotope_model/
├── v3.1_optimized_3x3.py              ← Optimized 3×3 model
├── v3.2_bb_fixed_2x2.py               ← BB-fixed 2×2 model
├── Output_v3.1_3x3/
│   ├── v3.1_results_smoothed.csv
│   ├── quality_report.json
│   ├── v3.1_hemispheric_sources_smoothed.png
│   └── v3.1_trend_histograms.png
├── Output_v3.2_2x2/
│   ├── v3.2_results_smoothed.csv
│   ├── quality_report.json
│   ├── v3.2_d13C_vs_dD_comparison.png
│   └── v3.2_trend_histograms.png
├── BEN_MODEL_ANALYSIS.md              ← Full analysis of Ben's model
└── CHANGELOG_v3.1_v3.2.md            ← This file
```
