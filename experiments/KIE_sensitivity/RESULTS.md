# KIE Sensitivity Experiment — Complete Results

**Date:** 2026-05-11  
**Repository:** Ilovecodinghhh/upgrade_two_isotope_model

---

## Executive Summary

This experiment tested whether combining δ¹³C and δD reduces sensitivity to the OH-¹³C KIE controversy (Saueressig 1.0039 vs Cantrell 1.0054). The answer is **definitively no** — adding δD as a WLS constraint *increases* KIE sensitivity in all tested configurations.

### Key Numbers

| Phase | Model | KSR (FF) | KSR (Mic) | Verdict |
|-------|-------|---------|----------|---------|
| 1+3 | 1-box, δ¹³C-only baseline | — | — | Spread = 2.0 Tg/yr |
| 2+3 | 1-box, dual-isotope WLS | 0.20 | 0.32 | **5× worse** |
| 4b | 2-box (fixed exchange), dual | 0.22 | 0.35 | **Same pattern** |
| 5 | Weight sweep (w_dD=0 to 1) | monotonically worsens | | **No optimal weight** |

**KSR < 1 everywhere** means dual-isotope always increases KIE sensitivity. There is no sweet spot.

---

## Phase 1: δ¹³C-Only Baseline (1-box)

| Run | OH-¹³C KIE | FF Trend ± σ (Tg/yr) | Mic Trend ± σ (Tg/yr) |
|-----|------------|----------------------|------------------------|
| A (Saueressig) | 1.0039 | +11.4 ± 4.1 | +62.5 ± 4.1 |
| B (Cantrell) | 1.0054 | +9.4 ± 4.2 | +64.4 ± 4.2 |
| C (Sampled) | U[1.0039–1.0054] | +10.4 ± 4.2 | +63.5 ± 4.2 |

**KIE spread:** |B − A| = **2.0 Tg/yr** for both FF and Mic  
**Key feature:** Per-run uncertainty is small (σ ≈ 4 Tg/yr) — the system is well-conditioned.

---

## Phase 2: Dual-Isotope WLS (1-box)

| Run | OH-¹³C KIE | FF Trend ± σ (Tg/yr) | Mic Trend ± σ (Tg/yr) |
|-----|------------|----------------------|------------------------|
| A (Saueressig) | 1.0039 | +11.1 ± 32.4 | +60.5 ± 22.6 |
| B (Cantrell) | 1.0054 | +21.2 ± 39.2 | +54.3 ± 27.1 |
| C (Sampled) | U[1.0039–1.0054] | +15.9 ± 36.1 | +57.6 ± 25.1 |

**KIE spread:** |B − A| = **10.1 Tg/yr** (FF), **6.2 Tg/yr** (Mic)  
**Total uncertainty:** σ ≈ 32–39 Tg/yr — **8× worse than δ¹³C-only!**

---

## Phase 4b: Two-Box (Fixed Exchange Isotopes)

With corrected isotopic exchange treatment, emissions are realistic:

| Run | Method | FF Mean (Tg/yr) | Mic Mean (Tg/yr) | FF Trend | Mic Trend |
|-----|--------|-----------------|-------------------|----------|-----------|
| A (Saueressig) | Dual | 145 | 414 | +4.2 ± 23.1 | +64.9 ± 16.4 |
| B (Cantrell) | Dual | 162 | 408 | +13.0 ± 23.3 | +59.6 ± 16.5 |
| A (Saueressig) | δ¹³C-only | 202 | 354 | +16.2 ± 3.8 | +57.9 ± 3.8 |
| B (Cantrell) | δ¹³C-only | 165 | 392 | +14.3 ± 3.7 | +59.8 ± 3.7 |

**KSR (2-box):** FF = 0.22, Mic = 0.35 — **same as 1-box**

**Sanity check:** Dual-isotope gives FF ≈ 145–162 Tg/yr, Mic ≈ 408–414 Tg/yr → Total ≈ 560 Tg/yr. This is consistent with He 2026 Science (575 Tg/yr total, ~52 Tg FF from oil/gas alone + coal + other → ~130–150 total FF).

---

## Phase 5: δD Weight Sweep + Cl Fraction Sensitivity

### Test 1: δD Weight vs KIE Spread (default Cl = 3.5%)

| δD Weight | KIE Spread FF (Tg/yr) | Total σ FF (Tg/yr) | KIE Spread Mic | Total σ Mic |
|-----------|----------------------|-------------------|----------------|-------------|
| **0.00** (pure δ¹³C) | **1.98** | **4.3** | **1.98** | **4.3** |
| 0.01 | 8.38 | 24.9 | 4.99 | 17.2 |
| 0.05 | 9.39 | 33.7 | 5.72 | 23.4 |
| 0.10 | 9.72 | 39.1 | 5.96 | 27.2 |
| 0.20 | 9.89 | 40.0 | 6.08 | 27.9 |
| 0.50 | 9.97 | 41.2 | 6.14 | 28.7 |
| 1.00 | 9.99 | 42.5 | 6.17 | 29.5 |

**Critical insight:** Even a tiny δD weight (w=0.01) causes KIE spread to jump from 2 → 8 Tg/yr (4× worse) and total uncertainty from 4 → 25 Tg/yr (6× worse). The function is essentially a step function — **there is no gradual trade-off**.

### Test 2: Cl Fraction × δD Weight Interaction

| Cl Config | Cl Fraction | w_dD=0 | w_dD=0.1 | w_dD=0.5 |
|-----------|------------|--------|----------|----------|
| Thanwerdas (low) | 0.6% | 2.14 | 6.49 | 6.49 |
| Default | 3.5% | 1.98 | 9.72 | 9.97 |
| High Cl | 6.5% | 1.84 | 14.45 | 15.10 |

**Critical interaction:** Higher Cl fraction *amplifies* the δD-induced KIE degradation (from 6.5 to 15 Tg/yr spread). This is because Cl has the largest δD KIE (α=1.52) — when OH-¹³C KIE changes, the δD budget propagates this through the Cl term more strongly.

---

## Root Cause Analysis

### Why does δD make things worse?

The fundamental mechanism:

1. **OH-¹³C KIE directly affects the δ¹³C source δ calculation** — shifting the inferred δ¹³C_src by ~1‰
2. **This shifts the δ¹³C equation** in the WLS system
3. **The δD equation is unaffected by OH-¹³C** (it depends on OH-D KIE, which we don't perturb)
4. **But the WLS solver must simultaneously satisfy both equations** — so the shifted δ¹³C row creates a contradiction with the unshifted δD row
5. **The solver resolves this contradiction by moving FF/Mic further** than it would with δ¹³C alone

In other words: δ¹³C-only gives a clean analytic solution where the KIE shift maps linearly to a small FF/Mic shift. Adding δD creates an over-determined system where the same KIE perturbation produces a larger least-squares residual that gets distributed across both unknowns.

### Mathematical explanation:

For δ¹³C-only:
```
FF = f(δ¹³C_src) = f(α_OH)  →  ∂FF/∂α_OH ≈ 2 Tg/yr per 0.0015 shift
```

For dual-isotope WLS:
```
[FF, Mic] = argmin ||W(Ax - b)||²
When row 2 (δ¹³C) shifts but row 3 (δD) doesn't:
→ larger residual → solver redistributes over both unknowns
→ ∂FF/∂α_OH ≈ 10 Tg/yr per 0.0015 shift (5× amplification)
```

---

## Implications

### 1. Simple WLS is the wrong approach for combining isotopes

The naive "more data = better constraint" intuition fails here. Adding δD to a WLS system:
- ✗ Does NOT reduce KIE sensitivity
- ✗ Does NOT reduce total uncertainty
- ✗ Gets WORSE with higher Cl fraction
- ✗ Has no optimal weighting (even w=0.01 breaks it)

### 2. What WOULD work instead?

Based on the literature (Thanwerdas 2024, He 2026 JGR):
- **Bayesian/MCMC:** Rather than WLS, use δD as a *prior constraint* that penalizes unphysical solutions without fully coupling the systems
- **Sequential filtering:** Solve δ¹³C first, then use δD to reject outliers rather than to influence the mean
- **Agreement metric:** Use the residual norm of the δD equation as a *diagnostic* (quality indicator) rather than a *constraint*
- **Separate inversions with consistency check:** Solve δ¹³C → FF/Mic, solve δD → FF/Mic, report overlap zone

### 3. The literature conflict resolved

| Study | Approach | Conclusion | Why |
|-------|----------|-----------|-----|
| Riddell-Young 2025 | Separate 2×2 (not WLS) | δD helps | They solve each isotope independently and check consistency |
| He 2026 JGR | 3D CTM + cost function | δ¹³C constrains, δD validates | They don't combine into WLS |
| Thanwerdas 2024 | Bayesian framework | δD adds "minor influence" | They use it as weak prior, not hard constraint |
| **This work** | WLS (hard coupling) | δD makes things worse | Over-determined system amplifies perturbations |

**Resolution:** δD's value is as an *independent validation* or *soft prior*, NOT as a hard algebraic constraint in a coupled system.

---

## Next Experiment Suggested

**Phase 6: Bayesian Agreement Framework**

Instead of WLS, implement:
1. Solve δ¹³C-only → get FF/Mic posterior (narrow, σ≈4 Tg/yr)
2. Independently solve δD-only → get FF/Mic posterior (wide, σ≈30+ Tg/yr)
3. Compute overlap zone = intersection of both posteriors
4. Test whether the overlap zone is *less* sensitive to KIE than δ¹³C alone

This would test the Riddell-Young (2025) "consistency" approach directly.

---

## Files Produced

```
results/
├── phase1_d13C_only/          (run_A/B/C.npz + summary.json)
├── phase2_dual_isotope/        (run_A/B/C.npz + summary.json)
├── phase3_comparison/          (summary.json)
├── phase4_two_box/             (original buggy — for reference only)
├── phase4b_two_box_fixed/      (corrected exchange — run_*.npz + summary.json)
└── phase5_weight_Cl_sweep/     (summary.json)

figures/
├── phase1_d13C_only_trends.png
├── phase2_dual_isotope_trends.png
├── fig1_KSR_summary.png
├── fig2_uncertainty_timeseries.png
├── fig3_emission_timeseries.png
├── fig4_KSR_1box_vs_2box.png      (original — misleading due to exchange bug)
├── fig5_2box_fixed.png
├── fig6_weight_sweep.png
└── fig7_Cl_weight_interaction.png
```
