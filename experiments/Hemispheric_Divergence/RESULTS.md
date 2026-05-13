# Hemispheric Divergence — Results

**Date:** 2026-05-13
**Status:** Phases 1–2 complete, Phase 4 running
**Branch:** three-box

---

## Research Question

Do NH and SH methane source trends diverge in a 2-box dual-isotope model, and does this divergence reconcile the disagreement between 3D inversions (FF increasing) and 1-box models (FF stable/declining)?

---

## Key Finding: Source Aliasing, Not Spatial Aliasing

The original hypothesis was that hemispheric aliasing (1-box averaging over NH↑ and SH↓) explains the literature split. **The actual result is more nuanced:**

### What the models show:

| Model | FF trend | Mic trend | BB trend | Total ΔS |
|-------|:---:|:---:|:---:|:---:|
| **2-box hemi** (NH) | **+1.04** (75% pos) | **+3.61** (100% sig ✓) | −0.41 | +4.24 |
| **2-box hemi** (SH) | +0.82 (83% pos) | +0.52 (70% pos) | −0.00 | +1.34 |
| **2-box hemi** (Global) | **+2.10** (89% pos) | **+3.99** (100% sig ✓) | −0.43 | +5.66 |
| **1-box** (Global) | **−0.82** (25% pos) | **+5.42** (100% sig ✓) | +0.99 | +5.59 |

### The Reconciliation

Both models see the **same total source growth** (~5.6 Tg/yr²), but they **partition it completely differently**:

- **2-box:** FF+2.1, Mic+4.0, BB−0.4 → FF is increasing, Mic dominates
- **1-box:** FF−0.8, Mic+5.4, BB+1.0 → FF is declining, Mic+BB absorb all growth

This is **source aliasing**: the 1-box cannot separate 3 sources (FF, Mic, BB) using 2 isotopes in a spatially homogeneous framework. It suffers from a 3×3 degeneracy that collapses FF+BB into a single "non-microbial" category. The 2-box breaks this degeneracy through hemispheric source-signature differences.

### Why This Matters More Than Hemispheric Aliasing

The original hypothesis (NH-FF↑, SH-FF↓, Global-FF≈0) is **not supported**: both NH and SH show slightly positive FF trends. Instead, the key finding is:

> **Adding spatial resolution (2-box vs 1-box) changes the FF/BB/Mic partition, not just the hemispheric split.** The 3×3 overdetermined system in the 2-box is better conditioned because hemispheric signatures provide independent constraints.

This is consistent with:
- **Basu 2022 (3D, FF↑):** Their spatial resolution also breaks the degeneracy
- **Riddell-Young 2025 (1-box, FF stable):** Their global framework inherits the aliasing
- **He 2026 (TROPOMI, FF stable):** Different method, but also spatially resolved

---

## Phase 2.2: Robustness

Pattern tested: NH_FF positive (>50%) AND SH_Mic positive (>60%)

| Config | NH_FF slope | NH_FF %pos | SH_Mic slope | SH_Mic %pos | Pattern? |
|--------|:---:|:---:|:---:|:---:|:---:|
| Default | +0.98 | 75% | +0.52 | 72% | ✓ |
| Fixed τ=9 | +0.95 | 75% | +0.29 | 60% | ✓ |
| Cantrell (1.0054) | +0.79 | 72% | +0.83 | 85% | ✓ |
| Saueressig (1.0039) | +1.26 | 81% | +0.18 | 56% | ✗ |
| Low Cl (0.6%) | +1.49 | 84% | −0.06 | 48% | ✗ |
| High Cl (6.5%) | +0.61 | 69% | +1.31 | 100% | ✓ |
| Fast τ_ex (0.8) | +1.46 | 82% | +0.56 | 80% | ✓ |
| Slow τ_ex (1.3) | +0.50 | 68% | +0.52 | 68% | ✓ |

**6/8 configurations** support the pattern. The two failures (Saueressig, Low Cl) suppress SH_Mic rather than reversing NH_FF.

Key observations:
- **NH_FF is positive in all 8 configurations** (68-84% probability)
- **SH_Mic varies most with Cl fraction:** High Cl amplifies the δD constraint, making SH_Mic strongly positive
- **KIE matters:** Cantrell shifts growth toward Mic, Saueressig toward FF (consistent with KIE_immunity findings)

---

## Phase 4: Exchange Rate Sensitivity

Observed IH δ¹³C gradient (NH−SH): **−0.237 ± 0.026‰**

| τ_ex (yr) | NH_FF slope | NH_FF %pos | SH_Mic slope | SH_Mic %pos | Global FF slope |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 0.7 | +1.74 | 85% | +0.51 | 86% | +2.16 |
| 0.8 | +1.46 | 82% | +0.56 | 80% | +2.08 |
| 0.9 | +1.19 | 79% | +0.53 | 72% | +2.10 |
| **1.0** | **+1.01** | **76%** | **+0.54** | **70%** | **+2.03** |
| 1.1 | +0.82 | 74% | +0.50 | 70% | +1.92 |
| 1.2 | +0.66 | 71% | +0.48 | 68% | +1.81 |
| 1.5 | +0.30 | 61% | +0.62 | 70% | +1.57 |
| 2.0 | +0.15 | 56% | +0.81 | 81% | +1.23 |

**NH_FF is positive at all τ_ex values.** The signal strength scales inversely with exchange time: faster exchange (more hemispheric independence) amplifies the NH_FF trend. Even at τ_ex = 2.0 yr (unphysically slow), NH_FF is still 56% positive.

SH_Mic is stable and positive across all values (68-86% positive), with a slight increase at extreme slow exchange.

Global FF trend is always positive (1.23–2.16 Tg/yr²), declining monotonically with τ_ex.

---

## Implications

1. **The 1-box vs 3D disagreement is a degeneracy problem, not a spatial aliasing problem.** Adding even one spatial dimension (2 boxes) breaks the FF/BB degeneracy that plagues 1-box dual-isotope inversions.

2. **Both Basu (2022) and Riddell-Young (2025) could be correct — but for the wrong reason.** Basu's FF increase may be real (consistent with our 2-box), but Riddell-Young's "stable FF" may be an artifact of 1-box degeneracy rather than a genuine signal.

3. **The key discriminator is BB.** 2-box finds BB declining (−0.4 Tg/yr²); 1-box finds BB increasing (+1.0). This is where the degeneracy manifests. GFED fire data can potentially arbitrate.

---

## File Inventory

```
experiments/Hemispheric_Divergence/
├── SUMMARY.md
├── PLAN.md
├── RESULTS.md
├── analysis/
│   ├── run_models.py          (Phase 1: 2-box hemi, 2-box global, 1-box)
│   ├── hemispheric_trends.py  (Phase 2: trend analysis + aliasing test)
│   ├── divergence_robustness.py (Phase 2.2: 8-config robustness)
│   └── exchange_rate_sensitivity.py (Phase 4: τ_ex sweep)
└── results/
    ├── twobox_hemi/
    ├── twobox_global_sigs/
    ├── onebox_reference/
    ├── trend_analysis.csv
    ├── aliasing_test.json
    └── robustness_table.csv
```
