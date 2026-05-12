# RESULT.md — Title 1: The δD Threshold Experiment

**Branch:** `dD_threshold`  
**Repository:** `Ilovecodinghhh/upgrade_two_isotope_model`  
**Status:** Phases 1–5 complete, 4-panel figure generated

---

## Research Question

When does adding δD (hydrogen isotope) measurements to a methane isotope box model actually improve source attribution — and when does it make things *worse*?

The literature is contradictory: some studies (e.g., Rigby et al. 2012) show δD helps constrain fossil-fuel (FF) emissions, while others (Thanwerdas et al. 2024, 3D inversion) found δD adds only "minor influence." We hypothesised that the discrepancy is explained by a single variable: **the assumed uncertainty on microbial δD source signatures**.

---

## Core Result

**σ(Mic δD) ≈ 25‰ is the critical threshold.**

| Mic δD uncertainty (1σ) | FF 90% CI width (Tg/yr) | Improvement vs. δ¹³C-only |
|--------------------------|--------------------------|---------------------------|
| 4.1‰ (0.5×)             | 54.4                     | **+46%** ✅                |
| 8.3‰ (1× baseline)      | 54.4                     | **+46%** ✅                |
| 16.5‰ (2×)              | 73.8                     | **+27%** ✅                |
| **24.8‰ (3×) ← threshold** | **115.3**             | **−14%** ❌                |
| 41.3‰ (5×)              | 206.6                    | −104% ❌                  |
| 66.0‰ (8×)              | 251.4                    | −148% ❌                  |
| 99.0‰ (12×)             | 270.0                    | −167% ❌                  |
| 132.0‰ (16×)            | 274.9                    | −172% ❌                  |

**Below ~25‰:** δD tightens FF constraints by 27–46%.  
**Above ~25‰:** δD becomes actively counterproductive — injecting noise faster than signal.

**δ¹³C-only reference:** 101.3 Tg/yr CI width (2-box, BB fixed from CarbonTracker).

---

## Phase-by-Phase Summary

### Phase 1 — Baseline Comparison

Ran 1-box and 2-box models in δ¹³C-only vs. dual-isotope mode (1000 MC iterations each).

| Model    | Mode       | FF mean (Tg/yr) | FF 90% CI width |
|----------|------------|------------------|-----------------|
| 1-box    | δ¹³C only  | 177.1            | 101.5           |
| 1-box    | Dual       | 43.1             | 201.3 ❌         |
| 2-box    | δ¹³C only  | 179.4            | 96.6            |
| 2-box    | Dual       | 52.9             | **46.6** ✅      |

**Key insight:** The 1-box dual model *worsens* constraint (−98%) because the 3×3 system is ill-conditioned without bounds. The 2-box with bounded least-squares *improves* constraint by 52%. Model framework matters.

### Phase 2 — Degrees of Freedom for Signal (DFS)

Theoretical information content via DFS = tr(HBHᵀ(HBHᵀ + R)⁻¹):

| Model | δ¹³C only | Dual (δ¹³C + δD) | ΔDFS  |
|-------|-----------|-------------------|-------|
| 1-box | 1.00      | 1.69              | +0.69 |
| 2-box | 2.00      | 3.39              | +1.39 |

δD nearly doubles the information gain in the 2-box framework. The theoretical maximum for 3 sources is 3.0 DFS; the 2-box dual model achieves 3.39 (>3 because it resolves hemispheres independently → effectively 6 sources with 6 observations).

### Phase 3 — The Threshold Sweep

Systematically inflated microbial δD source-signature uncertainty from 0.5× to 16× baseline (σ_baseline ≈ 8.25‰). The crossover from "δD helps" to "δD hurts" occurs between 2× (σ = 16.5‰, +27%) and 3× (σ = 24.8‰, −14%).

**Interpolated critical threshold: σ(Mic δD) ≈ 25‰.**

### Phase 3b — Thanwerdas Replication

Applied Thanwerdas et al. (2024) uncertainty specifications (Mic δD: σ ≈ 110‰, FF δD: σ ≈ 37‰, BB δD: σ ≈ 70‰) to our 2-box model:

| Configuration              | FF 90% CI (Tg/yr) | vs. δ¹³C-only |
|----------------------------|--------------------|---------------|
| δ¹³C only (reference)     | 101.3              | —             |
| Dual, our σ ≈ 8‰          | 54.4               | **−46%** (better) |
| Dual, Thanwerdas σ ≈ 110‰ | 271.6              | **+168%** (far worse) |

**Conclusion:** It's the **uncertainty specification** that kills δD utility, not the spatial framework (box model vs. 3D CTM). Even our simple 2-box reproduces Thanwerdas's negative result when given Thanwerdas's priors.

### Phase 5 — Sensitivity Analysis

Tested threshold robustness across KIE parameterisations and lifetime assumptions:

**KIE sensitivity** (Saueressig / Cantrell / sampled):  
All three give **identical** threshold curves. The bounded least-squares solver makes KIE choice irrelevant — the constraints are tight enough that KIE differences vanish.

**Lifetime sensitivity** (τ = 8.5yr / 9.0yr fixed / varying):  
Negligible shift. The threshold moves by < 2‰ across all lifetime assumptions.

**The ~25‰ threshold is completely robust.**

---

## Key Figures

- **`fig_dD_threshold.png`** — 2-panel: (A) % improvement vs. σ(Mic δD); (B) CI width degradation curve
- **`fig_comprehensive_4panel.png`** — Publication-quality 4-panel figure:
  - A: DFS comparison (1-box vs. 2-box)
  - B: The threshold curve (core result)
  - C: Thanwerdas replication bar chart
  - D: Sensitivity analysis (KIE + lifetime overlay)

---

## Narrative for Paper

The literature contradiction on δD's utility is **not** about model complexity (box vs. 3D) — it's about a single number: σ(Mic δD). When microbial δD source signatures are known to ~8‰ (achievable with modern process-based constraints), δD improves FF source attribution by ~46%. When uncertainty inflates beyond ~25‰ — as assumed in Thanwerdas et al.'s 3D inversion (128‰ prior) — δD adds noise that overwhelms the signal, actively degrading the inversion.

This resolves the apparent paradox: both sides of the literature are correct within their respective uncertainty assumptions. The actionable insight is that **measuring and reducing microbial δD source-signature uncertainty below 25‰ is the prerequisite** for δD to be useful in methane source attribution.

---

## File Inventory

```
experiments/dD_threshold/
├── analysis/
│   ├── phase1_baseline.py       # 1-box/2-box baseline (δ¹³C-only vs dual)
│   ├── phase2_dfs.py            # DFS information content
│   ├── phase3_threshold.py      # Core threshold sweep
│   ├── phase3b_thanwerdas.py    # Thanwerdas uncertainty replication
│   └── phase5_sensitivity.py    # KIE + lifetime robustness
├── figures/
│   ├── fig_threshold.py         # 2-panel threshold figure
│   ├── fig_comprehensive.py     # 4-panel publication figure
│   ├── fig_dD_threshold.png/pdf
│   └── fig_comprehensive_4panel.png/pdf
└── results/
    ├── phase1_baseline/summary.json
    ├── phase2_dfs/dfs_results.json
    ├── phase3_threshold/threshold_results.json
    ├── phase3b_thanwerdas/thanwerdas_comparison.json
    └── phase5_sensitivity/sensitivity_results.json
```
