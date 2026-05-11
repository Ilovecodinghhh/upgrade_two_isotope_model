# Experiment: Does Dual-Isotope (δ¹³C + δD) Reduce KIE Sensitivity?

## Hypothesis

> When both δ¹³C and δD are used simultaneously to constrain CH₄ source attribution, the resulting fossil-fuel vs. microbial emission estimates become **significantly less sensitive** to the choice of OH-¹³C KIE (Saueressig 1.0039 vs. Cantrell 1.0054) compared to using δ¹³C alone.

## Background

The OH-¹³C KIE controversy (Contradiction E from our gap analysis) is the single largest systematic uncertainty in isotope-based CH₄ budget studies:

- **Cantrell et al. (1990):** α_OH = 1.0054 ± 0.0009 (ε = −5.4‰)
- **Saueressig et al. (2001):** α_OH = 1.0039 ± 0.0004 (ε = −3.9‰)

This 1.5‰ difference propagates to a **20–40% shift** in the FF/microbial partition (Basu 2022, Lan 2021). Studies using only δ¹³C cannot break this degeneracy.

**The key insight:** δD is governed by a *different* KIE (OH-D α ≈ 1.294–1.327) and provides an independent constraint on source partitioning. If both isotopes are used simultaneously, the system becomes over-determined relative to the KIE uncertainty in either single isotope.

## Experimental Design

### Phase 1: Quantify Baseline KIE Sensitivity (δ¹³C-only)

Run the existing `2x2_one.py` model three times using only the δ¹³C-derived emissions:
1. `OH_13C = 'saueressig'` (fixed at 1.0039)
2. `OH_13C = 'cantrell'` (fixed at 1.0054)
3. `OH_13C = 'sampled'` (uniform across range)

**Metric:** Spread in FF and Microbial emission trends (Δ2020–2022 vs 2005–2007) between runs 1 and 2.

### Phase 2: Quantify Dual-Isotope KIE Sensitivity (δ¹³C + δD combined)

Create a new script `kie_sensitivity.py` that implements a **joint dual-isotope inversion** — using both δ¹³C and δD *simultaneously* as constraints in a single least-squares system (3 equations for 2 unknowns → over-determined):

```
Eq 1: S_total = FF + Mic + BB_fixed           (mass balance)
Eq 2: S × δ¹³C_src = FF×δ¹³C_FF + Mic×δ¹³C_Mic + BB×δ¹³C_BB  (¹³C balance)
Eq 3: S × δD_src = FF×δD_FF + Mic×δD_Mic + BB×δD_BB            (D balance)
```

Solve via weighted least squares (WLS) or bounded `lsq_linear`.

Run this joint system with the same three KIE settings and compute the same metric.

### Phase 3: Comparison & Visualization

**Key outputs:**
1. **KIE Sensitivity Ratio (KSR):**
   ```
   KSR = Spread(δ¹³C-only, Cantrell vs Saueressig) / Spread(dual-isotope, Cantrell vs Saueressig)
   ```
   If KSR > 1 (ideally >> 1), dual isotopes reduce KIE sensitivity.

2. **Histogram comparison:** Overlapping distributions of FF/Mic trends for all 6 runs (3 KIE × 2 methods)

3. **Uncertainty reduction plot:** 2σ uncertainty band width for FF and Mic emissions over time, comparing δ¹³C-only vs dual-isotope

### Phase 4: Extension to 2-Box (NH/SH)

Repeat Phases 1–3 using the hemispheric 2-box model (`3x3_two.py` logic) to test whether spatial disaggregation amplifies or reduces the dual-isotope KIE-dampening effect.

---

## KIE Reference Values (from research_methane_KIE)

### Atmospheric Sink KIEs Used in This Experiment

| Sink | ¹³C KIE (α) | D/H KIE (α) | Source |
|------|-------------|-------------|--------|
| OH (Saueressig) | 1.0039 ± 0.0004 | 1.294 ± 0.018 | Saueressig et al. (2001) |
| OH (Cantrell) | 1.0054 ± 0.0009 | — | Cantrell et al. (1990) |
| OH (D/H alternative) | — | 1.327 | Whitehill-Joelson average |
| Cl | 1.066 ± 0.002 | 1.508 ± 0.04 | Saueressig et al. (1995, 2001) |
| Stratosphere | 1.003 ± 0.001 | 1.179 ± 0.01 | Dyonisius et al. (2020) |
| Soil | 1.020 ± 0.003 | 1.083 ± 0.01 | Snover & Quay (2000) |

### Temperature Dependence (Arrhenius)

| Reaction | Isotope | Formula | Source |
|----------|---------|---------|--------|
| CH₄ + Cl | ¹³C | α(T) = 1.043·exp(6.455/T) | Saueressig (1995) |
| CH₄ + Cl | D/H | α(T) = 1.278·exp(51.31/T) | Saueressig (2001) |
| CH₄ + OH | D/H | α(T) = 1.097·exp(49/T) | Saueressig (2001) |

### Why This Matters

The 1.5‰ discrepancy in OH-¹³C KIE:
- Shifts modelled tropospheric δ¹³C by ~1‰
- Changes FF attribution by 33–65 Tg/yr (Lan et al. 2021)
- Neither value can be confirmed against the real atmosphere (both measured at unnaturally high [OH])
- A new independent measurement is urgently needed (cavity ring-down at realistic conditions)

**Our hypothesis:** The dual-isotope approach makes this uncertainty *less consequential* because δD provides an independent check that is insensitive to OH-¹³C KIE.

---

## Expected Outcomes

| Scenario | Expected Result | Implication |
|----------|----------------|-------------|
| KSR >> 1 (e.g., 3–5×) | δD dramatically dampens KIE sensitivity | Strong case for dual-isotope approaches; publication-ready finding |
| KSR ≈ 1.5–2 | Modest improvement | δD helps but doesn't fully resolve; useful methodological advance |
| KSR ≈ 1 | No improvement | δD and δ¹³C degenerate under KIE perturbation; would challenge Riddell-Young/He claims |
| KSR < 1 | Dual isotope is *worse* | Unlikely but would imply δD source-signature uncertainty dominates |

Based on the literature (Riddell-Young 2025 Fig. 3; He 2026), we expect KSR ≈ 2–4 for the global 1-box, potentially higher for the 2-box where NH/SH gradient provides additional constraint.

---

## File Structure

```
experiments/KIE_sensitivity/
├── README.md                  ← This file
├── PLAN.md                    ← Detailed coding plan with agent prompts
├── kie_sensitivity.py         ← Main experiment script (to be created)
├── results/                   ← Output directory
│   ├── phase1_d13C_only/
│   ├── phase2_dual_isotope/
│   ├── phase3_comparison/
│   └── phase4_two_box/
└── figures/                   ← Publication-quality figures
```

---

## Success Criteria

1. Demonstrate KSR > 1.5 with statistical significance (p < 0.05 via bootstrap)
2. Show that dual-isotope FF trend uncertainty is < 70% of δ¹³C-only uncertainty
3. Produce publication-quality figure suitable for a methods/results section

---

## References

- Cantrell, C. A. et al. (1990). J. Geophys. Res., 95, 22455–22462.
- Saueressig, G. et al. (2001). J. Geophys. Res., 106, 23127–23138.
- Lan, X. et al. (2021). Global Biogeochem. Cycles, 35, e2021GB007000.
- Basu, S. et al. (2022). Atmos. Chem. Phys., 22, 15351–15377.
- Riddell-Young, B. et al. (2025). PNAS, accepted.
- He, J. et al. (2026). Science, in press.
- Thanwerdas, J. et al. (2024). Atmos. Chem. Phys., 24, 2129–2167.
