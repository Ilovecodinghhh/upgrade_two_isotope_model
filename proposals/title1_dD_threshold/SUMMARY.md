# Title 1: The Spatial Threshold of Deuterium Constraint

**Full Title:** "The Spatial Threshold of Deuterium Constraint: Why δD Helps Box Models but Fails 3D Inversions in Methane Source Attribution"

**Target Journal:** Nature Communications

---

## Scientific Question

At what level of spatial disaggregation does δD-CH₄ lose its independent constraining power on methane source partitioning, and what source-signature uncertainty threshold governs this transition?

## Core Contradiction Addressed

Contradiction C from the literature:
- **δD is essential:** Riddell-Young 2025 (PNAS), He 2026 (Science), Fujita 2025 (JGR)
- **δD is useless:** Thanwerdas 2024 (ACP) — the only 3D dual-isotope inversion found δD added "only a minor influence"

Nobody has explained *why* these conclusions diverge.

## Novelty

Your model suite spans the methodological gap between 1-box models (where δD helps) and 3D inversions (where it doesn't):
- v2.0: Global 1-box (like Riddell-Young/He/Fujita)
- v3.0b/v3.1b: 2-box NH/SH (intermediate — nobody else has this)
- v3.0a/v3.1a: Bayesian framework for formal information quantification

This is the **only** dual-isotope multi-box model in existence.

## Key Datasets and Their Roles

| Dataset | Role |
|---------|------|
| `GlobMean_dD_iterations_UmezawaCal_noBUDS.xlsx` | Global δD MC observations — 1-box baseline |
| `ch4c13_nh_sh_mean.xlsx` | Hemispheric δ¹³C for 2-box; compare information gain |
| `Mic_dD_MC.csv`, `FF_dD_GlobMC_CTCH4.csv`, `BB_dD_annual.csv` | Source signatures — inflate uncertainty to find threshold |
| Riddell-Young 2025 PNAS harmonized δD record | Independent validation (2005–2023) |
| Thanwerdas 2024 Table 2 uncertainty specification | Replicate their setup in 2-box to isolate the failure mode |

## Falsifiable Prediction

δD constraining power collapses above a source-signature uncertainty threshold of approximately 25–30‰ for microbial δD. Below this threshold, even a 2-box model retains δD's discriminating power.

## Impact If Confirmed

Guides investment decisions (~$10M+) in global δD monitoring networks: tells the community exactly what measurement precision is needed before 3D inversions should bother assimilating δD.

---

## Analytical Strategy (Overview)

1. Run 1-box (v2.0) and 2-box (v3.1b) with δ¹³C-only vs. δ¹³C+δD
2. Compute information gain metric (DFS or Bayesian mutual information via PyMC v3.1a)
3. Progressively widen source-signature uncertainties until δD constraint vanishes
4. Identify critical threshold
5. Compare against Thanwerdas 2024's exact uncertainty prescription
6. Validate against Riddell-Young 2025 and He 2026 results
