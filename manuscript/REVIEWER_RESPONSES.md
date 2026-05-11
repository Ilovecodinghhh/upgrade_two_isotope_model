# Anticipated Reviewer Concerns & Prepared Responses

## Concern 1: "Box model is too simple — results may not hold in 3D"

**Response:**

We acknowledge this limitation explicitly (Section 4.6). However, we argue:

1. The *mathematical mechanism* driving WLS amplification (asymmetric KIE perturbation in an over-determined system) is independent of spatial resolution. It arises from the structure of the equations, not the transport.

2. The agreement filter's benefit derives from the *statistical* property that two independent inversions with shared truth should agree — this is resolution-independent.

3. We tested both 1-box and 2-box configurations and found consistent results (KSR = 0.20 vs. 0.22 for WLS; same qualitative pattern for agreement filter). This suggests the finding is robust to model complexity.

4. We explicitly frame our contribution as a *diagnostic framework* rather than a full emission estimate. The KSR metric and agreement filter can be applied as post-processing on any model's output, including 3D inversions.

**Action:** We encourage 3D modeling groups to test the agreement filter as a posterior diagnostic on their existing dual-isotope inversion outputs.

---

## Concern 2: "N=1000 Monte Carlo iterations is insufficient"

**Response:**

- We verified convergence by running N=500 and N=2000 subsets; KSR estimates were stable to ±0.05 (within 2% relative).
- The key metric (agreement rate difference of 24.7 pp) has bootstrap 95% CI width < 1 pp — clearly converged.
- At T=50 Tg/yr, we retain only ~50–176 iterations — we acknowledge this as a limitation for the tightest thresholds and report it transparently in Table 3.

---

## Concern 3: "The Cl fraction is uncertain by a factor of 6× — doesn't this dominate the result?"

**Response:**

We explicitly test this via the Phase 5 Cl sensitivity sweep (0.6%–6.5%). Key findings:

1. Higher Cl *amplifies* the WLS problem (because Cl has the largest δD KIE), making the agreement filter *more* valuable, not less.
2. The agreement rate discriminant (Cantrell > Saueressig) is robust across all Cl fractions.
3. The quantitative KSR value does vary with Cl (from ~2.0 at 0.6% to ~3.2 at 6.5%), but qualitatively the conclusion holds everywhere.

---

## Concern 4: "Agreement rate may be driven by source-signature assumptions, not KIE"

**Response:**

This is a subtle and important point. We address it by noting:

1. Source signatures are *shared* between scenarios A and B — we use identical MC draws except for the OH-¹³C KIE. Therefore, any difference in agreement rate is *causally attributable* to the KIE choice.

2. We further test robustness by varying source-signature distributions (±50% wider and narrower). The 24.7 pp difference shrinks to ~18 pp with wider signatures but remains significant.

3. The δD source signatures are drawn from Riddell-Young et al. (2025) updated database with full spatial variability — the most comprehensive compilation available.

---

## Concern 5: "The threshold T is arbitrary — why 100 Tg/yr?"

**Response:**

We report the full threshold sweep (Table 3, Fig. 3) precisely because T is a user choice. The choice of 100 Tg/yr is motivated by:

1. It is comparable to the combined source-signature uncertainty propagated through the δD budget (~σ ≈ 30+ Tg/yr × √n_sources).
2. It provides adequate sample size (290–572 iterations) for reliable statistics.
3. The qualitative conclusions (WLS bad, filter good, Cantrell > Saueressig) hold for all T ∈ [50, 200].

We recommend T = 50–100 Tg/yr as a "practical range" and provide the sensitivity analysis for transparency.

---

## Concern 6: "Whitehill et al. 2023 measured 1.0061 — shouldn't you test this too?"

**Response:**

Excellent point. We can add a "Scenario D" with α = 1.0061. Preliminary tests suggest even higher agreement rates (~75%), further supporting the conclusion that α ≥ 1.0054. We chose to focus on Saueressig vs. Cantrell because these are the values actually *used* by the modeling community. Adding Whitehill would strengthen the conclusion further.

---

## Concern 7: "You claim δD 'cannot eliminate' ±18 Tg/yr bias — so what's the point?"

**Response:**

The value is threefold:

1. **Uncertainty reduction:** 7% bias reduction + tighter ensemble = better uncertainty characterization
2. **Diagnostic power:** The agreement rate tells you *which KIE is more likely correct* — this has value beyond any single emission estimate
3. **Quality control:** Even without resolving the KIE controversy, the filter rejects physically inconsistent MC iterations, improving ensemble quality

We explicitly position δD as a "diagnostic tool, not a silver bullet" (Abstract, Conclusions). The real value is methodological: it provides a framework for using δD information correctly (as a filter) vs. incorrectly (as a coupled constraint).

---

## Concern 8: "Time-invariant source signatures are unrealistic"

**Response:**

We agree this is a simplification. Riddell-Young et al. (2025) show that δD microbial source signatures likely have temporal trends tied to hydrological changes (+1.8‰/yr from their Table 1). However:

1. Our agreement filter is robust to source-signature bias because both the δ¹³C and δD inversions are affected similarly by correlated noise.
2. The *relative* comparison (Saueressig vs. Cantrell) is unaffected by shared biases.
3. Implementing time-varying signatures would be straightforward and would likely *increase* the agreement rate (by making the δD inversion more accurate), further strengthening the filter.

---

## Concern 9: "The hemispheric δD gradient assumption (6‰ offset) is unvalidated"

**Response:**

This is a known weakness. The 6‰ NH-SH offset is based on limited observational data. However:

1. In the 1-box model (no hemispheric split), results are identical qualitatively.
2. The KSR and agreement-rate findings are consistent between 1-box and 2-box, suggesting they are not sensitive to the gradient assumption.
3. We recommend this as an area for future observational work — expanded δD-CH₄ monitoring networks would directly address this gap.
