# Editorial Assessment — JGR: Atmospheres (v1)

**Manuscript:** "Reducing Kinetic Isotope Effect Sensitivity in Methane Source Apportionment via a Dual-Isotope Agreement Filter"

**Role:** Senior Editor evaluation prior to submission

**Date:** 2026-05-12

---

## 1. The 'Big Picture' Critique

### The Introduction Is Too Narrow — and It Buries the Lede

Your Introduction does an excellent job on the *technical* setup (the Cantrell–Saueressig controversy, the KIE propagation math, the survey of 16 studies). But it reads like it's written for the 50 people who already care about OH-¹³C KIE values. A JGR editor's first question will be: **"So what?"**

**What's missing is the policy bridge.** The 20–40 Tg/yr uncertainty you're resolving is not just a methodological nuisance — it's comparable to the *entire* disputed signal that drives climate policy debates about whether to prioritize fossil fuel mitigation vs. agricultural reform. You should say this explicitly. Something like:

> *The difference between the Saueressig and Cantrell partitions — 30–40 Tg/yr shifted between fossil and microbial categories — is comparable to the total methane mitigation pledged under the Global Methane Pledge (2021). Whether post-2007 growth is primarily fossil or microbial directly informs whether policy should target oil/gas infrastructure or agricultural practices.*

Right now, Section 1.2 ("Consequences for Source Attribution") talks about Tg/yr shifts in an abstract way. Section 1.3 (the literature survey table) is useful but belongs in the Discussion — it reads as a literature review in the Introduction, which slows the narrative momentum. **Move the 16-study survey to Section 4.3 ("A Course Correction for the Community") where it actually supports your argument.**

### The "Methane Paradox" Framing Needs Sharpening

You use the phrase in the abstract ("methane paradox") but never define it in the body text. The paradox — rising CH₄ concurrent with declining δ¹³C — should be stated as a one-sentence definition at the top of Section 1.1, because it's the *hook*. Your paper's core claim is that this "paradox" may be partly an artifact of the wrong KIE default, not purely a real shift in source mix. That's a provocative and publishable framing. **Lean into it.**

### Recommendation

Restructure the Introduction as:
1. **§1.1** — The Methane Paradox (2 paragraphs: define it, state the policy stakes)
2. **§1.2** — The KIE Vulnerability (the Cantrell–Saueressig controversy — keep tight)
3. **§1.3** — The Promise of δD (what's been tried, what hasn't)
4. **§1.4** — This Study (your five-approach hierarchy — this is good as-is)

Delete §1.3 (literature survey table) from the Introduction entirely.

---

## 2. Data-to-Claim Mapping: Are You Over-Interpreting the 24.7 pp?

### The Short Answer: Partially Yes, But Fixably

Your headline claim — *"the real atmosphere is more internally consistent with Cantrell's KIE"* — is supported by the data, but the way you frame it conflates **statistical significance with physical proof**, and a careful reviewer will catch this.

### What the Data Actually Shows

Looking at your figures and results:

- **Fig. 13 (fine threshold sweep):** The discriminant is significant at *every* threshold from 30–220 Tg/yr, peaking at 25.4 pp at 90 Tg/yr. The bootstrap CIs are tight and non-overlapping. **This is clean.**
- **Fig. 14 (temporal stability):** 21.5–28.3 pp across three independent epochs. **This is your strongest figure** — it kills the "maybe it's just one atmospheric regime" objection.
- **Fig. 12 (time-varying KIE):** Even under symmetric drift, the discriminant survives at 12.8 pp. **Solid robustness check.**

### The Interpretive Overreach

The problem is in your language, not your numbers. You write (Section 4.2):

> *"This implies that the real atmosphere... is more internally consistent with α = 1.0054 than with 1.0039"*

This is defensible. But then in Section 5.2 you escalate to:

> *"We recommend 1.0054 as the preferred central value"*

That's a bigger leap. Your agreement filter shows that *within your box model framework, with your chosen source signatures, sink fractions, and lifetime parameterization*, Cantrell produces more self-consistent dual-isotope budgets. A reviewer will ask: **"What if the source signatures you chose for δD are systematically biased in a way that happens to favor Cantrell?"**

Your OSSE partially addresses this — but the OSSE uses the *same model* to generate and recover truth, so it can't detect structural model bias. This is a well-known limitation of OSSEs, and you should acknowledge it.

### The Fix

1. **Keep the 24.7 pp finding as-is** — it's novel and significant
2. **Soften the recommendation language** from "we recommend 1.0054 as preferred" to something like: *"Our results provide independent observational evidence consistent with α ≥ 1.0054, adding to the emerging consensus from He et al. (2026b), Dasgupta et al. (2025), and Whitehill et al. (2023). We recommend that the community treat Cantrell's value as at least equally plausible to Saueressig's in default model configurations."*
3. **Add an explicit caveat** in Limitations (§4.6) about the circular-model risk in the OSSE: *"Because the OSSE uses the same box-model structure for truth generation and recovery, it cannot detect systematic biases arising from model structural errors (e.g., neglected transport, time-varying source signatures)."*
4. **In the Discussion, explicitly address the δD source-signature sensitivity.** You partially do this in §4.5, but you should add a sentence like: *"If microbial δD signatures are systematically 20‰ heavier than our assumed range, the agreement-rate discriminant would weaken. However, the temporal stability across three epochs (Fig. 14) argues against a time-invariant signature bias as the explanation."*

---

## 3. Technical Red Flags for a Pedantic Reviewer

### 3a. Sink Fraction Assumptions

Your default sink fractions (OH=83.5%, Cl=3.5%, Strat=7%, Soil=6%) follow Saunois et al. (2020). However:

- **Thanwerdas et al. (2024)** use a lower Cl fraction (~0.6%). You tested this in Phase 5 and show it matters (KIE spread drops from 10 to 6.5 Tg/yr under WLS). **But you don't test the agreement-filter discriminant under different Cl fractions.** A reviewer will ask: "Does the 24.7 pp hold at Cl=0.6%?" This is a gap. Your Phase 5 only tested Cl sensitivity for the WLS approach, not for the agreement filter (Phases 6–8).
- **He et al. (2026b)** allow sink fractions to vary in their 3D CTM. Your fixed fractions are a limitation you should flag more prominently.

**Recommendation:** Run the agreement filter (Phase 6b) with Thanwerdas Cl=0.6% and high Cl=6.5% and report whether the discriminant survives. If you can't re-run, at minimum add a sentence in Limitations acknowledging this untested axis.

### 3b. The 80% Year-Agreement Threshold

In Algorithm 1, you keep iterations where ≥80% of years agree. This is a *second* free parameter (in addition to the Tg/yr threshold) that you never justify or sweep. A reviewer will ask why 80% and not 70% or 90%. You should either:
- Sweep this parameter and show insensitivity, or
- Justify it physically (e.g., "Allowing 20% of years to disagree accounts for years with anomalous biomass burning or δD measurement gaps")

### 3c. N=1000 Monte Carlo Iterations

For the strict 50 Tg/yr threshold, you retain only 53 (Saueressig) and 176 (Cantrell) iterations. The KSR=3.21 at this threshold is based on trend statistics from 53 iterations for the Saueressig case. **That's uncomfortably small.** A reviewer will question whether 53 iterations are sufficient for a robust trend estimate. Consider either:
- Increasing N to 5000 or 10000 (computationally cheap for a box model)
- Reporting confidence intervals on the KSR itself (not just the agreement rates)

### 3d. Lifetime Parameterization

You use τ(t) = 9.0 − 0.017(t − 2010), citing He et al. (2026a). This is a linear decline, giving τ ≈ 9.19 yr in 1999 and τ ≈ 8.80 yr in 2022. Several issues:

- This is a ~4.2% lifetime change over 23 years — significant but aggressive. Some studies (Rigby et al., 2017; Turner et al., 2019) argue for smaller or non-monotonic changes.
- You show (Phase 6b) that lifetime mode has "NO effect" on agreement rates — this is actually a *strength*, not a limitation. **Promote this finding more visibly.** It means your discriminant is insensitive to one of the biggest uncertainties in the field.

### 3e. The δD Data Gap

Your Methods cite Riddell-Young et al. (2025) and Rice et al. (2016) for δD observations. But δD-CH₄ has far sparser coverage than δ¹³C — especially before 2005. How do you handle years with missing or interpolated δD? This needs a sentence in §2.1.

---

## 4. Structural Recommendation: OSSE to Supplement?

### Yes — Move the OSSE Details to Supplementary

Here's the reasoning:

**The OSSE is a validation exercise, not a discovery.** Your paper has two discoveries:
1. WLS coupling amplifies KIE sensitivity (surprising negative result)
2. The agreement-rate discriminant favors Cantrell (novel positive result)

The OSSE confirms the filter gives modest improvement (~7%) and can't eliminate KIE bias. This is useful but expected — and it's not the finding that will get the paper cited. In its current position (§3.5), it deflates the narrative after the high point of §3.4 (the discriminant).

### Recommended Structure

**Main text:**
- §3.3 → Agreement filter results (KSR numbers)
- §3.4 → The discriminant finding (your headline result)
- §3.5 → **NEW: Time-varying robustness (Phase 7) + Temporal stability (Phase 8)** — this is currently *not in the manuscript at all* and it's arguably your strongest evidence
- §3.6 → One paragraph summarizing the OSSE conclusion: "An OSSE confirms the filter reduces bias by ~7% but cannot eliminate fundamental KIE uncertainty (see Supplementary §S3 for details)."

**Supplementary:**
- S1: WLS coupling details (Phase 1–3 numeric tables)
- S2: Weight sweep + Cl interaction details (Phase 5)
- S3: Full OSSE methodology, synthetic truth parameters, recovery tables
- S4: Fine threshold sweep full table (20 thresholds)

**Main text figures (revised):**
- Fig. 1: Schematic (keep)
- Fig. 2: Weight sweep step function (keep — the WLS failure is visually memorable)
- Fig. 3: Threshold sweep with KSR + discriminant (keep)
- Fig. 4: **Replace** the current per-year time series with Fig. 14 (temporal stability by epoch) — this is far more compelling
- Fig. 5: Fig. 12 (time-varying KIE robustness)
- Fig. 6: Summary KSR bar chart (keep)

This restructure gives you a tighter narrative arc: *"We tried coupling (it failed) → we tried filtering (it worked) → the filter reveals a KIE discriminant → the discriminant survives every robustness test we threw at it."*

---

## 5. Summary Scorecard

| Dimension | Current Grade | After Fixes |
|-----------|:---:|:---:|
| Novelty | **A** | A |
| Policy motivation | **C+** | A− |
| Statistical rigor | **B+** | A− |
| Claim calibration | **B−** | A− |
| Robustness testing | **A−** (data exists, not in paper) | A |
| Figure quality | **B+** | A− |
| Structure | **B** | A− |

The core science is strong. The agreement-rate discriminant is genuinely novel — I'm not aware of anyone using dual-isotope consistency as a KIE constraint before. The Phase 7/8 robustness results are your insurance policy against reviewers, and they *must* be in the main text. The main risks are (a) the over-claim on Cantrell, (b) the untested Cl sensitivity for the agreement filter, and (c) the small N at the strict threshold.

---

## Priority Action Items

1. **[HIGH]** Integrate Phase 7 + Phase 8 results into manuscript §3.5 (currently missing entirely)
2. **[HIGH]** Soften KIE recommendation language (§5.2) — "equally plausible" not "preferred"
3. **[HIGH]** Increase MC iterations to N=5000+ or add KSR confidence intervals
4. **[MEDIUM]** Restructure Introduction (move survey table to Discussion, add policy paragraph)
5. **[MEDIUM]** Test agreement-filter discriminant under Thanwerdas low-Cl scenario
6. **[MEDIUM]** Justify or sweep the 80% year-agreement parameter
7. **[LOW]** Move full OSSE to Supplementary; keep one-paragraph summary in main text
8. **[LOW]** Add sentence on δD data handling for pre-2005 years
9. **[LOW]** Add OSSE structural-bias caveat to Limitations
