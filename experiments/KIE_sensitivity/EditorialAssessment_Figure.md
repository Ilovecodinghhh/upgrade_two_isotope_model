# Editorial Assessment — Figure Review

**Manuscript:** "Reducing Kinetic Isotope Effect Sensitivity in Methane Source Apportionment via a Dual-Isotope Agreement Filter"

**Date:** 2026-05-12

---

## 1. Is the Agreement Filter Concept Clear From the Charts Alone?

### Fig. M1 (Schematic) — 7/10

The flowchart communicates the *mechanics* well: two independent solvers → consistency check → filtered ensemble. But it has two weaknesses:

- **Missing rejection path.** The "disagree" path just disappears. A reader scanning the figure won't see *why* filtering helps. Add a red "reject" arrow or a small ghost/faded box showing discarded iterations going to a trash bin or "rejected" pile. The asymmetry between keep/reject is the whole point.
- The σ ≈ 4 vs. σ ≈ 30+ annotation is great — it immediately tells the reader why δ¹³C is the "anchor" and δD is the "filter." Keep that.

### Fig. 8 (Agreement Framework, 4-panel) — 8/10

This is actually your most self-explanatory figure. Panels (a) vs. (b) — unfiltered vs. filtered histograms — make the KSR concept instantly visual: the blue/red overlap tightens after filtering. Panel (d) — per-year agreement rates — shows Cantrell (red) consistently above Saueressig (blue). A reader can *see* the discriminant without reading a word.

**This figure should be in the main text, not supplementary.** It communicates more intuitively than the threshold-sweep plot (Fig. 9).

### Fig. M6 (KSR Summary Bar Chart) — 9/10

This is your best "elevator pitch" figure. The red bars below 1.0 (WLS fails) vs. green bars above 1.0 (filter works) with the dashed KSR=1 line — you can understand the entire paper's conclusion in 3 seconds. The background shading (pink = amplifies, green = reduces) and the "THIS WORK" label are smart design choices.

**Only fix:** The orange T=200 bar is ambiguous — it's KSR=1.09, basically no improvement. Consider making it a more neutral gray to visually separate it from the green "success" bars.

### Fig. 9 (Threshold Sweep, 3-panel) — 6/10

This is information-dense but hard to parse quickly. Panel (a) is fine — two curves separating clearly. But panels (b) and (c) feel like they belong in a technical supplement. A reviewer will understand them; a skimming editor won't. The KSR curve (panel b) is redundant with Fig. M6. The discriminant bars (panel c) lack error bars or significance markers — which feeds directly into the significance question below.

### The Gap

None of your current figures show the **before/after on the actual emission time series**. A reader wants to see: "here's the FF estimate under Saueressig vs. Cantrell *before* filtering, and here it is *after*." Fig. 8 shows this for *trends* (histograms), but not for the year-by-year emission curves. Consider adding a panel with the median FF time series ± CI for unfiltered vs. filtered.

---

## 2. Making the Cantrell vs. Saueressig Significance More Visually Striking

### What's Wrong With Fig. 9(c) (Discriminant Bars)

The purple bars showing the pp difference are *unlabeled for significance*. There are no error bars, no asterisks, no CI brackets. A reader has to trust the text that says p < 0.05. For the figure that's supposed to be your "money shot," this is a missed opportunity.

### What's Wrong With Fig. 14 (Temporal Stability)

This figure is actually *close* to great. Panel (a) has error bars on the paired bars — you can see the Cantrell and Saueressig CIs don't overlap. But the error bars are tiny and hard to see at print scale. Panel (b) has the discriminant values (+28.3, +21.5, +24.1) but again no explicit significance markers.

### Specific Improvements

#### Option A — Redesign Fig. 9(c) as a "Forest Plot" Style Comparison

Instead of bars, show:

```
Cantrell:    ●——[67.5%———68.1%———68.7%]——●
Saueressig:  ●——[42.8%———43.5%———44.1%]——●
                         ← 24.7 pp →
                     *** p < 0.001 ***
```

Horizontal lines with 95% CI, point estimates as dots, and a bracket with the pp difference and significance stars between them. This is the standard way to visually communicate "these two things are different at p < 0.05" in biogeosciences journals. The non-overlap is dramatic at your CI widths (~1 pp each vs. ~25 pp gap).

#### Option B — Add Significance Annotation to Fig. 14

Panel (b) already has the discriminant bars. Add:
- Bracket above each bar pair connecting Saueressig ↔ Cantrell
- Asterisks: `***` (p < 0.001) — your bootstrap CIs are so non-overlapping that this is effectively p ≪ 0.001
- A horizontal dashed line at Δ = 0 (no discriminant)

#### Option C — Create a New Dedicated "Cantrell vs. Saueressig" Figure

A single, clean 2-panel figure:
- **Left panel:** Paired bar chart (like Fig. 14a) but for the *full record* (not split by epoch), with large readable error bars and a bracket showing Δ = 24.7 pp (p < 0.001)
- **Right panel:** The 3-epoch version (Fig. 14) showing the result holds across regimes

This becomes your **Fig. 4** in the restructured manuscript — the "headline figure" that reviewers will remember.

### Color Choice

Right now Saueressig = blue, Cantrell = red. In climate science, red often carries a "bad/warming" connotation. Consider swapping: Cantrell (the one you're arguing *for*) in blue, Saueressig in red/orange. Or use a more neutral palette — teal vs. coral — to avoid implicit bias.

---

## 3. Suggested Caption for the Primary Results Figure (KSR Summary)

> **Figure 3. KIE Sensitivity Ratio (KSR) across all dual-isotope methods tested.** Bars show the ratio of δ¹³C-only KIE sensitivity to dual-isotope KIE sensitivity for fossil fuel emission trends; KSR > 1 (green shading) indicates reduced sensitivity, KSR < 1 (red shading) indicates amplification. Weighted least-squares (WLS) coupling of δD with δ¹³C amplifies KIE sensitivity by 4–5× regardless of configuration (KSR = 0.20–0.24; red bars). In contrast, the agreement filter — which solves the two isotopic budgets independently and retains only Monte Carlo iterations with consistent fossil fuel estimates — progressively reduces KIE sensitivity as the agreement threshold tightens, reaching **KSR = 3.21 at T = 50 Tg yr⁻¹** (rightmost bar). This means the filtered ensemble is 3.2× less sensitive to the choice of OH-¹³C KIE (Saueressig 1.0039 vs. Cantrell 1.0054) than the standard δ¹³C-only inversion. The dashed line marks KSR = 1 (no improvement). The transition bar (T = 200; orange) shows that a lenient threshold provides negligible benefit (KSR = 1.09). See Table 1 for agreement rates and sample sizes at each threshold.

---

## 4. Summary: Recommended Figure Lineup for Main Text

| Position | Figure | Source | Status |
|----------|--------|--------|--------|
| Fig. 1 | Schematic of Agreement Filter | `figM1_schematic.png` | **Revise:** add reject path arrows |
| Fig. 2 | Weight sweep step function (WLS failure) | `fig6_weight_sweep.png` | Keep as-is |
| Fig. 3 | KSR summary bar chart | `figM6_KSR_summary.png` | **Revise:** gray out T=200 bar |
| Fig. 4 | **NEW** — Cantrell vs. Saueressig forest plot + epoch stability | Combine `fig14` + new forest plot | **Create new figure** |
| Fig. 5 | Time-varying KIE robustness | `fig12_timevarying_OH.png` | Keep; add to manuscript text |
| Fig. 6 | Agreement framework 4-panel | `fig8_agreement_framework.png` | **Promote from supplementary** |

### Move to Supplementary

| Figure | Source | Reason |
|--------|--------|--------|
| S1 | `phase1_d13C_only_trends.png` | Baseline reference only |
| S2 | `phase2_dual_isotope_trends.png` | WLS detail |
| S3 | `fig9_threshold_sweep.png` | Redundant with Fig. 3 (KSR summary) |
| S4 | `fig7_Cl_weight_interaction.png` | Phase 5 detail |
| S5 | `fig11_OSSE_recovery.png` | OSSE detail (moved per editorial assessment v1) |
| S6 | `fig13_fine_threshold.png` | Extended threshold table |

### Key Principle

Current lineup: 6 main figures, several underselling your results.
Proposed lineup: 6 main figures, each carrying a single clear message. No figure should require the caption to make its point — the visual alone should tell the story.
