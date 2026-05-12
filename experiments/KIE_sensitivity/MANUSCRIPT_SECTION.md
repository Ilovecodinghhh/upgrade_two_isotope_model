# Manuscript Section: The δ¹³C–δD Agreement Rate as a KIE Discriminant

**Status:** Draft (2026-05-12). Intended for the Results / Discussion sections
of the Phase 6+ manuscript (see `../../MANUSCRIPT_DRAFT.md`).

---

## Section X — A novel observational test of the OH-¹³C KIE controversy

The 25-year-old disagreement between the Cantrell (1990; α = 1.0054 ± 0.0009)
and Saueressig (2001; α = 1.0039 ± 0.0004) determinations of the OH-CH₄ ¹³C
kinetic isotope effect propagates to a ±20 Tg yr⁻¹ shift in fossil-fuel
versus microbial source attribution in every isotope-based budget study
(Lan et al., 2021; Basu et al., 2022). No subsequent laboratory measurement
has resolved the controversy. Here we show that **the consistency rate
between independent δ¹³C and δD inversions of the global atmospheric record
provides an observation-driven discriminant** between the two candidate
values.

### Method

For each Monte-Carlo iteration *k* and year *j* we solve the δ¹³C-only and
δD-only mass-balance equations independently (see Section X.X) to obtain
two estimates of fossil-fuel emissions, FF_C(*j*,*k*) and FF_D(*j*,*k*).
We define the *agreement indicator*

\[
A_{j,k}(\tau) =
\mathbb{1}\bigl[|FF_C(j,k) - FF_D(j,k)| < \tau\bigr],
\]

and the *agreement rate* as the mean of *A* over valid (years, iterations)
entries. The threshold τ has a clear physical interpretation: it is the
maximum acceptable disagreement between the two independent partitionings.
A small τ implies the two isotopes must give nearly the same FF estimate;
a large τ accepts any solution.

### Result 1 — The agreement rate depends sharply on the assumed OH-¹³C KIE

Figure 13a shows the overall agreement rate as a function of τ for
inversions that assume (i) the Saueressig KIE and (ii) the Cantrell KIE.
At every threshold tested (30 to 220 Tg yr⁻¹, ten-Tg/yr steps), Cantrell's
KIE yields a *higher* agreement rate than Saueressig's, with bootstrap 95 %
confidence intervals that do not overlap.

The agreement-rate difference Δ ≡ rate(Cantrell) − rate(Saueressig)
peaks at **Δ = +25.4 percentage points at τ = 90 Tg yr⁻¹**
(Cantrell: 61.9 %; Saueressig: 36.5 %; bootstrap CI on the difference:
[24.8, 26.0] pp). The KIE
sensitivity ratio KSR(FF) — defined as the ratio of the Cantrell–Saueressig
emission-trend spread under δ¹³C-only to that under agreement-filtering —
peaks at **3.21 at τ = 50 Tg yr⁻¹**.

### Result 2 — The discriminant is robust to time-varying KIE

A reasonable concern is that the bulk OH-¹³C KIE may have drifted over our
1999–2022 analysis window — for example, because He et al. (2026)
report a declining methane lifetime, implying changes in [OH]
which could in principle modulate the effective KIE. Phase 7 of our
analysis tested five trajectories (Table 3): constant Saueressig, constant
Cantrell, symmetric drift toward the 1.0046 midpoint from either end, and
a "convergent" trajectory in which Saueressig drifts to the midpoint while
Cantrell remains constant.

Even in the most aggressive damping case — **symmetric drift toward the
midpoint** — the discriminant retains **12.8 pp** (rate 49.9 % vs 62.7 %)
with statistically non-overlapping bootstrap CIs. In the convergent case
the discriminant remains at 18.7 pp.

Implication: even granting an unphysically large secular drift of either
endpoint toward the midpoint, the dual-isotope agreement rate remains
incompatible with a constant Saueressig KIE.

### Result 3 — The discriminant is stable across atmospheric regimes

Phase 8b splits the 1999–2022 record into three 8-year epochs spanning
qualitatively different atmospheric chemistry regimes (Table 4):

| Epoch | Regime | Δ (Cantrell − Saueressig) |
|-------|--------|---------------------------|
| 1999–2006 | Pre-renewed-growth plateau | **+28.3 pp** (sig.) |
| 2007–2014 | Renewed growth phase | **+21.5 pp** (sig.) |
| 2015–2022 | Post-2014 acceleration | **+24.1 pp** (sig.) |

The discriminant is significant in **every** epoch, despite the very
different growth rates, δ¹³C trends and δD trends that characterise these
periods. The signal is therefore not driven by any single atmospheric
event (e.g. the 2007 inflection, the post-2014 microbial surge, COVID-era
emission anomalies).

### Interpretation

The agreement-rate test rests on an observation-driven premise: that the
true atmospheric history is *internally consistent* — i.e. a single source
partition (FF, microbial, biomass burning) must simultaneously satisfy
the δ¹³C and δD budgets. If an inversion using the wrong OH-¹³C KIE
displaces the inferred δ¹³C source signature, that displacement *cannot*
be compensated by a corresponding shift in δD (which depends on the OH-D
KIE, an unrelated quantity). Iterations that recover a consistent
(FF, microbial) pair from both isotopes are therefore preferentially those
in which the assumed OH-¹³C KIE is closer to truth.

The 25.4 pp gap is large because the δ¹³C and δD source signature
distributions used in our Monte Carlo are independent: ¹³C source
endmembers are drawn from Sherwood et al. (2017) and Schwietzke et al.
(2016); D endmembers from Whiticar (1999) and Thanwerdas et al. (2024).
The probability of accidentally drawing matching FF estimates is small
unless the underlying transformation (i.e. the KIE) is correct.

### Limitations

1. **The test does not give an absolute value for the OH-¹³C KIE.** It
   only shows that Cantrell's value gives substantially higher internal
   consistency than Saueressig's; a value slightly above 1.0054 might
   give even higher agreement.
2. **The agreement rate is sensitive to the assumed δD source endmembers.**
   We use Sherwood/Whiticar/Thanwerdas distributions; sensitivity to
   alternative endmember sets (e.g. He 2026 JGR) is the subject of
   ongoing work.
3. **The 1-box global framework is geographically aggregated.** A 2-box
   hemispheric implementation (Phase 4b) is consistent with the 1-box
   result but has not yet been extended to the agreement-rate test.

### Conclusion

The δ¹³C–δD agreement rate is, to our knowledge, the first observation-
driven discriminant between the Cantrell and Saueressig OH-¹³C KIE
values. The discriminant favours Cantrell's α = 1.0054 over Saueressig's
α = 1.0039 with 25.4 pp peak power, is statistically significant across
a 30–220 Tg yr⁻¹ threshold range, survives plausible time-variation of
the bulk KIE, and is stable across three independent 8-year atmospheric
regimes. While this does not eliminate the need for an improved
laboratory measurement (which remains the gold standard), it provides
the first independent, atmospheric-record-based vote on the controversy.

---

## Suggested figures for this section

| In-manuscript | Source file | Caption headline |
|---------------|-------------|------------------|
| Fig. M-discriminant.a | `figures/fig13_fine_threshold.png` (panel b) | "Agreement-rate discriminant Δ vs threshold, with bootstrap 95% CIs" |
| Fig. M-discriminant.b | `figures/fig14_temporal_stability.png` (panel b) | "Discriminant Δ across three 8-year atmospheric regimes" |
| Fig. S-robustness   | `figures/fig12_timevarying_OH.png` | "Discriminant survival under time-varying OH-¹³C KIE" |

---

## Suggested tables

**Table 3 — Phase 7 scenarios:** copy from
`results/phase7_timevarying_OH/summary.json` (scenarios + discriminants
sections).

**Table 4 — Temporal stability:** copy from
`results/phase8_fine_thresholds/summary.json` (`temporal_stability` block).
