# Phase 9 — SH source-region sensitivity of the OH ¹³C KIE constraint

**Question (the hypothesis tested):** the KIE_sites headline result rests on two
Southern-Hemisphere sites, CGO and SPO, whose seasonal cycles are treated as
nearly source-free by assigning their wetland source to the **local SH band only**
(`SH_extra`, 2.9 Tg yr⁻¹). But these are background air monitors — the wetland CH₄
reaching them is a *transport-weighted mixture of upstream latitude bands*, not just
local austral wetlands. **Does the SH KIE constraint survive a physically realistic
source-region mixture?**

This phase is a **diagnostic**. It does not modify the Phase 6 production run, the
`SITE_BAND` assignment, or the published headline — it quantifies how fragile that
headline is to the one assumption it most depends on.

## Method

Each emitting band's seasonal wetland phasor, seen at an SH receptor after
interhemispheric transport with exchange time τ, is low-pass filtered:

```
H(τ) = 1 / (1 + i·ω·τ),   ω = 2π yr⁻¹
|H| = 1/√(1+(ωτ)²)   (amplitude attenuation)
arg H = −arctan(ωτ)  (phase lag)
```

At τ = 1.3 yr this gives |H| = 0.122 and a 2.77-month lag — **reproducing the exact
numbers quoted in `phase6_phasor.py`'s own comment block**, so the model is faithful
to the experiment's stated physics. Local (co-located) bands use τ = 0 → H = 1.

The effective source phasor at an SH site sums over emitting bands:

```
Z_src(iso) = Σ_b (δ_src,b − δ_atm) · [(B_Q,b + i·C_Q,b)/Q_total] · H(τ_b)
```

with band δD from Douglas (2021) zonal means (NH_high −374, NH_mid −324, Tropics/SH
−301 ‰) and the existing Phase 6 inversion R → α. δ¹³C is held at the Phase 6 global
value (−62 ‰) for every band so this test isolates the **source-region (δD-gradient +
seasonality)** effect, not a δ¹³C latitude sweep. Monte Carlo (N=50 000, seed 42)
propagates observed amplitude/phase, per-band δD, Q_total, wetland B_Q/C_Q spread,
τ_mix ~ N(1.3, 0.3), and the non-OH KIE/sink priors (reused from `phase6.mc_phasor`).

`analysis/phase9_sh_source_sensitivity.py` · `tests/test_phase9_sh_source.py` (7 tests)
· outputs in `results/phase9_sh_source_sensitivity/` and `figures/fig18,19`.

## Result

**Pooled SH (CGO+SPO) α¹³C_OH by source-region scenario:**

| Scenario | Source region for CGO/SPO | Pooled α | 95% CI | vs lab values |
|----------|---------------------------|---------|--------|---------------|
| **A** local only (Phase 6 baseline) | local SH wetlands, τ=0 | **1.0044** | [0.997, 1.015] | between Sau/Can ✓ |
| **E** local + transported **tropics** | + Tropics, τ=1.3 | 1.0040 | [0.997, 1.019] | between Sau/Can ✓ |
| **D** local + transported **Trop+NH** | + Trop+NH, τ=1.3 | 1.0128 | [0.9995, 1.038] | **above both** |
| **E** local + transported **NH** | + NH bands, τ=1.3 | 1.0114 | [0.999, 1.036] | **above both** |
| **B** tropics only, undamped | Tropics, τ=0 | 1.0184 | [1.001, 1.045] | **above both** |
| **C** global band, undamped | Global, τ=0 | 1.0183 | [1.002, 1.036] | **above both** |

(Saueressig = 1.0039, Cantrell = 1.0054. Scenario A reproduces the published Phase 6
pooled SH value of 1.0044 — a regression guard the test suite enforces.)

### What decides the answer

The result hinges on **one physical fact: which upstream wetlands are seasonal.**

| Band | Annual flux | Fractional seasonal amplitude |
|------|------------:|------------------------------:|
| Tropics | 114 Tg/yr | **0.07** (nearly aseasonal) |
| NH_mid | 30 Tg/yr | 0.98 |
| NH_high | 11 Tg/yr | **1.39** (strongly seasonal) |
| SH_extra | 3 Tg/yr | 0.32 |

- **Transporting tropical wetlands is harmless** (scenario E-tropics: α=1.004). They
  dominate the global *flux* but are almost aseasonal, so they add negligible seasonal
  source phasor regardless of how much reaches the SH.
- **Transporting NH boreal wetlands is what bites** (scenario E-NH / D: α≈1.011–1.013).
  They are strongly seasonal, so even after ~8× transport attenuation (|H|=0.12) their
  residual seasonal signal pulls the SH source phasor up enough to inflate α above both
  lab values.
- The τ_mix sweep (Fig 18b) shows the crossover: at the realistic τ≈1.3 yr the NH
  contribution is damped but **not** to zero, leaving CGO α≈1.009–1.012 and SPO
  α≈1.016 in the mixed scenario.

## Interpretation

**The SH KIE constraint is conditionally robust.** It holds (α≈1.004, favoring
Saueressig) *only if* the seasonal wetland signal reaching CGO/SPO is dominated by
local SH + transported **tropical** sources — both nearly aseasonal. It breaks
(α≈1.011–1.018, above both lab values) if a meaningful fraction of the **seasonal NH
boreal** wetland signal survives interhemispheric transport to the deep SH.

This means the headline result's robustness is **not** primarily a question of the δD
source *value* (the experiment showed, and Phase 9 confirms, that α moves only ~0.001
across −260→−440 ‰ at SH sites). It is a question of source-region **seasonality
geometry** — exactly the dimension Phase 6 collapsed into a single binary `SH_extra`
assignment. The honest SH constraint is therefore **wider than the published
[0.997, 1.015]** once source-region uncertainty is propagated: the physically
plausible scenarios (A, D, E) span pooled α ≈ **1.004 – 1.013**, which no longer sits
cleanly between the lab values and cannot favor Saueressig over Cantrell.

## Recommendation

Phase 6 assigns SH sites to `SH_extra` with a one-line τ_mix justification and the QA
notes record that the opposite extreme (`Global`, undamped) over-corrected by 7×.
Phase 9 shows the truth lives between: the right treatment is the **transport-mixed
scenario D** (local + low-pass-filtered remote bands), not either binary. The
remaining leverage is observational — whether deep-SH δD-CH₄ seasonal phase/amplitude
actually carries a damped NH-boreal signature — which the existing 2005–2010 record is
too short to settle. Until then, the SH constraint should be reported with the wider,
source-region-propagated CI (~1.004–1.013), not the source-region-fixed [0.997, 1.015].

## Files

- `analysis/phase9_sh_source_sensitivity.py` — transport model, 6 scenarios, τ_mix sweep, MC
- `tests/test_phase9_sh_source.py` — 7 tests (transport function, multiband reduction, phase6 reproduction)
- `results/phase9_sh_source_sensitivity/sh_source_sensitivity.json` — all numbers
- `figures/fig18_sh_source_scenarios.png` — α by scenario (a) and α vs τ_mix (b)
- `figures/fig19_sh_source_phasor.png` — source amplitude vs observed signal (a), R by scenario (b)
