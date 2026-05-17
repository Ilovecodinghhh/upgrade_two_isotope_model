# OSSE (Observing System Simulation Experiment) — Summary

## Purpose

Test whether the agreement filter can improve emission recovery when the true answer is known. A controlled end-to-end validation of the method.

## How It Works

### 1. Define Ground Truth

- **True KIE**: OH-¹³C = 1.0046 (midpoint of Saueressig/Cantrell)
- **True source signatures**: fixed constants (FF δ¹³C = −44‰, Mic δ¹³C = −62‰, BB δ¹³C = −22‰, plus corresponding δD values)
- **True emissions**: real observed CH₄ → compute total source from mass balance → partition as 24% FF, ~5% BB, ~71% Mic (fractions from He 2026 Science)

### 2. Forward-Model a Synthetic Atmosphere

Starting from real first-year δ¹³C/δD observations, integrate forward year by year:

```
n13[j+1] = n13[j] × (1 − α/τ[j]) + S[j] × f13_src[j]
```

This produces 23 years of synthetic δ¹³C and δD that are perfectly self-consistent with the true emissions and true KIE. Real CH₄ concentrations are used (only isotope values are synthetic).

### 3. Invert with Wrong KIEs

- Add realistic observational noise (σ = 0.05‰ for δ¹³C, σ = 3.0‰ for δD)
- Run 1000 MC iterations using **real source-signature distributions** (not the fixed true values)
- Test three KIE assumptions: true (1.0046), Saueressig (1.0039), Cantrell (1.0054)
- Compare: unfiltered δ¹³C-only results vs agreement-filtered results

## Key Design Choices

| Element | Forward model (truth) | Inversion (recovery) |
|---------|----------------------|---------------------|
| CH₄ concentrations | Real observed | Real observed |
| Source signatures | Fixed constants | MC-sampled, time-varying (real data) |
| KIE | 1.0046 (fixed) | Tested: 1.0046 / 1.0039 / 1.0054 |
| Atmospheric δ¹³C, δD | Synthetic (forward-modeled) | Synthetic + noise |

The deliberate mismatch in source signatures makes the test realistic — in the real world we don't know exact signatures either.

## Results

| Inversion KIE | Bias (unfiltered) | RMSE (unfiltered) | Bias (filtered) | RMSE (filtered) |
|---------------|:-:|:-:|:-:|:-:|
| True (1.0046) | +2.1 Tg/yr | 11.6 Tg/yr | +1.5 Tg/yr | 11.4 Tg/yr |
| Saueressig (1.0039) | +19.6 Tg/yr | 22.6 Tg/yr | +18.3 Tg/yr | 21.4 Tg/yr |
| Cantrell (1.0054) | −17.8 Tg/yr | 21.2 Tg/yr | −17.4 Tg/yr | 20.8 Tg/yr |

## Conclusions

1. **The agreement filter helps modestly** — ~7% bias reduction, ~5% RMSE reduction by trimming the worst MC iterations.
2. **It cannot fix the wrong KIE.** Using Saueressig when truth is 1.0046 gives ±18 Tg/yr bias regardless of filtering. The KIE bias is baked into the physics.
3. **The filter's value is as a diagnostic, not a corrector.** It identifies which KIE is more consistent with observations (via the agreement rate discriminant), but doesn't eliminate bias from choosing wrong.
4. **Source-signature uncertainty dominates the residual error.** Even with the correct KIE, RMSE ≈ 11 Tg/yr — driven by the mismatch between the fixed true signatures and the MC-sampled real distributions used in inversion.

This is why the recovered FF pattern in the figures looks different from the smooth true FF curve: the inversion uses time-varying, uncertain source signatures while the truth was generated with fixed constants.
