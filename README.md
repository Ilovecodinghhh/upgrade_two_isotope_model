# Methane Isotope Box Models — Clean Taxonomy

## Model Infra (2 × 3 = 6 models)

|  | **One Box (Global)** | **Two Boxes (NH + SH)** | **Three Boxes (NHext / Trop / SHext)** |
|---|---|---|---|
| **2×2** (BB fixed; δ¹³C and δD solved separately) | `22_one.py` | `22_two.py` | `2x2_three.py` |
| **3×3** (δ¹³C and δD used simultaneously) | `33_one.py` | `33_two.py` | `3x3_three.py` |

### What "2×2" means
Fix biomass-burning from CarbonTracker (GFED4 prior). Solve for FF and Mic using δ¹³C alone, then independently using δD alone. Cross-validate the two estimates.

### What "3×3" means
Solve a 3-equation system (mass + δ¹³C + δD) simultaneously for all three sources (BB, FF, Mic). No external BB constraint needed. Uses bounded least-squares for non-negativity.

### What "one" vs "two" vs "three" means
- **one**: Single global well-mixed box. τ = global lifetime.
- **two**: Two hemispheric boxes (NH/SH) connected by interhemispheric exchange (τ_ex ~ 1.0 ± 0.1 yr). Hemisphere-specific lifetimes, sink fractions, and real hemispheric δD.
- **three**: Three latitude-band boxes (NHext >30°N, Trop 30°S–30°N, SHext <30°S) with two exchange paths (NHext↔Trop τ_NT ~ 0.8 yr; Trop↔SHext τ_TS ~ 1.2 yr). Box-specific sink fractions, lifetimes, BB splits, and real per-box δD observations.

---

## Shared Infrastructure (`models/`)

| File | Contents |
|------|----------|
| `models/__init__.py` | Package exports |
| `models/inputs.py` | **Input choice catalog** — every configurable parameter with literature references and sensitivity presets |
| `models/core.py` | Data loading, KIE sampling, isotope math, lifetime, quality monitoring, smoothing |

---

## Input Choices (for sensitivity tests)

All documented with references in `models/inputs.py`. Run any model with a preset:

```bash
python 22_one.py default           # Default config
python 22_one.py cantrell_only     # OH KIE fixed to Cantrell (1990)
python 33_two.py fixed_lifetime    # τ = 9 yr constant
python 22_two.py CTCH4_FF          # CT-CH₄ fossil fuel signatures
```

### Available presets:

| Preset | What changes | Reference |
|--------|-------------|-----------|
| `default` | Full stochastic (sampled KIE + time-varying τ + EDGAR FF) | This work |
| `cantrell_only` | OH KIE = 1.0054 (¹³C), 1.327 (D) | Cantrell et al. (1990) |
| `saueressig_only` | OH KIE = 1.0039 (¹³C), 1.294 (D) | Saueressig et al. (2001) |
| `KIE_strat_soil_fixed` | Only OH and Cl sampled; Strat/Soil fixed | — |
| `fixed_lifetime` | τ = 9.0 yr for all years | IPCC AR5/AR6 |
| `CTCH4_FF` | CarbonTracker-CH₄ posterior FF signatures | Bruhwiler et al. (2014) |
| `mic_dD_7` | Microbial δD uncertainty = 7‰ (original) | Bao et al. |
| `BB_annual` | Time-varying BB from CT (not fixed mean) | GFED4 |
| `BB_declining` | Declining BB scenario | Worden et al. (2017) |
| `thanwerdas_sinks` | Thanwerdas sink fractions (low Cl) | Thanwerdas et al. (2024) |
| `tau_ex_fixed` | τ_ex = 1.0 yr (no sampling) | Patra et al. (2011) |
| `dD_offset_old` | δD NH/SH offset = ±1.5‰ (old estimate) | — |

---

## Not included here (set aside)

| Model | Reason |
|-------|--------|
| `v4.0_mic_vs_nonmic.py` | Different source taxonomy (2-source: Mic vs NonMic) |
| Bayesian MCMC/NUTS | Different statistical framework (PyMC) — separate branch |
| Frequentist scenarios (RedCl, IncOH, BBdrop) | Will be re-implemented as sensitivity presets |

---

## Data Lineage

All input data files in `TwoIsotopeBoxModel/rel/`:
- `data/` — atmospheric observations (CH₄, δ¹³C, δD, CarbonTracker)
- `output/` — pre-computed MC source signatures (1000 iterations × years)

See `models/inputs.py` → `ATM_OBS_CATALOG` and `FF_SIGNATURE_OPTIONS` for full documentation.

---

## Version History

| Date | Change |
|------|--------|
| 2026-05-12 | Three-box models (`2x2_three.py`, `3x3_three.py`) on `three-box` branch |
| 2026-05-12 | Real hemispheric δD upgrade in `common.py` (v2 KIE immunity experiment) |
| 2026-05-09 | Clean taxonomy: 4 models + shared `models/` infra + `inputs.py` catalog |
| 2026-05-05 | v3.1/v3.2 with Ben-model optimisations |
| 2026-05-04 | v2.0 upgrades (KIE sampling, time-varying τ, quality monitor) |
| 2026-05-04 | v1.0 original Bao code |
