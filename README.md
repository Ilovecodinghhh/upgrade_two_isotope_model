# Upgraded Two-Isotope Box Model

Dual-isotope (δ¹³C + δD) Monte Carlo mass-balance model for global methane source partitioning (BB/FF/Mic).

Upgraded from [TwoIsotopeBoxModel](https://github.com/Ilovecodinghhh/TwoIsotopeBoxModel) by Yufan Bao (2026).

## Upgrades

| # | Feature | Original | Upgraded |
|---|---------|----------|----------|
| 1 | **KIE sampling** | Fixed per scenario | Sampled from literature distributions per MC iteration |
| 2 | **Quality monitoring** | None | Condition number tracking, non-physical solution rates |
| 3 | **Lifetime** | Fixed τ = 9 yr | Time-varying τ(t) ≈ 9.0 − 0.017·(t − 2010) |

See [CHANGELOG.md](CHANGELOG.md) for detailed scientific justification and improvement suggestions.

## Quick Start

```bash
# Run the model
python upgraded_box_model.py

# Run with debug output (stops at first non-physical solution)
python upgraded_box_model.py --debug
```

## Output Files

| File | Description |
|------|-------------|
| `upgraded_base_results.csv` | Per-year statistics (mean ± std for BB/FF/Mic) |
| `quality_report.json` | Solution quality summary (rejection rates, condition numbers) |
| `quality_per_year.csv` | Per-year quality diagnostics |
| `KIE_samples.csv` | All 1000 KIE draws (for posterior analysis) |
| `BB/FF/Mic_upgraded_MC.csv` | Full 1000-iteration results per source |
| `upgraded_model_results.png` | Spaghetti + quality diagnostic plots |
| `lifetime_trajectory.png` | Time-varying lifetime visualization |

## Data Requirements

Input data in `rel/data/` and `rel/output/` — copied from the original repository.

### Simulated Data (needs replacement)

- **Time-varying lifetime**: Currently a linear parameterization. Needs actual values from He et al. (2026, Science) and Montzka et al. (2011) MCF-derived OH.
- **Microbial δD uncertainty**: Hardcoded at 7‰. Should be derived from EMID database (Menoud 2022).
- **δD 1999–2004**: Padded by repeating 2005 value. Needs real observations.

## References

- Cantrell et al. (1990) — OH KIE for ¹³C = 1.0054
- Saueressig et al. (2001) — OH KIE for ¹³C = 1.0039, D = 1.294
- Chandra et al. (2024, Comm. Earth Environ.) — KIE sensitivity analysis
- He et al. (2026, Science) — Time-varying methane lifetime
- Nguyen et al. (2020, GRL) — CH₄-OH feedback and perturbation lifetime
- Naus et al. (2019, ACP) — Two-box model framework
