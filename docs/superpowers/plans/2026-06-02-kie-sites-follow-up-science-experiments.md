# KIE_sites Follow-up Science Experiments Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add follow-up KIE_sites experiments that test whether the current seasonal-phasor method can recover known alpha values, quantify uncertainty drivers, and diagnose remaining NH bias without expanding into a full multi-source transport model.

**Architecture:** Keep Phase 1-8 immutable. Add standalone Phase 9-13 scripts that read existing KIE_sites results, write new outputs under `results/phase9_*` through `results/phase13_*`, and create optional diagnostic figures under `figures/`. Each phase exposes small pure functions for tests plus a `main()` entry point for reproducible outputs.

**Tech Stack:** Python, numpy, pandas, matplotlib, json, pathlib, pytest.

---

## File Structure

- Create `experiments/KIE_sites/analysis/phase9_osse.py`
  - Synthetic observation system simulation experiment (OSSE) for known alpha recovery.
- Create `experiments/KIE_sites/analysis/phase10_phasor_inversion.py`
  - Direct complex-phasor inversion and residual diagnostics.
- Create `experiments/KIE_sites/analysis/phase11_block_bootstrap.py`
  - Year-block bootstrap for alpha uncertainty.
- Create `experiments/KIE_sites/analysis/phase12_harmonic_sensitivity.py`
  - Annual, annual-plus-semiannual, monthly fixed-effect, and leave-one-year-out sensitivity checks.
- Create `experiments/KIE_sites/analysis/phase13_uncertainty_attribution.py`
  - One-at-a-time uncertainty attribution for alpha variance and NH residual diagnostics.
- Create tests:
  - `experiments/KIE_sites/tests/test_phase9_osse.py`
  - `experiments/KIE_sites/tests/test_phase10_phasor_inversion.py`
  - `experiments/KIE_sites/tests/test_phase11_block_bootstrap.py`
  - `experiments/KIE_sites/tests/test_phase12_harmonic_sensitivity.py`
  - `experiments/KIE_sites/tests/test_phase13_uncertainty_attribution.py`

Do not modify Phase 1-8 scripts or overwrite their result files.

## Task 1: OSSE Synthetic Recovery

- [ ] Add tests proving pure-OH synthetic data recovers the input alpha and wetland contamination biases the uncorrected scalar ratio high.
- [ ] Implement a compact monthly synthetic generator with explicit sink/source phasors, deterministic random seed support, and alpha inversion.
- [ ] Write `osse_results.json` and `fig18_osse_recovery.png`.

## Task 2: Direct Phasor Inversion

- [ ] Add tests proving a known observed/source phasor decomposes into the expected sink phasor and alpha.
- [ ] Implement per-site residual phasor diagnostics from existing Phase 6 and Phase 8 outputs.
- [ ] Write `phasor_inversion_results.json` and `fig19_phasor_inversion_diagnostics.png`.

## Task 3: Year-block Bootstrap

- [ ] Add tests proving fixed-seed year-block resampling is reproducible and preserves full years.
- [ ] Implement yearly phasor refits from Phase 1 monthly paired data.
- [ ] Write `block_bootstrap_results.json` and `fig20_block_bootstrap_alpha.png`.

## Task 4: Harmonic Model Sensitivity

- [ ] Add tests proving annual-plus-semiannual fitting preserves the annual amplitude on pure annual data.
- [ ] Implement annual, annual-plus-semiannual, monthly fixed-effect, and leave-one-year-out comparisons.
- [ ] Write `harmonic_sensitivity_results.json` and `fig21_harmonic_model_comparison.png`.

## Task 5: Uncertainty Attribution and NH Residual Diagnostics

- [ ] Add tests proving one-at-a-time variance attribution only flags the active uncertainty source.
- [ ] Implement grouped perturbation attribution for observation, wetland phasor, wetland isotope signatures, BB correction, sink fractions, alpha_D, and non-OH KIE.
- [ ] Write `uncertainty_attribution_results.json` and `fig22_uncertainty_attribution.png`.

## Test Plan

Run:

```bash
python -m pytest experiments/KIE_sites/tests -q
```

Then run:

```bash
python experiments/KIE_sites/analysis/phase9_osse.py
python experiments/KIE_sites/analysis/phase10_phasor_inversion.py
python experiments/KIE_sites/analysis/phase11_block_bootstrap.py
python experiments/KIE_sites/analysis/phase12_harmonic_sensitivity.py
python experiments/KIE_sites/analysis/phase13_uncertainty_attribution.py
```

Acceptance criteria:

- All tests pass.
- Each new phase writes a JSON result file.
- New phases do not modify Phase 1-8 scripts or result files.
- OSSE pure-OH recovery is within 0.001 of the input alpha.
- SH-only diagnostics remain close to the Phase 6 interpretation, while NH residual diagnostics are explicitly reported as still biased high.

## Assumptions

- Existing Phase 6 SH-only alpha remains the scientific baseline.
- Wetland and BB are the only source corrections used in this implementation.
- No transport lag, amplitude attenuation, or continuous clean-site weighting is introduced.
- New scripts are diagnostic/scientific experiments, not replacements for the existing published Phase 6 result.
