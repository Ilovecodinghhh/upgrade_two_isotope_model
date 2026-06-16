#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
phase9_sh_source_sensitivity.py — SH source-region sensitivity test (KIE_sites)
================================================================================

Tests the hypothesis that the wetland CH₄ reaching the two Southern-Hemisphere
constraint sites (CGO, SPO) does NOT originate only from local austral wetlands
(the "SH_extra" band assumed by Phase 6), but is a transport-weighted mixture
of upstream latitude bands (tropical wetlands dominate the *global* wetland
flux, and deep-SH air is ~1 yr old and globally mixed).

Why this matters
----------------
The headline KIE result, α¹³C_OH = 1.004 [SH-only], depends on CGO/SPO having a
negligible seasonal source phasor — which is true ONLY if their source is the
tiny local SH_extra band (2.9 Tg/yr). Phase 6 hard-codes that assignment
(`SITE_BAND`) and justifies it with a one-line interhemispheric-attenuation
estimate (τ_mix ≈ 1.3 yr → "<0.002 effect on R"). The project's own QA notes
record that an earlier *Global* assignment over-corrected by 7×. Nobody
bracketed the physically realistic middle ground. A throwaway calculation shows
flipping SH_extra → Global swings α by +0.014 (≈9× the Saueressig–Cantrell gap
of 0.0015), so this single discrete choice may dominate the whole result.

Transport model
---------------
A remote band's seasonal emission phasor, seen at an SH receptor after
interhemispheric transport with exchange timescale τ, is low-pass filtered:

    H(τ) = 1 / (1 + i·ω·τ),     ω = 2π yr⁻¹

    |H| = 1/√(1+(ωτ)²)   (amplitude attenuation)
    arg(H) = −arctan(ωτ) (phase lag, peak arrives later)

At τ = 1.3 yr: |H| = 0.122 and lag = 2.8 months — reproducing the exact numbers
quoted in phase6_phasor.py's comment block, confirming the model is faithful to
the experiment's own stated physics. Local (co-located) bands use τ = 0 → H = 1.

The effective source phasor at an SH site is the sum over emitting bands:

    Z_src(iso) = Σ_b (δ_src,b − δ_atm) · [(B_Q,b + i·C_Q,b)/Q_total] · H(τ_b)

and the corrected ratio / KIE follow the existing Phase 6 inversion:

    Z_sink = Z_obs − Z_src ;  R = |Z_sink₁₃C| / |Z_sink_D| ;  α = f(R)

Scenarios (run for CGO, SPO; appendix table for all clean sites)
----------------------------------------------------------------
  A_local_only        SH_extra only, τ=0            → reproduces Phase 6 (baseline)
  B_tropics_undamped  Tropics only, τ=0             → upper bound (tropics dominate flux)
  C_global_undamped   Global band only, τ=0         → the QA-notes "overcorrection"
  D_transport_mix     SH_extra(τ=0) + [Trop+NH](τ_mix) → physical estimate (τ_mix swept)
  E_local_plus_tropics  SH_extra + Tropics(τ_mix)   → isolates tropical contribution
  E_local_plus_NH     SH_extra + NH bands(τ_mix)    → isolates NH contribution

Outputs
-------
  results/phase9_sh_source_sensitivity/sh_source_sensitivity.json
  figures/fig18_sh_source_scenarios.png
  figures/fig19_sh_source_phasor.png

References
----------
  Li et al. (2026) ESSD — wetland emission seasonality per band
  Douglas et al. (2021) Biogeosciences — band δD source signatures
  Patra et al. (2011) ACP — interhemispheric exchange time (~1 yr)
"""

from pathlib import Path
import importlib.util
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ============================================================================
# PATHS
# ============================================================================
ANALYSIS_DIR = Path(__file__).resolve().parent
EXPT_DIR     = ANALYSIS_DIR.parent
PHASE2_JSON  = EXPT_DIR / "results" / "phase2_harmonics" / "harmonic_fits.json"
WETLAND_JSON = EXPT_DIR / "data" / "wetland_seasonality.json"
DD_SRC_JSON  = EXPT_DIR / "data" / "dD_source_database.json"
OUT_DIR      = EXPT_DIR / "results" / "phase9_sh_source_sensitivity"
FIG_DIR      = EXPT_DIR / "figures"

OUT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)
OUT_JSON = OUT_DIR / "sh_source_sensitivity.json"


# ============================================================================
# REUSE PHASE 6 HELPERS + CONSTANTS (load as a module; matches test house style)
# ============================================================================
def _load_phase6():
    spec = importlib.util.spec_from_file_location(
        "phase6_phasor", ANALYSIS_DIR / "phase6_phasor.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


phase6 = _load_phase6()

# Pull constants from phase6 so this stays in lock-step with the production run
D13C_ATM = phase6.D13C_ATM            # -47.3 ‰ VPDB
DD_ATM   = phase6.DD_ATM              # -86.0 ‰ VSMOW
Q_TOTAL_TG_MONTH        = phase6.Q_TOTAL_TG_MONTH
Q_TOTAL_TG_YR           = phase6.Q_TOTAL_TG_YR
Q_TOTAL_TG_YR_SIGMA     = phase6.Q_TOTAL_TG_YR_SIGMA
WETLAND_BC_FRAC_SIGMA   = phase6.WETLAND_BC_FRAC_SIGMA
D13C_WETLAND            = phase6.D13C_WETLAND        # -62 ‰ (global, all bands in main run)
D13C_WETLAND_SIGMA      = phase6.D13C_WETLAND_SIGMA  # 5 ‰
ALPHA_13C_SAUERESSIG    = phase6.ALPHA_13C_SAUERESSIG
ALPHA_13C_CANTRELL      = phase6.ALPHA_13C_CANTRELL

OMEGA = 2.0 * np.pi   # annual angular frequency (yr⁻¹); phasors use t in fractional years

# Band δD source signatures (Douglas 2021 zonal means). SH_extra uses the <30°
# default (−301), matching what Phase 6 assigns to CGO/SPO, so scenario A is a
# faithful reproduction of the production run.
BAND_dD = {
    "NH_high":  -374.0,
    "NH_mid":   -324.0,
    "Tropics":  -301.0,
    "SH_extra": -301.0,
    "Global":   -310.0,
}
BAND_dD_SIGMA = {
    "NH_high": 10.0, "NH_mid": 14.0, "Tropics": 15.0, "SH_extra": 15.0, "Global": 15.0,
}
# δ¹³C held at the Phase 6 global value for every band in the main scenarios, so
# this test isolates the *source-region* (δD-gradient + seasonality) effect the
# user asked about, without confounding it with a δ¹³C latitude sweep.
BAND_d13C = {b: D13C_WETLAND for b in BAND_dD}

# Saueressig central τ_mix and its spread for the MC (Patra 2011: IH exchange ~1 yr)
TAU_MIX_MEAN  = 1.3
TAU_MIX_SIGMA = 0.3
TAU_MIX_SWEEP = np.linspace(0.0, 3.0, 31)

SH_SITES = ["CGO", "SPO"]
CLEAN_SITES = list(phase6.CLEAN_SITES)
N_MC = 50_000
SEED = 42


# ============================================================================
# TRANSPORT TRANSFER FUNCTION
# ============================================================================
def transport_transfer(tau_years):
    """First-order interhemispheric low-pass H(τ) = 1/(1 + i·ω·τ).

    τ = 0  → H = 1 (local / co-located source, no attenuation or lag).
    Larger τ → stronger amplitude attenuation and larger phase lag.
    Accepts scalar or array. At τ=1.3 yr: |H|≈0.122, lag≈2.8 months.
    """
    return 1.0 / (1.0 + 1j * OMEGA * np.asarray(tau_years, dtype=float))


# ============================================================================
# SOURCE PHASOR CONSTRUCTION
# ============================================================================
def band_source_phasor(B_Q, C_Q, q_total, delta_src, delta_atm, tau, bc_scale=1.0):
    """Source phasor (‰) contributed by one emitting band, seen at the receptor.

    B_Q, C_Q : band seasonal-emission harmonic coefficients (Tg/month, month-index basis)
    q_total  : total global CH₄ source (Tg/month) — scalar or array
    delta_src, delta_atm : source and atmospheric δ (‰) for the isotope
    tau      : transport timescale (yr); 0 for local band
    bc_scale : multiplicative draw on (B_Q, C_Q) for MC ensemble spread
    Returns complex phasor (scalar or array, matching input shapes).
    """
    B_mid, C_mid = phase6.convert_wetland_to_phase2_phasor(B_Q * bc_scale, C_Q * bc_scale)
    z_frac = (B_mid + 1j * C_mid) / q_total
    return (delta_src - delta_atm) * z_frac * transport_transfer(tau)


def build_total_source(bands, wet_bands, q_total=Q_TOTAL_TG_MONTH,
                       band_dD=BAND_dD, band_d13C=BAND_d13C):
    """Sum band source phasors for a scenario (deterministic, central values).

    bands : list of (band_key, tau) tuples
    wet_bands : dict band_key → {B_Q_Tg_month, C_Q_Tg_month, ...} (from Li2026 JSON)
    Returns (Z_src_13c, Z_src_dD) as complex.
    """
    z13 = 0.0 + 0.0j
    zD = 0.0 + 0.0j
    for bkey, tau in bands:
        B_Q = wet_bands[bkey]["B_Q_Tg_month"]
        C_Q = wet_bands[bkey]["C_Q_Tg_month"]
        z13 += band_source_phasor(B_Q, C_Q, q_total, band_d13C[bkey], D13C_ATM, tau)
        zD  += band_source_phasor(B_Q, C_Q, q_total, band_dD[bkey],  DD_ATM,  tau)
    return z13, zD


def decompose(site_fit, z_src_13c, z_src_dD):
    """Vector-subtract the source phasor and invert R → α (Phase 6 inversion)."""
    z_obs_13 = complex(site_fit["d13C"]["B"], site_fit["d13C"]["C"])
    z_obs_dD = complex(site_fit["dD"]["B"], site_fit["dD"]["C"])
    z_sink_13 = z_obs_13 - z_src_13c
    z_sink_dD = z_obs_dD - z_src_dD
    R = abs(z_sink_13) / abs(z_sink_dD) if abs(z_sink_dD) > 0 else np.nan
    alpha = float(phase6.ratio_to_alpha_13c(R))
    return {
        "R_corrected": float(R),
        "alpha_13c_oh": alpha,
        "A_src_13c": float(abs(z_src_13c)),
        "A_src_dD": float(abs(z_src_dD)),
        "A_sink_13c": float(abs(z_sink_13)),
        "A_sink_dD": float(abs(z_sink_dD)),
        "peak_src_dD": float(phase6.phasor_peak_month(z_src_dD.real, z_src_dD.imag)),
    }


# ============================================================================
# SCENARIO DEFINITIONS
# ============================================================================
def scenario_bands(scenario, tau_mix):
    """Return the (band_key, tau) list for a named scenario."""
    if scenario == "A_local_only":
        return [("SH_extra", 0.0)]
    if scenario == "B_tropics_undamped":
        return [("Tropics", 0.0)]
    if scenario == "C_global_undamped":
        return [("Global", 0.0)]
    if scenario == "D_transport_mix":
        return [("SH_extra", 0.0), ("Tropics", tau_mix),
                ("NH_mid", tau_mix), ("NH_high", tau_mix)]
    if scenario == "E_local_plus_tropics":
        return [("SH_extra", 0.0), ("Tropics", tau_mix)]
    if scenario == "E_local_plus_NH":
        return [("SH_extra", 0.0), ("NH_mid", tau_mix), ("NH_high", tau_mix)]
    raise ValueError(f"Unknown scenario: {scenario}")


SCENARIOS = [
    ("A_local_only",        "SH_extra only (Phase 6 baseline)"),
    ("B_tropics_undamped",  "Tropics only, undamped (upper bound)"),
    ("C_global_undamped",   "Global band, undamped (QA overcorrection)"),
    ("D_transport_mix",     "SH local + transported Trop+NH (physical)"),
    ("E_local_plus_tropics","SH local + transported Tropics"),
    ("E_local_plus_NH",     "SH local + transported NH"),
]


# ============================================================================
# MONTE CARLO
# ============================================================================
def mc_scenario(site_fit, scenario, wet_bands, rng, tau_mix_mean=TAU_MIX_MEAN,
                tau_mix_sigma=TAU_MIX_SIGMA, n_mc=N_MC):
    """Monte-Carlo α samples for one site / scenario.

    Propagates: observed amplitude+phase (bootstrap CIs), per-band δD and δ¹³C,
    Q_total, wetland B_Q/C_Q ensemble spread, transport τ_mix, and the non-OH
    KIE / sink-fraction priors (reusing phase6.ratio_to_alpha_13c).
    Returns (R_samples, alpha_samples) arrays.
    """
    # --- observed phasor: draw amplitude and phase independently (as in phase6) ---
    def draw_obs(iso):
        B = site_fit[iso]["B"]; C = site_fit[iso]["C"]
        A = np.hypot(B, C)
        sig_A = phase6.ci_to_sigma(*site_fit[iso]["amplitude_ci95"])
        peak_ci = site_fit[iso].get(
            "peak_month_ci95",
            [site_fit[iso]["peak_month"] - 1.0, site_fit[iso]["peak_month"] + 1.0])
        sig_ph = phase6.ci_to_sigma(*peak_ci) * (2 * np.pi / 12)
        phase0 = np.arctan2(C, B)
        a = np.maximum(rng.normal(A, sig_A, n_mc), 1e-6)
        p = phase0 + rng.normal(0, sig_ph, n_mc)
        return a * np.cos(p) + 1j * a * np.sin(p)

    z_obs_13 = draw_obs("d13C")
    z_obs_dD = draw_obs("dD")

    # --- shared draws ---
    q_tot = np.maximum(rng.normal(Q_TOTAL_TG_YR, Q_TOTAL_TG_YR_SIGMA, n_mc), 300.0) / 12.0
    d13c_w = rng.normal(D13C_WETLAND, D13C_WETLAND_SIGMA, n_mc)  # one δ¹³C for all bands

    bands = scenario_bands(scenario, tau_mix_mean)

    z_src_13 = np.zeros(n_mc, dtype=complex)
    z_src_dD = np.zeros(n_mc, dtype=complex)
    for bkey, tau_central in bands:
        B_Q = wet_bands[bkey]["B_Q_Tg_month"]
        C_Q = wet_bands[bkey]["C_Q_Tg_month"]
        bc = 1.0 + rng.normal(0, WETLAND_BC_FRAC_SIGMA, n_mc)
        dD_b = rng.normal(BAND_dD[bkey], BAND_dD_SIGMA[bkey], n_mc)
        # local band τ fixed at 0; remote bands draw τ_mix (shared per draw)
        if tau_central == 0.0:
            tau = np.zeros(n_mc)
        else:
            tau = np.clip(rng.normal(tau_mix_mean, tau_mix_sigma, n_mc), 0.3, 3.5)
        z_src_13 += band_source_phasor(B_Q, C_Q, q_tot, d13c_w, D13C_ATM, tau, bc)
        z_src_dD += band_source_phasor(B_Q, C_Q, q_tot, dD_b, DD_ATM, tau, bc)

    z_sink_13 = z_obs_13 - z_src_13
    z_sink_dD = z_obs_dD - z_src_dD
    A_sink_dD = np.abs(z_sink_dD)
    R = np.divide(np.abs(z_sink_13), A_sink_dD,
                  out=np.full(n_mc, np.nan), where=A_sink_dD > 1e-6)

    # non-OH KIE / sink-fraction priors (identical to phase6.mc_phasor)
    f_oh = np.clip(rng.normal(phase6.F_OH, phase6.SIGMA_F_OH, n_mc), 0.5, 0.99)
    f_cl = np.clip(rng.normal(phase6.F_CL, phase6.SIGMA_F_CL, n_mc), 0.0, 0.1)
    f_soil = np.clip(rng.normal(phase6.F_SOIL, phase6.SIGMA_F_SOIL, n_mc), 0.0, 0.15)
    f_strat = 1.0 - f_oh - f_cl - f_soil
    alpha_d_oh = rng.normal(phase6.ALPHA_D_OH, phase6.SIGMA_ALPHA_D_OH, n_mc)
    alpha_13c_cl = rng.normal(phase6.ALPHA_13C_CL, phase6.SIGMA_ALPHA_13C_CL, n_mc)
    alpha_d_cl = rng.normal(phase6.ALPHA_D_CL, phase6.SIGMA_ALPHA_D_CL, n_mc)

    alpha = phase6.ratio_to_alpha_13c(
        R, alpha_d_oh, f_oh, f_cl, f_soil, f_strat, alpha_13c_cl, alpha_d_cl)
    return R, alpha


def mc_summary(alpha_samples):
    """Median + 95% CI over finite, physically-filtered α (wide filter, as phase6)."""
    a = np.asarray(alpha_samples)
    a = a[np.isfinite(a)]
    a = a[(a > 0.98) & (a < 1.05)]
    if len(a) < 100:
        return {"alpha_median": np.nan, "alpha_ci95": [np.nan, np.nan], "n": int(len(a))}
    return {
        "alpha_median": float(np.median(a)),
        "alpha_ci95": [float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))],
        "n": int(len(a)),
    }


# ============================================================================
# MAIN
# ============================================================================
def main():
    print("=" * 74)
    print("Phase 9: SH source-region sensitivity (do CGO/SPO see other bands?)")
    print("=" * 74)

    with open(PHASE2_JSON) as f:
        fits = json.load(f)
    with open(WETLAND_JSON) as f:
        wet = json.load(f)["bands"]

    results = {
        "metadata": {
            "purpose": "Test whether the SH KIE constraint (CGO, SPO) is sensitive "
                       "to the assumed wetland source region / latitude-band mixing.",
            "transport_model": "H(tau) = 1/(1 + i*omega*tau), omega = 2*pi/yr; "
                               "local band tau=0, remote bands tau=tau_mix.",
            "tau_mix_mean": TAU_MIX_MEAN, "tau_mix_sigma": TAU_MIX_SIGMA,
            "n_mc": N_MC, "seed": SEED,
            "band_dD": BAND_dD, "band_d13C": BAND_d13C,
            "saueressig": ALPHA_13C_SAUERESSIG, "cantrell": ALPHA_13C_CANTRELL,
            "note": "Diagnostic only — does not modify Phase 6 SITE_BAND or the "
                    "production headline. delta13C held at -62 per mil for all bands "
                    "to isolate the source-region (delta-D + seasonality) effect.",
        },
        "transport_check": {},
        "scenarios": {},
        "tau_sweep": {},
        "pooled_sh": {},
        "all_sites_appendix": {},
    }

    # --- transport sanity check (reproduce phase6's quoted numbers) ---
    H13 = transport_transfer(1.3)
    results["transport_check"] = {
        "tau_yr": 1.3, "abs_H": float(abs(H13)),
        "lag_months": float((-np.angle(H13)) / (2 * np.pi) * 12),
        "phase6_comment_says": "abs≈0.11-0.12, lag≈2.8 months",
    }
    print(f"\nTransport check  τ=1.3 yr:  |H|={abs(H13):.3f}  lag={(-np.angle(H13))/(2*np.pi)*12:.2f} mo "
          f"(phase6 comment: 0.11–0.12, 2.8 mo)")

    rng = np.random.default_rng(SEED)

    # --- per-scenario deterministic + MC for SH sites ---
    print(f"\n{'scenario':24} {'site':4} {'R_corr':>7} {'alpha':>7} {'95% CI':>20} "
          f"{'|Zsrc_dD|':>9}")
    print("-" * 80)
    for scen, label in SCENARIOS:
        results["scenarios"][scen] = {"label": label, "sites": {}}
        for site in SH_SITES:
            z13, zD = build_total_source(scenario_bands(scen, TAU_MIX_MEAN), wet)
            det = decompose(fits[site], z13, zD)
            _, alpha_mc = mc_scenario(fits[site], scen, wet, rng)
            summ = mc_summary(alpha_mc)
            results["scenarios"][scen]["sites"][site] = {**det, "mc": summ}
            print(f"{scen:24} {site:4} {det['R_corrected']:7.4f} "
                  f"{summ['alpha_median']:7.4f} "
                  f"[{summ['alpha_ci95'][0]:.4f},{summ['alpha_ci95'][1]:.4f}]  "
                  f"{det['A_src_dD']:9.3f}")

        # pooled SH (concatenate CGO+SPO α samples, wide filter)
        pooled = []
        rng_p = np.random.default_rng(SEED + 1)
        for site in SH_SITES:
            _, a = mc_scenario(fits[site], scen, wet, rng_p)
            pooled.append(a)
        results["scenarios"][scen]["pooled_sh"] = mc_summary(np.concatenate(pooled))

    # --- τ_mix sweep for the physical scenario D ---
    print(f"\nτ_mix sweep (scenario D_transport_mix), deterministic α:")
    print(f"{'tau_yr':>7} {'CGO_R':>7} {'CGO_a':>7} {'SPO_R':>7} {'SPO_a':>7}")
    for tau in TAU_MIX_SWEEP:
        row = {}
        for site in SH_SITES:
            z13, zD = build_total_source(scenario_bands("D_transport_mix", tau), wet)
            row[site] = decompose(fits[site], z13, zD)
        results["tau_sweep"][f"{tau:.2f}"] = {
            s: {"R_corrected": row[s]["R_corrected"],
                "alpha_13c_oh": row[s]["alpha_13c_oh"]} for s in SH_SITES}
        if abs(tau - round(tau, 1)) < 1e-9 and (round(tau * 10) % 5 == 0):
            print(f"{tau:7.2f} {row['CGO']['R_corrected']:7.4f} {row['CGO']['alpha_13c_oh']:7.4f} "
                  f"{row['SPO']['R_corrected']:7.4f} {row['SPO']['alpha_13c_oh']:7.4f}")

    # --- appendix: all clean sites under baseline (A) vs physical mix (D) ---
    for site in CLEAN_SITES:
        zA13, zAD = build_total_source(scenario_bands("A_local_only", TAU_MIX_MEAN), wet)
        # NH/tropical sites: "local" band is their own assigned band, not SH_extra.
        # For the appendix we only report SH sites' mix faithfully; others get A only.
        if site in SH_SITES:
            zD13, zDD = build_total_source(scenario_bands("D_transport_mix", TAU_MIX_MEAN), wet)
            results["all_sites_appendix"][site] = {
                "A_local_only": decompose(fits[site], zA13, zAD),
                "D_transport_mix": decompose(fits[site], zD13, zDD),
            }

    # --- pooled-SH summary across scenarios (the headline read) ---
    print(f"\n{'='*60}\nPOOLED SH (CGO+SPO) α BY SCENARIO\n{'='*60}")
    for scen, label in SCENARIOS:
        p = results["scenarios"][scen]["pooled_sh"]
        results["pooled_sh"][scen] = p
        print(f"  {scen:24} α = {p['alpha_median']:.4f} "
              f"[{p['alpha_ci95'][0]:.4f}, {p['alpha_ci95'][1]:.4f}]")
    print(f"\n  Saueressig = {ALPHA_13C_SAUERESSIG};  Cantrell = {ALPHA_13C_CANTRELL}")

    with open(OUT_JSON, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n✓ Saved: {OUT_JSON}")

    plot_scenarios(results)
    plot_phasor(results, fits)


# ============================================================================
# FIGURES
# ============================================================================
def plot_scenarios(results):
    """Fig 18: (a) α per scenario for CGO/SPO/pooled; (b) α vs τ_mix (scenario D)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    scen_keys = [s for s, _ in SCENARIOS]
    x = np.arange(len(scen_keys))

    for off, site, col in [(-0.18, "CGO", "C0"), (0.0, "SPO", "C1"), (0.18, "pooled", "k")]:
        meds, los, his = [], [], []
        for scen in scen_keys:
            if site == "pooled":
                m = results["scenarios"][scen]["pooled_sh"]
            else:
                m = results["scenarios"][scen]["sites"][site]["mc"]
            meds.append(m["alpha_median"])
            los.append(m["alpha_median"] - m["alpha_ci95"][0])
            his.append(m["alpha_ci95"][1] - m["alpha_median"])
        ax1.errorbar(x + off, meds, yerr=[los, his], fmt="o", color=col, capsize=3,
                     ms=6, lw=1.5, label=("pooled SH" if site == "pooled" else site))

    ax1.axhspan(ALPHA_13C_SAUERESSIG, ALPHA_13C_CANTRELL, alpha=0.15, color="green",
                label="Saueressig–Cantrell")
    ax1.axhline(ALPHA_13C_SAUERESSIG, color="green", ls="--", lw=0.8)
    ax1.axhline(ALPHA_13C_CANTRELL, color="green", ls="--", lw=0.8)
    ax1.set_xticks(x)
    ax1.set_xticklabels([s.replace("_", "\n") for s in scen_keys], fontsize=7)
    ax1.set_ylabel("α¹³C_OH", fontsize=11)
    ax1.set_title("(a) SH KIE constraint by source-region scenario", fontsize=11)
    ax1.legend(fontsize=8, loc="upper left")

    # (b) τ_mix sweep
    taus = np.array([float(k) for k in results["tau_sweep"].keys()])
    order = np.argsort(taus)
    taus = taus[order]
    for site, col in [("CGO", "C0"), ("SPO", "C1")]:
        a = np.array([results["tau_sweep"][f"{t:.2f}"][site]["alpha_13c_oh"] for t in taus])
        ax2.plot(taus, a, "-", color=col, lw=2, label=site)
    ax2.axhspan(ALPHA_13C_SAUERESSIG, ALPHA_13C_CANTRELL, alpha=0.15, color="green")
    ax2.axhline(ALPHA_13C_SAUERESSIG, color="green", ls="--", lw=0.8, label="Saueressig")
    ax2.axhline(ALPHA_13C_CANTRELL, color="green", ls="--", lw=0.8, label="Cantrell")
    ax2.axvline(TAU_MIX_MEAN, color="gray", ls=":", lw=1, label=f"τ_mix={TAU_MIX_MEAN} yr")
    ax2.set_xlabel("Interhemispheric transport time τ_mix (yr)", fontsize=11)
    ax2.set_ylabel("α¹³C_OH", fontsize=11)
    ax2.set_title("(b) Scenario D: α vs transport time", fontsize=11)
    ax2.legend(fontsize=8, loc="upper right")

    fig.suptitle("Fig 18: SH source-region sensitivity of the OH ¹³C KIE constraint",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = FIG_DIR / "fig18_sh_source_scenarios.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"✓ Figure saved: {out}")


def plot_phasor(results, fits):
    """Fig 19: (a) δD source-phasor amplitude per scenario; (b) R_corrected per scenario."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
    scen_keys = [s for s, _ in SCENARIOS]
    x = np.arange(len(scen_keys))
    width = 0.38

    for off, site, col in [(-width / 2, "CGO", "C0"), (width / 2, "SPO", "C1")]:
        a_src = [results["scenarios"][s]["sites"][site]["A_src_dD"] for s in scen_keys]
        ax1.bar(x + off, a_src, width, color=col, alpha=0.8, label=site)
    # observed δD amplitude reference lines
    for site, col in [("CGO", "C0"), ("SPO", "C1")]:
        A_obs = np.hypot(fits[site]["dD"]["B"], fits[site]["dD"]["C"])
        ax1.axhline(A_obs, color=col, ls=":", lw=1.2, alpha=0.7,
                    label=f"{site} |Z_obs(δD)|={A_obs:.2f}")
    ax1.set_xticks(x)
    ax1.set_xticklabels([s.replace("_", "\n") for s in scen_keys], fontsize=7)
    ax1.set_ylabel("Source phasor |Z_src(δD)|  (‰)", fontsize=11)
    ax1.set_title("(a) Wetland source amplitude vs observed δD signal", fontsize=11)
    ax1.legend(fontsize=7)

    for off, site, col in [(-width / 2, "CGO", "C0"), (width / 2, "SPO", "C1")]:
        r = [results["scenarios"][s]["sites"][site]["R_corrected"] for s in scen_keys]
        ax2.bar(x + off, r, width, color=col, alpha=0.8, label=site)
    # R bands for the two KIE values (bulk inversion)
    r_sau = _bulk_R(ALPHA_13C_SAUERESSIG)
    r_can = _bulk_R(ALPHA_13C_CANTRELL)
    ax2.axhspan(r_sau, r_can, alpha=0.15, color="green", label=f"Sau–Can R ({r_sau:.3f}–{r_can:.3f})")
    ax2.set_xticks(x)
    ax2.set_xticklabels([s.replace("_", "\n") for s in scen_keys], fontsize=7)
    ax2.set_ylabel("R_corrected = A(δ¹³C)/A(δD)", fontsize=11)
    ax2.set_title("(b) Corrected ratio by scenario", fontsize=11)
    ax2.legend(fontsize=7)

    fig.suptitle("Fig 19: How source-region choice reshapes the SH source phasor and R",
                 fontsize=12, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = FIG_DIR / "fig19_sh_source_phasor.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"✓ Figure saved: {out}")


def _bulk_R(alpha_13c_oh):
    """Bulk sink amplitude ratio R for a candidate OH ¹³C KIE (for figure bands)."""
    e13 = (phase6.F_OH * (alpha_13c_oh - 1) + phase6.F_CL * (phase6.ALPHA_13C_CL - 1)
           + phase6.F_SOIL * (phase6.ALPHA_13C_SOIL - 1)
           + phase6.F_STRAT * (phase6.ALPHA_13C_STRAT - 1)) * 1000
    eD = (phase6.F_OH * (phase6.ALPHA_D_OH - 1) + phase6.F_CL * (phase6.ALPHA_D_CL - 1)
          + phase6.F_SOIL * (phase6.ALPHA_D_SOIL - 1)
          + phase6.F_STRAT * (phase6.ALPHA_D_STRAT - 1)) * 1000
    return e13 / eD


if __name__ == "__main__":
    main()
