#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create manuscript-ready figures for the KIE sites draft.

The figures in this script are redrawn from existing derived JSON outputs so
that manuscript numbering does not conflict with exploratory phase figures.
"""

from __future__ import annotations

from pathlib import Path
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np


EXPT_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = EXPT_DIR / "results"
OUT_DIR = EXPT_DIR / "figures" / "manuscript"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SITE_ORDER = ["SPO", "CGO", "SMO", "ASC", "MLO", "KUM", "AZR", "CBA", "MHD", "BRW", "ZEP", "ALT"]
CLEAN_SITES = ["ALT", "ZEP", "BRW", "CBA", "MHD", "KUM", "CGO", "SPO"]
SH_SITES = ["CGO", "SPO"]
MONTH_NAMES = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

ALPHA_13C_SAUERESSIG = 1.0039
ALPHA_13C_CANTRELL = 1.0054
ALPHA_D_OH = 1.294
D13C_ATM = -47.3
D13C_WETLAND_BASE = -62.0

F_OH = 0.84
F_CL = 0.035
F_SOIL = 0.06
F_STRAT = 0.065
ALPHA_13C_CL = 1.066
ALPHA_D_CL = 1.508
ALPHA_13C_SOIL = 1.022
ALPHA_D_SOIL = 1.066
ALPHA_13C_STRAT = 1.013
ALPHA_D_STRAT = 1.16


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def ci_error(value: float, ci95: list[float]) -> list[list[float]]:
    return [[value - ci95[0]], [ci95[1] - value]]


def oh_only_ratio(alpha_13c: float) -> float:
    return (alpha_13c - 1.0) / (ALPHA_D_OH - 1.0)


def bulk_sink_ratio(alpha_13c_oh: float) -> float:
    eps_13c = (
        F_OH * (alpha_13c_oh - 1.0)
        + F_CL * (ALPHA_13C_CL - 1.0)
        + F_SOIL * (ALPHA_13C_SOIL - 1.0)
        + F_STRAT * (ALPHA_13C_STRAT - 1.0)
    ) * 1000.0
    eps_d = (
        F_OH * (ALPHA_D_OH - 1.0)
        + F_CL * (ALPHA_D_CL - 1.0)
        + F_SOIL * (ALPHA_D_SOIL - 1.0)
        + F_STRAT * (ALPHA_D_STRAT - 1.0)
    ) * 1000.0
    return eps_13c / eps_d


def banded_wetland_d13c_for_source_band(source_band: str) -> float:
    if source_band in {"NH_high", "NH_mid"}:
        return -67.8
    if source_band == "Tropics":
        return -56.7
    return D13C_WETLAND_BASE


def ganesan_banded_ratio(site_result: dict) -> float:
    banded_d13c = banded_wetland_d13c_for_source_band(site_result["source_band"])
    scale = (banded_d13c - D13C_ATM) / (D13C_WETLAND_BASE - D13C_ATM)
    z_obs_13c = complex(*site_result["Z_obs_13c"])
    z_src_13c = complex(*site_result["Z_src_13c"]) * scale
    z_sink_13c = z_obs_13c - z_src_13c
    z_sink_dD = complex(*site_result["Z_sink_dD"])
    return float(abs(z_sink_13c) / abs(z_sink_dD))


def set_common_style():
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "figure.dpi": 120,
            "savefig.dpi": 300,
        }
    )


def label_points(ax, xs, ys, labels, dx=1.5, dy=0.0):
    for x, y, label in zip(xs, ys, labels):
        ax.text(x + dx, y + dy, label, fontsize=8, va="center")


def shortest_phase_diff(a: float, b: float) -> float:
    diff = abs(a - b)
    return min(diff, 12.0 - diff)


def mass_conserving_scenario_name(f_nh: float, f_tropics: float) -> str:
    return f"mc_nh{int(round(100 * f_nh)):02d}_tr{int(round(100 * f_tropics)):02d}"


def mass_conserving_grid(phase14: dict, field: str) -> tuple[list[float], list[float], np.ndarray]:
    mass = phase14["mass_conserving"]
    nh_values = mass["nh_high_fractions"]
    tropics_values = mass["tropics_fractions"]
    grid = np.full((len(nh_values), len(tropics_values)), np.nan)
    for i, f_nh in enumerate(nh_values):
        for j, f_tropics in enumerate(tropics_values):
            name = mass_conserving_scenario_name(f_nh, f_tropics)
            grid[i, j] = mass["summary_by_scenario"][name][field]
    return nh_values, tropics_values, grid


def plot_mass_conserving_heatmap(
    ax,
    phase14: dict,
    field: str,
    title: str,
    mark_lab_range: bool = False,
):
    nh_values, tropics_values, grid = mass_conserving_grid(phase14, field)
    image = ax.imshow(grid, origin="lower", aspect="auto", cmap="viridis")
    ax.set_xticks(np.arange(len(tropics_values)))
    ax.set_xticklabels([f"{value:.2f}" for value in tropics_values], rotation=35, ha="right")
    ax.set_yticks(np.arange(len(nh_values)))
    ax.set_yticklabels([f"{value:.2f}" for value in nh_values])
    ax.set_xlabel("Tropics weight")
    ax.set_ylabel("Delayed NH high weight")
    ax.set_title(title)
    midpoint = 0.5 * (float(np.nanmin(grid)) + float(np.nanmax(grid)))
    for i in range(len(nh_values)):
        for j in range(len(tropics_values)):
            in_lab_range = (
                mark_lab_range
                and ALPHA_13C_SAUERESSIG <= grid[i, j] <= ALPHA_13C_CANTRELL
            )
            label = f"{grid[i, j]:.4f}" + ("*" if in_lab_range else "")
            color = "white" if grid[i, j] < midpoint else "black"
            ax.text(j, i, label, ha="center", va="center", fontsize=6.8, color=color)
            if in_lab_range:
                ax.add_patch(
                    Rectangle(
                        (j - 0.5, i - 0.5),
                        1.0,
                        1.0,
                        fill=False,
                        edgecolor="white",
                        linewidth=1.2,
                    )
                )
    return image


def mass_conserving_slice_values(phase14: dict, slice_key: str, field: str):
    mass = phase14["mass_conserving"][slice_key]
    fractions = mass["fractions"]
    summary = mass["summary_by_scenario"]
    if slice_key == "tropics_only":
        names = [f"tr_only_{int(round(100 * fraction)):02d}" for fraction in fractions]
    elif slice_key == "nh_high_only":
        names = [f"nh_high_only_{int(round(100 * fraction)):02d}" for fraction in fractions]
    else:
        raise ValueError(f"Unknown mass-conserving slice: {slice_key}")
    values = [summary[name][field] for name in names]
    return fractions, values


def plot_mass_conserving_slice(
    ax,
    phase14: dict,
    slice_key: str,
    title: str,
    xlabel: str,
    color: str,
    reference_alpha: float | None = None,
):
    fractions, values = mass_conserving_slice_values(
        phase14,
        slice_key,
        "mean_alpha_13C_OH",
    )
    ax.axhspan(
        ALPHA_13C_SAUERESSIG,
        ALPHA_13C_CANTRELL,
        color="#b7e1b4",
        alpha=0.35,
        label="Lab range",
    )
    ax.plot(fractions, values, marker="o", color=color, lw=1.8)
    if reference_alpha is not None:
        ax.axhline(
            reference_alpha,
            color="#222222",
            ls=":",
            lw=1.3,
            label="Main SH MC central",
        )
    ax.axhline(ALPHA_13C_SAUERESSIG, color="#2ca02c", ls="--", lw=1.0)
    ax.axhline(ALPHA_13C_CANTRELL, color="#d62728", ls="--", lw=1.0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Deterministic mean alpha13C_OH")
    ax.set_title(title)
    ax.set_xticks(fractions)
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=7)


def make_fig1_raw_ratio_phase():
    phase2 = load_json(RESULTS_DIR / "phase2_harmonics" / "harmonic_fits.json")
    summary = load_json(RESULTS_DIR / "phase1_data" / "site_summary.json")
    sites = sorted(summary.keys(), key=lambda s: summary[s]["latitude"])

    lat = np.array([summary[s]["latitude"] for s in sites])
    ratios = np.array([phase2[s]["ratio"]["value"] for s in sites])
    ratio_err = np.array(
        [
            [phase2[s]["ratio"]["value"] - phase2[s]["ratio"]["ci95"][0] for s in sites],
            [phase2[s]["ratio"]["ci95"][1] - phase2[s]["ratio"]["value"] for s in sites],
        ]
    )
    phase = np.array([phase2[s]["phase_diff_months"]["value"] for s in sites])
    phase_err = np.array(
        [
            [phase2[s]["phase_diff_months"]["value"] - phase2[s]["phase_diff_months"]["ci95"][0] for s in sites],
            [phase2[s]["phase_diff_months"]["ci95"][1] - phase2[s]["phase_diff_months"]["value"] for s in sites],
        ]
    )
    colors = ["#1f77b4" if s in CLEAN_SITES else "#9a9a9a" for s in sites]

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.1), constrained_layout=True)
    ax = axes[0]
    lo = bulk_sink_ratio(ALPHA_13C_SAUERESSIG)
    hi = bulk_sink_ratio(ALPHA_13C_CANTRELL)
    ax.axhspan(lo, hi, color="#b7e1b4", alpha=0.45, label="Bulk-sink lab range")
    ax.errorbar(lat, ratios, yerr=ratio_err, fmt="none", ecolor="#4a90c2", lw=1.4, capsize=3, zorder=1)
    ax.scatter(lat, ratios, s=46, c=colors, edgecolor="white", linewidth=0.8, zorder=2)
    label_points(ax, lat, ratios, sites)
    ax.set_xlabel("Latitude (deg)")
    ax.set_ylabel("A(d13C) / A(dD)")
    ax.set_title("(a) Raw seasonal amplitude ratio")
    ax.set_xlim(-95, 90)
    ax.set_ylim(0, max(0.16, float(np.nanmax(ratios + ratio_err[1])) * 1.05))
    ax.grid(alpha=0.25)
    ax.legend(loc="upper left", frameon=False)

    ax = axes[1]
    ax.axhline(0.0, color="#555555", ls="--", lw=1.1)
    ax.axhspan(-1.0, 1.0, color="#cccccc", alpha=0.25, label="+/- 1 month")
    ax.errorbar(lat, phase, yerr=phase_err, fmt="none", ecolor="#ff7f0e", lw=1.4, capsize=3, zorder=1)
    ax.scatter(lat, phase, s=44, c="#ff7f0e", edgecolor="white", linewidth=0.8, zorder=2)
    label_points(ax, lat, phase, sites)
    ax.set_xlabel("Latitude (deg)")
    ax.set_ylabel("Peak month: d13C minus dD")
    ax.set_title("(b) Raw isotope phase alignment")
    ax.set_xlim(-95, 90)
    ax.set_ylim(-6.5, 6.5)
    ax.grid(alpha=0.25)
    ax.legend(loc="upper left", frameon=False)

    out = OUT_DIR / "fig1_raw_ratio_phase.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def make_fig2_phase_convergence():
    phase6 = load_json(RESULTS_DIR / "phase6_phasor" / "phasor_results.json")
    summary = load_json(RESULTS_DIR / "phase1_data" / "site_summary.json")
    sites = ["ALT", "ZEP", "BRW", "CBA", "MHD", "KUM", "CGO", "SPO"]
    y = np.arange(len(sites))
    robs = np.array([phase6["sites"][s]["R_obs"] for s in sites])
    rcorr = np.array([phase6["sites"][s]["mc"]["R_corrected_median"] for s in sites])
    rcorr_ci = np.array([phase6["sites"][s]["mc"]["R_corrected_ci95"] for s in sites])
    rcorr_err = np.vstack([rcorr - rcorr_ci[:, 0], rcorr_ci[:, 1] - rcorr])
    pk_obs_13c = np.array([phase6["sites"][s]["peak_obs_13c"] for s in sites])
    pk_obs_dd = np.array([phase6["sites"][s]["peak_obs_dD"] for s in sites])
    pk_sink_13c = np.array([phase6["sites"][s]["peak_sink_13c"] for s in sites])
    pk_sink_dd = np.array([phase6["sites"][s]["peak_sink_dD"] for s in sites])
    diff_obs = np.array([shortest_phase_diff(a, b) for a, b in zip(pk_obs_13c, pk_obs_dd)])
    diff_sink = np.array([shortest_phase_diff(a, b) for a, b in zip(pk_sink_13c, pk_sink_dd)])

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(14.2, 4.4),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1.05, 1.35, 1.0]},
    )

    ax = axes[0]
    for i in range(len(sites)):
        ax.plot([robs[i], rcorr[i]], [y[i], y[i]], color="#9e9e9e", lw=1.3, zorder=1)
    ax.scatter(robs, y, marker="o", s=46, color="#bdbdbd", edgecolor="#555555",
               linewidth=0.8, label="Observed", zorder=3)
    ax.errorbar(rcorr, y, xerr=rcorr_err, fmt="none", ecolor="#1f77b4",
                elinewidth=1.2, capsize=2.5, alpha=0.45, zorder=2)
    ax.scatter(rcorr, y, marker="s", s=50, color="#1f77b4", edgecolor="white",
               linewidth=0.8, label="Corrected median", zorder=4)
    ax.axvline(bulk_sink_ratio(ALPHA_13C_SAUERESSIG), color="#2ca02c", ls="--", lw=1.2,
               label="Bulk sink, Saueressig")
    ax.axvline(bulk_sink_ratio(ALPHA_13C_CANTRELL), color="#d62728", ls="--", lw=1.2,
               label="Bulk sink, Cantrell")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{s} ({summary[s]['latitude']:+.0f} deg)" for s in sites])
    ax.set_xlabel("A(d13C) / A(dD)")
    ax.set_title("(a) Local wetland correction")
    ax.set_xlim(0.0, 0.15)
    ax.grid(axis="x", alpha=0.25)
    ax.invert_yaxis()
    ax.legend(loc="lower right", frameon=False)

    ax = axes[1]
    offset = 0.16
    for i in range(len(sites)):
        ax.plot([pk_obs_13c[i], pk_obs_dd[i]], [y[i] - offset, y[i] - offset],
                color="#d62728", lw=1.4, alpha=0.45, zorder=1)
        ax.plot([pk_sink_13c[i], pk_sink_dd[i]], [y[i] + offset, y[i] + offset],
                color="#1f77b4", lw=1.6, alpha=0.65, zorder=1)
    ax.scatter(pk_obs_13c, y - offset, marker="o", s=52, color="#d62728",
               edgecolor="white", linewidth=0.7, label="d13C before", zorder=3)
    ax.scatter(pk_obs_dd, y - offset, marker="D", s=48, color="#d62728",
               edgecolor="white", linewidth=0.7, label="dD before", zorder=3)
    ax.scatter(pk_sink_13c, y + offset, marker="o", s=52, color="#1f77b4",
               edgecolor="white", linewidth=0.7, label="d13C after", zorder=3)
    ax.scatter(pk_sink_dd, y + offset, marker="D", s=48, color="#1f77b4",
               edgecolor="white", linewidth=0.7, label="dD after", zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(sites)
    ax.set_xticks(range(12))
    ax.set_xticklabels(MONTH_NAMES)
    ax.set_xlim(-0.5, 11.5)
    ax.set_xlabel("Peak month")
    ax.set_title("(b) Peak months before and after correction")
    ax.grid(axis="x", alpha=0.25)
    ax.invert_yaxis()
    ax.legend(loc="lower right", ncol=2, frameon=False)

    ax = axes[2]
    x = np.arange(len(sites))
    width = 0.36
    ax.bar(x - width / 2, diff_obs, width, color="#d62728", alpha=0.60,
           label="Before correction")
    ax.bar(x + width / 2, diff_sink, width, color="#1f77b4", alpha=0.80,
           label="After correction")
    ax.axhline(1.0, color="#444444", ls="--", lw=1.1, label="1 month")
    ax.set_xticks(x)
    ax.set_xticklabels(sites)
    ax.set_ylabel("Abs. d13C-dD peak-month difference")
    ax.set_title("(c) Phase convergence")
    ax.set_ylim(0, max(2.2, float(np.nanmax(diff_obs)) * 1.2))
    ax.grid(axis="y", alpha=0.25)
    ax.legend(loc="upper left", frameon=False)

    out = OUT_DIR / "fig2_phase_convergence.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def make_fig3_alpha_sensitivity_uncertainty():
    phase6 = load_json(RESULTS_DIR / "phase6_phasor" / "phasor_results.json")
    phase14 = load_json(RESULTS_DIR / "phase14_sh_wetland_sensitivity" / "sh_wetland_sensitivity_results.json")
    phase13 = load_json(RESULTS_DIR / "phase13_uncertainty_attribution" / "uncertainty_attribution_results.json")

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2), constrained_layout=True)
    lab_lines = [
        (ALPHA_13C_SAUERESSIG, "Saueressig", "#2ca02c"),
        (ALPHA_13C_CANTRELL, "Cantrell", "#d62728"),
    ]

    ax = axes[0]
    items = [
        (
            "SH only",
            phase6["multi_site_result"]["sh_only"]["alpha_13c_oh_median"],
            phase6["multi_site_result"]["sh_only"]["alpha_13c_oh_ci95"],
            "#ff7f0e",
        ),
        (
            "All clean sites",
            phase6["multi_site_result"]["alpha_13c_oh_median"],
            phase6["multi_site_result"]["alpha_13c_oh_ci95"],
            "#222222",
        ),
    ]
    for value, label, color in lab_lines:
        ax.axvline(value, color=color, ls="--", lw=1.5, label=label)
    for i, (label, median, ci95, color) in enumerate(items):
        ax.errorbar(median, i, xerr=ci_error(median, ci95), fmt="D", ms=7, color=color,
                    ecolor=color, capsize=4, lw=2)
    ax.set_yticks(range(len(items)))
    ax.set_yticklabels([x[0] for x in items])
    ax.set_xlabel("alpha13C_OH")
    ax.set_title("(a) Atmospheric KIE constraint")
    ax.set_xlim(0.996, 1.027)
    ax.grid(axis="x", alpha=0.25)
    ax.legend(loc="lower right", frameon=False)

    ax = axes[1]
    image = plot_mass_conserving_heatmap(
        ax,
        phase14,
        "mean_alpha_13C_OH",
        "(b) Mass-conserving SH sensitivity",
    )
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)

    ax = axes[2]
    ranking = phase13["attribution"]["ranking"]
    pretty = {
        "observation": "Observations",
        "wetland_phasor": "Wetland phasor",
        "wetland_isotopes": "Wetland isotopes",
        "sink_fractions": "Sink fractions",
        "alpha_D_OH": "alphaD_OH",
        "bb_correction": "Biomass burning",
        "non_oh_kie": "Non-OH KIEs",
    }
    fractions = [
        phase13["attribution"][key]["fraction_of_oat_variance"]
        for key in ranking
    ]
    y = np.arange(len(ranking))
    ax.barh(y, fractions, color="#4c78a8")
    ax.set_yticks(y)
    ax.set_yticklabels([pretty.get(k, k) for k in ranking])
    ax.invert_yaxis()
    ax.set_xlabel("Fraction of diagnostic variance")
    ax.set_title("(c) Uncertainty attribution")
    ax.set_xlim(0, 0.42)
    ax.grid(axis="x", alpha=0.25)

    out = OUT_DIR / "fig3_alpha_sensitivity_uncertainty.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def make_figs5_ganesan_delta13c_sensitivity():
    phase6 = load_json(RESULTS_DIR / "phase6_phasor" / "phasor_results.json")
    labels = ["Uniform base", "NH high/mid", "Tropics", "SH/ambiguous"]
    values = [D13C_WETLAND_BASE, -67.8, -56.7, D13C_WETLAND_BASE]
    colors = ["#777777", "#1f77b4", "#2ca02c", "#777777"]

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 3.8), constrained_layout=True)
    ax = axes[0]
    y = np.arange(len(labels))
    for yy, value, color in zip(y, values, colors):
        ax.plot([D13C_WETLAND_BASE, value], [yy, yy], color=color, lw=2.6, alpha=0.85)
    ax.scatter(values, y, color=colors, s=58, edgecolor="white", linewidth=0.8, zorder=3)
    ax.axvline(D13C_WETLAND_BASE, color="#222222", lw=1.1, ls="--", label="Base value")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Wetland d13C-CH4 source signature (permil)")
    ax.set_title("(a) Source-signature substitution")
    ax.set_xlim(-72, -52)
    ax.grid(axis="x", alpha=0.25)
    ax.legend(frameon=False)

    ax = axes[1]
    sites = CLEAN_SITES
    y = np.arange(len(sites))
    base_r = np.array([phase6["sites"][site]["R_corrected"] for site in sites])
    banded_r = np.array([ganesan_banded_ratio(phase6["sites"][site]) for site in sites])
    for yy, x0, x1 in zip(y, base_r, banded_r):
        ax.plot([x0, x1], [yy, yy], color="#9a9a9a", lw=1.3, zorder=1)
    ax.scatter(base_r, y, marker="o", s=42, color="#777777", label="Uniform base", zorder=3)
    ax.scatter(banded_r, y, marker="s", s=46, color="#1f77b4", label="Banded d13C", zorder=4)
    ax.axvspan(
        bulk_sink_ratio(ALPHA_13C_SAUERESSIG),
        bulk_sink_ratio(ALPHA_13C_CANTRELL),
        color="#b7e1b4",
        alpha=0.42,
        label="Bulk-sink lab range",
        zorder=0,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(sites)
    ax.invert_yaxis()
    ax.set_xlabel("Wetland-corrected R")
    ax.set_title("(b) Resulting corrected-ratio shift")
    ax.set_xlim(0.02, max(0.09, float(np.nanmax(banded_r)) * 1.08))
    ax.grid(axis="x", alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    out = OUT_DIR / "figS5_ganesan_delta13c_sensitivity.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def make_figs7_sh_source_region_sensitivity():
    phase14 = load_json(RESULTS_DIR / "phase14_sh_wetland_sensitivity" / "sh_wetland_sensitivity_results.json")
    phase6 = load_json(RESULTS_DIR / "phase6_phasor" / "phasor_results.json")
    sh_mc_alpha = phase6["multi_site_result"]["sh_only"]["alpha_13c_oh_median"]

    fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.2), constrained_layout=True)
    image = plot_mass_conserving_heatmap(
        axes[0],
        phase14,
        "mean_alpha_13C_OH",
        "(a) Deterministic mass-conserving grid",
        mark_lab_range=True,
    )
    fig.colorbar(image, ax=axes[0], fraction=0.046, pad=0.04)
    axes[0].text(
        0.0,
        -0.24,
        "* within Saueressig-Cantrell lab range",
        transform=axes[0].transAxes,
        fontsize=8,
        ha="left",
        va="top",
    )

    plot_mass_conserving_slice(
        axes[1],
        phase14,
        "tropics_only",
        "(b) Tropics only",
        "Tropics response weight",
        "#2ca02c",
        reference_alpha=sh_mc_alpha,
    )
    plot_mass_conserving_slice(
        axes[2],
        phase14,
        "nh_high_only",
        "(c) Delayed NH high only",
        "Delayed NH high response weight",
        "#1f77b4",
        reference_alpha=sh_mc_alpha,
    )
    out = OUT_DIR / "figS7_sh_source_region_sensitivity.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def make_figs10_uncertainty_attribution():
    phase13 = load_json(RESULTS_DIR / "phase13_uncertainty_attribution" / "uncertainty_attribution_results.json")
    ranking = phase13["attribution"]["ranking"]
    pretty = {
        "observation": "Observations",
        "wetland_phasor": "Wetland phasor",
        "wetland_isotopes": "Wetland isotopes",
        "sink_fractions": "Sink fractions",
        "alpha_D_OH": "alphaD_OH",
        "bb_correction": "Biomass burning",
        "non_oh_kie": "Non-OH KIEs",
    }
    fractions = [phase13["attribution"][key]["fraction_of_oat_variance"] for key in ranking]

    fig, ax = plt.subplots(figsize=(6.2, 3.7), constrained_layout=True)
    y = np.arange(len(ranking))
    ax.barh(y, fractions, color="#4c78a8")
    ax.set_yticks(y)
    ax.set_yticklabels([pretty.get(k, k) for k in ranking])
    ax.invert_yaxis()
    ax.set_xlabel("Fraction of diagnostic one-at-a-time variance")
    ax.set_title("Grouped uncertainty attribution")
    ax.set_xlim(0, 0.42)
    ax.grid(axis="x", alpha=0.25)
    out = OUT_DIR / "figS10_uncertainty_attribution.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def main():
    set_common_style()
    outputs = [
        make_fig1_raw_ratio_phase(),
        make_fig2_phase_convergence(),
        make_fig3_alpha_sensitivity_uncertainty(),
        make_figs5_ganesan_delta13c_sensitivity(),
        make_figs7_sh_source_region_sensitivity(),
        make_figs10_uncertainty_attribution(),
    ]
    for path in outputs:
        print(path.relative_to(EXPT_DIR))


if __name__ == "__main__":
    main()
