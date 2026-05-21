#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig11_phase_convergence.py — Visualize phase convergence after phasor correction
================================================================================

Shows that after removing wetland source contamination, the peak months of
δ¹³C and δD seasonal cycles converge — both isotopes now reflect the same
OH-driven sink signal.

Three panels:
  (a) Before/after peak months for each site (paired dot plot)
  (b) Phase difference Δ(δ¹³C − δD) before vs after (bar chart)
  (c) Phasor clock diagrams for 4 representative sites
"""

from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import matplotlib.patches as mpatches

EXPT_DIR = Path(__file__).resolve().parent.parent
FIG_DIR  = EXPT_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

MONTH_NAMES = ["Jan","Feb","Mar","Apr","May","Jun",
               "Jul","Aug","Sep","Oct","Nov","Dec"]

SITE_ORDER = ["ALT","ZEP","BRW","CBA","MHD","KUM","CGO","SPO"]
SITE_LATS  = {"ALT":82.5, "ZEP":78.9, "BRW":71.3, "CBA":55.2,
              "MHD":53.3, "KUM":19.6, "CGO":-40.7, "SPO":-90.0}

def shortest_phase_diff(a, b):
    d = abs(a - b)
    return min(d, 12 - d)


def main():
    with open(EXPT_DIR / "results" / "phase6_phasor" / "phasor_results.json") as f:
        data = json.load(f)

    # Extract data
    sites = []
    for code in SITE_ORDER:
        s = data["sites"][code]
        d = {
            "code": code,
            "lat": SITE_LATS[code],
            "pk_obs_13c": s["peak_obs_13c"],
            "pk_obs_dD":  s["peak_obs_dD"],
            "pk_sink_13c": s["peak_sink_13c"],
            "pk_sink_dD":  s["peak_sink_dD"],
        }
        d["diff_obs"]  = shortest_phase_diff(d["pk_obs_13c"], d["pk_obs_dD"])
        d["diff_sink"] = shortest_phase_diff(d["pk_sink_13c"], d["pk_sink_dD"])
        sites.append(d)

    codes      = [s["code"] for s in sites]
    diff_obs   = [s["diff_obs"] for s in sites]
    diff_sink  = [s["diff_sink"] for s in sites]

    # =====================================================================
    # Figure
    # =====================================================================
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.30,
                          left=0.07, right=0.97, top=0.92, bottom=0.07)

    # ----- Panel (a): Peak month dot plot -----
    ax_a = fig.add_subplot(gs[0, 0])

    y = np.arange(len(codes))
    offset = 0.15  # vertical offset for before/after

    # Before correction
    ax_a.scatter([s["pk_obs_13c"] for s in sites], y - offset, 
                 marker="o", s=70, color="#d62728", edgecolors="k", linewidths=0.5,
                 zorder=5, label="δ¹³C (before)")
    ax_a.scatter([s["pk_obs_dD"] for s in sites], y - offset,
                 marker="D", s=70, color="#d62728", edgecolors="k", linewidths=0.5,
                 zorder=5, label="δD (before)")
    # Connect before pairs with gray lines
    for i, s in enumerate(sites):
        ax_a.plot([s["pk_obs_13c"], s["pk_obs_dD"]], [i - offset, i - offset],
                  "-", color="#d62728", lw=1.5, alpha=0.4)

    # After correction
    ax_a.scatter([s["pk_sink_13c"] for s in sites], y + offset,
                 marker="o", s=70, color="#1f77b4", edgecolors="k", linewidths=0.5,
                 zorder=5, label="δ¹³C (after)")
    ax_a.scatter([s["pk_sink_dD"] for s in sites], y + offset,
                 marker="D", s=70, color="#1f77b4", edgecolors="k", linewidths=0.5,
                 zorder=5, label="δD (after)")
    # Connect after pairs with blue lines
    for i, s in enumerate(sites):
        ax_a.plot([s["pk_sink_13c"], s["pk_sink_dD"]], [i + offset, i + offset],
                  "-", color="#1f77b4", lw=1.5, alpha=0.6)

    ax_a.set_yticks(y)
    ax_a.set_yticklabels([f"{s['code']} ({s['lat']:+.0f}°)" for s in sites], fontsize=9)
    ax_a.set_xticks(range(12))
    ax_a.set_xticklabels(MONTH_NAMES, fontsize=8)
    ax_a.set_xlim(-0.5, 11.5)
    ax_a.set_xlabel("Peak month", fontsize=10)
    ax_a.set_title("(a) Peak months: δ¹³C vs δD", fontsize=11, fontweight="bold")
    ax_a.legend(fontsize=7.5, loc="lower left", ncol=2, 
                framealpha=0.9, columnspacing=1.0)
    ax_a.grid(axis="x", alpha=0.2)
    # Shade NH summer and SH summer
    ax_a.axvspan(4.5, 7.5, alpha=0.04, color="orange", label="_NH summer")
    ax_a.axvspan(-0.5, 2.5, alpha=0.04, color="blue", label="_SH summer")
    ax_a.invert_yaxis()

    # ----- Panel (b): Phase difference bar chart -----
    ax_b = fig.add_subplot(gs[0, 1])

    x = np.arange(len(codes))
    w = 0.35

    bars_obs  = ax_b.bar(x - w/2, diff_obs,  w, color="#d62728", alpha=0.6,
                          edgecolor="k", linewidth=0.5, label="Before correction")
    bars_sink = ax_b.bar(x + w/2, diff_sink, w, color="#1f77b4", alpha=0.7,
                          edgecolor="k", linewidth=0.5, label="After correction")

    # Percentage improvement labels
    for i in range(len(codes)):
        if diff_obs[i] > 0.05:
            pct = (diff_obs[i] - diff_sink[i]) / diff_obs[i] * 100
            ax_b.annotate(f"{pct:+.0f}%",
                          xy=(x[i] + w/2, diff_sink[i]),
                          xytext=(0, 5), textcoords="offset points",
                          ha="center", fontsize=7.5, color="#1f77b4", fontweight="bold")

    ax_b.set_xticks(x)
    ax_b.set_xticklabels([f"{s['code']}" for s in sites], fontsize=9)
    ax_b.set_ylabel("Phase difference |Δ(δ¹³C − δD)| (months)", fontsize=10)
    ax_b.set_title("(b) Phase convergence after phasor correction", fontsize=11, fontweight="bold")
    ax_b.legend(fontsize=9, loc="upper left")
    ax_b.set_ylim(0, 2.2)
    ax_b.axhline(0, color="k", lw=0.5)
    ax_b.grid(axis="y", alpha=0.2)

    # Add mean lines
    mean_obs  = np.mean(diff_obs)
    mean_sink = np.mean(diff_sink)
    ax_b.axhline(mean_obs,  color="#d62728", ls="--", lw=1.2, alpha=0.5)
    ax_b.axhline(mean_sink, color="#1f77b4", ls="--", lw=1.2, alpha=0.5)
    ax_b.text(len(codes)-0.5, mean_obs + 0.05,  f"mean = {mean_obs:.2f} mo",
              fontsize=8, color="#d62728", ha="right")
    ax_b.text(len(codes)-0.5, mean_sink + 0.05, f"mean = {mean_sink:.2f} mo",
              fontsize=8, color="#1f77b4", ha="right")

    # ----- Panels (c–f): Clock diagrams for 4 sites -----
    gs_bottom = gs[1, :].subgridspec(1, 4, wspace=0.4)
    show_sites = ["BRW", "CBA", "CGO", "SPO"]
    for k, code in enumerate(show_sites):
        ax_c = fig.add_subplot(gs_bottom[0, k], projection="polar")
        s = [si for si in sites if si["code"] == code][0]

        # Convert month to angle (0=Jan at top, clockwise)
        # polar: 0 = right, counterclockwise. We want 0=top, clockwise.
        # angle = π/2 - 2π * month/12
        def month_to_angle(m):
            return np.pi/2 - 2*np.pi * m / 12

        # Clock face setup
        ax_c.set_theta_zero_location("N")
        ax_c.set_theta_direction(-1)  # clockwise
        ax_c.set_xticks(np.linspace(0, 2*np.pi, 12, endpoint=False))
        ax_c.set_xticklabels(MONTH_NAMES, fontsize=7)
        ax_c.set_ylim(0, 1.15)
        ax_c.set_yticks([])
        ax_c.grid(True, alpha=0.15)

        # Helper: draw arrow from center to month position
        def draw_hand(month, color, ls, lw, label, r=1.0):
            theta = 2 * np.pi * month / 12  # clockwise from top
            ax_c.annotate("", xy=(theta, r), xytext=(0, 0),
                          arrowprops=dict(arrowstyle="-|>", color=color,
                                          lw=lw, linestyle=ls, mutation_scale=14))
            ax_c.plot(theta, r, "o", color=color, ms=5, zorder=10)

        # Before correction
        draw_hand(s["pk_obs_13c"], "#d62728", "-",  2.0, "δ¹³C obs", r=0.85)
        draw_hand(s["pk_obs_dD"],  "#d62728", "--", 2.0, "δD obs",   r=0.85)

        # After correction
        draw_hand(s["pk_sink_13c"], "#1f77b4", "-",  2.5, "δ¹³C sink", r=1.05)
        draw_hand(s["pk_sink_dD"],  "#1f77b4", "--", 2.5, "δD sink",   r=1.05)

        # Title with phase diff improvement
        ax_c.set_title(f"{code} ({s['lat']:+.0f}°)\n"
                        f"Δ: {s['diff_obs']:.2f} → {s['diff_sink']:.2f} mo",
                        fontsize=9.5, fontweight="bold", pad=15)

    # Add a shared legend for the clock panels
    legend_elements = [
        mpatches.Patch(facecolor="#d62728", alpha=0.5, label="Before correction"),
        mpatches.Patch(facecolor="#1f77b4", alpha=0.7, label="After correction"),
        plt.Line2D([0],[0], color="gray", ls="-",  lw=1.5, label="δ¹³C (solid)"),
        plt.Line2D([0],[0], color="gray", ls="--", lw=1.5, label="δD (dashed)"),
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=4, fontsize=9,
               bbox_to_anchor=(0.5, -0.01), framealpha=0.9)

    fig.suptitle("Fig 11: Phase convergence — wetland source correction aligns δ¹³C and δD peak months",
                 fontsize=13, fontweight="bold")

    out = FIG_DIR / "fig11_phase_convergence.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"✓ Saved: {out}")


if __name__ == "__main__":
    main()
