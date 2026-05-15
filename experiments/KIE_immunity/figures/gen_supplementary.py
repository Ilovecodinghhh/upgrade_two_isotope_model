#!/usr/bin/env python3
"""Generate all 8 supplementary figures (S1–S8) for the KIE immunity manuscript."""

import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

RESULTS = Path(__file__).parent.parent / 'results'
OUTDIR = Path(__file__).parent

def load(name):
    with open(RESULTS / name) as f:
        return json.load(f)

# ── Color palette ──
C_KIE = '#e74c3c'
C_SIG = '#3498db'
C_TAU = '#2ecc71'
C_RES = '#95a5a6'
C_MAIN = '#2c3e50'

plt.rcParams.update({
    'font.size': 10, 'axes.titlesize': 12, 'axes.labelsize': 11,
    'figure.dpi': 200, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.15
})


# ═══════════════════════════════════════════
# Fig S1: Robustness Matrix (Phase 8)
# ═══════════════════════════════════════════
def fig_s1():
    d = {k: v for k, v in load('phase8_robustness_matrix.json').items() if not k.startswith('_')}
    taus = sorted(set(v['tau'] for v in d.values()))
    oh_ds = sorted(set(v['oh_d'] for v in d.values()))
    cls = sorted(set(v['cl_frac'] for v in d.values()))

    fig, axes = plt.subplots(1, len(taus), figsize=(5*len(taus), 4.5), sharey=True)
    if len(taus) == 1:
        axes = [axes]

    for ax, tau in zip(axes, taus):
        mat = np.full((len(oh_ds), len(cls)), np.nan)
        for i, oh in enumerate(oh_ds):
            for j, cl in enumerate(cls):
                key = f'τ={tau}_OHD={oh}_Cl={cl}'
                if key in d:
                    mat[i, j] = d[key]['trend_median']

        vmax = max(abs(np.nanmin(mat)), abs(np.nanmax(mat)))
        im = ax.imshow(mat, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto', origin='lower')
        ax.set_xticks(range(len(cls)))
        ax.set_xticklabels([f'{c*100:.1f}%' for c in cls], rotation=45, fontsize=8)
        ax.set_xlabel('Cl fraction')
        if ax == axes[0]:
            ax.set_yticks(range(len(oh_ds)))
            ax.set_yticklabels([f'{o:.3f}' for o in oh_ds])
            ax.set_ylabel('OH-D KIE (α)')
        ax.set_title(f'τ = {tau} yr')

        for i in range(len(oh_ds)):
            for j in range(len(cls)):
                if not np.isnan(mat[i, j]):
                    ax.text(j, i, f'{mat[i,j]:+.1f}', ha='center', va='center', fontsize=7,
                           color='white' if abs(mat[i,j]) > vmax*0.6 else 'black')

    fig.colorbar(im, ax=axes, label='ΔFF trend (Tg/yr)', shrink=0.8)
    fig.suptitle('Fig. S1: Robustness Matrix — ΔFF Across Parameter Space', fontsize=13, y=1.02)
    fig.savefig(OUTDIR / 'fig_S1_robustness_matrix.png')
    fig.savefig(OUTDIR / 'fig_S1_robustness_matrix.pdf')
    plt.close(fig)
    print('  S1 done')


# ═══════════════════════════════════════════
# Fig S2: Lifetime Sensitivity (Phase 5)
# ═══════════════════════════════════════════
def fig_s2():
    d = load('phase5_tau_sensitivity.json')
    labels = list(d.keys())
    medians = [d[k]['trend_median'] for k in labels]
    lo = [d[k]['trend_5pct'] for k in labels]
    hi = [d[k]['trend_95pct'] for k in labels]
    sigmas = [d[k]['sigma_ff'] for k in labels]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

    x = range(len(labels))
    errs = [[m - l for m, l in zip(medians, lo)], [h - m for h, m in zip(hi, medians)]]
    ax1.bar(x, medians, color=C_MAIN, alpha=0.7, edgecolor='white')
    ax1.errorbar(x, medians, yerr=errs, fmt='none', ecolor='black', capsize=4)
    ax1.axhline(0, color='grey', ls='--', lw=0.8)
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=30, ha='right', fontsize=8)
    ax1.set_ylabel('ΔFF trend (Tg/yr)')
    ax1.set_title('(a) Post-2007 FF Trend')

    ax2.bar(x, sigmas, color=C_SIG, alpha=0.7, edgecolor='white')
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, rotation=30, ha='right', fontsize=8)
    ax2.set_ylabel('σ(FF) (Tg/yr)')
    ax2.set_title('(b) FF Uncertainty')

    fig.suptitle('Fig. S2: Lifetime Sensitivity', fontsize=13)
    fig.tight_layout()
    fig.savefig(OUTDIR / 'fig_S2_tau_sensitivity.png')
    fig.savefig(OUTDIR / 'fig_S2_tau_sensitivity.pdf')
    plt.close(fig)
    print('  S2 done')


# ═══════════════════════════════════════════
# Fig S3: OH-D KIE Sensitivity (Phase 6)
# ═══════════════════════════════════════════
def fig_s3():
    d = load('phase6_OHD_sensitivity.json')
    labels = list(d.keys())
    medians = [d[k]['trend_median'] for k in labels]
    lo = [d[k]['trend_5pct'] for k in labels]
    hi = [d[k]['trend_95pct'] for k in labels]
    sigmas = [d[k]['sigma_ff'] for k in labels]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

    x = range(len(labels))
    errs = [[m - l for m, l in zip(medians, lo)], [h - m for h, m in zip(hi, medians)]]
    ax1.bar(x, medians, color=C_KIE, alpha=0.7, edgecolor='white')
    ax1.errorbar(x, medians, yerr=errs, fmt='none', ecolor='black', capsize=4)
    ax1.axhline(0, color='grey', ls='--', lw=0.8)
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=30, ha='right', fontsize=8)
    ax1.set_ylabel('ΔFF trend (Tg/yr)')
    ax1.set_title('(a) Post-2007 FF Trend')

    ax2.bar(x, sigmas, color=C_SIG, alpha=0.7, edgecolor='white')
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, rotation=30, ha='right', fontsize=8)
    ax2.set_ylabel('σ(FF) (Tg/yr)')
    ax2.set_title('(b) FF Uncertainty')

    fig.suptitle('Fig. S3: OH-D KIE Sensitivity', fontsize=13)
    fig.tight_layout()
    fig.savefig(OUTDIR / 'fig_S3_OHD_sensitivity.png')
    fig.savefig(OUTDIR / 'fig_S3_OHD_sensitivity.pdf')
    plt.close(fig)
    print('  S3 done')


# ═══════════════════════════════════════════
# Fig S4: Cl Fraction Sensitivity (Phase 7)
# ═══════════════════════════════════════════
def fig_s4():
    d = load('phase7_Cl_sensitivity.json')
    labels = list(d.keys())
    medians = [d[k]['trend_median'] for k in labels]
    lo = [d[k]['trend_5pct'] for k in labels]
    hi = [d[k]['trend_95pct'] for k in labels]
    kie_pcts = [d[k]['kie_pct'] for k in labels]
    sig_pcts = [d[k]['sig_pct'] for k in labels]

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    x = range(len(labels))
    errs = [[m - l for m, l in zip(medians, lo)], [h - m for h, m in zip(hi, medians)]]
    axes[0].bar(x, medians, color=C_MAIN, alpha=0.7, edgecolor='white')
    axes[0].errorbar(x, medians, yerr=errs, fmt='none', ecolor='black', capsize=4)
    axes[0].axhline(0, color='grey', ls='--', lw=0.8)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=35, ha='right', fontsize=7)
    axes[0].set_ylabel('ΔFF trend (Tg/yr)')
    axes[0].set_title('(a) Post-2007 FF Trend')

    axes[1].bar(x, kie_pcts, color=C_KIE, alpha=0.7, edgecolor='white')
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=35, ha='right', fontsize=7)
    axes[1].set_ylabel('KIE variance share (%)')
    axes[1].set_title('(b) KIE% vs Cl Fraction')

    axes[2].bar(x, sig_pcts, color=C_SIG, alpha=0.7, edgecolor='white')
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(labels, rotation=35, ha='right', fontsize=7)
    axes[2].set_ylabel('Source signature variance share (%)')
    axes[2].set_title('(c) Sig% vs Cl Fraction')

    fig.suptitle('Fig. S4: Cl Fraction Sensitivity with Variance Decomposition', fontsize=13)
    fig.tight_layout()
    fig.savefig(OUTDIR / 'fig_S4_Cl_sensitivity.png')
    fig.savefig(OUTDIR / 'fig_S4_Cl_sensitivity.pdf')
    plt.close(fig)
    print('  S4 done')


# ═══════════════════════════════════════════
# Fig S5: Interhemispheric Exchange (Phase 11)
# ═══════════════════════════════════════════
def fig_s5():
    d = load('phase11_tau_ex.json')
    labels = list(d.keys())
    medians = [d[k]['trend_median'] for k in labels]
    lo = [d[k]['trend_5pct'] for k in labels]
    hi = [d[k]['trend_95pct'] for k in labels]
    sigmas = [d[k]['sigma_ff'] for k in labels]
    kie_pcts = [d[k]['kie_pct'] for k in labels]

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    x = range(len(labels))
    errs = [[m - l for m, l in zip(medians, lo)], [h - m for h, m in zip(hi, medians)]]
    axes[0].bar(x, medians, color=C_MAIN, alpha=0.7, edgecolor='white')
    axes[0].errorbar(x, medians, yerr=errs, fmt='none', ecolor='black', capsize=4)
    axes[0].axhline(0, color='grey', ls='--', lw=0.8)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=30, ha='right', fontsize=8)
    axes[0].set_ylabel('ΔFF trend (Tg/yr)')
    axes[0].set_title('(a) Post-2007 FF Trend')

    axes[1].bar(x, sigmas, color=C_SIG, alpha=0.7, edgecolor='white')
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=30, ha='right', fontsize=8)
    axes[1].set_ylabel('σ(FF) (Tg/yr)')
    axes[1].set_title('(b) FF Uncertainty')

    axes[2].bar(x, kie_pcts, color=C_KIE, alpha=0.7, edgecolor='white')
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(labels, rotation=30, ha='right', fontsize=8)
    axes[2].set_ylabel('KIE variance share (%)')
    axes[2].set_title('(c) KIE%')

    fig.suptitle('Fig. S5: Interhemispheric Exchange Time Sensitivity', fontsize=13)
    fig.tight_layout()
    fig.savefig(OUTDIR / 'fig_S5_tau_ex_sensitivity.png')
    fig.savefig(OUTDIR / 'fig_S5_tau_ex_sensitivity.pdf')
    plt.close(fig)
    print('  S5 done')


# ═══════════════════════════════════════════
# Fig S6: Bootstrap CIs on Variance Decomposition (Phase 9)
# ═══════════════════════════════════════════
def fig_s6():
    d = load('phase9_bootstrap.json')
    configs = list(d.keys())
    components = ['kie_pct', 'sig_pct', 'tau_pct', 'resid_pct']
    labels_comp = ['KIE', 'Source Sig.', 'Lifetime', 'Residual']
    colors = [C_KIE, C_SIG, C_TAU, C_RES]

    fig, axes = plt.subplots(1, len(configs), figsize=(5*len(configs), 5), sharey=True)
    if len(configs) == 1:
        axes = [axes]

    for ax, cfg in zip(axes, configs):
        vals = d[cfg]
        meds = [vals[c][0] for c in components]
        los = [vals[c][1] for c in components]
        his = [vals[c][2] for c in components]

        x = range(len(components))
        errs = [[m - l for m, l in zip(meds, los)], [h - m for h, m in zip(his, meds)]]
        bars = ax.bar(x, meds, color=colors, alpha=0.8, edgecolor='white')
        ax.errorbar(x, meds, yerr=errs, fmt='none', ecolor='black', capsize=5, lw=1.5)
        ax.set_xticks(x)
        ax.set_xticklabels(labels_comp, rotation=30, ha='right')
        ax.set_title(cfg.replace('_', ' '), fontsize=10)
        if ax == axes[0]:
            ax.set_ylabel('Variance share (%)')

        for i, (m, lo, hi) in enumerate(zip(meds, los, his)):
            ax.text(i, hi + 1.5, f'{m:.1f}', ha='center', fontsize=7, fontweight='bold')

    fig.suptitle('Fig. S6: Bootstrap 95% CIs on Variance Decomposition', fontsize=13)
    fig.tight_layout()
    fig.savefig(OUTDIR / 'fig_S6_bootstrap_variance.png')
    fig.savefig(OUTDIR / 'fig_S6_bootstrap_variance.pdf')
    plt.close(fig)
    print('  S6 done')


# ═══════════════════════════════════════════
# Fig S7: W Matrix Sensitivity (Phase 14)
# ═══════════════════════════════════════════
def fig_s7():
    d = load('phase14_W_sensitivity.json')
    labels = list(d.keys())
    medians = [d[k]['trend_step_median'] for k in labels]
    lo = [d[k]['trend_step_5pct'] for k in labels]
    hi = [d[k]['trend_step_95pct'] for k in labels]
    kie_pcts = [d[k]['kie_pct'] for k in labels]
    sig_pcts = [d[k]['sig_pct'] for k in labels]

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))

    x = range(len(labels))
    errs = [[m - l for m, l in zip(medians, lo)], [h - m for h, m in zip(hi, medians)]]
    axes[0].bar(x, medians, color=C_MAIN, alpha=0.7, edgecolor='white')
    axes[0].errorbar(x, medians, yerr=errs, fmt='none', ecolor='black', capsize=4)
    axes[0].axhline(0, color='grey', ls='--', lw=0.8)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=35, ha='right', fontsize=7)
    axes[0].set_ylabel('ΔFF trend (Tg/yr)')
    axes[0].set_title('(a) Post-2007 FF Trend')

    axes[1].bar(x, kie_pcts, color=C_KIE, alpha=0.7, edgecolor='white')
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=35, ha='right', fontsize=7)
    axes[1].set_ylabel('KIE variance share (%)')
    axes[1].set_title('(b) KIE%')

    axes[2].bar(x, sig_pcts, color=C_SIG, alpha=0.7, edgecolor='white')
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(labels, rotation=35, ha='right', fontsize=7)
    axes[2].set_ylabel('Source sig. variance share (%)')
    axes[2].set_title('(c) Sig%')

    fig.suptitle('Fig. S7: W Matrix Sensitivity — ΔFF and Variance Components', fontsize=13)
    fig.tight_layout()
    fig.savefig(OUTDIR / 'fig_S7_W_sensitivity.png')
    fig.savefig(OUTDIR / 'fig_S7_W_sensitivity.pdf')
    plt.close(fig)
    print('  S7 done')


# ═══════════════════════════════════════════
# Fig S8: MC Convergence (Phase 16)
# ═══════════════════════════════════════════
def fig_s8():
    d = load('phase16_convergence.json')
    ns = sorted(d.keys(), key=lambda x: int(x))
    n_vals = [int(n) for n in ns]
    sigmas = [d[n]['sigma_ff'] for n in ns]
    medians = [d[n]['trend_median'] for n in ns]
    kie_pcts = [d[n]['kie_pct'] for n in ns]

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    for ax, vals, label, color, ylabel in zip(
        axes,
        [sigmas, medians, kie_pcts],
        ['σ(FF)', 'ΔFF trend', 'KIE%'],
        [C_SIG, C_MAIN, C_KIE],
        ['σ(FF) (Tg/yr)', 'ΔFF (Tg/yr)', 'KIE variance share (%)']
    ):
        ax.plot(n_vals, vals, 'o-', color=color, lw=2, markersize=6)
        ax.set_xlabel('N iterations')
        ax.set_ylabel(ylabel)
        ax.set_title(f'({"abc"[axes.tolist().index(ax)]}) {label}')
        ax.grid(True, alpha=0.3)
        # shade ±5% band around final value
        final = vals[-1]
        ax.axhspan(final*0.95, final*1.05, alpha=0.1, color=color)
        ax.axhline(final, color=color, ls='--', lw=0.8, alpha=0.5)

    fig.suptitle('Fig. S8: Monte Carlo Convergence', fontsize=13)
    fig.tight_layout()
    fig.savefig(OUTDIR / 'fig_S8_convergence.png')
    fig.savefig(OUTDIR / 'fig_S8_convergence.pdf')
    plt.close(fig)
    print('  S8 done')


if __name__ == '__main__':
    print('Generating supplementary figures...')
    fig_s1()
    fig_s2()
    fig_s3()
    fig_s4()
    fig_s5()
    fig_s6()
    fig_s7()
    fig_s8()
    print('All done!')
