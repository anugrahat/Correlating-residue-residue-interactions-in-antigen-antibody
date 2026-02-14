#!/usr/bin/env python3
"""
Publication figures for elastic net regression story (corrected regression only).
  Fig 1: Elastic net regression — fit, top contacts, residue ranking
  Fig 2: MM-GBSA validation — only mutations identified by corrected regression
  Fig 3: NetFavorability vs ΔΔG — does regression signal predict effect size?
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from sklearn.linear_model import ElasticNetCV
from sklearn.metrics import r2_score, mean_squared_error
import os

# ==============================================================================
# Config
# ==============================================================================
BASE_DIR = "/home/anugraha/c1_WT/analysis"
FIG_DIR = "/home/anugraha/antibody_optimization/figures"
FREQ_FILE = os.path.join(BASE_DIR, "average_frequency.csv")
ENERGY_FILE = os.path.join(BASE_DIR, "interaction_energy.csv")
GRO_PATH = "/home/anugraha/c1_WT/pull/md10_protein.gro"
COEFF_FILE = os.path.join(BASE_DIR, "elastic_net_coefficients_3000.csv")
NF_FILE = os.path.join(BASE_DIR, "net_favorability_3000.csv")

FRAME_END = 1500
REF_FRAME = 0
WT_DG = -25.31

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIE": "H", "HIS": "H",
    "HID": "H", "HIP": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T",
    "TRP": "W", "TYR": "Y", "VAL": "V",
}

# MM-GBSA data: ONLY mutations present in corrected regression
# (mutation, dG, SEM, NF, chemistry)
MMGBSA_DATA = [
    ('F287W',  -58.99, 1.44, 17.7, 'Aromatic extension'),
    ('S247G',  -56.87, 1.20,  3.4, 'Unfavorable removal'),
    ('D245N',  -43.47, 1.10, 23.2, 'Charge removal'),
    ('T330G',  -41.47, 0.65,  3.7, 'Unfavorable removal'),
    ('N248Y',  -36.20, 1.01,  6.2, 'Aromatic addition'),
    ('Y227W',  -35.99, 0.58, 67.0, 'Aromatic extension'),
    ('A246Y',  -31.41, 0.81,  7.8, 'Aromatic addition'),
    ('Y404W',  -28.00, 0.96, 23.4, 'Aromatic extension'),
    ('WT',     -25.31, 0.86,  0.0, 'Baseline'),
    ('Y404H',  -22.79, 0.77, 23.4, 'Non-conservative'),
    ('L306G',  -16.35, 0.60,  5.8, 'Unfavorable removal'),
    ('Y157A',  -13.73, 1.49,  0.0, 'Negative control'),
]

CHEM_COLORS = {
    'Charge removal':     '#1f77b4',
    'Aromatic extension': '#2ca02c',
    'Aromatic addition':  '#9467bd',
    'Unfavorable removal':'#ff7f0e',
    'Baseline':           '#555555',
    'Non-conservative':   '#d62728',
    'Negative control':   '#d62728',
}


def load_gro_resnames(gro_path):
    resnames = {}
    with open(gro_path) as f:
        lines = f.readlines()
    for line in lines[2:]:
        if len(line) < 20:
            continue
        try:
            resid = int(line[:5].strip())
            resname = line[5:10].strip()
            if resid not in resnames:
                resnames[resid] = resname
        except (ValueError, IndexError):
            continue
    return resnames


def gro_to_global(gro_resid):
    if 1 <= gro_resid <= 194:
        return gro_resid
    elif 195 <= gro_resid <= 314:
        return gro_resid - 194 + 302
    elif 315 <= gro_resid <= 421:
        return gro_resid - 314 + 195
    raise ValueError(f"GRO resid {gro_resid} out of range")


def get_label(gro_resid, resnames):
    glob = gro_to_global(gro_resid)
    rn = resnames.get(gro_resid, "UNK")
    aa1 = AA3_TO_1.get(rn, "X")
    return f"{aa1}{glob}"


# ==============================================================================
# FIGURE 1: Elastic Net Regression (3 panels)
# ==============================================================================
def make_figure1():
    print("=" * 70)
    print("FIGURE 1: Elastic Net Regression (corrected indexing)")
    print("=" * 70)

    resnames = load_gro_resnames(GRO_PATH)

    # Load and prepare data
    df_energy = pd.read_csv(ENERGY_FILE)
    df_frequency = pd.read_csv(FREQ_FILE)
    freq_col = 'MeanFrequency' if 'MeanFrequency' in df_frequency.columns else 'InteractionFrequency'
    df_pivot = df_frequency.pivot_table(
        index='Frame', columns='ResiduePair', values=freq_col, fill_value=0
    ).sort_index()
    df_energy_idx = df_energy.set_index('Frame').sort_index()
    common = df_pivot.index.intersection(df_energy_idx.index)
    df_E = df_energy_idx.loc[common, 'AverageInteractionEnergy']
    df_F = df_pivot.loc[common]
    E_ref = df_E.loc[REF_FRAME]
    F_ref = df_F.loc[REF_FRAME]
    dE = (df_E - E_ref).loc[:FRAME_END].drop(index=REF_FRAME)
    dF = df_F.subtract(F_ref, axis='columns').loc[:FRAME_END].drop(index=REF_FRAME)
    X = dF.values
    y = dE.values
    pair_names = dF.columns.tolist()

    # Fit regression
    print("Fitting ElasticNetCV...")
    model = ElasticNetCV(
        alphas=np.logspace(-2, 2, 10),
        l1_ratio=[0.3, 0.5, 0.7, 0.9],
        cv=5, fit_intercept=False, max_iter=200000,
        random_state=42, n_jobs=-1
    )
    model.fit(X, y)
    y_pred = model.predict(X)
    r2 = r2_score(y, y_pred)
    rmse = np.sqrt(mean_squared_error(y, y_pred))
    n_nz = (model.coef_ != 0).sum()

    # Load pre-computed net favorability
    df_nf = pd.read_csv(NF_FILE)

    # ---- Create figure: 3 panels ----
    fig = plt.figure(figsize=(18, 6))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1.2, 1.3], wspace=0.35)

    # --- Panel A: Predicted vs Actual ---
    ax1 = fig.add_subplot(gs[0])
    ax1.scatter(y, y_pred, s=3, alpha=0.3, c='steelblue', edgecolors='none', rasterized=True)
    lims = [min(y.min(), y_pred.min()) - 10, max(y.max(), y_pred.max()) + 10]
    ax1.plot(lims, lims, 'k--', lw=1, alpha=0.5)
    ax1.set_xlim(lims)
    ax1.set_ylim(lims)
    ax1.set_xlabel(r'Actual $\Delta E$ (kJ/mol)', fontsize=11)
    ax1.set_ylabel(r'Predicted $\Delta E$ (kJ/mol)', fontsize=11)
    ax1.set_title('A) Elastic Net Regression Fit', fontsize=12, fontweight='bold')
    stats_text = (
        f'$R^2$ = {r2:.4f}\n'
        f'RMSE = {rmse:.1f} kJ/mol\n'
        f'{n_nz}/{len(pair_names)} non-zero $\\beta$\n'
        f'$\\alpha$ = {model.alpha_:.2f}, $\\ell_1$ = {model.l1_ratio_:.2f}'
    )
    ax1.text(0.05, 0.95, stats_text, transform=ax1.transAxes, fontsize=9,
             verticalalignment='top', bbox=dict(boxstyle='round,pad=0.4', facecolor='wheat', alpha=0.8))
    eq_text = r'$\Delta E(t) = \sum_k \beta_k \cdot \Delta F_k(t)$'
    ax1.text(0.5, -0.18, eq_text, transform=ax1.transAxes, fontsize=13,
             ha='center', style='italic')

    # --- Panel B: Top 15 |beta x MeanFreq| pairs ---
    ax2 = fig.add_subplot(gs[1])
    df_coeff = pd.read_csv(COEFF_FILE)
    nz = df_coeff[df_coeff['Beta'] != 0].copy()
    nz['AbsBxF'] = nz['Beta_x_MeanFreq'].abs()
    top15 = nz.nlargest(15, 'AbsBxF').iloc[::-1]

    labels = []
    colors = []
    for _, r in top15.iterrows():
        pair = r['ResiduePair']
        parts = pair.split('-')
        g1, g2 = int(parts[0]), int(parts[1])
        l1 = get_label(g1, resnames)
        l2 = get_label(g2, resnames)
        labels.append(f"{l1}\u2013{l2}")
        colors.append('#2166ac' if r['Beta'] < 0 else '#b2182b')

    y_pos = np.arange(len(top15))
    ax2.barh(y_pos, top15['AbsBxF'].values, color=colors, edgecolor='white', linewidth=0.5, height=0.7)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(labels, fontsize=8.5, fontfamily='monospace')
    ax2.set_xlabel(r'$|\beta_k \times \bar{F}_k|$', fontsize=11)
    ax2.set_title(r'B) Top Contact Contributions', fontsize=12, fontweight='bold')
    legend_b = [
        Patch(facecolor='#2166ac', label='Favorable ($\\beta < 0$)'),
        Patch(facecolor='#b2182b', label='Unfavorable ($\\beta > 0$)')
    ]
    ax2.legend(handles=legend_b, loc='lower right', fontsize=8, framealpha=0.9)

    # --- Panel C: NetFavorability by antibody residue ---
    ax3 = fig.add_subplot(gs[2])
    df_nf_plot = df_nf[df_nf['NetFavorability'] > 0.5].copy()
    df_nf_plot = df_nf_plot.sort_values('NetFavorability', ascending=True)

    bar_colors = []
    for _, r in df_nf_plot.iterrows():
        clf = r['Classification']
        if clf == 'hot':
            bar_colors.append('#d62728')
        elif clf == 'warm':
            if r['DominantSign'] == 'unfavorable':
                bar_colors.append('#ff7f0e')
            else:
                bar_colors.append('#2ca02c')
        else:
            bar_colors.append('#7f7f7f')

    y_pos_nf = np.arange(len(df_nf_plot))
    ax3.barh(y_pos_nf, df_nf_plot['NetFavorability'].values, color=bar_colors,
             edgecolor='white', linewidth=0.5, height=0.7)
    ax3.set_yticks(y_pos_nf)
    ax3.set_yticklabels(df_nf_plot['Label'].values, fontsize=8.5, fontfamily='monospace')
    ax3.set_xlabel('NetFavorability  ' + r'$\left(\sum |\beta_k \times \bar{F}_k|\right)$', fontsize=10)
    ax3.set_title('C) Antibody Residue Ranking', fontsize=12, fontweight='bold')
    legend_c = [
        Patch(facecolor='#d62728', label='Hot (high leverage)'),
        Patch(facecolor='#2ca02c', label='Warm (favorable)'),
        Patch(facecolor='#ff7f0e', label='Warm (unfavorable)'),
    ]
    ax3.legend(handles=legend_c, loc='lower right', fontsize=8, framealpha=0.9)
    all_nf = df_nf['NetFavorability'].values
    p70 = np.percentile(all_nf[all_nf > 0], 70)
    ax3.axvline(p70, color='red', ls='--', lw=1, alpha=0.6)
    ax3.text(p70 + 0.5, len(df_nf_plot) - 1, '70th\npctl', fontsize=7, color='red', va='top')

    plt.savefig(os.path.join(FIG_DIR, 'fig_regression.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'fig_regression.pdf'), bbox_inches='tight')
    print("Saved fig_regression.png/pdf")
    plt.close()


# ==============================================================================
# FIGURE 2: MM-GBSA Bar Chart (only corrected regression targets)
# ==============================================================================
def make_figure2():
    print("=" * 70)
    print("FIGURE 2: MM-GBSA Validation (corrected regression only)")
    print("=" * 70)

    names = [d[0] for d in MMGBSA_DATA]
    dGs = np.array([d[1] for d in MMGBSA_DATA])
    SEMs = np.array([d[2] for d in MMGBSA_DATA])
    NFs = np.array([d[3] for d in MMGBSA_DATA])
    chems = [d[4] for d in MMGBSA_DATA]
    colors = [CHEM_COLORS[c] for c in chems]

    fig, ax = plt.subplots(figsize=(13, 7))
    x = np.arange(len(names))
    bars = ax.bar(x, dGs, yerr=SEMs, capsize=4, color=colors, edgecolor='white',
                  linewidth=0.8, width=0.7, error_kw=dict(lw=1.2, capthick=1.2))

    # WT baseline
    ax.axhline(WT_DG, color='#555555', ls='--', lw=1.5, alpha=0.7)
    ax.text(len(names) - 0.5, WT_DG + 0.8, f'WT = {WT_DG:.1f}', fontsize=9,
            ha='right', color='#555555', fontweight='bold')

    # Annotate ΔΔG and NF on each bar
    for i, (name, dg, sem, nf, chem) in enumerate(MMGBSA_DATA):
        if name in ('WT', 'Y157A'):
            continue
        ddg = dg - WT_DG
        sign = '+' if ddg > 0 else ''
        va = 'top' if dg < WT_DG else 'bottom'
        offset = -1.8 if dg < WT_DG else 1.8
        ax.text(i, dg + offset, f'{sign}{ddg:.1f}\n(NF={nf:.0f})', ha='center', va=va,
                fontsize=7.5, fontweight='bold', color='black')

    # Mark Y157A
    y157_idx = names.index('Y157A')
    ax.annotate('Antigen residue\n(negative control)',
                xy=(y157_idx, dGs[y157_idx]),
                xytext=(y157_idx - 4, -5),
                fontsize=9, ha='center',
                arrowprops=dict(arrowstyle='->', color='#d62728', lw=1.5),
                color='#d62728', fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=10, rotation=35, ha='right')
    ax.set_ylabel(r'$\Delta G_{bind}$ (kcal/mol)', fontsize=12)
    ax.set_title('MM-GBSA: Mutations Identified by Corrected Elastic Net Regression',
                 fontsize=13, fontweight='bold')

    # Arrow for binding direction
    ax.annotate('', xy=(-0.08, 0.15), xytext=(-0.08, 0.85),
                xycoords='axes fraction',
                arrowprops=dict(arrowstyle='->', color='green', lw=2))
    ax.text(-0.1, 0.5, 'Stronger\nbinding', transform=ax.transAxes,
            fontsize=9, color='green', va='center', ha='center', rotation=90)

    legend_elements = [
        Patch(facecolor='#1f77b4', label='Charge removal (D\u2192N)'),
        Patch(facecolor='#2ca02c', label='Aromatic extension (\u2192Trp)'),
        Patch(facecolor='#9467bd', label='Aromatic addition (\u2192Tyr)'),
        Patch(facecolor='#ff7f0e', label='Unfavorable removal (\u2192Gly)'),
        Patch(facecolor='#d62728', label='Failed / Negative control'),
        Patch(facecolor='#555555', label='WT baseline'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=9, framealpha=0.9,
              bbox_to_anchor=(0.02, 0.98))

    ax.set_xlim(-0.6, len(names) - 0.4)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.savefig(os.path.join(FIG_DIR, 'fig_mmgbsa_validation.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'fig_mmgbsa_validation.pdf'), bbox_inches='tight')
    print("Saved fig_mmgbsa_validation.png/pdf")
    plt.close()


# ==============================================================================
# FIGURE 3: NetFavorability vs ΔΔG scatter
# ==============================================================================
def make_figure3():
    print("=" * 70)
    print("FIGURE 3: Regression Signal vs Binding Improvement")
    print("=" * 70)

    # Only mutations with both NF and ΔΔG (exclude WT, Y157A)
    scatter_data = [d for d in MMGBSA_DATA if d[0] not in ('WT', 'Y157A')]

    fig, ax = plt.subplots(figsize=(8, 6))

    for name, dg, sem, nf, chem in scatter_data:
        ddg = dg - WT_DG
        color = CHEM_COLORS[chem]
        marker = 'o' if ddg < 0 else 'X'
        size = 120 if ddg < 0 else 100
        ax.scatter(nf, ddg, c=color, s=size, marker=marker, edgecolors='black',
                   linewidths=0.8, zorder=5)
        # Label each point
        x_off = 6
        y_off = 4
        if name == 'Y404W':
            y_off = 3
            x_off = 6
        elif name == 'Y404H':
            y_off = -8
            x_off = 6
        elif name == 'F287W':
            y_off = 6
            x_off = 4
        elif name == 'N248Y':
            x_off = 6
            y_off = -8
        elif name == 'A246Y':
            x_off = 6
            y_off = 4
        elif name == 'L306G':
            x_off = 6
            y_off = -4
        ax.annotate(name, (nf, ddg), xytext=(x_off, y_off),
                    textcoords='offset points', fontsize=8.5, fontweight='bold')

    # Zero line
    ax.axhline(0, color='#555555', ls='--', lw=1, alpha=0.5)
    ax.text(70, 1, 'WT baseline', fontsize=8, color='#555555')

    # Shade improvement zone
    ax.axhspan(ax.get_ylim()[0], 0, alpha=0.05, color='green')

    ax.set_xlabel('NetFavorability (regression signal)', fontsize=12)
    ax.set_ylabel(r'$\Delta\Delta G$ vs WT (kcal/mol)', fontsize=12)
    ax.set_title('Regression Signal vs Binding Improvement', fontsize=13, fontweight='bold')

    legend_elements = [
        Patch(facecolor='#1f77b4', label='Charge removal'),
        Patch(facecolor='#2ca02c', label='Aromatic extension'),
        Patch(facecolor='#9467bd', label='Aromatic addition'),
        Patch(facecolor='#ff7f0e', label='Unfavorable removal'),
        Patch(facecolor='#d62728', label='Non-conservative'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=9, framealpha=0.9)

    # Adjust limits for label visibility
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    ax.set_xlim(xlims[0] - 2, xlims[1] + 5)
    ax.set_ylim(min(ylims[0], -38), max(ylims[1], 14))

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.savefig(os.path.join(FIG_DIR, 'fig_nf_vs_ddg.png'), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(FIG_DIR, 'fig_nf_vs_ddg.pdf'), bbox_inches='tight')
    print("Saved fig_nf_vs_ddg.png/pdf")
    plt.close()


if __name__ == '__main__':
    os.makedirs(FIG_DIR, exist_ok=True)
    make_figure1()
    make_figure2()
    make_figure3()
    print("\nDone!")
