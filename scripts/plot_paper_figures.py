#!/usr/bin/env python3
"""
Publication-quality figures for elastic net regression + MM-GBSA validation.
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
# Global style
# ==============================================================================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10,
    'axes.linewidth': 1.0,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})

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

# All MM-GBSA results — ONLY residues present in the corrected regression
# Excluded: D265N (D265 SASA-filtered), S258G/S258A (S258 not in regression),
#           T403Y (T403 not in regression), Y157A (Y157 not in regression)
ALL_MMGBSA = [
    ('F287W',  -58.99, 1.44, 'Aromatic extension'),
    ('S247G',  -56.87, 1.20, 'Side chain removal'),
    ('D245N',  -43.47, 1.10, 'Charge removal'),
    ('T330G',  -41.47, 0.65, 'Side chain removal'),
    ('S262Y',  -40.28, 1.20, 'Aromatic addition'),
    ('N248Y',  -36.20, 1.01, 'Aromatic addition'),
    ('Y227W',  -35.99, 0.58, 'Aromatic extension'),
    ('A246Y',  -31.41, 0.81, 'Aromatic addition'),
    ('Y404W',  -28.00, 0.96, 'Aromatic extension'),
    ('WT',     -25.31, 0.86, 'Baseline'),
    ('Y404H',  -22.79, 0.77, 'Failed'),
    ('L306G',  -16.35, 0.60, 'Failed'),
    ('Y157A',  -13.73, 1.49, 'Failed'),
]

CHEM_COLORS = {
    'Charge removal':     '#1976D2',
    'Aromatic extension': '#388E3C',
    'Aromatic addition':  '#7B1FA2',
    'Side chain removal': '#F57C00',
    'Baseline':           '#616161',
    'Failed':             '#C62828',
    'Negative control':   '#C62828',
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
    print("FIGURE 1: Elastic Net Regression")
    print("=" * 70)

    resnames = load_gro_resnames(GRO_PATH)

    # Load data
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

    # Fit
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
    df_nf = pd.read_csv(NF_FILE)

    # ---- Figure ----
    fig = plt.figure(figsize=(18, 5.5))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1.2, 1.3], wspace=0.32)

    # --- Panel A ---
    ax1 = fig.add_subplot(gs[0])
    ax1.scatter(y, y_pred, s=4, alpha=0.35, c='#1565C0', edgecolors='none', rasterized=True)
    lims = [min(y.min(), y_pred.min()) - 10, max(y.max(), y_pred.max()) + 10]
    ax1.plot(lims, lims, 'k-', lw=0.8, alpha=0.4)
    ax1.set_xlim(lims); ax1.set_ylim(lims)
    ax1.set_xlabel(r'Actual $\Delta E$ (kJ/mol)')
    ax1.set_ylabel(r'Predicted $\Delta E$ (kJ/mol)')
    ax1.set_title('A', fontweight='bold', loc='left', fontsize=14)
    stats = (f'$R^2$ = {r2:.4f}\nRMSE = {rmse:.1f} kJ/mol\n'
             f'{n_nz}/{len(pair_names)} non-zero $\\beta$')
    ax1.text(0.05, 0.95, stats, transform=ax1.transAxes, fontsize=9,
             va='top', bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#ccc', alpha=0.9))
    eq = r'$\Delta E(t) = \sum_k \beta_k \cdot \Delta F_k(t)$'
    ax1.text(0.5, -0.15, eq, transform=ax1.transAxes, fontsize=12, ha='center', style='italic')

    # --- Panel B ---
    ax2 = fig.add_subplot(gs[1])
    df_coeff = pd.read_csv(COEFF_FILE)
    nz = df_coeff[df_coeff['Beta'] != 0].copy()
    nz['AbsBxF'] = nz['Beta_x_MeanFreq'].abs()
    top15 = nz.nlargest(15, 'AbsBxF').iloc[::-1]
    labels, colors = [], []
    for _, r in top15.iterrows():
        g1, g2 = [int(x) for x in r['ResiduePair'].split('-')]
        labels.append(f"{get_label(g1, resnames)}\u2013{get_label(g2, resnames)}")
        colors.append('#1565C0' if r['Beta'] < 0 else '#C62828')
    y_pos = np.arange(len(top15))
    ax2.barh(y_pos, top15['AbsBxF'].values, color=colors, edgecolor='white', lw=0.3, height=0.65)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(labels, fontsize=8, fontfamily='monospace')
    ax2.set_xlabel(r'$|\beta_k \times \bar{F}_k|$')
    ax2.set_title('B', fontweight='bold', loc='left', fontsize=14)
    ax2.legend(handles=[Patch(fc='#1565C0', label='Stabilizing ($\\beta<0$)'),
                        Patch(fc='#C62828', label='Destabilizing ($\\beta>0$)')],
               loc='lower right', fontsize=8, framealpha=0.9)

    # --- Panel C ---
    ax3 = fig.add_subplot(gs[2])
    df_nf_plot = df_nf[df_nf['NetFavorability'] > 0.5].sort_values('NetFavorability', ascending=True)
    bar_c = []
    for _, r in df_nf_plot.iterrows():
        if r['Classification'] == 'hot':
            bar_c.append('#C62828')
        elif r['DominantSign'] == 'unfavorable':
            bar_c.append('#F57C00')
        else:
            bar_c.append('#388E3C')
    y_nf = np.arange(len(df_nf_plot))
    ax3.barh(y_nf, df_nf_plot['NetFavorability'].values, color=bar_c, edgecolor='white', lw=0.3, height=0.65)
    ax3.set_yticks(y_nf)
    ax3.set_yticklabels(df_nf_plot['Label'].values, fontsize=8, fontfamily='monospace')
    ax3.set_xlabel(r'NetFavorability $\left(\sum |\beta_k \times \bar{F}_k|\right)$')
    ax3.set_title('C', fontweight='bold', loc='left', fontsize=14)
    ax3.legend(handles=[Patch(fc='#C62828', label='Hot spot'),
                        Patch(fc='#388E3C', label='Warm (favorable)'),
                        Patch(fc='#F57C00', label='Warm (unfavorable)')],
               loc='lower right', fontsize=8, framealpha=0.9)
    p70 = np.percentile(df_nf['NetFavorability'].values[df_nf['NetFavorability'].values > 0], 70)
    ax3.axvline(p70, color='#C62828', ls=':', lw=1, alpha=0.5)

    for ax in [ax1, ax2, ax3]:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.savefig(os.path.join(FIG_DIR, 'fig_regression.png'), dpi=300)
    plt.savefig(os.path.join(FIG_DIR, 'fig_regression.pdf'))
    print("Saved fig_regression.png/pdf")
    plt.close()


# ==============================================================================
# FIGURE 2: MM-GBSA — ALL mutations, grouped by chemistry
# ==============================================================================
def make_figure2():
    print("=" * 70)
    print("FIGURE 2: MM-GBSA Validation (all mutations)")
    print("=" * 70)

    names = [d[0] for d in ALL_MMGBSA]
    dGs = np.array([d[1] for d in ALL_MMGBSA])
    SEMs = np.array([d[2] for d in ALL_MMGBSA])
    chems = [d[3] for d in ALL_MMGBSA]
    colors = [CHEM_COLORS[c] for c in chems]

    fig, ax = plt.subplots(figsize=(12, 5.5))
    x = np.arange(len(names))

    bars = ax.bar(x, dGs, yerr=SEMs, capsize=3, color=colors, edgecolor='white',
                  linewidth=0.5, width=0.72, error_kw=dict(lw=1, capthick=1, color='#333'))

    # WT baseline
    ax.axhline(WT_DG, color='#616161', ls='--', lw=1.2, alpha=0.6)
    ax.text(len(names) - 0.3, WT_DG + 1.2, f'WT = {WT_DG:.1f}', fontsize=8,
            ha='right', color='#616161', fontweight='bold')

    # ΔΔG annotations
    for i, (name, dg, sem, chem) in enumerate(ALL_MMGBSA):
        if name == 'WT':
            continue
        ddg = dg - WT_DG
        sign = '+' if ddg > 0 else ''
        if dg < WT_DG:
            ax.text(i, dg - 1.5, f'{sign}{ddg:.0f}', ha='center', va='top',
                    fontsize=7, fontweight='bold', color='#333')
        else:
            ax.text(i, dg + 1.5, f'{sign}{ddg:.0f}', ha='center', va='bottom',
                    fontsize=7, fontweight='bold', color='#333')

    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=40, ha='right', fontsize=9)
    ax.set_ylabel(r'$\Delta G_{\mathrm{bind}}$ (kcal/mol)')

    legend_elements = [
        Patch(fc='#1976D2', label='Charge removal (D/E\u2192N/Q)'),
        Patch(fc='#388E3C', label='Aromatic extension (\u2192Trp)'),
        Patch(fc='#7B1FA2', label='Aromatic addition (\u2192Tyr)'),
        Patch(fc='#F57C00', label='Side chain removal (\u2192Gly)'),
        Patch(fc='#C62828', label='Failed / Control'),
        Patch(fc='#616161', label='WT'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=8,
              framealpha=0.95, edgecolor='#ccc')

    ax.set_xlim(-0.6, len(names) - 0.4)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.savefig(os.path.join(FIG_DIR, 'fig_mmgbsa_validation.png'), dpi=300)
    plt.savefig(os.path.join(FIG_DIR, 'fig_mmgbsa_validation.pdf'))
    print("Saved fig_mmgbsa_validation.png/pdf")
    plt.close()


# ==============================================================================
# FIGURE 3: Chemistry success rate summary
# ==============================================================================
def make_figure3():
    print("=" * 70)
    print("FIGURE 3: Chemistry Success Rates")
    print("=" * 70)

    # (chemistry, successes, failures) — only corrected-regression residues
    chem_data = [
        ('Charge removal\n(D/E\u2192N/Q)',  1, 0),
        ('Aromatic ext.\n(\u2192Trp)',       3, 0),
        ('Aromatic add.\n(\u2192Tyr)',       3, 0),
        ('Side chain rem.\n(\u2192Gly)',     2, 1),
        ('Non-conserv.\nat aromatic',        0, 1),
        ('Ala control',                      0, 1),
    ]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5),
                                    gridspec_kw={'width_ratios': [1.2, 1]})

    # --- Panel A: Stacked bar ---
    categories = [d[0] for d in chem_data]
    hits = [d[1] for d in chem_data]
    fails = [d[2] for d in chem_data]
    totals = [h + f for h, f in zip(hits, fails)]
    rates = [h / t * 100 if t > 0 else 0 for h, t in zip(hits, totals)]

    y_pos = np.arange(len(categories))
    ax1.barh(y_pos, hits, color='#388E3C', edgecolor='white', lw=0.5, height=0.55, label='Improved')
    ax1.barh(y_pos, fails, left=hits, color='#C62828', edgecolor='white', lw=0.5, height=0.55, label='Failed')

    for i, (h, f, rate) in enumerate(zip(hits, fails, rates)):
        ax1.text(h + f + 0.15, i, f'{rate:.0f}%', va='center', fontsize=10, fontweight='bold',
                 color='#388E3C' if rate > 50 else '#C62828')

    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(categories, fontsize=10)
    ax1.set_xlabel('Number of mutations tested')
    ax1.set_title('A', fontweight='bold', loc='left', fontsize=14)
    ax1.legend(loc='lower right', fontsize=9, framealpha=0.9)
    ax1.set_xlim(0, 5)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # --- Panel B: ΔΔG distribution by chemistry (corrected-regression residues only) ---
    chem_ddg = {
        'Charge removal':     [(-18.2, 'D245N')],
        'Aromatic extension': [(-33.7, 'F287W'), (-10.7, 'Y227W'), (-2.7, 'Y404W')],
        'Aromatic addition':  [(-15.0, 'S262Y'), (-10.9, 'N248Y'), (-6.1, 'A246Y')],
        'Side chain removal': [(-31.6, 'S247G'), (-16.2, 'T330G'), (+9.0, 'L306G')],
        'Non-conservative':   [(+2.5, 'Y404H')],
        'Ala control':        [(+11.6, 'Y157A')],
    }
    chem_colors_strip = {
        'Charge removal': '#1976D2',
        'Aromatic extension': '#388E3C',
        'Side chain removal': '#F57C00',
        'Aromatic addition': '#7B1FA2',
        'Non-conservative': '#C62828',
        'Ala control': '#795548',
    }

    ax2.axhline(0, color='#616161', ls='--', lw=1, alpha=0.5)
    ax2.axhspan(-45, 0, alpha=0.04, color='green')

    x_idx = 0
    xticks = []
    xticklabels = []
    for chem, vals in chem_ddg.items():
        color = chem_colors_strip[chem]
        for ddg, name in vals:
            marker = 'o' if ddg < 0 else 'X'
            ax2.scatter(x_idx, ddg, c=color, s=80, marker=marker, edgecolors='#333',
                        linewidths=0.6, zorder=5)
            ax2.text(x_idx, ddg - 2.2 if ddg < 0 else ddg + 2.2, name,
                     ha='center', va='top' if ddg < 0 else 'bottom',
                     fontsize=7, fontweight='bold', color='#333')
            xticks.append(x_idx)
            x_idx += 1
        x_idx += 0.5  # gap between groups

    ax2.set_ylabel(r'$\Delta\Delta G$ vs WT (kcal/mol)')
    ax2.set_title('B', fontweight='bold', loc='left', fontsize=14)
    ax2.set_xticks([])
    ax2.set_xlim(-1.5, x_idx + 0.5)
    ax2.set_ylim(-52, 16)
    ax2.text(0.5, 0, 'WT', fontsize=8, color='#616161', transform=ax2.get_yaxis_transform())

    # Group labels at bottom
    group_centers = []
    idx = 0
    for chem, vals in chem_ddg.items():
        center = idx + (len(vals) - 1) / 2
        group_centers.append((center, chem.split('(')[0].strip()))
        idx += len(vals) + 0.5
    for cx, label in group_centers:
        ax2.text(cx, -49, label, ha='center', va='top', fontsize=7, color='#555', style='italic',
                 rotation=0)

    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'fig_chemistry_summary.png'), dpi=300)
    plt.savefig(os.path.join(FIG_DIR, 'fig_chemistry_summary.pdf'))
    print("Saved fig_chemistry_summary.png/pdf")
    plt.close()


if __name__ == '__main__':
    os.makedirs(FIG_DIR, exist_ok=True)
    make_figure1()
    make_figure2()
    make_figure3()
    print("\nDone!")
