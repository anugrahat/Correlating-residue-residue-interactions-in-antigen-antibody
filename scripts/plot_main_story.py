#!/usr/bin/env python3
"""
Main story figures:
  Figure 1: NF ranking from replica CV (top 15), tested vs untested
  Figure 2: MM-GBSA validation — hits (blue), WT (gray), Y157A (red)
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os

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

FIG_DIR = "/home/anugraha/antibody_optimization/figures"
WT_DG = -25.31

# Tested residues (residue label -> mutation name)
TESTED = {
    'Y406': 'Y406W', 'Y227': 'Y227W', 'N248': 'N248Y',
    'D245': 'D245N', 'F287': 'F287W', 'S262': 'S262Y',
    'S247': 'S247G', 'A246': 'A246Y',
}

# ==============================================================================
# Figure 1: NF ranking (top 15)
# ==============================================================================
def fig_nf_ranking():
    nf_df = pd.read_csv("/home/anugraha/c1_WT/analysis/net_favorability_replicaCV.csv")
    top = nf_df.head(15)

    labels = top['Label'].values
    nfs = top['NetFavorability'].values

    color_tested = '#2196F3'
    color_untested = '#BDBDBD'

    colors = [color_tested if lab in TESTED else color_untested for lab in labels]

    fig, ax = plt.subplots(figsize=(8, 4.5))
    x_pos = np.arange(len(labels))
    ax.bar(x_pos, nfs, color=colors, edgecolor='white', linewidth=0.5, width=0.7)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, fontsize=10, rotation=45, ha='right')
    ax.set_ylabel('Net Favorability (NF)')
    ax.set_ylim(0, max(nfs) * 1.15)

    handles = [
        Patch(facecolor=color_tested, label='Tested'),
        Patch(facecolor=color_untested, label='Untested'),
    ]
    ax.legend(handles=handles, loc='upper right', framealpha=0.9, edgecolor='#cccccc')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    for fmt in ['png', 'pdf']:
        fig.savefig(os.path.join(FIG_DIR, f'fig_main_nf_ranking.{fmt}'))
    plt.close()
    print("Saved fig_main_nf_ranking")


# ==============================================================================
# Figure 2: MM-GBSA validation
# ==============================================================================
def fig_mmgbsa_validation():
    # (name, dG, sem, type)  type: 'hit', 'wt', 'control'
    data = [
        ('F287W',  -58.99, 1.44, 'hit'),
        ('S247G',  -56.87, 1.20, 'hit'),
        ('D245N',  -43.47, 1.10, 'hit'),
        ('Y406W',  -41.05, 1.22, 'hit'),
        ('S262Y',  -40.28, 1.20, 'hit'),
        ('N248Y',  -36.20, 1.01, 'hit'),
        ('Y227W',  -35.99, 0.58, 'hit'),
        ('A246Y',  -31.41, 0.81, 'hit'),
        ('WT',     -25.31, 0.86, 'wt'),
        ('Y157A',  -13.73, 1.49, 'control'),
    ]

    # Already sorted by dG (strongest first)
    names = [d[0] for d in data]
    dgs = [d[1] for d in data]
    sems = [d[2] for d in data]
    types = [d[3] for d in data]

    color_map = {'hit': '#2196F3', 'wt': '#757575', 'control': '#D32F2F'}
    colors = [color_map[t] for t in types]

    fig, ax = plt.subplots(figsize=(8, 5))
    x_pos = np.arange(len(names))

    bars = ax.bar(x_pos, dgs, yerr=sems, color=colors, edgecolor='white',
                  linewidth=0.5, width=0.7, capsize=4, error_kw={'linewidth': 1})

    # WT reference line
    ax.axhline(y=WT_DG, color='#757575', linestyle='--', linewidth=1.0, alpha=0.5)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(names, fontsize=10, rotation=45, ha='right')
    ax.set_ylabel('MM-GBSA ΔG (kcal/mol)')

    # y-axis: more negative = stronger binding at bottom
    ax.invert_yaxis()

    handles = [
        Patch(facecolor=color_map['hit'], label='Regression-guided mutation'),
        Patch(facecolor=color_map['wt'], label='Wild type'),
        Patch(facecolor=color_map['control'], label='Negative control (Y157A)'),
    ]
    ax.legend(handles=handles, loc='lower right', framealpha=0.9, edgecolor='#cccccc')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    for fmt in ['png', 'pdf']:
        fig.savefig(os.path.join(FIG_DIR, f'fig_main_mmgbsa.{fmt}'))
    plt.close()
    print("Saved fig_main_mmgbsa")


if __name__ == '__main__':
    fig_nf_ranking()
    fig_mmgbsa_validation()
    print("\nDone.")
