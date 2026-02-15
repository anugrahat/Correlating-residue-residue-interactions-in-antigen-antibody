#!/usr/bin/env python3
"""
Figures for NF-magnitude + chemistry-guided design story.
Uses replica CV regression results only.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
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

FIG_DIR = "/home/anugraha/antibody_optimization/figures"
os.makedirs(FIG_DIR, exist_ok=True)

WT_DG = -25.31

# ==============================================================================
# Data: all MM-GBSA results with NF from replica CV
# ==============================================================================
# (name, dG, sem, NF, strategy)
# strategy: 'improve' = chemistry substitution, 'remove' = side chain removal,
#           'control' = WT or alanine scan, 'remove_important' = removal at high NF
DATA = [
    ('F287W',  -58.99, 1.44, 17.2,  'improve'),
    ('S247G',  -56.87, 1.20,  7.7,  'remove'),
    ('D245N',  -43.47, 1.10, 29.6,  'improve'),
    ('Y406W',  -41.05, 1.22, 72.5,  'improve'),
    ('S262Y',  -40.28, 1.20,  9.6,  'improve'),
    ('N248Y',  -36.20, 1.01, 45.7,  'improve'),
    ('Y227W',  -35.99, 0.58, 65.9,  'improve'),
    ('A246Y',  -31.41, 0.81,  1.8,  'improve'),
    ('N288G',  -29.08, 1.25,  0.5,  'remove'),
    ('WT',     -25.31, 0.86,  None, 'control'),
    ('Y404W',  -28.00, 0.96, 30.0,  'improve'),
    ('Y404H',  -22.79, 0.77, 30.0,  'improve'),
    ('Y157A',  -13.73, 1.49,  None, 'control'),
    ('Y404G',   -9.97, 1.18, 30.0,  'remove_important'),
    ('H225G',   -5.31, 0.99, 43.4,  'remove_important'),
]

# Sort by dG for the bar chart
DATA_SORTED = sorted(DATA, key=lambda x: x[1])

# Colors
COLORS = {
    'improve':          '#2196F3',   # blue
    'remove':           '#4CAF50',   # green
    'remove_important': '#F44336',   # red
    'control':          '#9E9E9E',   # gray
}

LABELS = {
    'improve':          'Improve chemistry',
    'remove':           'Remove low-NF side chain',
    'remove_important': 'Remove high-NF side chain',
    'control':          'Control (WT / Ala scan)',
}


# ==============================================================================
# Figure 1: MM-GBSA bar chart colored by design strategy
# ==============================================================================
def fig_mmgbsa_strategy():
    fig, ax = plt.subplots(figsize=(8, 6))

    names = [d[0] for d in DATA_SORTED]
    dgs = [d[1] for d in DATA_SORTED]
    sems = [d[2] for d in DATA_SORTED]
    nfs = [d[3] for d in DATA_SORTED]
    strategies = [d[4] for d in DATA_SORTED]
    colors = [COLORS[s] for s in strategies]

    y_pos = np.arange(len(names))
    bars = ax.barh(y_pos, dgs, xerr=sems, color=colors, edgecolor='white',
                   linewidth=0.5, height=0.7, capsize=3, error_kw={'linewidth': 1})

    # WT reference line
    ax.axvline(x=WT_DG, color='#333333', linestyle='--', linewidth=1.0, alpha=0.7)
    ax.text(WT_DG - 0.5, len(names) - 0.5, 'WT', fontsize=9, ha='right',
            color='#333333', fontstyle='italic')

    # NF annotations on bars
    for i, (name, dg, sem, nf, strat) in enumerate(DATA_SORTED):
        if nf is not None:
            label = f'NF={nf:.0f}' if nf >= 1 else f'NF={nf:.1f}'
            if dg < WT_DG:
                ax.text(dg - sem - 1.0, i, label, va='center', ha='right',
                        fontsize=7.5, color='#555555')
            else:
                ax.text(dg + sem + 1.0, i, label, va='center', ha='left',
                        fontsize=7.5, color='#555555')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(names, fontsize=10)
    ax.set_xlabel('MM-GBSA ΔG (kcal/mol)')
    ax.invert_xaxis()

    # Legend
    handles = [Patch(facecolor=COLORS[k], label=LABELS[k])
               for k in ['improve', 'remove', 'remove_important', 'control']]
    ax.legend(handles=handles, loc='lower right', framealpha=0.9,
              edgecolor='#cccccc')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    for fmt in ['png', 'pdf']:
        fig.savefig(os.path.join(FIG_DIR, f'fig_strategy_validation.{fmt}'))
    plt.close()
    print("Saved fig_strategy_validation")


# ==============================================================================
# Figure 2: NF magnitude vs ΔΔG — substitution vs removal
# ==============================================================================
def fig_nf_vs_ddg():
    fig, ax = plt.subplots(figsize=(7, 5.5))

    for name, dg, sem, nf, strat in DATA:
        if nf is None:
            continue
        ddg = dg - WT_DG
        if strat == 'improve':
            marker = 'o'
            color = COLORS['improve']
            ms = 9
        elif strat == 'remove':
            marker = 's'
            color = COLORS['remove']
            ms = 9
        elif strat == 'remove_important':
            marker = 's'
            color = COLORS['remove_important']
            ms = 9
        else:
            continue

        ax.errorbar(nf, ddg, yerr=sem, fmt=marker, color=color, ms=ms,
                     markeredgecolor='white', markeredgewidth=0.8,
                     ecolor=color, elinewidth=1, capsize=3, zorder=5)
        # Label
        offset_x = 1.5
        offset_y = 0
        ha = 'left'
        if name == 'Y404W':
            offset_y = -2.0
        elif name == 'Y404H':
            offset_y = 2.0
        elif name == 'Y404G':
            offset_x = -1.5
            ha = 'right'
        elif name == 'A246Y':
            offset_x = 1.5
        elif name == 'N288G':
            offset_y = 1.5

        ax.text(nf + offset_x, ddg + offset_y, name, fontsize=8, ha=ha,
                va='center', color='#333333')

    # Reference line at ΔΔG = 0
    ax.axhline(y=0, color='#333333', linestyle='--', linewidth=1.0, alpha=0.5)
    ax.text(75, 1.5, 'WT baseline', fontsize=8, color='#555555', ha='right',
            fontstyle='italic')

    # Quadrant labels
    ax.text(60, -25, 'High NF\nsubstitution\n→ hits', fontsize=9,
            ha='center', va='center', color=COLORS['improve'], alpha=0.6,
            fontweight='bold')
    ax.text(60, 15, 'High NF\nremoval\n→ failure', fontsize=9,
            ha='center', va='center', color=COLORS['remove_important'], alpha=0.6,
            fontweight='bold')
    ax.text(5, -25, 'Low NF\nremoval\n→ hits', fontsize=9,
            ha='center', va='center', color=COLORS['remove'], alpha=0.6,
            fontweight='bold')

    ax.set_xlabel('Net Favorability (NF magnitude)')
    ax.set_ylabel('ΔΔG vs WT (kcal/mol)')

    # Legend
    handles = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['improve'],
                   markersize=9, label='Substitution (improve chemistry)'),
        plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=COLORS['remove'],
                   markersize=9, label='Removal at low NF'),
        plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=COLORS['remove_important'],
                   markersize=9, label='Removal at high NF'),
    ]
    ax.legend(handles=handles, loc='upper left', framealpha=0.9,
              edgecolor='#cccccc')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlim(-3, 80)

    for fmt in ['png', 'pdf']:
        fig.savefig(os.path.join(FIG_DIR, f'fig_nf_vs_ddg.{fmt}'))
    plt.close()
    print("Saved fig_nf_vs_ddg")


# ==============================================================================
# Figure 3: NF ranking from CV with tested mutations annotated
# ==============================================================================
def fig_nf_ranking_annotated():
    # Load CV NF data
    nf_df = pd.read_csv("/home/anugraha/c1_WT/analysis/net_favorability_replicaCV.csv")

    # Map of tested mutations and their outcomes
    tested = {
        'Y406': ('Y406W', -15.7, 'improve'),
        'Y227': ('Y227W', -10.7, 'improve'),
        'N248': ('N248Y', -10.9, 'improve'),
        'H225': ('H225G', +20.0, 'remove_important'),
        'Y404': ('Y404G', +15.3, 'remove_important'),
        'D245': ('D245N', -18.2, 'improve'),
        'A255': (None, None, None),  # untested
        'F287': ('F287W', -33.7, 'improve'),
        'S262': ('S262Y', -15.0, 'improve'),
        'S247': ('S247G', -31.6, 'remove'),
        'A246': ('A246Y', -6.1, 'improve'),
        'N288': ('N288G', -3.8, 'remove'),
    }

    fig, ax = plt.subplots(figsize=(10, 5))

    # Top 15 residues
    top = nf_df.head(15)
    labels = top['Label'].values
    nfs = top['NetFavorability'].values

    colors = []
    for lab in labels:
        if lab in tested:
            mut, ddg, strat = tested[lab]
            if strat is None:
                colors.append('#BBBBBB')  # untested
            else:
                colors.append(COLORS[strat])
        else:
            colors.append('#BBBBBB')  # untested

    x_pos = np.arange(len(labels))
    bars = ax.bar(x_pos, nfs, color=colors, edgecolor='white', linewidth=0.5,
                  width=0.7)

    # Annotate with mutation name and ΔΔG
    for i, lab in enumerate(labels):
        if lab in tested:
            mut, ddg, strat = tested[lab]
            if mut is not None:
                sign = '+' if ddg > 0 else ''
                ax.text(i, nfs[i] + 1.5, f'{mut}\n{sign}{ddg:.1f}',
                        ha='center', va='bottom', fontsize=7.5,
                        fontweight='bold', color='#333333')
            else:
                ax.text(i, nfs[i] + 1.5, 'untested',
                        ha='center', va='bottom', fontsize=7.5,
                        fontstyle='italic', color='#999999')

    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, fontsize=10, rotation=45, ha='right')
    ax.set_ylabel('Net Favorability (NF)')
    ax.set_ylim(0, max(nfs) * 1.35)

    # Legend
    handles = [
        Patch(facecolor=COLORS['improve'], label='Substitution → hit'),
        Patch(facecolor=COLORS['remove'], label='Low-NF removal → hit'),
        Patch(facecolor=COLORS['remove_important'], label='High-NF removal → failure'),
        Patch(facecolor='#BBBBBB', label='Untested'),
    ]
    ax.legend(handles=handles, loc='upper right', framealpha=0.9,
              edgecolor='#cccccc')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    for fmt in ['png', 'pdf']:
        fig.savefig(os.path.join(FIG_DIR, f'fig_nf_ranking_annotated.{fmt}'))
    plt.close()
    print("Saved fig_nf_ranking_annotated")


# ==============================================================================
# Figure 4: Design rules summary — grouped bar
# ==============================================================================
def fig_design_rules():
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [2, 1]})

    # Panel A: ΔΔG grouped by strategy
    ax = axes[0]

    groups = {
        'Improve chemistry\n(substitution at any NF)': [
            ('F287W',  -33.7, 17.2),
            ('D245N',  -18.2, 29.6),
            ('Y406W',  -15.7, 72.5),
            ('S262Y',  -15.0,  9.6),
            ('N248Y',  -10.9, 45.7),
            ('Y227W',  -10.7, 65.9),
            ('A246Y',   -6.1,  1.8),
        ],
        'Remove side chain\nat low NF': [
            ('S247G',  -31.6,  7.7),
            ('N288G',   -3.8,  0.5),
        ],
        'Remove side chain\nat high NF': [
            ('H225G',  +20.0, 43.4),
            ('Y404G',  +15.3, 30.0),
        ],
        'Controls': [
            ('Y157A',  +11.6, None),
        ],
    }

    group_colors = [COLORS['improve'], COLORS['remove'],
                    COLORS['remove_important'], COLORS['control']]

    x = 0
    xticks = []
    xticklabels = []
    group_centers = []
    group_names = list(groups.keys())

    for gi, (gname, muts) in enumerate(groups.items()):
        positions = []
        for name, ddg, nf in muts:
            bar = ax.bar(x, ddg, width=0.7, color=group_colors[gi],
                        edgecolor='white', linewidth=0.5)
            xticks.append(x)
            xticklabels.append(name)
            positions.append(x)
            x += 1
        group_centers.append(np.mean(positions))
        x += 0.8  # gap between groups

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, fontsize=8.5, rotation=45, ha='right')
    ax.axhline(y=0, color='#333333', linewidth=0.8)
    ax.set_ylabel('ΔΔG vs WT (kcal/mol)')

    # Legend for Panel A
    handles_a = [Patch(facecolor=COLORS[k], label=LABELS[k])
                 for k in ['improve', 'remove', 'remove_important', 'control']]
    ax.legend(handles=handles_a, loc='lower left', framealpha=0.9,
              edgecolor='#cccccc', fontsize=8)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(-40, 28)
    ax.text(-0.08, 1.02, 'A', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='bottom')

    # Panel B: hit rate summary
    ax2 = axes[1]

    categories = ['Substitution\n(any NF)', 'Removal\n(low NF)', 'Removal\n(high NF)']
    # hits: ΔΔG < -5. Substitution 7/7, low-NF removal 1/2 (S247G hit, N288G marginal)
    hits = [7, 1, 0]
    totals = [7, 2, 2]
    rates = [h/t * 100 for h, t in zip(hits, totals)]
    cat_colors = [COLORS['improve'], COLORS['remove'], COLORS['remove_important']]

    bars2 = ax2.bar(range(3), rates, color=cat_colors, edgecolor='white',
                    linewidth=0.5, width=0.6)

    for i, (h, t, r) in enumerate(zip(hits, totals, rates)):
        ax2.text(i, r + 2, f'{h}/{t}', ha='center', va='bottom',
                fontsize=11, fontweight='bold', color='#333333')

    ax2.set_xticks(range(3))
    ax2.set_xticklabels(categories, fontsize=9)
    ax2.set_ylabel('Hit rate (%)')
    ax2.set_ylim(0, 120)
    ax2.axhline(y=50, color='#cccccc', linestyle=':', linewidth=0.8)

    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.text(-0.12, 1.02, 'B', transform=ax2.transAxes, fontsize=14,
            fontweight='bold', va='bottom')

    plt.tight_layout()

    for fmt in ['png', 'pdf']:
        fig.savefig(os.path.join(FIG_DIR, f'fig_design_rules.{fmt}'))
    plt.close()
    print("Saved fig_design_rules")


# ==============================================================================
# Run all
# ==============================================================================
if __name__ == '__main__':
    fig_mmgbsa_strategy()
    fig_nf_vs_ddg()
    fig_nf_ranking_annotated()
    fig_design_rules()
    print("\nAll figures saved to", FIG_DIR)
