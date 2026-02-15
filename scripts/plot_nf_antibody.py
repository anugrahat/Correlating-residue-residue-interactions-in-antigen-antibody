#!/usr/bin/env python3
"""
Net Favorability (Σ|β_k × F_k|) ranking for antibody interface residues.
Clean figure — no favorable/unfavorable labels.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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

# Load NF data (already antibody-only)
nf = pd.read_csv("/home/anugraha/c1_WT/analysis/net_favorability_replicaCV.csv")
nf = nf.sort_values('NetFavorability', ascending=False)

labels = nf['Label'].values
nfs = nf['NetFavorability'].values

fig, ax = plt.subplots(figsize=(7, 7))
y_pos = np.arange(len(labels))

ax.barh(y_pos, nfs, color='#2196F3', edgecolor='white', linewidth=0.5, height=0.7)

ax.set_yticks(y_pos)
ax.set_yticklabels(labels, fontsize=10)
ax.set_xlabel(r'Net Favorability  $\Sigma |\beta_k \times \overline{F}_k|$')
ax.invert_yaxis()

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'fig_nf_antibody.{fmt}'))
plt.close()
print("Saved fig_nf_antibody")
