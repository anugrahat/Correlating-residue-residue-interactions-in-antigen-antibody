#!/usr/bin/env python3
"""
Figure: Top residue-pair β coefficients from replica CV regression.
Shows which antigen-antibody contacts the regression identified as important.
Y157 pairs appear naturally in the ranking.
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
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})

FIG_DIR = "/home/anugraha/antibody_optimization/figures"

# Read coefficients
df = pd.read_csv("/home/anugraha/c1_WT/analysis/elastic_net_coefficients_replicaCV.csv")
robust = df[df['Robust'] == True].copy()
robust['AbsBeta'] = robust['Beta_final'].abs()
robust = robust.sort_values('AbsBeta', ascending=False)

# Read GRO for residue names
resnames = {}
with open('/home/anugraha/c1_WT/pull/md10_protein.gro') as f:
    f.readline(); f.readline()
    for line in f:
        line = line.rstrip()
        if len(line) < 20:
            break
        try:
            resid = int(line[:5])
            resname = line[5:10].strip()
            if resid not in resnames:
                resnames[resid] = resname
        except ValueError:
            break

# 3-letter to 1-letter
AA3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'HIE': 'H',
    'HID': 'H', 'HIP': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T',
    'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

# GRO -> global numbering
def gro_to_global(gro_id):
    if gro_id <= 194:
        return gro_id + 1  # Chain A (antigen), GRO 0-indexed offset
    elif gro_id <= 314:
        return gro_id - 194 + 302  # Chain H
    else:
        return gro_id - 314 + 195  # Chain L

def gro_to_label(gro_id):
    rn = resnames.get(gro_id, 'UNK')
    aa1 = AA3TO1.get(rn, '?')
    glob = gro_to_global(gro_id)
    return f'{aa1}{glob}'

def pair_label(pair_str):
    g1, g2 = [int(x) for x in pair_str.split('-')]
    l1 = gro_to_label(g1)
    l2 = gro_to_label(g2)
    # First is antigen, second is antibody
    if g1 <= 194:
        return f'{l1} \u2013 {l2}'
    else:
        return f'{l2} \u2013 {l1}'

# Top 25 robust pairs
top = robust.head(25)

labels = [pair_label(p) for p in top['ResiduePair']]
betas = top['Beta_final'].values

# Color: highlight pairs involving Y157 (GRO 156)
colors = []
for i, pair in enumerate(top['ResiduePair']):
    g1, g2 = [int(x) for x in pair.split('-')]
    if g1 == 156 or g2 == 156:  # GRO 156 = Y157 (negative control target)
        colors.append('#D32F2F')  # red for Y157 pairs
    elif betas[i] < 0:
        colors.append('#2196F3')  # blue for negative beta
    else:
        colors.append('#FF9800')  # orange for positive beta

fig, ax = plt.subplots(figsize=(8, 7))

y_pos = np.arange(len(labels))
ax.barh(y_pos, betas, color=colors, edgecolor='white', linewidth=0.5, height=0.7)

ax.set_yticks(y_pos)
ax.set_yticklabels(labels, fontsize=9)
ax.set_xlabel('Regression coefficient (β)')
ax.axvline(x=0, color='#333333', linewidth=0.8)
ax.invert_yaxis()

from matplotlib.patches import Patch
handles = [
    Patch(facecolor='#2196F3', label='Negative β (more contact → stronger binding)'),
    Patch(facecolor='#FF9800', label='Positive β (more contact → weaker binding)'),
    Patch(facecolor='#D32F2F', label='Y157 pairs (negative control target)'),
]
ax.legend(handles=handles, loc='lower right', framealpha=0.9, edgecolor='#cccccc',
          fontsize=8)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'fig_pair_betas.{fmt}'))
plt.close()
print("Saved fig_pair_betas")
