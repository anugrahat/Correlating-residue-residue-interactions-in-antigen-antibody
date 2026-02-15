#!/usr/bin/env python3
"""
ΔΔG plot for mutations at the top 12 NF-ranked antibody interface
residues, plus Y157A antigen control.
"""

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
WT_DG = -25.31

# Top 12 NF residues: Y406, Y227, Y405, N248, H225, Y404, D245, A255,
#                     F244, F287, T403, W289
# Only mutations at these residues with GBSA data + Y157A control
# (Y405, A255, F244, W289 have no mutation GBSA data)
data = [
    ('F287W',   -58.99, 1.44, 'F287'),
    ('D245N',   -43.47, 1.10, 'D245'),
    ('Y406W',   -41.05, 1.22, 'Y406'),
    ('N248Y',   -36.20, 1.01, 'N248'),
    ('Y227W',   -35.99, 0.58, 'Y227'),
    ('Y404W',   -28.00, 0.96, 'Y404'),
    ('Y404H',   -22.79, 0.77, 'Y404'),
    ('T403Y',   -21.08, 1.57, 'T403'),
    ('Y157A',   -13.73, 1.49, 'Y157'),   # antigen control
    ('Y404G',    -9.97, 1.18, 'Y404'),
    ('H225G',    -5.31, 0.99, 'H225'),
]

# Compute ΔΔG
names = [d[0] for d in data]
ddgs = [d[1] - WT_DG for d in data]
sems = [d[2] for d in data]

# Sort by ΔΔG (most improved first)
order = np.argsort(ddgs)
names = [names[i] for i in order]
ddgs = [ddgs[i] for i in order]
sems = [sems[i] for i in order]

# Color: blue for improved, red for weakened, dark red for Y157A
colors = []
for i, name in enumerate(names):
    if name == 'Y157A':
        colors.append('#B71C1C')
    elif ddgs[i] < 0:
        colors.append('#2E7D32')
    else:
        colors.append('#D32F2F')

fig, ax = plt.subplots(figsize=(8, 6))
y_pos = np.arange(len(names))

ax.barh(y_pos, ddgs, xerr=sems, color=colors, edgecolor='white',
        linewidth=0.5, height=0.7, capsize=3, error_kw={'linewidth': 1})

ax.axvline(x=0, color='#333333', linewidth=1.0, linestyle='-')
ax.set_yticks(y_pos)
ax.set_yticklabels(names, fontsize=10)
ax.set_xlabel('\u0394\u0394G (kcal/mol)')
ax.invert_yaxis()

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Annotations
ax.text(0.02, 0.02, '\u2190 Stronger binding',
        transform=ax.transAxes, fontsize=9, color='#2E7D32', va='bottom')
ax.text(0.98, 0.02, 'Weaker binding \u2192',
        transform=ax.transAxes, fontsize=9, color='#D32F2F', va='bottom', ha='right')

for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'fig_ddg.{fmt}'))
plt.close()
print("Saved fig_ddg")
