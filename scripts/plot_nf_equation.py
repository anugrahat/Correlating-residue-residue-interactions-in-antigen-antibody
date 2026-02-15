#!/usr/bin/env python3
"""
Methods figure: Net Favorability equation and variable definitions.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os

FIG_DIR = "/home/anugraha/antibody_optimization/figures"

fig, ax = plt.subplots(figsize=(8, 4.5))
ax.axis('off')

# Equation
ax.text(0.50, 0.82,
        r'$\mathrm{NF}_j = \sum_{k \in \mathcal{P}_j} |\beta_k \times \bar{F}_k|$',
        transform=ax.transAxes, fontsize=28, ha='center', va='center',
        fontfamily='serif')

# Variable definitions
defs = [
    (r'$\mathrm{NF}_j$',
     'Net Favorability for antibody residue $j$\n'
     '(per-residue importance at the interface)'),
    (r'$\mathcal{P}_j$',
     'Set of all robust antigen\u2013antibody contact pairs\n'
     'involving residue $j$ (selected in >50% of CV folds)'),
    (r'$\beta_k$',
     'Elastic net coefficient for contact pair $k$'),
    (r'$\bar{F}_k$',
     'Mean contact frequency for pair $k$\n'
     '(averaged across all replicas and frames)'),
]

y_start = 0.58
y_step = 0.14
for i, (sym, desc) in enumerate(defs):
    y = y_start - i * y_step
    ax.text(0.08, y, sym, transform=ax.transAxes, fontsize=16,
            ha='left', va='top', fontfamily='serif')
    ax.text(0.22, y, desc, transform=ax.transAxes, fontsize=11,
            ha='left', va='top', fontfamily='sans-serif', color='#333333')

# Footer
ax.text(0.50, 0.02,
        'Higher NF = more important at the binding interface  |  '
        '26 antibody interface residues identified',
        transform=ax.transAxes, fontsize=9, ha='center', va='bottom',
        color='#666666', style='italic')

for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'fig_nf_equation.{fmt}'),
                dpi=300, bbox_inches='tight', pad_inches=0.3)
plt.close()
print("Saved fig_nf_equation")
