#!/usr/bin/env python3
"""
Methods figure: regression equation and variable definitions.
Single clean PDF showing what is being regressed.
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
        r'$\Delta E(t) = \sum_{k} \beta_k \, \Delta F_k(t)$',
        transform=ax.transAxes, fontsize=28, ha='center', va='center',
        fontfamily='serif')

# Variable definitions
defs = [
    (r'$\Delta E(t)$',
     'Change in antigen\u2013antibody interaction energy\n'
     '(Coulomb + LJ) at pulling frame $t$ relative to frame 0'),
    (r'$\Delta F_k(t)$',
     'Change in pairwise contact frequency for\n'
     'residue pair $k$ at frame $t$ relative to frame 0'),
    (r'$\beta_k$',
     'Elastic net regression coefficient for pair $k$\n'
     '(sparse, L1/L2 regularized)'),
]

y_start = 0.55
y_step = 0.18
for i, (sym, desc) in enumerate(defs):
    y = y_start - i * y_step
    ax.text(0.08, y, sym, transform=ax.transAxes, fontsize=16,
            ha='left', va='top', fontfamily='serif')
    ax.text(0.22, y, desc, transform=ax.transAxes, fontsize=11,
            ha='left', va='top', fontfamily='sans-serif', color='#333333')

# Footer: method summary
ax.text(0.50, 0.02,
        '10 replicas \u00d7 50 subsampled frames  |  '
        'Leave-2-replicas-out CV (45 folds)  |  '
        'Test R\u00b2 = 0.959',
        transform=ax.transAxes, fontsize=9, ha='center', va='bottom',
        color='#666666', style='italic')

for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'fig_regression_equation.{fmt}'),
                dpi=300, bbox_inches='tight', pad_inches=0.3)
plt.close()
print("Saved fig_regression_equation")
