#!/usr/bin/env python3
"""Plot pull force curves for all validated mutations + WT + Y157A control."""

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
    'legend.fontsize': 8,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})

FIG_DIR = "/home/anugraha/antibody_optimization/figures"

def read_pullf(path):
    """Read GROMACS pull force xvg file."""
    t, f = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith(('#', '@')):
                continue
            parts = line.split()
            t.append(float(parts[0]))
            f.append(float(parts[1]))
    return np.array(t), np.array(f)


def smooth(y, window=50):
    """Simple moving average."""
    kernel = np.ones(window) / window
    return np.convolve(y, kernel, mode='same')


# Mutations to plot: (label, directory, color, linestyle, linewidth)
HITS = [
    ('F287W',  '/home/anugraha/c1_F287W/pull/replica2_pullf.xvg',  '#1565C0', '-', 1.5),
    ('S247G',  '/home/anugraha/c1_S247G/pull/replica2_pullf.xvg',  '#0277BD', '-', 1.5),
    ('Y406W',  '/home/anugraha/c1_Y406W/pull/replica2_pullf.xvg',  '#00838F', '-', 1.5),
    ('S262Y',  '/home/anugraha/c1_S262Y/pull/replica2_pullf.xvg',  '#00695C', '-', 1.5),
    ('N248Y',  '/home/anugraha/c1_N248Y/pull/replica2_pullf.xvg',  '#2E7D32', '-', 1.5),
    ('Y227W',  '/home/anugraha/c1_Y227W/pull/replica2_pullf.xvg',  '#558B2F', '-', 1.5),
    ('A246Y',  '/home/anugraha/c1_A246Y/pull/replica2_pullf.xvg',  '#827717', '-', 1.5),
]

WT =    ('WT',     '/home/anugraha/c1_WT/pull/replica2_pullf.xvg',     '#555555', '--', 2.0)
CTRL =  ('Y157A',  '/home/anugraha/c1_Y157A/pull/replica2_pullf.xvg',  '#D32F2F', '--', 2.0)


# ==============================================================================
# Figure: All pull force curves overlaid
# ==============================================================================
fig, ax = plt.subplots(figsize=(10, 5))

# Plot hits
for label, path, color, ls, lw in HITS:
    t, f = read_pullf(path)
    # Convert ps to ns
    t_ns = t / 1000.0
    f_smooth = smooth(np.abs(f), window=100)
    ax.plot(t_ns, f_smooth, color=color, ls=ls, lw=lw, label=label, alpha=0.85)

# Plot WT
label, path, color, ls, lw = WT
t, f = read_pullf(path)
t_ns = t / 1000.0
f_smooth = smooth(np.abs(f), window=100)
ax.plot(t_ns, f_smooth, color=color, ls=ls, lw=lw, label=label, alpha=0.9)

# Plot Y157A control
label, path, color, ls, lw = CTRL
t, f = read_pullf(path)
t_ns = t / 1000.0
f_smooth = smooth(np.abs(f), window=100)
ax.plot(t_ns, f_smooth, color=color, ls=ls, lw=lw, label=label, alpha=0.9)

ax.set_xlabel('Time (ns)')
ax.set_ylabel('Pull Force (kJ/mol/nm)')
ax.legend(loc='upper right', ncol=3, framealpha=0.9, edgecolor='#cccccc')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlim(0, None)
ax.set_ylim(0, None)

for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'fig_pull_forces.{fmt}'))
plt.close()
print("Saved fig_pull_forces")


# ==============================================================================
# Figure: Peak rupture force bar chart
# ==============================================================================
fig2, ax2 = plt.subplots(figsize=(8, 5))

all_muts = HITS + [WT] + [CTRL]
names = []
peaks = []
colors = []

for label, path, color, ls, lw in all_muts:
    t, f = read_pullf(path)
    f_smooth = smooth(np.abs(f), window=100)
    peak = np.max(f_smooth)
    names.append(label)
    peaks.append(peak)
    colors.append(color)

# Sort by peak force descending
order = np.argsort(peaks)[::-1]
names = [names[i] for i in order]
peaks = [peaks[i] for i in order]
colors = [colors[i] for i in order]

y_pos = np.arange(len(names))
ax2.barh(y_pos, peaks, color=colors, edgecolor='white', linewidth=0.5, height=0.7)

# WT reference line
wt_idx = names.index('WT')
ax2.axvline(x=peaks[wt_idx], color='#555555', linestyle='--', linewidth=1.0, alpha=0.5)

ax2.set_yticks(y_pos)
ax2.set_yticklabels(names, fontsize=10)
ax2.set_xlabel('Peak Rupture Force (kJ/mol/nm)')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

for fmt in ['png', 'pdf']:
    fig2.savefig(os.path.join(FIG_DIR, f'fig_peak_force.{fmt}'))
plt.close()
print("Saved fig_peak_force")
