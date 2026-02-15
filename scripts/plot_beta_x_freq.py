#!/usr/bin/env python3
"""
Figures:
  1. Top 25 residue pairs ranked by |β × mean_freq|
  2. Top 25 residue pairs ranked by |β × max_freq|
  3. Temporal β × freq(t) for key pairs during pulling
Y157 pairs highlighted in red throughout.
"""

import pandas as pd
import numpy as np
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
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})

FIG_DIR = "/home/anugraha/antibody_optimization/figures"

# =========================================================================
# Residue labeling (GRO -> human-readable)
# =========================================================================
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

AA3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'HIE': 'H',
    'HID': 'H', 'HIP': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T',
    'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

def gro_to_global(gro_id):
    if gro_id <= 194:
        return gro_id + 1      # Chain A (antigen)
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

def is_y157_pair(pair_str):
    g1, g2 = [int(x) for x in pair_str.split('-')]
    return g1 == 156 or g2 == 156

# =========================================================================
# Load regression coefficients
# =========================================================================
coef_df = pd.read_csv("/home/anugraha/c1_WT/analysis/elastic_net_coefficients_replicaCV.csv")
robust = coef_df[coef_df['Robust'] == True].copy()

# =========================================================================
# Compute max frequency per pair from temporal data
# =========================================================================
print("Loading temporal contact frequencies...")
freq_df = pd.read_csv("/home/anugraha/c1_WT/analysis/average_frequency.csv")

max_freq = freq_df.groupby('ResiduePair')['InteractionFrequency'].max().reset_index()
max_freq.columns = ['ResiduePair', 'MaxFreq']

robust = robust.merge(max_freq, on='ResiduePair', how='left')
robust['MaxFreq'] = robust['MaxFreq'].fillna(0)
robust['Beta_x_MaxFreq'] = robust['Beta_final'] * robust['MaxFreq']

# =========================================================================
# Helper: bar chart of top N pairs
# =========================================================================
def plot_pair_bars(df, value_col, abs_col, xlabel, filename, n=25):
    top = df.nlargest(n, abs_col)

    labels = [pair_label(p) for p in top['ResiduePair']]
    values = top[value_col].values

    colors = []
    for i, pair in enumerate(top['ResiduePair']):
        if is_y157_pair(pair):
            colors.append('#D32F2F')
        elif values[i] < 0:
            colors.append('#2196F3')
        else:
            colors.append('#FF9800')

    fig, ax = plt.subplots(figsize=(8, 7))
    y_pos = np.arange(len(labels))
    ax.barh(y_pos, values, color=colors, edgecolor='white', linewidth=0.5, height=0.7)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel(xlabel)
    ax.axvline(x=0, color='#333333', linewidth=0.8)
    ax.invert_yaxis()

    handles = [
        Patch(facecolor='#2196F3', label='Negative (more contact \u2192 stronger binding)'),
        Patch(facecolor='#FF9800', label='Positive (more contact \u2192 weaker binding)'),
        Patch(facecolor='#D32F2F', label='Y157 pairs (negative control target)'),
    ]
    ax.legend(handles=handles, loc='lower right', framealpha=0.9, edgecolor='#cccccc',
              fontsize=8)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    for fmt in ['png', 'pdf']:
        fig.savefig(os.path.join(FIG_DIR, f'{filename}.{fmt}'))
    plt.close()
    print(f"Saved {filename}")

# =========================================================================
# Figure 1: |β × MeanFreq|
# =========================================================================
robust['AbsBetaMeanFreq'] = robust['Beta_x_MeanFreq'].abs()
plot_pair_bars(robust, 'Beta_x_MeanFreq', 'AbsBetaMeanFreq',
               '\u03b2 \u00d7 Mean Contact Frequency',
               'fig_beta_x_meanfreq')

# =========================================================================
# Figure 2: |β × MaxFreq|
# =========================================================================
robust['AbsBetaMaxFreq'] = robust['Beta_x_MaxFreq'].abs()
plot_pair_bars(robust, 'Beta_x_MaxFreq', 'AbsBetaMaxFreq',
               '\u03b2 \u00d7 Max Contact Frequency',
               'fig_beta_x_maxfreq')

# =========================================================================
# Figure 3: Temporal β × freq(t) for key pairs
# =========================================================================
print("\nPreparing temporal plot...")

# Select top pairs by |β × MeanFreq| plus any Y157 pairs in top 30
top_pairs_df = robust.nlargest(30, 'AbsBetaMeanFreq')
# Get the top 8 non-Y157 pairs + all Y157 pairs in top 30
non_y157 = [p for p in top_pairs_df['ResiduePair'] if not is_y157_pair(p)][:8]
y157_in_top = [p for p in top_pairs_df['ResiduePair'] if is_y157_pair(p)]

key_pairs = non_y157 + y157_in_top
beta_dict = dict(zip(robust['ResiduePair'], robust['Beta_final']))

# Filter temporal data to key pairs
temporal = freq_df[freq_df['ResiduePair'].isin(key_pairs)].copy()

# Convert frame to time (ps, 1 frame = 1 ps)
temporal['Time_ns'] = temporal['Frame'] / 1000.0

fig, ax = plt.subplots(figsize=(12, 6))

# Color palette for non-Y157
palette = ['#1565C0', '#0277BD', '#00838F', '#00695C', '#2E7D32',
           '#558B2F', '#827717', '#4527A0']

line_idx = 0
for pair in key_pairs:
    subset = temporal[temporal['ResiduePair'] == pair].sort_values('Frame')
    if len(subset) == 0:
        continue

    beta = beta_dict[pair]
    beta_x_freq = beta * subset['InteractionFrequency'].values

    # Smooth with moving average
    window = 100
    if len(beta_x_freq) > window:
        kernel = np.ones(window) / window
        beta_x_freq_smooth = np.convolve(beta_x_freq, kernel, mode='same')
    else:
        beta_x_freq_smooth = beta_x_freq

    t = subset['Time_ns'].values
    label = pair_label(pair)

    if is_y157_pair(pair):
        ax.plot(t, beta_x_freq_smooth, color='#D32F2F', lw=1.8,
                label=label, alpha=0.85, ls='--')
    else:
        color = palette[line_idx % len(palette)]
        ax.plot(t, beta_x_freq_smooth, color=color, lw=1.3,
                label=label, alpha=0.8)
        line_idx += 1

ax.set_xlabel('Time (ns)')
ax.set_ylabel('\u03b2 \u00d7 Contact Frequency')
ax.axhline(y=0, color='#333333', linewidth=0.8, alpha=0.5)
ax.legend(loc='best', ncol=2, framealpha=0.9, edgecolor='#cccccc', fontsize=7)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlim(0, 5.0)

for fmt in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'fig_temporal_beta_x_freq.{fmt}'))
plt.close()
print("Saved fig_temporal_beta_x_freq")

print("\nDone.")
