#!/usr/bin/env python3
"""
Publication-quality figure: C1-87g7 Wild Type Interaction Energy vs COM Distance.

Extracts Coul-SR:chauA-rest and LJ-SR:chauA-rest from energy rerun .edr files
for all 10 pulling replicates, computes total interaction energy (Coul + LJ),
and plots it against COM distance.

2-panel figure:
  Panel A: Individual interaction energy curves + bold mean
  Panel B: Mean ± SEM shaded region

Requires: energy rerun .edr files to already exist.
"""

import numpy as np
import subprocess
import os
import sys
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d

# ---------------------------------------------------------------------------
# 1. Paths
# ---------------------------------------------------------------------------
base_dir = "/home/anugraha/c1_WT"
rep_dirs = {1: os.path.join(base_dir, "pull")}
for i in range(2, 11):
    rep_dirs[i] = os.path.join(base_dir, f"pull_rep{i}")

out_path = "/home/anugraha/antibody_optimization/figures/fig2_interaction_energy.png"
os.makedirs(os.path.dirname(out_path), exist_ok=True)

# ---------------------------------------------------------------------------
# 2. Helper: read GROMACS .xvg
# ---------------------------------------------------------------------------
def read_xvg(filepath):
    """Return (time, values) arrays from a GROMACS .xvg file."""
    t, vals = [], []
    with open(filepath) as fh:
        for line in fh:
            if line.startswith(("#", "@")):
                continue
            parts = line.split()
            t.append(float(parts[0]))
            # May have multiple columns
            vals.append([float(x) for x in parts[1:]])
    return np.array(t), np.array(vals)


def extract_energy(edr_file, output_xvg):
    """Extract Coul-SR:chauA-rest and LJ-SR:chauA-rest from .edr file."""
    # Use gmx energy with echo to select terms 22 and 23
    cmd = f"echo '22 23 0' | gmx_mpi energy -f {edr_file} -o {output_xvg}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True,
                          env={**os.environ, 'GMX_MAXBACKUP': '-1'})
    if result.returncode != 0:
        print(f"  WARNING: gmx energy failed for {edr_file}")
        print(f"  stderr: {result.stderr[:500]}")
        return False
    return True


# ---------------------------------------------------------------------------
# 3. Extract energies and load data
# ---------------------------------------------------------------------------
data = []  # list of (com_distance, interaction_energy)
for rep_id in sorted(rep_dirs.keys()):
    d = rep_dirs[rep_id]
    edr_file = os.path.join(d, "energy_rerun.edr")

    if not os.path.exists(edr_file):
        print(f"  Rep {rep_id}: energy_rerun.edr not found, skipping")
        continue

    # Check if edr file is complete (not still being written)
    xvg_file = os.path.join(d, "energy_chauA_rest.xvg")

    if not os.path.exists(xvg_file):
        print(f"  Rep {rep_id}: extracting energies...")
        if not extract_energy(edr_file, xvg_file):
            continue
    else:
        print(f"  Rep {rep_id}: energy xvg already exists")

    # Load energy data (columns: Coul-SR:chauA-rest, LJ-SR:chauA-rest)
    t_en, en_vals = read_xvg(xvg_file)
    coul_sr = en_vals[:, 0]  # Coul-SR:chauA-rest
    lj_sr = en_vals[:, 1]    # LJ-SR:chauA-rest
    interaction_energy = coul_sr + lj_sr  # Total interaction energy

    # Load COM distance
    com_file = os.path.join(d, "replica2_pullx.xvg")
    _, com = read_xvg(com_file)
    com_dist = com[:, 0]  # First column is COM distance

    # Match lengths (energy rerun might have different frame count)
    min_len = min(len(t_en), len(com_dist))
    t_en = t_en[:min_len]
    interaction_energy = interaction_energy[:min_len]
    com_dist = com_dist[:min_len]

    data.append((com_dist, interaction_energy))
    print(f"  Rep {rep_id}: {min_len} frames, COM {com_dist.min():.3f}-{com_dist.max():.3f} nm, "
          f"E_int {interaction_energy.min():.1f} to {interaction_energy.max():.1f} kJ/mol")

n_reps = len(data)
print(f"\nLoaded {n_reps} replicates.")

if n_reps == 0:
    print("ERROR: No data loaded. Energy reruns may not be complete yet.")
    sys.exit(1)

# ---------------------------------------------------------------------------
# 4. Interpolate onto common COM grid
# ---------------------------------------------------------------------------
com_min = max(com.min() for com, _ in data)
com_max = min(com.max() for com, _ in data)
print(f"Common COM range: {com_min:.4f} - {com_max:.4f} nm")

n_grid = 5000
com_grid = np.linspace(com_min, com_max, n_grid)

energy_matrix = np.zeros((n_reps, n_grid))
for i, (com, energy) in enumerate(data):
    sort_idx = np.argsort(com)
    com_sorted = com[sort_idx]
    energy_sorted = energy[sort_idx]

    # Bin data
    bin_edges = np.linspace(com_min, com_max, n_grid + 1)
    bin_idx = np.digitize(com_sorted, bin_edges) - 1
    bin_idx = np.clip(bin_idx, 0, n_grid - 1)

    binned = np.full(n_grid, np.nan)
    for b in range(n_grid):
        mask = bin_idx == b
        if mask.any():
            binned[b] = energy_sorted[mask].mean()

    valid = ~np.isnan(binned)
    if valid.sum() < 100:
        f_interp = interp1d(com_sorted, energy_sorted, bounds_error=False,
                           fill_value="extrapolate")
        energy_matrix[i] = f_interp(com_grid)
    else:
        f_interp = interp1d(com_grid[valid], binned[valid],
                           bounds_error=False, fill_value="extrapolate")
        energy_matrix[i] = f_interp(com_grid)

# ---------------------------------------------------------------------------
# 5. Smoothing
# ---------------------------------------------------------------------------
sg_window = 101
sg_order = 3

energy_smooth = np.zeros_like(energy_matrix)
for i in range(n_reps):
    energy_smooth[i] = savgol_filter(energy_matrix[i], sg_window, sg_order)

avg_raw = energy_matrix.mean(axis=0)
avg_smooth = savgol_filter(avg_raw, sg_window, sg_order)

sem_raw = energy_matrix.std(axis=0, ddof=1) / np.sqrt(n_reps)
sem_smooth = savgol_filter(sem_raw, sg_window, sg_order)

# ---------------------------------------------------------------------------
# 6. Plotting
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 12,
    "axes.labelsize": 14,
    "axes.titlesize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 10,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.linewidth": 1.2,
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
    "xtick.major.size": 5,
    "ytick.major.size": 5,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
})

rep_colors = plt.cm.tab10(np.linspace(0, 1, 10))

fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True)

# --- Panel A: Individual curves + smoothed average ---
ax = axes[0]
for i in range(n_reps):
    ax.plot(com_grid, energy_smooth[i], color=rep_colors[i], alpha=0.35,
            linewidth=0.8, label=f"Rep {i+1}" if i < 10 else None)

ax.plot(com_grid, avg_smooth, color="black", linewidth=2.5, label="Mean (smoothed)",
        zorder=10)

ax.set_xlabel("COM Distance (nm)")
ax.set_ylabel("Interaction Energy (kJ/mol)")
ax.text(0.03, 0.95, "(A)", transform=ax.transAxes, fontsize=16,
        fontweight="bold", va="top")

leg = ax.legend(loc="lower right", ncol=2, frameon=True, framealpha=0.9,
                edgecolor="gray", fontsize=8.5)
leg.get_frame().set_linewidth(0.5)

# --- Panel B: Average +/- SEM ---
ax2 = axes[1]
ax2.fill_between(com_grid, avg_smooth - sem_smooth, avg_smooth + sem_smooth,
                 color="#C44E52", alpha=0.30, label="SEM")
ax2.plot(com_grid, avg_smooth, color="#C44E52", linewidth=2.2, label="Mean")

ax2.set_xlabel("COM Distance (nm)")
ax2.text(0.03, 0.95, "(B)", transform=ax2.transAxes, fontsize=16,
         fontweight="bold", va="top")

leg2 = ax2.legend(loc="lower right", frameon=True, framealpha=0.9,
                  edgecolor="gray")
leg2.get_frame().set_linewidth(0.5)

fig.suptitle("C1-87g7 Wild Type — Interaction Energy During SMD ($n$ = %d)" % n_reps,
             fontsize=15, fontweight="bold", y=1.02)

plt.tight_layout()
fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
print(f"\nFigure saved to: {out_path}")
plt.close()
