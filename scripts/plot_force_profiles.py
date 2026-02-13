#!/usr/bin/env python3
"""
Publication-quality figure: C1-87g7 Wild Type SMD Force Profiles (n=10).

Generates a 2-panel figure:
  Panel A: Individual force-vs-COM-distance curves (light colored lines)
           with a bold smoothed average line.
  Panel B: Average force profile +/- SEM shaded region.

Data: GROMACS pull force (.xvg) and COM distance (.xvg) files.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import ticker
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
import os

# ---------------------------------------------------------------------------
# 1. Paths
# ---------------------------------------------------------------------------
base_dir = "/home/anugraha/c1_WT"
rep_dirs = {
    1: os.path.join(base_dir, "pull"),
}
for i in range(2, 11):
    rep_dirs[i] = os.path.join(base_dir, f"pull_rep{i}")

out_path = "/home/anugraha/antibody_optimization/figures/fig1_force_profiles.png"

# ---------------------------------------------------------------------------
# 2. Helper: read GROMACS .xvg
# ---------------------------------------------------------------------------
def read_xvg(filepath):
    """Return (time, values) arrays from a GROMACS .xvg file."""
    t, v = [], []
    with open(filepath) as fh:
        for line in fh:
            if line.startswith(("#", "@")):
                continue
            parts = line.split()
            t.append(float(parts[0]))
            v.append(float(parts[1]))
    return np.array(t), np.array(v)

# ---------------------------------------------------------------------------
# 3. Load all replicates
# ---------------------------------------------------------------------------
forces_raw = []   # list of (com_distance, force) tuples
for rep_id in sorted(rep_dirs.keys()):
    d = rep_dirs[rep_id]
    _, force = read_xvg(os.path.join(d, "replica2_pullf.xvg"))
    _, com   = read_xvg(os.path.join(d, "replica2_pullx.xvg"))
    forces_raw.append((com, force))
    print(f"  Rep {rep_id:2d}: {len(com)} pts, COM range {com.min():.3f} - {com.max():.3f} nm, "
          f"Force range {force.min():.1f} - {force.max():.1f} kJ/mol/nm")

n_reps = len(forces_raw)
print(f"\nLoaded {n_reps} replicates.")

# ---------------------------------------------------------------------------
# 4. Build a common COM-distance grid and interpolate all replicates onto it
# ---------------------------------------------------------------------------
# Determine the common range (intersection of all replicate COM ranges).
com_min = max(com.min() for com, _ in forces_raw)
com_max = min(com.max() for com, _ in forces_raw)
print(f"Common COM range: {com_min:.4f} - {com_max:.4f} nm")

n_grid = 5000  # dense enough for smooth curves
com_grid = np.linspace(com_min, com_max, n_grid)

force_matrix = np.zeros((n_reps, n_grid))
for i, (com, force) in enumerate(forces_raw):
    # Sort by COM distance (should be monotonically increasing for SMD, but
    # small fluctuations can cause issues with interpolation).
    sort_idx = np.argsort(com)
    com_sorted = com[sort_idx]
    force_sorted = force[sort_idx]

    # Use a running average to pre-smooth before interpolation to handle
    # the case where COM is not strictly monotonic.
    # First, bin the data into the grid to handle duplicates.
    # Simple approach: use numpy digitize + binned means.
    bin_edges = np.linspace(com_min, com_max, n_grid + 1)
    bin_idx = np.digitize(com_sorted, bin_edges) - 1
    bin_idx = np.clip(bin_idx, 0, n_grid - 1)

    binned_force = np.full(n_grid, np.nan)
    for b in range(n_grid):
        mask = bin_idx == b
        if mask.any():
            binned_force[b] = force_sorted[mask].mean()

    # Fill any gaps via linear interpolation of the binned data.
    valid = ~np.isnan(binned_force)
    if valid.sum() < 100:
        # fallback: plain interp1d on sorted data
        f_interp = interp1d(com_sorted, force_sorted, bounds_error=False,
                            fill_value="extrapolate")
        force_matrix[i] = f_interp(com_grid)
    else:
        f_interp = interp1d(com_grid[valid], binned_force[valid],
                            bounds_error=False, fill_value="extrapolate")
        force_matrix[i] = f_interp(com_grid)

print("Interpolation onto common grid complete.")

# ---------------------------------------------------------------------------
# 5. Smoothing
# ---------------------------------------------------------------------------
# Savitzky-Golay: window must be odd
sg_window = 101  # ~101 grid points  (~0.5 nm window on 5000-point grid)
sg_order = 3

force_smooth = np.zeros_like(force_matrix)
for i in range(n_reps):
    force_smooth[i] = savgol_filter(force_matrix[i], sg_window, sg_order)

avg_raw = force_matrix.mean(axis=0)
avg_smooth = savgol_filter(avg_raw, sg_window, sg_order)

sem_raw = force_matrix.std(axis=0, ddof=1) / np.sqrt(n_reps)
sem_smooth = savgol_filter(sem_raw, sg_window, sg_order)

# ---------------------------------------------------------------------------
# 6. Plotting
# ---------------------------------------------------------------------------
# Publication style
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

# Colour palette for individual replicates (muted, distinguishable).
rep_colors = plt.cm.tab10(np.linspace(0, 1, 10))

fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True)

# --- Panel A: Individual curves + smoothed average ---
ax = axes[0]
for i in range(n_reps):
    ax.plot(com_grid, force_smooth[i], color=rep_colors[i], alpha=0.35,
            linewidth=0.8, label=f"Rep {i+1}" if i < 10 else None)

ax.plot(com_grid, avg_smooth, color="black", linewidth=2.5, label="Mean (smoothed)",
        zorder=10)

ax.set_xlabel("COM Distance (nm)")
ax.set_ylabel("Pull Force (kJ mol$^{-1}$ nm$^{-1}$)")
ax.text(0.03, 0.95, "(A)", transform=ax.transAxes, fontsize=16,
        fontweight="bold", va="top")

# Legend: two columns, outside bottom or inside
leg = ax.legend(loc="upper right", ncol=2, frameon=True, framealpha=0.9,
                edgecolor="gray", fontsize=8.5)
leg.get_frame().set_linewidth(0.5)

# --- Panel B: Average +/- SEM ---
ax2 = axes[1]
ax2.fill_between(com_grid, avg_smooth - sem_smooth, avg_smooth + sem_smooth,
                 color="#4C72B0", alpha=0.30, label="SEM")
ax2.plot(com_grid, avg_smooth, color="#4C72B0", linewidth=2.2, label="Mean")

ax2.set_xlabel("COM Distance (nm)")
# y-label already shared via sharey
ax2.text(0.03, 0.95, "(B)", transform=ax2.transAxes, fontsize=16,
         fontweight="bold", va="top")

leg2 = ax2.legend(loc="upper right", frameon=True, framealpha=0.9,
                  edgecolor="gray")
leg2.get_frame().set_linewidth(0.5)

# Suptitle
fig.suptitle("C1-87g7 Wild Type \u2014 SMD Force Profiles ($n$ = 10)",
             fontsize=15, fontweight="bold", y=1.02)

plt.tight_layout()
fig.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
print(f"\nFigure saved to: {out_path}")
plt.close()
