#!/usr/bin/env python3
"""Compare mutant pulling outputs against C1 baseline.

Expected mutant file locations:
  /home/anugraha/c1_<MUT>/pull/replica6_Ax6.xvg
  /home/anugraha/c1_<MUT>/pull/replica6_Af6.xvg
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASELINE_REPS = [
    (
        "c1_rep0",
        Path(
            "/home/anugraha/c1_87g7/rot_120_corrected/bigger_box/replica_0ns_pos/replica6_Ax6.xvg"
        ),
        Path(
            "/home/anugraha/c1_87g7/rot_120_corrected/bigger_box/replica_0ns_pos/replica6_Af6.xvg"
        ),
    ),
    (
        "c1_rep80",
        Path(
            "/home/anugraha/c1_87g7/rot_120_corrected/bigger_box/md2/eq80/replica_80ns_pos/replica6_Ax6.xvg"
        ),
        Path(
            "/home/anugraha/c1_87g7/rot_120_corrected/bigger_box/md2/eq80/replica_80ns_pos/replica6_Af6.xvg"
        ),
    ),
    (
        "c1_rep100",
        Path(
            "/home/anugraha/c1_87g7/rot_120_corrected/bigger_box/md2/eq100/replica_100ns/replica6_Ax6.xvg"
        ),
        Path(
            "/home/anugraha/c1_87g7/rot_120_corrected/bigger_box/md2/eq100/replica_100ns/replica6_Af6.xvg"
        ),
    ),
]


def read_xvg(path: Path) -> np.ndarray:
    vals = []
    with path.open() as f:
        for line in f:
            if line.startswith(("@", "#")):
                continue
            sp = line.split()
            if len(sp) < 2:
                continue
            try:
                vals.append(float(sp[1]))
            except ValueError:
                continue
    return np.asarray(vals, dtype=float)


def align_xy(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    n = min(len(x), len(y))
    if n == 0:
        return np.array([]), np.array([])
    x = x[:n] - x[0]
    y = y[:n]
    return x, y


def metrics(x: np.ndarray, y: np.ndarray) -> Dict[str, float]:
    if len(x) == 0:
        return {"peak_force": np.nan, "auc_force_dist": np.nan}
    return {
        "peak_force": float(np.nanmax(y)),
        "auc_force_dist": float(np.trapz(y, x)),
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--mutations",
        nargs="+",
        required=True,
        help="Mutation tags without prefix, e.g. Y405N Y406W",
    )
    p.add_argument(
        "--outdir",
        type=Path,
        default=Path("/home/anugraha/antibody_optimization/results"),
    )
    p.add_argument("--ylim", type=float, default=1200.0)
    args = p.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    rows = []

    plt.figure(figsize=(10, 6))

    # Plot baseline replicates
    c1_metrics = []
    for label, x_path, y_path in BASELINE_REPS:
        if not (x_path.exists() and y_path.exists()):
            continue
        x = read_xvg(x_path)
        y = read_xvg(y_path)
        x, y = align_xy(x, y)
        m = metrics(x, y)
        c1_metrics.append(m)
        rows.append(
            {
                "system": "c1_87g7",
                "replicate": label,
                "peak_force": m["peak_force"],
                "auc_force_dist": m["auc_force_dist"],
                "x_file": str(x_path),
                "y_file": str(y_path),
            }
        )
        plt.plot(x, y, color="purple", alpha=0.35, linewidth=0.7, label=label)

    # Plot mutants
    cmap = plt.get_cmap("tab10")
    for i, mut in enumerate(args.mutations):
        folder = Path(f"/home/anugraha/c1_{mut}/pull")
        x_path = folder / "replica6_Ax6.xvg"
        y_path = folder / "replica6_Af6.xvg"
        if not (x_path.exists() and y_path.exists()):
            rows.append(
                {
                    "system": f"c1_{mut}",
                    "replicate": "pull",
                    "peak_force": np.nan,
                    "auc_force_dist": np.nan,
                    "x_file": str(x_path),
                    "y_file": str(y_path),
                }
            )
            continue
        x = read_xvg(x_path)
        y = read_xvg(y_path)
        x, y = align_xy(x, y)
        m = metrics(x, y)
        rows.append(
            {
                "system": f"c1_{mut}",
                "replicate": "pull",
                "peak_force": m["peak_force"],
                "auc_force_dist": m["auc_force_dist"],
                "x_file": str(x_path),
                "y_file": str(y_path),
            }
        )
        plt.plot(x, y, color=cmap(i % 10), alpha=0.85, linewidth=1.2, label=f"c1_{mut}")

    # Add baseline mean row
    if c1_metrics:
        rows.append(
            {
                "system": "c1_87g7",
                "replicate": "mean",
                "peak_force": float(np.mean([x["peak_force"] for x in c1_metrics])),
                "auc_force_dist": float(np.mean([x["auc_force_dist"] for x in c1_metrics])),
                "x_file": "",
                "y_file": "",
            }
        )

    out_csv = args.outdir / "comparison_summary.csv"
    pd.DataFrame(rows).to_csv(out_csv, index=False)

    plt.xlabel("Normalized COM Distance (nm)")
    plt.ylabel("Force (kJ/mol/nm^2)")
    plt.ylim(0, args.ylim)
    plt.legend(loc="upper left", fontsize=8)
    plt.tight_layout()

    out_png = args.outdir / "c1_vs_mutants_overlay.png"
    plt.savefig(out_png, dpi=300)
    print(f"Wrote: {out_csv}")
    print(f"Wrote: {out_png}")


if __name__ == "__main__":
    main()
