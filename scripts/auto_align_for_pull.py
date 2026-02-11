#!/usr/bin/env python3
"""Auto-align protein interface for SMD pulling simulations.

Reads a vacuum .gro from ACPYPE (before solvation), identifies
antigen/antibody groups, rotates so the interface normal is parallel
to the pull axis (Y), builds an asymmetric box sized for the pull
displacement, and writes the corrected .gro.

Standalone script -- all helper functions are self-contained for
cluster portability (no imports from pull_cli_view.py).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Protein residue names (matches build_pull_index.py)
# ---------------------------------------------------------------------------
PROT_RESNAMES = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "CYX",
}

# ---------------------------------------------------------------------------
# GRO parser / writer
# ---------------------------------------------------------------------------

def parse_gro(path: Path) -> Tuple[Dict[str, np.ndarray], np.ndarray, str]:
    """Return (atoms dict, box[3], title string)."""
    lines = path.read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"Invalid .gro file: {path}")
    title = lines[0]
    n_atoms = int(lines[1].strip())
    atom_lines = lines[2 : 2 + n_atoms]
    box_line = lines[2 + n_atoms].split()

    if len(box_line) not in (3, 9):
        raise ValueError(f"Unsupported box format in {path}: {box_line}")
    box = np.array([float(box_line[0]), float(box_line[1]), float(box_line[2])], dtype=float)

    resid = np.zeros(n_atoms, dtype=int)
    resname: List[str] = []
    atomname: List[str] = []
    atomnr = np.zeros(n_atoms, dtype=int)
    xyz = np.zeros((n_atoms, 3), dtype=float)

    for i, line in enumerate(atom_lines):
        resid[i] = int(line[0:5])
        resname.append(line[5:10].strip())
        atomname.append(line[10:15].strip())
        atomnr[i] = int(line[15:20])
        xyz[i, 0] = float(line[20:28])
        xyz[i, 1] = float(line[28:36])
        xyz[i, 2] = float(line[36:44])

    atoms = {
        "resid": resid,
        "resname": np.array(resname, dtype=object),
        "atomname": np.array(atomname, dtype=object),
        "atomnr": atomnr,
        "xyz": xyz,
    }
    return atoms, box, title


def write_gro(path: Path, atoms: Dict[str, np.ndarray], box: np.ndarray, title: str) -> None:
    xyz = atoms["xyz"]
    resid = atoms["resid"]
    resname = atoms["resname"]
    atomname = atoms["atomname"]
    atomnr = atoms["atomnr"]
    with path.open("w") as f:
        f.write(title[:80] + "\n")
        f.write(f"{len(xyz):5d}\n")
        for i in range(len(xyz)):
            f.write(
                f"{int(resid[i]) % 100000:5d}"
                f"{str(resname[i])[:5]:>5s}"
                f"{str(atomname[i])[:5]:>5s}"
                f"{int(atomnr[i]) % 100000:5d}"
                f"{float(xyz[i, 0]):8.3f}"
                f"{float(xyz[i, 1]):8.3f}"
                f"{float(xyz[i, 2]):8.3f}\n"
            )
        f.write(f"{box[0]:10.5f}{box[1]:10.5f}{box[2]:10.5f}\n")


# ---------------------------------------------------------------------------
# Linear-algebra helpers
# ---------------------------------------------------------------------------

def _unit(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    if n < 1e-12:
        raise ValueError("Cannot normalize near-zero vector")
    return v / n


def _skew(v: np.ndarray) -> np.ndarray:
    return np.array(
        [
            [0.0, -v[2], v[1]],
            [v[2], 0.0, -v[0]],
            [-v[1], v[0], 0.0],
        ],
        dtype=float,
    )


def rotation_matrix_from_vectors(src: np.ndarray, dst: np.ndarray) -> np.ndarray:
    """Rodrigues rotation from *src* direction to *dst* direction."""
    a = _unit(src)
    b = _unit(dst)
    c = float(np.clip(np.dot(a, b), -1.0, 1.0))
    if c > 1.0 - 1e-10:
        return np.eye(3, dtype=float)
    if c < -1.0 + 1e-10:
        aux = np.array([1.0, 0.0, 0.0], dtype=float)
        if abs(np.dot(a, aux)) > 0.9:
            aux = np.array([0.0, 1.0, 0.0], dtype=float)
        u = _unit(np.cross(a, aux))
        return -np.eye(3, dtype=float) + 2.0 * np.outer(u, u)
    v = np.cross(a, b)
    s = np.linalg.norm(v)
    k = _skew(v)
    return np.eye(3, dtype=float) + k + (k @ k) * ((1.0 - c) / (s * s))


def rotate_about_point(coords: np.ndarray, R: np.ndarray, origin: np.ndarray) -> np.ndarray:
    return (coords - origin) @ R.T + origin


# ---------------------------------------------------------------------------
# Interface contact detection (chunked to keep memory bounded)
# ---------------------------------------------------------------------------

def interface_contact_masks(
    g1: np.ndarray,
    g2: np.ndarray,
    cutoff_nm: float,
    chunk_size: int = 2000,
) -> Tuple[np.ndarray, np.ndarray]:
    cutoff2 = float(cutoff_nm * cutoff_nm)
    n1 = len(g1)
    n2 = len(g2)
    mask1 = np.zeros(n1, dtype=bool)
    mask2 = np.zeros(n2, dtype=bool)
    if n1 == 0 or n2 == 0:
        return mask1, mask2
    b2 = np.einsum("ij,ij->i", g2, g2)
    bt = g2.T
    for i0 in range(0, n1, chunk_size):
        i1 = min(i0 + chunk_size, n1)
        a = g1[i0:i1]
        a2 = np.einsum("ij,ij->i", a, a)
        dist2 = a2[:, None] + b2[None, :] - 2.0 * (a @ bt)
        close = dist2 <= cutoff2
        mask1[i0:i1] = np.any(close, axis=1)
        mask2 |= np.any(close, axis=0)
    return mask1, mask2


def interface_normal_from_contacts(
    g1: np.ndarray,
    g2: np.ndarray,
    cutoff_nm: float,
) -> Tuple[np.ndarray, np.ndarray, int, int]:
    """SVD-based interface normal from contact atoms between two groups."""
    mask1, mask2 = interface_contact_masks(g1, g2, cutoff_nm=cutoff_nm)
    n1 = int(np.count_nonzero(mask1))
    n2 = int(np.count_nonzero(mask2))
    if n1 < 3 or n2 < 3:
        raise ValueError(f"Not enough interface contact atoms (group1={n1}, group2={n2})")
    pts = np.vstack([g1[mask1], g2[mask2]])
    center = np.mean(pts, axis=0)
    centered = pts - center
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    normal = vh[-1]
    return normal, center, n1, n2


# ---------------------------------------------------------------------------
# Group identification (by residue order, same logic as build_pull_index.py)
# ---------------------------------------------------------------------------

def identify_groups(
    atoms: Dict[str, np.ndarray], antigen_res_count: int
) -> Tuple[np.ndarray, np.ndarray]:
    """Split protein into chauA (antigen) and rest (antibody) by residue order.

    Returns (mask_chauA, mask_rest) -- boolean masks over all atoms.
    """
    resid = atoms["resid"]
    resname = atoms["resname"]
    n = len(resid)

    residue_order: List[Tuple[int, str]] = []
    residue_atom_indices: Dict[Tuple[int, str], List[int]] = {}
    last_key = None
    for i in range(n):
        rn = str(resname[i])
        if rn not in PROT_RESNAMES:
            continue
        key = (int(resid[i]), rn)
        if key != last_key:
            residue_order.append(key)
            residue_atom_indices[key] = []
            last_key = key
        residue_atom_indices[key].append(i)

    if len(residue_order) < antigen_res_count:
        raise ValueError(
            f"Protein residue count ({len(residue_order)}) < antigen count ({antigen_res_count})"
        )

    chau_res = residue_order[:antigen_res_count]
    rest_res = residue_order[antigen_res_count:]

    mask_chau = np.zeros(n, dtype=bool)
    mask_rest = np.zeros(n, dtype=bool)
    for key in chau_res:
        for idx in residue_atom_indices[key]:
            mask_chau[idx] = True
    for key in rest_res:
        for idx in residue_atom_indices[key]:
            mask_rest[idx] = True

    return mask_chau, mask_rest


# ---------------------------------------------------------------------------
# Alignment computation
# ---------------------------------------------------------------------------

def compute_alignment_rotation(
    g1_coords: np.ndarray,
    g2_coords: np.ndarray,
    pull_axis_vec: np.ndarray,
    cutoff_nm: float,
    fallback_to_com: bool,
) -> Tuple[np.ndarray, np.ndarray, str, float, int, int]:
    """Compute rotation matrix to align interface normal (or COM vector) to pull axis.

    Returns: (R, origin, method, angle_before_deg, n_contact_g1, n_contact_g2)
    """
    com1 = np.mean(g1_coords, axis=0)
    com2 = np.mean(g2_coords, axis=0)

    method = "interface-normal"
    n_contact_g1 = 0
    n_contact_g2 = 0

    try:
        normal, iface_center, n_contact_g1, n_contact_g2 = interface_normal_from_contacts(
            g1_coords, g2_coords, cutoff_nm=cutoff_nm
        )
        rot_vec = normal
        origin = iface_center
    except ValueError as exc:
        if not fallback_to_com:
            raise
        print(f"[WARN] Interface detection failed ({exc}), falling back to COM-COM vector")
        rot_vec = com2 - com1
        if np.linalg.norm(rot_vec) < 1e-8:
            raise ValueError("COM separation is ~0; cannot align to pull axis.")
        origin = 0.5 * (com1 + com2)
        method = "com-vector (fallback)"

    # Orient so that the vector points in the same direction as pull axis
    target = pull_axis_vec if np.dot(rot_vec, pull_axis_vec) >= 0.0 else -pull_axis_vec
    angle_before_deg = float(
        np.degrees(np.arccos(np.clip(abs(np.dot(_unit(rot_vec), _unit(pull_axis_vec))), 0.0, 1.0)))
    )

    R = rotation_matrix_from_vectors(rot_vec, target)
    return R, origin, method, angle_before_deg, n_contact_g1, n_contact_g2


# ---------------------------------------------------------------------------
# Box sizing
# ---------------------------------------------------------------------------

def compute_box(
    coords: np.ndarray,
    pull_axis_idx: int,
    pull_distance_nm: float,
    margin_pull_nm: float,
    margin_perp_nm: float,
    g1_coords: np.ndarray,
    g2_coords: np.ndarray,
) -> np.ndarray:
    """PBC-safe box for pulling simulations.

    Two constraints on the pull axis:
      1. Water padding:  box >= span + 2*margin
      2. Minimum image:  box/2 > COM_sep + pull_dist  (with safety margin)
         =>  box > 2*(COM_sep + pull_dist) + 2*margin
    We take the maximum of both.

    Perpendicular axes just need water padding: span + 2*margin.
    """
    span = np.ptp(coords, axis=0)
    com_sep = abs(np.mean(g1_coords[:, pull_axis_idx]) - np.mean(g2_coords[:, pull_axis_idx]))

    box = np.zeros(3, dtype=float)
    for ax in range(3):
        if ax == pull_axis_idx:
            padding_box = span[ax] + 2.0 * margin_pull_nm
            pbc_box = 2.0 * (com_sep + pull_distance_nm) + 2.0 * margin_pull_nm
            box[ax] = max(padding_box, pbc_box)
        else:
            box[ax] = span[ax] + 2.0 * margin_perp_nm

    # Report which constraint was binding
    if pbc_box > padding_box:
        print(f"[INFO] Pull-axis box set by PBC constraint: "
              f"2*({com_sep:.2f} + {pull_distance_nm:.2f}) + 2*{margin_pull_nm:.1f} = {pbc_box:.2f} nm")
    else:
        print(f"[INFO] Pull-axis box set by padding constraint: "
              f"span {span[pull_axis_idx]:.2f} + 2*{margin_pull_nm:.1f} = {padding_box:.2f} nm")
    print(f"[INFO] PBC check: box/2 = {box[pull_axis_idx]/2:.2f} nm > "
          f"COM_sep + pull = {com_sep:.2f} + {pull_distance_nm:.2f} = {com_sep + pull_distance_nm:.2f} nm  OK")

    return box


def center_protein(coords: np.ndarray, box: np.ndarray) -> np.ndarray:
    """Center geometric center at box/2."""
    geo_center = np.mean(coords, axis=0)
    shift = (box / 2.0) - geo_center
    return coords + shift


def position_for_pull(
    coords: np.ndarray,
    box: np.ndarray,
    pull_axis_idx: int,
    larger_mask: np.ndarray,
    margin_nm: float,
) -> np.ndarray:
    """Center on X/Z, offset along pull axis.

    Places the restrained (larger) group near -wall with just *margin_nm*
    of padding.  The PBC-safe box size is already set by compute_box(),
    so all the extra room ends up at +pull-axis where the smaller group
    will be pulled into.
    """
    out = coords.copy()
    # Perpendicular axes: center normally
    for ax in range(3):
        if ax == pull_axis_idx:
            continue
        geo_center = np.mean(out[:, ax])
        out[:, ax] += (box[ax] / 2.0) - geo_center
    # Pull axis: place larger (restrained) group's min at margin from -wall
    larger_min = np.min(out[larger_mask, pull_axis_idx])
    out[:, pull_axis_idx] += margin_nm - larger_min
    return out


# ---------------------------------------------------------------------------
# MDP parser for pull distance
# ---------------------------------------------------------------------------

def parse_pull_distance_from_mdp(mdp_path: Path) -> float:
    """Auto-detect pull distance from umb2.mdp (rate * nsteps * dt)."""
    params: Dict[str, str] = {}
    for line in mdp_path.read_text().splitlines():
        stripped = line.split(";")[0].strip()
        if "=" not in stripped:
            continue
        key, val = stripped.split("=", 1)
        params[key.strip()] = val.strip()

    dt = float(params.get("dt", "0.002"))
    nsteps = int(params.get("nsteps", "0"))
    rate = float(params.get("pull_coord1_rate", "0.0"))

    total_time_ps = nsteps * dt
    distance_nm = abs(rate) * total_time_ps
    return distance_nm


# ---------------------------------------------------------------------------
# Visual report
# ---------------------------------------------------------------------------

def _projection_subplot(
    ax: plt.Axes,
    coords_before: np.ndarray,
    coords_after: np.ndarray,
    g1_mask: np.ndarray,
    g2_mask: np.ndarray,
    box_before: np.ndarray,
    box_after: np.ndarray,
    ax_i: int,
    ax_j: int,
    label_x: str,
    label_y: str,
    is_before: bool,
) -> None:
    coords = coords_before if is_before else coords_after
    box = box_before if is_before else box_after
    g1 = coords[g1_mask]
    g2 = coords[g2_mask]
    other = coords[~(g1_mask | g2_mask)]

    if len(other) > 0:
        ax.scatter(other[:, ax_i], other[:, ax_j], s=0.5, alpha=0.15, c="#9ca3af")
    ax.scatter(g1[:, ax_i], g1[:, ax_j], s=2, alpha=0.7, c="#0ea5e9", label="chauA")
    ax.scatter(g2[:, ax_i], g2[:, ax_j], s=2, alpha=0.7, c="#ef4444", label="rest")

    # Draw box edges
    ax.plot(
        [0, box[ax_i], box[ax_i], 0, 0],
        [0, 0, box[ax_j], box[ax_j], 0],
        "k-", linewidth=0.8, alpha=0.5,
    )
    ax.set_xlabel(f"{label_x} (nm)", fontsize=8)
    ax.set_ylabel(f"{label_y} (nm)", fontsize=8)
    ax.set_aspect("equal", adjustable="box")
    ax.tick_params(labelsize=7)


def generate_report_png(
    out_png: Path,
    coords_before: np.ndarray,
    coords_after: np.ndarray,
    g1_mask: np.ndarray,
    g2_mask: np.ndarray,
    box_before: np.ndarray,
    box_after: np.ndarray,
    method: str,
    angle_before: float,
    angle_after: float,
) -> None:
    """2x3 grid: before/after x XY/YZ/XZ projections."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    planes = [(0, 1, "X", "Y"), (1, 2, "Y", "Z"), (0, 2, "X", "Z")]

    for col, (ai, aj, lx, ly) in enumerate(planes):
        _projection_subplot(
            axes[0, col], coords_before, coords_after,
            g1_mask, g2_mask, box_before, box_after,
            ai, aj, lx, ly, is_before=True,
        )
        axes[0, col].set_title(f"BEFORE  {lx}-{ly}", fontsize=9)

        _projection_subplot(
            axes[1, col], coords_before, coords_after,
            g1_mask, g2_mask, box_before, box_after,
            ai, aj, lx, ly, is_before=False,
        )
        axes[1, col].set_title(f"AFTER  {lx}-{ly}", fontsize=9)

    # Legend on first subplot only
    axes[0, 0].legend(loc="upper right", fontsize=7, markerscale=2)

    fig.suptitle(
        f"Auto-align for pull  |  method: {method}  |  "
        f"angle: {angle_before:.1f}° → {angle_after:.1f}°",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"[INFO] Wrote report: {out_png}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gro", type=Path, required=True, help="Input vacuum .gro")
    p.add_argument("--out", type=Path, required=True, help="Output aligned .gro")
    p.add_argument("--antigen-res-count", type=int, required=True,
                    help="Number of antigen residues (first N protein residues)")
    p.add_argument("--pull-dim", choices=["X", "Y", "Z"], default="Y",
                    help="Pull axis dimension (default Y)")
    p.add_argument("--mdp", type=Path, default=None,
                    help="umb2.mdp to auto-detect pull distance")
    p.add_argument("--pull-distance-nm", type=float, default=None,
                    help="Override pull distance in nm (used if --mdp not given)")
    p.add_argument("--margin-pull-nm", type=float, default=2.0,
                    help="Box margin along pull axis (default 2.0)")
    p.add_argument("--margin-perp-nm", type=float, default=1.5,
                    help="Box margin on perpendicular axes (default 1.5)")
    p.add_argument("--interface-cutoff-nm", type=float, default=0.6,
                    help="Contact cutoff for interface detection (default 0.6)")
    p.add_argument("--report-png", type=Path, default=None,
                    help="Output before/after visual report PNG")
    p.add_argument("--fallback-to-com", action="store_true",
                    help="Fall back to COM-COM vector if interface detection fails")
    p.add_argument("--skip-threshold-deg", type=float, default=5.0,
                    help="Skip rotation if already aligned within this angle (default 5)")
    p.add_argument("--orient-smaller-positive", action="store_true", default=True,
                    help="Flip so smaller group is at +pull-axis (pulled away). On by default.")
    p.add_argument("--no-orient-smaller-positive", dest="orient_smaller_positive",
                    action="store_false",
                    help="Disable automatic orientation of smaller group to +pull-axis.")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    pull_axis_idx = {"X": 0, "Y": 1, "Z": 2}[args.pull_dim]
    pull_axis_vec = np.zeros(3, dtype=float)
    pull_axis_vec[pull_axis_idx] = 1.0

    # --- Determine pull distance ---
    if args.mdp is not None:
        pull_distance_nm = parse_pull_distance_from_mdp(args.mdp)
        print(f"[INFO] Pull distance from {args.mdp}: {pull_distance_nm:.2f} nm")
    elif args.pull_distance_nm is not None:
        pull_distance_nm = args.pull_distance_nm
    else:
        pull_distance_nm = 5.0
        print(f"[WARN] No --mdp or --pull-distance-nm given; using default {pull_distance_nm} nm")

    # --- Parse input ---
    atoms, box_orig, title_orig = parse_gro(args.gro)
    xyz = atoms["xyz"].copy()
    n_atoms = len(xyz)
    print(f"[INFO] Read {n_atoms} atoms from {args.gro}")
    print(f"[INFO] Original box: {box_orig[0]:.3f} x {box_orig[1]:.3f} x {box_orig[2]:.3f} nm")

    # --- Identify groups ---
    mask_chau, mask_rest = identify_groups(atoms, args.antigen_res_count)
    n_chau = int(np.count_nonzero(mask_chau))
    n_rest = int(np.count_nonzero(mask_rest))
    n_other = n_atoms - n_chau - n_rest
    print(f"[INFO] Groups: chauA={n_chau} atoms, rest={n_rest} atoms, other={n_other} atoms")

    g1_coords = xyz[mask_chau]
    g2_coords = xyz[mask_rest]

    # Save before-state for report
    xyz_before = xyz.copy()

    # --- Compute alignment ---
    R, origin, method, angle_before_deg, nc1, nc2 = compute_alignment_rotation(
        g1_coords, g2_coords, pull_axis_vec,
        cutoff_nm=args.interface_cutoff_nm,
        fallback_to_com=args.fallback_to_com,
    )

    print(f"[INFO] Alignment method: {method}")
    if "interface" in method:
        print(f"[INFO] Interface contacts: chauA={nc1}, rest={nc2}")
    print(f"[INFO] Angle before: {angle_before_deg:.2f}°")

    # --- Apply rotation (or skip if already aligned) ---
    rotated = False
    if angle_before_deg > args.skip_threshold_deg:
        xyz = rotate_about_point(xyz, R, origin)
        rotated = True
        print(f"[INFO] Rotation applied")
    else:
        print(f"[INFO] Already aligned within {args.skip_threshold_deg}° -- skipping rotation")

    # Measure angle after
    g1_after = xyz[mask_chau]
    g2_after = xyz[mask_rest]
    try:
        normal_after, _, _, _ = interface_normal_from_contacts(
            g1_after, g2_after, cutoff_nm=args.interface_cutoff_nm
        )
        angle_after_deg = float(
            np.degrees(np.arccos(np.clip(abs(np.dot(_unit(normal_after), pull_axis_vec)), 0.0, 1.0)))
        )
    except ValueError:
        com_vec = np.mean(g2_after, axis=0) - np.mean(g1_after, axis=0)
        angle_after_deg = float(
            np.degrees(np.arccos(np.clip(abs(np.dot(_unit(com_vec), pull_axis_vec)), 0.0, 1.0)))
        )
    print(f"[INFO] Angle after: {angle_after_deg:.2f}°")

    # --- Second rotation: align COM-COM vector to pull axis ---
    # After interface alignment, the COM vector may still have off-axis components.
    # With dim = N Y N pulling, the COM vector should be along Y for clean separation.
    g1_now = xyz[mask_chau]
    g2_now = xyz[mask_rest]
    com_vec_now = np.mean(g1_now, axis=0) - np.mean(g2_now, axis=0)
    com_angle_before = float(
        np.degrees(np.arccos(np.clip(abs(np.dot(_unit(com_vec_now), pull_axis_vec)), 0.0, 1.0)))
    )
    print(f"[INFO] COM vector angle from {args.pull_dim}: {com_angle_before:.2f}°")
    if com_angle_before > args.skip_threshold_deg:
        com_target = pull_axis_vec if np.dot(com_vec_now, pull_axis_vec) >= 0.0 else -pull_axis_vec
        R_com = rotation_matrix_from_vectors(com_vec_now, com_target)
        com_origin = 0.5 * (np.mean(g1_now, axis=0) + np.mean(g2_now, axis=0))
        xyz = rotate_about_point(xyz, R_com, com_origin)
        com_vec_after = np.mean(xyz[mask_chau], axis=0) - np.mean(xyz[mask_rest], axis=0)
        com_angle_after = float(
            np.degrees(np.arccos(np.clip(abs(np.dot(_unit(com_vec_after), pull_axis_vec)), 0.0, 1.0)))
        )
        print(f"[INFO] COM rotation applied: {com_angle_before:.1f}° → {com_angle_after:.1f}°")
    else:
        print(f"[INFO] COM already aligned within {args.skip_threshold_deg}° -- skipping")

    # --- Orient smaller group to +pull-axis ---
    if args.orient_smaller_positive:
        g1_now = xyz[mask_chau]
        g2_now = xyz[mask_rest]
        com1_now = np.mean(g1_now, axis=0)
        com2_now = np.mean(g2_now, axis=0)
        n_g1 = int(np.count_nonzero(mask_chau))
        n_g2 = int(np.count_nonzero(mask_rest))

        # Determine which is the smaller group
        if n_g1 <= n_g2:
            smaller_com = com1_now[pull_axis_idx]
            larger_com = com2_now[pull_axis_idx]
            smaller_label = "chauA"
            larger_label = "rest"
        else:
            smaller_com = com2_now[pull_axis_idx]
            larger_com = com1_now[pull_axis_idx]
            smaller_label = "rest"
            larger_label = "chauA"

        if smaller_com < larger_com:
            # Flip: reflect all coords along pull axis about the midpoint
            midpoint = np.mean(xyz[:, pull_axis_idx])
            xyz[:, pull_axis_idx] = 2.0 * midpoint - xyz[:, pull_axis_idx]
            print(f"[INFO] Flipped along {args.pull_dim}: "
                  f"{smaller_label} (smaller, {min(n_g1, n_g2)} atoms) now at +{args.pull_dim}, "
                  f"{larger_label} (larger, {max(n_g1, n_g2)} atoms) at -{args.pull_dim}")
        else:
            print(f"[INFO] Orientation OK: {smaller_label} (smaller) already at +{args.pull_dim}")

    # --- Build PBC-safe asymmetric box ---
    g1_now = xyz[mask_chau]
    g2_now = xyz[mask_rest]
    new_box = compute_box(xyz, pull_axis_idx, pull_distance_nm,
                          args.margin_pull_nm, args.margin_perp_nm,
                          g1_now, g2_now)
    print(
        f"[INFO] New box: {new_box[0]:.3f} x {new_box[1]:.3f} x {new_box[2]:.3f} nm "
        f"(pull axis {args.pull_dim}={new_box[pull_axis_idx]:.3f})"
    )

    # --- Position: center on X/Z, offset along Y so restrained group is near -wall ---
    n_g1_pos = int(np.count_nonzero(mask_chau))
    n_g2_pos = int(np.count_nonzero(mask_rest))
    larger_mask = mask_rest if n_g2_pos >= n_g1_pos else mask_chau
    xyz = position_for_pull(xyz, new_box, pull_axis_idx, larger_mask,
                            margin_nm=args.margin_pull_nm)

    # --- Write output ---
    atoms_out = {**atoms, "xyz": xyz}
    out_title = f"Aligned for {args.pull_dim}-pull ({method}, {angle_before_deg:.1f} -> {angle_after_deg:.1f} deg)"
    write_gro(args.out, atoms_out, new_box, out_title)
    print(f"[INFO] Wrote {args.out}")

    # --- Generate visual report ---
    if args.report_png is not None:
        # Center before-state in original box for fair comparison
        xyz_before_centered = center_protein(xyz_before, box_orig)
        generate_report_png(
            args.report_png,
            coords_before=xyz_before_centered,
            coords_after=xyz,
            g1_mask=mask_chau,
            g2_mask=mask_rest,
            box_before=box_orig,
            box_after=new_box,
            method=method,
            angle_before=angle_before_deg,
            angle_after=angle_after_deg,
        )

    print("[DONE] Auto-alignment complete")


if __name__ == "__main__":
    main()
