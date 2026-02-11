#!/usr/bin/env python3
"""Terminal-first pull setup inspector.

Given a .gro and .ndx, this script:
1) strips solvent for inspection,
2) computes pull-group and box metrics in terminal output,
3) writes quick projection PNGs for optional visual checking.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np

try:
    import plotly.graph_objects as go
except Exception:
    go = None


SOLVENT_NAMES = {"SOL", "WAT", "HOH", "TIP3", "TIP3P", "SPC", "SPCE"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gro", type=Path, required=True, help="Input .gro structure")
    p.add_argument("--ndx", type=Path, required=True, help="Index file with pull groups")
    p.add_argument("--group1", default="chauA", help="Pull group 1 name in .ndx")
    p.add_argument("--group2", default="rest", help="Pull group 2 name in .ndx")
    p.add_argument(
        "--pull-dim",
        choices=["X", "Y", "Z"],
        default="Y",
        help="Pull axis dimension from mdp (default Y)",
    )
    p.add_argument(
        "--edge-warn-nm",
        type=float,
        default=1.0,
        help="Warn if group gets closer than this to a box wall",
    )
    p.add_argument(
        "--edge-loose-nm",
        type=float,
        default=2.5,
        help="Suggest box could be reduced if all groups are farther than this",
    )
    p.add_argument(
        "--target-margin-nm",
        type=float,
        default=1.5,
        help="Target margin to build suggested new box around protein",
    )
    p.add_argument(
        "--outdir",
        type=Path,
        default=Path("view_cli"),
        help="Output folder for generated files",
    )
    p.add_argument(
        "--interactive-html",
        dest="interactive_html",
        action="store_true",
        default=True,
        help="Write interactive 3D HTML viewer (default on, requires plotly).",
    )
    p.add_argument(
        "--no-interactive-html",
        dest="interactive_html",
        action="store_false",
        help="Disable interactive 3D HTML output.",
    )
    p.add_argument(
        "--open-html",
        action="store_true",
        help="Try to open the interactive HTML in your default browser.",
    )
    p.add_argument(
        "--auto-align-to-pull",
        action="store_true",
        help="Rotate so COM(group2-group1) aligns with pull axis before reporting/writing outputs.",
    )
    p.add_argument(
        "--align-interface-normal",
        action="store_true",
        help="With --auto-align-to-pull, align interface-plane normal (from contact atoms) to pull axis.",
    )
    p.add_argument(
        "--interface-cutoff-nm",
        type=float,
        default=0.6,
        help="Contact cutoff for interface atom detection (used with --align-interface-normal).",
    )
    p.add_argument(
        "--emit-gmx-rotate-cmd",
        action="store_true",
        help="Print and save a gmx editconf -rotate command for the computed alignment.",
    )
    p.add_argument(
        "--rebox-mode",
        choices=["none", "axis", "isotropic"],
        default="axis",
        help="Reboxing policy after alignment: none, axis-aware, or isotropic.",
    )
    p.add_argument(
        "--rebox-margin-nm",
        type=float,
        default=1.5,
        help="Margin used for non-pull axes in rebox suggestion/output.",
    )
    p.add_argument(
        "--rebox-pull-margin-nm",
        type=float,
        default=2.5,
        help="Margin used for pull axis in axis-aware rebox mode.",
    )
    p.add_argument(
        "--allow-shrink",
        action="store_true",
        help="Allow rebox dimensions to shrink if suggested box is smaller than input box.",
    )
    p.add_argument(
        "--write-aligned-gro",
        type=Path,
        default=None,
        help="Output full-system .gro after transform/rebox (default: <outdir>/aligned_for_pull.gro).",
    )
    p.add_argument(
        "--write-groups-ndx",
        type=Path,
        default=None,
        help="Output .ndx containing group1/group2 atom lists (default: <outdir>/pull_groups.ndx).",
    )
    p.add_argument(
        "--report-txt",
        type=Path,
        default=None,
        help="Write terminal summary to text file (default: <outdir>/summary.txt).",
    )
    p.add_argument(
        "--roll-around-pull-deg",
        type=float,
        default=0.0,
        help="Extra rigid rotation (deg) around pull axis after alignment; keeps perpendicular condition.",
    )
    return p.parse_args()


def parse_gro(path: Path) -> Tuple[Dict[str, np.ndarray], np.ndarray]:
    lines = path.read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"Invalid .gro file: {path}")
    n_atoms = int(lines[1].strip())
    atom_lines = lines[2 : 2 + n_atoms]
    box_line = lines[2 + n_atoms].split()

    if len(box_line) not in (3, 9):
        raise ValueError(f"Unsupported box format in {path}: {box_line}")
    box = np.array([float(box_line[0]), float(box_line[1]), float(box_line[2])], dtype=float)

    resid = np.zeros(n_atoms, dtype=int)
    resname = []
    atomname = []
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
    return atoms, box


def parse_ndx(path: Path) -> Dict[str, np.ndarray]:
    groups: Dict[str, List[int]] = {}
    current = None
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            current = line[1:-1].strip()
            groups[current] = []
            continue
        if current is None:
            continue
        groups[current].extend(int(x) for x in line.split())
    return {k: np.array(v, dtype=int) for k, v in groups.items()}


def write_ndx(path: Path, groups: Dict[str, np.ndarray]) -> None:
    with path.open("w") as f:
        for name, atom_ids in groups.items():
            f.write(f"[ {name} ]\n")
            vals = [str(int(x)) for x in atom_ids]
            for i in range(0, len(vals), 15):
                f.write(" ".join(vals[i : i + 15]) + "\n")
            f.write("\n")


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


def make_whole(coords: np.ndarray, box: np.ndarray) -> np.ndarray:
    if len(coords) == 0:
        return coords.copy()
    out = np.empty_like(coords)
    ref = coords[0].copy()
    out[0] = ref
    for i in range(1, len(coords)):
        d = coords[i] - ref
        d -= box * np.round(d / box)
        out[i] = ref + d
    return out


def center_in_box(coords: np.ndarray, box: np.ndarray, center_of: np.ndarray) -> np.ndarray:
    shift = (box / 2.0) - np.mean(center_of, axis=0)
    return np.mod(coords + shift, box)


def projection_plot(
    out_png: Path,
    box: np.ndarray,
    coords_all: np.ndarray,
    coords_g1: np.ndarray,
    coords_g2: np.ndarray,
    title: str,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    pairs = [(0, 1, "X", "Y"), (1, 2, "Y", "Z"), (0, 2, "X", "Z")]
    for ax, (i, j, lx, ly) in zip(axes, pairs):
        ax.scatter(coords_all[:, i], coords_all[:, j], s=1, alpha=0.15, c="#6b7280")
        ax.scatter(coords_g1[:, i], coords_g1[:, j], s=6, alpha=0.8, c="#0ea5e9", label="group1")
        ax.scatter(coords_g2[:, i], coords_g2[:, j], s=6, alpha=0.8, c="#ef4444", label="group2")
        ax.set_xlim(0, box[i])
        ax.set_ylim(0, box[j])
        ax.set_xlabel(f"{lx} (nm)")
        ax.set_ylabel(f"{ly} (nm)")
        ax.set_aspect("equal", adjustable="box")
    axes[0].legend(loc="upper right", fontsize=8)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


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
        # 180-degree rotation: R = -I + 2*u*u^T
        return -np.eye(3, dtype=float) + 2.0 * np.outer(u, u)

    v = np.cross(a, b)
    s = np.linalg.norm(v)
    k = _skew(v)
    return np.eye(3, dtype=float) + k + (k @ k) * ((1.0 - c) / (s * s))


def rotation_matrix_axis_angle(axis: np.ndarray, angle_rad: float) -> np.ndarray:
    u = _unit(axis)
    k = _skew(u)
    c = float(np.cos(angle_rad))
    s = float(np.sin(angle_rad))
    return np.eye(3, dtype=float) + s * k + (1.0 - c) * (k @ k)


def rotate_about_point(coords: np.ndarray, R: np.ndarray, origin: np.ndarray) -> np.ndarray:
    return (coords - origin) @ R.T + origin


def rot_x(a: float) -> np.ndarray:
    c, s = float(np.cos(a)), float(np.sin(a))
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]], dtype=float)


def rot_y(a: float) -> np.ndarray:
    c, s = float(np.cos(a)), float(np.sin(a))
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]], dtype=float)


def rot_z(a: float) -> np.ndarray:
    c, s = float(np.cos(a)), float(np.sin(a))
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=float)


def euler_xyz_deg_from_matrix(R: np.ndarray) -> Tuple[float, float, float]:
    """Extract x,y,z angles (deg) for composition R = Rz(z) * Ry(y) * Rx(x)."""
    sy = float(np.sqrt(R[0, 0] * R[0, 0] + R[1, 0] * R[1, 0]))
    singular = sy < 1e-10
    if not singular:
        x = float(np.arctan2(R[2, 1], R[2, 2]))
        y = float(np.arctan2(-R[2, 0], sy))
        z = float(np.arctan2(R[1, 0], R[0, 0]))
    else:
        x = float(np.arctan2(-R[1, 2], R[1, 1]))
        y = float(np.arctan2(-R[2, 0], sy))
        z = 0.0
    return float(np.degrees(x)), float(np.degrees(y)), float(np.degrees(z))


def matrix_from_euler_xyz_deg(x_deg: float, y_deg: float, z_deg: float) -> np.ndarray:
    xr = float(np.radians(x_deg))
    yr = float(np.radians(y_deg))
    zr = float(np.radians(z_deg))
    return rot_z(zr) @ rot_y(yr) @ rot_x(xr)


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
    mask1, mask2 = interface_contact_masks(g1, g2, cutoff_nm=cutoff_nm)
    n1 = int(np.count_nonzero(mask1))
    n2 = int(np.count_nonzero(mask2))
    if n1 < 3 or n2 < 3:
        raise ValueError(f"Not enough interface contact atoms (group1={n1}, group2={n2})")

    g1_if = g1[mask1]
    g2_if = g2[mask2]
    pts = np.vstack([g1_if, g2_if])
    center = np.mean(pts, axis=0)
    centered = pts - center
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    normal = vh[-1]
    return normal, center, n1, n2


def _box_wireframe_xyz(box: np.ndarray) -> Tuple[List[float], List[float], List[float]]:
    corners = np.array(
        [
            [0.0, 0.0, 0.0],
            [box[0], 0.0, 0.0],
            [box[0], box[1], 0.0],
            [0.0, box[1], 0.0],
            [0.0, 0.0, box[2]],
            [box[0], 0.0, box[2]],
            [box[0], box[1], box[2]],
            [0.0, box[1], box[2]],
        ],
        dtype=float,
    )
    edges = [
        (0, 1), (1, 2), (2, 3), (3, 0),
        (4, 5), (5, 6), (6, 7), (7, 4),
        (0, 4), (1, 5), (2, 6), (3, 7),
    ]
    xs: List[float] = []
    ys: List[float] = []
    zs: List[float] = []
    for i, j in edges:
        xs.extend([corners[i, 0], corners[j, 0], None])
        ys.extend([corners[i, 1], corners[j, 1], None])
        zs.extend([corners[i, 2], corners[j, 2], None])
    return xs, ys, zs


def _downsample(coords: np.ndarray, max_points: int, seed: int = 42) -> np.ndarray:
    if len(coords) <= max_points:
        return coords
    rng = np.random.default_rng(seed)
    keep = rng.choice(len(coords), size=max_points, replace=False)
    keep.sort()
    return coords[keep]


def write_interactive_html(
    out_html: Path,
    box: np.ndarray,
    coords_all: np.ndarray,
    coords_g1: np.ndarray,
    coords_g2: np.ndarray,
    title: str,
    pull_axis_vec: np.ndarray | None = None,
    interface_normal: np.ndarray | None = None,
    interface_center: np.ndarray | None = None,
) -> bool:
    if go is None:
        return False

    bg = _downsample(coords_all, max_points=15000)
    bx, by, bz = _box_wireframe_xyz(box)

    fig = go.Figure()
    fig.add_trace(
        go.Scatter3d(
            x=bx,
            y=by,
            z=bz,
            mode="lines",
            line=dict(color="black", width=4),
            name="box",
        )
    )
    fig.add_trace(
        go.Scatter3d(
            x=bg[:, 0],
            y=bg[:, 1],
            z=bg[:, 2],
            mode="markers",
            marker=dict(size=1.8, color="gray", opacity=0.10),
            name="protein (no solvent)",
        )
    )
    fig.add_trace(
        go.Scatter3d(
            x=coords_g1[:, 0],
            y=coords_g1[:, 1],
            z=coords_g1[:, 2],
            mode="markers",
            marker=dict(size=2.6, color="#0ea5e9", opacity=0.85),
            name="group1",
        )
    )
    fig.add_trace(
        go.Scatter3d(
            x=coords_g2[:, 0],
            y=coords_g2[:, 1],
            z=coords_g2[:, 2],
            mode="markers",
            marker=dict(size=2.6, color="#ef4444", opacity=0.85),
            name="group2",
        )
    )
    # Optional orientation guides for unambiguous visual validation.
    center = box / 2.0
    guide_len = float(np.min(box) * 0.45)
    if pull_axis_vec is not None:
        a = _unit(pull_axis_vec)
        p0 = center - a * guide_len
        p1 = center + a * guide_len
        fig.add_trace(
            go.Scatter3d(
                x=[p0[0], p1[0]],
                y=[p0[1], p1[1]],
                z=[p0[2], p1[2]],
                mode="lines",
                line=dict(color="#10b981", width=8),
                name="pull-axis",
            )
        )
    if interface_normal is not None:
        n = _unit(interface_normal)
        c = interface_center if interface_center is not None else center
        q0 = c - n * (guide_len * 0.55)
        q1 = c + n * (guide_len * 0.55)
        fig.add_trace(
            go.Scatter3d(
                x=[q0[0], q1[0]],
                y=[q0[1], q1[1]],
                z=[q0[2], q1[2]],
                mode="lines",
                line=dict(color="#f59e0b", width=8),
                name="interface-normal",
            )
        )
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis=dict(title="X (nm)", range=[0.0, float(box[0])]),
            yaxis=dict(title="Y (nm)", range=[0.0, float(box[1])]),
            zaxis=dict(title="Z (nm)", range=[0.0, float(box[2])]),
            aspectmode="data",
        ),
        legend=dict(x=0.01, y=0.99),
        margin=dict(l=0, r=0, b=0, t=40),
    )
    # Self-contained HTML for offline use on cluster/login nodes.
    fig.write_html(str(out_html), include_plotlyjs=True, full_html=True)
    return True


def write_solvent_free_pdb(path: Path, atoms: Dict[str, np.ndarray]) -> None:
    xyz = atoms["xyz"]
    resid = atoms["resid"]
    resname = atoms["resname"]
    atomname = atoms["atomname"]
    atomnr = atoms["atomnr"]
    with path.open("w") as f:
        for i in range(len(xyz)):
            # Chain is omitted in GRO, keep blank.
            f.write(
                "ATOM  {serial:5d} {name:<4s} {resn:>3s} {chain:1s}{resid:4d}"
                "    {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bf:6.2f}          {el:>2s}\n".format(
                    serial=int(atomnr[i]),
                    name=str(atomname[i])[:4],
                    resn=str(resname[i])[:3],
                    chain=" ",
                    resid=int(resid[i]),
                    x=float(xyz[i, 0] * 10.0),  # nm -> A
                    y=float(xyz[i, 1] * 10.0),
                    z=float(xyz[i, 2] * 10.0),
                    occ=1.00,
                    bf=0.00,
                    el=str(atomname[i])[0],
                )
            )
        f.write("END\n")


def min_wall_gap(coords: np.ndarray, box: np.ndarray) -> np.ndarray:
    lower = np.min(coords, axis=0)
    upper = box - np.max(coords, axis=0)
    return np.minimum(lower, upper)


def fmt3(v: np.ndarray) -> str:
    return f"({v[0]:.3f}, {v[1]:.3f}, {v[2]:.3f})"


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    atoms, box = parse_gro(args.gro)
    groups = parse_ndx(args.ndx)
    if args.group1 not in groups:
        raise ValueError(f"Group '{args.group1}' not found in {args.ndx}")
    if args.group2 not in groups:
        raise ValueError(f"Group '{args.group2}' not found in {args.ndx}")

    report_lines: List[str] = []

    def emit(msg: str) -> None:
        print(msg)
        report_lines.append(msg)

    idx1 = groups[args.group1] - 1
    idx2 = groups[args.group2] - 1
    pull_axis = {"X": 0, "Y": 1, "Z": 2}[args.pull_dim]
    axis_vec = {
        "X": np.array([1.0, 0.0, 0.0], dtype=float),
        "Y": np.array([0.0, 1.0, 0.0], dtype=float),
        "Z": np.array([0.0, 0.0, 1.0], dtype=float),
    }[args.pull_dim]

    xyz_all = atoms["xyz"].copy()
    keep_mask = ~np.isin(atoms["resname"], list(SOLVENT_NAMES))
    atoms_keep = {
        "resid": atoms["resid"][keep_mask],
        "resname": atoms["resname"][keep_mask],
        "atomname": atoms["atomname"][keep_mask],
        "atomnr": atoms["atomnr"][keep_mask],
        "xyz": xyz_all[keep_mask].copy(),
    }

    # Build robust whole-group references and center whole system accordingly.
    g1_whole = make_whole(xyz_all[idx1], box)
    g2_whole = make_whole(xyz_all[idx2], box)
    all_whole = np.vstack([g1_whole, g2_whole])
    shift0 = (box / 2.0) - np.mean(all_whole, axis=0)
    xyz_all_view = np.mod(xyz_all + shift0, box)
    g1_view = np.mod(g1_whole + shift0, box)
    g2_view = np.mod(g2_whole + shift0, box)
    xyz_view = xyz_all_view[keep_mask]

    aligned = False
    rolled = False
    rotation_used = np.eye(3, dtype=float)
    rotate_vec_before = None
    rotate_target = None
    angle_before_deg = None
    angle_after_deg = None
    plane_before_deg = None
    plane_after_deg = None
    align_target = "com-vector"
    align_fallback = None
    interface_n1 = None
    interface_n2 = None

    if args.auto_align_to_pull:
        com1_pre = np.mean(g1_view, axis=0)
        com2_pre = np.mean(g2_view, axis=0)
        dvec_pre = com2_pre - com1_pre
        if np.linalg.norm(dvec_pre) < 1e-8:
            raise ValueError("COM separation is ~0; cannot align to pull axis.")

        rot_vec = dvec_pre
        rot_origin = 0.5 * (com1_pre + com2_pre)
        if args.align_interface_normal:
            try:
                normal_pre, iface_center, n1c, n2c = interface_normal_from_contacts(
                    g1_view, g2_view, cutoff_nm=args.interface_cutoff_nm
                )
                rot_vec = normal_pre
                rot_origin = iface_center
                align_target = "interface-normal"
                interface_n1 = n1c
                interface_n2 = n2c
            except ValueError as exc:
                align_fallback = str(exc)
                align_target = "com-vector (fallback)"

        target = axis_vec if np.dot(rot_vec, axis_vec) >= 0.0 else -axis_vec
        R = rotation_matrix_from_vectors(rot_vec, target)
        rotation_used = R.copy()
        rotate_vec_before = _unit(rot_vec)
        rotate_target = _unit(target)

        xyz_all_rot = rotate_about_point(xyz_all_view, R, rot_origin)
        g1_rot = rotate_about_point(g1_view, R, rot_origin)
        g2_rot = rotate_about_point(g2_view, R, rot_origin)

        # Recenter after rotation for stable box-wall metrics.
        all_rot = np.vstack([g1_rot, g2_rot])
        shift1 = (box / 2.0) - np.mean(all_rot, axis=0)
        xyz_all_view = np.mod(xyz_all_rot + shift1, box)
        g1_view = np.mod(g1_rot + shift1, box)
        g2_view = np.mod(g2_rot + shift1, box)
        xyz_view = xyz_all_view[keep_mask]
        aligned = True

        vec_after = R @ _unit(rot_vec)
        angle_before_deg = float(
            np.degrees(np.arccos(np.clip(abs(np.dot(_unit(rot_vec), axis_vec)), -1.0, 1.0)))
        )
        angle_after_deg = float(
            np.degrees(np.arccos(np.clip(abs(np.dot(_unit(vec_after), axis_vec)), -1.0, 1.0)))
        )
        plane_before_deg = 90.0 - angle_before_deg
        plane_after_deg = 90.0 - angle_after_deg

    # Optional extra roll around pull axis (for clearly visible orientation changes).
    if abs(args.roll_around_pull_deg) > 1e-10:
        R_roll = rotation_matrix_axis_angle(axis_vec, float(np.radians(args.roll_around_pull_deg)))
        roll_origin = 0.5 * (np.mean(g1_view, axis=0) + np.mean(g2_view, axis=0))
        xyz_all_roll = rotate_about_point(xyz_all_view, R_roll, roll_origin)
        g1_roll = rotate_about_point(g1_view, R_roll, roll_origin)
        g2_roll = rotate_about_point(g2_view, R_roll, roll_origin)
        xyz_all_view = np.mod(xyz_all_roll, box)
        g1_view = np.mod(g1_roll, box)
        g2_view = np.mod(g2_roll, box)
        xyz_view = xyz_all_view[keep_mask]
        rolled = True

    # Rebox suggestion / application.
    if args.rebox_mode == "axis":
        margin_vec = np.array([args.rebox_margin_nm] * 3, dtype=float)
        margin_vec[pull_axis] = args.rebox_pull_margin_nm
    elif args.rebox_mode == "isotropic":
        m = max(args.rebox_margin_nm, args.rebox_pull_margin_nm)
        margin_vec = np.array([m, m, m], dtype=float)
    else:
        margin_vec = np.array([args.target_margin_nm] * 3, dtype=float)

    protein_span = np.ptp(xyz_view, axis=0)
    suggested_box = protein_span + 2.0 * margin_vec
    if args.rebox_mode == "none":
        final_box = box.copy()
    else:
        final_box = suggested_box.copy() if args.allow_shrink else np.maximum(box, suggested_box)
        center_ref = np.vstack([g1_view, g2_view])
        shiftb = (final_box / 2.0) - np.mean(center_ref, axis=0)
        xyz_all_view = np.mod(xyz_all_view + shiftb, final_box)
        g1_view = np.mod(g1_view + shiftb, final_box)
        g2_view = np.mod(g2_view + shiftb, final_box)
        xyz_view = xyz_all_view[keep_mask]

    com1 = np.mean(g1_view, axis=0)
    com2 = np.mean(g2_view, axis=0)
    dvec = com2 - com1
    gap1 = min_wall_gap(g1_view, final_box)
    gap2 = min_wall_gap(g2_view, final_box)

    aligned_gro = args.write_aligned_gro or (args.outdir / "aligned_for_pull.gro")
    out_ndx = args.write_groups_ndx or (args.outdir / "pull_groups.ndx")
    report_txt = args.report_txt or (args.outdir / "summary.txt")
    solvent_free_pdb = args.outdir / "solvent_free_view.pdb"
    projection_png = args.outdir / "pull_groups_box_projections.png"
    interactive_html = args.outdir / "pull_groups_box_interactive.html"

    write_gro(aligned_gro, {**atoms, "xyz": xyz_all_view}, final_box, title="Aligned for pull")
    write_ndx(out_ndx, {args.group1: groups[args.group1], args.group2: groups[args.group2]})
    atoms_keep["xyz"] = xyz_view
    write_solvent_free_pdb(solvent_free_pdb, atoms_keep)

    projection_plot(
        projection_png,
        box=final_box,
        coords_all=xyz_view,
        coords_g1=g1_view,
        coords_g2=g2_view,
        title=f"{args.gro.name}: {args.group1} vs {args.group2} (solvent-free)",
    )
    html_ok = False
    iface_n_final = None
    iface_c_final = None
    iface_angle_final = None
    try:
        iface_n_final, iface_c_final, _, _ = interface_normal_from_contacts(
            g1_view, g2_view, cutoff_nm=args.interface_cutoff_nm
        )
        iface_angle_final = float(
            np.degrees(np.arccos(np.clip(abs(np.dot(_unit(iface_n_final), axis_vec)), -1.0, 1.0)))
        )
    except ValueError:
        iface_n_final = None
        iface_c_final = None
        iface_angle_final = None

    if args.interactive_html:
        angle_tag = (
            f" | interface-normal vs pull: {iface_angle_final:.2f} deg"
            if iface_angle_final is not None
            else ""
        )
        html_ok = write_interactive_html(
            interactive_html,
            box=final_box,
            coords_all=xyz_view,
            coords_g1=g1_view,
            coords_g2=g2_view,
            title=f"{args.gro.name}: {args.group1} vs {args.group2} (interactive, solvent-free){angle_tag}",
            pull_axis_vec=axis_vec,
            interface_normal=iface_n_final,
            interface_center=iface_c_final,
        )

    emit(f"Input GRO           : {args.gro}")
    emit(f"Input NDX           : {args.ndx}")
    emit(f"Input box (nm)      : {fmt3(box)}")
    emit(f"Final box (nm)      : {fmt3(final_box)}")
    emit(f"Rebox mode          : {args.rebox_mode} (allow_shrink={'yes' if args.allow_shrink else 'no'})")
    emit(f"Group1              : {args.group1} ({len(idx1)} atoms)")
    emit(f"Group2              : {args.group2} ({len(idx2)} atoms)")
    emit(f"Auto align to pull  : {'ON' if aligned else 'OFF'}")
    emit(f"Extra roll (deg)    : {args.roll_around_pull_deg:.3f} ({'applied' if rolled else 'none'})")
    if aligned:
        emit(f"Alignment target    : {align_target}")
    if align_target.startswith("interface-normal") and interface_n1 is not None and interface_n2 is not None:
        emit(f"Interface contacts  : {args.group1}={interface_n1} atoms, {args.group2}={interface_n2} atoms")
    if align_fallback:
        emit(f"Align fallback      : {align_fallback}")
    if aligned and angle_before_deg is not None and angle_after_deg is not None:
        emit(f"Normal vs pull axis : {angle_before_deg:.2f} deg -> {angle_after_deg:.2f} deg")
    if aligned and plane_before_deg is not None and plane_after_deg is not None:
        emit(f"Plane vs pull axis  : {plane_before_deg:.2f} deg -> {plane_after_deg:.2f} deg (90 deg = perpendicular)")
    emit(f"COM group1 (nm)     : {fmt3(com1)}")
    emit(f"COM group2 (nm)     : {fmt3(com2)}")
    emit(f"COM delta (nm)      : {fmt3(dvec)}")
    emit(f"Pull-axis ({args.pull_dim}) delta (nm): {dvec[pull_axis]:.3f}")
    emit(f"Min wall gap {args.group1} (nm): {fmt3(gap1)}")
    emit(f"Min wall gap {args.group2} (nm): {fmt3(gap2)}")
    emit(f"Protein span (nm)   : {fmt3(protein_span)}")
    emit(f"Suggested box (nm)  : {fmt3(suggested_box)}  (margins {fmt3(margin_vec)})")

    tight = min(np.min(gap1), np.min(gap2))
    if tight < args.edge_warn_nm:
        emit(f"Recommendation      : INCREASE box (closest wall gap {tight:.3f} nm < {args.edge_warn_nm:.3f} nm)")
    elif tight > args.edge_loose_nm:
        emit(f"Recommendation      : Box may be reduced (all groups > {args.edge_loose_nm:.3f} nm from walls)")
    else:
        emit("Recommendation      : Box size looks reasonable for current pull setup")

    emit(f"Wrote               : {aligned_gro}")
    emit(f"Wrote               : {out_ndx}")
    emit(f"Wrote               : {solvent_free_pdb}")
    emit(f"Wrote               : {projection_png}")
    if args.interactive_html and html_ok:
        emit(f"Wrote               : {interactive_html}")
        if args.open_html:
            opener = shutil.which("xdg-open")
            if opener is None:
                emit("Open action         : skipped (xdg-open not found)")
            else:
                try:
                    subprocess.Popen(
                        [opener, str(interactive_html.resolve())],
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL,
                    )
                    emit("Open action         : requested via xdg-open")
                except Exception as exc:
                    emit(f"Open action         : failed ({exc})")
    elif args.interactive_html and not html_ok:
        emit("Interactive viewer  : plotly not available, skipped HTML output")

    if args.emit_gmx_rotate_cmd and aligned and rotate_vec_before is not None and rotate_target is not None:
        rx_deg, ry_deg, rz_deg = euler_xyz_deg_from_matrix(rotation_used)
        R_chk = matrix_from_euler_xyz_deg(rx_deg, ry_deg, rz_deg)
        v_chk = _unit(R_chk @ rotate_vec_before)
        err_deg = float(
            np.degrees(
                np.arccos(np.clip(np.dot(v_chk, rotate_target), -1.0, 1.0))
            )
        )
        gmx_cmd = (
            f"gmx editconf -f {args.gro} -o {args.outdir / 'aligned_gmx_editconf.gro'} "
            f"-rotate {rx_deg:.6f} {ry_deg:.6f} {rz_deg:.6f}"
        )
        rotate_file = args.outdir / "gmx_editconf_rotate_command.txt"
        rotate_file.write_text(gmx_cmd + "\n")
        emit(f"GMX rotate (deg)    : X={rx_deg:.6f} Y={ry_deg:.6f} Z={rz_deg:.6f}")
        emit(f"GMX rotate error    : {err_deg:.6f} deg (vector-target check)")
        emit(f"GMX command         : {gmx_cmd}")
        emit(f"Wrote               : {rotate_file}")

    report_txt.write_text("\n".join(report_lines) + "\n")


if __name__ == "__main__":
    main()
