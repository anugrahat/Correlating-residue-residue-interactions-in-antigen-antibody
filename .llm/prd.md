# PRD: .gro File Quick Viewer

## Overview

Single-file CLI tool at `scripts/view_gro.py` that converts `.gro` files into interactive 3D HTML viewers.

---

## Component 1: CLI Entry Point

**File:** `scripts/view_gro.py` — top-level `main()` with argparse.

```python
#!/usr/bin/env python3
"""Quick 3D viewer for GROMACS .gro files.

Usage:  python view_gro.py structure.gro
        python view_gro.py structure.gro --color-by chain --dark
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

try:
    import plotly.graph_objects as go
except ImportError:
    go = None


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("gro", type=Path, help="Input .gro file")
    p.add_argument("--color-by", choices=["element", "chain"], default="element",
                   help="Atom coloring mode (default: element)")
    p.add_argument("--no-open", action="store_true",
                   help="Don't auto-open HTML in browser")
    p.add_argument("--no-box", action="store_true",
                   help="Hide simulation box wireframe")
    p.add_argument("--max-atoms", type=int, default=20000,
                   help="Downsample limit for display (default: 20000)")
    p.add_argument("--output", type=Path, default=None,
                   help="Output HTML path (default: <stem>_view.html next to input)")
    p.add_argument("--dark", action="store_true",
                   help="Dark background theme")
    return p.parse_args()
```

**Acceptance criteria:**
- `python scripts/view_gro.py file.gro` works with no other args
- All flags are optional with sane defaults
- Exits with clear error if plotly not installed or file not found

---

## Component 2: GRO Parser

Reuse the fixed-width parsing pattern from `pull_cli_view.py:155-189`.

```python
SOLVENT_NAMES = {"SOL", "WAT", "HOH", "TIP3", "TIP3P", "SPC", "SPCE"}
ION_NAMES = {"NA", "CL", "K", "MG", "CA", "ZN", "NA+", "CL-"}


def parse_gro(path: Path) -> Tuple[Dict[str, np.ndarray], np.ndarray, str]:
    """Parse .gro file. Returns (atoms_dict, box_vec, title)."""
    lines = path.read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"Invalid .gro file: {path}")
    title = lines[0].strip()
    n_atoms = int(lines[1].strip())
    atom_lines = lines[2 : 2 + n_atoms]
    box_tokens = lines[2 + n_atoms].split()
    box = np.array([float(box_tokens[0]), float(box_tokens[1]), float(box_tokens[2])])

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
    return atoms, box, title
```

**Acceptance criteria:**
- Handles standard GRO format (3-component and 9-component box lines)
- Returns numpy arrays matching `pull_cli_view.py` conventions

---

## Component 3: Solvent/Ion Stripping

```python
def strip_solvent(atoms: Dict[str, np.ndarray]) -> Tuple[Dict[str, np.ndarray], int]:
    """Remove water and ion residues. Returns (filtered_atoms, n_removed)."""
    exclude = SOLVENT_NAMES | ION_NAMES
    keep = ~np.isin(atoms["resname"], list(exclude))
    n_removed = int(np.count_nonzero(~keep))
    filtered = {k: v[keep] for k, v in atoms.items()}
    return filtered, n_removed
```

**Acceptance criteria:**
- Removes all standard water models and common ions
- Reports count of removed atoms

---

## Component 4: Atom Coloring

### Element-based (CPK)

```python
CPK_COLORS = {
    "C": "#909090", "N": "#3050F8", "O": "#FF0D0D", "S": "#FFFF30",
    "H": "#FFFFFF", "P": "#FF8000", "FE": "#E06633", "ZN": "#7D80B0",
    "CA": "#3DFF00", "MG": "#8AFF00",
}
DEFAULT_COLOR = "#FF1493"


def element_from_atomname(name: str) -> str:
    """Derive element symbol from GROMACS atom name."""
    stripped = name.lstrip("0123456789")
    if len(stripped) >= 2 and stripped[:2].upper() in CPK_COLORS:
        return stripped[:2].upper()
    if stripped and stripped[0].upper() in CPK_COLORS:
        return stripped[0].upper()
    return stripped[0].upper() if stripped else "X"


def color_by_element(atomnames: np.ndarray) -> List[str]:
    return [CPK_COLORS.get(element_from_atomname(str(n)), DEFAULT_COLOR) for n in atomnames]
```

### Chain-based

```python
CHAIN_PALETTE = [
    "#0ea5e9", "#ef4444", "#10b981", "#f59e0b",
    "#8b5cf6", "#ec4899", "#06b6d4", "#84cc16",
]


def detect_chains(resid: np.ndarray) -> np.ndarray:
    """Assign chain index by detecting residue ID discontinuities."""
    chain_idx = np.zeros(len(resid), dtype=int)
    current_chain = 0
    for i in range(1, len(resid)):
        if resid[i] < resid[i - 1]:  # residue numbering reset = new chain
            current_chain += 1
        chain_idx[i] = current_chain
    return chain_idx


def color_by_chain(resid: np.ndarray) -> List[str]:
    chains = detect_chains(resid)
    return [CHAIN_PALETTE[c % len(CHAIN_PALETTE)] for c in chains]
```

**Acceptance criteria:**
- Element detection handles GROMACS naming (CA, CB, OG1, NZ, etc.)
- Chain detection catches residue numbering resets
- Fallback color for unknown elements

---

## Component 5: Downsampling

```python
def downsample(atoms: Dict[str, np.ndarray], colors: List[str],
               max_atoms: int, seed: int = 42) -> Tuple[Dict[str, np.ndarray], List[str]]:
    n = len(atoms["xyz"])
    if n <= max_atoms:
        return atoms, colors
    rng = np.random.default_rng(seed)
    keep = np.sort(rng.choice(n, size=max_atoms, replace=False))
    filtered = {k: v[keep] for k, v in atoms.items()}
    colors_filtered = [colors[i] for i in keep]
    return filtered, colors_filtered
```

**Acceptance criteria:**
- Deterministic (fixed seed)
- Preserves atom-color correspondence

---

## Component 6: HTML Generation

```python
def box_wireframe(box: np.ndarray) -> Tuple[List, List, List]:
    """12 edges of rectangular box as Scatter3d line coords."""
    c = np.array([
        [0,0,0],[box[0],0,0],[box[0],box[1],0],[0,box[1],0],
        [0,0,box[2]],[box[0],0,box[2]],[box[0],box[1],box[2]],[0,box[1],box[2]],
    ], dtype=float)
    edges = [(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
    xs, ys, zs = [], [], []
    for i, j in edges:
        xs.extend([c[i,0], c[j,0], None])
        ys.extend([c[i,1], c[j,1], None])
        zs.extend([c[i,2], c[j,2], None])
    return xs, ys, zs


def build_hover_text(atoms: Dict[str, np.ndarray]) -> List[str]:
    return [
        f"{rn} {ri} : {an}"
        for rn, ri, an in zip(atoms["resname"], atoms["resid"], atoms["atomname"])
    ]


def generate_html(
    atoms: Dict[str, np.ndarray],
    box: np.ndarray,
    colors: List[str],
    title: str,
    output: Path,
    show_box: bool = True,
    dark: bool = False,
) -> None:
    if go is None:
        print("ERROR: plotly is required. Install with: pip install plotly", file=sys.stderr)
        sys.exit(1)

    fig = go.Figure()

    # Box wireframe
    if show_box:
        bx, by, bz = box_wireframe(box)
        fig.add_trace(go.Scatter3d(
            x=bx, y=by, z=bz, mode="lines",
            line=dict(color="white" if dark else "black", width=3),
            name="box", hoverinfo="skip",
        ))

    # Atoms
    xyz = atoms["xyz"]
    hover = build_hover_text(atoms)
    fig.add_trace(go.Scatter3d(
        x=xyz[:, 0], y=xyz[:, 1], z=xyz[:, 2],
        mode="markers",
        marker=dict(size=2.5, color=colors, opacity=0.9),
        text=hover, hoverinfo="text",
        name=f"atoms ({len(xyz)})",
    ))

    bg = "#1a1a2e" if dark else "#ffffff"
    fg = "#ffffff" if dark else "#000000"
    grid = "#333355" if dark else "#cccccc"

    fig.update_layout(
        title=dict(text=title, font=dict(color=fg)),
        paper_bgcolor=bg,
        scene=dict(
            xaxis=dict(title="X (nm)", range=[0, float(box[0])],
                       backgroundcolor=bg, gridcolor=grid, color=fg),
            yaxis=dict(title="Y (nm)", range=[0, float(box[1])],
                       backgroundcolor=bg, gridcolor=grid, color=fg),
            zaxis=dict(title="Z (nm)", range=[0, float(box[2])],
                       backgroundcolor=bg, gridcolor=grid, color=fg),
            aspectmode="data",
            bgcolor=bg,
        ),
        legend=dict(x=0.01, y=0.99, font=dict(color=fg)),
        margin=dict(l=0, r=0, b=0, t=40),
    )

    fig.write_html(str(output), include_plotlyjs=True, full_html=True)
```

**Acceptance criteria:**
- Self-contained HTML (no CDN dependency)
- Hover shows residue info
- Dark mode works with readable labels
- Box wireframe adapts color to background

---

## Component 7: Auto-Open & Main Flow

```python
def auto_open(path: Path) -> None:
    opener = shutil.which("xdg-open")
    if opener:
        subprocess.Popen([opener, str(path.resolve())],
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print(f"Opened in browser via xdg-open")
    else:
        print(f"Auto-open skipped (xdg-open not found). Open manually:")
        print(f"  file://{path.resolve()}")


def main() -> None:
    args = parse_args()

    if not args.gro.exists():
        print(f"ERROR: file not found: {args.gro}", file=sys.stderr)
        sys.exit(1)

    atoms, box, gro_title = parse_gro(args.gro)
    total = len(atoms["xyz"])

    atoms, n_stripped = strip_solvent(atoms)
    protein_count = len(atoms["xyz"])

    if args.color_by == "chain":
        colors = color_by_chain(atoms["resid"])
    else:
        colors = color_by_element(atoms["atomname"])

    atoms, colors = downsample(atoms, colors, args.max_atoms)
    display_count = len(atoms["xyz"])

    output = args.output or args.gro.with_name(args.gro.stem + "_view.html")

    title = f"{args.gro.name} — {protein_count} atoms (solvent stripped)"

    print(f"Input       : {args.gro}")
    print(f"Total atoms : {total}")
    print(f"Stripped    : {n_stripped} (solvent + ions)")
    print(f"Displaying  : {display_count}" +
          (f" (downsampled from {protein_count})" if display_count < protein_count else ""))
    print(f"Color mode  : {args.color_by}")

    generate_html(atoms, box, colors, title, output,
                  show_box=not args.no_box, dark=args.dark)

    print(f"Output      : {output}")

    if not args.no_open:
        auto_open(output)


if __name__ == "__main__":
    main()
```

**Acceptance criteria:**
- Runs end-to-end with `python scripts/view_gro.py file.gro`
- Terminal output summarizes what happened
- HTML opens in browser automatically
- All error paths print clear messages

---

## Implementation Order

1. Write `scripts/view_gro.py` with all components above assembled into a single file
2. Test with a real `.gro` file from the simulation output directories
3. Verify HTML opens and renders correctly

## Non-Goals

- No bond rendering
- No trajectory support
- No ribbon/surface modes
- No editing or measurement tools
