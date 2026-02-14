# Design: .gro File Quick Viewer

## Executive Summary

A lightweight CLI tool that converts any `.gro` file into a self-contained, interactive 3D HTML viewer and auto-opens it in the browser. One command, instant visual sanity check — no VMD, no PyMol.

## Architecture

### Single-file script: `scripts/view_gro.py`

Follows project conventions (standalone, no cross-module imports, cluster-portable).

**Pipeline:** `.gro` file → parse → strip solvent → generate Plotly HTML → auto-open in browser

### Core Components

#### 1. GRO Parser (reuse pattern from `pull_cli_view.py`)
- Fixed-width column parsing (GRO format: resid, resname, atomname, atomnr, x, y, z)
- Extract box vectors from last line
- Return atoms dict + box array

#### 2. Solvent Stripper
- Remove residues matching: `SOL, WAT, HOH, TIP3, TIP3P, SPC, SPCE` (same set as `pull_cli_view.py`)
- Also strip common ions: `NA, CL, K, MG, CA, ZN` (residue names)
- Report count of stripped atoms in terminal

#### 3. Element Color Mapper (CPK scheme)
- Derive element from first character of atom name (standard GROMACS convention)
- Color map: C=`#909090`, N=`#3050F8`, O=`#FF0D0D`, S=`#FFFF30`, H=`#FFFFFF`, P=`#FF8000`, Fe=`#E06633`, default=`#FF1493`
- Per-atom color array passed to Plotly

#### 4. Chain/Residue Detector
- Detect chain breaks by residue ID gaps (>1 jump = new chain)
- Assign distinct colors per chain as an alternative coloring mode
- Default: color by element. Flag `--color-by chain` switches mode.

#### 5. HTML Generator (Plotly, self-contained)
- `Scatter3d` with per-atom coloring, marker size ~2.5
- Box wireframe overlay (12 edges, black lines)
- Downsample to 20,000 atoms max for browser performance (random subset, deterministic seed)
- `include_plotlyjs=True` for fully offline HTML
- Dark background option for better contrast
- Hover text: `resname resid : atomname`

#### 6. Auto-Open
- Use `xdg-open` (Linux) to launch HTML in default browser
- `--no-open` flag to suppress auto-open (just generate file)
- Print file path to terminal regardless

### CLI Interface

```
python scripts/view_gro.py <file.gro> [options]

Positional:
  file.gro              Path to .gro file

Options:
  --color-by {element,chain}   Coloring mode (default: element)
  --no-open                    Don't auto-open in browser
  --no-box                     Hide simulation box wireframe
  --max-atoms N                Downsample limit (default: 20000)
  --output PATH                Output HTML path (default: <input_stem>_view.html in same dir)
  --dark                       Dark background theme
```

### Output

- Single self-contained HTML file: `<input_stem>_view.html`
- Placed next to the input `.gro` file by default (or `--output` override)
- Terminal prints: atom count, stripped count, output path

### Integration with Claude Code

When the user asks me to display a `.gro` file:
1. Run `python scripts/view_gro.py <path>` via Bash
2. Open the generated HTML in Chrome using MCP browser tools (`navigate` to `file:///path/to/output.html`)
3. User sees the interactive 3D view in their browser

### Dependencies

- `numpy` (already in project)
- `plotly` (already optional in project, used by `pull_cli_view.py`)
- No new dependencies

### What This Does NOT Do

- No bond calculation (too slow for quick sanity check, not needed for scatter)
- No ribbon/surface rendering (keep it simple)
- No trajectory playback (single frame only)
- No editing capabilities
