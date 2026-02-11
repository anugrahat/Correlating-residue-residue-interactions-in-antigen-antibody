# Antibody Optimization Notes

## Goal
- Baseline system is `C1-87g7`.
- Objective is to increase antibody affinity (stronger binding) relative to this baseline.
- Mutations should be on antibody residues, not antigen residues.

## Chain / Residue Mapping (Current Omicron System)
Source files:
- `home/anugraha/omicron_pulling/bigger_box/md_om.pdb`
- `home/anugraha/omicron_pulling/bigger_box/50ns/replica_50ns_pos/heat_om.py`

Residue partitions:
- `1-195` = antigen (RBD)
- `196-302` = antibody chain 1 (light-like)
- `303-422` = antibody chain 2 (heavy-like)

Why this is confirmed:
- `group1_resids = (1, 195)` and `group2_resids = (196, 422)` in `heat_om.py`.
- Sequence check:
  - `1-195` starts with `TNLCPF...` (RBD-like motif).
  - `196-302` starts with `EIVLTQSP...` (light-chain-like motif).
  - `303-422` starts with `EVQLLESG...` (heavy-chain-like motif).

Important implication:
- Positions like `489` are antigen-side (RBD numbering), not antibody-side.
- Example: `489R` is antigen mutation, not antibody optimization.

## Regression Interpretation Rules
From:
- `home/anugraha/omicron_pulling/bigger_box/50ns/replica_50ns_pos/output_data/new_analysis/average_11_avg_f_rbd_2.py`

Pair remapping logic used in analysis:
- First residue in pair is shifted by `+332` to map into RBD numbering.
- Second residue is remapped from antibody internal index blocks:
  - `196..302 -> subtract 195`
  - `>=303 -> subtract 302`

Note:
- This remap introduces collisions (multiple original pairs mapping to the same shifted pair). Keep this in mind when ranking contacts.

Confirmed example (do not swap display unless explicitly needed):
- Original pair `157-405` means:
  - `157` = antigen-side residue (group1)
  - `405` = antibody-side residue (group2)
- With the chain-wise antibody remap above, `405 -> 103`, so:
  - `157-405` becomes `157-103`

## Mutation Scope for This Project
Only mutate antibody residues:
- Allowed range: `196-422`
- Prefer interface-focused candidates from regression/contact hotspots.

Do not mutate for antibody optimization:
- Antigen range `1-195` (RBD)
- Any RBD-numbered suggestion like 417/449/474/489 etc. unless intentionally doing antigen engineering.

## Folder Workflow Convention
For each new design, create:
- `/home/anugraha/c1_<mutation_name>/`

Recommended contents:
- system prep files (`.gro`, `.top`, `.ndx`, `.mdp`)
- pulling replica folders
- analysis outputs (`interaction_energy.csv`, `average_frequency.csv`, regression outputs)
- short run log/README for reproducibility

Examples:
- `/home/anugraha/c1_H331R/`
- `/home/anugraha/c1_S355W/`
- `/home/anugraha/c1_Y103F/`

## Immediate Next Step
- Build an antibody-only mutation shortlist from existing Elastic Net outputs by restricting candidate residues to `196-422`.

## C1 Mutation Suggestions (Engineering-Guided Panel)
Input used for prioritization:
- `elastic_net_coefficients_and_meanFreq.csv` (`Beta_x_MeanFreq`)
- Residue `103` (global `405`) is the strongest favorable hotspot by `sum(Beta_x_MeanFreq)`.
- Your observed outcomes: `Y103W` no clear gain, `Y103F` weaker binding.

Engineering interpretation:
- Position `103` is real and important, but Trp/Phe substitutions were the wrong chemistries for this site.
- Keep the aromatic platform around `102-104` active, but tune polarity and directionality at `103`.
- Build combinations only from single-mutant winners.

Numbering reminder:
- Local antibody numbering from analysis: `196..302 -> -195`, `303..422 -> -302`.
- Example: global `405 -> local 103`, global `406 -> local 104`, global `404 -> local 102`.

### GRO numbering caution (important)
- `mutant_GMX.gro` does not keep PDB chain IDs; residues are renumbered sequentially by file order.
- For current C1 build order:
  - chain `A` (antigen, residues `2..195`) becomes GRO residues `1..194`
  - chain `H` (antibody, residues `1..120`) becomes GRO residues `195..314`
  - chain `L` (antibody, residues `1..107`) becomes GRO residues `315..421`
- Therefore:
  - mutation `Y405N` (global `405` = chain `H` residue `103`) appears at `GRO residue 297` as `ASN`
  - `GRO residue 405` is a different site (`chain L residue 91`, `ARG`)
- Conclusion: checking `TYR` at `GRO 405` does not mean the `Y405N` mutation failed.

### Single mutants to test first
Tier 1 (highest priority):
1. `Y103H` (global `Y405H`)
2. `Y103N` (global `Y405N`)
3. `Y103Q` (global `Y405Q`)
4. `Y104W` (global `Y406W`)
5. `Y102W` (global `Y404W`)
6. `F92W` (global `F287W`)
7. `D70N` (global `D265N`)
8. `S67Y` (global `S262Y`)

Tier 2 (if needed after Tier 1):
1. `A51Y` (global `A246Y`)
2. `Y32W` (global `Y227W`)
3. `W94Y` (global `W289Y`)
4. `N57Y` (global `N359Y`)

Controls to keep in the same campaign:
1. WT (`Y103`)
2. `Y103W`
3. `Y103F`

### Multi-mutant panel (build only from single winners)
Start with double mutants:
1. `Y103H + Y104W`
2. `Y103N + Y102W`
3. `Y103Q + F92W`
4. `Y102W + Y104W`
5. `Y104W + F92W`

Then triple mutants:
1. `Y103H + Y104W + F92W`
2. `Y103N + Y102W + F92W`
3. `Y103Q + Y102W + Y104W`

### Execution workflow
1. Run Tier 1 singles first with identical pulling settings and replicate count vs C1 baseline.
2. Rank by consistent improvement in rupture-force profile and separation work (not one-frame peaks only).
3. Promote top 2-3 singles into double mutants; do not jump directly to 4-way combinations.
4. If all `Y103` variants fail, pivot to `102/104/92` only and continue combination design there.

## Automated Workflow

### Prerequisites

**Cluster modules:**
- `amber/24` — AmberTools25 (tleap) for rebuilding missing sidechain atoms after mutation
- `gromacs/2024` — MD engine (pdb2gmx, grompp, mdrun, trjconv, etc.)

**Force field (must exist before running):**
- `/home/anugraha/amber14sb.ff/` — native GROMACS port of Amber ff14SB (Maier et al. 2015)
- Source: `git clone https://github.com/intbio/gromacs_ff.git` then `cp -r amber14sb_parmbsc1.ff /home/anugraha/amber14sb.ff`
- Includes Na LJ parameter bugfix from intbio maintainers

**GMXLIB (set automatically by the sbatch script):**
- `export GMXLIB="/home/anugraha:/software/gromacs/2024/ucdhpc-20.04/share/gromacs/top"`
- First path: parent directory of `amber14sb.ff/`
- Second path: GROMACS default topologies (water boxes, standard ions)

**Input files:**
- Base PDB: `/home/anugraha/antibody_optimization/c1_87g7.pdb` (wildtype antibody-antigen complex)
- MDP files: placed in `inputs/` subfolder of each mutation directory
  - Copied from templates: `ions.mdp`, `step4.0_minimization.mdp`, `step4.1_equilibration.mdp`, `umb2.mdp`
  - Generated by `setup_smd_campaign.py`: `npt1.mdp` (Berendsen), `pull_npt.mdp` (Berendsen), `md10.mdp` (V-rescale)

**Python dependencies:**
- `numpy` — used by `auto_align_for_pull.py` for SVD/rotation
- `matplotlib` — used for alignment report PNG generation

**SLURM requirements:**
- Partition: `gpu-ahn`, Account: `ahnlab`
- 2 GPUs (`--gres=gpu:2`), 2 MPI ranks, 16 CPUs/task, 120h walltime

**System-specific assumptions baked into the pipeline:**
- Antigen = chain A (first 194 residues in GRO numbering)
- Antibody = chains H + L (residues 195–421 in GRO numbering)
- 6 disulfide bonds (auto-detected by pdb2gmx from S-S distance < 0.3 nm)
- Mutations restricted to antibody range (global residues 196–422)
- Pull axis = Y, pull distance = 5.0 nm (from `umb2.mdp`: rate=0.001 nm/ps, 2.5M steps, dt=0.002)
- Water model: TIP3P (consistent with ff14SB parameterization)
- Force field: Amber ff14SB via native GROMACS `pdb2gmx` (not acpype)

**Quick start:**
```bash
# Generate mutation folder(s)
python3 /home/anugraha/antibody_optimization/scripts/setup_smd_campaign.py \
    --mutations Y405N Y405Q Y406W

# Submit a job
sbatch /home/anugraha/c1_Y405N/run_pipeline.sbatch
```

### Scripts
- `/home/anugraha/antibody_optimization/scripts/setup_smd_campaign.py` -- create per-mutation folders + sbatch scripts
- `/home/anugraha/antibody_optimization/scripts/auto_align_for_pull.py` -- auto-align interface + asymmetric box sizing (called by sbatch)
- `/home/anugraha/antibody_optimization/scripts/pull_cli_view.py` -- inspect/visualize pull setups
- `/home/anugraha/antibody_optimization/scripts/compare_vs_c1.py` -- compare results vs C1 baseline
- `/home/anugraha/antibody_optimization/scripts/run_top_campaign.sh` -- batch submission helper

### Pipeline (what `run_pipeline.sbatch` does)

**Phase 1 — Normal MD in small box** (fast, fewer water molecules):
1. `tleap` (ff14SB) builds complete structure from mutant PDB (rebuilds missing sidechain atoms) -> `work/mutant_complete.pdb`
2. `gmx pdb2gmx` with native **amber14sb.ff** (`/home/anugraha/amber14sb.ff`) -> GROMACS `.gro`/`.top` with proper `#include`-based topology
   - Flags: `-ff amber14sb -water tip3p -ignh -ss` (no `-merge all`)
   - `-ignh`: strips input H, lets pdb2gmx rebuild with correct naming
   - **No `-merge all`**: chains kept as separate molecules (Protein, Protein2, Protein3). This is critical for correct PBC re-imaging — merged chains cannot be properly clustered by `trjconv -pbc cluster` since they appear as one molecule with no covalent bonds between chains.
   - `-ss`: interactive disulfide detection (piped `printf 'y\n'` for 6 SS bonds, all intra-chain)
   - `-i posre.itp`: position restraint file (bare name, then `mv posre_Protein*.itp work/`)
   - Generates separate chain ITPs: `mutant_GMX_Protein.itp` (antigen), `mutant_GMX_Protein2.itp` (heavy), `mutant_GMX_Protein3.itp` (light)
   - `GMXLIB` set to include `/home/anugraha` (where `amber14sb.ff/` lives)
3. Save clean topology (`mutant_GMX_clean.top`) before solvation modifies `[ molecules ]`
4. `editconf -d 1.2 -c -bt triclinic` -> small box for equilibration
5. `gmx solvate` + `gmx genion` (neutralize + 0.15 M NaCl)
6. EM -> NVT -> NPT -> 10 ns production MD
7. `trjconv -pbc cluster -center` (Protein/Protein/System) -> `md10_clustered.gro`, then extract protein-only -> `pull/md10_protein.gro`
   - Uses `-pbc cluster` (not `-pbc whole`) because chains are separate molecules — this correctly re-images all chains into the same periodic box near each other

**Phase 2 — Pulling in PBC-safe big box** (only for the pull step):
8. **Auto-align for pulling** (`auto_align_for_pull.py` on protein-only structure):
   - **Rotation 1 (interface normal):** Detects interface via SVD on contact atoms, rotates so interface normal is parallel to pull axis (Y)
   - **Rotation 2 (COM-COM vector):** After interface rotation, aligns the COM-COM vector between antigen and antibody to the pull axis (Y). This ensures clean separation along the pull direction even if the interface normal and COM vector differ (typically ~20° apart). Uses `rotation_matrix_from_vectors` with rotation center at midpoint of the two COMs.
   - Orients smaller group at +Y (pulled), larger at -Y (restrained)
   - Builds PBC-safe box: `max(span+2*margin, 2*(COM_sep+pull_dist)+2*margin)` along Y
   - Offsets restrained group near -Y wall, centered on X/Z
   - Falls back to COM-COM vector if interface detection fails
   - Writes `pull/aligned.gro` + `pull/alignment_report.png`
9. Copy clean topology + chain ITPs + posre files to `pull/` -> re-solvate in big box + re-add ions
10. Short re-equilibration: EM + 100ps NPT with position restraints (relax new water shell)
11. Build pull index (`chauA`/`rest`) + pulling run with `umb2.mdp`

**Why two phases**: The big PBC-safe box (~21 nm along Y) is only needed for pulling. Running EM/NVT/NPT/md10 in a small box (~12 nm) uses far fewer water molecules and runs much faster. The protein coordinates from md10 are preserved — only orientation, position, and box change.

**Checkpointing**: Every step checks if its output `.gro` already exists before running. If a job is resubmitted (e.g., after a crash on the pull step), it skips all completed steps automatically.

**GPU flags**:
- `GPU_FLAGS="-nb gpu -bonded gpu -pme gpu"` — used for equilibration steps (NVT, NPT, md10, pull NPT) which use V-rescale/Berendsen thermostats
- `PULL_GPU_FLAGS="-nb gpu -bonded gpu -pme gpu"` — used for the final pull mdrun. Does **not** include `-update gpu` because `umb2.mdp` uses Nose-Hoover thermostat, which is incompatible with GPU-accelerated update in GROMACS 2024.
- `EM_FLAGS="-nb gpu"` — EM only offloads nonbonded
- Multi-GPU: 2 MPI ranks (`srun $GMX mdrun -npme 1`), 1 rank dedicated to PME. Serial GROMACS tools use `srun --ntasks=1`.

### Box sizing logic (PBC-safe)
From `umb2.mdp`: rate=0.001 nm/ps, 2.5M steps, dt=0.002 -> **5.0 nm total pull distance**.

**Critical PBC constraint** (Lemkul GROMACS pulling tutorial):
- GROMACS calculates distances with periodicity: if you pull over > box/2, the periodic image distance becomes the reference
- Therefore: `box_Y / 2 > initial_COM_separation + pull_distance` at ALL times
- This is non-negotiable — violating it gives wrong forces and wrong results

The script takes the **maximum** of two constraints on the pull axis:
1. Water padding: `span + 2*margin`
2. PBC minimum image: `2*(COM_sep + pull_dist) + 2*margin`

| Axis | Formula | Example (Y405N) |
|------|---------|-----------------|
| Y (pull) | max(span + 2*2.0, 2*(COM_sep + 5.0) + 2*2.0) | 21.0 nm |
| X (perp) | span + 2*1.5 | 8.8 nm |
| Z (perp) | span + 2*1.5 | 7.4 nm |

### Positioning within the box
- **X, Z (perpendicular)**: protein centered (equal water padding both sides — PBC requires this)
- **Y (pull axis)**: protein offset toward -Y wall
  - Restrained (larger) group: placed near -Y wall with `margin_pull_nm` padding (2.0 nm)
  - Pulled (smaller) group: gets all remaining room at +Y
  - This is correct because the restrained group doesn't move — it doesn't need equal padding on both sides
  - The PBC box SIZE still satisfies minimum image convention regardless of offset

### Verification checklist (run after auto-align)
Always verify these before submitting a pull job:
1. `box_Y / 2 > COM_sep + pull_distance` (PBC minimum image)
2. Restrained group has >= 1.5 nm gap to -Y wall (water padding)
3. Pulled group has >= `pull_distance + margin` gap to +Y wall after pull
4. X/Z walls have >= 1.0 nm gap (water padding)
5. Smaller group is at +Y, larger at -Y (orientation for pull direction)
6. Interface normal rotated toward Y axis (rotation 1)
7. COM-COM vector aligned to Y axis (rotation 2, angle after ~0°)

Example verified output (Y405N):
```
Box Y:           21.03 nm   |  box/2 = 10.52 nm
rest  (restrained): Y = 2.00 - 6.83  COM=4.62   (2.0 nm to -Y wall)
chauA (pulled):     Y = 5.65 - 11.12  COM=8.13
After 5nm pull:     pulled max → 16.12  (4.9 nm still to +Y wall)
PBC:  COM_sep+pull = 8.52 < box/2 = 10.52  ✓
```

### Force field
- **Native GROMACS amber14sb** force field: `/home/anugraha/amber14sb.ff/`
  - Source: [intbio/gromacs_ff](https://github.com/intbio/gromacs_ff) (amber14sb_parmbsc1 port)
  - Includes known Na LJ parameter bugfix
  - ff14SB protein parameters (Maier et al. 2015) + parmbsc1 nucleic acid corrections
- **Previous workflow** (deprecated): `tleap` -> `acpype` -> inline parameters + `amber99sb-ildn.ff/spce.itp`
- **Current workflow**: `tleap` (structure building only) -> `gmx pdb2gmx -ff amber14sb` (topology)
  - Advantages: native `#include`-based topology, proper posre.itp, no acpype/fix_acpype_top.py needed
  - `GMXLIB` must include `/home/anugraha` so GROMACS can find `amber14sb.ff/`
  - The sbatch script sets: `export GMXLIB="/home/anugraha:/software/gromacs/2024/ucdhpc-20.04/share/gromacs/top"`
- Water model: TIP3P (consistent with ff14SB parameterization)
- Note: baseline CHARMM-GUI systems use CHARMM36; keep this in mind when comparing absolute values across mixed force-field campaigns
- Pull groups: `chauA` (antigen, first 194 residues) vs `rest` (antibody, remaining residues)

### Standalone auto-align usage (any system)
```bash
python3 /home/anugraha/antibody_optimization/scripts/auto_align_for_pull.py \
    --gro <vacuum.gro> --out <aligned.gro> \
    --antigen-res-count <N> --pull-dim Y \
    --mdp <umb2.mdp> \
    --margin-pull-nm 2.0 --margin-perp-nm 1.5 \
    --report-png <report.png> \
    --fallback-to-com --orient-smaller-positive
```

Key flags:
- `--antigen-res-count N`: first N protein residues = group 1 (chauA)
- `--orient-smaller-positive`: smaller group at +pull-axis (pulled), larger at - (restrained). On by default.
- `--fallback-to-com`: use COM-COM vector if <3 interface contacts found
- `--mdp`: auto-detect pull distance from MDP (rate * nsteps * dt)
- `--pull-distance-nm`: override pull distance manually
- `--skip-threshold-deg 5`: skip rotation if already within 5 degrees (applies to each rotation independently)

**Two-step rotation logic:**
1. **Interface-normal rotation**: SVD on interface contact atoms detects the plane where antigen and antibody meet. The normal of this plane is rotated to align with the pull axis (Y). This ensures the interface is perpendicular to the pull direction — forces are applied evenly across the binding surface rather than peeling from one edge.
2. **COM-COM rotation**: After the interface rotation, the center-of-mass vector between antigen and antibody may still be off-axis (typically ~20° for globular complexes where the interface plane normal and COM vector don't coincide). A second rotation aligns this COM vector to Y. This ensures clean separation along the pull axis.
   - Rotation center: midpoint of the two group COMs
   - If COM is already within `--skip-threshold-deg` of Y, the second rotation is skipped
   - Trade-off: the interface normal ends up ~20° off Y, but the separation trajectory is along Y. In practice, this produces cleaner pull force profiles because the groups separate along the box's long axis without lateral drift.

Top-mutation folders already created:
- `/home/anugraha/c1_Y405N`
- `/home/anugraha/c1_Y405Q`
- `/home/anugraha/c1_Y406W`
- `/home/anugraha/c1_Y404W`
- `/home/anugraha/c1_F287W`
- `/home/anugraha/c1_D265N`
- `/home/anugraha/c1_S262Y`

Submit jobs:
1. `sbatch /home/anugraha/c1_Y405N/run_pipeline.sbatch`
2. `sbatch /home/anugraha/c1_Y405Q/run_pipeline.sbatch`
3. `sbatch /home/anugraha/c1_Y406W/run_pipeline.sbatch`
4. `sbatch /home/anugraha/c1_Y404W/run_pipeline.sbatch`
5. `sbatch /home/anugraha/c1_F287W/run_pipeline.sbatch`
6. `sbatch /home/anugraha/c1_D265N/run_pipeline.sbatch`
7. `sbatch /home/anugraha/c1_S262Y/run_pipeline.sbatch`

After runs complete, compare vs C1 baseline:
1. `python /home/anugraha/antibody_optimization/scripts/compare_vs_c1.py --mutations Y405N Y405Q Y406W Y404W F287W D265N S262Y`

## CLI Viewer (No-Solvent, Terminal-First)
Script:
- `/home/anugraha/antibody_optimization/scripts/pull_cli_view.py`

Purpose:
1. Read `.gro` + `.ndx`
2. Strip solvent for inspection
3. Report box/pull-group geometry in terminal
4. Write quick visual artifacts for fast review

Run command:
1. `python3 /home/anugraha/antibody_optimization/scripts/pull_cli_view.py --gro /home/anugraha/c1_Y405N/work/md10.gro --ndx /home/anugraha/c1_Y405N/work/rbd_2_reg.ndx --group1 chauA --group2 rest --pull-dim Y --outdir /home/anugraha/c1_Y405N/view_cli`
2. Optional popup/open attempt: add `--open-html`
3. Optional orientation enforcement: add `--auto-align-to-pull`
4. Interface-perpendicular mode: add `--auto-align-to-pull --align-interface-normal --interface-cutoff-nm 0.6`
5. Emit equivalent `gmx editconf -rotate` command: add `--emit-gmx-rotate-cmd`

Outputs:
- `/home/anugraha/c1_Y405N/view_cli/solvent_free_view.pdb`
- `/home/anugraha/c1_Y405N/view_cli/pull_groups_box_projections.png`
- `/home/anugraha/c1_Y405N/view_cli/pull_groups_box_interactive.html` (rotatable 3D viewer)

## Full Auto-Prep For Any Two Proteins
This now works as a one-command autoprep for any pair of proteins, regardless of starting orientation, as long as you provide:
- a `.gro` with both proteins
- an `.ndx` containing the two pull groups (`group1`, `group2`)

What it does in one run:
1. Centers groups in box
2. Finds interface contacts and fits interface plane
3. Rotates so interface normal is parallel to pull axis (interface perpendicular to pull)
4. Reboxes with configurable margins
5. Writes transformed full-system `.gro` for next gromacs steps
6. Writes `pull_groups.ndx`, summary report, static and interactive viewers
7. Emits equivalent `gmx editconf -rotate` command

Generic command template:
1. `python3 /home/anugraha/antibody_optimization/scripts/pull_cli_view.py --gro <system.gro> --ndx <groups.ndx> --group1 <proteinA_group> --group2 <proteinB_group> --pull-dim Y --auto-align-to-pull --align-interface-normal --interface-cutoff-nm 0.6 --rebox-mode axis --rebox-margin-nm 1.5 --rebox-pull-margin-nm 2.5 --emit-gmx-rotate-cmd --outdir <autoprep_outdir>`

Main output files in `<autoprep_outdir>`:
- `aligned_for_pull.gro`
- `pull_groups.ndx`
- `summary.txt`
- `pull_groups_box_interactive.html`
- `gmx_editconf_rotate_command.txt`

Current `c1_Y405N` result snapshot:
- Box (nm): `(7.967, 9.274, 10.346)`
- Pull-axis `Y` COM delta (nm): `2.735`
- Minimum wall gap:
  - `chauA`: `(1.574, 0.836, 1.036)`
  - `rest`:  `(1.578, 0.972, 1.291)`
- Recommendation from CLI tool: increase box (tightest gap `0.836 nm` along `Y`, below `1.0 nm` warning threshold)

Interpretation for quick decision:
- The limiting axis is `Y` (pull axis), not `X` or `Z`.
- For safer pulling, prioritize adding space along `Y` (or reorient before boxing so pull direction has more margin).

Auto-align check (`--auto-align-to-pull`) on current `Y405N`:
- COM-vector angle to pull axis: `43.80 deg -> 0.00 deg` (aligned)
- After alignment, closest wall gap becomes ~`0.004 nm` in `Y`, which confirms the current box is too tight for that orientation.
- Practical implication: if you enforce perpendicular interface / pull-axis alignment, increase `Y` box size before pull setup.

Interface-normal align check (`--auto-align-to-pull --align-interface-normal`) on current `Y405N`:
- Interface contacts used (`0.6 nm` cutoff): `chauA=315 atoms`, `rest=291 atoms`
- Interface-normal angle to pull axis: `11.87 deg -> 0.00 deg` (perpendicular interface enforced)
- COM vector after this mode is not forced to axis, which is expected.
- Closest wall gap still tight in `Y` (`0.624 nm`), so increase `Y` box margin before production pull.

Equivalent `gmx editconf` rotation for the same interface-perpendicular alignment:
- `-rotate 5.280957 -0.490218 10.599453`
- Vector-target decomposition error check: `0.000000 deg`
- Generated command file:
  - `/home/anugraha/c1_Y405N/view_cli_interface_align/gmx_editconf_rotate_command.txt`
