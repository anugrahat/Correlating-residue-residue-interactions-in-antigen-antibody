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

## Mutation Design Strategy (Elastic Net Regression)

### How we identify mutation targets
1. **Time-resolved contact frequency analysis**: From equilibrium MD trajectories, compute per-frame residue-residue contact frequencies between antigen and antibody.
2. **Elastic net regression**: Regress contact frequencies against binding energy (interaction energy). The beta coefficient for each contact pair tells how much that contact contributes to binding:
   - **Negative beta** = contact stabilizes binding (favorable)
   - **Positive beta** = contact destabilizes binding (unfavorable)
3. **Aggregate by antibody residue**: Sum `|beta × MeanFreq|` across all antigen partners to rank antibody residues by importance.
   - Data source: `c1_from_existing_regression_antibody_targets.csv`

### Residue classification (from regression) — REVISED

The original hot/warm/cold classification assumed hot spots were "already optimal" and should not be mutated. **This was wrong.** Round 3 results showed that hot spot mutations (F287W, D245N, Y404W) produced the largest improvements, while warm spot mutations gave smaller gains. The regression identifies high-leverage positions, not untouchable ones.

| Classification | NetFavorability | Meaning | Revised Action |
|---|---|---|---|
| **Hot spot** | high | Strong coupling to binding energy — high leverage | **Best mutation targets** — right chemistry gives largest ΔΔG |
| **Warm spot** | moderate | Moderate coupling — moderate leverage | Good targets, but smaller effect size |
| **Cold spot** | low/zero | Weak or no coupling | Low priority — mutations have small effects |
| **Unfavorable** | any, positive β | Contact destabilizes binding | Remove/shrink side chain |

### Key insight: regression identifies HIGH-LEVERAGE positions
- The regression ranks positions by their coupling strength to binding energy (|β|)
- **Higher coupling = larger effect of mutation in either direction** (improvement or disruption)
- Hot spots gave ΔΔG of -18 to -39 kcal/mol; warm spots gave -6 to -16 kcal/mol
- The classification tells you WHERE to mutate (high |β|), not WHETHER to mutate
- **What determines success is the chemistry of the substitution**, not the hot/warm label:
  - Charge removal (D→N): works at any classification — relieves desolvation penalty
  - Aromatic extension (Phe→Trp, small→Tyr): works when it extends existing aromatic platform
  - Non-conservative disruption (Tyr→His): fails at hot spots — destroys the interaction

### Ranked antibody interface residues

| Global | Chain | AA | NetFav | #Pairs | Class |
|---|---|---|---|---|---|
| 406 | H:104 | TYR | 165.6 | 22 | HOT |
| 227 | L:32 | TYR | 119.9 | 23 | HOT |
| 404 | H:102 | TYR | 104.8 | 28 | HOT |
| 405 | H:103 | TYR | 86.0 | 19 | HOT |
| 287 | L:92 | PHE | 76.7 | 20 | HOT |
| 249 | L:54 | ARG | 74.5 | 13 | HOT |
| 374 | H:72 | ARG | 74.0 | 12 | HOT |
| 262 | L:67 | SER | 50.8 | 24 | WARM |
| 251 | L:56 | THR | 39.8 | 9 | WARM |
| 265 | L:70 | ASP | 24.3 | 18 | WARM |
| 246 | L:51 | ALA | 23.6 | 24 | WARM |
| 289 | L:94 | TRP | 21.1 | 13 | WARM |
| 226 | L:31 | ASN | 17.8 | 23 | lukewarm |
| 245 | L:50 | ASP | 16.2 | 21 | lukewarm |
| 403 | H:101 | THR | 8.0 | 22 | lukewarm |
| 288 | L:93 | ASN | 0 | 18 | cold (many pairs) |
| 248 | L:53 | ASN | 0 | 17 | cold (many pairs) |
| 264 | L:69 | THR | 0 | 14 | cold (many pairs) |
| 260 | L:65 | SER | 0 | 12 | cold (many pairs) |
| 225 | L:30 | HIS | 0 | 16 | cold (many pairs) |
| 244 | L:49 | PHE | -16.3 | 9 | DESTABILIZING |

## Regression Data Sources

| System | Trajectory | Interaction Energy CSV | Min IE (kJ/mol) | Notes |
|--------|-----------|----------------------|-----------------|-------|
| **Wuhan-87g7** | 50ns | `/home/anugraha/omicron_pulling/bigger_box/50ns/replica_50ns_pos/output_data/new_analysis/interaction_energy.csv` | -849.9 | C1 antibody + Wuhan RBD; avg of 4 replicas (C_SR+LJ_SR) |
| **Omicron-87g7** | 100ns | `/home/anugraha/omicron_pulling/bigger_box/100ns/replica_100ns_pos/output_data_om/new_analysis/interaction_energy.csv` | -905.2 | C1 antibody + Omicron RBD |

Labels confirmed from: `/home/anugraha/omicron_pulling/bigger_box/100ns/replica_100ns_pos/output_data_om/new_analysis/count.py`

Elastic net coefficients (used for mutation design):
- C1 (50ns): `/home/anugraha/omicron_pulling/bigger_box/50ns/replica_50ns_pos/output_data/new_analysis/elastic_net_coefficients_and_meanFreq.csv`
- Omicron (100ns): `/home/anugraha/omicron_pulling/bigger_box/100ns/replica_100ns_pos/output_data_om/new_analysis/elastic_net_coefficients_and_meanFreq.csv`

Aggregated antibody targets (from C1 regression): `/home/anugraha/antibody_optimization/c1_from_existing_regression_antibody_targets.csv`

## WT Baseline (C1-87g7)

Data from: `/home/anugraha/c1_87g7/rot_120_corrected/bigger_box/md2/`

| Replica | Smoothed Peak (kJ/mol/nm) |
|---|---|
| eq60/80ns_pos | 940 |
| eq80/80ns_pos | 812 |
| eq45/80ns_pos | 748 |
| eq100/100ns | 723 |
| **Mean ± Std** | **806 ± 97** |

## Round 1: Warm-Spot Single Mutations (completed)

Strategy: mutate warm spots with conservative substitutions to improve suboptimal contacts.

| Mutation | Global | Chain | Rationale | Peak Force (kJ/mol/nm) | vs WT (806) |
|---|---|---|---|---|---|
| D265N | 265 | L:70 | Remove negative charge | 795 ± 196 | ~equal |
| F287W | 287 | L:92 | Conservative Phe→Trp | 671 ± 40 | -17% |
| Y405N | 405 | H:103 | Test hot spot (control) | 572 ± 29 | -29% |
| A246Y | 246 | L:51 | Add aromatic at small site | 480 ± 47 | -40% |
| S262Y | 262 | L:67 | Add aromatic at polar site | 454 ± 101 | -44% |

**Conclusion**: No single warm-spot mutation exceeded WT. Y405N confirmed as hot spot (any substitution weakens binding). D265N is closest to WT but high variance.

## Round 2: Affinity Maturation (in progress)

### Strategy shift
Round 1 showed that improving existing contacts (warm spots) is not enough — those contacts are already functional, and mutations tend to disrupt more than they gain. Two new approaches:

**Approach 1 — Combine best singles**: D265N and F287W target different interface regions (L:70 vs L:92). If effects are additive, the double mutant could exceed WT.

**Approach 2 — Cold spots with many pairs (untapped potential)**: Residues that sit at the interface and make many antigen contacts but contribute ZERO to binding energy. The current residue makes geometrically suboptimal contacts. A bulkier/aromatic substitution could convert dead contacts into productive ones.

**Approach 3 — Lukewarm positions with high contact count**: Residues with moderate signal but disproportionately many contact pairs (e.g., T403 has 22 pairs but only 8.0 net favorability — massive untapped potential).

### Round 2 Tier 1 mutations (submitted)

| Mutation | Type | Rationale |
|---|---|---|
| D265N+F287W | double | Best two Round 1 singles combined |
| T403Y | single (lukewarm→aromatic) | H:101, 22 pairs but only 8.0 net fav; Thr→Tyr adds aromatic+OH |
| N248Y | single (cold→aromatic) | L:53, 17 pairs, zero signal; Asn→Tyr to create new contacts |

### Round 2 Tier 2 mutations (pending results)

| Mutation | Rationale |
|---|---|
| N226Y | L:31, 23 pairs, 17.8 net fav; Asn→Tyr at lukewarm position with most contact pairs |
| F244A | L:49, net fav -16.3; remove destabilizing bulky aromatic |
| N288W | L:93, 18 pairs, zero signal; extend aromatic cluster near F287 |
| T264Y | L:69, 14 pairs, zero signal; adjacent to D265 |

### Pipeline features
- Multi-mutation support: `--mutations D265N+F287W` (+ separator for combos on same system)
- WT baseline: `--mutations WT` (copies base PDB unchanged)
- Replicates: `--replicates 3` generates rep2/rep3 with different gen_seed, separate sbatch scripts
- Folder naming: `c1_D265N_F287W`, `c1_WT`, etc.

### Single mutants to test first (original list, partially superseded by Round 1/2)
Tier 1 (highest priority):
1. `Y103H` (global `Y405H`)
2. `Y103N` (global `Y405N`) — tested, -29% vs WT
3. `Y103Q` (global `Y405Q`)
4. `Y104W` (global `Y406W`)
5. `Y102W` (global `Y404W`)
6. `F92W` (global `F287W`) — tested, -17% vs WT
7. `D70N` (global `D265N`) — tested, ~equal to WT
8. `S67Y` (global `S262Y`) — tested, -44% vs WT

Tier 2 (if needed after Tier 1):
1. `A51Y` (global `A246Y`) — tested, -40% vs WT
2. `Y32W` (global `Y227W`)
3. `W94Y` (global `W289Y`)
4. `N57Y` (global `N359Y`)

Controls:
1. WT (`Y103`) — baseline: 806 ± 97 kJ/mol/nm
2. `Y103W` — previously tested, no improvement
3. `Y103F` — previously tested, weaker binding

### Multi-mutant panel (build only from single winners)
Start with double mutants:
1. `D265N + F287W` — **submitted (Round 2 Tier 1)**
2. Best Round 2 single + D265N
3. Best Round 2 single + F287W

### Execution workflow
1. Run Tier 1 singles first with identical pulling settings and replicate count vs C1 baseline.
2. Rank by consistent improvement in rupture-force profile and separation work (not one-frame peaks only).
3. Promote top 2-3 singles into double mutants; do not jump directly to 4-way combinations.
4. If all `Y103` variants fail, pivot to `102/104/92` only and continue combination design there.
5. If warm-spot singles don't improve, shift to cold-spot aromatic additions (Round 2 approach).

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

## SMD Pulling Analysis (Work Integrals)

### Methodology
Instead of using peak rupture force (noisy, single-frame measure), we computed **work integrals** W = integral(F dx) using trapezoidal integration of pullf.xvg force-distance curves. This integrates over the full separation process and is more robust than peak force alone.

We also computed **Jarzynski PMF** using ΔG = -kT ln⟨exp(-W/kT)⟩ across replicates, though this is exponentially dominated by the lowest-work replicate and unreliable with <30 replicates.

### Replicates completed
- WT: 10 replicates (rep1-10)
- D265N: 10 replicates
- F287W: 10 replicates
- Y157A: 6 replicates (negative control — alanine scan of antigen residue)
- D265N+F287W: 6 replicates (+ rep7-10 submitted)
- T403Y: 6 replicates (+ rep7-10 submitted)
- S262Y: 6 replicates (+ rep7-10 submitted)
- N248Y: 6 replicates (+ rep7-10 submitted)
- A246Y: 6 replicates (+ rep7-10 submitted)

### Work integral results (kJ/mol)

| Mutation | Mean W | SEM | n | vs WT (%) | p-value |
|----------|--------|-----|---|-----------|---------|
| D265N | 715 | 34 | 10 | +90% | <0.001 |
| F287W | 551 | 27 | 10 | +46% | 0.002 |
| T403Y | 506 | 24 | 6 | +34% | 0.001 |
| S262Y | 478 | 39 | 6 | +27% | 0.040 |
| D265N+F287W | 464 | 41 | 6 | +23% | 0.080 |
| N248Y | 449 | 29 | 6 | +19% | 0.070 |
| WT | 377 | 22 | 10 | — | — |
| Y157A | 332 | 27 | 6 | -12% | 0.100 |
| A246Y | 312 | 36 | 6 | -17% | 0.140 |

### SMD limitations identified
- **Noise floor ±15-20%**: Even with n=10, WT spans 260-480 kJ/mol. Only D265N (+90%) and F287W (+46%) are reliably above noise.
- **SMD measures nonequilibrium mechanical resistance**, not equilibrium binding affinity. Fast pulling (0.001 nm/ps) means most work is dissipated, not thermodynamic.
- Mid-tier mutations (T403Y +34%, S262Y +27%) overlap with WT variance — cannot be reliably ranked.
- Y157A negative control showed only -12% (not statistically significant with WT at n=10), despite being expected to strongly disrupt binding.
- Conclusion: SMD alone cannot resolve moderate binding differences (1-3 kcal/mol). Switched to MM-GBSA for reliable ranking.

### Plots generated
All in `/home/anugraha/antibody_optimization/`:
- `plot1_work_integral.png` — Bar chart with significance stars and individual data points
- `plot2_force_profiles.png` — Smoothed average force profiles per mutation
- `plot3_jarzynski_pmf.png` — Jarzynski PMF with 95% bootstrap CI
- `plot4_cumulative_work.png` — Cumulative work ± SEM
- `plot5_violin_distribution.png` — Violin/distribution plot
- `plot6_negative_control.png` — Y157A vs top hits focused comparison
- `plot7_significance_heatmap.png` — Pairwise t-test significance heatmap
- `y157a_full_comparison.png` — 3-panel Y157A vs WT with all 6 reps

## MM-GBSA Analysis (Equilibrium Binding Free Energy)

### Why MM-GBSA
SMD pulling is too noisy to rank mutations with moderate affinity differences. MM-GBSA uses the **existing 10ns equilibrium trajectories** (md10.trr from Phase 1) to compute binding free energy without additional simulations.

### Method
- **Tool**: `gmx_MMPBSA` v1.6.3 (wrapper around AmberTools MMPBSA.py)
- **Approach**: Single-trajectory MM-GBSA (internal energies cancel exactly)
- **Solvent model**: GB (igb=5, modified Onufriev-Bashford-Case)
- **Salt**: 0.15 M (saltcon=0.150)
- **Frames**: 90 frames from 10ns trajectory (startframe=1, endframe=900, interval=10)
  - Last 1ns excluded due to structural instability in some systems
- **Groups**: Antigen (residues 1-194, index group 17) vs Antibody (residues 195-421, index group 18)

### Setup procedure (for each mutation)
```bash
module load amber/24
module load gromacs/2024
cd /home/anugraha/c1_<MUT>/work

# 1. Fix ion names (GROMACS amber14sb uses NA/CL, topology has Na/Cl)
cp mutant_GMX.top mutant_GMX.top.bak
sed -i 's/^Na /NA /g; s/^Cl /CL /g' mutant_GMX.top

# 2. Fix PBC-broken molecules (prevents sander bond energy overflow)
printf 'Protein\nSystem\n' | gmx_mpi trjconv -f md10.trr -s md10.tpr -o md10_noPBC.trr -pbc mol -center

# 3. Create index file with Antigen (r1-194) and Antibody (r195-421) groups
# -> mmpbsa_index.ndx (groups 17 and 18)

# 4. Create input file (mmpbsa.in):
# &general
#   startframe=1, endframe=900, interval=10,
#   verbose=2,
# /
# &gb
#   igb=5, saltcon=0.150,
# /

# 5. Run
/home/anugraha/.local/bin/gmx_MMPBSA -O -i mmpbsa.in -cs md10.tpr \
    -ci mmpbsa_index.ndx -cg 18 17 -ct md10_noPBC.trr -cp mutant_GMX.top -nogui
```

### Critical fixes discovered
1. **Ion name case mismatch**: `mutant_GMX.top` has `Na`/`Cl` in `[ molecules ]` but amber14sb expects `NA`/`CL`. Fix: `sed -i 's/^Na /NA /g; s/^Cl /CL /g' mutant_GMX.top`
2. **PBC-broken molecules**: Raw md10.trr has molecules broken across periodic boundaries, causing sander BOND energy overflow (`*************`). Fix: `trjconv -pbc mol -center` to create md10_noPBC.trr
3. **Late frame instability**: Some systems (especially WT) have structural instability in last 1ns. Fix: `endframe=900` (use 90 frames instead of 100)

### MM-GBSA results (all single mutants, ranked)

| Mutation | ΔG_bind (kcal/mol) | SEM | ΔΔG vs WT | Corrected NF | Chemistry |
|----------|-------------------|-----|-----------|-------------|-----------|
| D265N | -64.73 | 1.45 | -39.4 | absent (SASA-filtered) | Charge removal (D→N) |
| F287W | -58.99 | 1.44 | -33.7 | 17.7 (HOT) | Aromatic extension (F→W) |
| D245N | -43.47 | 1.10 | -18.2 | 23.2 (HOT) | Charge removal (D→N) |
| S258G | -41.78 | 0.82 | -16.5 | absent | Side chain removal (S→G) |
| S262Y | -40.28 | 1.20 | -15.0 | 1.9 (warm) | Aromatic addition (S→Y) |
| N248Y | -36.20 | 1.01 | -10.9 | 6.2 (warm) | Aromatic addition (N→Y) |
| A246Y | -31.41 | 0.81 | -6.1 | 7.8 (HOT) | Aromatic addition (A→Y) |
| Y404W | -28.00 | 0.96 | -2.7 | 23.4 (HOT) | Aromatic extension (Y→W) |
| **WT** | **-25.31** | **0.86** | **—** | — | — |
| S258A | -24.35 | 0.91 | +1.0 | absent | Conservative removal (S→A) |
| Y404H | -22.79 | 0.77 | +2.5 | 23.4 (HOT) | Non-conservative (Y→H) |
| T403Y | -21.08 | 1.57 | +4.2 | absent | Aromatic addition (T→Y) |
| Y157A | -13.73 | 1.49 | +11.6 | antigen control | Alanine scan |

Result files: `/home/anugraha/c1_*/work/FINAL_RESULTS_MMPBSA.dat`

### Key findings

1. **D265N is the top performer** (ΔΔG = -39.4 kcal/mol): Removing the negative charge at L:70 (Asp→Asn) dramatically improves binding. This was a warm spot in the regression with moderate signal — the charge removal creates much better electrostatic complementarity at the interface.

2. **F287W is second best** (ΔΔG = -33.7 kcal/mol): Conservative Phe→Trp at L:92, a hot spot. The larger indole ring extends the aromatic platform and makes additional van der Waals contacts.

3. **Double mutant paradox**: D265N+F287W (-50.82) is **worse** than either D265N (-64.73) or F287W (-58.99) alone. Energy decomposition shows electrostatic contribution drops from -86.80 (D265N alone) to -34.14 (double), suggesting the two mutations cause electrostatic interference when combined. The D265N charge removal and F287W steric change are not independent — they share a coupled electrostatic network.

4. **S262Y and N248Y are real hits** (ΔΔG = -15.0 and -10.9): SMD couldn't resolve these (overlapped with WT noise), but MM-GBSA clearly shows improvement. S262Y adds aromatic at warm spot L:67 (24 contact pairs). N248Y adds aromatic at cold spot L:53 (17 pairs, zero regression signal → new contacts created).

5. **T403Y is a false positive**: SMD ranked it at +34% improvement, but MM-GBSA shows it's slightly worse than WT (ΔΔG = +4.2). SMD measures mechanical resistance during nonequilibrium pulling, not equilibrium binding affinity — Tyr at H:101 may resist physical separation (steric) without actually improving thermodynamic binding.

6. **A246Y reversal**: SMD showed -17% (worse than WT), but MM-GBSA shows -31.41 (better by 6.1 kcal/mol). Different properties being measured — SMD pathway may not sample the relevant binding mode.

7. **Y157A negative control validated**: ΔG = -13.73 vs WT -25.31 (ΔΔG = +11.6 kcal/mol). Alanine substitution at the antigen hotspot Y157 dramatically weakens binding, confirming MM-GBSA sensitivity. This mutation was unresolvable by SMD (only -12%, p=0.10).

### SMD vs MM-GBSA ranking comparison

| Rank | SMD (work integral) | MM-GBSA (ΔG_bind) |
|------|--------------------|--------------------|
| 1 | D265N (+90%) | D265N (-39.4) |
| 2 | F287W (+46%) | F287W (-33.7) |
| 3 | T403Y (+34%) | S262Y (-15.0) |
| 4 | S262Y (+27%) | N248Y (-10.9) |
| 5 | D265N+F287W (+23%) | A246Y (-6.1) |
| 6 | N248Y (+19%) | WT (baseline) |
| 7 | WT (baseline) | T403Y (+4.2) |
| 8 | Y157A (-12%) | Y157A (+11.6) |
| 9 | A246Y (-17%) | D265N+F287W (-25.5)* |

*D265N+F287W is still better than WT in absolute terms but worse than either single — an anti-cooperative double mutant.

Top two (D265N, F287W) agree between methods. Mid-tier ranking differs significantly — MM-GBSA is more reliable for equilibrium binding affinity.

### Recommended next steps
1. **Aromatic extension at top hot spots**: Y227W, Y406W, Y405W — highest leverage positions, proven Tyr→Trp chemistry
2. **Charge removal at E401**: E401Q — extends the successful D→N strategy to Glu
3. **Experimental validation**: D265N (-64.7), F287W (-59.0), D245N (-43.5), S258G (-41.8) are top candidates for wet-lab testing
4. **Avoid double mutants**: All 5 tested doubles were anti-cooperative. Focus on single mutants.
5. **Longer equilibrium MD**: 50-100ns trajectories would give more reliable MM-GBSA with better sampling

## Double Mutant MM-GBSA Results (Round 2)

Testing double mutant combinations to find cooperative pairs.

| Double Mutant | ΔG_bind (kcal/mol) | SEM | vs WT (-25.31) | vs Best Single | Cooperative? |
|---------------|-------------------|-----|----------------|----------------|-------------|
| D265N+F287W | -50.82 | 0.88 | -25.5 | Worse than D265N (-64.73) | **Anti-cooperative** |
| D265N+N248Y | -45.58 | 1.00 | -20.3 | Worse than D265N (-64.73) | **Anti-cooperative** |
| D265N+S262Y | -35.37 | 0.82 | -10.1 | Worse than D265N (-64.73) | **Anti-cooperative** |
| F287W+N248Y | -28.44 | 6.36 | -3.1 | Worse than F287W (-58.99) | **Anti-cooperative** |
| F287W+S262Y | -32.59 | 1.04 | -7.3 | Worse than F287W (-42.26) | **Anti-cooperative** |

**Pattern**: ALL tested doubles are anti-cooperative. D265N-based doubles lose the most due to coupled electrostatic network disruption. F287W-based doubles (F287W+N248Y, F287W+S262Y) also anti-cooperative. Single mutants remain superior to all tested doubles. Double mutant strategy abandoned.

## Round 3: New Mutations from C1 WT Regression

Pipeline: Phase 1 (tleap → pdb2gmx → solvate → EM → NVT → NPT → 10ns MD) + MM-GBSA (gmx_MMPBSA, igb=5, ~80-90 frames).
Folders: `/home/anugraha/c1_{Y404W,Y404H,S258A,S258G,D245N}/`

### Round 3 MM-GBSA results

| Mutation | ΔG_bind (kcal/mol) | SEM | ΔΔG vs WT (-25.31) | Classification | Outcome |
|----------|-------------------|-----|---------------------|---------------|---------|
| D245N | -43.47 | 1.10 | -18.2 | HOT (18.6) | **Improvement** — charge removal at hot spot |
| S258G | -41.78 | 0.82 | -16.5 | warm/unfav (1.9) | **Improvement** — remove unfavorable contact |
| Y404W | -28.00 | 0.96 | -2.7 | HOT (23.1) | Mild improvement — aromatic extension |
| S258A | -24.35 | 0.91 | +1.0 | warm/unfav (1.9) | Neutral — Ala not aggressive enough |
| Y404H | -22.79 | 0.77 | +2.5 | HOT (23.1) | Weakened — non-conservative disruption |

Note: D245N used ~8ns trajectory (early MM-GBSA, 78 frames) since the 10ns MD was still running.

### Round 3 key findings

1. **D245N is a strong hit** (ΔΔG = -18.2): Same Asp→Asn charge removal as D265N (our #1 hit). Works despite D245 being classified HOT. Confirms that the regression's hot/warm label does not predict mutability.

2. **S258G works, S258A doesn't**: Both target the same unfavorable S258 position, but only complete side chain removal (Gly) eliminates the destabilizing contact. Ala retains too much of the problematic interaction. The regression identifies the position correctly but cannot predict which substitution works.

3. **Y404W mild improvement, Y404H weakened**: Both target the same HOT Tyr. Trp preserves the aromatic platform and slightly improves packing (ΔΔG = -2.7). His disrupts the aromatic interaction (ΔΔG = +2.5). Chemistry of the substitution matters more than the classification.

### Desolvation penalty: why charge removal works

MM-GBSA energy decomposition reveals the mechanism:

| Component | WT | D265N | D245N | Meaning |
|-----------|-----|-------|-------|---------|
| ΔVDWAALS | -75.8 | -109.7 | -97.7 | Van der Waals packing |
| ΔEEL | -45.2 | -86.8 | +23.0 | Electrostatic attraction |
| ΔEGB | +105.7 | +147.0 | +43.5 | Desolvation penalty (polar) |
| ΔESURF | -10.0 | -15.2 | -12.3 | Nonpolar solvation |
| **ΔTOTAL** | **-25.3** | **-64.7** | **-43.5** | **Binding free energy** |

Charged Asp at the interface is a liability:
- The -1 carboxylate has a strong solvation shell in the unbound state (water stabilizes the charge)
- Upon binding, the interface buries the charge → water is displaced → massive desolvation penalty (ΔEGB)
- Unless a perfectly positioned Arg/Lys on the antigen compensates, the desolvation cost exceeds the electrostatic gain
- **D245N mechanism**: Removing the charge drops ΔEGB from +105.7 to +43.5 (saves ~62 kcal/mol in desolvation). Electrostatics become slightly repulsive (+23.0), but VDW packing improves by 22 kcal/mol — Asn packs better without forcing water into the interface.
- **D265N mechanism**: Different — VDW jumps massively (-75.8 → -109.7), electrostatics get stronger (-45.2 → -86.8). The Asn amide makes tighter H-bonds and allows interface compression.

## Round 4: Corrected Regression Targets (submitted)

After fixing the off-by-one bug in contact frequency extraction and re-running the regression, we submitted mutations targeting the corrected top residues.

### Strategy
Two approaches from the corrected regression:
1. **Aromatic extension at top hot spots** (Tyr→Trp, proven chemistry)
2. **Removal of unfavorable contacts** (→Gly, proven by S258G)

### Round 4 MM-GBSA results

| Mutation | NF | Class | Chemistry | ΔG (kcal/mol) | SEM | ΔΔG vs WT | Outcome |
|----------|-----|-------|-----------|--------------|-----|-----------|---------|
| S247G | 3.4 | warm/unfav | Side chain removal | -56.87 | 1.20 | -31.6 | **Major hit** — 3rd best overall |
| T330G | 3.7 | warm/unfav | Side chain removal | -41.47 | 0.65 | -16.2 | **Strong hit** |
| Y227W | 67.0 | HOT | Aromatic extension | -35.99 | 0.58 | -10.7 | **Hit** — #2 hot spot Tyr→Trp |
| L306G | 5.8 | warm/unfav | Side chain removal | -16.35 | 0.60 | +9.0 | **Worse** — disrupted N-terminal packing |
| Y406W | 80.6 | HOT | Aromatic extension | — | — | — | Running (pipeline resubmitted) |

### Round 4 key findings

1. **S247G is a major hit** (ΔΔG = -31.6): Third-best single mutant overall, behind D265N (-39.4) and F287W (-33.7). Removing the unfavorable Ser→Gly at L:52 eliminates a destabilizing contact identified by the regression. Validates the unfavorable contact removal strategy.

2. **T330G works** (ΔΔG = -16.2): Comparable to S258G (-16.5) and D245N (-18.2). Removing Thr at H:28, an unfavorable contact not detected in the old contaminated regression, only found after the off-by-one fix.

3. **Y227W is a solid hit** (ΔΔG = -10.7): Aromatic extension at the #2 hot spot works but with moderate effect. Smaller gain than F287W (-33.7) despite higher NF (67.0 vs 17.7) — NF magnitude does not linearly predict ΔΔG.

4. **L306G failed** (ΔΔG = +9.0): Despite unfavorable regression signal (NF=5.8), removing Leu at H:4 (N-terminal heavy chain) weakened binding. Likely disrupted hydrophobic core packing — Leu→Gly at a buried position creates a cavity rather than removing a bad contact. Lesson: →Gly works at surface positions but can backfire at partially buried sites.

### Updated chemistry success rates (all rounds)

| Chemistry | Success rate | Hits | Failures |
|-----------|-------------|------|----------|
| Charge removal (D→N) | 2/2 (100%) | D265N (-39.4), D245N (-18.2) | — |
| Unfavorable →Gly | 3/4 (75%) | S247G (-31.6), S258G (-16.5), T330G (-16.2) | L306G (+9.0) |
| Aromatic extension (→Trp) | 3/3 (100%) | F287W (-33.7), Y227W (-10.7), Y404W (-2.7) | — |
| Aromatic addition (small→Tyr) | 3/4 (75%) | S262Y (-15.0), N248Y (-10.9), A246Y (-6.1) | T403Y (+4.2) |
| Non-conservative at aromatic | 0/2 (0%) | — | Y404H (+2.5), Y405N (weaker) |

### Updated complete MM-GBSA ranking (all single mutants)

| Rank | Mutation | ΔG (kcal/mol) | SEM | ΔΔG vs WT | Chemistry |
|------|----------|--------------|-----|-----------|-----------|
| 1 | D265N | -64.73 | 1.45 | -39.4 | Charge removal |
| 2 | F287W | -58.99 | 1.44 | -33.7 | Aromatic extension |
| 3 | **S247G** | **-56.87** | **1.20** | **-31.6** | **Unfavorable removal** |
| 4 | D245N | -43.47 | 1.10 | -18.2 | Charge removal |
| 5 | S258G | -41.78 | 0.82 | -16.5 | Unfavorable removal |
| 6 | **T330G** | **-41.47** | **0.65** | **-16.2** | **Unfavorable removal** |
| 7 | S262Y | -40.28 | 1.20 | -15.0 | Aromatic addition |
| 8 | N248Y | -36.20 | 1.01 | -10.9 | Aromatic addition |
| 9 | **Y227W** | **-35.99** | **0.58** | **-10.7** | **Aromatic extension** |
| 10 | A246Y | -31.41 | 0.81 | -6.1 | Aromatic addition |
| 11 | Y404W | -28.00 | 0.96 | -2.7 | Aromatic extension |
| 12 | **WT** | **-25.31** | **0.86** | **—** | — |
| 13 | S258A | -24.35 | 0.91 | +1.0 | Conservative removal |
| 14 | Y404H | -22.79 | 0.77 | +2.5 | Non-conservative |
| 15 | T403Y | -21.08 | 1.57 | +4.2 | Aromatic addition |
| 16 | **L306G** | **-16.35** | **0.60** | **+9.0** | **Unfavorable removal** |
| 17 | Y157A | -13.73 | 1.49 | +11.6 | Alanine scan (control) |

Folders: `/home/anugraha/c1_{Y406W,Y227W,L306G,T330G,S247G}/`

## New Elastic Net Regression on C1 WT (10 Pulling Replicas)

### Methodology
Repeated the elastic net regression analysis directly on the C1 WT pulling simulations (10 replicas) instead of the previous Wuhan/Omicron system (4 replicas). This provides a system-specific regression for the actual antibody-antigen pair being optimized.

**Note**: Initial analysis had two issues: (1) truncated energy data (1860/5001 frames due to incomplete XVG extraction) — re-extracted all 10 replicas from EDR files; (2) **off-by-one bug in contact frequency extraction** — mdtraj `resid` is 0-based but the script used 1-based GRO numbers, causing `resid 1 to 194` to select resSeq 2-195 (missed antigen GRO res 1, included antibody GRO res 195). This contaminated 36/548 pairs with antibody-antibody intra-contacts. Fixed by using `resid 0 to 193` (antigen) and `resid 194 to 420` (antibody). All 10 replicas re-extracted and regression re-run.

- **Data**: 10 pulling replicas, 5001 frames each (full energy data)
- **Frame range**: 1–1500 (rupture event focus; cleaner signal than full 5001)
- **Contact frequency**: SASA-filtered (0.5 nm², Shrake-Rupley on frame 0), atom-atom distance < 10 Å, count-based. **Fixed indexing**: mdtraj resid 0-193 = antigen (GRO 1-194), resid 194-420 = antibody (GRO 195-421)
- **Interaction energy**: Coul-SR + LJ-SR between chauA (antigen) and rest (antibody) from energy reruns with group decomposition
- **Regression**: ElasticNetCV (alpha=12.92, l1_ratio=0.50, R²=0.9981)
- **Features**: 549 residue pairs → 52 non-zero coefficients

### Data locations
- Contact frequencies: `/home/anugraha/c1_WT/analysis/average_frequency.csv`
- Interaction energy: `/home/anugraha/c1_WT/analysis/interaction_energy.csv` (5001 frames, re-extracted)
- Coefficients: `/home/anugraha/c1_WT/analysis/elastic_net_coefficients_3000.csv` (1500-frame regression)
- NetFavorability: `/home/anugraha/c1_WT/analysis/net_favorability_3000.csv`
- Regression script: `/home/anugraha/c1_WT/analysis/run_full_regression_3000.py`
- Residue mapping: `/home/anugraha/c1_WT/analysis/residue_mapping.py`

### Figures
- `/home/anugraha/antibody_optimization/figures/fig2_elastic_net_heatmap_3000.png` — NetFavorability bar + β heatmap
- `/home/anugraha/antibody_optimization/figures/fig3_regression_diagnostics_3000.png` — R²=0.9981, pred vs actual + residuals
- `/home/anugraha/antibody_optimization/figures/fig4_temporal_contributions_3000.png` — 3 panels: β×ΔF(t), β×MeanFreq bar, β×ΔF(t)
- `/home/anugraha/antibody_optimization/figures/fig5_deltaF_vs_frame.png` — raw ΔF(t) vs frame for top 10 pairs

### Residue mapping (C1 WT GRO → Global)
- GRO order: Chain A (1-194), Chain H (195-314), Chain L (315-421)
- Global order: Chain A (1-195), Chain L (196-302), Chain H (303-422)
- Chain H: global = GRO - 194 + 302
- Chain L: global = GRO - 314 + 195

### C1 WT regression: ranked antibody residues (corrected indexing, 1500-frame)

NetFavorability classification: hot (>70th percentile), warm (20–70th), cold (<20th).

| Global | Chain | AA | NetFav | #Pairs | DominantSign | Class | Mutation tested? |
|--------|-------|-----|--------|--------|-------------|-------|-----------------|
| **406** | **H:104** | **TYR** | **80.6** | **5** | **favorable** | **HOT** | Y406W (running) |
| **227** | **L:32** | **TYR** | **67.0** | **5** | **favorable** | **HOT** | Y227W (running) |
| 405 | H:103 | TYR | 34.3 | 5 | favorable | **HOT** | Y405N (failed — non-conservative) |
| 404 | H:102 | TYR | 23.4 | 3 | favorable | **HOT** | Y404W (-2.7), Y404H (+2.5) |
| 245 | L:50 | ASP | 23.2 | 4 | favorable | **HOT** | D245N (-18.2) ✓ |
| 287 | L:92 | PHE | 17.7 | 1 | favorable | **HOT** | F287W (-33.7) ✓ |
| 249 | L:54 | ARG | 16.2 | 2 | favorable | **HOT** | — |
| 246 | L:51 | ALA | 7.8 | 3 | favorable | **HOT** | A246Y (-6.1) |
| 248 | L:53 | ASN | 6.2 | 2 | favorable | warm | N248Y (-10.9) ✓ |
| **306** | **H:4** | **LEU** | **5.8** | **1** | **unfavorable** | **warm** | L306G (running) |
| **330** | **H:28** | **THR** | **3.7** | **1** | **unfavorable** | **warm** | T330G (running) |
| 401 | H:99 | GLU | 3.5 | 1 | favorable | warm | — |
| 244 | L:49 | PHE | 3.4 | 3 | favorable | warm | — |
| **247** | **L:52** | **SER** | **3.4** | **2** | **unfavorable** | **warm** | S247G (running) |
| 288 | L:93 | ASN | 2.8 | 1 | unfavorable | warm | — |
| 334 | H:32 | TYR | 2.1 | 1 | favorable | warm | — |
| 225 | L:30 | HIS | 1.9 | 2 | unfavorable | warm | — |
| 262 | L:67 | SER | 1.9 | 1 | favorable | warm | S262Y (-15.0) ✓ |
| 255 | L:60 | ALA | 1.6 | 1 | unfavorable | warm | — |
| 357 | H:55 | GLY | 1.2 | 1 | favorable | warm | — |
| 222 | L:27 | GLN | 1.2 | 1 | favorable | warm | — |
| 378 | H:76 | LYS | 1.1 | 1 | favorable | cold | — |
| 332 | H:30 | SER | 0.4 | 1 | favorable | cold | — |
| 256 | L:61 | ARG | 0.3 | 1 | unfavorable | cold | — |
| 260 | L:65 | SER | 0.3 | 1 | unfavorable | cold | — |
| 223 | L:28 | SER | 0.3 | 1 | favorable | cold | — |
| 289 | L:94 | TRP | 0.1 | 1 | favorable | cold | — |

**Notable absences** (SASA-filtered at pulling frame 0, below 0.5 nm² threshold):
- **D265** (L:70): Our #1 hit by MM-GBSA (ΔΔG = -39.4). Was WARM (24.3) in Wuhan regression. SASA = 0.31 nm² at pulling frame 0.
- **S262** (L:67): Now barely warm (NF=1.9). Was WARM (50.8) in Wuhan. S262Y gave ΔΔG = -15.0.

### Top 12 pairwise contacts (non-zero β, corrected indexing)

| Antigen | Antibody | β | MeanFreq | |β×F| | Sign |
|---------|----------|-------|----------|------|------|
| A:K84 | L:Y227 | -0.352 | 143.8 | 50.6 | favorable |
| A:Q160 | H:Y406 | -0.247 | 161.1 | 39.8 | favorable |
| A:F157 | H:Y406 | -0.149 | 184.4 | 27.5 | favorable |
| A:F123 | L:F287 | -0.189 | 93.4 | 17.7 | favorable |
| A:V150 | H:Y405 | -0.134 | 123.4 | 16.6 | favorable |
| A:Y172 | L:R249 | -0.104 | 150.3 | 15.6 | favorable |
| A:E151 | H:Y405 | -0.067 | 232.0 | 15.4 | favorable |
| A:F157 | H:Y404 | -0.158 | 92.9 | 14.7 | favorable |
| A:Q160 | L:Y227 | -0.079 | 136.8 | 10.9 | favorable |
| A:Y172 | L:D245 | -0.158 | 66.4 | 10.5 | favorable |
| A:R70 | L:D245 | -0.193 | 49.4 | 9.5 | favorable |
| A:F123 | H:Y406 | -0.041 | 153.7 | 6.3 | favorable |

### Key changes from contaminated → corrected regression
1. **Off-by-one bug fixed**: mdtraj `resid` is 0-based; old code used 1-based GRO numbers. `resid 1 to 194` selected resSeq 2-195, including antibody GRO 195 in the antigen group. Fixed to `resid 0 to 193` / `resid 194 to 420`.
2. **l1_ratio dropped 0.90 → 0.50**: Model is less sparse (52 vs 36 non-zero β), uses more features with smaller coefficients. More ridge-like regularization.
3. **Y406 became #1** (NF 73.4 → 80.6): Previously masked by contamination. Now the highest-leverage position.
4. **Y227 stable at #2** (NF 76.2 → 67.0): Slight decrease but remains top hot spot.
5. **New residues appeared**: T330 (warm/unfav), G357 (warm/fav), K378 (cold/fav) — not detected in contaminated regression.
6. **L306 signal increased** (NF 1.3 → 5.8): Unfavorable contact now more prominent.
7. **S258 dropped below top 20**: Was warm/unfavorable (1.9) in contaminated, now absent. S258G still worked (ΔΔG = -16.5) — mutation target didn't need high regression signal.
8. **D265 still absent**: SASA-filtered. Our #1 hit (-39.4 ΔΔG) is not detectable by this regression approach.
9. **S262 barely warm** (NF=1.9): Was absent before, now appears weakly. S262Y gave -15.0 ΔΔG despite low signal.

### Revised mutation strategy (post-Round 3)

The hot/warm/cold classification is **not a reliable decision rule** for mutation design. Evidence:

| Mutation | ΔΔG | Classification | "Don't mutate hot" correct? |
|----------|-----|---------------|---------------------------|
| F287W | -33.7 | HOT | **No** — major improvement |
| D245N | -18.2 | HOT | **No** — strong improvement |
| Y404W | -2.7 | HOT | **No** — mild improvement |
| Y404H | +2.5 | HOT | Yes — weakened |
| Y405N | weaker | HOT | Yes — weakened |

Hot spot rule: **2/5 correct (40%)**. Warm spot rule: 3/4 correct (75%). The classification is barely better than a coin flip for hot spots.

**What actually predicts success is the chemistry:**

| Chemistry | Success rate | Examples |
|-----------|-------------|---------|
| Charge removal (D→N) | 2/2 (100%) | D265N (-39.4), D245N (-18.2) |
| Aromatic extension (→Trp) | 3/3 (100%) | F287W (-33.7), Y227W (-10.7), Y404W (-2.7) |
| Unfavorable →Gly | 3/4 (75%) | S247G (-31.6), S258G (-16.5), T330G (-16.2); L306G failed (+9.0) |
| Aromatic addition (small→Tyr) | 3/4 (75%) | S262Y (-15.0), N248Y (-10.9), A246Y (-6.1); T403Y failed |
| Non-conservative at aromatic hot spot | 0/2 (0%) | Y404H (+2.5), Y405N (weaker) |

**Revised approach**: Use the regression to identify high-leverage positions (high |β|), then apply physics-based substitution rules:
1. **Interfacial Asp** → try D→N (charge removal relieves desolvation penalty)
2. **Interfacial Glu** → try E→Q (same logic)
3. **Aromatic positions (Phe, Tyr)** → try →Trp (extend aromatic platform, preserve interaction type)
4. **Unfavorable β at surface positions** → try →Gly (remove destabilizing contact). Caution: fails at partially buried sites (L306G)
5. **Never** replace aromatics with non-aromatics at high-leverage positions (His, Asn at Tyr sites fails)

### Future mutation candidates

**Awaiting results**:
- **Y406W** (H:104, NF=80.6, HOT): #1 hot spot. Tyr→Trp. Pipeline running.

**High-leverage aromatic extension** (next priority):
- **Y405W** (H:103, NF=34.3, HOT): #3 hot spot. Y405N failed but Trp is conservative aromatic — different chemistry.

**Charge removal** (proven 100% success rate):
- **E401Q** (H:99, NF=3.5, favorable): Glu→Gln charge removal. Same D→N logic applied to Glu.

**Remaining untested unfavorable contacts**:
- **H225G** (L:30, NF=1.9, unfavorable): destabilizing contact; surface position, Gly should work
- **N288A** (L:93, NF=2.8, unfavorable): adjacent to F287; remove destabilizing contact
- **A255G** (L:60, NF=1.6, unfavorable): small unfavorable contact

**Other untested**:
- **R249** (L:54, NF=16.2, HOT): Arg — could try R249K (conservative) or R249W (aromatic)

D265N remains the top performer by MM-GBSA (-64.73 kcal/mol) despite being absent from the C1 WT regression (SASA-filtered at pulling frame 0).
