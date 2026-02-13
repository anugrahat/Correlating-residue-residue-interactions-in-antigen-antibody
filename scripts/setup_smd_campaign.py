#!/usr/bin/env python3
"""Create reproducible C1 mutant SMD folders and job scripts.

This script does not submit jobs. It creates per-mutation folders at
/home/anugraha/c1_<MUT> and writes:
  - mutated input PDB
  - shared MDP templates copied from your existing workflow
  - a Slurm script that builds, equilibrates, and runs pulling
"""

from __future__ import annotations

import argparse
import random
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


AA1_TO_3: Dict[str, str] = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}

AA3_TO_1: Dict[str, str] = {v: k for k, v in AA1_TO_3.items()}

TEMPLATE_PULL_DIR = Path(
    "/home/anugraha/Wuhan_103_PHE/charmm-gui-5253081658/gromacs/pulling"
)
DEFAULT_BASE_PDB = Path("/home/anugraha/antibody_optimization/c1_87g7.pdb")
DEFAULT_OUT_ROOT = Path("/home/anugraha")
DEFAULT_MUTATIONS = ["Y405N", "Y405Q", "Y406W", "Y404W", "F287W", "D265N", "S262Y"]
DEFAULT_FF_DIR = Path("/home/anugraha/amber14sb.ff")
DEFAULT_GMXLIB = "/home/anugraha:/software/gromacs/2024/ucdhpc-20.04/share/gromacs/top"

BACKBONE_ALLOWED = {
    "N",
    "CA",
    "C",
    "O",
    "OXT",
    "H",
    "H1",
    "H2",
    "H3",
    "HN",
    "HT1",
    "HT2",
    "HT3",
    "HA",
    "HA1",
    "HA2",
    "HA3",
    "CB",
    "HB",
    "HB1",
    "HB2",
    "HB3",
}

COPY_FILES = [
    "ions.mdp",
    "step4.0_minimization.mdp",
    "step4.1_equilibration.mdp",
    "umb2.mdp",
]


@dataclass(frozen=True)
class Mutation:
    wt: str
    global_resid: int
    mut: str

    @property
    def tag(self) -> str:
        return f"{self.wt}{self.global_resid}{self.mut}"


def parse_mutation(spec: str) -> Mutation:
    m = re.fullmatch(r"([A-Z])(\d+)([A-Z])", spec.strip().upper())
    if not m:
        raise ValueError(f"Invalid mutation '{spec}'. Expected form like Y405N.")
    wt, resid, mut = m.groups()
    if wt not in AA1_TO_3 or mut not in AA1_TO_3:
        raise ValueError(f"Unsupported amino acid in '{spec}'.")
    return Mutation(wt=wt, global_resid=int(resid), mut=mut)


def parse_mutation_spec(spec: str) -> List[Mutation] | None:
    """Parse a mutation specification that may contain multiple mutations.

    Supports:
      - "WT"           -> None (wild-type, no mutations)
      - "Y405N"        -> [Mutation(Y,405,N)]
      - "D265N+F287W"  -> [Mutation(D,265,N), Mutation(F,287,W)]
    """
    spec = spec.strip().upper()
    if spec == "WT":
        return None
    parts = spec.split("+")
    return [parse_mutation(p) for p in parts]


def mutation_spec_tag(mutations: List[Mutation] | None) -> str:
    """Return a folder-safe tag for a mutation spec: 'WT' or 'D265N_F287W'."""
    if mutations is None:
        return "WT"
    return "_".join(m.tag for m in mutations)


def global_to_chain_resid(global_resid: int) -> Tuple[str, int]:
    """Map regression/global index to PDB chain+resid.

    Mapping used in your workflow:
      1..195   = chain A residues 1..195  (antigen / RBD)
      196..302 = chain L residues 1..107
      303..422 = chain H residues 1..120
    """
    if 1 <= global_resid <= 195:
        return "A", global_resid
    if 196 <= global_resid <= 302:
        return "L", global_resid - 195
    if 303 <= global_resid <= 422:
        return "H", global_resid - 302
    raise ValueError(
        f"Global residue {global_resid} is outside range 1..422."
    )


def _residue_key(line: str) -> Tuple[str, int, str]:
    chain = line[21]
    resid = int(line[22:26])
    icode = line[26]
    return chain, resid, icode


def mutate_pdb_single(
    lines: List[str],
    mutation: Mutation,
    force: bool = False,
) -> List[str]:
    """Apply a single mutation to PDB lines in memory. Returns new lines."""
    chain, local_resid = global_to_chain_resid(mutation.global_resid)
    target_key = (chain, local_resid, " ")
    mut_res3 = AA1_TO_3[mutation.mut]

    found_resname = None
    out_lines: List[str] = []

    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            out_lines.append(line)
            continue

        key = _residue_key(line)
        if key[0] == target_key[0] and key[1] == target_key[1]:
            atom_name = line[12:16].strip()
            if mut_res3 == "GLY" and atom_name in {"CB", "HB", "HB1", "HB2", "HB3"}:
                continue
            if atom_name not in BACKBONE_ALLOWED:
                continue

            current_res3 = line[17:20].strip()
            if found_resname is None:
                found_resname = current_res3
            new_line = f"{line[:17]}{mut_res3:>3}{line[20:]}"
            out_lines.append(new_line)
        else:
            out_lines.append(line)

    if found_resname is None:
        raise ValueError(
            f"Did not find residue chain {chain} resid {local_resid}."
        )

    found_wt = AA3_TO_1.get(found_resname, "?")
    if found_wt != mutation.wt and not force:
        raise ValueError(
            f"Mutation {mutation.tag} expected WT {mutation.wt} but found "
            f"{found_resname} ({found_wt}) at chain {chain} resid {local_resid}. "
            "Use --force to override."
        )

    return out_lines


def mutate_pdb(
    pdb_in: Path,
    pdb_out: Path,
    mutations: List[Mutation],
    force: bool = False,
) -> List[Tuple[str, int, str, str]]:
    """Apply one or more mutations to a PDB file.

    Returns a list of (chain, local_resid, found_res3, mut_res3) per mutation.
    """
    lines = pdb_in.read_text().splitlines(keepends=True)
    results = []

    for mutation in mutations:
        chain, local_resid = global_to_chain_resid(mutation.global_resid)
        mut_res3 = AA1_TO_3[mutation.mut]

        # Find the WT residue before mutating
        found_resname = None
        for line in lines:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            key = _residue_key(line)
            if key[0] == chain and key[1] == local_resid:
                found_resname = line[17:20].strip()
                break

        lines = mutate_pdb_single(lines, mutation, force=force)
        results.append((chain, local_resid, found_resname or "???", mut_res3))

    pdb_out.write_text("".join(lines))
    return results


def write_sbatch_script(
    folder: Path,
    antigen_res_count: int,
    auto_align_script: Path | None = None,
    ff_dir: Path = DEFAULT_FF_DIR,
    gmxlib: str = DEFAULT_GMXLIB,
) -> None:
    if auto_align_script is None:
        auto_align_script = Path(__file__).resolve().parent / "auto_align_for_pull.py"
    script = f"""#!/bin/bash
#SBATCH --job-name={folder.name}_pull
#SBATCH --account=ahnlab
#SBATCH --partition=gpu-ahn
#SBATCH --chdir={folder}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:2
#SBATCH --time=120:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=anuthyagatur@ucdavis.edu

set -euo pipefail

module load amber/24
module load gromacs/2024

if command -v gmx_mpi >/dev/null 2>&1; then
  GMX=gmx_mpi
else
  GMX=gmx
fi

# Serial wrapper for non-mdrun GROMACS commands (1 rank)
RUN1="srun --ntasks=1"
# Parallel mdrun (2 ranks, 2 GPUs, 1 dedicated PME rank)
MDRUN="srun $GMX mdrun -npme 1"

GPU_FLAGS="-nb gpu -bonded gpu -pme gpu"
EM_FLAGS="-nb gpu"
PULL_GPU_FLAGS="-nb gpu -bonded gpu -pme gpu"  # no -update gpu (Nose-Hoover incompatible)
export GMXLIB="{gmxlib}"

echo "[INFO] Using GROMACS command: $GMX"
echo "[INFO] Working dir: $(pwd)"

mkdir -p work pull logs

# ============================================================
# PHASE 1: Normal MD in small box
# ============================================================

# 1) Build complete structure with tleap (rebuilds missing sidechain atoms)
if [ -f work/mutant_complete.pdb ]; then
  echo "[SKIP] Step 1: tleap already done"
else
  echo "[RUN]  Step 1: tleap"
  cat > work/tleap.in <<'EOF'
source leaprc.protein.ff14SB
mol = loadpdb mutant_input.pdb
check mol
savepdb mol work/mutant_complete.pdb
quit
EOF
  tleap -f work/tleap.in > logs/tleap.log 2>&1
fi

# 2) Generate GROMACS topology (separate chains, no -merge all)
if [ -f work/mutant_GMX.gro ] && [ -f work/mutant_GMX_clean.top ]; then
  echo "[SKIP] Step 2: pdb2gmx already done"
else
  echo "[RUN]  Step 2: pdb2gmx"
  printf 'y\\ny\\ny\\ny\\ny\\ny\\ny\\ny\\ny\\ny\\n' | $RUN1 $GMX pdb2gmx -f work/mutant_complete.pdb -o work/mutant_GMX.gro \\
      -p work/mutant_GMX.top -i posre.itp \\
      -ff amber14sb -water tip3p -ignh -ss > logs/pdb2gmx.log 2>&1
  # Move posre files into work/ alongside the chain ITPs
  mv posre_Protein*.itp work/
  # The chain ITPs include "posre_Protein*.itp" (bare name, same dir) -- correct
  # Save clean protein-only topology before solvation modifies it
  cp work/mutant_GMX.top work/mutant_GMX_clean.top
fi

# 3) Small box + solvate + ions (normal size for equilibration MD)
if [ -f work/ions.gro ]; then
  echo "[SKIP] Step 3: solvation + ions already done"
else
  echo "[RUN]  Step 3: solvation + ions"
  $RUN1 $GMX editconf -f work/mutant_GMX.gro -o work/newbox.gro -c -d 1.2 -bt triclinic
  $RUN1 $GMX solvate -cp work/newbox.gro -cs spc216.gro -o work/solv.gro -p work/mutant_GMX.top
  $RUN1 $GMX grompp -f inputs/ions.mdp -c work/solv.gro -p work/mutant_GMX.top -o work/ions.tpr -maxwarn 10
  echo SOL | $RUN1 $GMX genion -s work/ions.tpr -o work/ions.gro -p work/mutant_GMX.top -pname Na -nname Cl -neutral -conc 0.15
fi

# 4a) Energy minimization
if [ -f work/em.gro ]; then
  echo "[SKIP] Step 4a: EM already done"
else
  echo "[RUN]  Step 4a: EM"
  $RUN1 $GMX grompp -f inputs/step4.0_minimization.mdp -c work/ions.gro -p work/mutant_GMX.top -o work/em.tpr -maxwarn 10
  $MDRUN -deffnm work/em $EM_FLAGS
fi

# 4b) NVT equilibration
if [ -f work/nvt.gro ]; then
  echo "[SKIP] Step 4b: NVT already done"
else
  echo "[RUN]  Step 4b: NVT"
  $RUN1 $GMX grompp -f inputs/step4.1_equilibration.mdp -c work/em.gro -r work/em.gro -p work/mutant_GMX.top -o work/nvt.tpr -maxwarn 10
  $MDRUN -deffnm work/nvt $GPU_FLAGS
fi

# 4c) NPT equilibration
if [ -f work/npt.gro ]; then
  echo "[SKIP] Step 4c: NPT already done"
else
  echo "[RUN]  Step 4c: NPT"
  $RUN1 $GMX grompp -f inputs/npt1.mdp -c work/nvt.gro -r work/nvt.gro -t work/nvt.cpt -p work/mutant_GMX.top -o work/npt.tpr -maxwarn 10
  $MDRUN -deffnm work/npt $GPU_FLAGS
fi

# 5) 10 ns production MD (2 fs timestep -> 5,000,000 steps)
if [ -f work/md10.gro ]; then
  echo "[SKIP] Step 5: 10ns MD already done"
else
  echo "[RUN]  Step 5: 10ns MD"
  $RUN1 $GMX grompp -f inputs/md10.mdp -c work/npt.gro -r work/npt.gro -t work/npt.cpt -p work/mutant_GMX.top -o work/md10.tpr -maxwarn 10
  $MDRUN -deffnm work/md10 $GPU_FLAGS
fi

# 6) Cluster protein chains into same periodic image and extract protein
if [ -f pull/md10_protein.gro ]; then
  echo "[SKIP] Step 6: PBC cluster + protein extraction already done"
else
  echo "[RUN]  Step 6: PBC cluster + protein extraction"
  printf 'Protein\\nProtein\\nSystem\\n' | $RUN1 $GMX trjconv -f work/md10.gro -s work/md10.tpr -o work/md10_clustered.gro -pbc cluster -center
  echo Protein | $RUN1 $GMX trjconv -f work/md10_clustered.gro -s work/md10.tpr -o pull/md10_protein.gro
fi

echo "[INFO] Phase 1 complete: 10 ns MD done in small box"

# ============================================================
# PHASE 2: Pulling in PBC-safe big box
# ============================================================

# 7) Auto-align interface for pulling + build PBC-safe asymmetric box
if [ -f pull/aligned.gro ]; then
  echo "[SKIP] Step 7: auto-align already done"
else
  echo "[RUN]  Step 7: auto-align"
  python3 {auto_align_script} \\
      --gro pull/md10_protein.gro --out pull/aligned.gro \\
      --antigen-res-count {antigen_res_count} --pull-dim Y --mdp inputs/umb2.mdp \\
      --margin-pull-nm 2.0 --margin-perp-nm 1.5 \\
      --report-png pull/alignment_report.png --fallback-to-com \\
      --orient-smaller-positive
fi

# 8) Re-solvate in big box with clean topology
if [ -f pull/ions.gro ]; then
  echo "[SKIP] Step 8: re-solvation already done"
else
  echo "[RUN]  Step 8: re-solvation"
  cp work/mutant_GMX_clean.top pull/pull.top
  # Copy chain ITP and posre files for pull topology
  cp work/mutant_GMX_Protein*.itp pull/
  cp work/posre_Protein*.itp pull/
  $RUN1 $GMX solvate -cp pull/aligned.gro -cs spc216.gro -o pull/solv.gro -p pull/pull.top
  $RUN1 $GMX grompp -f inputs/ions.mdp -c pull/solv.gro -p pull/pull.top -o pull/ions.tpr -maxwarn 10
  echo SOL | $RUN1 $GMX genion -s pull/ions.tpr -o pull/ions.gro -p pull/pull.top -pname Na -nname Cl -neutral -conc 0.15
fi

# 9a) Pull-box EM
if [ -f pull/em.gro ]; then
  echo "[SKIP] Step 9a: pull EM already done"
else
  echo "[RUN]  Step 9a: pull EM"
  $RUN1 $GMX grompp -f inputs/step4.0_minimization.mdp -c pull/ions.gro -p pull/pull.top -o pull/em.tpr -maxwarn 10
  $MDRUN -deffnm pull/em $EM_FLAGS
fi

# 9b) Pull-box NPT re-equilibration
if [ -f pull/npt.gro ]; then
  echo "[SKIP] Step 9b: pull NPT already done"
else
  echo "[RUN]  Step 9b: pull NPT"
  $RUN1 $GMX grompp -f inputs/pull_npt.mdp -c pull/em.gro -r pull/em.gro -p pull/pull.top -o pull/npt.tpr -maxwarn 10
  $MDRUN -deffnm pull/npt $GPU_FLAGS
fi

# 10) Build pull index groups compatible with umb2.mdp names: chauA/rest
if [ -f pull/rbd_2_reg.ndx ]; then
  echo "[SKIP] Step 10: pull index already built"
else
  echo "[RUN]  Step 10: build pull index"
  python3 build_pull_index.py pull/npt.gro pull/rbd_2_reg.ndx {antigen_res_count}
fi

# 11) Pulling run from re-equilibrated structure
if [ -f pull/replica2.gro ]; then
  echo "[SKIP] Step 11: pull mdrun already done"
else
  echo "[RUN]  Step 11: pull mdrun"
  $RUN1 $GMX grompp -f inputs/umb2.mdp -r pull/npt.gro -c pull/npt.gro -n pull/rbd_2_reg.ndx -t pull/npt.cpt -p pull/pull.top -o pull/replica2.tpr -maxwarn 20
  $MDRUN -v -s pull/replica2.tpr -deffnm pull/replica2 -pf pull/replica2_pullf.xvg -px pull/replica2_pullx.xvg $PULL_GPU_FLAGS
fi

echo "[DONE] Pulling outputs:"
echo "  pull/replica2_pullf.xvg"
echo "  pull/replica2_pullx.xvg"
echo "  pull/alignment_report.png"
"""
    out = folder / "run_pipeline.sbatch"
    out.write_text(script)
    out.chmod(0o755)


def write_index_builder(folder: Path) -> None:
    code = """#!/usr/bin/env python3
import sys

PROT = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","HID","HIE","HIP",
    "ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL","CYX"
}

if len(sys.argv) != 4:
    print("Usage: build_pull_index.py <in.gro> <out.ndx> <antigen_res_count>", file=sys.stderr)
    sys.exit(1)

gro_path, ndx_path, antigen_res_count = sys.argv[1], sys.argv[2], int(sys.argv[3])

atoms = []
with open(gro_path) as f:
    f.readline()
    n_atoms = int(f.readline().strip())
    for idx in range(1, n_atoms + 1):
        line = f.readline()
        resid = int(line[0:5])
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        atoms.append((idx, resid, resname, atomname))

# Track residues in appearance order for protein atoms only.
residue_order = []
residue_to_atoms = {}
protein_atom_set = set()
last_key = None
for idx, resid, resname, atomname in atoms:
    if resname not in PROT:
        continue
    protein_atom_set.add(idx)
    key = (resid, resname)
    if key != last_key:
        residue_order.append(key)
        residue_to_atoms[key] = []
        last_key = key
    residue_to_atoms[key].append(idx)

if len(residue_order) < antigen_res_count:
    raise ValueError(
        f"Protein residue count ({len(residue_order)}) is smaller than antigen count ({antigen_res_count})."
    )

chau_res = residue_order[:antigen_res_count]
rest_res = residue_order[antigen_res_count:]
chau_atoms = [a for r in chau_res for a in residue_to_atoms[r]]
rest_atoms = [a for r in rest_res for a in residue_to_atoms[r]]
all_atoms = [a[0] for a in atoms]
protein_atoms = sorted(protein_atom_set)
non_protein_atoms = [a for a in all_atoms if a not in protein_atom_set]
water_atoms = [idx for idx, _resid, resname, _atom in atoms if resname in {"SOL", "WAT"}]

def write_group(fh, name, arr, chunk=15):
    fh.write(f"[ {name} ]\\n")
    for i in range(0, len(arr), chunk):
        fh.write(" ".join(str(x) for x in arr[i:i+chunk]) + "\\n")
    fh.write("\\n")

with open(ndx_path, "w") as out:
    write_group(out, "System", all_atoms)
    write_group(out, "Protein", protein_atoms)
    write_group(out, "Non-Protein", non_protein_atoms)
    write_group(out, "Water", water_atoms)
    write_group(out, "chauA", chau_atoms)
    write_group(out, "rest", rest_atoms)

print(
    f"Wrote {ndx_path}: Protein={len(protein_atoms)} "
    f"Non-Protein={len(non_protein_atoms)} chauA={len(chau_atoms)} rest={len(rest_atoms)}"
)
"""
    out = folder / "build_pull_index.py"
    out.write_text(code)
    out.chmod(0o755)


def count_antigen_residues(base_pdb: Path, antigen_chain: str = "A") -> int:
    seen = set()
    for line in base_pdb.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[21] != antigen_chain:
            continue
        resid = int(line[22:26])
        icode = line[26]
        seen.add((resid, icode))
    if not seen:
        raise ValueError(f"No residues found for antigen chain '{antigen_chain}' in {base_pdb}")
    return len(seen)


def write_npt_mdp(folder: Path) -> None:
    npt_text = """title       = NPT Equilibration
define      = -DPOSRES
integrator  = md
nsteps      = 50000
dt          = 0.002
nstxout     = 1000
nstvout     = 1000
nstenergy   = 1000
nstlog      = 1000
continuation         = yes
constraint_algorithm = lincs
constraints          = h-bonds
lincs_iter           = 1
lincs_order          = 4
ns_type     = grid
nstlist     = 5
rlist       = 1.4
rcoulomb    = 1.4
rvdw        = 1.4
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.16
tcoupl      = Berendsen
tc-grps     = Protein   Non-Protein
tau_t       = 0.1       0.1
ref_t       = 310       310
pcoupl              = Berendsen
pcoupltype          = isotropic
tau_p               = 2.0
ref_p               = 1.0
compressibility     = 4.5e-5
refcoord_scaling    = com
pbc     = xyz
DispCorr    = EnerPres
nstcomm         = 10
comm-mode       = Linear
comm-grps       = System
"""
    (folder / "npt1.mdp").write_text(npt_text)


def write_pull_npt_mdp(folder: Path) -> None:
    """Short 100ps NPT to relax water shell after re-solvation for pulling."""
    pull_npt_text = """title       = Pull-box NPT Re-equilibration
define      = -DPOSRES
integrator  = md
nsteps      = 50000
dt          = 0.002
nstxout     = 1000
nstvout     = 1000
nstenergy   = 1000
nstlog      = 1000
continuation         = no
constraint_algorithm = lincs
constraints          = h-bonds
lincs_iter           = 1
lincs_order          = 4
ns_type     = grid
nstlist     = 5
rlist       = 1.4
rcoulomb    = 1.4
rvdw        = 1.4
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.16
tcoupl      = Berendsen
tc-grps     = Protein   Non-Protein
tau_t       = 0.1       0.1
ref_t       = 310       310
pcoupl              = Berendsen
pcoupltype          = isotropic
tau_p               = 2.0
ref_p               = 1.0
compressibility     = 4.5e-5
refcoord_scaling    = com
pbc     = xyz
DispCorr    = EnerPres
gen_vel     = yes
gen_temp    = 310
gen_seed    = -1
nstcomm         = 10
comm-mode       = Linear
comm-grps       = System
"""
    (folder / "pull_npt.mdp").write_text(pull_npt_text)


def write_md10_mdp(folder: Path) -> None:
    md_text = """title       = 10 ns Production MD
integrator  = md
nsteps      = 5000000
dt          = 0.002
nstxout     = 5000
nstvout     = 5000
nstenergy   = 1000
nstlog      = 1000
continuation         = yes
constraint_algorithm = lincs
constraints          = h-bonds
lincs_iter           = 1
lincs_order          = 4
nstlist     = 20
rlist       = 1.2
rcoulomb    = 1.2
rvdw        = 1.2
coulombtype     = PME
pme_order       = 4
fourierspacing  = 0.16
tcoupl      = V-rescale
tc-grps     = Protein   Non-Protein
tau_t       = 1.0       1.0
ref_t       = 310       310
pcoupl              = Parrinello-Rahman
pcoupltype          = isotropic
tau_p               = 5.0
ref_p               = 1.0
compressibility     = 4.5e-5
pbc     = xyz
DispCorr    = EnerPres
nstcomm         = 100
comm-mode       = Linear
comm-grps       = System
"""
    (folder / "md10.mdp").write_text(md_text)


def write_replicate_mdp(folder: Path, rep: int, seed: int) -> None:
    """Write a pull MDP with gen_vel=yes and a unique seed for this replicate."""
    rep_dir = folder / f"pull_rep{rep}"
    rep_dir.mkdir(parents=True, exist_ok=True)
    mdp_text = f"""title       = Umbrella pulling simulation - replicate {rep}
define      = -DPOSRES_B
integrator  = md
dt          = 0.002
tinit       = 0
nsteps      = 2500000
nstcomm     = 10
nstxout     = 5000
nstvout     = 5000
nstfout     = 500
nstxtcout   = 500
nstenergy   = 500
constraint_algorithm    = lincs
constraints             = h-bonds
continuation            = no
cutoff-scheme   = Verlet
nstlist         = 20
ns_type         = grid
rlist           = 1.4
rcoulomb        = 1.4
rvdw            = 1.4
coulombtype     = PME
fourierspacing  = 0.12
fourier_nx      = 0
fourier_ny      = 0
fourier_nz      = 0
pme_order       = 4
ewald_rtol      = 1e-5
optimize_fft    = yes
Tcoupl      = Nose-Hoover
tc_grps     = Protein   Non-Protein
tau_t       = 1.0       1.0
ref_t       = 310       310
Pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau_p           = 2.0
compressibility = 4.5e-5
ref_p           = 1.0
refcoord_scaling = com
gen_vel     = yes
gen_temp    = 310
gen_seed    = {seed}
pbc     = xyz
DispCorr    = EnerPres
pull                    = yes
pull_ncoords            = 1
pull_ngroups            = 2
pull_group1_name        = chauA
pull_group2_name        = rest
pull_coord1_type        = umbrella
pull_coord1_geometry    = distance
pull_coord1_dim         = N Y N
pull_coord1_groups      = 1 2
pull_coord1_start       = yes
pull_coord1_rate        = 0.001
pull_coord1_k           = 1350
"""
    (rep_dir / f"umb2_rep{rep}.mdp").write_text(mdp_text)


def write_replicate_sbatch(
    folder: Path,
    rep: int,
    gmxlib: str = DEFAULT_GMXLIB,
) -> None:
    """Write a standalone sbatch script for a pull replicate (1 GPU)."""
    script = f"""#!/bin/bash
#SBATCH --job-name={folder.name}_rep{rep}
#SBATCH --account=ahnlab
#SBATCH --partition=gpu-ahn
#SBATCH --chdir={folder}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --time=120:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=anuthyagatur@ucdavis.edu

set -euo pipefail
module load gromacs/2024
if command -v gmx_mpi >/dev/null 2>&1; then GMX=gmx_mpi; else GMX=gmx; fi
RUN1="srun --ntasks=1"
PULL_GPU_FLAGS="-nb gpu -bonded gpu -pme gpu"
export GMXLIB="{gmxlib}"

echo "[INFO] Replicate {rep} pull run for {folder.name}"

if [ -f pull_rep{rep}/replica2.gro ]; then
  echo "[SKIP] Pull mdrun already done"
else
  echo "[RUN]  grompp + mdrun for replicate {rep}"
  $RUN1 $GMX grompp -f pull_rep{rep}/umb2_rep{rep}.mdp -r pull/npt.gro -c pull/npt.gro \\
      -n pull/rbd_2_reg.ndx -p pull/pull.top -o pull_rep{rep}/replica2.tpr -maxwarn 20
  srun $GMX mdrun -v -s pull_rep{rep}/replica2.tpr -deffnm pull_rep{rep}/replica2 \\
      -pf pull_rep{rep}/replica2_pullf.xvg -px pull_rep{rep}/replica2_pullx.xvg $PULL_GPU_FLAGS
fi

echo "[DONE] Replicate {rep} outputs:"
echo "  pull_rep{rep}/replica2_pullf.xvg"
echo "  pull_rep{rep}/replica2_pullx.xvg"
"""
    out = folder / f"run_pull_rep{rep}.sbatch"
    out.write_text(script)
    out.chmod(0o755)


def setup_one(
    out_root: Path,
    base_pdb: Path,
    mutations: List[Mutation] | None,
    template_dir: Path,
    antigen_res_count: int,
    force: bool,
    auto_align_script: Path | None = None,
    replicates: int = 1,
) -> Path:
    tag = mutation_spec_tag(mutations)
    folder = out_root / f"c1_{tag}"
    folder.mkdir(parents=True, exist_ok=True)

    inputs = folder / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    for name in COPY_FILES:
        shutil.copy2(template_dir / name, inputs / name)
    write_npt_mdp(inputs)
    write_pull_npt_mdp(inputs)
    write_md10_mdp(inputs)

    mut_pdb = folder / "mutant_input.pdb"
    if mutations is None:
        # Wild-type: just copy the base PDB unchanged
        shutil.copy2(base_pdb, mut_pdb)
        report = (
            "mutation=WT (no mutations)\n"
            f"antigen_residue_count={antigen_res_count}\n"
            "pull_groups=chauA(restraint group1, antigen),rest(group2, antibody)\n"
        )
    else:
        results = mutate_pdb(
            pdb_in=base_pdb, pdb_out=mut_pdb, mutations=mutations, force=force
        )
        report_lines = []
        for i, (mutation, (chain, local_resid, found_res3, mut_res3)) in enumerate(
            zip(mutations, results)
        ):
            prefix = f"mutation{i+1}" if len(mutations) > 1 else "mutation"
            report_lines.append(f"{prefix}={mutation.tag}")
            report_lines.append(f"{prefix}_global_resid={mutation.global_resid}")
            report_lines.append(f"{prefix}_mapped_chain={chain}")
            report_lines.append(f"{prefix}_mapped_chain_resid={local_resid}")
            report_lines.append(f"{prefix}_found_res3={found_res3}")
            report_lines.append(f"{prefix}_target_res3={mut_res3}")
        report_lines.append(f"antigen_residue_count={antigen_res_count}")
        report_lines.append(
            "pull_groups=chauA(restraint group1, antigen),rest(group2, antibody)"
        )
        report = "\n".join(report_lines) + "\n"

    (folder / "mutation_mapping.txt").write_text(report)
    write_index_builder(folder)
    write_sbatch_script(folder, antigen_res_count, auto_align_script)

    # Generate replicate pull scripts (rep2, rep3, ...) with different seeds
    for rep in range(2, replicates + 1):
        seed = random.randint(10000, 99999)
        write_replicate_mdp(folder, rep, seed)
        write_replicate_sbatch(folder, rep)

    return folder


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--mutations",
        nargs="+",
        default=DEFAULT_MUTATIONS,
        help="Mutation specs: single (Y405N), multi (D265N+F287W), or WT. "
             "Use + to combine mutations on the same system.",
    )
    p.add_argument("--base-pdb", type=Path, default=DEFAULT_BASE_PDB)
    p.add_argument("--out-root", type=Path, default=DEFAULT_OUT_ROOT)
    p.add_argument("--template-dir", type=Path, default=TEMPLATE_PULL_DIR)
    p.add_argument(
        "--auto-align-script",
        type=Path,
        default=Path(__file__).resolve().parent / "auto_align_for_pull.py",
        help="Path to auto_align_for_pull.py (default: alongside this script)",
    )
    p.add_argument("--antigen-chain", default="A")
    p.add_argument(
        "--replicates",
        type=int,
        default=1,
        help="Number of pull replicates (default 1). Rep1 is the main pipeline. "
             "Rep2+ get separate sbatch scripts with different initial velocities.",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="Allow WT mismatch in input PDB without aborting.",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if not args.base_pdb.exists():
        raise FileNotFoundError(f"Missing base PDB: {args.base_pdb}")
    if not args.template_dir.exists():
        raise FileNotFoundError(f"Missing template directory: {args.template_dir}")

    antigen_res_count = count_antigen_residues(args.base_pdb, args.antigen_chain)
    mutation_specs = [parse_mutation_spec(m) for m in args.mutations]

    auto_align = args.auto_align_script.resolve()
    if not auto_align.exists():
        raise FileNotFoundError(f"Missing auto-align script: {auto_align}")

    created = []
    for mutations in mutation_specs:
        folder = setup_one(
            out_root=args.out_root,
            base_pdb=args.base_pdb,
            mutations=mutations,
            template_dir=args.template_dir,
            antigen_res_count=antigen_res_count,
            force=args.force,
            auto_align_script=auto_align,
            replicates=args.replicates,
        )
        created.append(folder)

    print("Created mutation folders:")
    for folder in created:
        print(f"  {folder}")
        print(f"    rep1: sbatch {folder / 'run_pipeline.sbatch'}")
        for rep in range(2, args.replicates + 1):
            print(f"    rep{rep}: sbatch --dependency=afterok:$MAIN_JOB {folder / f'run_pull_rep{rep}.sbatch'}")
    if args.replicates > 1:
        print(f"\nTo submit all replicates automatically:")
        print(f"  for MUT_DIR in {' '.join(str(f) for f in created)}; do")
        print(f"    MAIN_JOB=$(sbatch --parsable $MUT_DIR/run_pipeline.sbatch)")
        for rep in range(2, args.replicates + 1):
            print(f"    sbatch --dependency=afterok:$MAIN_JOB $MUT_DIR/run_pull_rep{rep}.sbatch")
        print(f"  done")


if __name__ == "__main__":
    main()
