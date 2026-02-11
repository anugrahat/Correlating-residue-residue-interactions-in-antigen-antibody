#!/usr/bin/env bash
set -euo pipefail

BASE_PDB="/home/anugraha/antibody_optimization/c1_87g7.pdb"
SCRIPT="/home/anugraha/antibody_optimization/scripts/setup_smd_campaign.py"

MUTS=(
  "Y405N"
  "Y405Q"
  "Y406W"
  "Y404W"
  "F287W"
  "D265N"
  "S262Y"
)

python "$SCRIPT" --base-pdb "$BASE_PDB" --mutations "${MUTS[@]}"

echo
echo "Folders created. Submit jobs with:"
for m in "${MUTS[@]}"; do
  echo "  sbatch /home/anugraha/c1_${m}/run_pipeline.sbatch"
done
