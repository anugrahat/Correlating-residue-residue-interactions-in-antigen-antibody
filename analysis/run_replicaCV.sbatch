#!/bin/bash
#SBATCH --job-name=repCV
#SBATCH --partition=high
#SBATCH --account=fallerlab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --output=/home/anugraha/c1_WT/analysis/replicaCV.out
#SBATCH --chdir=/home/anugraha/c1_WT/analysis

python3 -u run_regression_replicaCV.py --no-figs

echo "[DONE] Replica CV regression complete"
