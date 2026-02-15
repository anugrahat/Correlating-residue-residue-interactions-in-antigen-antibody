# Contact Frequency Regression on SMD Ensembles for Rational Antibody Interface Engineering

## Overview

We developed a computational framework for antibody affinity maturation that uses elastic net regression on steered molecular dynamics (SMD) pulling simulations to rank antibody interface residues by importance. The regression decomposes binding into per-residue-pair contributions, and a per-residue metric — Net Favorability (NF) — tells us which positions matter most. Chemistry-guided mutations at NF-ranked positions were validated by MM-GBSA, with 6 out of 11 mutations at top-12 residues improving binding, including a best hit of ΔΔG = -33.7 kcal/mol.

## Method

### 1. Steered MD pulling simulations

The C1-87g7 antibody–Wuhan RBD complex was pulled apart in 10 independent replicas (5 ns each, 0.001 nm/ps pull rate). During pulling, two quantities were recorded per frame:

- **Interaction energy** E(t) = Coulomb_SR + LJ_SR between antigen and antibody
- **Contact frequencies** F_k(t) = atom-atom pair counts within 10 Å for each antigen–antibody residue pair k

Surface residues identified by SASA > 0.5 nm² on the bound structure.

### 2. Elastic net regression

The regression asks: which contacts predict interaction energy during unbinding?

**Equation**: ΔE(t) = Σ_k β_k × ΔF_k(t)

- ΔE(t) = change in interaction energy at frame t relative to frame 0
- ΔF_k(t) = change in contact frequency for pair k at frame t relative to frame 0
- β_k = elastic net coefficient (sparse, L1/L2 regularized)
- No intercept — when all contacts are at reference values, ΔE = 0

**Cross-validation**: Leave-2-replicas-out (45 folds). Train on 8 replicas, test on 2. Subsampled every 30 frames to reduce autocorrelation (50 frames per replica, 500 total). Fixed hyperparameters: α = 12.92, l1_ratio = 0.90.

- Train R² = 0.989 ± 0.001
- **Test R² = 0.959 ± 0.027**
- Test RMSE = 37.9 ± 15.8 kJ/mol
- 75 robust features (selected in >50% of folds) out of 516 total pairs

### 3. Net Favorability (NF)

Per-residue importance score:

**NF_j = Σ_{k ∈ P_j} |β_k × F̄_k|**

- P_j = set of all robust contact pairs involving antibody residue j
- F̄_k = mean contact frequency for pair k across all replicas and frames
- Higher NF = more important at the interface

This collapses pair-level coefficients into one number per antibody residue. 26 interface residues identified, ranked from Y406 (NF=72.5) to N288 (NF=0.5).

### 4. NF Ranking (top 12)

| Rank | Residue | NF |
|------|---------|-----|
| 1 | Y406 | 72.5 |
| 2 | Y227 | 65.9 |
| 3 | Y405 | 61.0 |
| 4 | N248 | 45.7 |
| 5 | H225 | 43.4 |
| 6 | Y404 | 30.0 |
| 7 | D245 | 29.6 |
| 8 | A255 | 27.2 |
| 9 | F244 | 23.4 |
| 10 | F287 | 17.2 |
| 11 | T403 | 14.1 |
| 12 | W289 | 13.3 |

## Design Rules

NF tells you WHERE a residue sits in the importance hierarchy. Chemistry tells you WHAT to do about it. Two rules:

**Rule 1: High-NF residues — improve their chemistry, never remove them.**
Extend aromatics (Y→W), add aromatics (S→Y, N→Y), neutralize buried charges (D→N). Stripping a high-NF residue to Gly is catastrophic.

**Rule 2: Low-NF residues — safe to remove.**
These contribute little. Stripping to Gly clears steric clutter without losing meaningful contacts.

## Validation: MM-GBSA

WT baseline: ΔG = -25.31 ± 0.86 kcal/mol. Method: gmx_MMPBSA v1.6.3, igb=5, 0.15 M salt, 90 frames from 10 ns equilibrium MD.

### Mutations at top-12 NF residues

| Mutation | NF Rank | NF | ΔΔG (kcal/mol) | Result |
|----------|---------|-----|-----------------|--------|
| F287W | 10 | 17.2 | **-33.7** | Hit |
| D245N | 7 | 29.6 | **-18.2** | Hit |
| Y406W | 1 | 72.5 | **-15.7** | Hit |
| N248Y | 4 | 45.7 | **-10.9** | Hit |
| Y227W | 2 | 65.9 | **-10.7** | Hit |
| Y404W | 6 | 30.0 | -2.7 | Marginal |
| Y404H | 6 | 30.0 | +2.5 | Failure |
| T403Y | 11 | 14.1 | +4.2 | Failure |
| Y157A | — | — | +11.6 | Neg. control |
| Y404G | 6 | 30.0 | +15.3 | Failure |
| H225G | 5 | 43.4 | +20.0 | Failure |

### Key observations

**Biggest gains come from suboptimal chemistry at important positions.** F287 (NF=17.2, #10) gives the best ΔΔG (-33.7) because Phe→Trp extends the aromatic contact surface. Y406 (NF=72.5, #1) only gives -15.7 because Tyr is already a good aromatic. The improvement ceiling depends on how suboptimal the current chemistry is.

**Stripping high-NF residues is catastrophic.** H225G (NF=43.4): ΔΔG = +20.0. Y404G (NF=30.0): ΔΔG = +15.3. NF magnitude predicts how badly removal will fail.

**Y404 is structurally coupled.** Y404 sits in the aromatic triad (Y404/Y405/Y406). Even conservative substitutions (Y404W: -2.7, Y404H: +2.5) fail because any modification at this position disrupts the neighboring hot spots.

**Y157A validates the framework.** Alanine substitution at a known antigen-side binding residue (Y157 on RBD) gives ΔΔG = +11.6. Y157 appears naturally in the regression — the pair Y157–Y405 has β × freq = -14.7, one of the largest contributions. Disrupting the antigen side of this contact confirms the regression identified real binding determinants.

## Figures

| Figure | File | Description |
|--------|------|-------------|
| Regression equation | `fig_regression_equation.pdf` | ΔE = Σ β_k × ΔF_k with variable definitions |
| NF equation | `fig_nf_equation.pdf` | NF_j = Σ \|β_k × F̄_k\| with variable definitions |
| NF ranking | `fig_nf_antibody.pdf` | All 26 antibody interface residues ranked by NF |
| ΔΔG validation | `fig_ddg.pdf` | Mutations at top-12 NF residues + Y157A control |
| β × mean freq | `fig_beta_x_meanfreq.pdf` | Top 25 residue pairs by \|β × mean freq\|, Y157 in red |
| Temporal β × freq | `fig_temporal_beta_x_freq.pdf` | β × freq(t) during 5 ns pull for key pairs |
| Pull force curves | `fig_pull_forces.pdf` | Smoothed pull force profiles for all validated mutations |
| Peak rupture force | `fig_peak_force.pdf` | Peak force bar chart |

## File Paths

### Data
- Per-replica contact frequencies: `/home/anugraha/c1_WT/analysis/rep{1-10}/frame_residuepair_interaction_frequencies.dat`
- Averaged contact frequencies: `/home/anugraha/c1_WT/analysis/average_frequency.csv`
- Regression coefficients: `/home/anugraha/c1_WT/analysis/elastic_net_coefficients_replicaCV.csv`
- NF ranking: `/home/anugraha/c1_WT/analysis/net_favorability_replicaCV.csv`
- CV fold stats: `/home/anugraha/c1_WT/analysis/replicaCV_fold_stats.csv`
- MM-GBSA results: `/home/anugraha/c1_*/work/FINAL_RESULTS_MMPBSA.dat`

### Scripts
- Regression: `/home/anugraha/c1_WT/analysis/run_regression_replicaCV.py`
- Contact frequency: `/home/anugraha/c1_WT/run_contact_freq.py`
- Averaging: `/home/anugraha/c1_WT/analysis/average_contact_freq.py`
- Mutation pipeline: `/home/anugraha/antibody_optimization/scripts/setup_smd_campaign.py`
- Auto-align: `/home/anugraha/antibody_optimization/scripts/auto_align_for_pull.py`
- Plotting: `/home/anugraha/antibody_optimization/scripts/plot_*.py`
- Figures: `/home/anugraha/antibody_optimization/figures/`
