# Regression-Guided Antibody Affinity Maturation via SMD Contact Analysis

## Overview

We developed a computational pipeline for antibody affinity maturation that uses elastic net regression on steered molecular dynamics (SMD) pulling simulations to identify mutation targets. The regression decomposes the binding interface into residue-pair contributions, classifying each antibody residue as a hot spot (don't touch), warm spot (improve existing contacts), or cold spot (create new contacts). Mutations designed from this classification were validated using MM-GBSA on independent equilibrium MD trajectories, yielding a top single mutant (D265N) with a predicted ΔΔG of -39.4 kcal/mol improvement over wild type.

## Step 1: SMD pulling simulations

The C1-87g7 antibody-antigen complex was subjected to steered molecular dynamics (SMD) pulling simulations in GROMACS. The antibody was restrained while the antigen (RBD) was pulled along the Y axis at a constant velocity of 0.001 nm/ps over 5 nm of separation. Four independent replicas were run from different equilibration starting points (0ns, 60ns, 80ns, 100ns pre-equilibrated structures) to sample different starting configurations of the bound complex.

During each pulling simulation, two quantities were recorded per frame:

**Interaction energy** (dependent variable for regression):
- Short-range Coulomb electrostatic energy (C_SR) between antigen and antibody groups
- Short-range Lennard-Jones energy (LJ_SR) between antigen and antibody groups
- Total interaction energy E(t) = C_SR(t) + LJ_SR(t), in kJ/mol
- Extracted from GROMACS energy output (.xvg files)

**Contact frequencies** (independent variables for regression):
- For each frame, all atom-atom pairwise distances were computed between surface residues of the antigen (group 1, residues 1-195) and antibody (group 2, residues 196-422)
- Surface residues were identified using Shrake-Rupley solvent accessible surface area (SASA) on the initial bound structure (frame 0 only), with a threshold of 0.5 nm² — only residues exposed to solvent were considered interface candidates
- Any atom-atom pair within 10 Angstroms was counted as an interaction
- For each residue pair (res_i on antigen, res_j on antibody), the contact frequency for that frame = total number of atom-atom pairs within the distance cutoff between those two residues
- This means contact frequency is NOT binary (yes/no contact) but a count reflecting how many atomic contacts exist — a large tyrosine touching many atoms on the partner scores higher than a small alanine touching one atom
- Script: `heat_count_surf.py` (MDAnalysis + mdtraj)

Both quantities were then **averaged frame-by-frame across the 4 replicas** to reduce noise:
- `interaction_energy.csv`: Frame, AverageInteractionEnergy (mean of 4 replicas)
- `average_frequency.csv`: Frame, ResiduePair, InteractionFrequency (mean of 4 replicas)
- Averaging script: `average3.py`
- Certain residue pairs involving boundary residue 196 were excluded as artifacts

## Step 2: Elastic net regression

The key idea: as the antibody is pulled away from the antigen, contacts break sequentially. Some contacts breaking causes a large energy loss (important contacts), others break with no energy change (irrelevant contacts). The regression identifies which contact pairs matter.

**Setup** (script: `average_11_avg_f_rbd_2.py`):

1. **Reference frame**: Frame 0 (fully bound state)
2. **Delta computation**:
   - ΔE(t) = E(t) - E(frame 0) — change in interaction energy relative to bound state
   - ΔF_k(t) = F_k(t) - F_k(frame 0) — change in contact frequency for pair k relative to bound state
3. **Regression model**: ΔE(t) ≈ Σ_k β_k × ΔF_k(t)
   - Each β_k tells how much energy changes when contact pair k changes by one unit
   - No intercept (fit_intercept=False) — when all contacts are at their reference values, ΔE = 0

**ElasticNetCV parameters**:
- alphas: logspace(-2, 2, 10) — regularization strength range
- l1_ratio: [0.3, 0.5, 0.7, 0.9] — balance between L1 (sparsity) and L2 (ridge) penalties
- cv: 5-fold cross-validation
- max_iter: 200,000
- Frame range: 0-3100 (50ns system) or 0-5000 (100ns system)

**Why elastic net?** The L1 penalty drives most coefficients to exactly zero (feature selection), keeping only the truly important contact pairs. The L2 penalty handles collinearity between neighboring residue pairs. The result is a sparse model where non-zero β values identify the critical contacts.

**Model performance**: R² ≈ 0.65, Pearson r ≈ 0.80 — the selected contact pairs explain ~65% of the variance in interaction energy during pulling.

## Step 3: Interpreting β coefficients

A β coefficient for contact pair (antigen_res_i, antibody_res_j) means:
- **Negative β**: When this contact frequency decreases (contact breaks during pulling), energy increases (becomes less negative = weaker binding). This is a FAVORABLE contact — it stabilizes binding.
- **Positive β**: When this contact breaks, energy decreases. This is an UNFAVORABLE contact — it destabilizes binding.
- **Zero β** (driven to zero by L1 penalty): This contact is not important for binding energy, even if it exists geometrically.

To rank antibody residues (not pairs), we aggregate across all antigen partners:
- **β × MeanFreq**: The product of the regression coefficient and the average contact frequency across frames. This weights both the coefficient strength and how often the contact exists.
- **Sum by antibody residue**: For each antibody residue j, sum |β_k × MeanFreq_k| across all antigen partners i. This gives the total energetic contribution of that antibody residue = "NetFavorability."

## Step 4: Hot / Warm / Cold spot classification

From the per-residue NetFavorability and number of contact pairs:

| Classification | NetFavorability | #Pairs | Meaning | Mutation strategy |
|---|---|---|---|---|
| **Hot spot** | > 70 | any | Already optimal, dominant binding contributor | DO NOT MUTATE — any change destroys binding |
| **Warm spot** | 20-70 | moderate | Makes contacts but suboptimally | IMPROVE existing contacts — charge removal, conservative aromatic substitutions |
| **Cold spot** | ≈ 0 | ≥ 8 | Geometrically at interface but zero energetic contribution | CREATE new contacts — add bulky/aromatic side chains to convert dead contacts to productive ones |
| **Lukewarm** | 5-20 | high | Weak signal despite many geometric contacts | Risky — geometry may be wrong for this position |
| **Destabilizing** | < 0 | any | Actively hurts binding | Remove/shrink the offending side chain |

**The critical insight**: Traditional alanine scanning identifies hot spots (residues where mutation to alanine causes large ΔΔG). But it cannot identify cold spots — residues that contribute nothing but have untapped potential. Only the combination of zero regression signal + many geometric contact pairs reveals these targets. This is the novel contribution of the regression approach.

## Step 5: Mutation design

Based on the classification, we designed a panel of mutations:

**Hot spot controls (don't mutate):**
- Y405N (H:103, NetFav=86.0) — confirmed: -29% peak force vs WT. Any substitution at this position weakens binding.

**Warm spot mutations (improve suboptimal contacts):**
- D265N (L:70, NetFav=24.3, 18 pairs): Remove negative charge (Asp→Asn) to improve electrostatic complementarity
- S262Y (L:67, NetFav=50.8, 24 pairs): Add aromatic (Ser→Tyr) to extend contact surface
- A246Y (L:51, NetFav=23.6, 24 pairs): Add aromatic at small warm site

**Hot spot conservative substitution:**
- F287W (L:92, NetFav=76.7, 20 pairs): Keep aromatic but extend to indole (Phe→Trp). Not a new chemistry — just more of the same productive interaction.

**Cold spot mutation (create new contacts):**
- N248Y (L:53, NetFav=0, 17 pairs): Zero regression signal but 17 contact pairs. Asparagine is too small and polar to form productive hydrophobic contacts. Tyrosine should create new interactions.

**Lukewarm position:**
- T403Y (H:101, NetFav=8.0, 22 pairs): 22 contact pairs but only 8.0 net favorability — massive untapped potential, or wrong geometry?

**Negative control:**
- Y157A: Alanine substitution at a known antigen-side binding hotspot (Tyr157 on RBD). Should dramatically weaken binding if our assay is sensitive.

**Double mutant:**
- D265N+F287W: Combine the two best-performing singles from different interface regions (L:70 and L:92)

## Step 6: Validation by SMD pulling (work integrals)

All mutations were run through the same SMD pulling pipeline with 6-10 replicates each. Work integrals W = ∫F·dx were computed by trapezoidal integration of force-distance curves.

| Mutation | Mean W (kJ/mol) | n | vs WT (%) | Classification |
|----------|-----------------|---|-----------|----------------|
| D265N | 715 ± 34 | 10 | +90% | warm |
| F287W | 551 ± 27 | 10 | +46% | hot (conservative) |
| T403Y | 506 ± 24 | 6 | +34% | lukewarm |
| S262Y | 478 ± 39 | 6 | +27% | warm |
| WT | 377 ± 22 | 10 | — | baseline |
| Y157A | 332 ± 27 | 6 | -12% | neg. control |
| A246Y | 312 ± 36 | 6 | -17% | warm |

**SMD limitation**: The noise floor is ±15-20%. Only D265N (+90%) and F287W (+46%) are reliably resolved. Mid-tier mutations overlap with WT variance. Y157A (negative control) showed only -12%, not statistically significant — SMD cannot resolve moderate binding differences from fast nonequilibrium pulling.

## Step 7: Validation by MM-GBSA (equilibrium binding free energy)

To overcome SMD noise limitations, we computed MM-GBSA binding free energies on the existing 10ns equilibrium MD trajectories (Phase 1 of the pipeline). This uses a completely independent method — no pulling involved.

**Method**: gmx_MMPBSA v1.6.3, igb=5 (Generalized Born), 0.15 M salt, 90 frames from 10ns trajectories, single-trajectory approach.

| Mutation | ΔG_bind (kcal/mol) | SEM | ΔΔG vs WT | Classification |
|----------|-------------------|-----|-----------|----------------|
| D265N | -64.73 | 1.45 | -39.4 | warm |
| F287W | -58.99 | 1.44 | -33.7 | hot (conservative) |
| D265N+F287W | -50.82 | 0.88 | -25.5 | double |
| S262Y | -40.28 | 1.20 | -15.0 | warm |
| N248Y | -36.20 | 1.01 | -10.9 | cold |
| A246Y | -31.41 | 0.81 | -6.1 | warm |
| **WT** | **-25.31** | **0.86** | **—** | **baseline** |
| T403Y | -21.08 | 1.57 | +4.2 | lukewarm |
| Y157A | -13.73 | 1.49 | +11.6 | neg. control |

Additional double mutants tested:
| Mutation | ΔG_bind (kcal/mol) | SEM | ΔΔG vs WT |
|----------|-------------------|-----|-----------|
| D265N+N248Y | -45.58 | 1.00 | -20.3 |
| D265N+S262Y | -35.37 | 0.82 | -10.1 |

F287W+S262Y and F287W+N248Y are currently running (jobs 754578, 754579).

## Step 8: Key findings

**1. Warm spots are the highest-yield targets.**
D265N (ΔΔG = -39.4 kcal/mol) is the best mutation found. The regression identified D265 as a warm spot — making contacts but suboptimally. Removing the aspartate's negative charge (Asp→Asn) dramatically improved electrostatic complementarity at the interface. This is the classic affinity maturation result: the residue is in the right place, just the wrong chemistry.

**2. Cold spots with many contact pairs are a novel target class.**
N248Y (ΔΔG = -10.9) validated the cold spot strategy. N248 had zero regression signal despite 17 geometric contact pairs — it sits at the interface doing nothing. Traditional alanine scanning would never flag this residue because there is nothing to lose. But the regression + geometric analysis reveals the opportunity: adding a tyrosine converted dead contacts into productive ones.

**3. Conservative substitutions can improve even hot spots.**
F287W (ΔΔG = -33.7) shows that hot spots are not untouchable if the substitution extends rather than replaces the existing interaction. Phe→Trp keeps the aromatic platform but adds a larger indole ring with more van der Waals surface area.

**4. Lukewarm positions are risky.**
T403Y (ΔΔG = +4.2, slightly worse than WT) failed despite 22 contact pairs. Low NetFavorability + high contact count may indicate wrong geometry for productive contacts at this position. The SMD result (+34%) was misleading — T403Y resists mechanical separation (steric friction from bulky Tyr) but does not improve equilibrium binding.

**5. The negative control validates sensitivity.**
Y157A (ΔΔG = +11.6 kcal/mol, strongly reduced binding) confirms that MM-GBSA can detect binding disruption. This alanine substitution at a key antigen hotspot (Tyr157 on RBD) was unresolvable by SMD (only -12%, p=0.10) but clearly detected by MM-GBSA.

**6. Double mutants show systematic anti-cooperativity.**
D265N+F287W (-50.82) is worse than either D265N (-64.73) or F287W (-58.99) alone. Energy decomposition shows the electrostatic contribution drops from -86.80 (D265N) to -34.14 (double), indicating electrostatic interference. Subsequent doubles (D265N+S262Y: -35.37, D265N+N248Y: -45.58) showed the same pattern. The binding interface is a coupled electrostatic/steric network — local improvements don't combine additively.

## Summary: the story in one paragraph

Elastic net regression on SMD pulling trajectories decomposes antibody-antigen binding into per-residue-pair contributions by correlating frame-by-frame contact frequency changes with interaction energy changes during forced separation. This classifies interface residues into hot spots (already optimal), warm spots (suboptimal contacts improvable by targeted substitutions), and cold spots (geometrically positioned but energetically silent, representing untapped potential for new contacts). Applied to the C1-87g7 antibody, the approach identified D265N (warm spot charge removal, ΔΔG = -39.4 kcal/mol) and F287W (conservative hot spot extension, ΔΔG = -33.7 kcal/mol) as strong affinity-improving mutations, validated the cold spot concept through N248Y (ΔΔG = -10.9 kcal/mol from a position with zero energetic signal), and revealed that double mutants exhibit systematic anti-cooperativity due to the coupled nature of the binding interface. The negative control Y157A confirmed assay sensitivity with +11.6 kcal/mol binding loss. Critically, the cold spot target class — invisible to traditional alanine scanning because there is nothing to lose — represents a novel avenue for affinity maturation that only emerges from the combination of regression-based energetic decomposition and geometric contact analysis.

## File Paths and Scripts

### Contact frequency computation
- Script: `/home/anugraha/c1_87g7/rot_120_corrected/bigger_box/md2/eq80/replica_80ns_pos/heat_count_surf.py`
- Also at: `/home/anugraha/omicron_pulling/bigger_box/50ns/replica_50ns_pos/heat_count_surf.py`

### Pulling simulation replicas (used for regression)
- Replica 1 (eq0ns): `/home/anugraha/omicron_pulling/bigger_box/50ns/replica_50ns_pos/`
  - Contact data: `frame_residuepair_interaction_frequencies.dat`
  - Energy: `/home/anugraha/c1_87g7/rot_120_corrected/bigger_box/replica_0ns_pos/C_SR.xvg`, `LJ_SR.xvg`
- Replica 2 (eq80ns): `/home/anugraha/c1_87g7/rot_120_corrected/bigger_box/md2/eq80/replica_80ns_pos/`
  - Contact data: `output_data/frame_residuepair_interaction_frequencies.dat`
  - Energy: `C_SR.xvg`, `LJ_SR.xvg`
- Replica 3 (eq60ns): `/home/anugraha/c1_87g7/rot_120_corrected/bigger_box/md2/eq60/replica_80ns_pos/`
  - Contact data: `output_data/frame_residuepair_interaction_frequencies.dat`
  - Energy: `C_SR.xvg`, `LJ_SR.xvg`
- Replica 4 (eq100ns): `/home/anugraha/c1_87g7/rot_120_corrected/bigger_box/md2/eq100/replica_100ns/`
  - Contact data: `output_data/frame_residuepair_interaction_frequencies.dat`
  - Energy: `C_SR.xvg`, `LJ_SR.xvg`

### Averaging script (4 replicas → averaged CSVs)
- Script: `/home/anugraha/omicron_pulling/bigger_box/50ns/replica_50ns_pos/output_data/average3.py`
- Output: `interaction_energy.csv`, `average_frequency.csv`

### Elastic net regression
- Script: `/home/anugraha/omicron_pulling/bigger_box/50ns/replica_50ns_pos/output_data/new_analysis/average_11_avg_f_rbd_2.py`
- Input: `interaction_energy.csv`, `average_frequency.csv` (in same directory)
- Output: `elastic_net_coefficients_and_meanFreq.csv`
- Per-residue sums: `sum_of_beta_meanFreq_by_firstResidue.csv`
- Aggregated antibody targets: `/home/anugraha/antibody_optimization/c1_from_existing_regression_antibody_targets.csv`

### New C1 WT pulling simulations (10 replicates)
- Rep 1: `/home/anugraha/c1_WT/pull/` (replica2_pullf.xvg, replica2_pullx.xvg, replica2.edr)
- Rep 2-10: `/home/anugraha/c1_WT/pull_rep{2..10}/` (same files)
- Pipeline sbatch: `/home/anugraha/c1_WT/run_pipeline.sbatch`

### MM-GBSA results (per mutation)
- Results: `/home/anugraha/c1_*/work/FINAL_RESULTS_MMPBSA.dat`
- PBC-corrected trajectories: `/home/anugraha/c1_*/work/md10_noPBC.trr`
- Index files: `/home/anugraha/c1_*/work/mmpbsa_index.ndx`
- Input: `/home/anugraha/c1_*/work/mmpbsa.in`

### Mutation pipeline
- Setup script: `/home/anugraha/antibody_optimization/scripts/setup_smd_campaign.py`
- Auto-align: `/home/anugraha/antibody_optimization/scripts/auto_align_for_pull.py`
- Base PDB: `/home/anugraha/antibody_optimization/c1_87g7.pdb`

### Analysis plots
- SMD work integral plots: `/home/anugraha/antibody_optimization/plot1-7_*.png`
- Y157A comparison: `/home/anugraha/antibody_optimization/y157a_full_comparison.png`
- Exit seminar figures: `/home/anugraha/antibody_optimization/figures/`
