# C1 Antibody Mutation Shortlist (From Existing Elastic Net Outputs)

## Data source used
- `/home/anugraha/omicron_pulling/bigger_box/50ns/replica_50ns_pos/output_data/new_analysis/elastic_net_coefficients_and_Fmax.csv`
- `/home/anugraha/omicron_pulling/bigger_box/50ns/replica_50ns_pos/output_data/new_analysis/max_contribution_beta_x_fmax.csv`
- Aggregated table saved to:
  - `/home/anugraha/antibody_optimization/c1_from_existing_regression_antibody_targets.csv`

No new regression was fit. This uses your existing C1 regression output files in the Omicron project folder.

## Interpretation rule used
- Negative `Beta_x_Fmax` = contact increase is associated with more favorable interaction energy.
- Antibody residues ranked by summed favorable contributions across all antigen-contact pairs.

## Top antibody hotspot residues (ranked)
1. 406 (TYR)
2. 227 (TYR)
3. 404 (TYR)
4. 405 (TYR)
5. 287 (PHE)
6. 249 (ARG)
7. 374 (ARG)
8. 262 (SER)
9. 251 (THR)
10. 265 (ASP)
11. 246 (ALA)
12. 289 (TRP)

## First mutation panel to test
Primary (highest priority):
- `Y406W`
- `Y405W`
- `Y404W`
- `F287W`
- `S262Y`
- `A246Y`

Secondary (conservative/alternative chemistry):
- `Y406F`
- `Y405F`
- `Y404F`
- `F287Y`
- `T251Y`
- `D265E`

## Why these
- Strongest regression signal is concentrated at an antibody aromatic cluster around `Y404/Y405/Y406` with large negative pair contributions.
- Additional favorable signal appears at `F287`, `S262`, `A246`, suggesting potential to improve packing/contact persistence by adding aromatic or larger side chains.
- Conservative alternatives are included to reduce steric-risk while probing mechanism.

## Constraints
- Antibody optimization only:
  - mutate residues in `196-422` only.
- Do not mutate antigen-side RBD residues (`1-195`) for this objective.

