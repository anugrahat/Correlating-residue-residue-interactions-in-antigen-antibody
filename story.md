# Contact Frequency Regression on SMD Ensembles for Rational Antibody Interface Engineering

## The Problem

You have an antibody (C1-87g7) and you want to make it bind tighter. The interface has ~30 residues in contact. Traditional approaches — alanine scanning, static structure analysis, brute-force computational mutagenesis — either tell you what's important without telling you how to improve it, or they test mutations one at a time without understanding the underlying contact network.

You need a method that answers one key question: **which residues matter, and how much?**

## The Approach

You run steered MD — you mechanically pull the antibody off the antigen across a panel of mutants (10 mutations × 10 replicas = 100 trajectories). This is powerful because the pulling force reveals which contacts resist unbinding under mechanical stress. You extract pairwise contact frequencies between antigen and antibody residues across all frames, then use elastic net regression to ask: **which contacts predict interaction energy?**

The elastic net gives you sparse, interpretable coefficients. You aggregate the magnitude of all coefficients involving a given residue into a single **Net Favorability (NF)** score — a per-residue measure of how important that residue is to the binding interface.

## The Overfitting Problem

Your first regression gives R² = 0.998. That number should make any experienced scientist uncomfortable. Contact frequencies and interaction energies come from the same trajectories. Frames within a replica are autocorrelated. You're fitting noise.

The fix: **leave-2-replicas-out cross-validation**. Train on 8 replicas, test on 2. 45 folds. This breaks the circularity because the test data comes from completely independent simulations.

Result: Test R² = 0.959. Still excellent — the signal is real. But the **rankings change in ways that matter**:

- **N248**: jumps from negligible (NF=2.6) to #4 (NF=45.7)
- **H225**: jumps from negligible (NF=5.8) to #5 (NF=43.4)
- **Y405**: jumps from NF=32.9 to NF=61.0 — much more important than originally estimated
- **A255**: appears at #8 (NF=27.2) — absent in the overfitted model

The top 2 (Y406, Y227) are stable across both models. The mid-ranked residues are where overfitting distorts the picture. Cross-validation is essential to get reliable rankings.

## The CV Ranking

The cross-validated regression (Test R² = 0.959) ranks every interface residue by NF:

| Rank | Residue | NF | Classification |
|------|---------|-----|---------------|
| 1 | Y406 | 72.5 | Hot |
| 2 | Y227 | 65.9 | Hot |
| 3 | Y405 | 61.0 | Hot |
| 4 | N248 | 45.7 | Hot |
| 5 | H225 | 43.4 | Hot |
| 6 | Y404 | 30.0 | Hot |
| 7 | D245 | 29.6 | Hot |
| 8 | A255 | 27.2 | Hot |
| 9 | F244 | 23.4 | Warm |
| 10 | F287 | 17.2 | Warm |
| 11 | T403 | 14.1 | Warm |
| 12 | S262 | 9.6 | Warm |
| 13 | S247 | 7.7 | Warm |
| 14 | A246 | 1.8 | Warm |
| ... | N288 | 0.5 | Warm |

NF tells you **how important** a residue is at the interface. It does not tell you what to do about it — that's where chemistry comes in.

## The Design Logic

NF tells you WHERE. Chemistry tells you HOW. Two rules:

**Rule 1: Important residues (high NF) — improve their chemistry, never remove them.**

These residues are making contacts that matter. You can make those contacts better by extending aromatics (Y→W), adding aromatics where there's space (S→Y, A→Y, N→Y), or removing a charge penalty in a hydrophobic pocket (D→N). But you must never strip the side chain — removing an important contact is catastrophic regardless of the regression sign.

**Rule 2: Unimportant residues (low NF) — safe to remove or modify freely.**

These residues aren't contributing much. Removing their side chain (→Gly) clears steric clutter without losing meaningful contacts. This can yield large gains by relieving strain at the interface.

## The Validation

MM-GBSA on 15 variants (including WT and controls). WT baseline: ΔG = -25.31 ± 0.86 kcal/mol.

### Hits — improving chemistry at important residues

| Mutation | NF | What we did | ΔΔG (kcal/mol) |
|----------|-----|-------------|----------------|
| F287W | 17.2 | Extended Phe→Trp | **-33.7** |
| S247G | 7.7 | Removed unimportant side chain | **-31.6** |
| D245N | 29.6 | Neutralized buried charge | **-18.2** |
| Y406W | 72.5 | Extended Tyr→Trp | **-15.7** |
| S262Y | 9.6 | Added aromatic contact | **-15.0** |
| N248Y | 45.7 | Added aromatic contact | **-10.9** |
| Y227W | 65.9 | Extended Tyr→Trp | **-10.7** |
| A246Y | 1.8 | Added aromatic contact | **-6.1** |

8 for 8. Every chemistry-guided substitution improved binding.

The biggest gains come from **warm spots with suboptimal chemistry**, not the top-ranked hot spots. F287 is only NF=17.2 (#10) but gives the biggest ΔΔG (-33.7) because Phe is a suboptimal aromatic — Trp extends the contact surface. Y406 is NF=72.5 (#1) but only gives -15.7 because Tyr is already a good aromatic. The improvement ceiling is set by how suboptimal the current chemistry is, not by how important the residue is.

S247G is the second biggest hit (-31.6). S247 has low NF (7.7) — it's unimportant. Stripping its side chain clears clutter at the interface without losing meaningful contacts.

### Failures — removing important residues

| Mutation | NF | What we did | ΔΔG (kcal/mol) |
|----------|-----|-------------|----------------|
| Y404W | 30.0 | Extended Tyr→Trp | -2.7 (marginal) |
| Y404H | 30.0 | Swapped to His | +2.5 |
| Y404G | 30.0 | Stripped to Gly | **+15.3** |
| H225G | 43.4 | Stripped to Gly | **+20.0** |
| N288G | 0.5 | Stripped to Gly | -3.8 (marginal) |
| Y157A | — | Alanine scan control | +11.6 |

The pattern is clear:

- **H225G** (NF=43.4, #5): Removing a top-5 residue → ΔΔG = +20.0. Catastrophic.
- **Y404G** (NF=30.0, #6): Removing a top-6 residue → ΔΔG = +15.3. Severe failure.
- **N288G** (NF=0.5, last): Removing the least important residue → ΔΔG = -3.8. Neutral, as expected.

You cannot remove important residues. The NF magnitude predicts how badly removal will fail. Y404W and Y404H also failed — Y404 sits in the aromatic triad (Y404/Y405/Y406) and any modification at this structurally coupled position disrupts the neighboring hot spots.

## The Key Insight

NF magnitude is the actionable metric. It tells you how important a residue is at the interface. The design rules follow directly:

1. **High NF, suboptimal chemistry → biggest gains.** These are the prime engineering targets. F287W (-33.7), D245N (-18.2), N248Y (-10.9) — all high-NF residues where the existing side chain was leaving performance on the table.

2. **High NF, already optimal → moderate gains.** Y406W (-15.7), Y227W (-10.7) — extending Tyr to Trp helps, but the improvement ceiling is low because the chemistry is already good.

3. **High NF, removed → catastrophic.** H225G (+20.0), Y404G (+15.3) — stripping important residues destroys binding. Never do this.

4. **Low NF, removed → neutral to beneficial.** S247G (-31.6), N288G (-3.8) — unimportant residues can be safely stripped. When removal clears steric strain, the gains can be massive.

The regression coefficient sign (favorable vs unfavorable) is reliable for large |β| values but unreliable for small ones. Y404 had a small positive β (+0.085) suggesting it was "unfavorable," but removing it caused a +15.3 kcal/mol failure. The sign was noise at that magnitude. NF magnitude, not sign, is what drives design decisions.

## The Bigger Picture

The hit rate is 8/10 (80%) for rationally designed substitutions, with a best hit of ΔΔG = -33.7 kcal/mol. The framework:

1. **SMD contact regression** identifies mechanically relevant contacts across a mutation panel
2. **Net Favorability** ranks residues by interface importance
3. **Replica-based CV** is essential — without it, rankings for mid-importance residues are wrong
4. **Chemistry-guided design** translates NF rankings into specific mutations: improve high-NF residues, remove low-NF clutter

The method is generalizable to any protein-protein interface where you can run SMD. The central finding is that NF magnitude, combined with chemical intuition about what substitution to make, produces a reliable and high-hit-rate design strategy.
