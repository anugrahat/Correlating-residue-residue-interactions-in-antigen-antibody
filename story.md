# Contact Frequency Regression on SMD Ensembles for Rational Antibody Interface Engineering

## The Problem

You have an antibody (C1-87g7) and you want to make it bind tighter. The interface has ~30 residues in contact. Traditional approaches — alanine scanning, static structure analysis, brute-force computational mutagenesis — either tell you what's important without telling you how to improve it, or they test mutations one at a time without understanding the underlying contact network.

You need a method that answers two questions simultaneously: **which residues matter**, and **are they helping or hurting?**

## The Approach

You run steered MD — you mechanically pull the antibody off the antigen across a panel of mutants (10 mutations × 10 replicas = 100 trajectories). This is powerful because the pulling force reveals which contacts resist unbinding under mechanical stress. You extract pairwise contact frequencies between antigen and antibody residues across all frames, then use elastic net regression to ask: **which contacts predict rupture force?**

The elastic net gives you sparse, interpretable coefficients. A negative β means that contact forming more frequently correlates with stronger binding (favorable). Positive β means it correlates with weaker binding (unfavorable). You then aggregate all coefficients involving a given residue into a single **Net Favorability (NF)** score — a per-residue summary of how much it contributes to holding the complex together.

## The Overfitting Problem

Your first regression gives R² = 0.998. That number should make any experienced scientist uncomfortable. Contact frequencies and rupture forces come from the same trajectories. Frames within a replica are autocorrelated. You're fitting noise.

The fix: **leave-2-replicas-out cross-validation**. Train on 8 replicas, test on 2. 45 folds. This breaks the circularity because the test data comes from completely independent simulations.

Result: Test R² = 0.959. Still excellent — the signal is real. But the **rankings change in ways that matter**:

- **N248**: jumps from negligible (NF=2.6) to #4 hotspot (NF=45.7)
- **H225**: flips from unfavorable to favorable #5 (NF=43.4)
- **Y404**: flips from favorable #4 to unfavorable
- **Y405**: jumps from 32.9 to 61.0 — much more important than originally estimated

The top 3 (Y406, Y227, Y405) are stable across both models. The mid-ranked residues are where overfitting distorts the picture.

## The Design Logic

NF tells you WHERE. Chemistry tells you HOW. The regression classifies the interface into three tiers:

**Tier 1 — Hot spots** (Y406, Y227, Y405): High NF, aromatic, already making strong contacts. These are load-bearing. You don't redesign load-bearing walls — you either leave them alone or make conservative extensions. Y→W adds one ring to the indole, slightly extending the aromatic surface without changing the fundamental interaction geometry.

**Tier 2 — Warm spots** (N248, D245, F287, S262, A246, S247): Moderate NF, making real contributions but with suboptimal chemistry. These are your engineering targets. A serine making a contact could be a glycine (remove steric clash) or a tyrosine (add aromatic stacking). An aspartate could be an asparagine (remove charge penalty in a hydrophobic pocket). A phenylalanine could be a tryptophan (extend the aromatic surface).

**Tier 3 — Unfavorable/weak** (Y404 in CV): Contacts that actually detract from binding or contribute noise. Leave alone or consider removal.

## The Validation

You run MM-GBSA on 14 variants (12 designed + WT + alanine control). The results:

### Hits (ΔΔG < -5 kcal/mol)

| Mutation | ΔΔG (kcal/mol) | Logic |
|----------|----------------|-------|
| F287W | -33.7 | Warm spot, Phe→Trp aromatic extension |
| S247G | -31.6 | Warm spot, remove suboptimal side chain |
| D245N | -18.2 | Warm spot, neutralize buried charge |
| Y406W | -15.7 | Hot spot, conservative Tyr→Trp extension |
| S262Y | -15.0 | Warm spot, add aromatic contact |
| N248Y | -10.9 | Warm spot, add aromatic where CV says it matters |
| Y227W | -10.7 | Hot spot, conservative extension |
| A246Y | -6.1 | Warm spot, add aromatic surface |

### Failures

| Mutation | ΔΔG (kcal/mol) | Why it failed |
|----------|----------------|---------------|
| Y404W | -2.7 | CV says Y404 is unfavorable — extending it doesn't help |
| Y404H | +2.5 | Replacing aromatic with imidazole at unfavorable position |
| Y157A | +11.6 | Alanine control — removing a real contact |
| H225G | +20.0 | **The smoking gun** |

## The Smoking Gun: H225G

This is the single most important data point in the study.

The overfitted model says H225 is unfavorable (positive β, NF=5.8). If you trust that model, H225→Gly should help — you're removing a bad contact. The CV model says H225 is favorable and #5 most important (NF=43.4). If you trust that model, H225→Gly should be catastrophic.

**Result: ΔΔG = +20.0 kcal/mol.** Binding virtually abolished. The CV model was right. The overfitted model would have sent you in exactly the wrong direction.

This isn't a subtle effect. This is a 20 kcal/mol error in prediction. For context, your best hit (F287W) improved binding by 33.7. H225G destroyed it by 20. If you'd trusted the naive regression and spent wet-lab resources on H225G, you'd have wasted months on a mutation that obliterates binding.

## The Bigger Picture

The hit rate is 8/12 (67%) for rationally designed mutations, with a best hit of ΔΔG = -33.7 kcal/mol. But the real contribution isn't the hit rate — it's the framework:

1. **SMD contact regression** identifies mechanically relevant contacts across a mutation panel
2. **Net Favorability** converts regression coefficients into actionable per-residue scores
3. **Replica-based CV** is essential — without it, you get qualitatively wrong predictions for mid-ranked residues
4. **Chemistry-guided design** translates NF rankings into specific mutations based on the local structural context

The method is generalizable to any protein-protein interface where you can run SMD. And the H225G result is the definitive proof that the CV step isn't optional — it's the difference between a predictive tool and an expensive way to generate wrong answers.
