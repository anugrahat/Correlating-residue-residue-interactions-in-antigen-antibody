#!/usr/bin/env python3
"""
Elastic net regression with replica-based cross-validation.

Key differences from run_full_regression_3000.py:
  1. Loads all 10 replicas SEPARATELY (not averaged)
  2. Subsamples every N frames to reduce autocorrelation
  3. Leave-2-replicas-out CV: train on 8, test on 2
  4. Reports TEST-SET R² (honest generalization metric)
  5. Only reports coefficients consistently selected across CV folds
"""

import os
import sys
import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNetCV, ElasticNet
from sklearn.metrics import r2_score, mean_squared_error
from itertools import combinations

MAKE_FIGS = '--no-figs' not in sys.argv

# ==============================================================================
# Config
# ==============================================================================
BASE_DIR = "/home/anugraha/c1_WT/analysis"
FIG_DIR = "/home/anugraha/antibody_optimization/figures"
GRO_PATH = "/home/anugraha/c1_WT/pull/md10_protein.gro"

FRAME_END = 1500
REF_FRAME = 0
SUBSAMPLE_STEP = 30  # take every 30th frame to reduce autocorrelation
N_REPLICAS = 10

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIE": "H", "HIS": "H",
    "HID": "H", "HIP": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T",
    "TRP": "W", "TYR": "Y", "VAL": "V",
}


def load_gro_resnames(gro_path):
    resnames = {}
    with open(gro_path) as f:
        lines = f.readlines()
    for line in lines[2:]:
        if len(line) < 20:
            continue
        try:
            resid = int(line[:5].strip())
            resname = line[5:10].strip()
            if resid not in resnames:
                resnames[resid] = resname
        except (ValueError, IndexError):
            continue
    return resnames


def gro_to_global(gro_resid):
    if 1 <= gro_resid <= 194:
        return gro_resid
    elif 195 <= gro_resid <= 314:
        return gro_resid - 194 + 302
    elif 315 <= gro_resid <= 421:
        return gro_resid - 314 + 195
    raise ValueError(f"GRO resid {gro_resid} out of range")


def get_label(gro_resid, resnames):
    glob = gro_to_global(gro_resid)
    rn = resnames.get(gro_resid, "UNK")
    aa1 = AA3_TO_1.get(rn, "X")
    return f"{aa1}{glob}"


def is_antibody(gro_resid):
    return 195 <= gro_resid <= 421


def is_inter_contact(pair_str):
    g1, g2 = [int(x) for x in pair_str.split('-')]
    return not (g1 >= 195 and g2 >= 195)


def load_replica(rep_num):
    """Load energy and contact frequency for one replica."""
    rep_dir = os.path.join(BASE_DIR, f"rep{rep_num}")

    # Energy
    energy_file = os.path.join(rep_dir, "energy_chauA_rest.xvg")
    energy_data = []
    with open(energy_file) as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                continue
            parts = line.split()
            frame = int(float(parts[0]))
            # Sum Coulomb + LJ
            e_coul = float(parts[1])
            e_lj = float(parts[2])
            energy_data.append((frame, e_coul + e_lj))
    df_energy = pd.DataFrame(energy_data, columns=['Frame', 'Energy']).set_index('Frame')

    # Contact frequencies
    freq_file = os.path.join(rep_dir, "frame_residuepair_interaction_frequencies.dat")
    freq_data = []
    with open(freq_file) as f:
        header = f.readline()
        for line in f:
            parts = line.split()
            frame = int(parts[0])
            pair = parts[1]
            freq = float(parts[2])
            freq_data.append((frame, pair, freq))
    df_freq = pd.DataFrame(freq_data, columns=['Frame', 'ResiduePair', 'InteractionFrequency'])

    # Filter intra-antibody
    df_freq = df_freq[df_freq['ResiduePair'].apply(is_inter_contact)]

    # Pivot
    df_pivot = df_freq.pivot_table(
        index='Frame', columns='ResiduePair',
        values='InteractionFrequency', fill_value=0
    ).sort_index()

    return df_energy, df_pivot


def main():
    print("=" * 70)
    print("ELASTIC NET WITH REPLICA-BASED CROSS-VALIDATION")
    print(f"Subsampling every {SUBSAMPLE_STEP} frames, frames 1-{FRAME_END}")
    print("=" * 70)

    resnames = load_gro_resnames(GRO_PATH)

    # ------------------------------------------------------------------
    # 1) Load all replicas
    # ------------------------------------------------------------------
    print("\nLoading replicas...")
    replica_data = {}
    all_pairs = set()

    for rep in range(1, N_REPLICAS + 1):
        print(f"  Loading replica {rep}...")
        df_e, df_f = load_replica(rep)
        replica_data[rep] = (df_e, df_f)
        all_pairs.update(df_f.columns.tolist())

    # Use union of all pairs (fill missing with 0)
    all_pairs = sorted(all_pairs)
    print(f"  Total unique inter-contact pairs across replicas: {len(all_pairs)}")

    # ------------------------------------------------------------------
    # 2) Build per-replica ΔE and ΔF matrices (subsampled)
    # ------------------------------------------------------------------
    print("\nBuilding subsampled ΔE/ΔF matrices...")
    rep_X = {}
    rep_y = {}

    for rep in range(1, N_REPLICAS + 1):
        df_e, df_f = replica_data[rep]

        # Align columns
        df_f = df_f.reindex(columns=all_pairs, fill_value=0)

        # Common frames
        common = df_f.index.intersection(df_e.index)
        E = df_e.loc[common, 'Energy']
        F = df_f.loc[common]

        # Delta from ref
        if REF_FRAME not in E.index:
            print(f"  WARNING: rep {rep} missing ref frame, skipping")
            continue
        E_ref = E.loc[REF_FRAME]
        F_ref = F.loc[REF_FRAME]
        dE = (E - E_ref).loc[1:FRAME_END]
        dF = F.subtract(F_ref, axis='columns').loc[1:FRAME_END]

        # Subsample
        frames = dE.index[::SUBSAMPLE_STEP]
        dE_sub = dE.loc[frames]
        dF_sub = dF.loc[frames]

        rep_X[rep] = dF_sub.values
        rep_y[rep] = dE_sub.values
        print(f"  Rep {rep}: {len(frames)} subsampled frames")

    # ------------------------------------------------------------------
    # 3) Replica-based CV: leave-2-out
    # ------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("LEAVE-2-REPLICAS-OUT CROSS-VALIDATION")
    print("=" * 70)

    reps = sorted(rep_X.keys())
    test_combos = list(combinations(reps, 2))
    n_folds = len(test_combos)
    print(f"  {n_folds} folds (C(10,2) = 45)")

    # Use fixed hyperparameters (from previous full-data fit) for speed
    FIXED_ALPHA = 12.92
    FIXED_L1 = 0.90

    fold_r2_train = []
    fold_r2_test = []
    fold_rmse_test = []
    fold_coefs = []

    for fold_i, test_reps in enumerate(test_combos):
        train_reps = [r for r in reps if r not in test_reps]

        X_train = np.vstack([rep_X[r] for r in train_reps])
        y_train = np.concatenate([rep_y[r] for r in train_reps])
        X_test = np.vstack([rep_X[r] for r in test_reps])
        y_test = np.concatenate([rep_y[r] for r in test_reps])

        model = ElasticNet(
            alpha=FIXED_ALPHA, l1_ratio=FIXED_L1,
            fit_intercept=False, max_iter=200000,
        )
        model.fit(X_train, y_train)

        y_pred_train = model.predict(X_train)
        y_pred_test = model.predict(X_test)

        r2_train = r2_score(y_train, y_pred_train)
        r2_test = r2_score(y_test, y_pred_test)
        rmse_test = np.sqrt(mean_squared_error(y_test, y_pred_test))

        fold_r2_train.append(r2_train)
        fold_r2_test.append(r2_test)
        fold_rmse_test.append(rmse_test)
        fold_coefs.append(model.coef_)

        if (fold_i + 1) % 10 == 0:
            print(f"    Fold {fold_i+1}/45 done (test R²={r2_test:.4f})")

    fold_r2_train = np.array(fold_r2_train)
    fold_r2_test = np.array(fold_r2_test)
    fold_rmse_test = np.array(fold_rmse_test)
    coef_matrix = np.array(fold_coefs)  # (n_folds, n_features)

    print(f"\n  Train R²: {fold_r2_train.mean():.4f} ± {fold_r2_train.std():.4f}")
    print(f"  Test  R²: {fold_r2_test.mean():.4f} ± {fold_r2_test.std():.4f}")
    print(f"  Test RMSE: {fold_rmse_test.mean():.1f} ± {fold_rmse_test.std():.1f} kJ/mol")
    print(f"  Fixed alpha: {FIXED_ALPHA}, l1_ratio: {FIXED_L1}")

    # ------------------------------------------------------------------
    # 4) Selection frequency: how often each feature is non-zero
    # ------------------------------------------------------------------
    selection_freq = (np.abs(coef_matrix) > 1e-10).mean(axis=0)  # fraction of folds
    mean_coef = coef_matrix.mean(axis=0)
    median_coef = np.median(coef_matrix, axis=0)

    # Features selected in >50% of folds = robust
    robust_mask = selection_freq > 0.5
    n_robust = robust_mask.sum()
    print(f"\n  Features selected in >50% of folds: {n_robust}/{len(all_pairs)}")

    # ------------------------------------------------------------------
    # 5) Final model: train on ALL data
    # ------------------------------------------------------------------
    print("\nFitting final model on all replicas...")
    X_all = np.vstack([rep_X[r] for r in reps])
    y_all = np.concatenate([rep_y[r] for r in reps])

    model_final = ElasticNetCV(
        alphas=np.logspace(-1, 3, 15),
        l1_ratio=[0.5, 0.7, 0.9, 0.95],
        cv=5, fit_intercept=False, max_iter=200000,
        random_state=42, n_jobs=-1
    )
    model_final.fit(X_all, y_all)

    y_pred_all = model_final.predict(X_all)
    r2_all = r2_score(y_all, y_pred_all)
    beta_final = pd.Series(model_final.coef_, index=all_pairs)
    n_nz = (beta_final.abs() > 1e-10).sum()

    print(f"  alpha={model_final.alpha_:.4f}, l1_ratio={model_final.l1_ratio_:.2f}")
    print(f"  Full-data R²={r2_all:.4f}")
    print(f"  Non-zero β: {n_nz}/{len(all_pairs)}")

    # ------------------------------------------------------------------
    # 6) Compute NF using only robust features
    # ------------------------------------------------------------------
    print("\nComputing NetFavorability (robust features only)...")
    # Mean frequency across all replicas
    F_all_frames = np.vstack([
        replica_data[r][1].reindex(columns=all_pairs, fill_value=0).loc[1:FRAME_END].values
        for r in reps
    ])
    mean_freq = pd.Series(F_all_frames.mean(axis=0), index=all_pairs)
    bxf = beta_final * mean_freq

    # NF per antibody residue (robust only)
    residue_contrib = {}
    for i, pair_str in enumerate(all_pairs):
        if not robust_mask[i]:
            continue
        bxf_val = bxf[pair_str]
        if pd.isna(bxf_val) or abs(bxf_val) < 1e-10:
            continue
        parts = pair_str.split('-')
        gro1, gro2 = int(parts[0]), int(parts[1])
        ab_gro = gro2 if is_antibody(gro2) else (gro1 if is_antibody(gro1) else None)
        if ab_gro is None:
            continue
        residue_contrib.setdefault(ab_gro, []).append((bxf_val, beta_final[pair_str], selection_freq[i]))

    nf_rows = []
    for gro_resid, contribs in residue_contrib.items():
        nf = sum(abs(c[0]) for c in contribs)
        dom_idx = np.argmax([abs(c[0]) for c in contribs])
        dom_beta = contribs[dom_idx][1]
        dom_sign = "favorable" if dom_beta < 0 else "unfavorable"
        avg_sel = np.mean([c[2] for c in contribs])
        nf_rows.append({
            'GRO_resid': gro_resid,
            'Global_resid': gro_to_global(gro_resid),
            'Label': get_label(gro_resid, resnames),
            'NetFavorability': nf,
            'DominantSign': dom_sign,
            'DominantBeta': dom_beta,
            'NumPairs': len(contribs),
            'AvgSelectionFreq': avg_sel,
        })

    nf_df = pd.DataFrame(nf_rows).sort_values('NetFavorability', ascending=False).reset_index(drop=True)
    p70 = nf_df['NetFavorability'].quantile(0.70) if len(nf_df) > 0 else 0
    nf_df['Classification'] = nf_df['NetFavorability'].apply(
        lambda x: "hot" if x > p70 else "warm")

    nf_out = os.path.join(BASE_DIR, "net_favorability_replicaCV.csv")
    nf_df.to_csv(nf_out, index=False)
    print(f"Saved: {nf_out}")

    print("\nRobust residues (selected >50% of folds):")
    for _, row in nf_df.iterrows():
        print(f"  {row['Label']:>6s} (Global {row['Global_resid']:>3d}): "
              f"NF={row['NetFavorability']:.2f} [{row['DominantSign']}] "
              f"({row['NumPairs']} pairs, sel={row['AvgSelectionFreq']:.0%}) - {row['Classification']}")

    # ------------------------------------------------------------------
    # 7) Save coefficients
    # ------------------------------------------------------------------
    coeff_df = pd.DataFrame({
        'ResiduePair': all_pairs,
        'Beta_final': beta_final.values,
        'MeanFreq': mean_freq.values,
        'Beta_x_MeanFreq': bxf.values,
        'SelectionFreq': selection_freq,
        'Robust': robust_mask,
    })
    coeff_out = os.path.join(BASE_DIR, "elastic_net_coefficients_replicaCV.csv")
    coeff_df.to_csv(coeff_out, index=False)
    print(f"Saved: {coeff_out}")

    # ------------------------------------------------------------------
    # 8) Save CV fold stats for figure generation
    # ------------------------------------------------------------------
    cv_stats = pd.DataFrame({
        'fold': range(n_folds),
        'r2_train': fold_r2_train,
        'r2_test': fold_r2_test,
        'rmse_test': fold_rmse_test,
    })
    cv_stats.to_csv(os.path.join(BASE_DIR, "replicaCV_fold_stats.csv"), index=False)
    np.save(os.path.join(BASE_DIR, "replicaCV_selection_freq.npy"), selection_freq)
    print(f"Saved CV fold stats and selection frequencies")

    if not MAKE_FIGS:
        print("\nSkipping figures (--no-figs). Run locally to generate plots.")
        print("\n" + "=" * 70)
        print("ALL DONE")
        print("=" * 70)
        return

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    # ------------------------------------------------------------------
    # 8b) Figure: CV summary
    # ------------------------------------------------------------------
    print("\nGenerating figures...")
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    # A: Train vs Test R²
    ax = axes[0]
    ax.scatter(fold_r2_train, fold_r2_test, s=30, alpha=0.6, c='#1565C0', edgecolors='#333', lw=0.5)
    ax.plot([0, 1], [0, 1], 'k--', lw=0.8, alpha=0.4)
    ax.set_xlabel(r'Train $R^2$')
    ax.set_ylabel(r'Test $R^2$')
    ax.set_title('A', fontweight='bold', loc='left', fontsize=14)
    stats = (f'Train R² = {fold_r2_train.mean():.4f} ± {fold_r2_train.std():.4f}\n'
             f'Test R² = {fold_r2_test.mean():.4f} ± {fold_r2_test.std():.4f}\n'
             f'Test RMSE = {fold_rmse_test.mean():.1f} ± {fold_rmse_test.std():.1f} kJ/mol\n'
             f'{n_folds} folds (leave-2-replicas-out)')
    ax.text(0.05, 0.95, stats, transform=ax.transAxes, fontsize=9, va='top',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#ccc', alpha=0.9))
    ax.set_xlim(0.8, 1.0)
    ax.set_ylim(0.0, 1.0)

    # B: Selection frequency histogram
    ax = axes[1]
    sel_nonzero = selection_freq[selection_freq > 0]
    ax.hist(sel_nonzero, bins=20, color='#1565C0', edgecolor='white', lw=0.5, alpha=0.8)
    ax.axvline(0.5, color='#C62828', ls='--', lw=1.5, label='>50% threshold')
    ax.set_xlabel('Selection frequency (fraction of 45 folds)')
    ax.set_ylabel('Number of features')
    ax.set_title('B', fontweight='bold', loc='left', fontsize=14)
    ax.legend(fontsize=9)
    ax.text(0.95, 0.95, f'{n_robust} robust features\n(selected >50% of folds)',
            transform=ax.transAxes, fontsize=9, va='top', ha='right',
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#ccc', alpha=0.9))

    # C: NF bar chart (robust only)
    ax = axes[2]
    if len(nf_df) > 0:
        nf_plot = nf_df.sort_values('NetFavorability', ascending=True)
        colors = []
        for _, r in nf_plot.iterrows():
            if r['Classification'] == 'hot':
                colors.append('#00838F')
            elif r['DominantSign'] == 'unfavorable':
                colors.append('#E65100')
            else:
                colors.append('#455A64')
        y_pos = np.arange(len(nf_plot))
        ax.barh(y_pos, nf_plot['NetFavorability'].values, color=colors,
                edgecolor='white', lw=0.3, height=0.65)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(nf_plot['Label'].values, fontsize=9, fontfamily='monospace')
        ax.set_xlabel(r'Net Favorability $\left(\sum |\beta_k \times \bar{F}_k|\right)$')
    ax.set_title('C', fontweight='bold', loc='left', fontsize=14)
    ax.legend(handles=[Patch(fc='#00838F', label='Hot spot'),
                       Patch(fc='#455A64', label='Warm (favorable)'),
                       Patch(fc='#E65100', label='Warm (unfavorable)')],
              loc='lower right', fontsize=9, framealpha=0.9)

    for a in axes:
        a.spines['top'].set_visible(False)
        a.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'fig_replicaCV_regression.png'), dpi=300)
    plt.savefig(os.path.join(FIG_DIR, 'fig_replicaCV_regression.pdf'))
    print("Saved fig_replicaCV_regression.png/pdf")
    plt.close()

    print("\n" + "=" * 70)
    print("ALL DONE")
    print("=" * 70)


if __name__ == "__main__":
    main()
