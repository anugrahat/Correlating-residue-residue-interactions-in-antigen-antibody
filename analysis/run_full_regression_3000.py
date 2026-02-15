#!/usr/bin/env python3
"""
Full elastic net regression on C1 WT pulling data — first 3000 frames only.
Generates all figures:
  fig2: NetFavorability bar chart + β heatmap
  fig3: Regression diagnostics (pred vs actual, residuals)
  fig4: Temporal contributions (β×ΔF(t), β×F(t)) for top 10 pairs
Saves CSVs: coefficients, net favorability.
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from sklearn.linear_model import ElasticNetCV
from sklearn.metrics import r2_score
from scipy.stats import pearsonr

# ==============================================================================
# Config
# ==============================================================================
BASE_DIR = "/home/anugraha/c1_WT/analysis"
FIG_DIR = "/home/anugraha/antibody_optimization/figures"
FREQ_FILE = os.path.join(BASE_DIR, "average_frequency.csv")
ENERGY_FILE = os.path.join(BASE_DIR, "interaction_energy.csv")
GRO_PATH = "/home/anugraha/c1_WT/pull/md10_protein.gro"

FRAME_END = 1500
REF_FRAME = 0

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

def pair_label(pair_str, resnames):
    parts = pair_str.split('-')
    labels = []
    for p in parts:
        gro = int(p)
        labels.append(get_label(gro, resnames))
    return f"{labels[0]}–{labels[1]}"

def is_antibody(gro_resid):
    return 195 <= gro_resid <= 421

# ==============================================================================
# Main
# ==============================================================================
def main():
    print("=" * 70)
    print(f"ELASTIC NET REGRESSION — C1 WT, FRAME_END={FRAME_END}")
    print("=" * 70)

    resnames = load_gro_resnames(GRO_PATH)
    print(f"Loaded {len(resnames)} residue names from GRO.")

    # ------------------------------------------------------------------
    # 1) Load data
    # ------------------------------------------------------------------
    print("\nLoading data...")
    df_energy = pd.read_csv(ENERGY_FILE)
    df_frequency = pd.read_csv(FREQ_FILE)
    print(f"  Energy: {df_energy.shape[0]} frames")
    print(f"  Frequency: {df_frequency['ResiduePair'].nunique()} unique pairs")

    # Pivot
    # Handle both column names (MeanFrequency from average_contact_freq.py, InteractionFrequency from average_data.py)
    freq_col = 'MeanFrequency' if 'MeanFrequency' in df_frequency.columns else 'InteractionFrequency'
    df_pivot = df_frequency.pivot_table(
        index='Frame', columns='ResiduePair',
        values=freq_col, fill_value=0
    ).sort_index()

    # Filter out intra-antibody pairs (both residues >= 195 in GRO numbering)
    def is_inter_contact(pair_str):
        g1, g2 = [int(x) for x in pair_str.split('-')]
        return not (g1 >= 195 and g2 >= 195)

    inter_cols = [c for c in df_pivot.columns if is_inter_contact(c)]
    n_removed = len(df_pivot.columns) - len(inter_cols)
    df_pivot = df_pivot[inter_cols]
    print(f"  Removed {n_removed} intra-antibody pairs, {len(inter_cols)} inter-contact pairs remain")

    # Align
    df_energy_idx = df_energy.set_index('Frame').sort_index()
    common = df_pivot.index.intersection(df_energy_idx.index)
    df_E = df_energy_idx.loc[common, 'AverageInteractionEnergy']
    df_F = df_pivot.loc[common]

    # Delta from ref frame
    E_ref = df_E.loc[REF_FRAME]
    F_ref = df_F.loc[REF_FRAME]
    dE = (df_E - E_ref).loc[:FRAME_END]
    dF = df_F.subtract(F_ref, axis='columns').loc[:FRAME_END]

    # Drop ref frame
    dE = dE.drop(index=REF_FRAME)
    dF = dF.drop(index=REF_FRAME)

    X = dF.values
    y = dE.values
    pair_names = dF.columns.tolist()
    print(f"  Regression: {len(y)} frames, {X.shape[1]} features (frames 1–{FRAME_END})")

    # ------------------------------------------------------------------
    # 2) Fit ElasticNetCV
    # ------------------------------------------------------------------
    print("\nFitting ElasticNetCV (n_jobs=-1)...")
    model = ElasticNetCV(
        alphas=np.logspace(-2, 2, 10),
        l1_ratio=[0.3, 0.5, 0.7, 0.9],
        cv=5, fit_intercept=False, max_iter=200000,
        random_state=42, n_jobs=-1
    )
    model.fit(X, y)

    y_pred = model.predict(X)
    r2 = r2_score(y, y_pred)
    pr, _ = pearsonr(y, y_pred)
    beta = pd.Series(model.coef_, index=pair_names, name='Beta')
    n_nz = (beta != 0).sum()

    print(f"  alpha={model.alpha_:.4f}, l1_ratio={model.l1_ratio_:.2f}")
    print(f"  R²={r2:.4f}, Pearson r={pr:.4f}")
    print(f"  Non-zero β: {n_nz}/{len(beta)}")

    # ------------------------------------------------------------------
    # 3) β × mean frequency
    # ------------------------------------------------------------------
    df_F_reg = df_F.loc[1:FRAME_END]
    mean_freq = df_F_reg.mean(axis=0)
    bxf = beta * mean_freq

    # ------------------------------------------------------------------
    # 4) Save coefficients CSV
    # ------------------------------------------------------------------
    coeff_df = pd.DataFrame({
        'ResiduePair': pair_names,
        'Beta': beta.values,
        'MeanFreq': mean_freq.values,
        'Beta_x_MeanFreq': bxf.values,
    })
    coeff_out = os.path.join(BASE_DIR, "elastic_net_coefficients_3000.csv")
    coeff_df.to_csv(coeff_out, index=False)
    print(f"\nSaved: {coeff_out}")

    # ------------------------------------------------------------------
    # 5) NetFavorability per antibody residue
    # ------------------------------------------------------------------
    residue_contrib = {}
    for pair_str, bxf_val in bxf.items():
        if pd.isna(bxf_val) or bxf_val == 0.0:
            continue
        parts = pair_str.split('-')
        gro1, gro2 = int(parts[0]), int(parts[1])
        ab_gro = gro2 if is_antibody(gro2) else (gro1 if is_antibody(gro1) else None)
        if ab_gro is None:
            continue
        residue_contrib.setdefault(ab_gro, []).append((bxf_val, beta[pair_str]))

    nf_rows = []
    for gro_resid, contribs in residue_contrib.items():
        nf = sum(abs(c[0]) for c in contribs)
        dom_idx = np.argmax([abs(c[0]) for c in contribs])
        dom_beta = contribs[dom_idx][1]
        dom_sign = "favorable" if dom_beta < 0 else "unfavorable"
        nf_rows.append({
            'GRO_resid': gro_resid,
            'Global_resid': gro_to_global(gro_resid),
            'Label': get_label(gro_resid, resnames),
            'NetFavorability': nf,
            'DominantSign': dom_sign,
            'DominantBeta': dom_beta,
            'NumPairs': len(contribs),
        })

    nf_df = pd.DataFrame(nf_rows).sort_values('NetFavorability', ascending=False).reset_index(drop=True)
    p70 = nf_df['NetFavorability'].quantile(0.70)
    p20 = nf_df['NetFavorability'].quantile(0.20)
    nf_df['Classification'] = nf_df['NetFavorability'].apply(
        lambda x: "hot" if x > p70 else ("warm" if x > p20 else "cold"))

    nf_out = os.path.join(BASE_DIR, "net_favorability_3000.csv")
    nf_df.to_csv(nf_out, index=False)
    print(f"Saved: {nf_out}")

    print("\nTop residues:")
    for _, row in nf_df.head(15).iterrows():
        print(f"  {row['Label']:>6s} (Global {row['Global_resid']:>3d}): "
              f"NF={row['NetFavorability']:.2f} [{row['DominantSign']}] "
              f"({row['NumPairs']} pairs) - {row['Classification']}")

    # ==================================================================
    # FIGURE 2: NetFavorability bar + β heatmap
    # ==================================================================
    print("\n\nGenerating Fig 2 (bar + heatmap)...")

    top20 = nf_df.head(20).copy().iloc[::-1]
    bar_colors = ['#2166AC' if s == 'favorable' else '#B2182B' for s in top20['DominantSign']]

    # Top 30 pairs for heatmap
    abs_beta = beta.abs()
    top30_pairs = abs_beta[abs_beta > 0].sort_values(ascending=False).head(30).index.tolist()

    heatmap_data = {}
    for pair_str in top30_pairs:
        parts = pair_str.split('-')
        gro_ant, gro_ab = int(parts[0]), int(parts[1])
        ant_label = get_label(gro_ant, resnames)
        ab_label = get_label(gro_ab, resnames)
        heatmap_data.setdefault(ant_label, {})[ab_label] = beta[pair_str]

    hm_df = pd.DataFrame(heatmap_data).T.fillna(0)
    # Sort columns by global resid
    ab_order = sorted(hm_df.columns, key=lambda x: int(x[1:]))
    ant_order = sorted(hm_df.index, key=lambda x: int(x[1:]))
    hm_df = hm_df.reindex(columns=ab_order, index=ant_order, fill_value=0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7),
                                     gridspec_kw={'width_ratios': [1, 1.5]})

    # Panel A: bar chart
    ax1.barh(range(len(top20)), top20['NetFavorability'].values, color=bar_colors)
    ax1.set_yticks(range(len(top20)))
    ax1.set_yticklabels(top20['Label'].values, fontsize=9)
    ax1.set_xlabel('NetFavorability (|β × MeanFreq|)', fontsize=11)
    ax1.set_title(f'A. NetFavorability (frames 1–{FRAME_END})', fontsize=13, fontweight='bold')

    # Annotate hot/warm spots
    hot_globals = {404, 405, 406}
    warm_globals = {287, 265, 262, 246, 248}
    for i, (_, row) in enumerate(top20.iterrows()):
        g = row['Global_resid']
        if g in hot_globals:
            ax1.text(row['NetFavorability'] + 0.5, i, 'hot spot', fontsize=8, color='red', va='center')
        elif g in warm_globals:
            ax1.text(row['NetFavorability'] + 0.5, i, 'warm spot', fontsize=8, color='orange', va='center')

    ax1.legend(handles=[
        plt.Rectangle((0, 0), 1, 1, fc='#2166AC', label='Favorable (β<0, stabilizing)'),
        plt.Rectangle((0, 0), 1, 1, fc='#B2182B', label='Unfavorable (β>0, destabilizing)'),
    ], fontsize=8, loc='lower right')

    # Panel B: heatmap
    vals = hm_df.values
    vmax = max(abs(vals.min()), abs(vals.max())) or 1
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    im = ax2.imshow(vals, cmap='RdBu_r', norm=norm, aspect='auto')
    ax2.set_xticks(range(len(hm_df.columns)))
    ax2.set_xticklabels(hm_df.columns, rotation=45, ha='right', fontsize=8)
    ax2.set_yticks(range(len(hm_df.index)))
    ax2.set_yticklabels(hm_df.index, fontsize=8)
    ax2.set_xlabel('Antibody Residue (Global)', fontsize=11)
    ax2.set_ylabel('Antigen Residue (Global)', fontsize=11)
    ax2.set_title(f'B. β coefficients (top 30 pairs)', fontsize=13, fontweight='bold')

    # Annotate non-zero cells
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            if vals[i, j] != 0:
                ax2.text(j, i, f"{vals[i, j]:.2f}", ha='center', va='center', fontsize=7,
                         color='white' if abs(vals[i, j]) > vmax * 0.6 else 'black')

    plt.colorbar(im, ax=ax2, label='β coefficient', shrink=0.8)
    plt.tight_layout()
    fig2_path = os.path.join(FIG_DIR, "fig2_elastic_net_heatmap_3000.png")
    fig.savefig(fig2_path, dpi=200, bbox_inches='tight')
    print(f"Saved: {fig2_path}")
    plt.close()

    # ==================================================================
    # FIGURE 3: Regression diagnostics
    # ==================================================================
    print("Generating Fig 3 (regression diagnostics)...")
    residuals = y - y_pred

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # A: Predicted vs Actual
    ax = axes[0]
    ax.scatter(y, y_pred, alpha=0.15, s=8, c='#2166AC', edgecolors='none')
    lims = [min(y.min(), y_pred.min()), max(y.max(), y_pred.max())]
    ax.plot(lims, lims, 'k--', lw=1, alpha=0.7)
    ax.set_xlabel('Actual ΔE (kcal/mol)', fontsize=12)
    ax.set_ylabel('Predicted ΔE (kcal/mol)', fontsize=12)
    ax.set_title('A. Predicted vs Actual', fontsize=13, fontweight='bold')
    stats = (f"R² = {r2:.4f}\nPearson r = {pr:.4f}\n"
             f"α = {model.alpha_:.2f}\nl₁ ratio = {model.l1_ratio_:.2f}\n"
             f"Non-zero β: {n_nz}/{len(beta)}\nN frames = {len(y)}")
    ax.text(0.05, 0.95, stats, transform=ax.transAxes, fontsize=10,
            va='top', bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.8))
    ax.set_aspect('equal', adjustable='box')

    # B: Residuals
    ax = axes[1]
    ax.scatter(y, residuals, alpha=0.15, s=8, c='#B2182B', edgecolors='none')
    ax.axhline(0, color='k', ls='--', lw=1, alpha=0.7)
    ax.set_xlabel('Actual ΔE (kcal/mol)', fontsize=12)
    ax.set_ylabel('Residual (kcal/mol)', fontsize=12)
    ax.set_title('B. Residuals', fontsize=13, fontweight='bold')
    rmse = np.sqrt(np.mean(residuals**2))
    mae = np.mean(np.abs(residuals))
    ax.text(0.05, 0.95, f"RMSE = {rmse:.2f}\nMAE = {mae:.2f}",
            transform=ax.transAxes, fontsize=10, va='top',
            bbox=dict(boxstyle='round,pad=0.5', fc='lightyellow', alpha=0.8))

    # C: Residual distribution
    ax = axes[2]
    ax.hist(residuals, bins=60, color='#2166AC', alpha=0.7, edgecolor='white', lw=0.5)
    ax.axvline(0, color='k', ls='--', lw=1, alpha=0.7)
    ax.set_xlabel('Residual (kcal/mol)', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('C. Residual Distribution', fontsize=13, fontweight='bold')
    ax.text(0.05, 0.95, f"μ = {np.mean(residuals):.3f}\nσ = {np.std(residuals):.2f}",
            transform=ax.transAxes, fontsize=10, va='top',
            bbox=dict(boxstyle='round,pad=0.5', fc='lightcyan', alpha=0.8))

    plt.tight_layout()
    fig3_path = os.path.join(FIG_DIR, "fig3_regression_diagnostics_3000.png")
    fig.savefig(fig3_path, dpi=200, bbox_inches='tight')
    print(f"Saved: {fig3_path}")
    plt.close()

    # ==================================================================
    # FIGURE 4: Temporal contributions for top 10 pairs
    # ==================================================================
    print("Generating Fig 4 (temporal contributions)...")

    # Top 10 by |β × MeanFreq|
    coeff_nz = coeff_df[coeff_df['Beta'] != 0].copy()
    coeff_nz['AbsBxF'] = coeff_nz['Beta_x_MeanFreq'].abs()
    top10 = coeff_nz.nlargest(10, 'AbsBxF')
    top10_pairs = top10['ResiduePair'].tolist()
    top10_betas = dict(zip(top10['ResiduePair'], top10['Beta']))

    frames = df_F_reg.index.values
    colors = plt.cm.tab10(np.linspace(0, 1, 10))

    fig, axes = plt.subplots(3, 1, figsize=(14, 14))

    # A: β × ΔF(t) time series
    ax = axes[0]
    for i, pair in enumerate(top10_pairs):
        b = top10_betas[pair]
        if pair in df_F.columns:
            delta_f = df_F.loc[1:FRAME_END, pair].values - F_ref[pair]
            contrib = b * delta_f
            ax.plot(df_F.loc[1:FRAME_END].index, contrib,
                    label=pair_label(pair, resnames), color=colors[i], alpha=0.8, lw=1.2)
    ax.axhline(0, color='k', ls='--', lw=0.8, alpha=0.5)
    ax.set_xlabel('Frame index', fontsize=12)
    ax.set_ylabel('β × ΔF(t)  (contribution to ΔE)', fontsize=12)
    ax.set_title(f'A. Per-pair contribution to predicted ΔE (frames 1–{FRAME_END})',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=9, ncol=2, loc='upper left', framealpha=0.9)
    ax.grid(alpha=0.3)

    # B: β × Mean Frequency (bar chart)
    ax = axes[1]
    top10_labels = [pair_label(p, resnames) for p in top10_pairs]
    top10_bxf = [top10.loc[top10['ResiduePair'] == p, 'Beta_x_MeanFreq'].values[0] for p in top10_pairs]
    bar_colors_b = ['#2166AC' if v < 0 else '#B2182B' for v in top10_bxf]
    bars = ax.barh(range(len(top10_labels)), top10_bxf, color=bar_colors_b)
    ax.set_yticks(range(len(top10_labels)))
    ax.set_yticklabels(top10_labels, fontsize=9)
    ax.axvline(0, color='k', ls='--', lw=0.8, alpha=0.5)
    ax.set_xlabel('β × Mean Frequency', fontsize=12)
    ax.set_title(f'B. β × Mean Frequency (top 10 pairs, frames 1–{FRAME_END})',
                 fontsize=13, fontweight='bold')
    for i, (bar, val) in enumerate(zip(bars, top10_bxf)):
        ax.text(val - 1 if val < 0 else val + 0.5, i, f'{val:.1f}', va='center', fontsize=8)
    ax.invert_yaxis()
    ax.grid(alpha=0.3, axis='x')

    # C: β × ΔF(t) time series (same pairs, delta frequency from ref frame)
    ax = axes[2]
    for i, pair in enumerate(top10_pairs):
        b = top10_betas[pair]
        if pair in dF.columns:
            bxdf_t = b * dF[pair].values
            ax.plot(dF.index, bxdf_t,
                    label=pair_label(pair, resnames), color=colors[i], alpha=0.8, lw=1.2)
    ax.axhline(0, color='k', ls='--', lw=0.8, alpha=0.5)
    ax.set_xlabel('Frame index', fontsize=12)
    ax.set_ylabel('β × ΔF(t)', fontsize=12)
    ax.set_title(f'C. β × ΔFrequency over frames (top 10 pairs, frames 1–{FRAME_END})',
                 fontsize=13, fontweight='bold')
    ax.legend(fontsize=9, ncol=2, loc='upper left', framealpha=0.9)
    ax.grid(alpha=0.3)

    plt.tight_layout()
    fig4_path = os.path.join(FIG_DIR, "fig4_temporal_contributions_3000.png")
    fig.savefig(fig4_path, dpi=200, bbox_inches='tight')
    print(f"Saved: {fig4_path}")
    plt.close()

    print("\n" + "=" * 70)
    print("ALL DONE")
    print("=" * 70)

if __name__ == "__main__":
    main()
