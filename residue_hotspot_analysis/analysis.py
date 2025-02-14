import numpy as np
import pandas as pd
import MDAnalysis as mda
import mdtraj as md
from tqdm import tqdm
import os
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge

# --------------------Class for MD Interaction Analysis (Within 6 angstrom cutoff, computes interaction freequencies within the cutoff.) --------------------

class InteractionAnalyzer:
    """
    Extracts interaction frequencies between surface residues from an MD trajectory.
    """
    def __init__(self, topology_file, traj_file, group1_resids, group2_resids, distance_threshold=6.0):
        self.topology_file = topology_file
        self.traj_file = traj_file
        self.group1_resids = group1_resids
        self.group2_resids = group2_resids
        self.distance_threshold = distance_threshold
        self.output_file = "output_data/interaction_frequencies.csv"
        os.makedirs("output_data", exist_ok=True)

    def calculate_interaction_frequencies(self):
        """
        Compute interaction frequencies between two residue groups within a given distance threshold.
        """
        print("\n[INFO] Starting MD trajectory processing...")
        universe = mda.Universe(self.topology_file, self.traj_file)
        group1 = universe.select_atoms(f"resid {self.group1_resids[0]}:{self.group1_resids[1]} and protein")
        group2 = universe.select_atoms(f"resid {self.group2_resids[0]}:{self.group2_resids[1]} and protein")

        interaction_data = []
        
        for ts in tqdm(universe.trajectory, desc="Processing Frames"):
            distances = np.linalg.norm(group1.positions[:, np.newaxis] - group2.positions, axis=2)
            interactions = distances < self.distance_threshold

            for i, residue1 in enumerate(group1.residues):
                for j, residue2 in enumerate(group2.residues):
                    interaction_freq = np.sum(interactions[group1.resids == residue1.resid, :][:, group2.resids == residue2.resid])
                    if interaction_freq > 0:
                        interaction_data.append([ts.frame, f"{residue1.resid}-{residue2.resid}", interaction_freq])

        df_interactions = pd.DataFrame(interaction_data, columns=["Frame", "ResiduePair", "InteractionFrequency"])
        df_interactions.to_csv(self.output_file, index=False)
        print(f"[INFO] Interaction frequencies saved to {self.output_file}")
        return self.output_file

# -------------------- Class for Ridge Regression Analysis (Get regression coefficients/ compares the plot of change in interaction energy with change in frequency - inversly proprtional in general) --------------------

class RegressionAnalyzer:
    """
    Performs Ridge regression on interaction frequency data.
    """
    def __init__(self, interaction_file, energy_file, mutation_residues=None, interaction_threshold=70, alpha=10):
        self.interaction_file = interaction_file
        self.energy_file = energy_file
        self.mutation_residues = mutation_residues if mutation_residues else []
        self.interaction_threshold = interaction_threshold
        self.alpha = alpha
        self.output_csv = "ridge_regression_coefficients.csv"

    def run_analysis(self):
        """
        Run Ridge regression to identify important residue interactions affecting energy changes.
        """
        print("\n[INFO] Loading data for regression analysis...")
        df_energy = pd.read_csv(self.energy_file)
        df_frequency = pd.read_csv(self.interaction_file)

        df_residue_pairs = df_frequency.pivot_table(index='Frame', columns='ResiduePair', values='InteractionFrequency', fill_value=0)
        
        filtered_columns = [col for col in df_residue_pairs.columns if not (int(col.split('-')[0]) in self.mutation_residues or int(col.split('-')[1]) in self.mutation_residues)]
        df_residue_pairs_filtered = df_residue_pairs[filtered_columns]

        frequent_pairs = df_residue_pairs_filtered.ge(self.interaction_threshold).any()
        df_residue_pairs_filtered = df_residue_pairs_filtered.loc[:, frequent_pairs]

        common_frames = df_residue_pairs_filtered.index.intersection(df_energy.index)
        df_residue_pairs_aligned = df_residue_pairs_filtered.loc[common_frames]
        df_energy_aligned = df_energy.loc[common_frames, 'AverageInteractionEnergy']

        df_residue_pairs_diff = df_residue_pairs_aligned.diff().dropna()
        df_energy_change = df_energy_aligned.diff().dropna()

        common_frames_diff = df_residue_pairs_diff.index.intersection(df_energy_change.index)
        df_residue_pairs_diff = df_residue_pairs_diff.loc[common_frames_diff]
        df_energy_change = df_energy_change.loc[common_frames_diff]

        scaler = StandardScaler()
        residue_pairs_scaled = scaler.fit_transform(df_residue_pairs_diff)

        ridge = Ridge(alpha=self.alpha)
        ridge.fit(residue_pairs_scaled, df_energy_change)
        ridge_coef = ridge.coef_

        selected_features = np.where(ridge_coef != 0)[0]
        ranked_indices = np.argsort(np.abs(ridge_coef[selected_features]))[::-1]
        ranked_features = df_residue_pairs_diff.columns[selected_features][ranked_indices]
        ranked_coefficients = ridge_coef[selected_features][ranked_indices]

        df_ridge_coefficients = pd.DataFrame({"Residue Pair": ranked_features, "Coefficient": ranked_coefficients})
        df_ridge_coefficients.to_csv(self.output_csv, index=False)
        print(f"[INFO] Regression results saved to {self.output_csv}")

        self.generate_plots(df_energy_change, ridge.predict(residue_pairs_scaled))

    def generate_plots(self, actual_energy_change, predicted_energy_change):
        """
        Generate regression plots for analysis.
        """
        plt.figure(figsize=(8, 6))
        plt.scatter(actual_energy_change, predicted_energy_change, color='b', alpha=0.6, label="Predicted vs Actual")
        plt.plot(actual_energy_change, actual_energy_change, 'r--', label="Ideal Fit")
        plt.xlabel("Actual Change in Interaction Energy")
        plt.ylabel("Predicted Change in Interaction Energy")
        plt.title("Predicted vs Actual Energy Change (Ridge Regression)")
        plt.legend()
        plt.grid(True)
        plt.savefig("regression_plot.png")
        print("[INFO] Regression plot saved as 'regression_plot.png'")
        plt.show()


class HotspotAnalysisPipeline:
    """
    Manages the full workflow: MD interaction analysis + regression analysis.
    """
    def __init__(self, topology, trajectory, energy_file, group1_resids, group2_resids, mutation_residues):
        self.topology = topology
        self.trajectory = trajectory
        self.energy_file = energy_file
        self.group1_resids = group1_resids
        self.group2_resids = group2_resids
        self.mutation_residues = mutation_residues

    def run(self):
        """
        Run the full pipeline: compute interactions & perform regression.
        """
        print("\n### Step 1: Computing Interaction Frequencies ###")
        interaction_analyzer = InteractionAnalyzer(self.topology, self.trajectory, self.group1_resids, self.group2_resids)
        interaction_file = interaction_analyzer.calculate_interaction_frequencies()

        print("\n### Step 2: Running Regression Analysis ###")
        regression_analyzer = RegressionAnalyzer(interaction_file, self.energy_file, self.mutation_residues)
        regression_analyzer.run_analysis()


if __name__ == "__main__":
    pipeline = HotspotAnalysisPipeline(
        topology="your_topology.gro",
        trajectory="your_trajectory.xtc",
        energy_file="interaction_energy.csv",
        group1_resids=[1, 195],
        group2_resids=[196, 422],
        mutation_residues=[7, 39, 41, 43, 85]
    )
    pipeline.run()
