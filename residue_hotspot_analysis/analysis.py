import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge, LinearRegression
import matplotlib.pyplot as plt

import seaborn as sns

df_energy = pd.read_csv('interaction_energy.csv')
df_frequency = pd.read_csv('average_frequency.csv')

mutation_residues = []

#mutation_residues = [7, 39, 41, 43, 85, 108, 114, 145, 146, 152, 161, 164, 166, 169, 173]

df_residue_pairs = df_frequency.pivot_table(index='Frame', columns='ResiduePair', values='InteractionFrequency', fill_value=0)

def exclude_mutation_pairs(df_residue_pairs, mutation_residues):
    filtered_columns = [col for col in df_residue_pairs.columns 
                        if not (int(col.split('-')[0]) in mutation_residues or int(col.split('-')[1]) in mutation_residues)]
    
    return df_residue_pairs[filtered_columns]

df_residue_pairs_filtered = exclude_mutation_pairs(df_residue_pairs, mutation_residues)

INTERACTION_THRESHOLD = 70
frequent_pairs = df_residue_pairs_filtered.ge(INTERACTION_THRESHOLD).any()
df_residue_pairs_filtered = df_residue_pairs_filtered.loc[:, frequent_pairs]

common_frames = df_residue_pairs_filtered.index.intersection(df_energy.index)
df_residue_pairs_aligned = df_residue_pairs_filtered.loc[common_frames]
df_energy_aligned = df_energy.loc[common_frames, 'AverageInteractionEnergy']

fixed_window_size = 100
df_residue_pairs_smoothed = df_residue_pairs_aligned.rolling(window=fixed_window_size, min_periods=fixed_window_size).mean()
df_energy_smoothed = df_energy_aligned.rolling(window=fixed_window_size, min_periods=fixed_window_size).mean()

df_residue_pairs_smoothed = df_residue_pairs_smoothed.loc[1010:1300]
df_energy_smoothed = df_energy_smoothed.loc[1010:1300]

df_residue_pairs_diff = df_residue_pairs_smoothed.diff().dropna()
df_energy_change = df_energy_smoothed.diff().dropna()

common_frames_diff = df_residue_pairs_diff.index.intersection(df_energy_change.index)
df_residue_pairs_diff = df_residue_pairs_diff.loc[common_frames_diff]
df_energy_change = df_energy_change.loc[common_frames_diff]

scaler = StandardScaler()
residue_pairs_scaled = scaler.fit_transform(df_residue_pairs_diff)

ridge = Ridge(alpha=10)  # The alpha parameter controls the regularization strength
ridge.fit(residue_pairs_scaled, df_energy_change)
ridge_coef = ridge.coef_

selected_features = np.where(ridge_coef != 0)[0]

ridge_coef_selected = ridge_coef[selected_features]  
ranked_indices = np.argsort(np.abs(ridge_coef_selected))[::-1]  
ranked_features = df_residue_pairs_diff.columns[selected_features][ranked_indices]  
ranked_coefficients = ridge_coef_selected[ranked_indices]  


print("Ranked Ridge Coefficients for Selected Features (Residue Pairs):")
for i, feature in enumerate(ranked_features):
    print(f"Rank {i+1}: {feature} -> Coefficient: {ranked_coefficients[i]}")

ridge_predictions = ridge.predict(residue_pairs_scaled)

r_squared = ridge.score(residue_pairs_scaled, df_energy_change)
print(f"R-squared value for Ridge regression: {r_squared:.4f}")

plt.figure(figsize=(10, 6))
plt.scatter(df_energy_change, ridge_predictions, label='Predicted vs Actual', color='b')
plt.plot(df_energy_change, df_energy_change, color='r', linestyle='--', label=f'Ideal Fit (R² = {r_squared:.4f})')
plt.xlabel('Actual Change in Interaction Energy')
plt.ylabel('Predicted Change in Interaction Energy')
plt.title('Predicted vs Actual Interaction Energy Change (Ridge Regression)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

df_ridge_coefficients = pd.DataFrame({
    'Residue Pair': df_residue_pairs_diff.columns[selected_features],  # Residue pairs
    'Coefficient': ridge_coef[selected_features]  # Corresponding Ridge coefficients
})

df_ridge_coefficients_sorted = df_ridge_coefficients.reindex(np.argsort(np.abs(df_ridge_coefficients['Coefficient']))[::-1])

df_ridge_coefficients_sorted.to_csv('ridge_regression_coefficients.csv', index=False)

print("Ridge regression coefficients have been saved to 'ridge_regression_coefficients.csv'.")


top_10_positive = ranked_features[:10]
top_10_negative = ranked_features[-10:]

top_20_pairs = list(top_10_positive) + list(top_10_negative)
top_20_coefficients = np.concatenate((ranked_coefficients[:10], ranked_coefficients[-10:]))

print("Top 10 Positive and Top 10 Negative Residue Pairs with Coefficients:")
for i, (pair, coef) in enumerate(zip(top_20_pairs, top_20_coefficients)):
    sign = "Positive" if coef > 0 else "Negative"
    print(f"{sign} Coefficient Rank {i+1}: Residue Pair {pair} -> Coefficient: {coef:.4f}")

fig, axs = plt.subplots(5, 4, figsize=(20, 15))  # 5 rows, 4 columns for 20 residue pairs

for i, residue_pair in enumerate(top_20_pairs):
    row = i // 4  # Determine the row index for the subplot
    col = i % 4   # Determine the column index for the subplot

    
    if residue_pair in df_residue_pairs_diff.columns:
        interaction_frequencies = df_residue_pairs_diff[residue_pair]

   
        axs[row, col].plot(interaction_frequencies.index, interaction_frequencies, label=f"Frequency ({residue_pair})", color='b')
        axs[row, col].set_xlabel("Frames")
        axs[row, col].set_ylabel("Interaction Frequency")
        axs[row, col].set_title(f"Freq: Residue Pair {residue_pair} (Coef: {top_20_coefficients[i]:.4f})")
        axs[row, col].grid(True)
        axs[row, col].legend()


        ax2 = axs[row, col].twinx()
        ax2.plot(df_energy_change.index, df_energy_change, label="Energy Change", color='r', linestyle='--')
        ax2.set_ylabel("Energy Change")
        ax2.legend(loc='upper right')

plt.tight_layout()
plt.show()
fig, axs = plt.subplots(5, 4, figsize=(20, 15)) 


norm = plt.Normalize(df_energy_change.index.min(), df_energy_change.index.max())
cmap = plt.get_cmap('viridis')  

for i, residue_pair in enumerate(top_20_pairs):
    row = i // 4  
    col = i % 4   


    if residue_pair in df_residue_pairs_diff.columns:
        interaction_frequencies = df_residue_pairs_diff[residue_pair]


        sc = axs[row, col].scatter(interaction_frequencies, df_energy_change, 
                                   c=df_energy_change.index, cmap=cmap, norm=norm, s=20)
        axs[row, col].set_xlabel(f'Freq ({residue_pair})')
        axs[row, col].set_ylabel('Energy Change')
        axs[row, col].set_title(f"Freq vs Energy: {residue_pair} (Coef: {top_20_coefficients[i]:.4f})")
        axs[row, col].grid(True)


plt.subplots_adjust(bottom=0.15)


cbar_ax = fig.add_axes([0.15, 0.05, 0.7, 0.03])  # [left, bottom, width, height] for the color bar
cbar = fig.colorbar(sc, cax=cbar_ax, orientation='horizontal')
cbar.set_label('Frame Progression')

plt.tight_layout(rect=[0, 0.1, 1, 1]) 
plt.show()


def identify_peak_frames(df_residue_pairs_diff, top_20_pairs, threshold_factor=1):
    peak_frames = {}
    
    for pair in top_20_pairs:
        
        mean_change = df_residue_pairs_diff[pair].mean()
        std_change = df_residue_pairs_diff[pair].std()
        
       
        positive_threshold = mean_change + threshold_factor * std_change
        negative_threshold = mean_change - threshold_factor * std_change
        
        
        peak_frames_pos = df_residue_pairs_diff.index[df_residue_pairs_diff[pair] >= positive_threshold].tolist()
        peak_frames_neg = df_residue_pairs_diff.index[df_residue_pairs_diff[pair] <= negative_threshold].tolist()
        
        
        peak_frames[pair] = sorted(set(peak_frames_pos + peak_frames_neg))
    
    return peak_frames


peak_frames_dict = identify_peak_frames(df_residue_pairs_diff, top_20_pairs)

 

fig, axs = plt.subplots(5, 4, figsize=(20, 15))  # 5 rows, 4 columns for 20 residue pairs



predicted_energy_changes = ridge.predict(scaler.transform(df_residue_pairs_diff)) 
fig, axs = plt.subplots(5, 4, figsize=(20, 15)) 

# Plot each residue pair's interaction frequency, actual energy change, and predicted energy change
for i, residue_pair in enumerate(top_20_pairs):
    row = i // 4  # Determine the row index for the subplot
    col = i % 4   # Determine the column index for the subplot

    # Extract interaction frequencies for the residue pair
    if residue_pair in df_residue_pairs_diff.columns:
        interaction_frequencies = df_residue_pairs_diff[residue_pair]

        # Plot interaction frequency on the left y-axis
        axs[row, col].plot(interaction_frequencies.index, interaction_frequencies, label=f"Frequency ({residue_pair})", color='b')
        axs[row, col].set_xlabel("Frames")
        axs[row, col].set_ylabel("Interaction Frequency")
        axs[row, col].set_title(f"Freq: Residue Pair {residue_pair} (Coef: {top_20_coefficients[i]:.4f})")
        axs[row, col].grid(True)
        axs[row, col].legend()

        # Create twin axis for actual and predicted energy change (right y-axis)
        ax2 = axs[row, col].twinx()
        ax2.plot(df_energy_change.index, df_energy_change, label="Actual Energy Change", color='r', linestyle='--')
        ax2.plot(df_energy_change.index, predicted_energy_changes[df_energy_change.index], label="Predicted Energy Change", color='purple', linestyle=':')
        ax2.set_ylabel("Energy Change")
        ax2.legend(loc='upper right')

plt.tight_layout()
plt.show()


beta_times_deltaF = {}
for i, residue_pair in enumerate(top_20_pairs):
    delta_F = df_residue_pairs_diff[residue_pair]
    beta_times_deltaF[residue_pair] = top_20_coefficients[i] * delta_F


fig, axs = plt.subplots(5, 4, figsize=(20, 15))

for i, residue_pair in enumerate(top_20_pairs):
    row = i // 4
    col = i % 4
    beta_times_deltaF_values = beta_times_deltaF[residue_pair]
    axs[row, col].plot(df_residue_pairs_diff.index, beta_times_deltaF_values, color='green', label=f"Beta * Delta F ({residue_pair})")
    axs[row, col].set_xlabel("Frames")
    axs[row, col].set_ylabel("Beta * Delta F")
    axs[row, col].set_title(f"Beta * Delta F: {residue_pair}")
    axs[row, col].grid(True)
    axs[row, col].legend()

plt.tight_layout()
plt.savefig('beta_deltaF_vs_frames.png', dpi=300)  
plt.show()







