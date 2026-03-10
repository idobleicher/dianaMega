#!/usr/bin/env python3
"""
Create heatmaps showing amino acid enrichment/depletion in unstable vs stable peptides
Similar to Timms et al. 2019 Figure 2A,B
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

OUTPUT_DIR = "logo_results_motif_analysis"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("="*80)
print("PSI-BASED AMINO ACID ENRICHMENT HEATMAPS")
print("="*80)

# Load data
print("\nLoading data...")
screen = pd.read_excel('UBR3 Nt screen.xlsx')

# Find PSI column
psi_col = None
for col in screen.columns:
    if 'PSI' in col.upper() or 'SCORE' in col.upper():
        psi_col = col
        break

if not psi_col:
    print("ERROR: No PSI column found!")
    exit(1)

print(f"Using PSI column: {psi_col}")

# Filter out rows with missing PSI values
screen_with_psi = screen.dropna(subset=[psi_col])
print(f"Total sequences with PSI data: {len(screen_with_psi)}")

# Define PSI threshold
psi_threshold = 3.6

# Split into unstable and stable groups
unstable = screen_with_psi[screen_with_psi[psi_col] <= psi_threshold]
stable = screen_with_psi[screen_with_psi[psi_col] > psi_threshold]

print(f"\nUnstable peptides (PSI <= {psi_threshold}): {len(unstable)}")
print(f"Stable peptides (PSI > {psi_threshold}): {len(stable)}")

# Get sequences
unstable_sequences = unstable['AA_seq'].tolist()
stable_sequences = stable['AA_seq'].tolist()
all_sequences = screen_with_psi['AA_seq'].tolist()

# Determine sequence length to analyze (let's do 2-24 like the paper)
max_position = 24

# Define amino acids in order
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

print(f"\nAnalyzing positions 2-{max_position}")

# ============================================================
# Calculate frequencies for each group
# ============================================================

def calculate_frequency_matrix(sequences, positions, amino_acids):
    """Calculate amino acid frequencies at each position"""
    matrix = np.zeros((len(amino_acids), len(positions)))
    
    for pos_idx, pos in enumerate(positions):
        # Count amino acids at this position
        aa_counts = {aa: 0 for aa in amino_acids}
        total = 0
        
        for seq in sequences:
            if len(seq) > pos:
                aa = seq[pos]
                if aa in amino_acids:
                    aa_counts[aa] += 1
                    total += 1
        
        # Calculate frequencies
        for aa_idx, aa in enumerate(amino_acids):
            if total > 0:
                matrix[aa_idx, pos_idx] = aa_counts[aa] / total
            else:
                matrix[aa_idx, pos_idx] = 0
    
    return matrix

# Positions to analyze (0-indexed, so position 1 is index 1, which is the 2nd amino acid)
positions = list(range(1, max_position))  # positions 2-24 (indices 1-23)
position_labels = [str(i+1) for i in positions]  # labels "2" through "24"

print("\nCalculating amino acid frequencies...")

# Calculate frequencies
unstable_freq = calculate_frequency_matrix(unstable_sequences, positions, amino_acids)
stable_freq = calculate_frequency_matrix(stable_sequences, positions, amino_acids)
overall_freq = calculate_frequency_matrix(all_sequences, positions, amino_acids)

# ============================================================
# Calculate enrichment/depletion relative to overall frequencies
# ============================================================

print("Calculating enrichment ratios...")

# Calculate log2 enrichment (relative to overall library)
# Use pseudocount to avoid division by zero
pseudocount = 0.0001

unstable_enrichment = np.log2((unstable_freq + pseudocount) / (overall_freq + pseudocount))
stable_enrichment = np.log2((stable_freq + pseudocount) / (overall_freq + pseudocount))

# Alternative: Calculate relative to expected uniform distribution (1/20 = 0.05)
# expected_freq = 0.05
# unstable_enrichment = np.log2((unstable_freq + pseudocount) / (expected_freq + pseudocount))
# stable_enrichment = np.log2((stable_freq + pseudocount) / (expected_freq + pseudocount))

print(f"Enrichment range (unstable): {unstable_enrichment.min():.2f} to {unstable_enrichment.max():.2f}")
print(f"Enrichment range (stable): {stable_enrichment.min():.2f} to {stable_enrichment.max():.2f}")

# ============================================================
# Create side-by-side heatmaps
# ============================================================

print("\nCreating heatmaps...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))

# Determine color scale limits (symmetric around 0)
vmin = min(unstable_enrichment.min(), stable_enrichment.min())
vmax = max(unstable_enrichment.max(), stable_enrichment.max())
abs_max = max(abs(vmin), abs(vmax))
vmin, vmax = -abs_max, abs_max

# Clip to reasonable range for visualization
vmin = max(vmin, -0.5)
vmax = min(vmax, 0.5)

# Panel A: Unstable peptides
df_unstable = pd.DataFrame(unstable_enrichment, 
                           index=amino_acids, 
                           columns=position_labels)

sns.heatmap(df_unstable, 
            cmap='RdBu_r', 
            center=0,
            vmin=vmin, 
            vmax=vmax,
            cbar_kws={'label': 'Log2(Frequency Ratio)', 'shrink': 0.8},
            linewidths=0.5,
            linecolor='white',
            ax=ax1,
            square=False,
            xticklabels=True,
            yticklabels=True)

ax1.set_xlabel('Residue position', fontsize=14, weight='bold')
ax1.set_ylabel('Amino acid', fontsize=14, weight='bold')
ax1.set_title(f'A    Unstable peptide-GFP fusions (PSI <= {psi_threshold})', 
              fontsize=16, weight='bold', loc='left', pad=20)

# Panel B: Stable peptides
df_stable = pd.DataFrame(stable_enrichment,
                         index=amino_acids,
                         columns=position_labels)

sns.heatmap(df_stable,
            cmap='RdBu_r',
            center=0,
            vmin=vmin,
            vmax=vmax,
            cbar_kws={'label': 'Log2(Frequency Ratio)', 'shrink': 0.8},
            linewidths=0.5,
            linecolor='white',
            ax=ax2,
            square=False,
            xticklabels=True,
            yticklabels=True)

ax2.set_xlabel('Residue position', fontsize=14, weight='bold')
ax2.set_ylabel('Amino acid', fontsize=14, weight='bold')
ax2.set_title(f'B    Stable peptide-GFP fusions (PSI > {psi_threshold})',
              fontsize=16, weight='bold', loc='left', pad=20)

# Main title
fig.suptitle('The effect of peptide composition on protein stability\n' +
             'Heatmaps showing relative depletion (blue) or enrichment (red) of each amino acid',
             fontsize=18, weight='bold', y=0.98)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'psi_heatmaps_unstable_vs_stable.png'),
            dpi=300, bbox_inches='tight')
plt.close()

print(f"Saved: psi_heatmaps_unstable_vs_stable.png")

# ============================================================
# Create separate high-resolution heatmaps
# ============================================================

# Unstable only (larger)
fig, ax = plt.subplots(1, 1, figsize=(14, 10))

sns.heatmap(df_unstable,
            cmap='RdBu_r',
            center=0,
            vmin=vmin,
            vmax=vmax,
            cbar_kws={'label': 'Log2(Frequency Ratio)', 'shrink': 0.8},
            linewidths=0.5,
            linecolor='white',
            ax=ax,
            square=False,
            xticklabels=True,
            yticklabels=True)

ax.set_xlabel('Residue position', fontsize=14, weight='bold')
ax.set_ylabel('Amino acid', fontsize=14, weight='bold')
ax.set_title(f'Unstable peptide-GFP fusions (PSI <= {psi_threshold})\n' +
             'Relative depletion (blue) or enrichment (red) of each amino acid',
             fontsize=16, weight='bold', pad=20)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'psi_heatmap_unstable_only.png'),
            dpi=300, bbox_inches='tight')
plt.close()

print(f"Saved: psi_heatmap_unstable_only.png")

# Stable only (larger)
fig, ax = plt.subplots(1, 1, figsize=(14, 10))

sns.heatmap(df_stable,
            cmap='RdBu_r',
            center=0,
            vmin=vmin,
            vmax=vmax,
            cbar_kws={'label': 'Log2(Frequency Ratio)', 'shrink': 0.8},
            linewidths=0.5,
            linecolor='white',
            ax=ax,
            square=False,
            xticklabels=True,
            yticklabels=True)

ax.set_xlabel('Residue position', fontsize=14, weight='bold')
ax.set_ylabel('Amino acid', fontsize=14, weight='bold')
ax.set_title(f'Stable peptide-GFP fusions (PSI > {psi_threshold})\n' +
             'Relative depletion (blue) or enrichment (red) of each amino acid',
             fontsize=16, weight='bold', pad=20)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'psi_heatmap_stable_only.png'),
            dpi=300, bbox_inches='tight')
plt.close()

print(f"Saved: psi_heatmap_stable_only.png")

# ============================================================
# Save data to CSV
# ============================================================

print("\nSaving frequency data...")

df_unstable.to_csv(os.path.join(OUTPUT_DIR, 'unstable_aa_enrichment.csv'))
df_stable.to_csv(os.path.join(OUTPUT_DIR, 'stable_aa_enrichment.csv'))

print(f"Saved: unstable_aa_enrichment.csv")
print(f"Saved: stable_aa_enrichment.csv")

# ============================================================
# Summary statistics
# ============================================================

print("\n" + "="*80)
print("SUMMARY - Position 2 Analysis")
print("="*80)

print("\n--- Most ENRICHED at Position 2 in UNSTABLE peptides ---")
pos2_unstable = df_unstable['2'].sort_values(ascending=False)
for aa, val in pos2_unstable.head(5).items():
    freq = unstable_freq[amino_acids.index(aa), 0]
    print(f"  {aa}: {val:+.3f} (frequency: {freq*100:.1f}%)")

print("\n--- Most DEPLETED at Position 2 in UNSTABLE peptides ---")
for aa, val in pos2_unstable.tail(5).items():
    freq = unstable_freq[amino_acids.index(aa), 0]
    print(f"  {aa}: {val:+.3f} (frequency: {freq*100:.1f}%)")

print("\n--- Most ENRICHED at Position 2 in STABLE peptides ---")
pos2_stable = df_stable['2'].sort_values(ascending=False)
for aa, val in pos2_stable.head(5).items():
    freq = stable_freq[amino_acids.index(aa), 0]
    print(f"  {aa}: {val:+.3f} (frequency: {freq*100:.1f}%)")

print("\n--- Most DEPLETED at Position 2 in STABLE peptides ---")
for aa, val in pos2_stable.tail(5).items():
    freq = stable_freq[amino_acids.index(aa), 0]
    print(f"  {aa}: {val:+.3f} (frequency: {freq*100:.1f}%)")

print("\n" + "="*80)
print("ANALYSIS COMPLETE!")
print("="*80)
print(f"\nGenerated files in {OUTPUT_DIR}/:")
print("  1. psi_heatmaps_unstable_vs_stable.png - Side-by-side comparison")
print("  2. psi_heatmap_unstable_only.png - Unstable peptides only")
print("  3. psi_heatmap_stable_only.png - Stable peptides only")
print("  4. unstable_aa_enrichment.csv - Unstable enrichment data")
print("  5. stable_aa_enrichment.csv - Stable enrichment data")
print("\nDone!")
