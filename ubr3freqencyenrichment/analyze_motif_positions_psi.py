#!/usr/bin/env python3
"""
Analyze the occurrence of specific motifs (PD, PE, GE, GD, PT) 
across all positions in unstable vs stable peptides
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

OUTPUT_DIR = "logo_results_motif_analysis"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("="*80)
print("MOTIF POSITION ANALYSIS - UNSTABLE VS STABLE PEPTIDES")
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

# Target motifs
target_motifs = ['PD', 'PE', 'GE', 'GD', 'PT']

print(f"\nTarget motifs: {', '.join(target_motifs)}")

# ============================================================
# ANALYSIS 1: Count motifs anywhere in sequences
# ============================================================

print("\n" + "="*80)
print("ANALYSIS 1: MOTIF FREQUENCY IN SEQUENCES")
print("="*80)

def count_motif_occurrences(sequences, motif):
    """Count how many sequences contain the motif (at any position)"""
    count = 0
    for seq in sequences:
        if motif in seq:
            count += 1
    return count

def count_motif_total_occurrences(sequences, motif):
    """Count total number of times motif appears across all sequences"""
    total = 0
    for seq in sequences:
        total += seq.count(motif)
    return total

motif_results = []

for motif in target_motifs:
    # Count in unstable
    unstable_seqs_with_motif = count_motif_occurrences(unstable_sequences, motif)
    unstable_total_occurrences = count_motif_total_occurrences(unstable_sequences, motif)
    unstable_freq = unstable_seqs_with_motif / len(unstable_sequences)
    
    # Count in stable
    stable_seqs_with_motif = count_motif_occurrences(stable_sequences, motif)
    stable_total_occurrences = count_motif_total_occurrences(stable_sequences, motif)
    stable_freq = stable_seqs_with_motif / len(stable_sequences)
    
    # Enrichment
    enrichment = unstable_freq / stable_freq if stable_freq > 0 else float('inf')
    log2_enrichment = np.log2(enrichment) if enrichment > 0 and enrichment != float('inf') else 0
    
    # Fisher's exact test
    contingency_table = [
        [unstable_seqs_with_motif, len(unstable_sequences) - unstable_seqs_with_motif],
        [stable_seqs_with_motif, len(stable_sequences) - stable_seqs_with_motif]
    ]
    odds_ratio, p_value = stats.fisher_exact(contingency_table)
    
    # Adjusted p-value
    p_adjusted = min(p_value * len(target_motifs), 1.0)
    
    # Significance
    if p_adjusted < 0.001:
        significance = '***'
    elif p_adjusted < 0.01:
        significance = '**'
    elif p_adjusted < 0.05:
        significance = '*'
    else:
        significance = 'ns'
    
    motif_results.append({
        'Motif': motif,
        'Unstable_Seqs_With_Motif': unstable_seqs_with_motif,
        'Unstable_Freq': unstable_freq,
        'Unstable_Total_Occurrences': unstable_total_occurrences,
        'Stable_Seqs_With_Motif': stable_seqs_with_motif,
        'Stable_Freq': stable_freq,
        'Stable_Total_Occurrences': stable_total_occurrences,
        'Enrichment_Unstable_vs_Stable': enrichment,
        'Log2_Enrichment': log2_enrichment,
        'P_value': p_value,
        'P_adjusted': p_adjusted,
        'Significance': significance
    })

results_df = pd.DataFrame(motif_results)

print("\n--- MOTIF OCCURRENCE IN SEQUENCES ---")
print(results_df.to_string(index=False))

results_df.to_csv(os.path.join(OUTPUT_DIR, 'motif_occurrence_unstable_vs_stable.csv'), index=False)
print(f"\nSaved: motif_occurrence_unstable_vs_stable.csv")

# ============================================================
# ANALYSIS 2: Position-specific motif occurrence
# ============================================================

print("\n" + "="*80)
print("ANALYSIS 2: MOTIF OCCURRENCE BY POSITION")
print("="*80)

def count_motif_at_position(sequences, motif, position):
    """Count motif occurrences starting at a specific position (0-indexed)"""
    count = 0
    for seq in sequences:
        if position + len(motif) <= len(seq):
            if seq[position:position + len(motif)] == motif:
                count += 1
    return count

# Analyze positions 1-23 (positions 2-24, allowing for 2-character motifs)
max_position = 23
positions = list(range(1, max_position))  # positions 2-23 (0-indexed: 1-22)

# Create position-specific frequency matrix
unstable_position_freq = np.zeros((len(target_motifs), len(positions)))
stable_position_freq = np.zeros((len(target_motifs), len(positions)))

for motif_idx, motif in enumerate(target_motifs):
    for pos_idx, pos in enumerate(positions):
        unstable_count = count_motif_at_position(unstable_sequences, motif, pos)
        stable_count = count_motif_at_position(stable_sequences, motif, pos)
        
        unstable_position_freq[motif_idx, pos_idx] = unstable_count / len(unstable_sequences)
        stable_position_freq[motif_idx, pos_idx] = stable_count / len(stable_sequences)

# Calculate enrichment (log2 ratio)
pseudocount = 0.0001
enrichment_matrix = np.log2((unstable_position_freq + pseudocount) / (stable_position_freq + pseudocount))

print(f"\nAnalyzed positions 2-{max_position}")

# Save position data
position_labels = [str(i+1) for i in positions]
df_enrichment = pd.DataFrame(enrichment_matrix, 
                              index=target_motifs, 
                              columns=position_labels)
df_enrichment.to_csv(os.path.join(OUTPUT_DIR, 'motif_position_enrichment.csv'))
print(f"Saved: motif_position_enrichment.csv")

# ============================================================
# VISUALIZATION 1: Overall motif frequency comparison
# ============================================================

fig, axes = plt.subplots(2, 2, figsize=(16, 12))

# Panel A: Percentage of sequences containing motif
ax1 = axes[0, 0]
x = np.arange(len(target_motifs))
width = 0.35

bars1 = ax1.bar(x - width/2, results_df['Unstable_Freq']*100, width,
                label='Unstable (PSI <= 3.6)', color='#3498db', edgecolor='black', linewidth=1.5)
bars2 = ax1.bar(x + width/2, results_df['Stable_Freq']*100, width,
                label='Stable (PSI > 3.6)', color='#e74c3c', edgecolor='black', linewidth=1.5)

ax1.set_xlabel('Motif', fontsize=13, weight='bold')
ax1.set_ylabel('% of Sequences Containing Motif', fontsize=13, weight='bold')
ax1.set_title('A. Motif Frequency in Sequences', fontsize=14, weight='bold', pad=15)
ax1.set_xticks(x)
ax1.set_xticklabels(target_motifs, fontsize=12, weight='bold')
ax1.legend(fontsize=11)
ax1.grid(axis='y', alpha=0.3)
ax1.set_facecolor('#F5F5F5')

# Add percentage labels on bars
for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.1f}%', ha='center', va='bottom', fontsize=9)

# Panel B: Log2 enrichment
ax2 = axes[0, 1]
colors = ['#CC3333' if sig != 'ns' else '#6699CC' for sig in results_df['Significance']]
bars = ax2.bar(results_df['Motif'], results_df['Log2_Enrichment'],
               color=colors, edgecolor='black', linewidth=2)

# Add significance stars
for i, (motif, log2_enr, sig) in enumerate(zip(results_df['Motif'],
                                                  results_df['Log2_Enrichment'],
                                                  results_df['Significance'])):
    if sig != 'ns':
        y_pos = log2_enr + 0.05 if log2_enr > 0 else log2_enr - 0.1
        ax2.text(i, y_pos, sig, ha='center', va='bottom' if log2_enr > 0 else 'top',
                fontsize=16, weight='bold')

ax2.axhline(y=0, color='black', linestyle='-', linewidth=1)
ax2.set_xlabel('Motif', fontsize=13, weight='bold')
ax2.set_ylabel('Log2 Enrichment\n(Unstable vs Stable)', fontsize=13, weight='bold')
ax2.set_title('B. Motif Enrichment in Unstable Peptides', fontsize=14, weight='bold', pad=15)
ax2.set_xticklabels(target_motifs, fontsize=12, weight='bold')
ax2.grid(axis='y', alpha=0.3)
ax2.set_facecolor('#F5F5F5')

# Panel C: Total occurrences
ax3 = axes[1, 0]
x = np.arange(len(target_motifs))

bars1 = ax3.bar(x - width/2, results_df['Unstable_Total_Occurrences'], width,
                label='Unstable (PSI <= 3.6)', color='#3498db', edgecolor='black', linewidth=1.5)
bars2 = ax3.bar(x + width/2, results_df['Stable_Total_Occurrences'], width,
                label='Stable (PSI > 3.6)', color='#e74c3c', edgecolor='black', linewidth=1.5)

ax3.set_xlabel('Motif', fontsize=13, weight='bold')
ax3.set_ylabel('Total Number of Occurrences', fontsize=13, weight='bold')
ax3.set_title('C. Total Motif Occurrences Across All Sequences', fontsize=14, weight='bold', pad=15)
ax3.set_xticks(x)
ax3.set_xticklabels(target_motifs, fontsize=12, weight='bold')
ax3.legend(fontsize=11)
ax3.grid(axis='y', alpha=0.3)
ax3.set_facecolor('#F5F5F5')

# Panel D: Summary table
ax4 = axes[1, 1]
ax4.axis('tight')
ax4.axis('off')

table_data = []
table_data.append(['Motif', 'Unstable\n(%)', 'Stable\n(%)', 'Enrichment', 'P-value', 'Sig'])

for _, row in results_df.iterrows():
    table_data.append([
        row['Motif'],
        f"{row['Unstable_Freq']*100:.1f}%",
        f"{row['Stable_Freq']*100:.1f}%",
        f"{row['Enrichment_Unstable_vs_Stable']:.2f}x",
        f"{row['P_adjusted']:.2e}",
        row['Significance']
    ])

table = ax4.table(cellText=table_data, cellLoc='center', loc='center',
                  colWidths=[0.12, 0.15, 0.15, 0.15, 0.18, 0.1])
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 2.5)

for i in range(6):
    cell = table[(0, i)]
    cell.set_facecolor('#4A90E2')
    cell.set_text_props(weight='bold', color='white')

for i in range(1, len(table_data)):
    if table_data[i][5] != 'ns':
        for j in range(6):
            table[(i, j)].set_facecolor('#FFE6E6')

ax4.set_title('D. Summary Statistics', fontsize=14, weight='bold', pad=20)

plt.suptitle('Motif Occurrence in Unstable vs Stable Peptides',
             fontsize=16, weight='bold', y=0.995)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'motif_occurrence_comparison.png'),
            dpi=300, bbox_inches='tight')
plt.close()

print(f"Saved: motif_occurrence_comparison.png")

# ============================================================
# VISUALIZATION 2: Position-specific heatmap
# ============================================================

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(16, 12))

# Panel 1: Unstable peptides
df_unstable_pos = pd.DataFrame(unstable_position_freq * 100,
                                index=target_motifs,
                                columns=position_labels)

sns.heatmap(df_unstable_pos, cmap='Blues', vmin=0, vmax=2,
            cbar_kws={'label': '% of Sequences', 'shrink': 0.8},
            linewidths=0.5, linecolor='white', ax=ax1,
            annot=False, fmt='.1f')

ax1.set_xlabel('Position (start of motif)', fontsize=13, weight='bold')
ax1.set_ylabel('Motif', fontsize=13, weight='bold')
ax1.set_title('A. Motif Position Frequency in UNSTABLE Peptides (PSI <= 3.6)',
              fontsize=14, weight='bold', pad=15)

# Panel 2: Stable peptides
df_stable_pos = pd.DataFrame(stable_position_freq * 100,
                              index=target_motifs,
                              columns=position_labels)

sns.heatmap(df_stable_pos, cmap='Reds', vmin=0, vmax=2,
            cbar_kws={'label': '% of Sequences', 'shrink': 0.8},
            linewidths=0.5, linecolor='white', ax=ax2,
            annot=False, fmt='.1f')

ax2.set_xlabel('Position (start of motif)', fontsize=13, weight='bold')
ax2.set_ylabel('Motif', fontsize=13, weight='bold')
ax2.set_title('B. Motif Position Frequency in STABLE Peptides (PSI > 3.6)',
              fontsize=14, weight='bold', pad=15)

# Panel 3: Enrichment (log2 ratio)
sns.heatmap(df_enrichment, cmap='RdBu_r', center=0, vmin=-1, vmax=1,
            cbar_kws={'label': 'Log2(Unstable/Stable)', 'shrink': 0.8},
            linewidths=0.5, linecolor='white', ax=ax3,
            annot=False, fmt='.2f')

ax3.set_xlabel('Position (start of motif)', fontsize=13, weight='bold')
ax3.set_ylabel('Motif', fontsize=13, weight='bold')
ax3.set_title('C. Position-Specific Enrichment (Red = Enriched in Unstable, Blue = Enriched in Stable)',
              fontsize=14, weight='bold', pad=15)

plt.suptitle('Position-Specific Motif Analysis',
             fontsize=16, weight='bold', y=0.995)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'motif_position_heatmaps.png'),
            dpi=300, bbox_inches='tight')
plt.close()

print(f"Saved: motif_position_heatmaps.png")

# ============================================================
# SUMMARY
# ============================================================

print("\n" + "="*80)
print("SUMMARY - KEY FINDINGS")
print("="*80)

print("\n--- Motif Frequency (% of sequences containing motif) ---")
for _, row in results_df.iterrows():
    print(f"\n{row['Motif']}:")
    print(f"  Unstable: {row['Unstable_Seqs_With_Motif']:,} seqs ({row['Unstable_Freq']*100:.2f}%)")
    print(f"  Stable:   {row['Stable_Seqs_With_Motif']:,} seqs ({row['Stable_Freq']*100:.2f}%)")
    print(f"  Enrichment: {row['Enrichment_Unstable_vs_Stable']:.2f}x")
    print(f"  P-value: {row['P_adjusted']:.2e} {row['Significance']}")

print("\n--- Most common positions for each motif in UNSTABLE peptides ---")
for motif_idx, motif in enumerate(target_motifs):
    top_positions = df_unstable_pos.iloc[motif_idx].nlargest(3)
    print(f"\n{motif}:")
    for pos, freq in top_positions.items():
        if freq > 0:
            print(f"  Position {pos}: {freq:.2f}%")

print("\n" + "="*80)
print("ANALYSIS COMPLETE!")
print("="*80)
print(f"\nGenerated files in {OUTPUT_DIR}/:")
print("  1. motif_occurrence_unstable_vs_stable.csv")
print("  2. motif_position_enrichment.csv")
print("  3. motif_occurrence_comparison.png")
print("  4. motif_position_heatmaps.png")
print("\nDone!")
