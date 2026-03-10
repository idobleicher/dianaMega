#!/usr/bin/env python3
"""
Comprehensive analysis of FIRST 5 RESIDUES ONLY for 54 hits
Generates all plot types focused on N-terminal recognition motif
"""

import pandas as pd
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Set output directory
OUTPUT_DIR = "logo_results_54hits_first5"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("="*80)
print("ANALYZING FIRST 5 RESIDUES ONLY (54 HITS)")
print("="*80)

# Load data
print("\nLoading data...")
screen = pd.read_excel('UBR3 Nt screen.xlsx')
hits_file = pd.read_excel('ubr3_best (1).xlsx')

# Get gene list (53 genes)
gene_col = hits_file.columns[2]
genes = hits_file[gene_col].unique()[:54]
print(f"Using {len(genes)} genes")

# Get ONE sequence per gene
hits_data = screen[screen['Gene_ID'].isin(genes)].groupby('Gene_ID').first().reset_index()
hits_sequences = hits_data['AA_seq'].tolist()

# Extract only first 5 residues
hits_sequences_first5 = [seq[:5] for seq in hits_sequences if len(seq) >= 5]
screen_sequences = screen['AA_seq'].tolist()
screen_sequences_first5 = [seq[:5] for seq in screen_sequences if len(seq) >= 5]

print(f"Hits sequences (first 5): {len(hits_sequences_first5)}")
print(f"Screen sequences (first 5): {len(screen_sequences_first5)}")

# Save sequences
pd.DataFrame({
    'gene': hits_data['Gene_ID'],
    'full_sequence': hits_sequences,
    'first_5_residues': hits_sequences_first5
}).to_csv(os.path.join(OUTPUT_DIR, 'hits_first5_sequences.csv'), index=False)

print("\n" + "="*80)
print("CALCULATING ENRICHMENT (POSITIONS 1-5)")
print("="*80)

# Calculate position-specific frequencies
def calc_position_frequencies(sequences, max_pos=5):
    """Calculate frequency of each AA at each position"""
    freq_dict = {}
    
    for pos in range(max_pos):
        aa_counts = Counter()
        for seq in sequences:
            if pos < len(seq):
                aa_counts[seq[pos]] += 1
        
        total = sum(aa_counts.values())
        freq_dict[pos] = {aa: count/total for aa, count in aa_counts.items()}
    
    return freq_dict

hits_freq = calc_position_frequencies(hits_sequences_first5, 5)
screen_freq = calc_position_frequencies(screen_sequences_first5, 5)

# Calculate enrichment ratios
enrichment_data = []
for pos in range(5):  # Positions 0-4 (displayed as 1-5)
    for aa in 'ACDEFGHIKLMNPQRSTVWY':
        hits_f = hits_freq[pos].get(aa, 0.0001)  # Pseudocount
        screen_f = screen_freq[pos].get(aa, 0.0001)
        
        enrichment = hits_f / screen_f
        log2_enrichment = np.log2(enrichment)
        
        enrichment_data.append({
            'position': pos + 1,  # Display as 1-5
            'amino_acid': aa,
            'hits_freq': hits_f,
            'screen_freq': screen_f,
            'hits_count': int(hits_f * len(hits_sequences_first5)),
            'screen_count': int(screen_f * len(screen_sequences_first5)),
            'enrichment': enrichment,
            'log2_enrichment': log2_enrichment
        })

enrichment_df = pd.DataFrame(enrichment_data)
enrichment_df.to_csv(os.path.join(OUTPUT_DIR, 'enrichment_first5_data.csv'), index=False)

print("\nTop enriched amino acids per position:")
for pos in range(1, 6):
    pos_data = enrichment_df[enrichment_df['position'] == pos].nlargest(5, 'enrichment')
    print(f"\nPosition {pos}:")
    for _, row in pos_data.iterrows():
        print(f"  {row['amino_acid']}: {row['enrichment']:.2f}x ({row['hits_freq']:.1%} in hits, {row['screen_freq']:.1%} in screen)")

print("\n" + "="*80)
print("GENERATING PLOTS")
print("="*80)

# Set style
sns.set_style("whitegrid")
plt.rcParams['font.size'] = 10

# 1. ENRICHMENT HEATMAP
print("\n1. Enrichment heatmap...")
fig, ax = plt.subplots(figsize=(8, 10))
pivot_enrichment = enrichment_df.pivot(index='amino_acid', columns='position', values='log2_enrichment')

sns.heatmap(pivot_enrichment, cmap='RdBu_r', center=0, vmin=-2, vmax=2,
            cbar_kws={'label': 'Log₂(Enrichment)'}, linewidths=0.5, linecolor='white', ax=ax)

ax.set_title('Position-Specific Amino Acid Enrichment\nFirst 5 Residues (54 Hits vs Full Screen)', 
             fontsize=14, weight='bold', pad=15)
ax.set_xlabel('Position', fontsize=12, weight='bold')
ax.set_ylabel('Amino Acid', fontsize=12, weight='bold')
ax.set_xticklabels(['1 (M)', '2', '3', '4', '5'], fontsize=11)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'enrichment_heatmap_first5.png'), dpi=300, bbox_inches='tight')
plt.close()

# 2. ANNOTATED ENRICHMENT HEATMAP
print("2. Annotated enrichment heatmap...")
fig, ax = plt.subplots(figsize=(8, 10))

sns.heatmap(pivot_enrichment, cmap='RdBu_r', center=0, vmin=-2.5, vmax=2.5,
            cbar_kws={'label': 'Log₂(Enrichment)'}, linewidths=0.8, linecolor='white', 
            ax=ax, annot=False)

# Annotate highly enriched cells
for i, aa in enumerate(pivot_enrichment.index):
    for j, pos in enumerate(pivot_enrichment.columns):
        val = pivot_enrichment.loc[aa, pos]
        if abs(val) > 1.0:  # Annotate if |log2| > 1 (>2x or <0.5x)
            color = 'white' if abs(val) > 1.5 else 'black'
            ax.text(j + 0.5, i + 0.5, f'{val:.1f}', 
                   ha='center', va='center', fontsize=9, 
                   weight='bold', color=color)

ax.set_title('Annotated Enrichment Heatmap (First 5 Residues)\nValues shown for |Log₂| > 1.0', 
             fontsize=14, weight='bold', pad=15)
ax.set_xlabel('Position', fontsize=12, weight='bold')
ax.set_ylabel('Amino Acid', fontsize=12, weight='bold')
ax.set_xticklabels(['1 (M)', '2', '3', '4', '5'], fontsize=11)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'enrichment_heatmap_annotated_first5.png'), dpi=300, bbox_inches='tight')
plt.close()

# 3. DUAL FREQUENCY HEATMAP (Hits vs Screen)
print("3. Dual frequency heatmap...")
fig = plt.figure(figsize=(10, 12))
gs = fig.add_gridspec(2, 1, hspace=0.3)

# Top: Hits
ax1 = fig.add_subplot(gs[0])
pivot_hits = enrichment_df.pivot(index='amino_acid', columns='position', values='hits_freq')
sns.heatmap(pivot_hits, cmap='YlOrRd', ax=ax1, vmin=0, vmax=0.4,
            cbar_kws={'label': 'Frequency'}, linewidths=0.5, linecolor='gray')
ax1.set_title('A. Hits (53 sequences) - First 5 Residues', fontsize=14, weight='bold', loc='left', pad=10)
ax1.set_xlabel('')
ax1.set_ylabel('Amino Acid', fontsize=12, weight='bold')
ax1.set_xticklabels(['1 (M)', '2', '3', '4', '5'], fontsize=11)
ax1.tick_params(axis='x', labelbottom=False)

# Bottom: Screen
ax2 = fig.add_subplot(gs[1])
pivot_screen = enrichment_df.pivot(index='amino_acid', columns='position', values='screen_freq')
sns.heatmap(pivot_screen, cmap='YlOrRd', ax=ax2, vmin=0, vmax=0.4,
            cbar_kws={'label': 'Frequency'}, linewidths=0.5, linecolor='gray')
ax2.set_title('B. Full Screen (16,514 sequences) - First 5 Residues', fontsize=14, weight='bold', loc='left', pad=10)
ax2.set_xlabel('Position', fontsize=12, weight='bold')
ax2.set_ylabel('Amino Acid', fontsize=12, weight='bold')
ax2.set_xticklabels(['1 (M)', '2', '3', '4', '5'], fontsize=11)

fig.suptitle('Amino Acid Frequency Comparison: First 5 Residues', 
             fontsize=16, weight='bold', y=0.995)
plt.savefig(os.path.join(OUTPUT_DIR, 'frequency_dual_heatmap_first5.png'), 
            dpi=300, bbox_inches='tight')
plt.close()

# 4. POSITION PANELS (5 panels for 5 positions)
print("4. Position-specific panels...")
fig, axes = plt.subplots(1, 5, figsize=(20, 8))

for idx, pos in enumerate(range(1, 6)):
    ax = axes[idx]
    pos_data = enrichment_df[enrichment_df['position'] == pos].sort_values('enrichment', ascending=True)
    
    # Color by enrichment
    colors = []
    for e in pos_data['enrichment']:
        if e > 3:
            colors.append('darkred')
        elif e > 2:
            colors.append('red')
        elif e > 1.5:
            colors.append('orange')
        elif e > 1:
            colors.append('gold')
        elif e > 0.5:
            colors.append('lightblue')
        else:
            colors.append('blue')
    
    ax.barh(pos_data['amino_acid'], pos_data['enrichment'], color=colors, 
            edgecolor='black', linewidth=0.5)
    ax.axvline(x=1, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
    ax.axvline(x=2, color='gray', linestyle=':', linewidth=1, alpha=0.5)
    ax.set_xlabel('Enrichment', fontsize=11, weight='bold')
    
    if pos == 1:
        ax.set_title(f'Position {pos}\n(Methionine)', fontsize=12, weight='bold')
    else:
        ax.set_title(f'Position {pos}', fontsize=12, weight='bold')
    
    ax.grid(axis='x', alpha=0.3)
    ax.set_xlim(0, max(5, pos_data['enrichment'].max() * 1.1))
    
    # Highlight top 3
    if pos > 1:  # Skip position 1 (all Met)
        top3 = pos_data.nlargest(3, 'enrichment')
        for aa in top3['amino_acid']:
            idx_aa = list(pos_data['amino_acid']).index(aa)
            ax.get_yticklabels()[idx_aa].set_weight('bold')
            ax.get_yticklabels()[idx_aa].set_color('darkred')

plt.suptitle('Amino Acid Enrichment: First 5 Residues (54 Hits)', 
             fontsize=16, weight='bold')
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'position_panels_first5.png'), 
            dpi=300, bbox_inches='tight')
plt.close()

# 5. TOP ENRICHED BARS (positions 2-5 only, skip Met)
print("5. Top enriched bars...")
fig, axes = plt.subplots(2, 2, figsize=(14, 12))
axes = axes.flatten()

for idx, pos in enumerate(range(2, 6)):
    ax = axes[idx]
    pos_data = enrichment_df[enrichment_df['position'] == pos].nlargest(10, 'enrichment')
    
    colors = ['darkred' if e > 3 else 'red' if e > 2 else 'orange' if e > 1.5 else 'gold' 
              for e in pos_data['enrichment']]
    
    bars = ax.barh(range(len(pos_data)), pos_data['enrichment'], color=colors, 
                   edgecolor='black', linewidth=1)
    ax.set_yticks(range(len(pos_data)))
    ax.set_yticklabels(pos_data['amino_acid'], fontsize=11, weight='bold')
    ax.set_xlabel('Enrichment Ratio', fontsize=11, weight='bold')
    ax.set_title(f'Position {pos}', fontsize=13, weight='bold')
    ax.axvline(x=1, color='black', linestyle='--', linewidth=1.5)
    ax.axvline(x=2, color='gray', linestyle=':', linewidth=1)
    ax.grid(axis='x', alpha=0.3)
    
    # Add enrichment values
    for i, (aa, e, hf, sf) in enumerate(zip(pos_data['amino_acid'], pos_data['enrichment'],
                                             pos_data['hits_freq'], pos_data['screen_freq'])):
        ax.text(e + 0.15, i, f'{e:.2f}x\n({hf:.1%} vs {sf:.1%})', 
               va='center', fontsize=8, weight='bold')

plt.suptitle('Top 10 Enriched Amino Acids per Position (First 5 Residues)', 
             fontsize=16, weight='bold')
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'top_enriched_bars_first5.png'), 
            dpi=300, bbox_inches='tight')
plt.close()

# 6. COMBINED ENRICHMENT PLOT (Bar graph style)
print("6. Combined enrichment plot...")
fig, ax = plt.subplots(figsize=(14, 8))

positions = range(2, 6)  # Skip position 1
bar_width = 0.15
x_offset = 0

for pos in positions:
    pos_data = enrichment_df[enrichment_df['position'] == pos].nlargest(5, 'enrichment')
    
    for rank, (idx, row) in enumerate(pos_data.iterrows()):
        x = pos + (rank - 2) * bar_width
        height = row['enrichment']
        
        if height > 3:
            color = 'darkred'
        elif height > 2:
            color = 'red'
        elif height > 1.5:
            color = 'orange'
        else:
            color = 'gold'
        
        bar = ax.bar(x, height, width=bar_width, color=color, 
                    edgecolor='black', linewidth=0.8)
        
        # Add AA label
        if height > 1.3:
            ax.text(x, height + 0.15, row['amino_acid'], 
                   ha='center', va='bottom', fontsize=10, weight='bold')

ax.axhline(y=1, color='black', linestyle='-', linewidth=2, alpha=0.6, label='No enrichment')
ax.axhline(y=2, color='gray', linestyle='--', linewidth=1.5, alpha=0.5, label='2x')
ax.axhline(y=3, color='darkgray', linestyle=':', linewidth=1.5, alpha=0.5, label='3x')

ax.set_xlabel('Position', fontsize=14, weight='bold')
ax.set_ylabel('Enrichment Ratio', fontsize=14, weight='bold')
ax.set_title('Top 5 Enriched Amino Acids: First 5 Residues (54 Hits)', 
             fontsize=16, weight='bold')
ax.set_xticks(range(2, 6))
ax.set_xticklabels(['2', '3', '4', '5'], fontsize=12)
ax.set_xlim(1.5, 5.5)
ax.set_ylim(0, max(enrichment_df[enrichment_df['position'] > 1]['enrichment']) * 1.2)
ax.legend(loc='upper right', fontsize=12)
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'enrichment_combined_first5.png'), 
            dpi=300, bbox_inches='tight')
plt.close()

# 7. OVERALL AA ENRICHMENT (averaged across positions 2-5)
print("7. Overall amino acid enrichment...")
overall_data = enrichment_df[enrichment_df['position'] > 1].groupby('amino_acid').agg({
    'hits_freq': 'mean',
    'screen_freq': 'mean',
    'enrichment': 'mean'
}).reset_index()
overall_data = overall_data.sort_values('enrichment', ascending=True)

fig, ax = plt.subplots(figsize=(10, 10))

colors = []
for e in overall_data['enrichment']:
    if e > 1.3:
        colors.append('darkred')
    elif e > 1.15:
        colors.append('red')
    elif e > 1:
        colors.append('orange')
    elif e > 0.85:
        colors.append('lightblue')
    else:
        colors.append('blue')

bars = ax.barh(overall_data['amino_acid'], overall_data['enrichment'], 
               color=colors, edgecolor='black', linewidth=1.5)

ax.axvline(x=1, color='black', linestyle='-', linewidth=2.5, label='No enrichment', alpha=0.7)
ax.set_xlabel('Mean Enrichment Ratio', fontsize=13, weight='bold')
ax.set_ylabel('Amino Acid', fontsize=13, weight='bold')
ax.set_title('Overall Amino Acid Enrichment (Average, Positions 2-5)\n54 Hits vs Full Screen', 
             fontsize=14, weight='bold')
ax.grid(axis='x', alpha=0.3)
ax.legend(fontsize=11)

# Add values
for i, (aa, enrich, hf, sf) in enumerate(zip(overall_data['amino_acid'], 
                                              overall_data['enrichment'],
                                              overall_data['hits_freq'],
                                              overall_data['screen_freq'])):
    x_pos = enrich + 0.03 if enrich > 1 else enrich - 0.03
    ha = 'left' if enrich > 1 else 'right'
    ax.text(x_pos, i, f'{enrich:.2f}x', va='center', ha=ha, fontsize=9, weight='bold')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'overall_enrichment_first5.png'), 
            dpi=300, bbox_inches='tight')
plt.close()

# 8. STACKED FREQUENCY PLOT
print("8. Stacked frequency plot...")
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

positions = [1, 2, 3, 4, 5]
amino_acids = sorted(enrichment_df['amino_acid'].unique())

# Hits stacked
for aa in amino_acids:
    aa_data = enrichment_df[enrichment_df['amino_acid'] == aa].sort_values('position')
    bottom_vals = [enrichment_df[(enrichment_df['position'] == p) & 
                                (enrichment_df['amino_acid'] < aa)]['hits_freq'].sum() 
                   for p in aa_data['position']]
    
    label = aa if aa in ['M', 'P', 'G', 'D', 'E', 'Y', 'I', 'K'] else ''
    ax1.bar(aa_data['position'], aa_data['hits_freq'], bottom=bottom_vals, label=label)

ax1.set_ylabel('Frequency (Cumulative)', fontsize=12, weight='bold')
ax1.set_title('A. Hits (53 sequences)', fontsize=13, weight='bold', loc='left')
ax1.legend(ncol=4, loc='upper right', title='Key AAs', fontsize=9)
ax1.grid(axis='y', alpha=0.3)
ax1.set_ylim(0, 1.05)

# Screen stacked
for aa in amino_acids:
    aa_data = enrichment_df[enrichment_df['amino_acid'] == aa].sort_values('position')
    bottom_vals = [enrichment_df[(enrichment_df['position'] == p) & 
                                (enrichment_df['amino_acid'] < aa)]['screen_freq'].sum() 
                   for p in aa_data['position']]
    ax2.bar(aa_data['position'], aa_data['screen_freq'], bottom=bottom_vals)

ax2.set_xlabel('Position', fontsize=12, weight='bold')
ax2.set_ylabel('Frequency (Cumulative)', fontsize=12, weight='bold')
ax2.set_title('B. Full Screen (16,514 sequences)', fontsize=13, weight='bold', loc='left')
ax2.grid(axis='y', alpha=0.3)
ax2.set_xticks(positions)
ax2.set_xticklabels(['1 (M)', '2', '3', '4', '5'], fontsize=11)
ax2.set_ylim(0, 1.05)

fig.suptitle('Stacked Amino Acid Frequencies: First 5 Residues', 
             fontsize=16, weight='bold', y=0.995)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'stacked_frequencies_first5.png'), 
            dpi=300, bbox_inches='tight')
plt.close()

# 9. ENRICHMENT LOGO STYLE
print("9. Enrichment logo style...")
fig, ax = plt.subplots(figsize=(10, 8))

for pos in range(2, 6):  # Positions 2-5
    pos_data = enrichment_df[enrichment_df['position'] == pos].nlargest(5, 'enrichment')
    
    for rank, (idx, row) in enumerate(pos_data.iterrows()):
        x = pos + (rank - 2) * 0.12
        height = row['enrichment']
        
        if height > 3:
            color = 'darkred'
        elif height > 2:
            color = 'red'
        elif height > 1.5:
            color = 'orange'
        else:
            color = 'gold'
        
        ax.bar(x, height, width=0.12, color=color, edgecolor='black', linewidth=0.8)
        
        if height > 1.3:
            ax.text(x, height + 0.12, row['amino_acid'], 
                   ha='center', va='bottom', fontsize=11, weight='bold',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                           edgecolor='gray', alpha=0.7))

ax.axhline(y=1, color='black', linestyle='-', linewidth=2, alpha=0.5)
ax.axhline(y=2, color='gray', linestyle='--', linewidth=1.5, alpha=0.5)
ax.axhline(y=3, color='darkgray', linestyle=':', linewidth=1.5, alpha=0.5)

ax.set_xlabel('Position', fontsize=14, weight='bold')
ax.set_ylabel('Enrichment Ratio', fontsize=14, weight='bold')
ax.set_title('Enrichment Logo: First 5 Residues (54 Hits vs Full Screen)', 
             fontsize=16, weight='bold')
ax.set_xticks(range(1, 6))
ax.set_xticklabels(['1\n(M)', '2', '3', '4', '5'], fontsize=12)
ax.set_xlim(1.5, 5.5)
ax.set_ylim(0, 5.5)
ax.grid(axis='y', alpha=0.3)

# Add legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='darkred', edgecolor='black', label='>3x'),
    Patch(facecolor='red', edgecolor='black', label='2-3x'),
    Patch(facecolor='orange', edgecolor='black', label='1.5-2x'),
    Patch(facecolor='gold', edgecolor='black', label='1-1.5x')
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=11, title='Enrichment')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'enrichment_logo_first5.png'), 
            dpi=300, bbox_inches='tight')
plt.close()

# 10. SIMPLE BAR COMPARISON (positions 2-5, top 3 each)
print("10. Simple bar comparison...")
fig, axes = plt.subplots(1, 4, figsize=(16, 6), sharey=True)

for idx, pos in enumerate(range(2, 6)):
    ax = axes[idx]
    pos_data = enrichment_df[enrichment_df['position'] == pos].nlargest(8, 'enrichment')
    
    colors = ['darkred' if e > 3 else 'red' if e > 2 else 'orange' for e in pos_data['enrichment']]
    
    ax.bar(range(len(pos_data)), pos_data['enrichment'], color=colors, 
           edgecolor='black', linewidth=1)
    ax.set_xticks(range(len(pos_data)))
    ax.set_xticklabels(pos_data['amino_acid'], fontsize=11, weight='bold')
    ax.set_ylabel('Enrichment' if idx == 0 else '', fontsize=12, weight='bold')
    ax.set_title(f'Position {pos}', fontsize=13, weight='bold')
    ax.axhline(y=1, color='black', linestyle='--', linewidth=1.5)
    ax.axhline(y=2, color='gray', linestyle=':', linewidth=1)
    ax.grid(axis='y', alpha=0.3)
    ax.set_ylim(0, 5)

plt.suptitle('Top 8 Enriched Amino Acids: Positions 2-5 (54 Hits)', 
             fontsize=16, weight='bold')
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'simple_bars_first5.png'), 
            dpi=300, bbox_inches='tight')
plt.close()

print("\n" + "="*80)
print("ANALYSIS COMPLETE!")
print("="*80)
print(f"\nAll results saved to: {OUTPUT_DIR}/")
print(f"\nGenerated files:")
print(f"  1. hits_first5_sequences.csv - Sequences data")
print(f"  2. enrichment_first5_data.csv - Enrichment data")
print(f"  3. enrichment_heatmap_first5.png - Basic heatmap")
print(f"  4. enrichment_heatmap_annotated_first5.png - Annotated heatmap")
print(f"  5. frequency_dual_heatmap_first5.png - Dual frequency heatmap")
print(f"  6. position_panels_first5.png - 5 position panels")
print(f"  7. top_enriched_bars_first5.png - Top enriched bars")
print(f"  8. enrichment_combined_first5.png - Combined enrichment")
print(f"  9. overall_enrichment_first5.png - Overall AA enrichment")
print(f"  10. stacked_frequencies_first5.png - Stacked frequencies")
print(f"  11. enrichment_logo_first5.png - Logo-style plot")
print(f"  12. simple_bars_first5.png - Simple bar comparison")
print("\nTotal: 12 plots + 2 data files")
print("\nDone!")
