#!/usr/bin/env python3
"""
Create publication-quality bar graph showing enrichment of first 5 residues
(positions 2-6) with statistical significance
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import fisher_exact
import os

print("\n" + "="*80)
print("ENRICHMENT ANALYSIS: FIRST 5 RESIDUES (Positions 2-6)")
print("="*80 + "\n")

# Load PWMs
pwm_hits = pd.read_csv('logo_results_hits/position_weight_matrix.csv', index_col=0)
pwm_screen = pd.read_csv('logo_results_full_screen/position_weight_matrix.csv', index_col=0)

print(f"Loaded PWM data for positions 1-{pwm_hits.index.max()}")

# Total sequences
n_hits = 91
n_screen = 16514

print(f"\nSequences: {n_hits} hits vs {n_screen} screen")

# Analyze positions 2-6
positions = [2, 3, 4, 5, 6]
position_labels = ['Position 2', 'Position 3', 'Position 4', 'Position 5', 'Position 6']

enrichment_data = []

for pos in positions:
    print(f"\n{'='*60}")
    print(f"POSITION {pos}")
    print(f"{'='*60}")
    
    # Get frequencies
    hits_freq = pwm_hits.loc[pos]
    screen_freq = pwm_screen.loc[pos]
    
    # Get counts (frequency * total sequences)
    hits_counts = (hits_freq * n_hits).round().astype(int)
    screen_counts = (screen_freq * n_screen).round().astype(int)
    
    # Calculate enrichment for ALL amino acids
    pseudocount = 0.001
    enrichments = {}
    for aa in hits_freq.index:
        hits_f = hits_freq[aa]
        screen_f = screen_freq[aa]
        enrichment = (hits_f + pseudocount) / (screen_f + pseudocount)
        enrichments[aa] = enrichment
    
    # Find top 5 amino acids by ENRICHMENT (not frequency!)
    enrichment_series = pd.Series(enrichments)
    top_aas = enrichment_series.nlargest(5).index.tolist()
    
    print(f"\nTop 5 amino acids by ENRICHMENT at position {pos}:")
    
    for aa in top_aas:
        hits_count = hits_counts[aa]
        screen_count = screen_counts[aa]
        hits_f = hits_freq[aa]
        screen_f = screen_freq[aa]
        
        # Use pre-calculated enrichment
        enrichment = enrichments[aa]
        log2_enr = np.log2(enrichment)
        
        # Fisher's exact test for significance
        # Contingency table: [[hits with AA, hits without AA], [screen with AA, screen without AA]]
        contingency = [
            [hits_count, n_hits - hits_count],
            [screen_count, n_screen - screen_count]
        ]
        
        try:
            oddsratio, p_value = fisher_exact(contingency, alternative='two-sided')
        except:
            p_value = 1.0
            oddsratio = enrichment
        
        # Significance level
        if p_value < 0.001:
            sig = '***'
        elif p_value < 0.01:
            sig = '**'
        elif p_value < 0.05:
            sig = '*'
        else:
            sig = 'ns'
        
        print(f"  {aa}: {hits_f:.3f} ({hits_count}) vs {screen_f:.3f} ({screen_count}) | "
              f"Enrichment: {enrichment:.2f}x | p={p_value:.2e} {sig}")
        
        enrichment_data.append({
            'position': pos,
            'amino_acid': aa,
            'hits_freq': hits_f,
            'screen_freq': screen_f,
            'hits_count': hits_count,
            'screen_count': screen_count,
            'enrichment': enrichment,
            'log2_enrichment': log2_enr,
            'p_value': p_value,
            'odds_ratio': oddsratio,
            'significance': sig
        })

# Convert to DataFrame
df_enrichment = pd.DataFrame(enrichment_data)

# Save data
df_enrichment.to_csv('logo_results_enrichment/first_five_residues_enrichment.csv', index=False)
print(f"\n\nSaved enrichment data to: logo_results_enrichment/first_five_residues_enrichment.csv")

# Create visualization
print("\nCreating visualization...")

# Figure 1: Comprehensive multi-panel plot
fig = plt.figure(figsize=(20, 12))
gs = fig.add_gridspec(3, 2, hspace=0.35, wspace=0.3)

# Color scheme for amino acids (by properties)
aa_colors = {
    # Hydrophobic
    'A': '#8DD3C7', 'V': '#8DD3C7', 'I': '#8DD3C7', 'L': '#8DD3C7', 'M': '#8DD3C7',
    'F': '#80B1D3', 'W': '#80B1D3', 'Y': '#80B1D3',  # Aromatic
    'P': '#FDB462', 'G': '#FDB462',  # Special
    'S': '#B3DE69', 'T': '#B3DE69', 'C': '#B3DE69', 'N': '#B3DE69', 'Q': '#B3DE69',  # Polar
    'K': '#FB8072', 'R': '#FB8072', 'H': '#FB8072',  # Basic
    'D': '#BEBADA', 'E': '#BEBADA'  # Acidic
}

# Panel 1: Enrichment by position (grouped bar)
ax1 = fig.add_subplot(gs[0, :])

for i, pos in enumerate(positions):
    pos_data = df_enrichment[df_enrichment['position'] == pos].nlargest(5, 'enrichment')
    
    x_positions = np.arange(len(pos_data)) + i * (len(pos_data) + 1)
    colors = [aa_colors.get(aa, '#CCCCCC') for aa in pos_data['amino_acid']]
    
    bars = ax1.bar(x_positions, pos_data['enrichment'], 
                   color=colors, edgecolor='black', linewidth=1.5, alpha=0.8,
                   label=f'Position {pos}' if i == 0 else '')
    
    # Add significance stars
    for j, (idx, row) in enumerate(pos_data.iterrows()):
        if row['significance'] != 'ns':
            ax1.text(x_positions[j], row['enrichment'] + 0.15,
                    row['significance'], ha='center', va='bottom',
                    fontsize=12, fontweight='bold', color='red')
        
        # Add amino acid labels
        ax1.text(x_positions[j], -0.3, row['amino_acid'],
                ha='center', va='top', fontsize=11, fontweight='bold')
    
    # Add position label
    mid_x = x_positions[len(pos_data)//2]
    ax1.text(mid_x, -0.8, f'Position {pos}',
            ha='center', va='top', fontsize=12, fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.5))

ax1.axhline(y=1, color='black', linestyle='--', linewidth=2, alpha=0.5, label='No enrichment')
ax1.set_ylabel('Enrichment Ratio (Hits/Screen)', fontsize=14, fontweight='bold')
ax1.set_title('Top 5 Enriched Amino Acids at Positions 2-6\n*** p<0.001, ** p<0.01, * p<0.05', 
             fontsize=16, fontweight='bold', pad=20)
ax1.set_ylim(-1, max(df_enrichment['enrichment']) * 1.2)
ax1.set_xlim(-1, max(x_positions) + 2)
ax1.set_xticks([])
ax1.grid(axis='y', alpha=0.3)

# Panel 2: Heatmap of log2 enrichment
ax2 = fig.add_subplot(gs[1, 0])

# Create pivot table for heatmap
pivot_data = df_enrichment.pivot_table(values='log2_enrichment', 
                                       index='amino_acid', 
                                       columns='position', 
                                       fill_value=0)

sns.heatmap(pivot_data, cmap='RdBu_r', center=0, 
           cbar_kws={'label': 'Log2 Enrichment'}, 
           ax=ax2, linewidths=2, linecolor='white',
           vmin=-2, vmax=2.5, annot=True, fmt='.2f',
           annot_kws={'fontsize': 10, 'fontweight': 'bold'})
ax2.set_xlabel('Position', fontsize=12, fontweight='bold')
ax2.set_ylabel('Amino Acid', fontsize=12, fontweight='bold')
ax2.set_title('Log2 Enrichment Heatmap (Positions 2-6)', fontsize=14, fontweight='bold')

# Panel 3: P-value significance plot
ax3 = fig.add_subplot(gs[1, 1])

for pos in positions:
    pos_data = df_enrichment[df_enrichment['position'] == pos].nlargest(5, 'enrichment')
    
    # Plot -log10(p-value)
    neg_log_p = -np.log10(pos_data['p_value'])
    
    x_pos = pos + np.linspace(-0.15, 0.15, len(pos_data))
    
    for x, y, aa, sig in zip(x_pos, neg_log_p, pos_data['amino_acid'], pos_data['significance']):
        color = 'red' if sig != 'ns' else 'gray'
        ax3.scatter(x, y, s=200, alpha=0.7, color=color, edgecolors='black', linewidth=2)
        ax3.text(x, y, aa, ha='center', va='center', fontsize=10, fontweight='bold')

# Significance thresholds
ax3.axhline(y=-np.log10(0.05), color='orange', linestyle='--', linewidth=2, label='p=0.05')
ax3.axhline(y=-np.log10(0.01), color='red', linestyle='--', linewidth=2, label='p=0.01')
ax3.axhline(y=-np.log10(0.001), color='darkred', linestyle='--', linewidth=2, label='p=0.001')

ax3.set_xlabel('Position', fontsize=12, fontweight='bold')
ax3.set_ylabel('-log10(p-value)', fontsize=12, fontweight='bold')
ax3.set_title('Statistical Significance (Fisher\'s Exact Test)', fontsize=14, fontweight='bold')
ax3.set_xticks(positions)
ax3.set_xticklabels([f'Pos {p}' for p in positions])
ax3.legend(loc='upper right')
ax3.grid(alpha=0.3)

# Panel 4: Frequency comparison
ax4 = fig.add_subplot(gs[2, :])

x = np.arange(len(positions))
width = 0.15

for i, (pos_idx, pos) in enumerate(zip(x, positions)):
    pos_data = df_enrichment[df_enrichment['position'] == pos].nlargest(3, 'enrichment')
    
    for j, (idx, row) in enumerate(pos_data.iterrows()):
        offset = (j - 1) * width
        
        # Hits bar
        ax4.bar(pos_idx + offset - 0.2, row['hits_freq'], width, 
               color=aa_colors.get(row['amino_acid'], '#CCCCCC'),
               edgecolor='black', linewidth=1, alpha=0.8)
        
        # Screen bar (lighter)
        ax4.bar(pos_idx + offset + 0.2, row['screen_freq'], width,
               color=aa_colors.get(row['amino_acid'], '#CCCCCC'),
               edgecolor='black', linewidth=1, alpha=0.3)
        
        # Add AA label
        if j == 1:  # Middle AA
            ax4.text(pos_idx + offset, max(row['hits_freq'], row['screen_freq']) + 0.02,
                    row['amino_acid'], ha='center', va='bottom',
                    fontsize=10, fontweight='bold')

ax4.set_xlabel('Position', fontsize=12, fontweight='bold')
ax4.set_ylabel('Frequency', fontsize=12, fontweight='bold')
ax4.set_title('Top 3 Amino Acids: Frequency Comparison (Dark=Hits, Light=Screen)', 
             fontsize=14, fontweight='bold')
ax4.set_xticks(x)
ax4.set_xticklabels([f'Position {p}' for p in positions])
ax4.grid(axis='y', alpha=0.3)

plt.suptitle('UBR3 Recognition Motif: First 5 Residues Enrichment Analysis',
            fontsize=18, fontweight='bold', y=0.995)

plt.savefig('logo_results_enrichment/first_five_residues_enrichment_comprehensive.png',
           dpi=300, bbox_inches='tight')
print("[OK] Saved comprehensive plot: first_five_residues_enrichment_comprehensive.png")
plt.close()

# Figure 2: Clean publication-quality bar graph
fig, ax = plt.subplots(figsize=(16, 8))

# Get top enriched at each position
x_offset = 0
x_labels = []
x_positions = []

for i, pos in enumerate(positions):
    pos_data = df_enrichment[df_enrichment['position'] == pos].nlargest(5, 'enrichment')
    
    x_pos = np.arange(len(pos_data)) + x_offset
    colors = [aa_colors.get(aa, '#CCCCCC') for aa in pos_data['amino_acid']]
    
    bars = ax.bar(x_pos, pos_data['enrichment'], 
                  color=colors, edgecolor='black', linewidth=2, alpha=0.85,
                  width=0.8)
    
    # Add significance stars
    for j, (idx, row) in enumerate(pos_data.iterrows()):
        # Stars above bars
        y_pos = row['enrichment'] + 0.1
        if row['significance'] != 'ns':
            ax.text(x_pos[j], y_pos, row['significance'],
                   ha='center', va='bottom', fontsize=16,
                   fontweight='bold', color='darkred')
        
        # Amino acid labels on bars
        ax.text(x_pos[j], row['enrichment']/2, row['amino_acid'],
               ha='center', va='center', fontsize=14,
               fontweight='bold', color='white',
               bbox=dict(boxstyle='circle', facecolor='black', alpha=0.6))
        
        # Enrichment value below bar
        ax.text(x_pos[j], -0.15, f'{row["enrichment"]:.2f}x',
               ha='center', va='top', fontsize=10,
               fontweight='bold', color='black')
    
    # Position label
    mid_x = x_pos[len(pos_data)//2]
    x_labels.append(f'Position {pos}')
    x_positions.append(mid_x)
    
    x_offset += len(pos_data) + 2

# Reference line
ax.axhline(y=1, color='gray', linestyle='--', linewidth=2, alpha=0.7, 
          label='No enrichment (ratio=1)', zorder=0)

# Styling
ax.set_ylabel('Enrichment Ratio (Hits/Screen)', fontsize=16, fontweight='bold')
ax.set_xlabel('Position in Sequence', fontsize=16, fontweight='bold')
ax.set_title('UBR3 Recognition Motif: First 5 Residues Enrichment\nTop 5 Amino Acids per Position (*** p<0.001, ** p<0.01, * p<0.05, ns=not significant)',
            fontsize=18, fontweight='bold', pad=20)

ax.set_xticks(x_positions)
ax.set_xticklabels(x_labels, fontsize=14, fontweight='bold')
ax.set_ylim(-0.5, max(df_enrichment['enrichment']) * 1.25)
ax.legend(loc='upper right', fontsize=12)
ax.grid(axis='y', alpha=0.3, linestyle=':', linewidth=1)

# Add legend for amino acid colors
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#BEBADA', edgecolor='black', label='Acidic (D/E)'),
    Patch(facecolor='#FB8072', edgecolor='black', label='Basic (K/R/H)'),
    Patch(facecolor='#8DD3C7', edgecolor='black', label='Hydrophobic'),
    Patch(facecolor='#80B1D3', edgecolor='black', label='Aromatic (F/W/Y)'),
    Patch(facecolor='#B3DE69', edgecolor='black', label='Polar (S/T/N/Q/C)'),
    Patch(facecolor='#FDB462', edgecolor='black', label='Special (P/G)')
]
ax.legend(handles=legend_elements, loc='upper left', fontsize=11, ncol=2,
         title='Amino Acid Properties', title_fontsize=12)

plt.tight_layout()
plt.savefig('logo_results_enrichment/first_five_residues_enrichment_publication.png',
           dpi=300, bbox_inches='tight')
print("[OK] Saved publication plot: first_five_residues_enrichment_publication.png")
plt.close()

# Generate summary statistics
print("\n" + "="*80)
print("SUMMARY STATISTICS")
print("="*80)

for pos in positions:
    pos_data = df_enrichment[df_enrichment['position'] == pos]
    sig_count = len(pos_data[pos_data['p_value'] < 0.05])
    
    print(f"\nPosition {pos}:")
    print(f"  Significant enrichments (p<0.05): {sig_count}/5")
    
    top = pos_data.nlargest(1, 'enrichment').iloc[0]
    print(f"  Top enriched: {top['amino_acid']} ({top['enrichment']:.2f}x, p={top['p_value']:.2e})")

# Overall statistics
print("\n" + "-"*80)
highly_sig = len(df_enrichment[df_enrichment['p_value'] < 0.001])
sig = len(df_enrichment[df_enrichment['p_value'] < 0.05])
total = len(df_enrichment)

print(f"\nOverall significance:")
print(f"  Highly significant (p<0.001): {highly_sig}/{total} ({highly_sig/total*100:.1f}%)")
print(f"  Significant (p<0.05): {sig}/{total} ({sig/total*100:.1f}%)")

print("\n" + "="*80)
print("ANALYSIS COMPLETE!")
print("="*80)
print("\nGenerated files:")
print("  * first_five_residues_enrichment.csv - Complete data with statistics")
print("  * first_five_residues_enrichment_comprehensive.png - Multi-panel analysis")
print("  * first_five_residues_enrichment_publication.png - Clean bar graph")
print("="*80 + "\n")
