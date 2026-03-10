#!/usr/bin/env python3
"""
Create a single combined bar graph showing enrichment for all 5 first residues
(positions 2-6) in one figure - paper style
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

print("\n" + "="*80)
print("CREATING COMBINED BAR GRAPH - ALL 5 POSITIONS")
print("="*80 + "\n")

# Load PWMs
pwm_hits = pd.read_csv('logo_results_hits/position_weight_matrix.csv', index_col=0)
pwm_screen = pd.read_csv('logo_results_full_screen/position_weight_matrix.csv', index_col=0)

n_hits = 91
n_screen = 16514

positions = [2, 3, 4, 5, 6]

# Color scheme by amino acid properties
aa_colors = {
    # Hydrophobic - Blue
    'A': '#4BACC6', 'V': '#4BACC6', 'I': '#4BACC6', 'L': '#4BACC6', 'M': '#4BACC6',
    # Aromatic - Green
    'F': '#9BBB59', 'W': '#9BBB59', 'Y': '#9BBB59',
    # Positively charged - Orange
    'K': '#F79646', 'R': '#F79646', 'H': '#F79646',
    # Negatively charged - Red
    'D': '#C0504D', 'E': '#C0504D',
    # Polar uncharged - Purple
    'S': '#8064A2', 'T': '#8064A2', 'N': '#8064A2', 'Q': '#8064A2', 'C': '#8064A2',
    # Special - Yellow
    'G': '#F2C80F', 'P': '#F2C80F'
}

# Calculate enrichment for all positions
all_data = []

for pos in positions:
    hits_freq = pwm_hits.loc[pos]
    screen_freq = pwm_screen.loc[pos]
    
    hits_counts = (hits_freq * n_hits).round().astype(int)
    screen_counts = (screen_freq * n_screen).round().astype(int)
    
    # Calculate enrichment for ALL amino acids
    pseudocount = 0.001
    enrichments = {}
    for aa in hits_freq.index:
        enrichment = (hits_freq[aa] + pseudocount) / (screen_freq[aa] + pseudocount)
        enrichments[aa] = enrichment
    
    # Get top 5 enriched amino acids by ENRICHMENT (not frequency!)
    enrichment_series = pd.Series(enrichments)
    top_aas = enrichment_series.nlargest(5).index.tolist()
    
    for aa in top_aas:
        hits_count = hits_counts[aa]
        screen_count = screen_counts[aa]
        
        # Use pre-calculated enrichment
        enrichment = enrichments[aa]
        
        # Fisher's exact test
        contingency = [
            [hits_count, n_hits - hits_count],
            [screen_count, n_screen - screen_count]
        ]
        
        try:
            oddsratio, p_value = fisher_exact(contingency)
        except:
            p_value = 1.0
        
        all_data.append({
            'position': pos,
            'amino_acid': aa,
            'enrichment': enrichment,
            'p_value': p_value,
            'hits_freq': hits_freq[aa],
            'screen_freq': screen_freq[aa]
        })

df = pd.DataFrame(all_data)

# Create the combined figure
fig, ax = plt.subplots(figsize=(18, 8))

# Create grouped bars
bar_width = 0.15
positions_array = np.arange(len(positions))

# For each position, plot top 5 amino acids
x_offset = 0
all_bars = []
x_labels = []
x_ticks = []

for i, pos in enumerate(positions):
    pos_data = df[df['position'] == pos].nlargest(5, 'enrichment')
    
    x_positions = np.arange(len(pos_data)) * bar_width + x_offset
    colors = [aa_colors.get(aa, '#CCCCCC') for aa in pos_data['amino_acid']]
    
    bars = ax.bar(x_positions, pos_data['enrichment'],
                  width=bar_width * 0.9,
                  color=colors, 
                  edgecolor='black', 
                  linewidth=1.5,
                  alpha=0.85)
    
    all_bars.extend(bars)
    
    # Add significance stars above bars
    for j, (_, row) in enumerate(pos_data.iterrows()):
        if row['p_value'] < 0.001:
            sig = '***'
            y_offset = 0.2
        elif row['p_value'] < 0.01:
            sig = '**'
            y_offset = 0.15
        elif row['p_value'] < 0.05:
            sig = '*'
            y_offset = 0.12
        else:
            sig = ''
            y_offset = 0
        
        if sig:
            ax.text(x_positions[j], row['enrichment'] + y_offset, sig,
                   ha='center', va='bottom', fontsize=14,
                   fontweight='bold', color='black')
        
        # Add amino acid labels below bars
        ax.text(x_positions[j], -0.3, row['amino_acid'],
               ha='center', va='top', fontsize=11,
               fontweight='bold', color='black')
    
    # Store position label location
    mid_x = x_positions[len(pos_data)//2]
    x_ticks.append(mid_x)
    x_labels.append(f'Position +{pos-1}\n(Pos {pos})')
    
    # Add separator line
    if i < len(positions) - 1:
        separator_x = x_offset + len(pos_data) * bar_width + bar_width/2
        ax.axvline(x=separator_x, color='gray', linestyle=':', 
                  linewidth=2, alpha=0.5)
    
    x_offset += len(pos_data) * bar_width + bar_width

# Styling
ax.set_ylabel('Enrichment ratio (Hits / Screen)', fontsize=16, fontweight='bold')
ax.set_xlabel('N-terminal position', fontsize=16, fontweight='bold')
ax.set_title('Amino acid enrichment at N-terminal positions 2-6\nTop 5 enriched amino acids per position (*** p<0.001, ** p<0.01, * p<0.05)',
            fontsize=18, fontweight='bold', pad=20)

# X-axis
ax.set_xticks(x_ticks)
ax.set_xticklabels(x_labels, fontsize=13, fontweight='bold')

# Y-axis
ax.set_ylim(-0.6, max(df['enrichment']) * 1.3)

# Reference line
ax.axhline(y=1, color='gray', linestyle='--', linewidth=2, alpha=0.7,
          label='No enrichment (ratio = 1)', zorder=0)

# Grid
ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.8)
ax.set_axisbelow(True)

# Spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)

# Legend for amino acid types
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#C0504D', edgecolor='black', label='Acidic (D/E)', linewidth=1.5),
    Patch(facecolor='#F79646', edgecolor='black', label='Basic (K/R/H)', linewidth=1.5),
    Patch(facecolor='#4BACC6', edgecolor='black', label='Hydrophobic (A/V/I/L/M)', linewidth=1.5),
    Patch(facecolor='#9BBB59', edgecolor='black', label='Aromatic (F/W/Y)', linewidth=1.5),
    Patch(facecolor='#8064A2', edgecolor='black', label='Polar (S/T/N/Q/C)', linewidth=1.5),
    Patch(facecolor='#F2C80F', edgecolor='black', label='Special (G/P)', linewidth=1.5)
]
ax.legend(handles=legend_elements, loc='upper right', ncol=2,
         fontsize=12, frameon=True, title='Amino acid properties',
         title_fontsize=13, edgecolor='black', fancybox=False)

plt.tight_layout()
plt.savefig('logo_results_enrichment/enrichment_combined_all_positions.png',
           dpi=300, bbox_inches='tight')
print("[OK] Saved: enrichment_combined_all_positions.png")
plt.close()

# Create a cleaner version with only top 3 per position
fig, ax = plt.subplots(figsize=(16, 8))

bar_width = 0.2
x_offset = 0
x_ticks = []
x_labels = []

for i, pos in enumerate(positions):
    pos_data = df[df['position'] == pos].nlargest(3, 'enrichment')  # Top 3 only
    
    x_positions = np.arange(len(pos_data)) * bar_width + x_offset
    colors = [aa_colors.get(aa, '#CCCCCC') for aa in pos_data['amino_acid']]
    
    bars = ax.bar(x_positions, pos_data['enrichment'],
                  width=bar_width * 0.85,
                  color=colors, 
                  edgecolor='black', 
                  linewidth=2,
                  alpha=0.9)
    
    # Add significance and values
    for j, (_, row) in enumerate(pos_data.iterrows()):
        # Significance stars
        if row['p_value'] < 0.001:
            sig = '***'
        elif row['p_value'] < 0.01:
            sig = '**'
        elif row['p_value'] < 0.05:
            sig = '*'
        else:
            sig = 'ns'
        
        if sig != 'ns':
            ax.text(x_positions[j], row['enrichment'] + 0.2, sig,
                   ha='center', va='bottom', fontsize=16,
                   fontweight='bold', color='darkred')
        
        # Amino acid label ON the bar
        ax.text(x_positions[j], row['enrichment']/2, row['amino_acid'],
               ha='center', va='center', fontsize=14,
               fontweight='bold', color='white',
               bbox=dict(boxstyle='circle', facecolor='black', 
                        edgecolor='white', linewidth=2, alpha=0.7, pad=0.3))
        
        # Enrichment value below bar
        ax.text(x_positions[j], -0.25, f'{row["enrichment"]:.2f}×',
               ha='center', va='top', fontsize=10,
               fontweight='bold', color='black')
    
    # Position label
    mid_x = x_positions[1]  # Middle of 3 bars
    x_ticks.append(mid_x)
    x_labels.append(f'Position +{pos-1}')
    
    # Separator
    if i < len(positions) - 1:
        separator_x = x_offset + 3 * bar_width + bar_width/2
        ax.axvline(x=separator_x, color='lightgray', linestyle='-', 
                  linewidth=3, alpha=0.6)
    
    x_offset += 3 * bar_width + bar_width * 1.5

# Styling
ax.set_ylabel('Enrichment ratio (Hits / Screen)', fontsize=16, fontweight='bold')
ax.set_xlabel('N-terminal position', fontsize=16, fontweight='bold')
ax.set_title('UBR3 N-degron motif: Top 3 enriched amino acids at positions 2-6\n(*** p<0.001, ** p<0.01, * p<0.05, ns=not significant)',
            fontsize=18, fontweight='bold', pad=20)

ax.set_xticks(x_ticks)
ax.set_xticklabels(x_labels, fontsize=14, fontweight='bold')
ax.set_ylim(-0.5, max(df['enrichment']) * 1.35)

ax.axhline(y=1, color='gray', linestyle='--', linewidth=2.5, alpha=0.7,
          label='No enrichment', zorder=0)

ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=1)
ax.set_axisbelow(True)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_linewidth(2)
ax.spines['bottom'].set_linewidth(2)

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#C0504D', edgecolor='black', label='Acidic', linewidth=1.5),
    Patch(facecolor='#F79646', edgecolor='black', label='Basic', linewidth=1.5),
    Patch(facecolor='#4BACC6', edgecolor='black', label='Hydrophobic', linewidth=1.5),
    Patch(facecolor='#9BBB59', edgecolor='black', label='Aromatic', linewidth=1.5),
    Patch(facecolor='#8064A2', edgecolor='black', label='Polar', linewidth=1.5),
    Patch(facecolor='#F2C80F', edgecolor='black', label='Special', linewidth=1.5)
]
ax.legend(handles=legend_elements, loc='upper left', ncol=3,
         fontsize=11, frameon=True, title='Amino acid type',
         title_fontsize=12, edgecolor='black')

plt.tight_layout()
plt.savefig('logo_results_enrichment/enrichment_combined_top3.png',
           dpi=300, bbox_inches='tight')
print("[OK] Saved: enrichment_combined_top3.png")
plt.close()

# Summary table
print("\n" + "="*80)
print("ENRICHMENT SUMMARY - ALL POSITIONS")
print("="*80)

for pos in positions:
    print(f"\nPosition +{pos-1} (Position {pos}):")
    pos_data = df[df['position'] == pos].nlargest(5, 'enrichment')
    for rank, (_, row) in enumerate(pos_data.iterrows(), 1):
        sig = '***' if row['p_value'] < 0.001 else '**' if row['p_value'] < 0.01 else '*' if row['p_value'] < 0.05 else 'ns'
        print(f"  {rank}. {row['amino_acid']}: {row['enrichment']:.2f}× "
              f"({row['hits_freq']:.1%} hits vs {row['screen_freq']:.1%} screen) "
              f"p={row['p_value']:.2e} {sig}")

print("\n" + "="*80)
print("FILES GENERATED")
print("="*80)
print("\n1. enrichment_combined_all_positions.png")
print("   - All top 5 amino acids per position in one graph")
print("\n2. enrichment_combined_top3.png")
print("   - Cleaner version with top 3 per position")
print("\n" + "="*80 + "\n")
