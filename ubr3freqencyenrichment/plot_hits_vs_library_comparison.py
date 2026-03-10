#!/usr/bin/env python3
"""
Create publication-quality bar graphs comparing best hits vs full library
for the first 5 residues (positions 2-6)
Multiple visualization styles for different publication needs
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
import seaborn as sns

print("\n" + "="*80)
print("PUBLICATION FIGURE: HITS vs LIBRARY COMPARISON")
print("First 5 residues (Positions 2-6)")
print("="*80 + "\n")

# Load PWMs
pwm_hits = pd.read_csv('logo_results_hits/position_weight_matrix.csv', index_col=0)
pwm_screen = pd.read_csv('logo_results_full_screen/position_weight_matrix.csv', index_col=0)

n_hits = 91
n_screen = 16514

positions = [2, 3, 4, 5, 6]

# Color scheme
aa_colors = {
    'D': '#C0504D', 'E': '#C0504D',  # Acidic - Red
    'K': '#F79646', 'R': '#F79646', 'H': '#F79646',  # Basic - Orange
    'A': '#4BACC6', 'V': '#4BACC6', 'I': '#4BACC6', 'L': '#4BACC6', 'M': '#4BACC6',  # Hydrophobic - Blue
    'F': '#9BBB59', 'W': '#9BBB59', 'Y': '#9BBB59',  # Aromatic - Green
    'S': '#8064A2', 'T': '#8064A2', 'N': '#8064A2', 'Q': '#8064A2', 'C': '#8064A2',  # Polar - Purple
    'G': '#F2C80F', 'P': '#F2C80F'  # Special - Yellow
}

# Calculate enrichment data
enrichment_data = []
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
    
    # Get top 3 enriched amino acids by ENRICHMENT (not frequency!)
    enrichment_series = pd.Series(enrichments)
    top_aas = enrichment_series.nlargest(3).index.tolist()
    
    for aa in top_aas:
        enrichment = enrichments[aa]
        
        contingency = [[hits_counts[aa], n_hits - hits_counts[aa]],
                      [screen_counts[aa], n_screen - screen_counts[aa]]]
        try:
            _, p_value = fisher_exact(contingency)
        except:
            p_value = 1.0
        
        enrichment_data.append({
            'position': pos,
            'amino_acid': aa,
            'hits_freq': hits_freq[aa],
            'screen_freq': screen_freq[aa],
            'enrichment': enrichment,
            'p_value': p_value
        })

df = pd.DataFrame(enrichment_data)

# ============================================================================
# Figure 1: Side-by-side bars (Hits vs Library)
# ============================================================================
print("Creating Figure 1: Side-by-side comparison...")

fig, axes = plt.subplots(1, 5, figsize=(20, 6), sharey=True)

for idx, (ax, pos) in enumerate(zip(axes, positions)):
    pos_data = df[df['position'] == pos]
    
    x = np.arange(len(pos_data))
    width = 0.35
    
    # Bars for hits (darker) and screen (lighter)
    colors = [aa_colors.get(aa, '#CCCCCC') for aa in pos_data['amino_acid']]
    
    bars1 = ax.bar(x - width/2, pos_data['hits_freq'] * 100, width,
                   label='Hits', color=colors, edgecolor='black',
                   linewidth=1.5, alpha=0.9)
    
    bars2 = ax.bar(x + width/2, pos_data['screen_freq'] * 100, width,
                   label='Library', color=colors, edgecolor='black',
                   linewidth=1.5, alpha=0.4, hatch='//')
    
    # Add significance stars
    for i, (_, row) in enumerate(pos_data.iterrows()):
        if row['p_value'] < 0.001:
            sig = '***'
        elif row['p_value'] < 0.01:
            sig = '**'
        elif row['p_value'] < 0.05:
            sig = '*'
        else:
            sig = ''
        
        if sig:
            max_height = max(row['hits_freq'], row['screen_freq']) * 100
            ax.text(i, max_height + 1, sig, ha='center', va='bottom',
                   fontsize=14, fontweight='bold', color='darkred')
    
    # Styling
    ax.set_xlabel(f'Position +{pos-1}', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(pos_data['amino_acid'], fontsize=13, fontweight='bold')
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    if idx == 0:
        ax.legend(loc='upper left', fontsize=11, frameon=True)

axes[0].set_ylabel('Frequency (%)', fontsize=14, fontweight='bold')

fig.suptitle('Amino acid frequency comparison: Best hits vs Full library\nFirst 5 N-terminal positions (*** p<0.001, ** p<0.01, * p<0.05)',
            fontsize=16, fontweight='bold', y=1.02)

plt.tight_layout()
plt.savefig('logo_results_enrichment/publication_hits_vs_library_sidebyside.png',
           dpi=300, bbox_inches='tight')
print("[OK] Saved: publication_hits_vs_library_sidebyside.png")
plt.close()

# ============================================================================
# Figure 2: Grouped comparison (all amino acids, all positions)
# ============================================================================
print("Creating Figure 2: Grouped comparison...")

fig, ax = plt.subplots(figsize=(18, 7))

# Create grouped bars
x_offset = 0
x_ticks = []
x_labels = []
bar_width = 0.35

for pos in positions:
    pos_data = df[df['position'] == pos]
    
    n_aas = len(pos_data)
    x_positions = np.arange(n_aas) + x_offset
    
    colors = [aa_colors.get(aa, '#CCCCCC') for aa in pos_data['amino_acid']]
    
    # Hits bars
    bars1 = ax.bar(x_positions - bar_width/2, pos_data['hits_freq'] * 100,
                   width=bar_width, color=colors, edgecolor='black',
                   linewidth=1.5, alpha=0.95, label='Hits' if pos == 2 else '')
    
    # Library bars
    bars2 = ax.bar(x_positions + bar_width/2, pos_data['screen_freq'] * 100,
                   width=bar_width, color=colors, edgecolor='black',
                   linewidth=1.5, alpha=0.4, hatch='//',
                   label='Library' if pos == 2 else '')
    
    # Add significance
    for i, (_, row) in enumerate(pos_data.iterrows()):
        if row['p_value'] < 0.001:
            sig = '***'
        elif row['p_value'] < 0.01:
            sig = '**'
        elif row['p_value'] < 0.05:
            sig = '*'
        else:
            sig = ''
        
        if sig:
            max_h = max(row['hits_freq'], row['screen_freq']) * 100
            ax.text(x_positions[i], max_h + 1.5, sig,
                   ha='center', va='bottom', fontsize=12,
                   fontweight='bold', color='darkred')
        
        # AA label below
        ax.text(x_positions[i], -2, row['amino_acid'],
               ha='center', va='top', fontsize=11, fontweight='bold')
    
    # Position label
    mid_x = x_positions[len(pos_data)//2]
    x_ticks.append(mid_x)
    x_labels.append(f'Pos +{pos-1}')
    
    # Separator
    if pos < 6:
        sep_x = x_offset + n_aas + 0.5
        ax.axvline(x=sep_x, color='gray', linestyle=':', linewidth=2, alpha=0.5)
    
    x_offset += n_aas + 1

ax.set_ylabel('Frequency (%)', fontsize=15, fontweight='bold')
ax.set_title('Amino acid enrichment in best hits vs full library (positions 2-6)\nTop 3 enriched amino acids per position',
            fontsize=17, fontweight='bold', pad=20)

ax.set_xticks(x_ticks)
ax.set_xticklabels(x_labels, fontsize=14, fontweight='bold')
ax.set_ylim(-3, max(df['hits_freq'].max(), df['screen_freq'].max()) * 110)

ax.legend(loc='upper right', fontsize=13, frameon=True, edgecolor='black')
ax.grid(axis='y', alpha=0.3, linestyle='--')
ax.set_axisbelow(True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('logo_results_enrichment/publication_hits_vs_library_grouped.png',
           dpi=300, bbox_inches='tight')
print("[OK] Saved: publication_hits_vs_library_grouped.png")
plt.close()

# ============================================================================
# Figure 3: Enrichment-focused (fold change bars)
# ============================================================================
print("Creating Figure 3: Enrichment fold-change...")

fig, ax = plt.subplots(figsize=(16, 7))

x_offset = 0
x_ticks = []
x_labels = []

for pos in positions:
    pos_data = df[df['position'] == pos]
    
    x_positions = np.arange(len(pos_data)) + x_offset
    colors = [aa_colors.get(aa, '#CCCCCC') for aa in pos_data['amino_acid']]
    
    bars = ax.bar(x_positions, pos_data['enrichment'], width=0.7,
                  color=colors, edgecolor='black', linewidth=2, alpha=0.85)
    
    # Add significance and labels
    for i, (_, row) in enumerate(pos_data.iterrows()):
        # Significance
        if row['p_value'] < 0.001:
            sig = '***'
        elif row['p_value'] < 0.01:
            sig = '**'
        elif row['p_value'] < 0.05:
            sig = '*'
        else:
            sig = 'ns'
        
        if sig != 'ns':
            ax.text(x_positions[i], row['enrichment'] + 0.15, sig,
                   ha='center', va='bottom', fontsize=14,
                   fontweight='bold', color='darkred')
        
        # AA on bar
        ax.text(x_positions[i], row['enrichment']/2, row['amino_acid'],
               ha='center', va='center', fontsize=13, fontweight='bold',
               color='white',
               bbox=dict(boxstyle='circle', facecolor='black',
                        edgecolor='white', linewidth=2, alpha=0.7))
        
        # Frequency info below
        ax.text(x_positions[i], -0.3,
               f'{row["hits_freq"]*100:.1f}% / {row["screen_freq"]*100:.1f}%',
               ha='center', va='top', fontsize=8, color='black')
    
    # Position label
    mid_x = x_positions[len(pos_data)//2]
    x_ticks.append(mid_x)
    x_labels.append(f'Position +{pos-1}')
    
    # Separator
    if pos < 6:
        sep_x = x_offset + len(pos_data) + 0.5
        ax.axvline(x=sep_x, color='lightgray', linestyle='-',
                  linewidth=3, alpha=0.6)
    
    x_offset += len(pos_data) + 1

ax.axhline(y=1, color='gray', linestyle='--', linewidth=2, alpha=0.7,
          label='No enrichment (ratio = 1)')

ax.set_ylabel('Enrichment ratio (Hits / Library)', fontsize=15, fontweight='bold')
ax.set_xlabel('N-terminal position', fontsize=15, fontweight='bold')
ax.set_title('Amino acid enrichment in best hits relative to library\n(*** p<0.001, ** p<0.01, * p<0.05)',
            fontsize=17, fontweight='bold', pad=20)

ax.set_xticks(x_ticks)
ax.set_xticklabels(x_labels, fontsize=14, fontweight='bold')
ax.set_ylim(-0.5, max(df['enrichment']) * 1.3)

# Legend for colors
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#C0504D', edgecolor='black', label='Acidic (D/E)'),
    Patch(facecolor='#F79646', edgecolor='black', label='Basic (K/R/H)'),
    Patch(facecolor='#4BACC6', edgecolor='black', label='Hydrophobic'),
    Patch(facecolor='#9BBB59', edgecolor='black', label='Aromatic (F/W/Y)'),
    Patch(facecolor='#8064A2', edgecolor='black', label='Polar (S/T/N/Q/C)'),
    Patch(facecolor='#F2C80F', edgecolor='black', label='Special (G/P)')
]
ax.legend(handles=legend_elements, loc='upper left', ncol=2,
         fontsize=11, frameon=True, title='Amino acid type')

ax.grid(axis='y', alpha=0.3, linestyle='--')
ax.set_axisbelow(True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('logo_results_enrichment/publication_hits_vs_library_enrichment.png',
           dpi=300, bbox_inches='tight')
print("[OK] Saved: publication_hits_vs_library_enrichment.png")
plt.close()

# ============================================================================
# Figure 4: Clean minimal version (main figure candidate)
# ============================================================================
print("Creating Figure 4: Clean minimal version...")

fig, axes = plt.subplots(1, 5, figsize=(18, 5), sharey=True)

for idx, (ax, pos) in enumerate(zip(axes, positions)):
    pos_data = df[df['position'] == pos]
    
    # Only show top amino acid
    top_aa = pos_data.iloc[0]
    
    x = [0, 1]
    heights = [top_aa['hits_freq'] * 100, top_aa['screen_freq'] * 100]
    colors_pair = [aa_colors.get(top_aa['amino_acid'], '#CCCCCC')] * 2
    
    bars = ax.bar(x, heights, width=0.6, color=colors_pair,
                  edgecolor='black', linewidth=2)
    
    # Set alpha individually
    bars[0].set_alpha(0.95)
    bars[1].set_alpha(0.4)
    
    # Add hatch to library bar
    bars[1].set_hatch('//')
    
    # Significance
    if top_aa['p_value'] < 0.001:
        sig = '***'
    elif top_aa['p_value'] < 0.01:
        sig = '**'
    elif top_aa['p_value'] < 0.05:
        sig = '*'
    else:
        sig = 'ns'
    
    if sig != 'ns':
        ax.text(0.5, max(heights) + 2, sig, ha='center', va='bottom',
               fontsize=16, fontweight='bold', color='darkred')
    
    # Enrichment value
    ax.text(0.5, max(heights) * 1.15, f'{top_aa["enrichment"]:.2f}×',
           ha='center', va='bottom', fontsize=12, fontweight='bold',
           bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))
    
    # Labels
    ax.set_xticks(x)
    ax.set_xticklabels(['Hits', 'Library'], fontsize=11, fontweight='bold')
    ax.set_title(f'Position +{pos-1}\n{top_aa["amino_acid"]}',
                fontsize=14, fontweight='bold', pad=10)
    
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

axes[0].set_ylabel('Frequency (%)', fontsize=14, fontweight='bold')

fig.suptitle('Most enriched amino acid at each position: Hits vs Library',
            fontsize=16, fontweight='bold', y=1.02)

plt.tight_layout()
plt.savefig('logo_results_enrichment/publication_hits_vs_library_minimal.png',
           dpi=300, bbox_inches='tight')
print("[OK] Saved: publication_hits_vs_library_minimal.png")
plt.close()

# Print summary
print("\n" + "="*80)
print("SUMMARY TABLE")
print("="*80)

for pos in positions:
    print(f"\nPosition +{pos-1} (Position {pos}):")
    print(f"{'AA':<5} {'Hits %':<10} {'Library %':<12} {'Enrichment':<12} {'P-value':<12} {'Sig'}")
    print("-" * 60)
    pos_data = df[df['position'] == pos]
    for _, row in pos_data.iterrows():
        sig = '***' if row['p_value'] < 0.001 else '**' if row['p_value'] < 0.01 else '*' if row['p_value'] < 0.05 else 'ns'
        print(f"{row['amino_acid']:<5} {row['hits_freq']*100:>8.1f}% {row['screen_freq']*100:>10.1f}% "
              f"{row['enrichment']:>10.2f}x {row['p_value']:>10.2e}  {sig}")

print("\n" + "="*80)
print("FILES GENERATED - PUBLICATION READY")
print("="*80)
print("\n1. publication_hits_vs_library_sidebyside.png")
print("   - Side-by-side comparison for each position")
print("\n2. publication_hits_vs_library_grouped.png")
print("   - All positions grouped together")
print("\n3. publication_hits_vs_library_enrichment.png")
print("   - Enrichment fold-change focused")
print("\n4. publication_hits_vs_library_minimal.png")
print("   - Clean minimal version (recommended for main figure)")
print("\n" + "="*80 + "\n")
