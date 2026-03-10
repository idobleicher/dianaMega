#!/usr/bin/env python3
"""
Create graphs matching the exact style from the N-terminal motif enrichment paper.
Style: Clean bar graphs and heatmaps similar to ATE1 substrate analysis.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
import sys

print("\n" + "="*80)
print("CREATING GRAPHS - PAPER EXACT STYLE")
print("N-terminal motif enrichment (UBR3 substrates/Proteome)")
print("="*80 + "\n")

# Load PWMs
pwm_hits = pd.read_csv('logo_results_hits/position_weight_matrix.csv', index_col=0)
pwm_screen = pd.read_csv('logo_results_full_screen/position_weight_matrix.csv', index_col=0)

n_hits = 91
n_screen = 16514

positions = [2, 3, 4, 5, 6]

# Calculate enrichment for all amino acids at each position
all_data = []

for pos in positions:
    hits_freq = pwm_hits.loc[pos]
    screen_freq = pwm_screen.loc[pos]
    
    hits_counts = (hits_freq * n_hits).round().astype(int)
    screen_counts = (screen_freq * n_screen).round().astype(int)
    
    for aa in hits_freq.index:
        hits_count = hits_counts[aa]
        screen_count = screen_counts[aa]
        
        # Enrichment
        pseudocount = 0.001
        enrichment = (hits_freq[aa] + pseudocount) / (screen_freq[aa] + pseudocount)
        log2_enr = np.log2(enrichment)
        
        # Fisher's exact test
        contingency = [
            [hits_count, n_hits - hits_count],
            [screen_count, n_screen - screen_count]
        ]
        
        try:
            oddsratio, p_value = fisher_exact(contingency)
        except:
            p_value = 1.0
        
        # Significance
        if p_value < 0.001:
            sig = '***'
        elif p_value < 0.01:
            sig = '**'
        elif p_value < 0.05:
            sig = '*'
        else:
            sig = ''
        
        all_data.append({
            'position': pos,
            'position_label': f'+{pos-1}',  # Position labeling like in paper (+1, +2, etc.)
            'amino_acid': aa,
            'motif': f'P{pos}{aa}',  # Like "MCH", "MCF", etc.
            'hits_freq': hits_freq[aa],
            'screen_freq': screen_freq[aa],
            'enrichment': enrichment,
            'log2_enrichment': log2_enr,
            'p_value': p_value,
            'significance': sig
        })

df = pd.DataFrame(all_data)

# ============================================================================
# FIGURE 1: Bar graph like Panel B - showing top enriched motifs
# ============================================================================

print("Creating Figure 1: Bar graph (like Panel B)...")

# Get top 20 enriched amino acids across all positions
top_enriched = df.nlargest(20, 'enrichment').copy()

# Create motif labels with position
top_enriched['motif_label'] = top_enriched['amino_acid'] + ' (P' + top_enriched['position'].astype(str) + ')'

fig, ax = plt.subplots(figsize=(14, 6))

# Create bar positions
x_pos = np.arange(len(top_enriched))

# Colors based on significance
colors = []
for sig in top_enriched['significance']:
    if sig == '***':
        colors.append('#4472C4')  # Dark blue
    elif sig == '**':
        colors.append('#5B9BD5')  # Medium blue
    elif sig == '*':
        colors.append('#9DC3E6')  # Light blue
    else:
        colors.append('#DEEBF7')  # Very light blue/gray

# Create bars
bars = ax.bar(x_pos, top_enriched['enrichment'], color=colors, 
              edgecolor='black', linewidth=0.5, width=0.8)

# Reference line at y=1 (no enrichment baseline) - like in paper
ax.axhline(y=1, color='black', linestyle='--', linewidth=2, alpha=0.7, zorder=0)

# Add significance markers above bars
for i, (idx, row) in enumerate(top_enriched.iterrows()):
    if row['significance']:
        y_pos = row['enrichment'] + 0.15
        ax.text(i, y_pos, row['significance'], ha='center', va='bottom',
               fontsize=12, fontweight='bold')

# Styling - clean and minimal like the paper
ax.set_ylabel('Fold change', fontsize=14, fontweight='bold')
ax.set_xlabel('N-terminal motif', fontsize=14, fontweight='bold')
ax.set_title('N-terminal motif enrichment\n(UBR3 substrates/Proteome)', 
            fontsize=14, fontweight='bold', pad=15)

# X-axis labels
ax.set_xticks(x_pos)
ax.set_xticklabels(top_enriched['motif_label'], rotation=45, ha='right', fontsize=11)

# Y-axis
ax.set_ylim(0, max(top_enriched['enrichment']) * 1.15)
ax.tick_params(axis='both', which='major', labelsize=11)

# Grid - subtle like in paper
ax.grid(axis='y', alpha=0.3, linestyle='-', linewidth=0.5)
ax.set_axisbelow(True)

# Remove top and right spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('logo_results_enrichment/paper_style_bar_graph.png', dpi=300, bbox_inches='tight')
sys.stdout.buffer.write(b"[OK] Saved: paper_style_bar_graph.png\n")
sys.stdout.flush()
plt.close()

# ============================================================================
# FIGURE 2: Heatmap like Panel C - Log2 enrichment scores
# ============================================================================

print("Creating Figure 2: Heatmap (like Panel C)...")

# Create matrix for heatmap: positions x amino acids
pivot_data = df.pivot(index='amino_acid', columns='position_label', values='log2_enrichment')

# Sort by overall enrichment
row_order = df.groupby('amino_acid')['log2_enrichment'].mean().sort_values(ascending=False).index
pivot_data = pivot_data.loc[row_order]

# Reorder columns
pivot_data = pivot_data[sorted(pivot_data.columns, key=lambda x: int(x[1:]))]

fig, ax = plt.subplots(figsize=(8, 10))

# Create heatmap with similar color scheme to paper
# Using a diverging colormap centered at 0
cmap = sns.diverging_palette(220, 10, as_cmap=True, center='light')

# Plot heatmap
sns.heatmap(pivot_data, cmap=cmap, center=0, 
            vmin=-3, vmax=3,
            annot=True, fmt='.1f', 
            linewidths=0.5, linecolor='gray',
            cbar_kws={'label': 'Log₂ enrichment score', 'shrink': 0.8},
            ax=ax, annot_kws={'fontsize': 9})

# Styling
ax.set_xlabel('N-terminal position', fontsize=12, fontweight='bold')
ax.set_ylabel('Amino acid', fontsize=12, fontweight='bold')
ax.set_title('Motif Log₂ enrichment score\n(UBR3 Positives/UBR3 Negatives)', 
            fontsize=13, fontweight='bold', pad=15)

# Rotate labels
ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize=11)
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=10)

plt.tight_layout()
plt.savefig('logo_results_enrichment/paper_style_heatmap.png', dpi=300, bbox_inches='tight')
sys.stdout.buffer.write(b"[OK] Saved: paper_style_heatmap.png\n")
sys.stdout.flush()
plt.close()

# ============================================================================
# FIGURE 3: Position-specific bar graphs (5 panels like Figure B)
# ============================================================================

print("Creating Figure 3: Position-specific panels...")

fig, axes = plt.subplots(1, 5, figsize=(18, 5), sharey=True)

for i, (ax, pos) in enumerate(zip(axes, positions)):
    pos_data = df[df['position'] == pos].nlargest(5, 'enrichment')
    
    x_pos = np.arange(len(pos_data))
    
    # Colors based on significance
    colors = []
    for sig in pos_data['significance']:
        if sig == '***':
            colors.append('#4472C4')
        elif sig == '**':
            colors.append('#5B9BD5')
        elif sig == '*':
            colors.append('#9DC3E6')
        else:
            colors.append('#DEEBF7')
    
    bars = ax.bar(x_pos, pos_data['enrichment'], color=colors,
                  edgecolor='black', linewidth=0.5, width=0.7)
    
    # Reference line at y=1 (no enrichment baseline)
    ax.axhline(y=1, color='black', linestyle='--', linewidth=2, alpha=0.7, zorder=0)
    
    # Add significance markers
    for j, (idx, row) in enumerate(pos_data.iterrows()):
        if row['significance']:
            y_pos = row['enrichment'] + 0.15
            ax.text(j, y_pos, row['significance'], ha='center', va='bottom',
                   fontsize=11, fontweight='bold')
    
    # Styling
    ax.set_title(f'Position {pos}\n(+{pos-1})', fontsize=12, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(pos_data['amino_acid'], fontsize=11, fontweight='bold')
    ax.grid(axis='y', alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    if i == 0:
        ax.set_ylabel('Fold change', fontsize=13, fontweight='bold')
    
    ax.set_xlabel('Amino acid', fontsize=11)

plt.suptitle('N-terminal motif enrichment by position (UBR3 substrates/Proteome)',
            fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('logo_results_enrichment/paper_style_position_panels.png', dpi=300, bbox_inches='tight')
sys.stdout.buffer.write(b"[OK] Saved: paper_style_position_panels.png\n")
sys.stdout.flush()
plt.close()

# ============================================================================
# FIGURE 4: Simplified heatmap (fewer amino acids, like panel C style)
# ============================================================================

print("Creating Figure 4: Simplified heatmap...")

# Get top 15 amino acids by overall enrichment
top_aas = df.groupby('amino_acid')['log2_enrichment'].mean().nlargest(15).index
df_top = df[df['amino_acid'].isin(top_aas)]

# Create pivot
pivot_simple = df_top.pivot(index='amino_acid', columns='position_label', values='log2_enrichment')
pivot_simple = pivot_simple.loc[top_aas]
pivot_simple = pivot_simple[sorted(pivot_simple.columns, key=lambda x: int(x[1:]))]

fig, ax = plt.subplots(figsize=(7, 8))

# Heatmap with annotations showing enrichment values
sns.heatmap(pivot_simple, cmap=cmap, center=0,
            vmin=-3, vmax=3,
            annot=True, fmt='.1f',
            linewidths=1, linecolor='white',
            cbar_kws={'label': '', 'shrink': 0.8},
            ax=ax, annot_kws={'fontsize': 10, 'fontweight': 'bold'})

# Styling
ax.set_xlabel('N-terminal position', fontsize=12, fontweight='bold')
ax.set_ylabel('', fontsize=12, fontweight='bold')
ax.set_title('Motif Log₂ enrichment score (Top 15 amino acids)\n(UBR3 substrates/Proteome)',
            fontsize=12, fontweight='bold', pad=15)

ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize=11, fontweight='bold')
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=11, fontweight='bold')

# Add text annotations for highly significant entries
for i, aa in enumerate(pivot_simple.index):
    for j, pos_label in enumerate(pivot_simple.columns):
        pos = int(pos_label[1:]) + 1
        sig = df[(df['amino_acid'] == aa) & (df['position'] == pos)]['significance'].values
        if len(sig) > 0 and sig[0]:
            ax.text(j + 0.5, i + 0.8, sig[0], ha='center', va='bottom',
                   fontsize=8, fontweight='bold', color='darkred')

plt.tight_layout()
plt.savefig('logo_results_enrichment/paper_style_heatmap_top15.png', dpi=300, bbox_inches='tight')
sys.stdout.buffer.write(b"[OK] Saved: paper_style_heatmap_top15.png\n")
sys.stdout.flush()
plt.close()

# ============================================================================
# Print Summary
# ============================================================================

print("\n" + "="*80)
print("SUMMARY")
print("="*80)

print("\nTop 10 most enriched motifs:")
print("-" * 60)
for i, (idx, row) in enumerate(df.nlargest(10, 'enrichment').iterrows(), 1):
    sig_label = f"({row['significance']})" if row['significance'] else "(ns)"
    print(f"{i:2d}. {row['amino_acid']} at Position {row['position']} (+{row['position']-1}): "
          f"{row['enrichment']:.2f}× {sig_label:6s} p={row['p_value']:.2e}")

print("\n" + "="*80)
print("FILES GENERATED (Paper-exact style)")
print("="*80)
print("\n1. paper_style_bar_graph.png")
print("   - Bar graph showing top 20 enriched motifs (like Panel B)")
print("\n2. paper_style_heatmap.png")
print("   - Full heatmap with log2 enrichment scores (like Panel C)")
print("\n3. paper_style_position_panels.png")
print("   - 5-panel view showing enrichment by position")
print("\n4. paper_style_heatmap_top15.png")
print("   - Simplified heatmap with top 15 amino acids")
print("\n" + "="*80 + "\n")
