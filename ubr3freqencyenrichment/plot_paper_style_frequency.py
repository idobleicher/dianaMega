#!/usr/bin/env python3
"""
Create graphs in paper style showing FREQUENCY instead of enrichment.
Shows amino acid composition at each position in UBR3 substrates vs Proteome.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
import sys

print("\n" + "="*80)
print("CREATING FREQUENCY GRAPHS - PAPER STYLE")
print("Amino acid frequency at N-terminal positions")
print("="*80 + "\n")

# Load PWMs
pwm_hits = pd.read_csv('logo_results_hits/position_weight_matrix.csv', index_col=0)
pwm_screen = pd.read_csv('logo_results_full_screen/position_weight_matrix.csv', index_col=0)

n_hits = 91
n_screen = 16514

positions = [2, 3, 4, 5, 6]

# Calculate data for all amino acids at each position
all_data = []

for pos in positions:
    hits_freq = pwm_hits.loc[pos]
    screen_freq = pwm_screen.loc[pos]
    
    hits_counts = (hits_freq * n_hits).round().astype(int)
    screen_counts = (screen_freq * n_screen).round().astype(int)
    
    for aa in hits_freq.index:
        hits_count = hits_counts[aa]
        screen_count = screen_counts[aa]
        
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
            'position_label': f'+{pos-1}',
            'amino_acid': aa,
            'hits_freq': hits_freq[aa] * 100,  # Convert to percentage
            'screen_freq': screen_freq[aa] * 100,  # Convert to percentage
            'hits_count': hits_count,
            'screen_count': screen_count,
            'p_value': p_value,
            'significance': sig
        })

df = pd.DataFrame(all_data)

# ============================================================================
# FIGURE 1: Bar graph showing frequency in HITS
# ============================================================================

print("Creating Figure 1: Frequency bar graph (Hits)...")

# Get top amino acids by frequency in hits
top_hits = df.nlargest(20, 'hits_freq').copy()
top_hits['motif_label'] = top_hits['amino_acid'] + ' (P' + top_hits['position'].astype(str) + ')'

fig, ax = plt.subplots(figsize=(14, 6))

x_pos = np.arange(len(top_hits))

# Colors based on significance
colors = []
for sig in top_hits['significance']:
    if sig == '***':
        colors.append('#4472C4')
    elif sig == '**':
        colors.append('#5B9BD5')
    elif sig == '*':
        colors.append('#9DC3E6')
    else:
        colors.append('#DEEBF7')

bars = ax.bar(x_pos, top_hits['hits_freq'], color=colors,
              edgecolor='black', linewidth=0.5, width=0.8)

# Add significance markers
for i, (idx, row) in enumerate(top_hits.iterrows()):
    if row['significance']:
        y_pos = row['hits_freq'] + 1
        ax.text(i, y_pos, row['significance'], ha='center', va='bottom',
               fontsize=12, fontweight='bold')

# Styling
ax.set_ylabel('Frequency (%)', fontsize=14, fontweight='bold')
ax.set_xlabel('N-terminal motif', fontsize=14, fontweight='bold')
ax.set_title('N-terminal amino acid frequency\n(UBR3 substrates)', 
            fontsize=14, fontweight='bold', pad=15)

ax.set_xticks(x_pos)
ax.set_xticklabels(top_hits['motif_label'], rotation=45, ha='right', fontsize=11)
ax.set_ylim(0, max(top_hits['hits_freq']) * 1.15)
ax.tick_params(axis='both', which='major', labelsize=11)
ax.grid(axis='y', alpha=0.3, linestyle='-', linewidth=0.5)
ax.set_axisbelow(True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('logo_results_enrichment/paper_style_frequency_hits.png', dpi=300, bbox_inches='tight')
sys.stdout.buffer.write(b"[OK] Saved: paper_style_frequency_hits.png\n")
sys.stdout.flush()
plt.close()

# ============================================================================
# FIGURE 2: Side-by-side comparison (Hits vs Library)
# ============================================================================

print("Creating Figure 2: Side-by-side frequency comparison...")

# Get top 15 motifs by hits frequency
top_motifs = df.nlargest(15, 'hits_freq').copy()
top_motifs['motif_label'] = top_motifs['amino_acid'] + '\n(P' + top_motifs['position'].astype(str) + ')'

fig, ax = plt.subplots(figsize=(15, 7))

x_pos = np.arange(len(top_motifs))
width = 0.35

# Hits bars
bars1 = ax.bar(x_pos - width/2, top_motifs['hits_freq'], width,
               label='UBR3 substrates', color='#4472C4', edgecolor='black', linewidth=0.5)

# Library bars
bars2 = ax.bar(x_pos + width/2, top_motifs['screen_freq'], width,
               label='Library/Proteome', color='#DEEBF7', edgecolor='black', linewidth=0.5)

# Add significance markers
for i, (idx, row) in enumerate(top_motifs.iterrows()):
    if row['significance']:
        y_pos = max(row['hits_freq'], row['screen_freq']) + 1.5
        ax.text(i, y_pos, row['significance'], ha='center', va='bottom',
               fontsize=11, fontweight='bold', color='darkred')

# Styling
ax.set_ylabel('Frequency (%)', fontsize=14, fontweight='bold')
ax.set_xlabel('N-terminal motif', fontsize=14, fontweight='bold')
ax.set_title('N-terminal amino acid frequency comparison\n(UBR3 substrates vs Library)',
            fontsize=14, fontweight='bold', pad=15)

ax.set_xticks(x_pos)
ax.set_xticklabels(top_motifs['motif_label'], fontsize=10)
ax.legend(loc='upper right', fontsize=12, frameon=True, edgecolor='black')
ax.grid(axis='y', alpha=0.3, linestyle='-', linewidth=0.5)
ax.set_axisbelow(True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('logo_results_enrichment/paper_style_frequency_comparison.png', dpi=300, bbox_inches='tight')
sys.stdout.buffer.write(b"[OK] Saved: paper_style_frequency_comparison.png\n")
sys.stdout.flush()
plt.close()

# ============================================================================
# FIGURE 3: Position-specific frequency panels (5 panels)
# ============================================================================

print("Creating Figure 3: Position-specific frequency panels...")

fig, axes = plt.subplots(1, 5, figsize=(18, 5), sharey=True)

for i, (ax, pos) in enumerate(zip(axes, positions)):
    pos_data = df[df['position'] == pos].nlargest(5, 'hits_freq')
    
    x_pos = np.arange(len(pos_data))
    width = 0.35
    
    # Hits bars
    bars1 = ax.bar(x_pos - width/2, pos_data['hits_freq'], width,
                   label='Hits' if i == 0 else '', color='#4472C4',
                   edgecolor='black', linewidth=0.5)
    
    # Library bars
    bars2 = ax.bar(x_pos + width/2, pos_data['screen_freq'], width,
                   label='Library' if i == 0 else '', color='#DEEBF7',
                   edgecolor='black', linewidth=0.5)
    
    # Add significance markers
    for j, (idx, row) in enumerate(pos_data.iterrows()):
        if row['significance']:
            y_pos = max(row['hits_freq'], row['screen_freq']) + 1.5
            ax.text(j, y_pos, row['significance'], ha='center', va='bottom',
                   fontsize=10, fontweight='bold')
    
    # Styling
    ax.set_title(f'Position {pos}\n(+{pos-1})', fontsize=12, fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(pos_data['amino_acid'], fontsize=11, fontweight='bold')
    ax.grid(axis='y', alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    if i == 0:
        ax.set_ylabel('Frequency (%)', fontsize=13, fontweight='bold')
        ax.legend(loc='upper right', fontsize=10)
    
    ax.set_xlabel('Amino acid', fontsize=11)

plt.suptitle('Amino acid frequency by position (Hits vs Library)',
            fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('logo_results_enrichment/paper_style_frequency_panels.png', dpi=300, bbox_inches='tight')
sys.stdout.buffer.write(b"[OK] Saved: paper_style_frequency_panels.png\n")
sys.stdout.flush()
plt.close()

# ============================================================================
# FIGURE 4: Heatmap - Frequency in Hits
# ============================================================================

print("Creating Figure 4: Frequency heatmap (Hits)...")

# Get top 15 amino acids by average frequency
top_aas = df.groupby('amino_acid')['hits_freq'].mean().nlargest(15).index
df_top = df[df['amino_acid'].isin(top_aas)]

# Create pivot
pivot_hits = df_top.pivot(index='amino_acid', columns='position_label', values='hits_freq')
pivot_hits = pivot_hits.loc[top_aas]
pivot_hits = pivot_hits[sorted(pivot_hits.columns, key=lambda x: int(x[1:]))]

fig, ax = plt.subplots(figsize=(7, 8))

# Heatmap with frequency values
sns.heatmap(pivot_hits, cmap='YlOrRd', vmin=0, vmax=20,
            annot=True, fmt='.1f',
            linewidths=1, linecolor='white',
            cbar_kws={'label': 'Frequency (%)', 'shrink': 0.8},
            ax=ax, annot_kws={'fontsize': 10, 'fontweight': 'bold'})

# Styling
ax.set_xlabel('N-terminal position', fontsize=12, fontweight='bold')
ax.set_ylabel('', fontsize=12, fontweight='bold')
ax.set_title('Amino acid frequency (%) - UBR3 substrates\n(Top 15 amino acids)',
            fontsize=12, fontweight='bold', pad=15)

ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize=11, fontweight='bold')
ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=11, fontweight='bold')

# Add significance markers
for i, aa in enumerate(pivot_hits.index):
    for j, pos_label in enumerate(pivot_hits.columns):
        pos = int(pos_label[1:]) + 1
        sig = df[(df['amino_acid'] == aa) & (df['position'] == pos)]['significance'].values
        if len(sig) > 0 and sig[0]:
            ax.text(j + 0.5, i + 0.8, sig[0], ha='center', va='bottom',
                   fontsize=8, fontweight='bold', color='darkred')

plt.tight_layout()
plt.savefig('logo_results_enrichment/paper_style_frequency_heatmap.png', dpi=300, bbox_inches='tight')
sys.stdout.buffer.write(b"[OK] Saved: paper_style_frequency_heatmap.png\n")
sys.stdout.flush()
plt.close()

# ============================================================================
# FIGURE 5: Dual heatmaps (Hits and Library side by side)
# ============================================================================

print("Creating Figure 5: Dual heatmaps...")

# Get top 12 amino acids
top_aas = df.groupby('amino_acid')['hits_freq'].mean().nlargest(12).index
df_top = df[df['amino_acid'].isin(top_aas)]

# Create pivots
pivot_hits = df_top.pivot(index='amino_acid', columns='position_label', values='hits_freq')
pivot_screen = df_top.pivot(index='amino_acid', columns='position_label', values='screen_freq')

pivot_hits = pivot_hits.loc[top_aas]
pivot_screen = pivot_screen.loc[top_aas]

pivot_hits = pivot_hits[sorted(pivot_hits.columns, key=lambda x: int(x[1:]))]
pivot_screen = pivot_screen[sorted(pivot_screen.columns, key=lambda x: int(x[1:]))]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 8))

# Hits heatmap
sns.heatmap(pivot_hits, cmap='Blues', vmin=0, vmax=20,
            annot=True, fmt='.1f',
            linewidths=0.5, linecolor='gray',
            cbar_kws={'label': 'Frequency (%)'},
            ax=ax1, annot_kws={'fontsize': 9})

ax1.set_title('UBR3 substrates', fontsize=13, fontweight='bold', pad=10)
ax1.set_xlabel('N-terminal position', fontsize=11, fontweight='bold')
ax1.set_ylabel('Amino acid', fontsize=11, fontweight='bold')
ax1.set_xticklabels(ax1.get_xticklabels(), rotation=0, fontsize=10)
ax1.set_yticklabels(ax1.get_yticklabels(), rotation=0, fontsize=10)

# Library heatmap
sns.heatmap(pivot_screen, cmap='Greens', vmin=0, vmax=20,
            annot=True, fmt='.1f',
            linewidths=0.5, linecolor='gray',
            cbar_kws={'label': 'Frequency (%)'},
            ax=ax2, annot_kws={'fontsize': 9})

ax2.set_title('Library/Proteome', fontsize=13, fontweight='bold', pad=10)
ax2.set_xlabel('N-terminal position', fontsize=11, fontweight='bold')
ax2.set_ylabel('', fontsize=11, fontweight='bold')
ax2.set_xticklabels(ax2.get_xticklabels(), rotation=0, fontsize=10)
ax2.set_yticklabels(ax2.get_yticklabels(), rotation=0, fontsize=10)

plt.suptitle('Amino acid frequency comparison (Top 12 amino acids)',
            fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout()
plt.savefig('logo_results_enrichment/paper_style_frequency_dual_heatmap.png', dpi=300, bbox_inches='tight')
sys.stdout.buffer.write(b"[OK] Saved: paper_style_frequency_dual_heatmap.png\n")
sys.stdout.flush()
plt.close()

# ============================================================================
# Print Summary
# ============================================================================

print("\n" + "="*80)
print("SUMMARY")
print("="*80)

print("\nTop 10 most frequent amino acids in UBR3 substrates:")
print("-" * 60)
for i, (idx, row) in enumerate(df.nlargest(10, 'hits_freq').iterrows(), 1):
    sig_label = f"({row['significance']})" if row['significance'] else "(ns)"
    print(f"{i:2d}. {row['amino_acid']} at Position {row['position']} (+{row['position']-1}): "
          f"{row['hits_freq']:.1f}% in hits vs {row['screen_freq']:.1f}% in library "
          f"{sig_label:6s} p={row['p_value']:.2e}")

print("\n" + "="*80)
print("FILES GENERATED (Frequency-based, Paper style)")
print("="*80)
print("\n1. paper_style_frequency_hits.png")
print("   - Bar graph showing top 20 most frequent amino acids in hits")
print("\n2. paper_style_frequency_comparison.png")
print("   - Side-by-side comparison: Hits vs Library")
print("\n3. paper_style_frequency_panels.png")
print("   - 5-panel view showing frequency by position")
print("\n4. paper_style_frequency_heatmap.png")
print("   - Heatmap showing frequency in hits (top 15 amino acids)")
print("\n5. paper_style_frequency_dual_heatmap.png")
print("   - Dual heatmaps: Hits and Library side-by-side")
print("\n" + "="*80 + "\n")
