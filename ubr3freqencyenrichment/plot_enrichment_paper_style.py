#!/usr/bin/env python3
"""
Create bar graphs in the style of:
"Global profiling of N-terminal cysteine-dependent degradation mechanisms"
Bekturova et al., PNAS 2025

Clean, publication-quality style with grouped amino acids by properties
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact
import os

print("\n" + "="*80)
print("CREATING BAR GRAPHS - PAPER STYLE")
print("Positions 2-6 Enrichment Analysis")
print("="*80 + "\n")

# Load PWMs
pwm_hits = pd.read_csv('logo_results_hits/position_weight_matrix.csv', index_col=0)
pwm_screen = pd.read_csv('logo_results_full_screen/position_weight_matrix.csv', index_col=0)

n_hits = 91
n_screen = 16514

positions = [2, 3, 4, 5, 6]

# Define amino acid groups (like in the paper)
aa_groups = {
    'Hydrophobic': ['L', 'V', 'I', 'M', 'A'],
    'Aromatic': ['F', 'W', 'Y'],
    'Positively charged': ['K', 'R', 'H'],
    'Negatively charged': ['D', 'E'],
    'Polar uncharged': ['S', 'T', 'N', 'Q', 'C'],
    'Special': ['G', 'P']
}

# Color scheme (clean, professional)
group_colors = {
    'Hydrophobic': '#4BACC6',  # Blue
    'Aromatic': '#9BBB59',  # Green
    'Positively charged': '#F79646',  # Orange
    'Negatively charged': '#C0504D',  # Red
    'Polar uncharged': '#8064A2',  # Purple
    'Special': '#F2C80F'  # Yellow
}

# Calculate enrichment for all positions
all_enrichment = []

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
            oddsratio = enrichment
        
        # Find group
        aa_group = None
        for group, aas in aa_groups.items():
            if aa in aas:
                aa_group = group
                break
        
        all_enrichment.append({
            'position': pos,
            'amino_acid': aa,
            'group': aa_group,
            'hits_freq': hits_freq[aa],
            'screen_freq': screen_freq[aa],
            'enrichment': enrichment,
            'log2_enrichment': log2_enr,
            'p_value': p_value
        })

df = pd.DataFrame(all_enrichment)

# Figure 1: Main paper-style figure - Top enriched per position
fig, axes = plt.subplots(1, 5, figsize=(20, 5), sharey=True)

for idx, (ax, pos) in enumerate(zip(axes, positions)):
    # Get top 5 enriched amino acids at this position
    pos_data = df[df['position'] == pos].nlargest(5, 'enrichment')
    
    # Create bars
    x_positions = np.arange(len(pos_data))
    colors = [group_colors.get(g, '#CCCCCC') for g in pos_data['group']]
    
    bars = ax.bar(x_positions, pos_data['enrichment'], 
                  color=colors, edgecolor='black', linewidth=1.5,
                  alpha=0.85, width=0.7)
    
    # Add significance stars
    for i, (_, row) in enumerate(pos_data.iterrows()):
        if row['p_value'] < 0.001:
            sig = '***'
            y_offset = 0.15
        elif row['p_value'] < 0.01:
            sig = '**'
            y_offset = 0.12
        elif row['p_value'] < 0.05:
            sig = '*'
            y_offset = 0.10
        else:
            sig = ''
            y_offset = 0
        
        if sig:
            ax.text(i, row['enrichment'] + y_offset, sig,
                   ha='center', va='bottom', fontsize=14,
                   fontweight='bold', color='black')
    
    # X-axis labels (amino acids)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(pos_data['amino_acid'], fontsize=13, fontweight='bold')
    
    # Title
    ax.set_title(f'Position +{pos-1}', fontsize=14, fontweight='bold', pad=10)
    
    # Grid
    ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)
    
    # Reference line
    ax.axhline(y=1, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
    
    # Spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# Y-axis label (only on first subplot)
axes[0].set_ylabel('Enrichment ratio\n(Hits / Screen)', fontsize=14, fontweight='bold')

# Overall title
fig.suptitle('Amino acid enrichment at N-terminal positions 2-6',
            fontsize=16, fontweight='bold', y=1.02)

# Legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=color, edgecolor='black', label=group)
                  for group, color in group_colors.items()]
fig.legend(handles=legend_elements, loc='center', bbox_to_anchor=(0.5, -0.05),
          ncol=6, fontsize=11, frameon=False)

plt.tight_layout(rect=[0, 0.05, 1, 0.98])
plt.savefig('logo_results_enrichment/enrichment_paper_style_main.png',
           dpi=300, bbox_inches='tight')
print("[OK] Saved: enrichment_paper_style_main.png")
plt.close()

# Figure 2: Grouped by amino acid properties (like the paper does)
fig, axes = plt.subplots(1, 5, figsize=(20, 5), sharey=True)

for idx, (ax, pos) in enumerate(zip(axes, positions)):
    pos_data = df[df['position'] == pos]
    
    # For each group, get the most enriched member
    group_data = []
    for group in ['Positively charged', 'Negatively charged', 'Hydrophobic', 
                  'Aromatic', 'Polar uncharged', 'Special']:
        group_members = pos_data[pos_data['group'] == group]
        if len(group_members) > 0:
            # Get the most enriched in this group
            top_member = group_members.nlargest(1, 'enrichment').iloc[0]
            group_data.append({
                'group': group,
                'amino_acid': top_member['amino_acid'],
                'enrichment': top_member['enrichment'],
                'p_value': top_member['p_value']
            })
    
    group_df = pd.DataFrame(group_data)
    
    # Create bars
    x_positions = np.arange(len(group_df))
    colors = [group_colors[g] for g in group_df['group']]
    
    bars = ax.bar(x_positions, group_df['enrichment'],
                  color=colors, edgecolor='black', linewidth=1.5,
                  alpha=0.85, width=0.7)
    
    # Add significance
    for i, row in group_df.iterrows():
        if row['p_value'] < 0.001:
            sig = '***'
        elif row['p_value'] < 0.01:
            sig = '**'
        elif row['p_value'] < 0.05:
            sig = '*'
        else:
            sig = ''
        
        if sig:
            ax.text(i, row['enrichment'] + 0.15, sig,
                   ha='center', va='bottom', fontsize=14,
                   fontweight='bold', color='black')
        
        # Add amino acid label on bar
        ax.text(i, row['enrichment']/2, row['amino_acid'],
               ha='center', va='center', fontsize=12,
               fontweight='bold', color='white')
    
    # X-axis (group names, rotated)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(group_df['group'], rotation=45, ha='right', fontsize=10)
    
    # Title
    ax.set_title(f'Position +{pos-1}', fontsize=14, fontweight='bold', pad=10)
    
    # Grid
    ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)
    
    # Reference line
    ax.axhline(y=1, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
    
    # Spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

axes[0].set_ylabel('Enrichment ratio\n(Hits / Screen)', fontsize=14, fontweight='bold')

fig.suptitle('Amino acid enrichment by biochemical property (positions 2-6)',
            fontsize=16, fontweight='bold', y=1.02)

plt.tight_layout(rect=[0, 0, 1, 0.98])
plt.savefig('logo_results_enrichment/enrichment_paper_style_grouped.png',
           dpi=300, bbox_inches='tight')
print("[OK] Saved: enrichment_paper_style_grouped.png")
plt.close()

# Figure 3: Combined view - all significant enrichments across positions
# Similar to how the paper shows key results
sig_enrichments = df[df['p_value'] < 0.05].sort_values('enrichment', ascending=False)

if len(sig_enrichments) > 0:
    # Limit to top 20 for clarity
    sig_enrichments = sig_enrichments.head(20)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Create labels
    labels = [f"Pos {row['position']}: {row['amino_acid']}" 
             for _, row in sig_enrichments.iterrows()]
    y_pos = np.arange(len(labels))
    
    colors = [group_colors.get(row['group'], '#CCCCCC') 
             for _, row in sig_enrichments.iterrows()]
    
    bars = ax.barh(y_pos, sig_enrichments['enrichment'],
                   color=colors, edgecolor='black', linewidth=1.5,
                   alpha=0.85)
    
    # Add significance markers
    for i, (_, row) in enumerate(sig_enrichments.iterrows()):
        if row['p_value'] < 0.001:
            sig = '***'
        elif row['p_value'] < 0.01:
            sig = '**'
        else:
            sig = '*'
        
        ax.text(row['enrichment'] + 0.1, i, sig,
               va='center', ha='left', fontsize=12,
               fontweight='bold', color='black')
        
        # Add enrichment value
        ax.text(row['enrichment'] - 0.05, i, f"{row['enrichment']:.2f}×",
               va='center', ha='right', fontsize=9,
               fontweight='bold', color='white')
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=11)
    ax.set_xlabel('Enrichment ratio (Hits / Screen)', fontsize=13, fontweight='bold')
    ax.set_title('Top significant enrichments (p < 0.05) across positions 2-6',
                fontsize=14, fontweight='bold', pad=15)
    
    ax.axvline(x=1, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
    ax.grid(axis='x', alpha=0.3, linestyle='--', linewidth=0.5)
    ax.set_axisbelow(True)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=color, edgecolor='black', label=group)
                      for group, color in group_colors.items()]
    ax.legend(handles=legend_elements, loc='lower right', ncol=2, 
             fontsize=10, frameon=True, title='Amino acid type',
             title_fontsize=11)
    
    plt.tight_layout()
    plt.savefig('logo_results_enrichment/enrichment_paper_style_significant.png',
               dpi=300, bbox_inches='tight')
    print("[OK] Saved: enrichment_paper_style_significant.png")
    plt.close()

# Figure 4: Simple clean version - one bar per position showing top enriched
fig, ax = plt.subplots(figsize=(10, 6))

top_per_position = []
for pos in positions:
    pos_data = df[df['position'] == pos].nlargest(1, 'enrichment').iloc[0]
    top_per_position.append(pos_data)

top_df = pd.DataFrame(top_per_position)

x_pos = np.arange(len(positions))
colors = [group_colors.get(row['group'], '#CCCCCC') for _, row in top_df.iterrows()]

bars = ax.bar(x_pos, top_df['enrichment'],
              color=colors, edgecolor='black', linewidth=2,
              alpha=0.85, width=0.6)

# Add significance
for i, (_, row) in enumerate(top_df.iterrows()):
    if row['p_value'] < 0.001:
        sig = '***'
    elif row['p_value'] < 0.01:
        sig = '**'
    elif row['p_value'] < 0.05:
        sig = '*'
    else:
        sig = 'ns'
    
    # Amino acid on bar
    ax.text(i, row['enrichment']/2, row['amino_acid'],
           ha='center', va='center', fontsize=16,
           fontweight='bold', color='white')
    
    # Significance above
    if sig != 'ns':
        ax.text(i, row['enrichment'] + 0.15, sig,
               ha='center', va='bottom', fontsize=16,
               fontweight='bold', color='black')
    
    # Enrichment value below
    ax.text(i, -0.2, f"{row['enrichment']:.2f}×",
           ha='center', va='top', fontsize=11,
           fontweight='bold', color='black')

ax.set_xticks(x_pos)
ax.set_xticklabels([f'Position +{p-1}' for p in positions], 
                   fontsize=13, fontweight='bold')
ax.set_ylabel('Enrichment ratio (Hits / Screen)', fontsize=14, fontweight='bold')
ax.set_title('Most enriched amino acid at each N-terminal position',
            fontsize=15, fontweight='bold', pad=15)

ax.axhline(y=1, color='gray', linestyle='--', linewidth=2, alpha=0.7,
          label='No enrichment')
ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
ax.set_axisbelow(True)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.set_ylim(-0.5, max(top_df['enrichment']) * 1.3)

# Legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=color, edgecolor='black', label=group)
                  for group, color in group_colors.items()]
ax.legend(handles=legend_elements, loc='upper right', ncol=2,
         fontsize=10, frameon=True)

plt.tight_layout()
plt.savefig('logo_results_enrichment/enrichment_paper_style_simple.png',
           dpi=300, bbox_inches='tight')
print("[OK] Saved: enrichment_paper_style_simple.png")
plt.close()

# Print summary
print("\n" + "="*80)
print("SUMMARY")
print("="*80)

for pos in positions:
    print(f"\nPosition +{pos-1} (Position {pos}):")
    pos_data = df[df['position'] == pos].nlargest(3, 'enrichment')
    for _, row in pos_data.iterrows():
        sig = '***' if row['p_value'] < 0.001 else '**' if row['p_value'] < 0.01 else '*' if row['p_value'] < 0.05 else 'ns'
        print(f"  {row['amino_acid']} ({row['group']}): {row['enrichment']:.2f}× (p={row['p_value']:.2e}) {sig}")

print("\n" + "="*80)
print("COMPLETE!")
print("="*80)
print("\nGenerated 4 paper-style figures:")
print("  1. enrichment_paper_style_main.png - Top 5 per position")
print("  2. enrichment_paper_style_grouped.png - By biochemical property")
print("  3. enrichment_paper_style_significant.png - All significant hits")
print("  4. enrichment_paper_style_simple.png - Clean single view")
print("="*80 + "\n")
