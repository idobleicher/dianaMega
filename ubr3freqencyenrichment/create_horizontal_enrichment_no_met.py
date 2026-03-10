#!/usr/bin/env python3
"""
Create horizontal enrichment plot for positions 2-5 (no methionine)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

OUTPUT_DIR = "logo_results_54hits_first5"

print("Creating horizontal enrichment plot (positions 2-5, no Met)...")

# Load enrichment data
enrichment_df = pd.read_csv(os.path.join(OUTPUT_DIR, 'enrichment_first5_data.csv'))

# Filter to positions 2-5 only
enrichment_df = enrichment_df[enrichment_df['position'].isin([2, 3, 4, 5])]

# Create horizontal plot
fig, axes = plt.subplots(1, 4, figsize=(20, 8), sharey=True)

for idx, pos in enumerate([2, 3, 4, 5]):
    ax = axes[idx]
    pos_data = enrichment_df[enrichment_df['position'] == pos].sort_values('enrichment', ascending=True)
    
    # Color by enrichment level
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
            edgecolor='black', linewidth=0.8)
    ax.axvline(x=1, color='black', linestyle='--', linewidth=2, alpha=0.7)
    ax.axvline(x=2, color='gray', linestyle=':', linewidth=1.5, alpha=0.5)
    ax.set_xlabel('Enrichment', fontsize=13, weight='bold')
    ax.set_title(f'Position {pos}', fontsize=15, weight='bold', pad=10)
    ax.grid(axis='x', alpha=0.3)
    ax.set_xlim(0, max(5.5, pos_data['enrichment'].max() * 1.1))
    
    # Make y-labels bold on first panel only
    if idx == 0:
        ax.set_ylabel('Amino Acid', fontsize=13, weight='bold')
    
    # Highlight top 3 enriched
    top3_aas = pos_data.nlargest(3, 'enrichment')['amino_acid'].tolist()
    for label in ax.get_yticklabels():
        if label.get_text() in top3_aas:
            label.set_weight('bold')
            label.set_color('darkred')
            label.set_fontsize(11)

# Overall title
fig.suptitle('Amino Acid Enrichment: Positions 2-5 (54 Hits vs 16,514 Screen)', 
             fontsize=18, weight='bold', y=0.98)

# Add legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='darkred', edgecolor='black', label='>3x enrichment', linewidth=1.5),
    Patch(facecolor='red', edgecolor='black', label='2-3x enrichment', linewidth=1.5),
    Patch(facecolor='orange', edgecolor='black', label='1.5-2x enrichment', linewidth=1.5),
    Patch(facecolor='gold', edgecolor='black', label='1-1.5x enrichment', linewidth=1.5),
    Patch(facecolor='lightblue', edgecolor='black', label='0.5-1x (depleted)', linewidth=1.5),
    Patch(facecolor='blue', edgecolor='black', label='<0.5x (strongly depleted)', linewidth=1.5)
]
fig.legend(handles=legend_elements, loc='lower center', ncol=3, fontsize=11, 
          bbox_to_anchor=(0.5, -0.05), framealpha=0.95)

plt.tight_layout(rect=[0, 0.02, 1, 0.96])
plt.savefig(os.path.join(OUTPUT_DIR, 'enrichment_horizontal_no_met.png'), 
            dpi=300, bbox_inches='tight')
plt.close()

print("Saved: enrichment_horizontal_no_met.png")

# Also create a cleaner version with only enriched (>1.0)
fig, axes = plt.subplots(1, 4, figsize=(20, 8), sharey=True)

for idx, pos in enumerate([2, 3, 4, 5]):
    ax = axes[idx]
    pos_data = enrichment_df[enrichment_df['position'] == pos].copy()
    
    # Only show enriched (>1.0)
    pos_data = pos_data[pos_data['enrichment'] > 1.0].sort_values('enrichment', ascending=True)
    
    # Color by enrichment level
    colors = []
    for e in pos_data['enrichment']:
        if e > 3:
            colors.append('darkred')
        elif e > 2:
            colors.append('red')
        elif e > 1.5:
            colors.append('orange')
        else:
            colors.append('gold')
    
    bars = ax.barh(pos_data['amino_acid'], pos_data['enrichment'], color=colors, 
                   edgecolor='black', linewidth=0.8)
    ax.axvline(x=1, color='black', linestyle='--', linewidth=2, alpha=0.7, label='No enrichment')
    ax.axvline(x=2, color='gray', linestyle=':', linewidth=1.5, alpha=0.5, label='2x')
    ax.set_xlabel('Enrichment Ratio', fontsize=13, weight='bold')
    ax.set_title(f'Position {pos}', fontsize=15, weight='bold', pad=10)
    ax.grid(axis='x', alpha=0.3, linewidth=0.8)
    ax.set_xlim(0.8, max(5, pos_data['enrichment'].max() * 1.1))
    
    if idx == 0:
        ax.set_ylabel('Amino Acid', fontsize=13, weight='bold')
        ax.legend(loc='lower right', fontsize=9)
    
    # Bold all labels
    for label in ax.get_yticklabels():
        label.set_weight('bold')
        label.set_fontsize(11)

fig.suptitle('Enriched Amino Acids (>1.0x): Positions 2-5\nN-Terminal Recognition Motif (54 Hits)', 
             fontsize=18, weight='bold', y=0.98)

legend_elements_clean = [
    Patch(facecolor='darkred', edgecolor='black', label='>3x', linewidth=1.5),
    Patch(facecolor='red', edgecolor='black', label='2-3x', linewidth=1.5),
    Patch(facecolor='orange', edgecolor='black', label='1.5-2x', linewidth=1.5),
    Patch(facecolor='gold', edgecolor='black', label='1-1.5x', linewidth=1.5)
]
fig.legend(handles=legend_elements_clean, loc='lower center', ncol=4, fontsize=12, 
          bbox_to_anchor=(0.5, -0.03), framealpha=0.95, title='Enrichment Level',
          title_fontsize=12)

plt.tight_layout(rect=[0, 0.03, 1, 0.96])
plt.savefig(os.path.join(OUTPUT_DIR, 'enrichment_horizontal_enriched_only.png'), 
            dpi=300, bbox_inches='tight')
plt.close()

print("Saved: enrichment_horizontal_enriched_only.png")

# Create publication version with annotations
fig, axes = plt.subplots(1, 4, figsize=(18, 7), sharey=True)

for idx, pos in enumerate([2, 3, 4, 5]):
    ax = axes[idx]
    pos_data = enrichment_df[enrichment_df['position'] == pos].copy()
    pos_data = pos_data[pos_data['enrichment'] > 1.15].sort_values('enrichment', ascending=True)
    
    colors = []
    for e in pos_data['enrichment']:
        if e > 3:
            colors.append('#8B0000')  # Dark red
        elif e > 2:
            colors.append('#DC143C')  # Crimson
        elif e > 1.5:
            colors.append('#FF8C00')  # Dark orange
        else:
            colors.append('#FFD700')  # Gold
    
    bars = ax.barh(pos_data['amino_acid'], pos_data['enrichment'], color=colors, 
                   edgecolor='black', linewidth=1)
    
    # Add enrichment values on bars
    for i, (aa, enrich) in enumerate(zip(pos_data['amino_acid'], pos_data['enrichment'])):
        if enrich > 1.3:
            ax.text(enrich + 0.15, i, f'{enrich:.2f}x', 
                   va='center', fontsize=9, weight='bold')
    
    ax.axvline(x=1, color='black', linestyle='-', linewidth=2, alpha=0.7)
    ax.axvline(x=2, color='gray', linestyle='--', linewidth=1.5, alpha=0.5)
    ax.axvline(x=3, color='gray', linestyle=':', linewidth=1.5, alpha=0.5)
    
    ax.set_xlabel('Enrichment', fontsize=12, weight='bold')
    ax.set_title(f'Position {pos}', fontsize=14, weight='bold', pad=10)
    ax.grid(axis='x', alpha=0.25, linewidth=0.5)
    ax.set_xlim(0.9, max(5, pos_data['enrichment'].max() * 1.15))
    
    if idx == 0:
        ax.set_ylabel('Amino Acid', fontsize=12, weight='bold')
    
    # Bold labels
    for label in ax.get_yticklabels():
        label.set_weight('bold')
        label.set_fontsize(10)

fig.suptitle('UBR3 N-Terminal Enrichment Profile (Positions 2-5)', 
             fontsize=17, weight='bold', y=0.98)

# Key findings annotation
key_text = (
    "Key Enrichments:\n"
    "Pos 2: P(4.39x), G(3.25x)\n"
    "Pos 3: D(4.57x) ★HIGHEST\n"
    "Pos 4: Y(4.37x) ★KEY\n"
    "Pos 5: I(2.96x), K(2.49x)"
)
fig.text(0.98, 0.5, key_text, transform=fig.transFigure,
         fontsize=10, family='monospace', weight='bold',
         bbox=dict(boxstyle='round,pad=0.8', facecolor='lightyellow', 
                  edgecolor='black', linewidth=2, alpha=0.95),
         ha='right', va='center')

plt.tight_layout(rect=[0, 0, 0.97, 0.96])
plt.savefig(os.path.join(OUTPUT_DIR, 'enrichment_horizontal_publication.png'), 
            dpi=300, bbox_inches='tight')
plt.close()

print("Saved: enrichment_horizontal_publication.png")

print("\nCreated 3 horizontal layouts (positions 2-5, no Met):")
print("  1. enrichment_horizontal_no_met.png - All amino acids")
print("  2. enrichment_horizontal_enriched_only.png - Only enriched (>1.0x)")
print("  3. enrichment_horizontal_publication.png - Publication version with annotations")
print("\nDone!")
