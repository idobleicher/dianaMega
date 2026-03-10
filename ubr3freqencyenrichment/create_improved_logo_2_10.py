#!/usr/bin/env python3
"""
Create improved, highly readable sequence logo for positions 2-10
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.transforms import Affine2D
import os

OUTPUT_DIR = "logo_results_54hits_first5"

print("Creating improved sequence logo for positions 2-10...")

# Load data
screen = pd.read_excel('UBR3 Nt screen.xlsx')
hits_file = pd.read_excel('ubr3_best (1).xlsx')

gene_col = hits_file.columns[2]
genes = hits_file[gene_col].unique()[:54]

hits_data = screen[screen['Gene_ID'].isin(genes)].groupby('Gene_ID').first().reset_index()
hits_sequences = hits_data['AA_seq'].tolist()
hits_sequences_10 = [seq[:10] for seq in hits_sequences if len(seq) >= 10]

print(f"Total sequences: {len(hits_sequences_10)}")

AA_COLORS = {
    'A': '#000000', 'V': '#000000', 'I': '#000000', 'L': '#000000', 
    'M': '#000000', 'F': '#000000', 'W': '#000000',
    'S': '#33AA33', 'T': '#33AA33', 'N': '#33AA33', 'Q': '#33AA33',
    'C': '#33AA33', 'Y': '#33AA33',
    'D': '#CC0000', 'E': '#CC0000',
    'K': '#0000CC', 'R': '#0000CC', 'H': '#6699CC',
    'G': '#FF9900', 'P': '#FF9900'
}

def calculate_information_content(sequences, position):
    from collections import Counter
    aas = [seq[position] for seq in sequences if position < len(seq)]
    total = len(aas)
    
    if total == 0:
        return {}, 0
    
    counts = Counter(aas)
    frequencies = {aa: count/total for aa, count in counts.items()}
    
    entropy = 0
    for freq in frequencies.values():
        if freq > 0:
            entropy -= freq * np.log2(freq)
    
    max_entropy = np.log2(20)
    information_content = max_entropy - entropy
    
    letter_heights = {aa: freq * information_content for aa, freq in frequencies.items()}
    return letter_heights, information_content

def draw_letter(ax, letter, x, y, height, width, color):
    if height <= 0.01:
        return
    
    fp = FontProperties(family='Arial', weight='bold', size=100)
    text_path = TextPath((0, 0), letter, size=1, prop=fp)
    
    bbox = text_path.get_extents()
    letter_width = bbox.width
    letter_height = bbox.height
    
    scale_x = width / letter_width
    scale_y = height / letter_height
    
    t = Affine2D().scale(scale_x, scale_y).translate(x, y)
    text_path = text_path.transformed(t)
    
    # Thicker border for better visibility
    patch = PathPatch(text_path, facecolor=color, edgecolor='black', linewidth=2)
    ax.add_patch(patch)

# 1. TOP 3 VERSION (CLEAREST)
print("\n1. Creating top 3 version (clearest)...")
fig, ax = plt.subplots(figsize=(20, 8))

for pos_idx in range(1, 10):
    letter_heights_all, ic = calculate_information_content(hits_sequences_10, pos_idx)
    
    # Keep only top 3
    top3 = dict(sorted(letter_heights_all.items(), key=lambda x: x[1], reverse=True)[:3])
    sorted_letters = sorted(top3.items(), key=lambda x: x[1])
    
    y_position = 0
    for letter, height in sorted_letters:
        color = AA_COLORS.get(letter, '#666666')
        draw_letter(ax, letter, pos_idx - 0.45, y_position, height, 0.9, color)
        y_position += height

ax.set_xlim(0.3, 9.7)
ax.set_ylim(0, 1.5)
ax.set_xticks(range(1, 10))
ax.set_xticklabels(['2', '3', '4', '5', '6', '7', '8', '9', '10'], fontsize=14, weight='bold')
ax.set_xlabel('Position', fontsize=16, weight='bold')
ax.set_ylabel('Information Content (bits)', fontsize=16, weight='bold')
ax.set_title('UBR3 N-Terminal Recognition Motif (Positions 2-10)\nTop 3 Amino Acids per Position', 
             fontsize=18, weight='bold', pad=20)

ax.set_yticks([0, 0.5, 1.0, 1.5])
ax.grid(axis='y', alpha=0.25, linestyle='-', linewidth=1, color='gray')
ax.set_facecolor('#F8F8F8')

from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#CC0000', edgecolor='black', label='Acidic (D,E)', linewidth=1.5),
    Patch(facecolor='#0000CC', edgecolor='black', label='Basic (K,R)', linewidth=1.5),
    Patch(facecolor='#6699CC', edgecolor='black', label='His (H)', linewidth=1.5),
    Patch(facecolor='#000000', edgecolor='black', label='Hydrophobic', linewidth=1.5),
    Patch(facecolor='#33AA33', edgecolor='black', label='Polar', linewidth=1.5),
    Patch(facecolor='#FF9900', edgecolor='black', label='G, P', linewidth=1.5)
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=13, framealpha=0.95,
         title='AA Properties', title_fontsize=13, edgecolor='black')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'sequence_logo_pos2_10_improved_top3.png'), dpi=300, bbox_inches='tight')
plt.close()
print("Saved: sequence_logo_pos2_10_improved_top3.png")

# 2. SPLIT PANEL VERSION (2-6 and 7-10)
print("\n2. Creating split panel version...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

# Panel 1: Positions 2-6
for pos_idx in range(1, 6):
    letter_heights_all, ic = calculate_information_content(hits_sequences_10, pos_idx)
    top4 = dict(sorted(letter_heights_all.items(), key=lambda x: x[1], reverse=True)[:4])
    sorted_letters = sorted(top4.items(), key=lambda x: x[1])
    
    y_position = 0
    for letter, height in sorted_letters:
        color = AA_COLORS.get(letter, '#666666')
        draw_letter(ax1, letter, pos_idx - 0.45, y_position, height, 0.9, color)
        y_position += height

ax1.set_xlim(0.3, 5.7)
ax1.set_ylim(0, 1.5)
ax1.set_xticks(range(1, 6))
ax1.set_xticklabels(['2', '3', '4', '5', '6'], fontsize=14, weight='bold')
ax1.set_xlabel('Position', fontsize=16, weight='bold')
ax1.set_ylabel('Information Content (bits)', fontsize=16, weight='bold')
ax1.set_title('A. N-Terminal Core (Positions 2-6)', fontsize=16, weight='bold', loc='left', pad=15)
ax1.set_yticks([0, 0.5, 1.0, 1.5])
ax1.grid(axis='y', alpha=0.25, linestyle='-', linewidth=1, color='gray')
ax1.set_facecolor('#F8F8F8')

# Panel 2: Positions 7-10
for pos_idx in range(6, 10):
    letter_heights_all, ic = calculate_information_content(hits_sequences_10, pos_idx)
    top4 = dict(sorted(letter_heights_all.items(), key=lambda x: x[1], reverse=True)[:4])
    sorted_letters = sorted(top4.items(), key=lambda x: x[1])
    
    # Adjust x position to start from 1
    y_position = 0
    for letter, height in sorted_letters:
        color = AA_COLORS.get(letter, '#666666')
        draw_letter(ax2, letter, (pos_idx - 5) - 0.45, y_position, height, 0.9, color)
        y_position += height

ax2.set_xlim(0.3, 4.7)
ax2.set_ylim(0, 1.5)
ax2.set_xticks(range(1, 5))
ax2.set_xticklabels(['7', '8', '9', '10'], fontsize=14, weight='bold')
ax2.set_xlabel('Position', fontsize=16, weight='bold')
ax2.set_ylabel('Information Content (bits)', fontsize=16, weight='bold')
ax2.set_title('B. Extended Region (Positions 7-10)', fontsize=16, weight='bold', loc='left', pad=15)
ax2.set_yticks([0, 0.5, 1.0, 1.5])
ax2.grid(axis='y', alpha=0.25, linestyle='-', linewidth=1, color='gray')
ax2.set_facecolor('#F8F8F8')

# Shared legend
ax2.legend(handles=legend_elements, loc='upper right', fontsize=12, framealpha=0.95,
          title='AA Properties', title_fontsize=12, edgecolor='black')

fig.suptitle('UBR3 N-Terminal Recognition Motif (Positions 2-10)\nTop 4 Amino Acids per Position', 
             fontsize=20, weight='bold', y=0.98)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(os.path.join(OUTPUT_DIR, 'sequence_logo_pos2_10_split_panels.png'), dpi=300, bbox_inches='tight')
plt.close()
print("Saved: sequence_logo_pos2_10_split_panels.png")

# 3. EXTRA WIDE VERSION (top 4)
print("\n3. Creating extra wide version...")
fig, ax = plt.subplots(figsize=(24, 7))

for pos_idx in range(1, 10):
    letter_heights_all, ic = calculate_information_content(hits_sequences_10, pos_idx)
    top4 = dict(sorted(letter_heights_all.items(), key=lambda x: x[1], reverse=True)[:4])
    sorted_letters = sorted(top4.items(), key=lambda x: x[1])
    
    y_position = 0
    for letter, height in sorted_letters:
        color = AA_COLORS.get(letter, '#666666')
        draw_letter(ax, letter, pos_idx - 0.45, y_position, height, 0.9, color)
        y_position += height

ax.set_xlim(0.3, 9.7)
ax.set_ylim(0, 1.5)
ax.set_xticks(range(1, 10))
ax.set_xticklabels(['Pos 2', 'Pos 3', 'Pos 4', 'Pos 5', 'Pos 6', 'Pos 7', 'Pos 8', 'Pos 9', 'Pos 10'], 
                    fontsize=14, weight='bold')
ax.set_xlabel('Position', fontsize=18, weight='bold')
ax.set_ylabel('Information Content (bits)', fontsize=18, weight='bold')
ax.set_title('UBR3 N-Terminal Recognition Motif - Extended Sequence Logo (Positions 2-10)\nTop 4 Amino Acids per Position', 
             fontsize=20, weight='bold', pad=20)

ax.set_yticks([0, 0.5, 1.0, 1.5])
ax.grid(axis='y', alpha=0.25, linestyle='-', linewidth=1, color='gray')
ax.set_facecolor('#F8F8F8')

ax.legend(handles=legend_elements, loc='upper right', fontsize=14, framealpha=0.95,
         title='AA Properties', title_fontsize=14, edgecolor='black')

# Add key regions annotation
key_regions = (
    "Key Regions:\n"
    "• Pos 2-3: Core motif (G/P, D/E)\n"
    "• Pos 7-8: Arg-rich region\n"
    "• Pos 4-6, 9-10: Variable"
)
ax.text(0.02, 0.98, key_regions, transform=ax.transAxes, ha='left', va='top',
        fontsize=11, family='monospace',
        bbox=dict(boxstyle='round,pad=0.7', facecolor='lightyellow', 
                 edgecolor='black', linewidth=1.5, alpha=0.95))

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'sequence_logo_pos2_10_extra_wide.png'), dpi=300, bbox_inches='tight')
plt.close()
print("Saved: sequence_logo_pos2_10_extra_wide.png")

print("\n" + "="*80)
print("IMPROVED SEQUENCE LOGOS CREATED!")
print("="*80)
print(f"\nCreated 3 improved versions in {OUTPUT_DIR}/:")
print("  1. sequence_logo_pos2_10_improved_top3.png - Top 3 per position (CLEAREST)")
print("  2. sequence_logo_pos2_10_split_panels.png - Split into 2 panels (BEST ORGANIZATION)")
print("  3. sequence_logo_pos2_10_extra_wide.png - Extra wide format (FOR PRESENTATIONS)")
print("\nAll versions have:")
print("  - Thicker letter borders (2px)")
print("  - Better spacing")
print("  - Only top 3-4 AAs per position")
print("  - Light gray background")
print("  - Clear, readable letters")
print("\nDone!")
