#!/usr/bin/env python3
"""
Create sequence logos for positions 2-10 (extended N-terminal motif)
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

print("="*80)
print("CREATING SEQUENCE LOGOS FOR POSITIONS 2-10")
print("="*80)

# Load data - need to get first 10 positions
print("\nLoading sequences...")
screen = pd.read_excel('UBR3 Nt screen.xlsx')
hits_file = pd.read_excel('ubr3_best (1).xlsx')

# Get gene list (53 genes)
gene_col = hits_file.columns[2]
genes = hits_file[gene_col].unique()[:54]
print(f"Using {len(genes)} genes")

# Get ONE sequence per gene
hits_data = screen[screen['Gene_ID'].isin(genes)].groupby('Gene_ID').first().reset_index()
hits_sequences = hits_data['AA_seq'].tolist()

# Extract first 10 positions
hits_sequences_10 = [seq[:10] for seq in hits_sequences if len(seq) >= 10]

print(f"Total sequences: {len(hits_sequences_10)}")
print(f"Sample: {hits_sequences_10[:3]}")

# Amino acid colors
AA_COLORS = {
    'A': '#000000', 'V': '#000000', 'I': '#000000', 'L': '#000000', 
    'M': '#000000', 'F': '#000000', 'W': '#000000',
    'S': '#33AA33', 'T': '#33AA33', 'N': '#33AA33', 'Q': '#33AA33',
    'C': '#33AA33', 'Y': '#33AA33',
    'D': '#CC0000', 'E': '#CC0000',
    'K': '#0000CC', 'R': '#0000CC', 'H': '#6699CC',
    'G': '#FF9900', 'P': '#FF9900'
}

def calculate_information_content(sequences, position, min_freq=0.0):
    """Calculate information content at a position"""
    from collections import Counter
    
    aas = [seq[position] for seq in sequences if position < len(seq)]
    total = len(aas)
    
    if total == 0:
        return {}, 0
    
    counts = Counter(aas)
    all_frequencies = {aa: count/total for aa, count in counts.items()}
    
    # Filter by frequency
    if min_freq > 0:
        all_frequencies = {aa: freq for aa, freq in all_frequencies.items() if freq >= min_freq}
    
    if not all_frequencies:
        return {}, 0
    
    # Calculate entropy
    entropy = 0
    for freq in all_frequencies.values():
        if freq > 0:
            entropy -= freq * np.log2(freq)
    
    max_entropy = np.log2(20)
    information_content = max_entropy - entropy
    
    # Letter heights
    letter_heights = {aa: freq * information_content 
                     for aa, freq in all_frequencies.items()}
    
    return letter_heights, information_content

def draw_letter(ax, letter, x, y, height, width, color):
    """Draw a letter"""
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
    
    patch = PathPatch(text_path, facecolor=color, edgecolor='black', linewidth=1)
    ax.add_patch(patch)

print("\nCreating sequence logos...")

# 1. STANDARD LOGO (all AAs)
fig, ax = plt.subplots(figsize=(18, 7))

for pos_idx in range(1, 10):  # Positions 1-9 (which are 2-10 in sequence)
    letter_heights, ic = calculate_information_content(hits_sequences_10, pos_idx, min_freq=0.0)
    sorted_letters = sorted(letter_heights.items(), key=lambda x: x[1])
    
    y_position = 0
    for letter, height in sorted_letters:
        color = AA_COLORS.get(letter, '#666666')
        draw_letter(ax, letter, pos_idx - 0.4, y_position, height, 0.8, color)
        y_position += height

ax.set_xlim(0.5, 9.5)
ax.set_ylim(0, 1.5)
ax.set_xticks(range(1, 10))
ax.set_xticklabels(['2', '3', '4', '5', '6', '7', '8', '9', '10'], fontsize=12, weight='bold')
ax.set_xlabel('Position', fontsize=14, weight='bold')
ax.set_ylabel('Information Content (bits)', fontsize=14, weight='bold')
ax.set_title('UBR3 N-Terminal Recognition Motif (Positions 2-10)\nSequence Logo - 53 Hits', 
             fontsize=16, weight='bold', pad=15)

ax.set_yticks([0, 0.5, 1.0, 1.5])
ax.grid(axis='y', alpha=0.25, linestyle='-', linewidth=0.8, color='gray')
ax.set_facecolor('#F8F8F8')

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#CC0000', edgecolor='black', label='Acidic (D,E)'),
    Patch(facecolor='#0000CC', edgecolor='black', label='Basic (K,R)'),
    Patch(facecolor='#6699CC', edgecolor='black', label='His (H)'),
    Patch(facecolor='#000000', edgecolor='black', label='Hydrophobic'),
    Patch(facecolor='#33AA33', edgecolor='black', label='Polar'),
    Patch(facecolor='#FF9900', edgecolor='black', label='G, P')
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=11, framealpha=0.95,
         title='AA Properties', title_fontsize=11, edgecolor='black')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'sequence_logo_pos2_10.png'), dpi=300, bbox_inches='tight')
plt.close()

print("Saved: sequence_logo_pos2_10.png")

# 2. CLEAN VERSION (>5% frequency)
fig, ax = plt.subplots(figsize=(18, 7))

for pos_idx in range(1, 10):
    letter_heights, ic = calculate_information_content(hits_sequences_10, pos_idx, min_freq=0.05)
    sorted_letters = sorted(letter_heights.items(), key=lambda x: x[1])
    
    y_position = 0
    for letter, height in sorted_letters:
        color = AA_COLORS.get(letter, '#666666')
        draw_letter(ax, letter, pos_idx - 0.4, y_position, height, 0.8, color)
        y_position += height

ax.set_xlim(0.5, 9.5)
ax.set_ylim(0, 1.5)
ax.set_xticks(range(1, 10))
ax.set_xticklabels(['2', '3', '4', '5', '6', '7', '8', '9', '10'], fontsize=12, weight='bold')
ax.set_xlabel('Position', fontsize=14, weight='bold')
ax.set_ylabel('Information Content (bits)', fontsize=14, weight='bold')
ax.set_title('UBR3 N-Terminal Recognition Motif (Positions 2-10)\nSequence Logo - Amino Acids >5% Frequency', 
             fontsize=16, weight='bold', pad=15)

ax.set_yticks([0, 0.5, 1.0, 1.5])
ax.grid(axis='y', alpha=0.25, linestyle='-', linewidth=0.8, color='gray')
ax.set_facecolor('#F8F8F8')

ax.legend(handles=legend_elements, loc='upper right', fontsize=11, framealpha=0.95,
         title='AA Properties', title_fontsize=11, edgecolor='black')

note = "Showing amino acids with >5% frequency\nLetter height = frequency × information content"
ax.text(0.02, 0.98, note, transform=ax.transAxes, ha='left', va='top',
        fontsize=9, bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', 
                               edgecolor='black', linewidth=1.5, alpha=0.9))

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'sequence_logo_pos2_10_clean.png'), dpi=300, bbox_inches='tight')
plt.close()

print("Saved: sequence_logo_pos2_10_clean.png")

# 3. PUBLICATION VERSION (top 5 per position)
fig, ax = plt.subplots(figsize=(18, 7))

for pos_idx in range(1, 10):
    letter_heights_all, ic = calculate_information_content(hits_sequences_10, pos_idx, min_freq=0.0)
    
    # Keep only top 5
    top5 = dict(sorted(letter_heights_all.items(), key=lambda x: x[1], reverse=True)[:5])
    sorted_letters = sorted(top5.items(), key=lambda x: x[1])
    
    y_position = 0
    for letter, height in sorted_letters:
        color = AA_COLORS.get(letter, '#666666')
        draw_letter(ax, letter, pos_idx - 0.4, y_position, height, 0.8, color)
        y_position += height

ax.set_xlim(0.5, 9.5)
ax.set_ylim(0, 1.5)
ax.set_xticks(range(1, 10))
ax.set_xticklabels(['2', '3', '4', '5', '6', '7', '8', '9', '10'], fontsize=12, weight='bold')
ax.set_xlabel('Position', fontsize=14, weight='bold')
ax.set_ylabel('Information Content (bits)', fontsize=14, weight='bold')
ax.set_title('UBR3 N-Terminal Recognition Motif (Positions 2-10)\nSequence Logo - Top 5 Amino Acids per Position', 
             fontsize=16, weight='bold', pad=15)

ax.set_yticks([0, 0.5, 1.0, 1.5])
ax.grid(axis='y', alpha=0.25, linestyle='-', linewidth=0.8, color='gray')
ax.set_facecolor('#F8F8F8')

ax.legend(handles=legend_elements, loc='upper right', fontsize=11, framealpha=0.95,
         title='AA Properties', title_fontsize=11, edgecolor='black')

note = "Showing top 5 amino acids per position\nLetter height = frequency × information content"
ax.text(0.02, 0.98, note, transform=ax.transAxes, ha='left', va='top',
        fontsize=9, bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', 
                               edgecolor='black', linewidth=1.5, alpha=0.9))

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'sequence_logo_pos2_10_top5.png'), dpi=300, bbox_inches='tight')
plt.close()

print("Saved: sequence_logo_pos2_10_top5.png")

# Calculate and display statistics
print("\n" + "="*80)
print("STATISTICS (POSITIONS 2-10)")
print("="*80)

for pos_idx in range(1, 10):
    letter_heights, ic = calculate_information_content(hits_sequences_10, pos_idx, min_freq=0.0)
    
    print(f"\nPosition {pos_idx + 1}:")
    print(f"  Information Content: {ic:.2f} bits")
    print(f"  Conservation: {ic/4.322*100:.1f}%")
    
    sorted_letters = sorted(letter_heights.items(), key=lambda x: x[1], reverse=True)[:5]
    print(f"  Top 5 amino acids:")
    for letter, height in sorted_letters:
        freq = height / ic if ic > 0 else 0
        print(f"    {letter}: {freq*100:.1f}%")

print("\n" + "="*80)
print("SEQUENCE LOGOS CREATED (POSITIONS 2-10)!")
print("="*80)
print(f"\nCreated files in {OUTPUT_DIR}/:")
print("  1. sequence_logo_pos2_10.png - All amino acids")
print("  2. sequence_logo_pos2_10_clean.png - AAs >5% frequency (cleaner)")
print("  3. sequence_logo_pos2_10_top5.png - Top 5 per position (clearest)")
print("\nAll logos show extended N-terminal motif (positions 2-10)")
print("\nDone!")
