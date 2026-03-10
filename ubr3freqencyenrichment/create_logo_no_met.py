#!/usr/bin/env python3
"""
Create traditional sequence logo starting from position 2 (no methionine)
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
print("CREATING SEQUENCE LOGO (POSITIONS 2-5, NO METHIONINE)")
print("="*80)

# Load data
print("\nLoading sequences...")
hits_data = pd.read_csv(os.path.join(OUTPUT_DIR, 'hits_first5_sequences.csv'))
sequences = hits_data['first_5_residues'].tolist()

print(f"Total sequences: {len(sequences)}")

# Amino acid color scheme (standard WebLogo coloring)
AA_COLORS = {
    # Hydrophobic (black)
    'A': '#000000', 'V': '#000000', 'I': '#000000', 'L': '#000000', 
    'M': '#000000', 'F': '#000000', 'W': '#000000',
    # Polar (green)
    'S': '#33AA33', 'T': '#33AA33', 'N': '#33AA33', 'Q': '#33AA33',
    'C': '#33AA33', 'Y': '#33AA33',
    # Acidic (red)
    'D': '#CC0000', 'E': '#CC0000',
    # Basic (blue)
    'K': '#0000CC', 'R': '#0000CC', 'H': '#6699CC',
    # Special (orange/brown)
    'G': '#FF9900', 'P': '#FF9900'
}

def calculate_information_content(sequences, position):
    """Calculate information content (bits) at a position"""
    from collections import Counter
    
    aas = [seq[position] for seq in sequences if position < len(seq)]
    total = len(aas)
    
    if total == 0:
        return {}, 0
    
    counts = Counter(aas)
    frequencies = {aa: count/total for aa, count in counts.items()}
    
    # Calculate Shannon entropy
    entropy = 0
    for freq in frequencies.values():
        if freq > 0:
            entropy -= freq * np.log2(freq)
    
    # Information content
    max_entropy = np.log2(20)
    information_content = max_entropy - entropy
    
    # Height of each letter = frequency * information_content
    letter_heights = {aa: freq * information_content 
                     for aa, freq in frequencies.items()}
    
    return letter_heights, information_content

def draw_letter(ax, letter, x, y, height, width, color):
    """Draw a letter using matplotlib TextPath"""
    if height <= 0:
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
    
    patch = PathPatch(text_path, facecolor=color, edgecolor='none')
    ax.add_patch(patch)

print("\nCreating sequence logos (positions 2-5)...")

# 1. STANDARD LOGO
fig, ax = plt.subplots(figsize=(10, 6))

# Start from position 1 (index 1 = position 2)
for pos_idx in range(1, 5):  # Positions 1-4 (which are positions 2-5 in sequence)
    letter_heights, ic = calculate_information_content(sequences, pos_idx)
    sorted_letters = sorted(letter_heights.items(), key=lambda x: x[1])
    
    # Plot at x = pos_idx (will be 1, 2, 3, 4)
    y_position = 0
    for letter, height in sorted_letters:
        color = AA_COLORS.get(letter, '#666666')
        draw_letter(ax, letter, pos_idx - 0.4, y_position, height, 0.8, color)
        y_position += height

ax.set_xlim(0.5, 4.5)
ax.set_ylim(0, 4.5)
ax.set_xticks([1, 2, 3, 4])
ax.set_xticklabels(['2', '3', '4', '5'], fontsize=13, weight='bold')
ax.set_xlabel('Position', fontsize=14, weight='bold')
ax.set_ylabel('Information Content (bits)', fontsize=14, weight='bold')
ax.set_title('Sequence Logo: N-Terminal Motif (Positions 2-5)\n54 Hits - Information Content Based', 
             fontsize=16, weight='bold', pad=15)

ax.set_yticks([0, 1, 2, 3, 4])
ax.grid(axis='y', alpha=0.3, linestyle='--')

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
ax.legend(handles=legend_elements, loc='upper right', fontsize=10, framealpha=0.9,
         title='AA Properties', title_fontsize=10)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'sequence_logo_pos2_5.png'), dpi=300, bbox_inches='tight')
plt.close()

print("Saved: sequence_logo_pos2_5.png")

# 2. PUBLICATION VERSION
fig, ax = plt.subplots(figsize=(12, 7))

for pos_idx in range(1, 5):
    letter_heights, ic = calculate_information_content(sequences, pos_idx)
    sorted_letters = sorted(letter_heights.items(), key=lambda x: x[1])
    
    y_position = 0
    for letter, height in sorted_letters:
        color = AA_COLORS.get(letter, '#666666')
        draw_letter(ax, letter, pos_idx - 0.45, y_position, height, 0.9, color)
        y_position += height

ax.set_xlim(0.3, 4.7)
ax.set_ylim(0, 4.5)
ax.set_xticks([1, 2, 3, 4])
ax.set_xticklabels(['Position 2', 'Position 3', 'Position 4', 'Position 5'], 
                    fontsize=12, weight='bold')
ax.set_xlabel('Sequence Position', fontsize=14, weight='bold')
ax.set_ylabel('Information Content (bits)', fontsize=14, weight='bold')
ax.set_title('UBR3 N-Terminal Recognition Motif (Positions 2-5)\nSequence Logo - 53 Selected Hits', 
             fontsize=16, weight='bold', pad=20)

ax.set_yticks([0, 1, 2, 3, 4])
ax.grid(axis='y', alpha=0.2, linestyle='-', linewidth=0.5)

ax.legend(handles=legend_elements, loc='upper right', fontsize=11, framealpha=0.95,
         title='Amino Acid Properties', title_fontsize=11)

note = "Letter height = frequency × information content\nPositions shown: 2-5 (Met at position 1 excluded)"
ax.text(0.02, 0.98, note, transform=ax.transAxes, ha='left', va='top',
        fontsize=9, bbox=dict(boxstyle='round,pad=0.5', facecolor='wheat', alpha=0.8))

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'sequence_logo_pos2_5_publication.png'), dpi=300, bbox_inches='tight')
plt.close()

print("Saved: sequence_logo_pos2_5_publication.png")

# 3. WIDE FORMAT
fig, ax = plt.subplots(figsize=(16, 6))

for pos_idx in range(1, 5):
    letter_heights, ic = calculate_information_content(sequences, pos_idx)
    sorted_letters = sorted(letter_heights.items(), key=lambda x: x[1])
    
    y_position = 0
    for letter, height in sorted_letters:
        color = AA_COLORS.get(letter, '#666666')
        draw_letter(ax, letter, pos_idx - 0.45, y_position, height, 0.9, color)
        y_position += height

ax.set_xlim(0.3, 4.7)
ax.set_ylim(0, 4.5)
ax.set_xticks([1, 2, 3, 4])
ax.set_xticklabels(['Pos 2', 'Pos 3', 'Pos 4', 'Pos 5'], 
                    fontsize=16, weight='bold')
ax.set_xlabel('Position', fontsize=18, weight='bold')
ax.set_ylabel('Information (bits)', fontsize=18, weight='bold')
ax.set_title('UBR3 Recognition Motif - Sequence Logo (Positions 2-5)', 
             fontsize=20, weight='bold', pad=15)

ax.set_yticks([0, 1, 2, 3, 4])
ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=1)

ax.legend(handles=legend_elements, loc='upper right', fontsize=13, framealpha=0.95,
         ncol=2, title='AA Properties', title_fontsize=13)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'sequence_logo_pos2_5_wide.png'), dpi=300, bbox_inches='tight')
plt.close()

print("Saved: sequence_logo_pos2_5_wide.png")

# Calculate statistics for positions 2-5
print("\n" + "="*80)
print("STATISTICS (POSITIONS 2-5)")
print("="*80)

for pos_idx in range(1, 5):
    letter_heights, ic = calculate_information_content(sequences, pos_idx)
    
    print(f"\nPosition {pos_idx + 1}:")
    print(f"  Information Content: {ic:.2f} bits")
    print(f"  Conservation: {ic/4.322*100:.1f}%")
    
    sorted_letters = sorted(letter_heights.items(), key=lambda x: x[1], reverse=True)[:5]
    print(f"  Top 5 amino acids:")
    for letter, height in sorted_letters:
        freq = height / ic if ic > 0 else 0
        print(f"    {letter}: {freq*100:.1f}%")

print("\n" + "="*80)
print("SEQUENCE LOGOS CREATED (POSITIONS 2-5)!")
print("="*80)
print(f"\nCreated files in {OUTPUT_DIR}/:")
print("  1. sequence_logo_pos2_5.png - Standard logo (positions 2-5)")
print("  2. sequence_logo_pos2_5_publication.png - Publication version")
print("  3. sequence_logo_pos2_5_wide.png - Wide format for presentations")
print("\nAll logos start from position 2 (methionine excluded)")
print("\nDone!")
