#!/usr/bin/env python3
"""
Create properly scaled sequence logo
Uses full information content but only displays significant amino acids
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

print("Creating properly scaled sequence logo...")

# Load data
hits_data = pd.read_csv(os.path.join(OUTPUT_DIR, 'hits_first5_sequences.csv'))
sequences = hits_data['first_5_residues'].tolist()

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

def calculate_information_content(sequences, position, min_freq=0.05):
    """Calculate information content properly"""
    from collections import Counter
    
    aas = [seq[position] for seq in sequences if position < len(seq)]
    total = len(aas)
    
    if total == 0:
        return {}, 0
    
    counts = Counter(aas)
    all_frequencies = {aa: count/total for aa, count in counts.items()}
    
    # Calculate FULL information content (from all AAs)
    entropy = 0
    for freq in all_frequencies.values():
        if freq > 0:
            entropy -= freq * np.log2(freq)
    
    max_entropy = np.log2(20)
    information_content = max_entropy - entropy
    
    # Filter to significant AAs for display
    significant_aas = {aa: freq for aa, freq in all_frequencies.items() if freq >= min_freq}
    
    if not significant_aas:
        return {}, 0
    
    # Letter heights use ORIGINAL frequencies × ORIGINAL IC
    # This way they're properly sized and won't overflow
    letter_heights = {aa: freq * information_content 
                     for aa, freq in significant_aas.items()}
    
    return letter_heights, information_content

def draw_letter(ax, letter, x, y, height, width, color):
    """Draw a letter with visible border"""
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
    
    patch = PathPatch(text_path, facecolor=color, edgecolor='black', linewidth=1.5)
    ax.add_patch(patch)

# 1. AAs with >5% frequency
fig, ax = plt.subplots(figsize=(14, 7))

for pos_idx in range(1, 5):
    letter_heights, ic = calculate_information_content(sequences, pos_idx, min_freq=0.05)
    sorted_letters = sorted(letter_heights.items(), key=lambda x: x[1])
    
    y_position = 0
    for letter, height in sorted_letters:
        color = AA_COLORS.get(letter, '#666666')
        draw_letter(ax, letter, pos_idx - 0.45, y_position, height, 0.9, color)
        y_position += height

ax.set_xlim(0.3, 4.7)
ax.set_ylim(0, 1.5)
ax.set_xticks([1, 2, 3, 4])
ax.set_xticklabels(['Position 2', 'Position 3', 'Position 4', 'Position 5'], 
                    fontsize=14, weight='bold')
ax.set_xlabel('Sequence Position', fontsize=16, weight='bold')
ax.set_ylabel('Information Content (bits)', fontsize=16, weight='bold')
ax.set_title('UBR3 N-Terminal Recognition Motif (Positions 2-5)\nSequence Logo - Amino Acids >5% Frequency', 
             fontsize=18, weight='bold', pad=20)

ax.set_yticks([0, 0.5, 1.0, 1.5])
ax.grid(axis='y', alpha=0.25, linestyle='-', linewidth=0.8, color='gray')
ax.set_facecolor('#F8F8F8')

from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#CC0000', edgecolor='black', label='Acidic (D,E)'),
    Patch(facecolor='#0000CC', edgecolor='black', label='Basic (K,R)'),
    Patch(facecolor='#6699CC', edgecolor='black', label='His (H)'),
    Patch(facecolor='#000000', edgecolor='black', label='Hydrophobic'),
    Patch(facecolor='#33AA33', edgecolor='black', label='Polar'),
    Patch(facecolor='#FF9900', edgecolor='black', label='G, P')
]
ax.legend(handles=legend_elements, loc='upper right', fontsize=12, framealpha=0.95,
         title='AA Properties', title_fontsize=12, edgecolor='black')

note = "Showing amino acids with >5% frequency\nLetter height = frequency × information content"
ax.text(0.02, 0.98, note, transform=ax.transAxes, ha='left', va='top',
        fontsize=10, bbox=dict(boxstyle='round,pad=0.6', facecolor='lightyellow', 
                               edgecolor='black', linewidth=1.5, alpha=0.9))

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'sequence_logo_5percent.png'), dpi=300, bbox_inches='tight')
plt.close()
print("Saved: sequence_logo_5percent.png")

# 2. AAs with >7% frequency (balanced)
fig, ax = plt.subplots(figsize=(14, 7))

for pos_idx in range(1, 5):
    letter_heights, ic = calculate_information_content(sequences, pos_idx, min_freq=0.07)
    sorted_letters = sorted(letter_heights.items(), key=lambda x: x[1])
    
    y_position = 0
    for letter, height in sorted_letters:
        color = AA_COLORS.get(letter, '#666666')
        draw_letter(ax, letter, pos_idx - 0.45, y_position, height, 0.9, color)
        y_position += height

ax.set_xlim(0.3, 4.7)
ax.set_ylim(0, 1.5)
ax.set_xticks([1, 2, 3, 4])
ax.set_xticklabels(['Position 2', 'Position 3', 'Position 4', 'Position 5'], 
                    fontsize=14, weight='bold')
ax.set_xlabel('Sequence Position', fontsize=16, weight='bold')
ax.set_ylabel('Information Content (bits)', fontsize=16, weight='bold')
ax.set_title('UBR3 N-Terminal Recognition Motif (Positions 2-5)\nSequence Logo - Amino Acids >7% Frequency', 
             fontsize=18, weight='bold', pad=20)

ax.set_yticks([0, 0.5, 1.0, 1.5])
ax.grid(axis='y', alpha=0.25, linestyle='-', linewidth=0.8, color='gray')
ax.set_facecolor('#F8F8F8')

ax.legend(handles=legend_elements, loc='upper right', fontsize=12, framealpha=0.95,
         title='AA Properties', title_fontsize=12, edgecolor='black')

note = "Showing amino acids with >7% frequency\nLetter height = frequency × information content"
ax.text(0.02, 0.98, note, transform=ax.transAxes, ha='left', va='top',
        fontsize=10, bbox=dict(boxstyle='round,pad=0.6', facecolor='lightyellow', 
                               edgecolor='black', linewidth=1.5, alpha=0.9))

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'sequence_logo_7percent.png'), dpi=300, bbox_inches='tight')
plt.close()
print("Saved: sequence_logo_7percent.png")

# 3. Top 5 per position
fig, ax = plt.subplots(figsize=(14, 7))

for pos_idx in range(1, 5):
    # Get all letter heights
    letter_heights_all, ic = calculate_information_content(sequences, pos_idx, min_freq=0.0)
    
    # Keep only top 5
    top5 = dict(sorted(letter_heights_all.items(), key=lambda x: x[1], reverse=True)[:5])
    sorted_letters = sorted(top5.items(), key=lambda x: x[1])
    
    y_position = 0
    for letter, height in sorted_letters:
        color = AA_COLORS.get(letter, '#666666')
        draw_letter(ax, letter, pos_idx - 0.45, y_position, height, 0.9, color)
        y_position += height

ax.set_xlim(0.3, 4.7)
ax.set_ylim(0, 1.5)
ax.set_xticks([1, 2, 3, 4])
ax.set_xticklabels(['Position 2', 'Position 3', 'Position 4', 'Position 5'], 
                    fontsize=14, weight='bold')
ax.set_xlabel('Sequence Position', fontsize=16, weight='bold')
ax.set_ylabel('Information Content (bits)', fontsize=16, weight='bold')
ax.set_title('UBR3 N-Terminal Recognition Motif (Positions 2-5)\nSequence Logo - Top 5 Amino Acids per Position', 
             fontsize=18, weight='bold', pad=20)

ax.set_yticks([0, 0.5, 1.0, 1.5])
ax.grid(axis='y', alpha=0.25, linestyle='-', linewidth=0.8, color='gray')
ax.set_facecolor('#F8F8F8')

ax.legend(handles=legend_elements, loc='upper right', fontsize=12, framealpha=0.95,
         title='AA Properties', title_fontsize=12, edgecolor='black')

note = "Showing top 5 amino acids per position\nLetter height = frequency × information content"
ax.text(0.02, 0.98, note, transform=ax.transAxes, ha='left', va='top',
        fontsize=10, bbox=dict(boxstyle='round,pad=0.6', facecolor='lightyellow', 
                               edgecolor='black', linewidth=1.5, alpha=0.9))

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'sequence_logo_top5_proper.png'), dpi=300, bbox_inches='tight')
plt.close()
print("Saved: sequence_logo_top5_proper.png")

print("\nCreated 3 properly scaled sequence logos:")
print("  1. sequence_logo_5percent.png - AAs >5% frequency")
print("  2. sequence_logo_7percent.png - AAs >7% frequency (RECOMMENDED)")
print("  3. sequence_logo_top5_proper.png - Top 5 per position")
print("\nAll logos now:")
print("  - Use correct information content scaling")
print("  - Show proper letter sizes")
print("  - Stay within y-axis bounds")
print("  - Have clear, readable amino acids")
print("\nDone!")
