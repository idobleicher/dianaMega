#!/usr/bin/env python3
"""
Create traditional sequence logo based on information content (bits)
Similar to WebLogo - the standard way sequence logos are shown in papers
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
print("CREATING TRADITIONAL SEQUENCE LOGO (INFORMATION CONTENT)")
print("="*80)

# Load data
print("\nLoading sequences...")
hits_data = pd.read_csv(os.path.join(OUTPUT_DIR, 'hits_first5_sequences.csv'))
sequences = hits_data['first_5_residues'].tolist()

print(f"Total sequences: {len(sequences)}")
print(f"Sample sequences: {sequences[:5]}")

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
    
    # Get amino acids at this position
    aas = [seq[position] for seq in sequences if position < len(seq)]
    total = len(aas)
    
    if total == 0:
        return {}, 0
    
    # Count frequencies
    counts = Counter(aas)
    frequencies = {aa: count/total for aa, count in counts.items()}
    
    # Calculate Shannon entropy
    entropy = 0
    for freq in frequencies.values():
        if freq > 0:
            entropy -= freq * np.log2(freq)
    
    # Information content (bits)
    # Maximum entropy for 20 amino acids = log2(20) = 4.322 bits
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

print("\n" + "="*80)
print("CREATING LOGOS")
print("="*80)

# 1. STANDARD SEQUENCE LOGO
print("\n1. Creating standard sequence logo (information content)...")

fig, ax = plt.subplots(figsize=(12, 6))

max_ic = 0  # Track maximum information content for y-axis

for pos in range(5):  # Positions 0-4 (displayed as 1-5)
    letter_heights, ic = calculate_information_content(sequences, pos)
    max_ic = max(max_ic, ic)
    
    # Sort letters by height (ascending for stacking)
    sorted_letters = sorted(letter_heights.items(), key=lambda x: x[1])
    
    # Stack letters
    y_position = 0
    for letter, height in sorted_letters:
        color = AA_COLORS.get(letter, '#666666')
        draw_letter(ax, letter, pos + 1 - 0.4, y_position, height, 0.8, color)
        y_position += height

# Formatting
ax.set_xlim(0.5, 5.5)
ax.set_ylim(0, 4.5)
ax.set_xticks(range(1, 6))
ax.set_xticklabels(['1', '2', '3', '4', '5'], fontsize=13, weight='bold')
ax.set_xlabel('Position', fontsize=14, weight='bold')
ax.set_ylabel('Information Content (bits)', fontsize=14, weight='bold')
ax.set_title('Sequence Logo: N-Terminal Motif (54 Hits)\nInformation Content Based', 
             fontsize=16, weight='bold', pad=15)

# Add y-axis grid
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
ax.legend(handles=legend_elements, loc='upper right', fontsize=10, framealpha=0.9)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'sequence_logo_traditional.png'), dpi=300, bbox_inches='tight')
plt.close()

print("Saved: sequence_logo_traditional.png")

# 2. PUBLICATION-QUALITY LOGO
print("\n2. Creating publication-quality sequence logo...")

fig, ax = plt.subplots(figsize=(14, 7))

for pos in range(5):
    letter_heights, ic = calculate_information_content(sequences, pos)
    sorted_letters = sorted(letter_heights.items(), key=lambda x: x[1])
    
    y_position = 0
    for letter, height in sorted_letters:
        color = AA_COLORS.get(letter, '#666666')
        draw_letter(ax, letter, pos + 1 - 0.45, y_position, height, 0.9, color)
        y_position += height

ax.set_xlim(0.3, 5.7)
ax.set_ylim(0, 4.5)
ax.set_xticks(range(1, 6))
ax.set_xticklabels(['Position 1\n(Start)', 'Position 2', 'Position 3', 
                    'Position 4', 'Position 5'], fontsize=12, weight='bold')
ax.set_xlabel('Sequence Position', fontsize=14, weight='bold')
ax.set_ylabel('Information Content (bits)', fontsize=14, weight='bold')
ax.set_title('UBR3 Recognition Motif: Sequence Logo\n53 Selected Hits from Screen', 
             fontsize=16, weight='bold', pad=20)

ax.set_yticks([0, 1, 2, 3, 4])
ax.grid(axis='y', alpha=0.2, linestyle='-', linewidth=0.5)

ax.legend(handles=legend_elements, loc='upper right', fontsize=11, framealpha=0.95,
         title='Amino Acid Properties', title_fontsize=11)

# Add note about information content
note = "Letter height = frequency × information content\nStack height = total information at position (0-4.32 bits)"
ax.text(0.02, 0.98, note, transform=ax.transAxes, ha='left', va='top',
        fontsize=9, bbox=dict(boxstyle='round,pad=0.5', facecolor='wheat', alpha=0.8))

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'sequence_logo_publication.png'), dpi=300, bbox_inches='tight')
plt.close()

print("Saved: sequence_logo_publication.png")

# 3. WIDE FORMAT FOR PRESENTATIONS
print("\n3. Creating wide format sequence logo...")

fig, ax = plt.subplots(figsize=(18, 6))

for pos in range(5):
    letter_heights, ic = calculate_information_content(sequences, pos)
    sorted_letters = sorted(letter_heights.items(), key=lambda x: x[1])
    
    y_position = 0
    for letter, height in sorted_letters:
        color = AA_COLORS.get(letter, '#666666')
        draw_letter(ax, letter, pos + 1 - 0.45, y_position, height, 0.9, color)
        y_position += height

ax.set_xlim(0.3, 5.7)
ax.set_ylim(0, 4.5)
ax.set_xticks(range(1, 6))
ax.set_xticklabels(['Pos 1', 'Pos 2', 'Pos 3', 'Pos 4', 'Pos 5'], 
                    fontsize=16, weight='bold')
ax.set_xlabel('Position', fontsize=18, weight='bold')
ax.set_ylabel('Information (bits)', fontsize=18, weight='bold')
ax.set_title('UBR3 N-Terminal Recognition Motif - Sequence Logo', 
             fontsize=20, weight='bold', pad=15)

ax.set_yticks([0, 1, 2, 3, 4])
ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=1)

ax.legend(handles=legend_elements, loc='upper right', fontsize=13, framealpha=0.95,
         ncol=2, title='AA Properties', title_fontsize=13)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'sequence_logo_wide.png'), dpi=300, bbox_inches='tight')
plt.close()

print("Saved: sequence_logo_wide.png")

# 4. CALCULATE AND DISPLAY STATISTICS
print("\n" + "="*80)
print("SEQUENCE STATISTICS")
print("="*80)

for pos in range(5):
    letter_heights, ic = calculate_information_content(sequences, pos)
    
    print(f"\nPosition {pos + 1}:")
    print(f"  Information Content: {ic:.2f} bits")
    
    # Get top 5 amino acids
    sorted_letters = sorted(letter_heights.items(), key=lambda x: x[1], reverse=True)[:5]
    print(f"  Top amino acids:")
    for letter, height in sorted_letters:
        freq = height / ic if ic > 0 else 0
        print(f"    {letter}: {freq*100:.1f}% (contributes {height:.2f} bits)")

# Save statistics to file
with open(os.path.join(OUTPUT_DIR, 'sequence_logo_statistics.txt'), 'w') as f:
    f.write("SEQUENCE LOGO STATISTICS\n")
    f.write("="*60 + "\n\n")
    f.write(f"Dataset: 53 sequences (54 hits, one per gene)\n")
    f.write(f"Positions analyzed: 1-5 (N-terminal motif)\n\n")
    
    for pos in range(5):
        letter_heights, ic = calculate_information_content(sequences, pos)
        f.write(f"\nPosition {pos + 1}:\n")
        f.write(f"  Information Content: {ic:.3f} bits\n")
        f.write(f"  Conservation: {ic/4.322*100:.1f}% (of maximum)\n")
        f.write(f"  Top amino acids:\n")
        
        sorted_letters = sorted(letter_heights.items(), key=lambda x: x[1], reverse=True)
        for letter, height in sorted_letters[:10]:
            freq = height / ic if ic > 0 else 0
            f.write(f"    {letter}: {freq*100:5.1f}% frequency, {height:.3f} bits\n")

print("\nSaved: sequence_logo_statistics.txt")

print("\n" + "="*80)
print("TRADITIONAL SEQUENCE LOGOS CREATED!")
print("="*80)
print(f"\nCreated files in {OUTPUT_DIR}/:")
print("  1. sequence_logo_traditional.png - Standard sequence logo")
print("  2. sequence_logo_publication.png - Publication-quality version")
print("  3. sequence_logo_wide.png - Wide format for presentations")
print("  4. sequence_logo_statistics.txt - Detailed statistics")
print("\nThese are TRADITIONAL sequence logos based on information content,")
print("as typically shown in scientific papers (like WebLogo).")
print("\nLetter heights = frequency × information content (bits)")
print("Stack heights = total information at each position (0-4.32 bits)")
print("\nDone!")
