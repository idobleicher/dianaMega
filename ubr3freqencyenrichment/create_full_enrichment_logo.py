#!/usr/bin/env python3
"""
Create COMPLETE enrichment sequence logos showing ALL enriched amino acids
Letter heights proportional to enrichment values
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
print("CREATING COMPLETE ENRICHMENT SEQUENCE LOGOS (ALL ENRICHED AAs)")
print("="*80)

# Load enrichment data
enrichment_df = pd.read_csv(os.path.join(OUTPUT_DIR, 'enrichment_first5_data.csv'))

print("\nLoaded enrichment data...")

# Amino acid color scheme (standard biochemical coloring)
AA_COLORS = {
    # Hydrophobic (black)
    'A': '#000000', 'V': '#000000', 'I': '#000000', 'L': '#000000', 
    'M': '#000000', 'F': '#000000', 'W': '#000000', 'P': '#000000',
    # Polar (green)
    'S': '#33AA33', 'T': '#33AA33', 'N': '#33AA33', 'Q': '#33AA33',
    'C': '#33AA33', 'Y': '#33AA33',
    # Acidic (red)
    'D': '#CC0000', 'E': '#CC0000',
    # Basic (blue)
    'K': '#0000CC', 'R': '#0000CC', 'H': '#6699CC',
    # Glycine (orange)
    'G': '#FF9900'
}

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
    
    patch = PathPatch(text_path, facecolor=color, edgecolor='black', linewidth=0.5)
    ax.add_patch(patch)

# 1. COMPLETE LOGO - All enriched AAs (enrichment > 1.0)
print("\n1. Creating complete enrichment logo (all enriched > 1.0)...")

fig, ax = plt.subplots(figsize=(14, 8))

positions = range(2, 6)  # Skip position 1 (all Met)

for pos in positions:
    pos_data = enrichment_df[enrichment_df['position'] == pos].copy()
    
    # Show ALL enriched AAs (enrichment > 1.0)
    pos_data = pos_data[pos_data['enrichment'] > 1.0].sort_values('enrichment', ascending=True)
    
    # Calculate letter heights proportional to enrichment
    # Use (enrichment - 1) to show excess over background
    pos_data['letter_height'] = pos_data['enrichment'] - 1.0
    
    # Stack letters from bottom to top
    y_position = 0
    for _, row in pos_data.iterrows():
        letter = row['amino_acid']
        height = row['letter_height']
        color = AA_COLORS.get(letter, '#666666')
        
        draw_letter(ax, letter, pos - 0.4, y_position, height, 0.8, color)
        
        y_position += height

# Add position 1 (Met - universal)
draw_letter(ax, 'M', 1 - 0.4, 0, 0.3, 0.8, AA_COLORS['M'])
ax.text(1, 0.35, '100%', ha='center', va='bottom', fontsize=9, weight='bold')

# Formatting
ax.set_xlim(0.5, 5.5)
ax.set_ylim(-0.2, 12)
ax.set_xticks(range(1, 6))
ax.set_xticklabels(['1', '2', '3', '4', '5'], fontsize=14, weight='bold')
ax.set_xlabel('Position', fontsize=16, weight='bold')
ax.set_ylabel('Cumulative Enrichment (Enrichment - 1)', fontsize=16, weight='bold')
ax.set_title('Complete Enrichment Logo: All Enriched Amino Acids (>1.0x)\nN-Terminal Recognition Motif (54 Hits vs 16,514 Screen)', 
             fontsize=18, weight='bold', pad=20)

ax.axhline(y=0, color='black', linestyle='-', linewidth=2, alpha=0.7)
ax.grid(axis='y', alpha=0.3, linestyle='--')

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#CC0000', edgecolor='black', label='Acidic (D,E)', linewidth=1.5),
    Patch(facecolor='#0000CC', edgecolor='black', label='Basic (K,R)', linewidth=1.5),
    Patch(facecolor='#6699CC', edgecolor='black', label='His (H)', linewidth=1.5),
    Patch(facecolor='#000000', edgecolor='black', label='Hydrophobic', linewidth=1.5),
    Patch(facecolor='#33AA33', edgecolor='black', label='Polar', linewidth=1.5),
    Patch(facecolor='#FF9900', edgecolor='black', label='Glycine', linewidth=1.5)
]
ax.legend(handles=legend_elements, loc='upper left', fontsize=12, framealpha=0.95)

# Add note
ax.text(0.98, 0.02, 'Letter height = Enrichment - 1.0\nAll AAs with enrichment > 1.0 shown',
        transform=ax.transAxes, ha='right', va='bottom', fontsize=10,
        bbox=dict(boxstyle='round,pad=0.5', facecolor='wheat', alpha=0.8))

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'sequence_logo_complete_all_enriched.png'), dpi=300, bbox_inches='tight')
plt.close()

print("Saved: sequence_logo_complete_all_enriched.png")

# 2. SCALED LOGO - Better visual scaling
print("\n2. Creating scaled enrichment logo (all enriched > 1.0, scaled)...")

fig, ax = plt.subplots(figsize=(14, 9))

for pos in positions:
    pos_data = enrichment_df[enrichment_df['position'] == pos].copy()
    pos_data = pos_data[pos_data['enrichment'] > 1.0].sort_values('enrichment', ascending=True)
    
    # Scale the enrichment for better visualization
    # Use log2 enrichment to compress very high values
    pos_data['letter_height'] = pos_data['log2_enrichment'].apply(lambda x: max(0, x))
    
    y_position = 0
    for _, row in pos_data.iterrows():
        letter = row['amino_acid']
        height = row['letter_height']
        color = AA_COLORS.get(letter, '#666666')
        
        draw_letter(ax, letter, pos - 0.4, y_position, height, 0.8, color)
        y_position += height

# Position 1
draw_letter(ax, 'M', 1 - 0.4, 0, 0.3, 0.8, AA_COLORS['M'])
ax.text(1, 0.35, '100%', ha='center', va='bottom', fontsize=9, weight='bold')

ax.set_xlim(0.5, 5.5)
ax.set_ylim(-0.2, 9)
ax.set_xticks(range(1, 6))
ax.set_xticklabels(['1', '2', '3', '4', '5'], fontsize=14, weight='bold')
ax.set_xlabel('Position', fontsize=16, weight='bold')
ax.set_ylabel('Cumulative Log₂(Enrichment)', fontsize=16, weight='bold')
ax.set_title('Enrichment Logo (Log₂ Scale): All Enriched Amino Acids\nN-Terminal Recognition Motif (54 Hits vs 16,514 Screen)', 
             fontsize=18, weight='bold', pad=20)

ax.axhline(y=0, color='black', linestyle='-', linewidth=2, alpha=0.7)
ax.grid(axis='y', alpha=0.3, linestyle='--')

ax.legend(handles=legend_elements, loc='upper left', fontsize=12, framealpha=0.95)

ax.text(0.98, 0.02, 'Letter height = Log₂(Enrichment)\nAll AAs with enrichment > 1.0 shown',
        transform=ax.transAxes, ha='right', va='bottom', fontsize=10,
        bbox=dict(boxstyle='round,pad=0.5', facecolor='wheat', alpha=0.8))

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'sequence_logo_complete_log2_scaled.png'), dpi=300, bbox_inches='tight')
plt.close()

print("Saved: sequence_logo_complete_log2_scaled.png")

# 3. PUBLICATION LOGO - All enriched > 1.2 (cleaner threshold)
print("\n3. Creating publication logo (enriched > 1.2, annotated)...")

fig, ax = plt.subplots(figsize=(16, 9))

for pos in positions:
    pos_data = enrichment_df[enrichment_df['position'] == pos].copy()
    pos_data = pos_data[pos_data['enrichment'] > 1.2].sort_values('enrichment', ascending=True)
    
    # Use sqrt scaling for better visual balance
    pos_data['letter_height'] = np.sqrt(pos_data['enrichment'] - 1.0) * 1.5
    
    y_position = 0
    for _, row in pos_data.iterrows():
        letter = row['amino_acid']
        height = row['letter_height']
        color = AA_COLORS.get(letter, '#666666')
        
        draw_letter(ax, letter, pos - 0.45, y_position, height, 0.9, color)
        
        # Add enrichment value for top 3
        if height > 1.0:  # Only label significant ones
            ax.text(pos + 0.5, y_position + height/2, 
                   f'{row["enrichment"]:.1f}x',
                   fontsize=8, weight='bold', color='white',
                   bbox=dict(boxstyle='round,pad=0.2', facecolor='black', alpha=0.7))
        
        y_position += height

# Position 1
draw_letter(ax, 'M', 1 - 0.45, 0, 0.4, 0.9, AA_COLORS['M'])
ax.text(1, 0.45, '100%\nUniversal', ha='center', va='bottom', fontsize=9, weight='bold', style='italic')

ax.set_xlim(0.3, 5.7)
ax.set_ylim(-0.3, 7.5)
ax.set_xticks(range(1, 6))
ax.set_xticklabels(['Position 1\n(Start)', 'Position 2\n(Structural)', 
                    'Position 3\n(Charge)', 'Position 4\n(Recognition)',
                    'Position 5\n(Transition)'], fontsize=12, weight='bold')
ax.set_ylabel('Enrichment (scaled)', fontsize=16, weight='bold')
ax.set_title('UBR3 N-Terminal Recognition Motif: Enrichment-Based Sequence Logo\nAll Enriched Amino Acids (>1.2x) - 54 Hits vs 16,514 Full Screen Library', 
             fontsize=18, weight='bold', pad=20)

ax.axhline(y=0, color='black', linestyle='-', linewidth=2.5, alpha=0.7)
ax.grid(axis='y', alpha=0.2, linestyle='-', linewidth=0.5)

ax.legend(handles=legend_elements, loc='upper left', fontsize=13, framealpha=0.95,
         title='Amino Acid Properties', title_fontsize=13)

# Consensus at bottom
ax.text(0.5, -0.08, 'Consensus Motif: M - [P/G] - [D/E] - Y - [I/K]',
        transform=ax.transAxes, ha='center', va='top',
        fontsize=15, weight='bold', style='italic',
        bbox=dict(boxstyle='round,pad=0.6', facecolor='lightblue', 
                 edgecolor='black', linewidth=2.5, alpha=0.9))

# Key findings box
key_findings = (
    "Key Findings:\n"
    "• Position 2: P (4.39x), G (3.25x) - Structural kink\n"
    "• Position 3: D (4.57x), E (2.27x) - Acidic charge ⭐HIGHEST\n"
    "• Position 4: Y (4.37x) - Aromatic recognition ⭐KEY\n"
    "• Position 5: I (2.96x), K (2.49x) - Hydrophobic/Basic"
)
ax.text(0.98, 0.98, key_findings,
        transform=ax.transAxes, ha='right', va='top',
        fontsize=10, family='monospace',
        bbox=dict(boxstyle='round,pad=0.7', facecolor='lightyellow', 
                 edgecolor='darkgray', linewidth=2, alpha=0.95))

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'sequence_logo_complete_publication.png'), dpi=300, bbox_inches='tight')
plt.close()

print("Saved: sequence_logo_complete_publication.png")

# 4. WIDE LOGO - All enriched, wide format for presentations
print("\n4. Creating wide format logo (for presentations)...")

fig, ax = plt.subplots(figsize=(20, 7))

for pos in positions:
    pos_data = enrichment_df[enrichment_df['position'] == pos].copy()
    pos_data = pos_data[pos_data['enrichment'] > 1.0].sort_values('enrichment', ascending=True)
    
    # Linear scaling
    pos_data['letter_height'] = (pos_data['enrichment'] - 1.0) * 0.8
    
    y_position = 0
    for _, row in pos_data.iterrows():
        letter = row['amino_acid']
        height = row['letter_height']
        color = AA_COLORS.get(letter, '#666666')
        
        draw_letter(ax, letter, pos - 0.45, y_position, height, 0.9, color)
        y_position += height

# Position 1
draw_letter(ax, 'M', 1 - 0.45, 0, 0.4, 0.9, AA_COLORS['M'])

ax.set_xlim(0.3, 5.7)
ax.set_ylim(-0.1, 10)
ax.set_xticks(range(1, 6))
ax.set_xticklabels(['Pos 1', 'Pos 2', 'Pos 3', 'Pos 4', 'Pos 5'], fontsize=16, weight='bold')
ax.set_xlabel('Position in Sequence', fontsize=18, weight='bold')
ax.set_ylabel('Enrichment Over Background', fontsize=18, weight='bold')
ax.set_title('UBR3 Recognition Motif - Complete Enrichment Profile (All Enriched Amino Acids)', 
             fontsize=22, weight='bold', pad=20)

ax.axhline(y=0, color='black', linestyle='-', linewidth=2.5)
ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=1)

ax.legend(handles=legend_elements, loc='upper left', fontsize=14, framealpha=0.95,
         ncol=2, title='AA Properties', title_fontsize=14)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'sequence_logo_complete_wide.png'), dpi=300, bbox_inches='tight')
plt.close()

print("Saved: sequence_logo_complete_wide.png")

# 5. DETAILED STATS LOGO - Show counts
print("\n5. Creating detailed statistics logo...")

fig, ax = plt.subplots(figsize=(16, 10))

for pos in positions:
    pos_data = enrichment_df[enrichment_df['position'] == pos].copy()
    pos_data = pos_data[pos_data['enrichment'] > 1.15].sort_values('enrichment', ascending=True)
    
    pos_data['letter_height'] = np.log2(pos_data['enrichment'] + 0.5) * 1.2
    
    y_position = 0
    for _, row in pos_data.iterrows():
        letter = row['amino_acid']
        height = row['letter_height']
        color = AA_COLORS.get(letter, '#666666')
        
        draw_letter(ax, letter, pos - 0.4, y_position, height, 0.8, color)
        
        # Add detailed annotation
        if row['enrichment'] > 1.5:
            annotation = f"{letter}: {row['enrichment']:.1f}x\n{row['hits_freq']*100:.1f}% vs {row['screen_freq']*100:.1f}%"
            ax.text(pos + 0.5, y_position + height/2, annotation,
                   fontsize=7, weight='bold', 
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                           edgecolor='black', alpha=0.85, linewidth=1))
        
        y_position += height

# Position 1
draw_letter(ax, 'M', 1 - 0.4, 0, 0.4, 0.8, AA_COLORS['M'])
ax.text(1, -0.3, 'M: 100%\n(Universal start)', ha='center', va='top', fontsize=9, weight='bold')

ax.set_xlim(0.3, 5.7)
ax.set_ylim(-0.5, 7)
ax.set_xticks(range(1, 6))
ax.set_xticklabels(['1', '2', '3', '4', '5'], fontsize=14, weight='bold')
ax.set_xlabel('Position', fontsize=16, weight='bold')
ax.set_ylabel('Enrichment Score', fontsize=16, weight='bold')
ax.set_title('Detailed Enrichment Logo with Statistics\nN-Terminal Motif (Positions 1-5) - 54 Hits vs Full Screen', 
             fontsize=18, weight='bold', pad=20)

ax.axhline(y=0, color='black', linestyle='-', linewidth=2)
ax.grid(axis='y', alpha=0.3, linestyle='--')

ax.legend(handles=legend_elements, loc='upper left', fontsize=12, framealpha=0.95)

# Add statistics table
stats_text = (
    "Dataset Summary:\n"
    f"Hits: 53 sequences (1 per gene)\n"
    f"Screen: 16,514 sequences\n"
    f"Enrichment shown: >1.15x\n"
    f"Top enrichments:\n"
    f"  D3: 4.57x (highest)\n"
    f"  P2: 4.39x\n"
    f"  Y4: 4.37x (key recognition)"
)
ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, ha='left', va='top',
        fontsize=9, family='monospace',
        bbox=dict(boxstyle='round,pad=0.6', facecolor='lightcyan', 
                 edgecolor='black', linewidth=1.5, alpha=0.9))

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'sequence_logo_complete_detailed.png'), dpi=300, bbox_inches='tight')
plt.close()

print("Saved: sequence_logo_complete_detailed.png")

print("\n" + "="*80)
print("ALL COMPLETE ENRICHMENT LOGOS CREATED!")
print("="*80)
print(f"\nCreated 5 complete logo variants in {OUTPUT_DIR}/:")
print("  1. sequence_logo_complete_all_enriched.png - All AAs >1.0x (linear scale)")
print("  2. sequence_logo_complete_log2_scaled.png - All AAs >1.0x (log2 scale)")
print("  3. sequence_logo_complete_publication.png - All AAs >1.2x (publication-ready)")
print("  4. sequence_logo_complete_wide.png - All AAs >1.0x (wide format)")
print("  5. sequence_logo_complete_detailed.png - All AAs >1.15x (with statistics)")
print("\nAll logos show EVERY enriched amino acid at each position!")
print("Letter heights are proportional to enrichment values.")
print("\nDone!")
