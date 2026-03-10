#!/usr/bin/env python3
"""
Create sequence logo-style plot with letter heights based on enrichment
Similar to WebLogo but using enrichment values instead of information content
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
import os

OUTPUT_DIR = "logo_results_54hits_first5"

print("="*80)
print("CREATING ENRICHMENT SEQUENCE LOGO")
print("="*80)

# Load enrichment data
enrichment_df = pd.read_csv(os.path.join(OUTPUT_DIR, 'enrichment_first5_data.csv'))

print("\nLoaded enrichment data...")

# Amino acid color scheme (standard biochemical coloring)
AA_COLORS = {
    # Hydrophobic (black/gray)
    'A': '#000000', 'V': '#000000', 'I': '#000000', 'L': '#000000', 
    'M': '#000000', 'F': '#000000', 'W': '#000000', 'P': '#000000',
    # Polar (green)
    'S': '#33AA33', 'T': '#33AA33', 'N': '#33AA33', 'Q': '#33AA33',
    'C': '#33AA33', 'Y': '#33AA33',
    # Acidic (red)
    'D': '#CC0000', 'E': '#CC0000',
    # Basic (blue)
    'K': '#0000CC', 'R': '#0000CC', 'H': '#6699CC',
    # Glycine (orange/brown)
    'G': '#FF9900'
}

def draw_letter(ax, letter, x, y, height, width, color):
    """Draw a letter using matplotlib TextPath"""
    # Create text path
    fp = FontProperties(family='Arial', weight='bold', size=100)
    text_path = TextPath((0, 0), letter, size=1, prop=fp)
    
    # Get bounds
    bbox = text_path.get_extents()
    letter_width = bbox.width
    letter_height = bbox.height
    
    # Scale and position
    scale_x = width / letter_width
    scale_y = height / letter_height
    
    # Create transformation matrix
    from matplotlib.transforms import Affine2D
    t = Affine2D().scale(scale_x, scale_y).translate(x, y)
    
    # Transform path
    text_path = text_path.transformed(t)
    
    # Draw
    patch = PathPatch(text_path, facecolor=color, edgecolor='black', linewidth=0.5)
    ax.add_patch(patch)

def create_enrichment_logo(enrichment_df, output_file, title="Enrichment Logo"):
    """Create sequence logo with letter heights based on enrichment"""
    
    print(f"\nCreating {output_file}...")
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    positions = sorted(enrichment_df['position'].unique())
    
    for pos in positions:
        if pos == 1:
            # Position 1 is all Met, skip or show differently
            continue
        
        # Get data for this position
        pos_data = enrichment_df[enrichment_df['position'] == pos].copy()
        
        # Only show enriched AAs (enrichment > 1.2)
        pos_data = pos_data[pos_data['enrichment'] > 1.2].sort_values('enrichment', ascending=True)
        
        if len(pos_data) == 0:
            continue
        
        # Calculate heights
        # Scale enrichment so that max ~= 1.0 for display
        # Use (enrichment - 1) to show relative enrichment
        pos_data['letter_height'] = (pos_data['enrichment'] - 1) * 0.8
        
        # Stack letters
        y_position = 0
        for _, row in pos_data.iterrows():
            letter = row['amino_acid']
            height = row['letter_height']
            color = AA_COLORS.get(letter, '#666666')
            
            # Draw letter
            draw_letter(ax, letter, pos - 0.4, y_position, height, 0.8, color)
            
            y_position += height
    
    # Add position 1 (Met) separately
    draw_letter(ax, 'M', 1 - 0.4, 0, 0.5, 0.8, AA_COLORS['M'])
    ax.text(1, -0.3, '(100%)', ha='center', va='top', fontsize=9, style='italic')
    
    # Set axis properties
    ax.set_xlim(0.5, 5.5)
    ax.set_ylim(-0.4, 3.5)
    ax.set_xticks(range(1, 6))
    ax.set_xticklabels(['1\n(M)', '2', '3', '4', '5'], fontsize=12)
    ax.set_xlabel('Position', fontsize=14, weight='bold')
    ax.set_ylabel('Enrichment (scaled)', fontsize=14, weight='bold')
    ax.set_title(title, fontsize=16, weight='bold', pad=20)
    
    # Add horizontal line at 0
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1.5, alpha=0.5)
    
    # Add grid
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Add legend for colors
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#000000', edgecolor='black', label='Hydrophobic'),
        Patch(facecolor='#33AA33', edgecolor='black', label='Polar'),
        Patch(facecolor='#CC0000', edgecolor='black', label='Acidic'),
        Patch(facecolor='#0000CC', edgecolor='black', label='Basic'),
        Patch(facecolor='#FF9900', edgecolor='black', label='Glycine')
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=10)
    
    # Add note
    ax.text(0.98, 0.02, 'Letter height = (Enrichment - 1) × 0.8\nOnly showing enrichment > 1.2x',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=8,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved: {output_file}")

# Create multiple logo variants

# 1. Standard enrichment logo
create_enrichment_logo(
    enrichment_df,
    os.path.join(OUTPUT_DIR, 'sequence_logo_enrichment_first5.png'),
    "N-Terminal Recognition Motif: Enrichment Logo (54 Hits)"
)

# 2. Create a simplified version showing only top 3 per position
def create_simplified_logo(enrichment_df, output_file):
    """Create simplified logo with only top 3 AAs per position"""
    
    print(f"\nCreating {output_file}...")
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    positions = range(2, 6)  # Skip position 1
    
    for pos in positions:
        pos_data = enrichment_df[enrichment_df['position'] == pos].nlargest(3, 'enrichment')
        
        # Scale heights
        pos_data = pos_data.copy()
        pos_data['letter_height'] = (pos_data['enrichment'] - 1) * 0.8
        
        # Stack letters
        y_position = 0
        for _, row in pos_data.iterrows():
            letter = row['amino_acid']
            height = row['letter_height']
            color = AA_COLORS.get(letter, '#666666')
            
            draw_letter(ax, letter, pos - 0.4, y_position, height, 0.8, color)
            
            # Add enrichment value
            ax.text(pos, y_position + height/2, f'{row["enrichment"]:.1f}x',
                   ha='center', va='center', fontsize=8, weight='bold',
                   color='white', 
                   bbox=dict(boxstyle='round,pad=0.2', facecolor='black', alpha=0.6))
            
            y_position += height
    
    # Add position 1
    draw_letter(ax, 'M', 1 - 0.4, 0, 0.5, 0.8, AA_COLORS['M'])
    ax.text(1, -0.3, 'Universal\nstart', ha='center', va='top', fontsize=8, style='italic')
    
    ax.set_xlim(0.5, 5.5)
    ax.set_ylim(-0.5, 3.5)
    ax.set_xticks(range(1, 6))
    ax.set_xticklabels(['1', '2', '3', '4', '5'], fontsize=13, weight='bold')
    ax.set_xlabel('Position', fontsize=14, weight='bold')
    ax.set_ylabel('Enrichment (scaled)', fontsize=14, weight='bold')
    ax.set_title('UBR3 Recognition Motif: Top 3 Enriched AAs per Position', 
                 fontsize=16, weight='bold', pad=20)
    
    ax.axhline(y=0, color='black', linestyle='-', linewidth=1.5, alpha=0.5)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#000000', edgecolor='black', label='Hydrophobic (I,L,M,V,P)'),
        Patch(facecolor='#33AA33', edgecolor='black', label='Polar (S,T,N,Q,C,Y)'),
        Patch(facecolor='#CC0000', edgecolor='black', label='Acidic (D,E)'),
        Patch(facecolor='#0000CC', edgecolor='black', label='Basic (K,R,H)'),
        Patch(facecolor='#FF9900', edgecolor='black', label='Glycine (G)')
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=9)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved: {output_file}")

create_simplified_logo(
    enrichment_df,
    os.path.join(OUTPUT_DIR, 'sequence_logo_top3_first5.png')
)

# 3. Create publication-quality logo
def create_publication_logo(enrichment_df, output_file):
    """Create clean publication-ready logo"""
    
    print(f"\nCreating {output_file}...")
    
    fig, ax = plt.subplots(figsize=(14, 7))
    
    positions = range(2, 6)
    
    # Get max enrichment for scaling
    max_enrichment = enrichment_df[enrichment_df['position'] > 1]['enrichment'].max()
    
    for pos in positions:
        # Top 4 per position
        pos_data = enrichment_df[enrichment_df['position'] == pos].nlargest(4, 'enrichment')
        pos_data = pos_data[pos_data['enrichment'] > 1.3].copy()
        
        # Scale to make display nice
        pos_data['letter_height'] = (pos_data['enrichment'] - 1) * 1.0
        
        y_position = 0
        for idx, row in pos_data.iterrows():
            letter = row['amino_acid']
            height = row['letter_height']
            color = AA_COLORS.get(letter, '#666666')
            
            draw_letter(ax, letter, pos - 0.45, y_position, height, 0.9, color)
            
            y_position += height
    
    # Position 1
    draw_letter(ax, 'M', 1 - 0.45, 0, 0.6, 0.9, AA_COLORS['M'])
    
    # Formatting
    ax.set_xlim(0.3, 5.7)
    ax.set_ylim(-0.1, 4.0)
    ax.set_xticks(range(1, 6))
    ax.set_xticklabels(['Position 1\n(Start)', 'Position 2\n(Kink)', 
                        'Position 3\n(Charge)', 'Position 4\n(Recognition)',
                        'Position 5\n(Core)'], fontsize=11, weight='bold')
    ax.set_ylabel('Relative Enrichment', fontsize=14, weight='bold')
    ax.set_title('UBR3 N-Terminal Recognition Motif\nEnrichment-Based Sequence Logo (54 Hits vs 16,514 Screen)', 
                 fontsize=16, weight='bold', pad=20)
    
    ax.axhline(y=0, color='black', linestyle='-', linewidth=2, alpha=0.7)
    ax.grid(axis='y', alpha=0.2, linestyle='-', linewidth=0.5)
    
    # Add key findings as annotations
    ax.annotate('4.39x', xy=(2, 2.7), fontsize=11, weight='bold', color='red',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))
    ax.annotate('4.57x\nHIGHEST', xy=(3, 2.8), fontsize=11, weight='bold', color='red',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))
    ax.annotate('4.37x\nKEY', xy=(4, 2.7), fontsize=11, weight='bold', color='red',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))
    
    # Color legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#CC0000', edgecolor='black', label='Acidic (D,E)', linewidth=1.5),
        Patch(facecolor='#0000CC', edgecolor='black', label='Basic (K,R,H)', linewidth=1.5),
        Patch(facecolor='#000000', edgecolor='black', label='Hydrophobic', linewidth=1.5),
        Patch(facecolor='#33AA33', edgecolor='black', label='Polar', linewidth=1.5),
        Patch(facecolor='#FF9900', edgecolor='black', label='Glycine', linewidth=1.5)
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=11, framealpha=0.9)
    
    # Add consensus sequence at bottom
    ax.text(0.5, -0.05, 'Consensus: M - [P/G] - [D/E] - Y - [I/K]',
            transform=ax.transAxes, ha='center', va='top',
            fontsize=13, weight='bold', style='italic',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='lightblue', 
                     edgecolor='black', linewidth=2))
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Saved: {output_file}")

create_publication_logo(
    enrichment_df,
    os.path.join(OUTPUT_DIR, 'sequence_logo_publication_first5.png')
)

print("\n" + "="*80)
print("SEQUENCE LOGOS CREATED!")
print("="*80)
print(f"\nCreated 3 logo variants in {OUTPUT_DIR}/:")
print("  1. sequence_logo_enrichment_first5.png - Standard enrichment logo")
print("  2. sequence_logo_top3_first5.png - Simplified (top 3 per position)")
print("  3. sequence_logo_publication_first5.png ⭐ - Publication-ready with annotations")
print("\nDone!")
