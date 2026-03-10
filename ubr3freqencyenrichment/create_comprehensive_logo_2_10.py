#!/usr/bin/env python3
"""
Create comprehensive sequence logo analysis for positions 2-10
Similar to the full analysis but shortened to positions 2-10
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath
from matplotlib.patches import PathPatch
from matplotlib.transforms import Affine2D
import seaborn as sns
import os

OUTPUT_DIR = "logo_results_54hits_first5"

print("Creating comprehensive sequence logo analysis (positions 2-10)...")

# Load data
screen = pd.read_excel('UBR3 Nt screen.xlsx')
hits_file = pd.read_excel('ubr3_best (1).xlsx')

gene_col = hits_file.columns[2]
genes = hits_file[gene_col].unique()[:54]

hits_data = screen[screen['Gene_ID'].isin(genes)].groupby('Gene_ID').first().reset_index()
hits_sequences = hits_data['AA_seq'].tolist()
screen_sequences = screen['AA_seq'].tolist()

# Truncate to position 10
hits_sequences_10 = [seq[:10] for seq in hits_sequences if len(seq) >= 10]
screen_sequences_10 = [seq[:10] for seq in screen_sequences if len(seq) >= 10]

print(f"Hits sequences: {len(hits_sequences_10)}")
print(f"Screen sequences: {len(screen_sequences_10)}")

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
    if height <= 0.001:
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

def get_consensus(sequences):
    from collections import Counter
    consensus = ""
    for pos in range(10):
        aas = [seq[pos] for seq in sequences if pos < len(seq)]
        if aas:
            most_common = Counter(aas).most_common(1)[0][0]
            consensus += most_common
    return consensus

# Create figure with 5 subplots
fig = plt.figure(figsize=(20, 16))
gs = fig.add_gridspec(5, 1, height_ratios=[1, 1, 1.2, 0.8, 0.8], hspace=0.35)

# ============================================================
# PANEL 1: HITS SEQUENCE LOGO
# ============================================================
ax1 = fig.add_subplot(gs[0])
print("\nPanel 1: Hits sequence logo...")

for pos_idx in range(1, 10):
    letter_heights_all, ic = calculate_information_content(hits_sequences_10, pos_idx)
    sorted_letters = sorted(letter_heights_all.items(), key=lambda x: x[1])
    
    y_position = 0
    for letter, height in sorted_letters:
        color = AA_COLORS.get(letter, '#666666')
        draw_letter(ax1, letter, pos_idx - 0.45, y_position, height, 0.9, color)
        y_position += height

ax1.set_xlim(0.3, 9.7)
ax1.set_ylim(0, 0.85)
ax1.set_xticks(range(1, 10))
ax1.set_xticklabels(['2', '3', '4', '5', '6', '7', '8', '9', '10'], fontsize=11)
ax1.set_ylabel('Information\nContent (bits)', fontsize=11, weight='bold')
ax1.set_title('UBR3 Hits - Sequence Logo (Positions 2-10)', fontsize=14, weight='bold', pad=10)
ax1.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
ax1.grid(axis='y', alpha=0.2)

# Add consensus
hits_consensus = "M" + get_consensus(hits_sequences_10)[1:]
ax1.text(0.98, 0.95, f'Consensus: {hits_consensus}', transform=ax1.transAxes, 
         ha='right', va='top', fontsize=10, 
         bbox=dict(boxstyle='round,pad=0.5', facecolor='wheat', edgecolor='black'))

# ============================================================
# PANEL 2: FULL SCREEN SEQUENCE LOGO
# ============================================================
ax2 = fig.add_subplot(gs[1])
print("Panel 2: Full screen sequence logo...")

for pos_idx in range(1, 10):
    letter_heights_all, ic = calculate_information_content(screen_sequences_10, pos_idx)
    sorted_letters = sorted(letter_heights_all.items(), key=lambda x: x[1])
    
    y_position = 0
    for letter, height in sorted_letters:
        color = AA_COLORS.get(letter, '#666666')
        draw_letter(ax2, letter, pos_idx - 0.45, y_position, height, 0.9, color)
        y_position += height

ax2.set_xlim(0.3, 9.7)
ax2.set_ylim(0, 0.45)
ax2.set_xticks(range(1, 10))
ax2.set_xticklabels(['2', '3', '4', '5', '6', '7', '8', '9', '10'], fontsize=11)
ax2.set_ylabel('Information\nContent (bits)', fontsize=11, weight='bold')
ax2.set_title('Full Screen - Sequence Logo (Positions 2-10)', fontsize=14, weight='bold', pad=10)
ax2.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
ax2.grid(axis='y', alpha=0.2)

# Add consensus
screen_consensus = "M" + get_consensus(screen_sequences_10)[1:]
ax2.text(0.98, 0.95, f'Consensus: {screen_consensus}', transform=ax2.transAxes, 
         ha='right', va='top', fontsize=10,
         bbox=dict(boxstyle='round,pad=0.5', facecolor='lightblue', edgecolor='black'))

# ============================================================
# PANEL 3: ENRICHMENT LOGO
# ============================================================
ax3 = fig.add_subplot(gs[2])
print("Panel 3: Enrichment logo...")

# Calculate enrichment
from collections import Counter
enrichment_data = {}

for pos_idx in range(1, 10):
    hits_aas = [seq[pos_idx] for seq in hits_sequences_10 if pos_idx < len(seq)]
    screen_aas = [seq[pos_idx] for seq in screen_sequences_10 if pos_idx < len(seq)]
    
    hits_counts = Counter(hits_aas)
    screen_counts = Counter(screen_aas)
    
    hits_total = len(hits_aas)
    screen_total = len(screen_aas)
    
    all_aas = set(list(hits_counts.keys()) + list(screen_counts.keys()))
    
    for aa in all_aas:
        hits_freq = hits_counts.get(aa, 0) / hits_total if hits_total > 0 else 0
        screen_freq = screen_counts.get(aa, 0) / screen_total if screen_total > 0 else 0
        
        if screen_freq > 0 and hits_freq > 0:
            enrichment = hits_freq / screen_freq
        elif hits_freq > 0:
            enrichment = 10
        else:
            enrichment = 0.1
        
        if pos_idx not in enrichment_data:
            enrichment_data[pos_idx] = {}
        enrichment_data[pos_idx][aa] = enrichment * hits_freq

# Draw enrichment logo
for pos_idx in range(1, 10):
    enrichments = enrichment_data[pos_idx]
    sorted_letters = sorted(enrichments.items(), key=lambda x: x[1])
    
    y_position = 0
    for letter, height in sorted_letters:
        color = AA_COLORS.get(letter, '#666666')
        draw_letter(ax3, letter, pos_idx - 0.45, y_position, height, 0.9, color)
        y_position += height

ax3.set_xlim(0.3, 9.7)
ax3.set_ylim(0, 2.2)
ax3.set_xticks(range(1, 10))
ax3.set_xticklabels(['2', '3', '4', '5', '6', '7', '8', '9', '10'], fontsize=11)
ax3.set_ylabel('Enrichment\nHeight', fontsize=11, weight='bold')
ax3.set_title('Enrichment Logo - Shows Discriminating Features (Hits vs Screen)', 
              fontsize=14, weight='bold', color='darkred', pad=10)
ax3.set_yticks([0, 0.5, 1.0, 1.5, 2.0])
ax3.grid(axis='y', alpha=0.2)

# Add yellow box explanation
ax3.text(0.98, 0.95, '*** KEY PLOT ***\nHeight = Frequency × Enrichment', 
         transform=ax3.transAxes, ha='right', va='top', fontsize=10, weight='bold',
         bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', edgecolor='black', linewidth=2))

# ============================================================
# PANEL 4: LOG2 ENRICHMENT HEATMAP
# ============================================================
ax4 = fig.add_subplot(gs[3])
print("Panel 4: Log2 enrichment heatmap...")

# Calculate log2 enrichment
amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
               'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
log2_matrix = []

for aa in amino_acids:
    row = []
    for pos_idx in range(1, 10):
        hits_aas = [seq[pos_idx] for seq in hits_sequences_10 if pos_idx < len(seq)]
        screen_aas = [seq[pos_idx] for seq in screen_sequences_10 if pos_idx < len(seq)]
        
        hits_freq = hits_aas.count(aa) / len(hits_aas) if hits_aas else 0
        screen_freq = screen_aas.count(aa) / len(screen_aas) if screen_aas else 0
        
        # Add pseudocount
        hits_freq = (hits_aas.count(aa) + 0.5) / (len(hits_aas) + 10)
        screen_freq = (screen_aas.count(aa) + 0.5) / (len(screen_aas) + 10)
        
        log2_enrichment = np.log2(hits_freq / screen_freq)
        row.append(log2_enrichment)
    log2_matrix.append(row)

log2_df = pd.DataFrame(log2_matrix, index=amino_acids, columns=['2', '3', '4', '5', '6', '7', '8', '9', '10'])

sns.heatmap(log2_df, cmap='RdBu_r', center=0, vmin=-1.5, vmax=1.5,
            cbar_kws={'label': 'Log2 Enrichment (Hits/Screen)', 'shrink': 0.8},
            linewidths=0.5, linecolor='white', ax=ax4)
ax4.set_xlabel('Position', fontsize=11, weight='bold')
ax4.set_ylabel('Amino Acid', fontsize=11, weight='bold')
ax4.set_title('Log2 Enrichment Heatmap (Red = Enriched in Hits, Blue = Depleted)', 
              fontsize=13, weight='bold', pad=10)

# ============================================================
# PANEL 5: POSITION DISCRIMINATION SCORE
# ============================================================
ax5 = fig.add_subplot(gs[4])
print("Panel 5: Position discrimination score...")

# Calculate discrimination score (KL divergence)
discrimination_scores = []
for pos_idx in range(1, 10):
    hits_aas = [seq[pos_idx] for seq in hits_sequences_10 if pos_idx < len(seq)]
    screen_aas = [seq[pos_idx] for seq in screen_sequences_10 if pos_idx < len(seq)]
    
    hits_counts = Counter(hits_aas)
    screen_counts = Counter(screen_aas)
    
    hits_total = len(hits_aas)
    screen_total = len(screen_aas)
    
    all_aas = set(amino_acids)
    kl_divergence = 0
    
    for aa in all_aas:
        p_hits = (hits_counts.get(aa, 0) + 0.5) / (hits_total + 10)
        p_screen = (screen_counts.get(aa, 0) + 0.5) / (screen_total + 10)
        kl_divergence += p_hits * np.log2(p_hits / p_screen)
    
    discrimination_scores.append(kl_divergence)

# Find top 5 discriminating positions
top5_positions = sorted(range(len(discrimination_scores)), 
                        key=lambda i: discrimination_scores[i], reverse=True)[:5]

colors = ['#CC3333' if i in top5_positions else '#6699CC' for i in range(len(discrimination_scores))]

bars = ax5.bar(range(1, 10), discrimination_scores, color=colors, edgecolor='black', linewidth=1.5)
ax5.set_xlabel('Position', fontsize=11, weight='bold')
ax5.set_ylabel('Enrichment Score', fontsize=11, weight='bold')
ax5.set_title('Position Discrimination Score (Red = Top 5 Most Discriminating)', 
              fontsize=13, weight='bold', pad=10)
ax5.set_xticks(range(1, 10))
ax5.set_xticklabels(['2', '3', '4', '5', '6', '7', '8', '9', '10'])
ax5.grid(axis='y', alpha=0.3)
ax5.set_facecolor('#F5F5F5')

# Annotate top positions
for i in top5_positions:
    ax5.text(i + 1, discrimination_scores[i] + 0.02, f'Pos {i+2}',
            ha='center', va='bottom', fontsize=9, weight='bold')

# Main title
fig.suptitle('UBR3 Comprehensive Sequence Analysis (Positions 2-10)', 
             fontsize=18, weight='bold', y=0.995)

plt.savefig(os.path.join(OUTPUT_DIR, 'comprehensive_analysis_pos2_10.png'), 
            dpi=300, bbox_inches='tight')
plt.close()

print("\n" + "="*80)
print("COMPREHENSIVE ANALYSIS COMPLETE!")
print("="*80)
print(f"\nSaved: comprehensive_analysis_pos2_10.png")
print("\nThis figure includes:")
print("  1. Hits sequence logo (positions 2-10)")
print("  2. Full screen sequence logo (positions 2-10)")
print("  3. Enrichment logo (positions 2-10)")
print("  4. Log2 enrichment heatmap (positions 2-10)")
print("  5. Position discrimination scores (positions 2-10)")
print("\nDone!")
