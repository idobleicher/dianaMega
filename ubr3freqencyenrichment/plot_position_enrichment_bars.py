#!/usr/bin/env python3
"""
Generate position-specific enrichment bar graphs for UBR3 analysis.
Shows enrichment organized by residue number (position).
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys

# Define amino acid color groups
AA_COLORS = {
    # Hydrophobic (orange/red)
    'A': '#E74C3C', 'V': '#E74C3C', 'L': '#E74C3C', 'I': '#E74C3C', 
    'M': '#E74C3C', 'F': '#E74C3C', 'W': '#E74C3C', 'P': '#E74C3C',
    # Polar (blue)
    'S': '#3498DB', 'T': '#3498DB', 'N': '#3498DB', 'Q': '#3498DB',
    'C': '#3498DB', 'Y': '#3498DB',
    # Charged positive (purple)
    'K': '#9B59B6', 'R': '#9B59B6', 'H': '#9B59B6',
    # Charged negative (green)
    'D': '#2ECC71', 'E': '#2ECC71',
    # Special (gray)
    'G': '#95A5A6'
}

def load_enrichment_data():
    """Load enrichment data for hits vs screen comparison."""
    # Load the position-specific enrichment data from logo results
    enrichment_file = Path('logo_results_enrichment/enrichment_matrix_pos2+.csv')
    
    if not enrichment_file.exists():
        print(f"Error: {enrichment_file} not found!")
        sys.exit(1)
    
    # Load enrichment ratios
    # Data format: rows = positions, columns = amino acids
    enrichment_df = pd.read_csv(enrichment_file, index_col=0)
    
    # Transpose so that columns are positions
    enrichment_df = enrichment_df.T
    
    # Rename columns to have 'Pos' prefix
    enrichment_df.columns = [f'Pos{col}' for col in enrichment_df.columns]
    
    # Also load log2 enrichment
    log2_file = Path('logo_results_enrichment/log2_enrichment_pos2+.csv')
    if log2_file.exists():
        log2_df = pd.read_csv(log2_file, index_col=0).T
        log2_df.columns = [f'Pos{col}' for col in log2_df.columns]
    else:
        log2_df = np.log2(enrichment_df.clip(lower=0.01))  # Avoid log(0)
    
    return enrichment_df, log2_df

def plot_enrichment_by_position(enrichment_df, log2_df, positions_to_show=5, top_n=10):
    """
    Create bar graphs showing enrichment organized by position.
    
    Parameters:
    -----------
    enrichment_df : DataFrame
        Enrichment ratios (hits/screen) for each AA at each position
    log2_df : DataFrame
        Log2 enrichment values
    positions_to_show : int
        Number of positions to display (default: 5 for positions 2-6)
    top_n : int
        Number of top enriched AAs to show per position
    """
    
    # Select positions (columns should be like 'Pos2', 'Pos3', etc.)
    available_positions = [col for col in enrichment_df.columns if col.startswith('Pos')]
    positions = available_positions[:positions_to_show]
    
    if len(positions) == 0:
        print("Error: No position columns found!")
        return
    
    # Create output directory
    output_dir = Path('logo_results_enrichment')
    output_dir.mkdir(exist_ok=True)
    
    # ========== FIGURE 1: Multi-panel Enrichment Ratios ==========
    fig, axes = plt.subplots(1, len(positions), figsize=(4*len(positions), 6))
    if len(positions) == 1:
        axes = [axes]
    
    for idx, pos in enumerate(positions):
        ax = axes[idx]
        
        # Get enrichment data for this position
        pos_data = enrichment_df[pos].dropna()
        
        # Sort and get top N
        pos_data_sorted = pos_data.sort_values(ascending=True)
        top_data = pos_data_sorted.tail(top_n)
        
        # Create colors based on AA properties
        colors = [AA_COLORS.get(aa, '#95A5A6') for aa in top_data.index]
        
        # Plot horizontal bars
        y_pos = np.arange(len(top_data))
        bars = ax.barh(y_pos, top_data.values, color=colors, alpha=0.85, edgecolor='black', linewidth=0.5)
        
        # Add vertical line at 1.0 (no enrichment)
        ax.axvline(x=1.0, color='black', linestyle='--', linewidth=1, alpha=0.5)
        
        # Customize axes
        ax.set_yticks(y_pos)
        ax.set_yticklabels(top_data.index, fontsize=11, fontweight='bold')
        ax.set_xlabel('Enrichment Ratio', fontsize=11, fontweight='bold')
        
        # Extract position number for title
        pos_num = pos.replace('Pos', '')
        ax.set_title(f'Position {pos_num}', fontsize=12, fontweight='bold')
        
        # Grid
        ax.grid(axis='x', alpha=0.3, linestyle=':')
        ax.set_axisbelow(True)
        
        # Set x-axis limits with some padding
        max_val = top_data.max()
        ax.set_xlim(0, max_val * 1.1)
        
        # Add value labels on bars
        for i, (aa, val) in enumerate(zip(top_data.index, top_data.values)):
            ax.text(val + 0.05, i, f'{val:.2f}', va='center', fontsize=9)
    
    plt.suptitle('UBR3 Enrichment Ratios by Position (Hits vs Library)', 
                 fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    output_file = output_dir / 'position_enrichment_ratios.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    sys.stdout.buffer.write(f"[OK] Saved: {output_file}\n".encode('utf-8'))
    sys.stdout.flush()
    plt.close()
    
    # ========== FIGURE 2: Multi-panel Log2 Enrichment ==========
    fig, axes = plt.subplots(1, len(positions), figsize=(4*len(positions), 6))
    if len(positions) == 1:
        axes = [axes]
    
    for idx, pos in enumerate(positions):
        ax = axes[idx]
        
        # Get log2 enrichment data for this position
        pos_data = log2_df[pos].dropna()
        
        # Sort and get top N
        pos_data_sorted = pos_data.sort_values(ascending=True)
        top_data = pos_data_sorted.tail(top_n)
        
        # Create colors: red for enriched (>0), blue for depleted (<0)
        colors = ['#E74C3C' if val > 0 else '#3498DB' for val in top_data.values]
        
        # Plot horizontal bars
        y_pos = np.arange(len(top_data))
        bars = ax.barh(y_pos, top_data.values, color=colors, alpha=0.85, edgecolor='black', linewidth=0.5)
        
        # Add vertical line at 0 (no enrichment)
        ax.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
        
        # Customize axes
        ax.set_yticks(y_pos)
        ax.set_yticklabels(top_data.index, fontsize=11, fontweight='bold')
        ax.set_xlabel('Log2 Enrichment', fontsize=11, fontweight='bold')
        
        # Extract position number for title
        pos_num = pos.replace('Pos', '')
        ax.set_title(f'Position {pos_num}', fontsize=12, fontweight='bold')
        
        # Grid
        ax.grid(axis='x', alpha=0.3, linestyle=':')
        ax.set_axisbelow(True)
        
        # Set x-axis limits symmetrically
        max_abs = abs(top_data).max()
        ax.set_xlim(-max_abs * 0.3, max_abs * 1.1)
        
        # Add value labels on bars
        for i, (aa, val) in enumerate(zip(top_data.index, top_data.values)):
            x_offset = 0.05 if val > 0 else -0.05
            ha = 'left' if val > 0 else 'right'
            ax.text(val + x_offset, i, f'{val:.2f}', va='center', ha=ha, fontsize=9)
    
    plt.suptitle('UBR3 Log2 Enrichment by Position (Hits vs Library)', 
                 fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    
    output_file = output_dir / 'position_log2_enrichment.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    sys.stdout.buffer.write(f"[OK] Saved: {output_file}\n".encode('utf-8'))
    sys.stdout.flush()
    plt.close()
    
    # ========== FIGURE 3: Stacked view - All positions in one figure ==========
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Combine data from all positions
    all_data = []
    all_labels = []
    all_colors = []
    
    y_offset = 0
    for pos in positions:
        pos_data = enrichment_df[pos].dropna()
        pos_data_sorted = pos_data.sort_values(ascending=False)
        top_data = pos_data_sorted.head(top_n)
        
        pos_num = pos.replace('Pos', '')
        
        for aa, val in top_data.items():
            all_data.append(val)
            all_labels.append(f"{aa} (P{pos_num})")
            all_colors.append(AA_COLORS.get(aa, '#95A5A6'))
    
    # Create horizontal bars
    y_pos = np.arange(len(all_data))
    ax.barh(y_pos, all_data, color=all_colors, alpha=0.85, edgecolor='black', linewidth=0.5)
    
    # Add vertical line at 1.0
    ax.axvline(x=1.0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
    
    # Customize axes
    ax.set_yticks(y_pos)
    ax.set_yticklabels(all_labels, fontsize=9)
    ax.set_xlabel('Enrichment Ratio (Hits/Library)', fontsize=12, fontweight='bold')
    ax.set_title('UBR3 Top Enriched Amino Acids by Position', fontsize=14, fontweight='bold')
    
    # Grid
    ax.grid(axis='x', alpha=0.3, linestyle=':')
    ax.set_axisbelow(True)
    
    plt.tight_layout()
    
    output_file = output_dir / 'position_enrichment_stacked.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    sys.stdout.buffer.write(f"[OK] Saved: {output_file}\n".encode('utf-8'))
    sys.stdout.flush()
    plt.close()
    
    # ========== FIGURE 4: Heatmap view of enrichment ==========
    fig, ax = plt.subplots(figsize=(8, 10))
    
    # Prepare data for heatmap: all amino acids x all positions
    heatmap_data = enrichment_df[positions].copy()
    
    # Sort amino acids by overall enrichment
    overall_enrichment = heatmap_data.mean(axis=1).sort_values(ascending=False)
    heatmap_data = heatmap_data.loc[overall_enrichment.index]
    
    # Create heatmap
    sns.heatmap(heatmap_data, cmap='RdBu_r', center=1.0, vmin=0, vmax=3,
                annot=True, fmt='.2f', cbar_kws={'label': 'Enrichment Ratio'},
                linewidths=0.5, linecolor='gray', ax=ax)
    
    # Customize
    pos_labels = [pos.replace('Pos', 'Position ') for pos in positions]
    ax.set_xticklabels(pos_labels, rotation=45, ha='right')
    ax.set_ylabel('Amino Acid', fontsize=12, fontweight='bold')
    ax.set_title('UBR3 Enrichment Heatmap (Hits vs Library)', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    
    output_file = output_dir / 'position_enrichment_heatmap_detailed.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    sys.stdout.buffer.write(f"[OK] Saved: {output_file}\n".encode('utf-8'))
    sys.stdout.flush()
    plt.close()

def main():
    """Main execution function."""
    print("\n" + "="*60)
    print("  UBR3 POSITION-SPECIFIC ENRICHMENT BAR GRAPHS")
    print("="*60 + "\n")
    
    # Load data
    print("Loading enrichment data...")
    enrichment_df, log2_df = load_enrichment_data()
    
    sys.stdout.buffer.write(f"[OK] Loaded enrichment data for {len(enrichment_df)} amino acids\n".encode('utf-8'))
    sys.stdout.buffer.write(f"[OK] Positions available: {', '.join(enrichment_df.columns)}\n\n".encode('utf-8'))
    sys.stdout.flush()
    
    # Generate plots
    print("Generating position-specific enrichment bar graphs...\n")
    plot_enrichment_by_position(enrichment_df, log2_df, positions_to_show=5, top_n=10)
    
    print("\n" + "="*60)
    sys.stdout.buffer.write("  [OK] ANALYSIS COMPLETE!\n".encode('utf-8'))
    sys.stdout.flush()
    print("="*60)
    print("\nGenerated files in 'logo_results_enrichment/':")
    print("  1. position_enrichment_ratios.png - Multi-panel enrichment ratios")
    print("  2. position_log2_enrichment.png - Multi-panel log2 enrichment")
    print("  3. position_enrichment_stacked.png - All positions stacked")
    print("  4. position_enrichment_heatmap_detailed.png - Heatmap overview")
    print("\n")

if __name__ == "__main__":
    main()
