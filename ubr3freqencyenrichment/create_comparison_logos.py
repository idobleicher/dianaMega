#!/usr/bin/env python3
"""
Create side-by-side comparison of logo plots between hits and full screen
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

try:
    import logomaker
    LOGOMAKER_AVAILABLE = True
except ImportError:
    LOGOMAKER_AVAILABLE = False
    print("Warning: logomaker not installed")


def create_comparison_plot():
    """Create comparison visualization"""
    
    print("Creating comparison visualization...")
    
    # Load PWMs
    pwm_hits = pd.read_csv('logo_results_hits/position_weight_matrix.csv', index_col=0)
    pwm_screen = pd.read_csv('logo_results_full_screen/position_weight_matrix.csv', index_col=0)
    
    # Load IC matrices
    ic_hits = pd.read_csv('logo_results_hits/information_content_matrix.csv', index_col=0)
    ic_screen = pd.read_csv('logo_results_full_screen/information_content_matrix.csv', index_col=0)
    
    if LOGOMAKER_AVAILABLE:
        # Create figure with logomaker
        fig, axes = plt.subplots(4, 1, figsize=(18, 16))
        
        # Plot 1: Hits IC logo
        ax1 = axes[0]
        logo1 = logomaker.Logo(ic_hits, ax=ax1, color_scheme='chemistry')
        ax1.set_ylabel('Information\nContent (bits)', fontsize=11)
        ax1.set_title('UBR3 Hits - Enriched Sequences (91 sequences)', 
                     fontsize=14, fontweight='bold', pad=10)
        ax1.set_xlim([0.5, len(ic_hits) + 0.5])
        ax1.text(0.98, 0.95, 'Consensus: MPDLSVLLLSLSGTGGGSSGGTRL', 
                transform=ax1.transAxes, fontsize=10, verticalalignment='top',
                horizontalalignment='right', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # Plot 2: Full screen IC logo
        ax2 = axes[1]
        logo2 = logomaker.Logo(ic_screen, ax=ax2, color_scheme='chemistry')
        ax2.set_ylabel('Information\nContent (bits)', fontsize=11)
        ax2.set_title('Full UBR3 Screen - All Library (16,514 sequences)', 
                     fontsize=14, fontweight='bold', pad=10)
        ax2.set_xlim([0.5, len(ic_screen) + 0.5])
        ax2.text(0.98, 0.95, 'Consensus: MAALLLLLLLLLLLLLLLLLLLLL', 
                transform=ax2.transAxes, fontsize=10, verticalalignment='top',
                horizontalalignment='right', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
        
        # Plot 3: Difference heatmap (Hits - Screen)
        ax3 = axes[2]
        
        # Calculate difference (enrichment in hits)
        pwm_diff = pwm_hits - pwm_screen
        
        # Transpose for visualization
        pwm_diff_t = pwm_diff.T
        
        sns.heatmap(pwm_diff_t, cmap='RdBu_r', center=0, 
                   cbar_kws={'label': 'Frequency Difference (Hits - Screen)'}, 
                   ax=ax3, linewidths=0.3, linecolor='gray', vmin=-0.15, vmax=0.15)
        ax3.set_xlabel('Position', fontsize=12)
        ax3.set_ylabel('Amino Acid', fontsize=12)
        ax3.set_title('Enrichment in Hits vs Full Screen (Red = Enriched in Hits, Blue = Depleted)', 
                     fontsize=14, fontweight='bold')
        
        # Plot 4: Key enriched positions bar plot
        ax4 = axes[3]
        
        # Calculate per-position enrichment score (sum of absolute differences)
        position_enrichment = pwm_diff.abs().sum(axis=1)
        
        positions = position_enrichment.index
        bars = ax4.bar(positions, position_enrichment.values, color='steelblue', 
                      edgecolor='black', alpha=0.8)
        
        # Highlight top positions
        top_positions = position_enrichment.nlargest(5)
        for pos in top_positions.index:
            idx = list(positions).index(pos)
            bars[idx].set_color('red')
            bars[idx].set_alpha(0.9)
        
        ax4.set_xlabel('Position', fontsize=12)
        ax4.set_ylabel('Total Enrichment Score', fontsize=12)
        ax4.set_title('Position-wise Enrichment Scores (Red = Top 5 Discriminating Positions)', 
                     fontsize=14, fontweight='bold')
        ax4.grid(axis='y', alpha=0.3)
        ax4.set_xticks(positions)
        
        # Add text labels for top positions
        for pos in top_positions.index:
            idx = list(positions).index(pos)
            height = position_enrichment.loc[pos]
            ax4.text(pos, height + 0.01, f'Pos {pos}', ha='center', va='bottom', 
                    fontsize=8, fontweight='bold', color='darkred')
        
        plt.tight_layout()
        
        # Save
        plt.savefig('logo_results/comparison_hits_vs_screen.png', dpi=300, bbox_inches='tight')
        print("Saved comparison plot to logo_results/comparison_hits_vs_screen.png")
        
        plt.close()
    
    # Create additional comparison - amino acid enrichment
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Calculate overall amino acid frequencies
    hits_overall = pwm_hits.mean()
    screen_overall = pwm_screen.mean()
    
    # Calculate enrichment ratios
    enrichment_ratios = hits_overall / screen_overall
    log2_enrichment = np.log2(enrichment_ratios)
    
    # Sort by enrichment
    log2_enrichment_sorted = log2_enrichment.sort_values(ascending=True)
    
    # Plot 1: Log2 enrichment
    ax1 = axes[0]
    colors = ['red' if x > 0 else 'blue' for x in log2_enrichment_sorted]
    ax1.barh(range(len(log2_enrichment_sorted)), log2_enrichment_sorted.values, 
            color=colors, alpha=0.7, edgecolor='black')
    ax1.set_yticks(range(len(log2_enrichment_sorted)))
    ax1.set_yticklabels(log2_enrichment_sorted.index, fontsize=11, family='monospace')
    ax1.axvline(x=0, color='black', linestyle='--', linewidth=2)
    ax1.set_xlabel('Log2 Enrichment (Hits/Screen)', fontsize=12)
    ax1.set_ylabel('Amino Acid', fontsize=12)
    ax1.set_title('Overall AA Enrichment in Hits vs Full Screen', fontsize=14, fontweight='bold')
    ax1.grid(axis='x', alpha=0.3)
    
    # Add value labels
    for i, (aa, val) in enumerate(log2_enrichment_sorted.items()):
        x_pos = val + (0.02 if val > 0 else -0.02)
        ha = 'left' if val > 0 else 'right'
        ax1.text(x_pos, i, f'{val:.2f}', va='center', ha=ha, fontsize=8)
    
    # Plot 2: Frequency comparison scatter
    ax2 = axes[1]
    
    ax2.scatter(screen_overall, hits_overall, s=200, alpha=0.6, c='purple', edgecolors='black')
    
    # Add AA labels
    for aa in screen_overall.index:
        ax2.annotate(aa, (screen_overall[aa], hits_overall[aa]), 
                    fontsize=10, ha='center', va='center', fontweight='bold')
    
    # Add diagonal line
    max_val = max(screen_overall.max(), hits_overall.max())
    ax2.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, linewidth=2)
    
    ax2.set_xlabel('Frequency in Full Screen', fontsize=12)
    ax2.set_ylabel('Frequency in Hits', fontsize=12)
    ax2.set_title('Amino Acid Frequency Comparison', fontsize=14, fontweight='bold')
    ax2.grid(alpha=0.3)
    ax2.set_xlim([0, max_val * 1.1])
    ax2.set_ylim([0, max_val * 1.1])
    
    # Add annotation
    ax2.text(0.05, 0.95, 'Above line = Enriched in Hits\nBelow line = Depleted in Hits', 
            transform=ax2.transAxes, fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.3))
    
    plt.tight_layout()
    plt.savefig('logo_results/aa_enrichment_comparison.png', dpi=300, bbox_inches='tight')
    print("Saved AA enrichment comparison to logo_results/aa_enrichment_comparison.png")
    plt.close()
    
    # Save enrichment data
    enrichment_df = pd.DataFrame({
        'amino_acid': hits_overall.index,
        'hits_frequency': hits_overall.values,
        'screen_frequency': screen_overall.values,
        'enrichment_ratio': enrichment_ratios.values,
        'log2_enrichment': log2_enrichment.values
    }).sort_values('log2_enrichment', ascending=False)
    
    enrichment_df.to_csv('logo_results/overall_aa_enrichment_comparison.csv', index=False)
    print("Saved enrichment data to logo_results/overall_aa_enrichment_comparison.csv")
    
    print("\n=== TOP ENRICHED AMINO ACIDS ===")
    print(enrichment_df.head(10).to_string(index=False))
    
    print("\n=== MOST DEPLETED AMINO ACIDS ===")
    print(enrichment_df.tail(10).to_string(index=False))


if __name__ == "__main__":
    create_comparison_plot()
    print("\nComparison analysis complete!")
