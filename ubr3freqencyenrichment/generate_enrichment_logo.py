#!/usr/bin/env python3
"""
Generate logo plots from enrichment data (hits vs screen)
Starting from position 2 to exclude the universal methionine
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

try:
    import logomaker
    LOGOMAKER_AVAILABLE = True
except ImportError:
    LOGOMAKER_AVAILABLE = False
    print("Warning: logomaker not installed")


def create_enrichment_logos():
    """Create logo plots based on enrichment ratios"""
    
    print("\n" + "="*80)
    print("CREATING ENRICHMENT-BASED LOGO PLOTS (Starting from Position 2)")
    print("="*80 + "\n")
    
    # Load PWMs from both datasets
    pwm_hits = pd.read_csv('logo_results_hits/position_weight_matrix.csv', index_col=0)
    pwm_screen = pd.read_csv('logo_results_full_screen/position_weight_matrix.csv', index_col=0)
    
    print(f"Loaded Hits PWM: {pwm_hits.shape}")
    print(f"Loaded Screen PWM: {pwm_screen.shape}")
    
    # Start from position 2 (skip position 1 which is always M)
    pwm_hits = pwm_hits.loc[2:]
    pwm_screen = pwm_screen.loc[2:]
    
    print(f"\nAnalyzing positions 2-{pwm_hits.index.max()}")
    
    # Calculate enrichment ratios with pseudocount
    pseudocount = 0.001
    enrichment_matrix = (pwm_hits + pseudocount) / (pwm_screen + pseudocount)
    
    # Calculate log2 enrichment
    log2_enrichment = np.log2(enrichment_matrix)
    
    # For logo plot, we want to show enrichment weighted by frequency
    # Higher frequency + higher enrichment = taller letter
    enrichment_weighted = pwm_hits * log2_enrichment
    
    # Clip extreme values for better visualization
    enrichment_weighted_clipped = enrichment_weighted.clip(-2, 2)
    
    # Save enrichment matrix
    os.makedirs('logo_results_enrichment', exist_ok=True)
    enrichment_matrix.to_csv('logo_results_enrichment/enrichment_matrix_pos2+.csv')
    log2_enrichment.to_csv('logo_results_enrichment/log2_enrichment_pos2+.csv')
    enrichment_weighted.to_csv('logo_results_enrichment/enrichment_weighted_pos2+.csv')
    
    print("\nTop enriched amino acids at each position (positions 2-10):")
    for pos in range(2, min(11, pwm_hits.index.max() + 1)):
        if pos in log2_enrichment.index:
            top_enriched = log2_enrichment.loc[pos].nlargest(3)
            print(f"\nPosition {pos}:")
            for aa, log2_enr in top_enriched.items():
                freq_hits = pwm_hits.loc[pos, aa]
                freq_screen = pwm_screen.loc[pos, aa]
                print(f"  {aa}: {freq_hits:.3f} (hits) vs {freq_screen:.3f} (screen) = {2**log2_enr:.2f}x (log2: {log2_enr:.2f})")
    
    if LOGOMAKER_AVAILABLE:
        # Create comprehensive figure
        fig = plt.figure(figsize=(20, 16))
        gs = fig.add_gridspec(5, 1, height_ratios=[1, 1, 1, 1, 0.8], hspace=0.35)
        
        # Plot 1: Hits frequency logo (for reference)
        ax1 = fig.add_subplot(gs[0])
        ic_hits = calculate_information_content(pwm_hits)
        logo1 = logomaker.Logo(ic_hits, ax=ax1, color_scheme='chemistry')
        ax1.set_ylabel('Information\nContent (bits)', fontsize=11)
        ax1.set_title('UBR3 Hits - Sequence Logo (Positions 2-24)', 
                     fontsize=14, fontweight='bold', pad=10)
        ax1.set_xlim([1.5, len(ic_hits) + 1.5])
        
        # Add consensus sequence
        consensus_hits = get_consensus(pwm_hits)
        ax1.text(0.98, 0.95, f'Consensus: M-{consensus_hits}', 
                transform=ax1.transAxes, fontsize=10, verticalalignment='top',
                horizontalalignment='right', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        # Plot 2: Screen frequency logo (for reference)
        ax2 = fig.add_subplot(gs[1])
        ic_screen = calculate_information_content(pwm_screen)
        logo2 = logomaker.Logo(ic_screen, ax=ax2, color_scheme='chemistry')
        ax2.set_ylabel('Information\nContent (bits)', fontsize=11)
        ax2.set_title('Full Screen - Sequence Logo (Positions 2-24)', 
                     fontsize=14, fontweight='bold', pad=10)
        ax2.set_xlim([1.5, len(ic_screen) + 1.5])
        
        consensus_screen = get_consensus(pwm_screen)
        ax2.text(0.98, 0.95, f'Consensus: M-{consensus_screen}', 
                transform=ax2.transAxes, fontsize=10, verticalalignment='top',
                horizontalalignment='right', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
        
        # Plot 3: Enrichment-weighted logo (KEY PLOT!)
        ax3 = fig.add_subplot(gs[2])
        
        # For positive enrichment, use hits frequency scaled by log2 enrichment
        # This shows: what's enriched AND abundant in hits
        enrichment_logo_matrix = pwm_hits.copy()
        for pos in enrichment_logo_matrix.index:
            for aa in enrichment_logo_matrix.columns:
                log2_enr = log2_enrichment.loc[pos, aa]
                if log2_enr > 0:
                    # Scale up enriched AAs
                    enrichment_logo_matrix.loc[pos, aa] *= (1 + log2_enr)
                else:
                    # Scale down depleted AAs
                    enrichment_logo_matrix.loc[pos, aa] *= max(0.1, 2**log2_enr)
        
        # Normalize each position to sum to reasonable height
        for pos in enrichment_logo_matrix.index:
            total = enrichment_logo_matrix.loc[pos].sum()
            if total > 0:
                enrichment_logo_matrix.loc[pos] = enrichment_logo_matrix.loc[pos] / total * 2
        
        logo3 = logomaker.Logo(enrichment_logo_matrix, ax=ax3, color_scheme='chemistry')
        ax3.set_ylabel('Enrichment-Weighted\nHeight', fontsize=11)
        ax3.set_title('Enrichment Logo - Shows Discriminating Features (Hits vs Screen)', 
                     fontsize=14, fontweight='bold', pad=10, color='darkred')
        ax3.set_xlim([1.5, len(enrichment_logo_matrix) + 1.5])
        
        # Highlight that this is the key plot
        ax3.text(0.98, 0.95, '*** KEY PLOT ***\nHeight = Frequency x Enrichment', 
                transform=ax3.transAxes, fontsize=10, verticalalignment='top',
                horizontalalignment='right', bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.9),
                fontweight='bold')
        
        # Plot 4: Log2 enrichment heatmap
        ax4 = fig.add_subplot(gs[3])
        
        # Transpose for better visualization
        log2_enr_t = log2_enrichment.T
        
        # Only show positions with interesting enrichment
        max_pos = min(24, log2_enr_t.shape[1])
        log2_enr_t_subset = log2_enr_t.iloc[:, :max_pos-1]
        
        sns.heatmap(log2_enr_t_subset, cmap='RdBu_r', center=0, 
                   cbar_kws={'label': 'Log2 Enrichment (Hits/Screen)'}, 
                   ax=ax4, linewidths=0.5, linecolor='white', 
                   vmin=-1.5, vmax=1.5, annot=False)
        ax4.set_xlabel('Position', fontsize=12)
        ax4.set_ylabel('Amino Acid', fontsize=12)
        ax4.set_title('Log2 Enrichment Heatmap (Red = Enriched in Hits, Blue = Depleted)', 
                     fontsize=14, fontweight='bold')
        
        # Plot 5: Position enrichment score
        ax5 = fig.add_subplot(gs[4])
        
        # Calculate position-wise enrichment score (variance of log2 enrichment)
        position_scores = []
        positions = []
        for pos in log2_enrichment.index:
            # Use weighted variance (more weight to frequent AAs)
            weights = pwm_hits.loc[pos].values
            log2_vals = log2_enrichment.loc[pos].values
            
            # Score = sum of (frequency * |log2_enrichment|)
            score = np.sum(weights * np.abs(log2_vals))
            position_scores.append(score)
            positions.append(pos)
        
        bars = ax5.bar(positions, position_scores, color='steelblue', edgecolor='black', alpha=0.8)
        
        # Highlight top 5 positions
        top_indices = np.argsort(position_scores)[-5:]
        for idx in top_indices:
            bars[idx].set_color('red')
            bars[idx].set_alpha(0.9)
            # Add label
            ax5.text(positions[idx], position_scores[idx] + 0.005, 
                    f'Pos {positions[idx]}', ha='center', va='bottom', 
                    fontsize=8, fontweight='bold', color='darkred')
        
        ax5.set_xlabel('Position', fontsize=12)
        ax5.set_ylabel('Enrichment Score', fontsize=12)
        ax5.set_title('Position Discrimination Score (Red = Top 5 Most Discriminating)', 
                     fontsize=13, fontweight='bold')
        ax5.grid(axis='y', alpha=0.3)
        ax5.set_xticks(positions)
        
        plt.savefig('logo_results_enrichment/enrichment_logo_comprehensive.png', 
                   dpi=300, bbox_inches='tight')
        print(f"\n[OK] Saved comprehensive enrichment logo to: logo_results_enrichment/enrichment_logo_comprehensive.png")
        plt.close()
        
        # Create a simplified version with just the enrichment logo
        fig, ax = plt.subplots(1, 1, figsize=(18, 6))
        
        logo = logomaker.Logo(enrichment_logo_matrix, ax=ax, color_scheme='chemistry')
        ax.set_ylabel('Enrichment-Weighted Height', fontsize=13, fontweight='bold')
        ax.set_xlabel('Position', fontsize=13, fontweight='bold')
        ax.set_title('UBR3 Recognition Motif - Enrichment Logo (Positions 2-24)\nHeight shows discriminating features (high = enriched in hits)', 
                    fontsize=15, fontweight='bold', pad=15)
        ax.set_xlim([1.5, len(enrichment_logo_matrix) + 1.5])
        
        # Add consensus and interpretation
        ax.text(0.02, 0.98, f'Consensus: M-{consensus_hits}', 
               transform=ax.transAxes, fontsize=11, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
               fontweight='bold')
        
        ax.text(0.98, 0.98, 'Tall letters = Enriched in hits\nShort letters = Depleted/Not discriminating', 
               transform=ax.transAxes, fontsize=10, verticalalignment='top',
               horizontalalignment='right', 
               bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
        
        plt.tight_layout()
        plt.savefig('logo_results_enrichment/enrichment_logo_simplified.png', 
                   dpi=300, bbox_inches='tight')
        print(f"[OK] Saved simplified enrichment logo to: logo_results_enrichment/enrichment_logo_simplified.png")
        plt.close()
        
    # Create enrichment bar plots for key positions
    create_position_enrichment_plots(pwm_hits, pwm_screen, log2_enrichment)
    
    # Generate summary statistics
    generate_enrichment_summary(log2_enrichment, pwm_hits, pwm_screen)
    
    print("\n" + "="*80)
    print("ENRICHMENT LOGO ANALYSIS COMPLETE!")
    print("="*80)
    print(f"\nResults saved to: logo_results_enrichment/")
    print("\nKey files:")
    print("  * enrichment_logo_comprehensive.png - Full comparison (*** BEST)")
    print("  * enrichment_logo_simplified.png - Clean enrichment logo")
    print("  * position_enrichment_plots.png - Detailed position analysis")
    print("  * enrichment_summary.txt - Statistical summary")
    print("="*80 + "\n")


def calculate_information_content(pwm):
    """Calculate information content matrix"""
    max_entropy = np.log2(20)
    ic_matrix = pd.DataFrame(0.0, index=pwm.index, columns=pwm.columns)
    
    for pos in pwm.index:
        freqs = pwm.loc[pos]
        entropy = 0
        for aa in freqs.index:
            if freqs[aa] > 0:
                entropy -= freqs[aa] * np.log2(freqs[aa])
        
        ic = max_entropy - entropy
        ic_matrix.loc[pos] = pwm.loc[pos] * ic
    
    return ic_matrix


def get_consensus(pwm):
    """Get consensus sequence from PWM"""
    consensus = ""
    for pos in pwm.index:
        max_aa = pwm.loc[pos].idxmax()
        consensus += max_aa
    return consensus


def create_position_enrichment_plots(pwm_hits, pwm_screen, log2_enrichment):
    """Create detailed enrichment plots for key positions"""
    
    print("\nCreating position-specific enrichment plots...")
    
    # Find top 8 discriminating positions
    position_scores = {}
    for pos in log2_enrichment.index:
        weights = pwm_hits.loc[pos].values
        log2_vals = log2_enrichment.loc[pos].values
        score = np.sum(weights * np.abs(log2_vals))
        position_scores[pos] = score
    
    top_positions = sorted(position_scores.items(), key=lambda x: x[1], reverse=True)[:8]
    top_pos_list = [p[0] for p in top_positions]
    
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()
    
    for idx, pos in enumerate(top_pos_list):
        ax = axes[idx]
        
        # Get data for this position
        hits_freq = pwm_hits.loc[pos].sort_values(ascending=False).head(10)
        screen_freq = pwm_screen.loc[pos].loc[hits_freq.index]
        log2_enr = log2_enrichment.loc[pos].loc[hits_freq.index]
        
        # Create grouped bar plot
        x = np.arange(len(hits_freq))
        width = 0.35
        
        bars1 = ax.bar(x - width/2, hits_freq.values, width, label='Hits', 
                      color='red', alpha=0.7, edgecolor='black')
        bars2 = ax.bar(x + width/2, screen_freq.values, width, label='Screen', 
                      color='blue', alpha=0.7, edgecolor='black')
        
        # Add enrichment as text
        for i, (aa, log2_val) in enumerate(zip(hits_freq.index, log2_enr.values)):
            if log2_val > 0.5:
                ax.text(i, max(hits_freq.values[i], screen_freq.values[i]) + 0.01, 
                       f'{2**log2_val:.1f}x', ha='center', va='bottom', 
                       fontsize=8, color='darkred', fontweight='bold')
        
        ax.set_xlabel('Amino Acid', fontsize=10)
        ax.set_ylabel('Frequency', fontsize=10)
        ax.set_title(f'Position {pos} (Score: {position_scores[pos]:.3f})', 
                    fontsize=11, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(hits_freq.index)
        ax.legend(fontsize=8)
        ax.grid(axis='y', alpha=0.3)
    
    plt.suptitle('Top 8 Discriminating Positions - Frequency Comparison (Hits vs Screen)', 
                fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig('logo_results_enrichment/position_enrichment_plots.png', 
               dpi=300, bbox_inches='tight')
    print("[OK] Saved position enrichment plots")
    plt.close()


def generate_enrichment_summary(log2_enrichment, pwm_hits, pwm_screen):
    """Generate summary statistics"""
    
    print("\nGenerating enrichment summary...")
    
    summary = []
    summary.append("="*80)
    summary.append("UBR3 ENRICHMENT ANALYSIS SUMMARY (Positions 2-24)")
    summary.append("="*80)
    summary.append("")
    
    # Top enriched features at each position
    summary.append("TOP ENRICHED AMINO ACIDS BY POSITION:")
    summary.append("-" * 80)
    
    for pos in range(2, min(25, log2_enrichment.index.max() + 1)):
        if pos not in log2_enrichment.index:
            continue
        
        top_3 = log2_enrichment.loc[pos].nlargest(3)
        summary.append(f"\nPosition {pos}:")
        
        for aa, log2_enr in top_3.items():
            freq_hits = pwm_hits.loc[pos, aa]
            freq_screen = pwm_screen.loc[pos, aa]
            enrichment = 2**log2_enr
            
            if enrichment > 1.2:  # Only show meaningful enrichments
                summary.append(f"  {aa}: {freq_hits:.3f} (hits) vs {freq_screen:.3f} (screen)")
                summary.append(f"      Enrichment: {enrichment:.2f}x (log2: {log2_enr:.2f})")
    
    summary.append("\n" + "="*80)
    summary.append("MOST DISCRIMINATING POSITIONS:")
    summary.append("-" * 80)
    
    # Calculate position discrimination scores
    position_scores = {}
    for pos in log2_enrichment.index:
        weights = pwm_hits.loc[pos].values
        log2_vals = log2_enrichment.loc[pos].values
        score = np.sum(weights * np.abs(log2_vals))
        position_scores[pos] = score
    
    top_positions = sorted(position_scores.items(), key=lambda x: x[1], reverse=True)[:10]
    
    for rank, (pos, score) in enumerate(top_positions, 1):
        # Get top enriched AA at this position
        top_aa = log2_enrichment.loc[pos].idxmax()
        top_log2 = log2_enrichment.loc[pos, top_aa]
        top_freq_hits = pwm_hits.loc[pos, top_aa]
        
        summary.append(f"\n{rank}. Position {pos} (Discrimination Score: {score:.3f})")
        summary.append(f"   Top enriched: {top_aa} ({top_freq_hits:.1%} in hits, {2**top_log2:.2f}x enrichment)")
    
    summary.append("\n" + "="*80)
    summary.append("KEY FINDINGS:")
    summary.append("-" * 80)
    
    # Identify key patterns
    summary.append("\n1. STRUCTURAL FEATURES:")
    
    # Check for proline enrichment at position 2
    if 2 in log2_enrichment.index:
        p_enrichment = 2**log2_enrichment.loc[2, 'P']
        p_freq = pwm_hits.loc[2, 'P']
        summary.append(f"   - Proline at position 2: {p_freq:.1%} ({p_enrichment:.2f}x)")
    
    # Check for glycine-rich region
    gly_positions = []
    for pos in range(10, min(22, log2_enrichment.index.max() + 1)):
        if pos in log2_enrichment.index and 'G' in log2_enrichment.columns:
            g_enrichment = 2**log2_enrichment.loc[pos, 'G']
            g_freq = pwm_hits.loc[pos, 'G']
            if g_enrichment > 1.5:
                gly_positions.append((pos, g_freq, g_enrichment))
    
    if gly_positions:
        summary.append(f"\n2. GLYCINE-RICH FLEXIBLE REGION:")
        for pos, freq, enr in gly_positions:
            summary.append(f"   - Position {pos}: {freq:.1%} glycine ({enr:.2f}x enrichment)")
    
    # Check for charged residues
    summary.append(f"\n3. CHARGED RESIDUES:")
    for aa in ['R', 'K', 'D', 'E']:
        if aa not in log2_enrichment.columns:
            continue
        avg_log2 = log2_enrichment[aa].mean()
        avg_freq_hits = pwm_hits[aa].mean()
        avg_freq_screen = pwm_screen[aa].mean()
        avg_enrichment = 2**avg_log2
        
        if abs(avg_log2) > 0.1:
            summary.append(f"   - {aa}: {avg_freq_hits:.2%} (hits) vs {avg_freq_screen:.2%} (screen) = {avg_enrichment:.2f}x")
    
    summary.append("\n" + "="*80)
    
    # Save summary
    with open('logo_results_enrichment/enrichment_summary.txt', 'w') as f:
        f.write('\n'.join(summary))
    
    print("[OK] Saved enrichment summary")
    
    # Also print to console
    print("\n" + '\n'.join(summary[:50]))  # Print first 50 lines


if __name__ == "__main__":
    create_enrichment_logos()
