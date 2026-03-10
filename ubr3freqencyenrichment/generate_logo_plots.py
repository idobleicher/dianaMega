#!/usr/bin/env python3
"""
Generate sequence logo plots from UBR3 screen data
Shows consensus sequence and information content
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import os
import argparse

try:
    import logomaker
    LOGOMAKER_AVAILABLE = True
except ImportError:
    LOGOMAKER_AVAILABLE = False
    print("Warning: logomaker not installed. Install with: pip install logomaker")


class SequenceLogoGenerator:
    """Generate sequence logo plots from amino acid sequences"""
    
    def __init__(self, input_file, output_dir="logo_results", start_position=1, max_positions=24):
        """
        Initialize the logo generator
        
        Args:
            input_file: Path to the translated sequences CSV or Excel file
            output_dir: Directory to save results
            start_position: Starting position for logo (1-based)
            max_positions: Maximum number of positions to show
        """
        self.input_file = input_file
        self.output_dir = output_dir
        self.start_position = start_position
        self.start_idx = start_position - 1  # Convert to 0-based
        self.max_positions = max_positions
        self.sequences = []
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
    
    def load_sequences(self):
        """Load amino acid sequences from file"""
        print(f"Loading sequences from {self.input_file}...")
        
        file_ext = os.path.splitext(self.input_file)[1].lower()
        
        if file_ext in ['.xlsx', '.xls']:
            data = pd.read_excel(self.input_file)
        elif file_ext in ['.csv']:
            data = pd.read_csv(self.input_file)
        else:
            raise ValueError(f"Unsupported file format: {file_ext}")
        
        # Find AA sequence column
        aa_col = None
        for col in data.columns:
            if 'aa_seq' in col.lower() or col.lower() == 'aa_sequence':
                aa_col = col
                break
        
        if aa_col is None:
            raise ValueError("Could not find amino acid sequence column in file")
        
        # Load sequences
        self.sequences = [str(seq) for seq in data[aa_col] if pd.notna(seq) and str(seq).strip() != '']
        
        print(f"Loaded {len(self.sequences)} sequences")
        
        # Print some example sequences
        print("\nExample sequences:")
        for i, seq in enumerate(self.sequences[:5]):
            print(f"  {i+1}. {seq}")
        
        return self.sequences
    
    def create_position_weight_matrix(self):
        """Create position weight matrix from sequences"""
        print("\nCreating position weight matrix...")
        
        # Determine sequence length range
        if not self.sequences:
            raise ValueError("No sequences loaded")
        
        # Find the effective max length
        max_len = max(len(seq) for seq in self.sequences)
        end_position = min(self.start_idx + self.max_positions, max_len)
        
        # Count amino acids at each position
        amino_acids = list('ACDEFGHIKLMNPQRSTVWY*')
        position_counts = {}
        
        for pos in range(self.start_idx, end_position):
            counts = Counter()
            for seq in self.sequences:
                if len(seq) > pos:
                    aa = seq[pos]
                    if aa in amino_acids:
                        counts[aa] += 1
            
            # Convert to frequencies
            total = sum(counts.values())
            if total > 0:
                position_counts[pos - self.start_idx + 1] = {aa: counts.get(aa, 0) / total for aa in amino_acids}
        
        # Convert to DataFrame
        pwm_df = pd.DataFrame(position_counts).T
        pwm_df = pwm_df.fillna(0)
        
        # Remove amino acids that never appear (for cleaner logo)
        pwm_df = pwm_df.loc[:, (pwm_df != 0).any(axis=0)]
        
        print(f"Created PWM with {len(pwm_df)} positions and {len(pwm_df.columns)} amino acids")
        
        self.pwm = pwm_df
        return pwm_df
    
    def calculate_information_content(self):
        """Calculate information content matrix for logo plot"""
        print("Calculating information content...")
        
        # Calculate information content using Shannon entropy
        # IC = log2(20) - entropy
        # For each position: IC_pos = sum(p * log2(p / background))
        
        background_freq = 1.0 / 20  # Uniform background
        max_entropy = np.log2(20)  # Maximum possible entropy for 20 amino acids
        
        ic_matrix = pd.DataFrame(0.0, index=self.pwm.index, columns=self.pwm.columns)
        
        for pos in self.pwm.index:
            # Calculate entropy for this position
            freqs = self.pwm.loc[pos]
            entropy = 0
            for aa in freqs.index:
                if freqs[aa] > 0:
                    entropy -= freqs[aa] * np.log2(freqs[aa])
            
            # Information content = max_entropy - entropy
            ic = max_entropy - entropy
            
            # Scale frequencies by information content
            ic_matrix.loc[pos] = self.pwm.loc[pos] * ic
        
        self.ic_matrix = ic_matrix
        return ic_matrix
    
    def get_consensus_sequence(self):
        """Extract consensus sequence from PWM"""
        print("Extracting consensus sequence...")
        
        consensus = ""
        consensus_scores = []
        
        for pos in self.pwm.index:
            # Get the most frequent amino acid at this position
            max_aa = self.pwm.loc[pos].idxmax()
            max_freq = self.pwm.loc[pos].max()
            consensus += max_aa
            consensus_scores.append(max_freq)
        
        print(f"\nConsensus sequence (positions {self.start_position}-{self.start_position + len(consensus) - 1}):")
        print(f"  {consensus}")
        print(f"  Average conservation: {np.mean(consensus_scores):.3f}")
        
        self.consensus = consensus
        self.consensus_scores = consensus_scores
        
        return consensus
    
    def create_logo_plot_logomaker(self):
        """Create sequence logo using logomaker library"""
        if not LOGOMAKER_AVAILABLE:
            print("Logomaker not available. Skipping logo plot.")
            return
        
        print("\nCreating sequence logo plot...")
        
        # Create figure with multiple subplots
        fig, axes = plt.subplots(3, 1, figsize=(16, 12))
        
        # Plot 1: Information content logo
        ax1 = axes[0]
        logo1 = logomaker.Logo(self.ic_matrix, ax=ax1, color_scheme='chemistry')
        ax1.set_ylabel('Information Content\n(bits)', fontsize=12)
        ax1.set_title(f'UBR3 Screen Consensus Sequence (Positions {self.start_position}-{self.start_position + len(self.ic_matrix) - 1})', 
                     fontsize=14, fontweight='bold', pad=20)
        ax1.set_xlim([0.5, len(self.ic_matrix) + 0.5])
        
        # Plot 2: Probability logo
        ax2 = axes[1]
        logo2 = logomaker.Logo(self.pwm, ax=ax2, color_scheme='hydrophobicity')
        ax2.set_ylabel('Probability', fontsize=12)
        ax2.set_title('Position-Specific Amino Acid Probabilities', fontsize=14, fontweight='bold')
        ax2.set_xlim([0.5, len(self.pwm) + 0.5])
        
        # Plot 3: Consensus bar plot
        ax3 = axes[2]
        positions = list(range(self.start_position, self.start_position + len(self.consensus)))
        colors = []
        for i, aa in enumerate(self.consensus):
            # Color by conservation score
            score = self.consensus_scores[i]
            if score > 0.5:
                colors.append('darkgreen')
            elif score > 0.3:
                colors.append('orange')
            else:
                colors.append('gray')
        
        bars = ax3.bar(positions, self.consensus_scores, color=colors, alpha=0.7, edgecolor='black')
        ax3.set_ylabel('Conservation Score', fontsize=12)
        ax3.set_xlabel('Position', fontsize=12)
        ax3.set_title('Consensus Sequence Conservation', fontsize=14, fontweight='bold')
        ax3.set_xticks(positions)
        ax3.set_ylim([0, 1.0])
        ax3.grid(axis='y', alpha=0.3)
        
        # Add amino acid labels on bars
        for i, (pos, aa, score) in enumerate(zip(positions, self.consensus, self.consensus_scores)):
            ax3.text(pos, score + 0.02, aa, ha='center', va='bottom', 
                    fontsize=10, fontweight='bold')
        
        plt.tight_layout()
        
        # Save plot
        output_path = os.path.join(self.output_dir, 'ubr3_sequence_logo.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved sequence logo to {output_path}")
        
        plt.close()
    
    def create_alternative_logo(self):
        """Create alternative visualization if logomaker is not available"""
        print("\nCreating alternative sequence visualization...")
        
        fig, axes = plt.subplots(2, 1, figsize=(18, 10))
        
        # Plot 1: Heatmap of position frequencies
        ax1 = axes[0]
        
        # Transpose for better visualization (amino acids as rows, positions as columns)
        pwm_t = self.pwm.T
        
        sns.heatmap(pwm_t, cmap='YlOrRd', annot=False, cbar_kws={'label': 'Frequency'}, 
                   ax=ax1, linewidths=0.5, linecolor='white')
        ax1.set_xlabel('Position', fontsize=12)
        ax1.set_ylabel('Amino Acid', fontsize=12)
        ax1.set_title(f'Position-Specific Amino Acid Frequencies (Positions {self.start_position}-{self.start_position + len(self.pwm) - 1})', 
                     fontsize=14, fontweight='bold')
        
        # Adjust x-axis labels to show actual positions
        positions = list(range(self.start_position, self.start_position + len(self.pwm)))
        ax1.set_xticklabels(positions, rotation=0)
        
        # Plot 2: Stacked bar chart showing amino acid composition at each position
        ax2 = axes[1]
        
        # Get top amino acids at each position for stacked bars
        bottom = np.zeros(len(self.pwm))
        
        # Get unique amino acids that appear
        all_aas = sorted(self.pwm.columns)
        
        # Use a color palette
        colors_dict = {
            'A': '#8DD3C7', 'C': '#FFFFB3', 'D': '#BEBADA', 'E': '#FB8072',
            'F': '#80B1D3', 'G': '#FDB462', 'H': '#B3DE69', 'I': '#FCCDE5',
            'K': '#D9D9D9', 'L': '#BC80BD', 'M': '#CCEBC5', 'N': '#FFED6F',
            'P': '#E5C494', 'Q': '#B3B3B3', 'R': '#8DD3C7', 'S': '#FFFFB3',
            'T': '#BEBADA', 'V': '#FB8072', 'W': '#80B1D3', 'Y': '#FDB462',
            '*': '#000000'
        }
        
        for aa in all_aas:
            if aa in self.pwm.columns:
                values = self.pwm[aa].values
                ax2.bar(positions, values, bottom=bottom, label=aa, 
                       color=colors_dict.get(aa, '#CCCCCC'), edgecolor='white', linewidth=0.5)
                bottom += values
        
        ax2.set_xlabel('Position', fontsize=12)
        ax2.set_ylabel('Frequency', fontsize=12)
        ax2.set_title('Amino Acid Composition by Position (Stacked)', fontsize=14, fontweight='bold')
        ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', ncol=2)
        ax2.set_xticks(positions)
        ax2.set_ylim([0, 1.0])
        ax2.grid(axis='y', alpha=0.3)
        
        plt.tight_layout()
        
        # Save plot
        output_path = os.path.join(self.output_dir, 'ubr3_sequence_heatmap.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved sequence heatmap to {output_path}")
        
        plt.close()
    
    def analyze_position_specific_enrichment(self):
        """Analyze which amino acids are enriched at each position"""
        print("\nAnalyzing position-specific enrichment...")
        
        # Expected background frequencies for amino acids
        background = {aa: 0.05 for aa in self.pwm.columns}  # Uniform for simplicity
        
        enrichment_data = []
        
        for pos in self.pwm.index:
            for aa in self.pwm.columns:
                observed = self.pwm.loc[pos, aa]
                expected = background.get(aa, 0.05)
                
                if expected > 0:
                    enrichment = observed / expected
                    log2_enrichment = np.log2(enrichment) if enrichment > 0 else 0
                else:
                    enrichment = 0
                    log2_enrichment = 0
                
                enrichment_data.append({
                    'position': pos + self.start_idx,
                    'amino_acid': aa,
                    'observed_freq': observed,
                    'expected_freq': expected,
                    'enrichment_ratio': enrichment,
                    'log2_enrichment': log2_enrichment
                })
        
        enrichment_df = pd.DataFrame(enrichment_data)
        
        # Save to CSV
        output_path = os.path.join(self.output_dir, 'position_specific_enrichment.csv')
        enrichment_df.to_csv(output_path, index=False)
        print(f"Saved position-specific enrichment to {output_path}")
        
        # Find top enriched amino acids at each position
        print("\nTop enriched amino acids by position:")
        for pos in sorted(enrichment_df['position'].unique())[:10]:  # Show first 10 positions
            pos_data = enrichment_df[enrichment_df['position'] == pos].nlargest(3, 'enrichment_ratio')
            print(f"\nPosition {pos}:")
            for _, row in pos_data.iterrows():
                print(f"  {row['amino_acid']}: {row['observed_freq']:.3f} (enrichment: {row['enrichment_ratio']:.2f}x)")
        
        return enrichment_df
    
    def generate_motif_patterns(self, window_size=3, top_n=20):
        """Extract and visualize common sequence motifs"""
        print(f"\nAnalyzing {window_size}-mer motifs...")
        
        motif_counts = Counter()
        for seq in self.sequences:
            # Only analyze from start position
            seq_region = seq[self.start_idx:self.start_idx + self.max_positions]
            for i in range(len(seq_region) - window_size + 1):
                motif = seq_region[i:i+window_size]
                if len(motif) == window_size:
                    motif_counts[motif] += 1
        
        # Get top motifs
        top_motifs = motif_counts.most_common(top_n)
        
        print(f"\nTop {top_n} enriched {window_size}-mer motifs:")
        for i, (motif, count) in enumerate(top_motifs, 1):
            freq = count / len(self.sequences)
            print(f"  {i}. {motif}: {count} occurrences ({freq:.3f})")
        
        # Create visualization
        fig, ax = plt.subplots(figsize=(12, 8))
        
        motifs = [m[0] for m in top_motifs]
        counts = [m[1] for m in top_motifs]
        
        bars = ax.barh(range(len(motifs)), counts, color='steelblue', alpha=0.8, edgecolor='black')
        ax.set_yticks(range(len(motifs)))
        ax.set_yticklabels(motifs, fontsize=11, family='monospace')
        ax.set_xlabel('Count', fontsize=12)
        ax.set_ylabel(f'{window_size}-mer Motif', fontsize=12)
        ax.set_title(f'Top {top_n} Enriched {window_size}-mer Motifs in UBR3 Screen', 
                    fontsize=14, fontweight='bold')
        ax.grid(axis='x', alpha=0.3)
        ax.invert_yaxis()
        
        # Add count labels
        for i, (bar, count) in enumerate(zip(bars, counts)):
            ax.text(count + max(counts)*0.01, i, str(count), 
                   va='center', fontsize=9)
        
        plt.tight_layout()
        
        output_path = os.path.join(self.output_dir, f'ubr3_top_motifs_{window_size}mer.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Saved motif plot to {output_path}")
        
        plt.close()
        
        return top_motifs
    
    def save_results(self):
        """Save all results to files"""
        print("\nSaving results...")
        
        # Save PWM
        pwm_path = os.path.join(self.output_dir, 'position_weight_matrix.csv')
        self.pwm.to_csv(pwm_path)
        print(f"Saved PWM to {pwm_path}")
        
        # Save information content matrix
        if hasattr(self, 'ic_matrix'):
            ic_path = os.path.join(self.output_dir, 'information_content_matrix.csv')
            self.ic_matrix.to_csv(ic_path)
            print(f"Saved IC matrix to {ic_path}")
        
        # Save consensus sequence
        consensus_path = os.path.join(self.output_dir, 'consensus_sequence.txt')
        with open(consensus_path, 'w') as f:
            f.write(f"Consensus sequence (positions {self.start_position}-{self.start_position + len(self.consensus) - 1}):\n")
            f.write(f"{self.consensus}\n\n")
            f.write("Position-by-position breakdown:\n")
            for i, (aa, score) in enumerate(zip(self.consensus, self.consensus_scores)):
                pos = self.start_position + i
                f.write(f"Position {pos}: {aa} ({score:.3f})\n")
        print(f"Saved consensus sequence to {consensus_path}")
    
    def run_full_analysis(self):
        """Run complete logo plot generation and analysis"""
        print("\n" + "="*80)
        print("STARTING UBR3 SEQUENCE LOGO ANALYSIS")
        print("="*80 + "\n")
        
        # Load sequences
        self.load_sequences()
        
        # Create PWM
        self.create_position_weight_matrix()
        
        # Calculate information content
        self.calculate_information_content()
        
        # Get consensus
        self.get_consensus_sequence()
        
        # Create logo plots
        if LOGOMAKER_AVAILABLE:
            self.create_logo_plot_logomaker()
        else:
            print("\nWarning: logomaker not available. Creating alternative visualizations.")
        
        # Always create alternative visualizations
        self.create_alternative_logo()
        
        # Analyze enrichment
        self.analyze_position_specific_enrichment()
        
        # Analyze motifs
        for window_size in [3, 4, 5]:
            self.generate_motif_patterns(window_size=window_size, top_n=20)
        
        # Save results
        self.save_results()
        
        print("\n" + "="*80)
        print("SEQUENCE LOGO ANALYSIS COMPLETE!")
        print("="*80)
        print(f"\nResults saved to: {self.output_dir}/")
        print(f"Consensus sequence: {self.consensus}")
        print("="*80 + "\n")


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description='Generate sequence logo plots from UBR3 screen data'
    )
    parser.add_argument(
        'input_file',
        nargs='?',
        default='results/ubr3_translated_sequences.csv',
        help='Path to translated sequences file (default: results/ubr3_translated_sequences.csv)'
    )
    parser.add_argument(
        '-o', '--output',
        default='logo_results',
        help='Output directory (default: logo_results)'
    )
    parser.add_argument(
        '-s', '--start-position',
        type=int,
        default=1,
        help='Starting position for logo (1-based, default: 1)'
    )
    parser.add_argument(
        '-m', '--max-positions',
        type=int,
        default=24,
        help='Maximum number of positions to show (default: 24)'
    )
    
    args = parser.parse_args()
    
    # Check if file exists
    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' not found!")
        return
    
    # Run analysis
    generator = SequenceLogoGenerator(
        input_file=args.input_file,
        output_dir=args.output,
        start_position=args.start_position,
        max_positions=args.max_positions
    )
    generator.run_full_analysis()


if __name__ == "__main__":
    main()
