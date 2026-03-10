#!/usr/bin/env python3
"""
UBR3 Enrichment Analysis Script
Performs enrichment analysis on amino acid sequences from UBR3 NT screen
"""

import pandas as pd
import numpy as np
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from Bio.Seq import Seq
from Bio import SeqIO
import argparse
import os

class UBR3EnrichmentAnalyzer:
    """Analyzer for UBR3 amino acid sequence enrichment"""
    
    def __init__(self, input_file, output_dir="results"):
        """
        Initialize the analyzer
        
        Args:
            input_file: Path to the UBR3 best file (CSV, TSV, or FASTA)
            output_dir: Directory to save results
        """
        self.input_file = input_file
        self.output_dir = output_dir
        self.data = None
        self.aa_sequences = []
        self.enrichment_results = {}
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
    def load_data(self):
        """Load the UBR3 best file"""
        print(f"Loading data from {self.input_file}...")
        
        file_ext = os.path.splitext(self.input_file)[1].lower()
        
        if file_ext in ['.fasta', '.fa', '.fna']:
            # Load FASTA file
            sequences = []
            for record in SeqIO.parse(self.input_file, "fasta"):
                sequences.append({
                    'id': record.id,
                    'nt_sequence': str(record.seq),
                    'description': record.description
                })
            self.data = pd.DataFrame(sequences)
            
        elif file_ext in ['.xlsx', '.xls']:
            # Load Excel file
            self.data = pd.read_excel(self.input_file)
            
            # Check if this is a gene list file (no sequence columns)
            has_sequences = any(col.lower() in ['aa_seq', 'aa_sequence', 'amino_acid', 'protein', 'nt_seq', 'nt_sequence', 'nucleotide', 'sequence'] 
                              for col in self.data.columns)
            
            if not has_sequences and len(self.data.columns) <= 4:
                # This looks like a gene list file, try to load sequences from screen file
                print("Detected gene list format. Looking up sequences from screen file...")
                screen_file = 'UBR3 Nt screen.xlsx'
                if os.path.exists(screen_file):
                    screen_data = pd.read_excel(screen_file)
                    
                    # Try to find gene column in the gene list
                    gene_col = None
                    
                    # Heuristic: find column with short names that vary (not all the same)
                    for col in self.data.columns:
                        col_values = self.data[col].dropna().astype(str).tolist()
                        unique_values = set(col_values)
                        
                        # Skip if all values are the same
                        if len(unique_values) <= 1:
                            continue
                        
                        # Check if values look like gene names
                        sample_values = list(unique_values)[:5]
                        is_gene_col = all(
                            not val.startswith('ENST') and 
                            not val.endswith('.pdf') and 
                            not val.lower() == 'ubn24withatg' and
                            3 < len(val) < 20 
                            for val in sample_values
                        )
                        
                        if is_gene_col and len(unique_values) > 10:  # Should have multiple unique genes
                            gene_col = col
                            break
                    
                    if gene_col is not None:
                        gene_list = self.data[gene_col].dropna().unique()
                        print(f"Found {len(gene_list)} unique genes in hits list: {gene_list[:5].tolist()}...")
                        
                        # Filter screen data to only include these genes
                        if 'Gene_ID' in screen_data.columns:
                            self.data = screen_data[screen_data['Gene_ID'].isin(gene_list)].copy()
                            print(f"Matched {len(self.data)} sequences from screen file")
                            if len(self.data) == 0:
                                print("Warning: No matching genes found in screen file")
                                print(f"Sample genes from hits: {gene_list[:5].tolist()}")
                                print(f"Sample genes from screen: {screen_data['Gene_ID'].head(5).tolist()}")
                        else:
                            print("Warning: Could not find 'Gene_ID' column in screen file")
                            print(f"Available columns: {screen_data.columns.tolist()}")
                    else:
                        print("Warning: Could not identify gene name column")
                        print(f"Available columns: {self.data.columns.tolist()}")
                else:
                    print(f"Warning: Screen file '{screen_file}' not found. Cannot lookup sequences.")
            
        elif file_ext in ['.csv']:
            self.data = pd.read_csv(self.input_file)
            
        elif file_ext in ['.tsv', '.txt']:
            self.data = pd.read_csv(self.input_file, sep='\t')
            
        else:
            raise ValueError(f"Unsupported file format: {file_ext}")
        
        print(f"Loaded {len(self.data)} sequences")
        return self.data
    
    def translate_sequences(self):
        """Translate nucleotide sequences to amino acid sequences or use existing AA sequences"""
        
        # First check if AA sequences already exist
        aa_col = None
        for col in self.data.columns:
            if col.lower() in ['aa_seq', 'aa_sequence', 'amino_acid', 'protein']:
                aa_col = col
                break
        
        if aa_col is not None:
            print(f"Using existing AA sequences from column '{aa_col}'...")
            self.data['aa_sequence'] = self.data[aa_col]
            self.aa_sequences = [str(seq) for seq in self.data[aa_col] if pd.notna(seq) and str(seq).strip() != '']
            print(f"Loaded {len(self.aa_sequences)} AA sequences")
            return self.aa_sequences
        
        # If no AA sequences, translate from NT
        print("Translating NT sequences to AA sequences...")
        
        # Find the column containing nucleotide sequences
        nt_col = None
        for col in self.data.columns:
            if any(keyword in col.lower() for keyword in ['nt_seq', 'nt_sequence', 'nucleotide']):
                nt_col = col
                break
        
        if nt_col is None and 'nt_sequence' in self.data.columns:
            nt_col = 'nt_sequence'
        elif nt_col is None:
            # Try the first column that looks like a sequence
            for col in self.data.columns:
                if self.data[col].dtype == object:
                    sample = str(self.data[col].iloc[0])
                    if any(nt in sample.upper() for nt in ['A', 'T', 'G', 'C']):
                        nt_col = col
                        break
        
        if nt_col is None:
            raise ValueError("Could not find nucleotide or amino acid sequence column. Please check your input file.")
        
        # Translate sequences
        aa_sequences = []
        for nt_seq in self.data[nt_col]:
            try:
                seq = Seq(str(nt_seq))
                aa_seq = str(seq.translate())
                aa_sequences.append(aa_seq)
            except Exception as e:
                print(f"Warning: Could not translate sequence {str(nt_seq)[:30]}... : {e}")
                aa_sequences.append("")
        
        self.data['aa_sequence'] = aa_sequences
        self.aa_sequences = [seq for seq in aa_sequences if seq]
        
        print(f"Successfully translated {len(self.aa_sequences)} sequences")
        return aa_sequences
    
    def calculate_aa_frequency(self):
        """Calculate amino acid frequency at each position"""
        print("Calculating amino acid frequencies...")
        
        # Check if we have sequences
        if not self.aa_sequences:
            print("Error: No amino acid sequences found!")
            return pd.DataFrame()
        
        # Find the maximum sequence length
        max_length = max(len(seq) for seq in self.aa_sequences)
        
        # Count AA at each position
        position_counts = defaultdict(Counter)
        for seq in self.aa_sequences:
            for pos, aa in enumerate(seq):
                position_counts[pos][aa] += 1
        
        # Convert to DataFrame
        freq_data = []
        for pos in range(max_length):
            total = sum(position_counts[pos].values())
            if total > 0:
                for aa, count in position_counts[pos].items():
                    freq_data.append({
                        'position': pos + 1,
                        'amino_acid': aa,
                        'count': count,
                        'frequency': count / total
                    })
        
        self.position_frequency = pd.DataFrame(freq_data)
        return self.position_frequency
    
    def calculate_overall_aa_enrichment(self):
        """Calculate overall amino acid enrichment across all sequences"""
        print("Calculating overall AA enrichment...")
        
        # Check if we have sequences
        if not self.aa_sequences:
            print("Error: No amino acid sequences found!")
            return pd.DataFrame()
        
        # Expected frequencies (based on standard codon usage)
        expected_freq = {
            'A': 0.0825, 'C': 0.0137, 'D': 0.0545, 'E': 0.0675,
            'F': 0.0386, 'G': 0.0707, 'H': 0.0227, 'I': 0.0596,
            'K': 0.0584, 'L': 0.0966, 'M': 0.0242, 'N': 0.0406,
            'P': 0.0470, 'Q': 0.0393, 'R': 0.0553, 'S': 0.0656,
            'T': 0.0534, 'V': 0.0687, 'W': 0.0108, 'Y': 0.0292,
            '*': 0.0010  # Stop codon
        }
        
        # Count all amino acids
        all_aa = ''.join(self.aa_sequences)
        observed_counts = Counter(all_aa)
        total_aa = len(all_aa)
        
        # Calculate enrichment
        enrichment_data = []
        for aa in observed_counts.keys():
            observed_freq = observed_counts[aa] / total_aa
            expected = expected_freq.get(aa, 0.05)  # Default if AA not in expected
            
            # Calculate enrichment ratio
            if expected > 0:
                enrichment_ratio = observed_freq / expected
            else:
                enrichment_ratio = float('inf')
            
            # Chi-square test
            expected_count = total_aa * expected
            chi2_stat = ((observed_counts[aa] - expected_count) ** 2) / expected_count if expected_count > 0 else 0
            p_value = 1 - stats.chi2.cdf(chi2_stat, df=1)
            
            enrichment_data.append({
                'amino_acid': aa,
                'observed_count': observed_counts[aa],
                'observed_frequency': observed_freq,
                'expected_frequency': expected,
                'enrichment_ratio': enrichment_ratio,
                'log2_enrichment': np.log2(enrichment_ratio) if enrichment_ratio > 0 and enrichment_ratio != float('inf') else 0,
                'chi2_statistic': chi2_stat,
                'p_value': p_value
            })
        
        self.overall_enrichment = pd.DataFrame(enrichment_data)
        self.overall_enrichment = self.overall_enrichment.sort_values('enrichment_ratio', ascending=False)
        
        return self.overall_enrichment
    
    def analyze_motifs(self, window_size=3):
        """Analyze enriched motifs/patterns in sequences"""
        print(f"Analyzing {window_size}-mer motifs...")
        
        # Check if we have sequences
        if not self.aa_sequences:
            print("Error: No amino acid sequences found!")
            return pd.DataFrame()
        
        motif_counts = Counter()
        for seq in self.aa_sequences:
            for i in range(len(seq) - window_size + 1):
                motif = seq[i:i+window_size]
                motif_counts[motif] += 1
        
        # Get top motifs
        total_motifs = sum(motif_counts.values())
        motif_data = []
        for motif, count in motif_counts.most_common(50):
            motif_data.append({
                'motif': motif,
                'count': count,
                'frequency': count / total_motifs
            })
        
        self.motif_enrichment = pd.DataFrame(motif_data)
        return self.motif_enrichment
    
    def plot_aa_enrichment(self):
        """Create visualization of AA enrichment"""
        print("Creating enrichment plots...")
        
        # Check if we have data to plot
        if not hasattr(self, 'overall_enrichment') or self.overall_enrichment.empty:
            print("Warning: No enrichment data to plot. Skipping plot generation.")
            return
        
        # Figure 1: Overall AA enrichment bar plot
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # Plot 1: Enrichment ratio
        ax1 = axes[0, 0]
        data = self.overall_enrichment.sort_values('enrichment_ratio', ascending=True)
        colors = ['red' if x > 1 else 'blue' for x in data['enrichment_ratio']]
        ax1.barh(data['amino_acid'], data['enrichment_ratio'], color=colors, alpha=0.7)
        ax1.axvline(x=1, color='black', linestyle='--', linewidth=1)
        ax1.set_xlabel('Enrichment Ratio')
        ax1.set_ylabel('Amino Acid')
        ax1.set_title('UBR3 AA Enrichment Ratios')
        ax1.grid(axis='x', alpha=0.3)
        
        # Plot 2: Log2 enrichment
        ax2 = axes[0, 1]
        data_log = self.overall_enrichment.sort_values('log2_enrichment', ascending=True)
        colors_log = ['red' if x > 0 else 'blue' for x in data_log['log2_enrichment']]
        ax2.barh(data_log['amino_acid'], data_log['log2_enrichment'], color=colors_log, alpha=0.7)
        ax2.axvline(x=0, color='black', linestyle='--', linewidth=1)
        ax2.set_xlabel('Log2 Enrichment')
        ax2.set_ylabel('Amino Acid')
        ax2.set_title('UBR3 AA Log2 Enrichment')
        ax2.grid(axis='x', alpha=0.3)
        
        # Plot 3: Observed vs Expected frequency
        ax3 = axes[1, 0]
        ax3.scatter(self.overall_enrichment['expected_frequency'], 
                   self.overall_enrichment['observed_frequency'],
                   s=100, alpha=0.6, c='purple')
        for idx, row in self.overall_enrichment.iterrows():
            ax3.annotate(row['amino_acid'], 
                        (row['expected_frequency'], row['observed_frequency']),
                        fontsize=8)
        max_val = max(self.overall_enrichment['expected_frequency'].max(),
                     self.overall_enrichment['observed_frequency'].max())
        ax3.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)
        ax3.set_xlabel('Expected Frequency')
        ax3.set_ylabel('Observed Frequency')
        ax3.set_title('Observed vs Expected AA Frequencies')
        ax3.grid(alpha=0.3)
        
        # Plot 4: Top enriched motifs
        ax4 = axes[1, 1]
        if hasattr(self, 'motif_enrichment') and not self.motif_enrichment.empty:
            top_motifs = self.motif_enrichment.head(20)
            ax4.barh(range(len(top_motifs)), top_motifs['count'], color='green', alpha=0.7)
            ax4.set_yticks(range(len(top_motifs)))
            ax4.set_yticklabels(top_motifs['motif'])
            ax4.set_xlabel('Count')
            ax4.set_ylabel('Motif')
            ax4.set_title('Top 20 Enriched 3-mer Motifs')
            ax4.grid(axis='x', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'ubr3_enrichment_analysis.png'), dpi=300)
        print(f"Saved plot to {os.path.join(self.output_dir, 'ubr3_enrichment_analysis.png')}")
        
        # Figure 2: Heatmap of position-specific frequencies (if available)
        if hasattr(self, 'position_frequency') and not self.position_frequency.empty:
            # Create pivot table for heatmap
            pivot_data = self.position_frequency.pivot_table(
                values='frequency',
                index='amino_acid',
                columns='position',
                fill_value=0
            )
            
            # Limit to first 50 positions if too long
            if pivot_data.shape[1] > 50:
                pivot_data = pivot_data.iloc[:, :50]
            
            plt.figure(figsize=(20, 8))
            sns.heatmap(pivot_data, cmap='YlOrRd', cbar_kws={'label': 'Frequency'})
            plt.title('Position-Specific AA Frequency in UBR3 Sequences')
            plt.xlabel('Position')
            plt.ylabel('Amino Acid')
            plt.tight_layout()
            plt.savefig(os.path.join(self.output_dir, 'ubr3_position_heatmap.png'), dpi=300)
            print(f"Saved heatmap to {os.path.join(self.output_dir, 'ubr3_position_heatmap.png')}")
        
        plt.close('all')
    
    def save_results(self):
        """Save all results to CSV files"""
        print("Saving results...")
        
        # Save overall enrichment
        if hasattr(self, 'overall_enrichment') and not self.overall_enrichment.empty:
            self.overall_enrichment.to_csv(
                os.path.join(self.output_dir, 'ubr3_overall_enrichment.csv'),
                index=False
            )
        else:
            print("Warning: No overall enrichment data to save")
        
        # Save position frequencies
        if hasattr(self, 'position_frequency'):
            self.position_frequency.to_csv(
                os.path.join(self.output_dir, 'ubr3_position_frequencies.csv'),
                index=False
            )
        
        # Save motif enrichment
        if hasattr(self, 'motif_enrichment'):
            self.motif_enrichment.to_csv(
                os.path.join(self.output_dir, 'ubr3_motif_enrichment.csv'),
                index=False
            )
        
        # Save translated sequences
        self.data.to_csv(
            os.path.join(self.output_dir, 'ubr3_translated_sequences.csv'),
            index=False
        )
        
        print(f"All results saved to {self.output_dir}/")
    
    def generate_summary_report(self):
        """Generate a summary report (console output only)"""
        print("\n" + "="*80)
        print("UBR3 ENRICHMENT ANALYSIS SUMMARY REPORT")
        print("="*80)
        print(f"\nInput file: {self.input_file}")
        print(f"Total sequences analyzed: {len(self.aa_sequences)}")
        print(f"Total amino acids: {sum(len(seq) for seq in self.aa_sequences)}")
        
        print("\n--- TOP 10 MOST ENRICHED AMINO ACIDS ---")
        print(self.overall_enrichment[['amino_acid', 'enrichment_ratio', 'log2_enrichment', 'p_value']].head(10).to_string(index=False))
        
        print("\n--- TOP 10 MOST DEPLETED AMINO ACIDS ---")
        print(self.overall_enrichment[['amino_acid', 'enrichment_ratio', 'log2_enrichment', 'p_value']].tail(10).to_string(index=False))
        
        if hasattr(self, 'motif_enrichment'):
            print("\n--- TOP 10 MOST COMMON MOTIFS ---")
            print(self.motif_enrichment[['motif', 'count', 'frequency']].head(10).to_string(index=False))
        
        print("\n" + "="*80)
    
    def run_full_analysis(self):
        """Run the complete enrichment analysis pipeline"""
        print("\n" + "="*80)
        print("STARTING UBR3 ENRICHMENT ANALYSIS")
        print("="*80 + "\n")
        
        # Load data
        self.load_data()
        
        # Translate sequences
        self.translate_sequences()
        
        # Calculate frequencies
        self.calculate_aa_frequency()
        
        # Calculate enrichment
        self.calculate_overall_aa_enrichment()
        
        # Analyze motifs
        self.analyze_motifs(window_size=3)
        
        # Create plots
        self.plot_aa_enrichment()
        
        # Save results
        self.save_results()
        
        # Generate report
        self.generate_summary_report()
        
        print("\nAnalysis complete!")


class UBR3ComparativeAnalyzer:
    """Comparative analyzer for UBR3 hits vs background screen"""
    
    def __init__(self, hits_file, screen_file, output_dir="comparative_results", start_position=2):
        """
        Initialize the comparative analyzer
        
        Args:
            hits_file: Path to the UBR3 best hits file
            screen_file: Path to the UBR3 NT screen file
            output_dir: Directory to save results
            start_position: Starting position for analysis (1-based, default=2)
        """
        self.hits_file = hits_file
        self.screen_file = screen_file
        self.output_dir = output_dir
        self.start_position = start_position
        self.start_idx = start_position - 1  # Convert to 0-based index
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        self.hits_data = None
        self.screen_data = None
        self.hits_aa_sequences = []
        self.screen_aa_sequences = []
        
    def load_and_translate(self, file_path, label):
        """Load data and translate sequences"""
        print(f"\n{'='*60}")
        print(f"Loading {label} from {file_path}...")
        print(f"{'='*60}")
        
        file_ext = os.path.splitext(file_path)[1].lower()
        
        if file_ext in ['.xlsx', '.xls']:
            data = pd.read_excel(file_path)
        elif file_ext in ['.csv']:
            data = pd.read_csv(file_path)
        elif file_ext in ['.tsv', '.txt']:
            data = pd.read_csv(file_path, sep='\t')
        else:
            raise ValueError(f"Unsupported file format: {file_ext}")
        
        print(f"Loaded {len(data)} sequences")
        
        # Check if this is a gene list file (no sequence columns)
        has_sequences = any(col.lower() in ['aa_seq', 'aa_sequence', 'amino_acid', 'protein', 'nt_seq', 'nt_sequence', 'nucleotide', 'sequence'] 
                          for col in data.columns)
        
        if not has_sequences and len(data.columns) <= 4 and label == "HITS":
            # This looks like a gene list file, try to load sequences from screen file
            print("Detected gene list format. Looking up sequences from screen file...")
            screen_file = 'UBR3 Nt screen.xlsx'
            if os.path.exists(screen_file):
                screen_data = pd.read_excel(screen_file)
                
                # Find gene column with varying values
                gene_col = None
                for col in data.columns:
                    col_values = data[col].dropna().astype(str).tolist()
                    unique_values = set(col_values)
                    
                    if len(unique_values) <= 1:
                        continue
                    
                    sample_values = list(unique_values)[:5]
                    is_gene_col = all(
                        not val.startswith('ENST') and 
                        not val.endswith('.pdf') and 
                        not val.lower() == 'ubn24withatg' and
                        3 < len(val) < 20 
                        for val in sample_values
                    )
                    
                    if is_gene_col and len(unique_values) > 10:
                        gene_col = col
                        break
                
                if gene_col is not None:
                    gene_list = data[gene_col].dropna().unique()
                    print(f"Found {len(gene_list)} unique genes in hits list")
                    
                    if 'Gene_ID' in screen_data.columns:
                        data = screen_data[screen_data['Gene_ID'].isin(gene_list)].copy()
                        print(f"Matched {len(data)} sequences from screen file")
        
        print(f"Final data size: {len(data)} sequences")
        
        # Check for existing AA sequences
        aa_col = None
        for col in data.columns:
            if col.lower() in ['aa_seq', 'aa_sequence', 'amino_acid', 'protein']:
                aa_col = col
                break
        
        if aa_col is not None:
            print(f"Using existing AA sequences from column '{aa_col}'...")
            aa_sequences = [str(seq) for seq in data[aa_col] if pd.notna(seq) and str(seq).strip() != '']
        else:
            # Translate from NT
            print("Translating NT sequences to AA sequences...")
            nt_col = None
            for col in data.columns:
                if any(keyword in col.lower() for keyword in ['nt_seq', 'nt_sequence', 'nucleotide', 'sequence']):
                    nt_col = col
                    break
            
            if nt_col is None:
                for col in data.columns:
                    if data[col].dtype == object:
                        sample = str(data[col].iloc[0])
                        if any(nt in sample.upper() for nt in ['A', 'T', 'G', 'C']):
                            nt_col = col
                            break
            
            if nt_col is None:
                raise ValueError(f"Could not find sequence column in {label}")
            
            aa_sequences = []
            for nt_seq in data[nt_col]:
                try:
                    seq = Seq(str(nt_seq))
                    aa_seq = str(seq.translate())
                    if len(aa_seq) >= self.start_position:  # Only keep sequences long enough
                        aa_sequences.append(aa_seq)
                except Exception as e:
                    pass
        
        # Filter sequences that are long enough to start from start_position
        aa_sequences = [seq for seq in aa_sequences if len(seq) >= self.start_position]
        
        print(f"Successfully loaded {len(aa_sequences)} AA sequences (filtered for length >= {self.start_position})")
        return data, aa_sequences
    
    def calculate_position_frequencies(self, sequences, label):
        """Calculate AA frequency at each position starting from start_position"""
        print(f"\nCalculating position-specific frequencies for {label}...")
        
        # Check if we have sequences
        if not sequences:
            print(f"Error: No sequences found for {label}!")
            return pd.DataFrame()
        
        # Find max length
        max_length = max(len(seq) for seq in sequences)
        
        # Count AA at each position (starting from start_idx)
        position_counts = defaultdict(Counter)
        for seq in sequences:
            for pos in range(self.start_idx, len(seq)):
                aa = seq[pos]
                position_counts[pos][aa] += 1
        
        # Convert to frequency matrix
        freq_data = []
        for pos in range(self.start_idx, max_length):
            total = sum(position_counts[pos].values())
            if total > 0:
                for aa in 'ACDEFGHIKLMNPQRSTVWY*':
                    count = position_counts[pos].get(aa, 0)
                    freq_data.append({
                        'position': pos + 1,  # 1-based
                        'amino_acid': aa,
                        'count': count,
                        'frequency': count / total if total > 0 else 0
                    })
        
        return pd.DataFrame(freq_data)
    
    def calculate_enrichment(self):
        """Calculate enrichment of hits vs screen"""
        print("\n" + "="*60)
        print("CALCULATING ENRICHMENT (HITS vs SCREEN)")
        print("="*60)
        
        # Calculate frequencies for both
        self.hits_freq = self.calculate_position_frequencies(self.hits_aa_sequences, "HITS")
        self.screen_freq = self.calculate_position_frequencies(self.screen_aa_sequences, "SCREEN")
        
        # Merge and calculate enrichment
        merged = pd.merge(
            self.hits_freq,
            self.screen_freq,
            on=['position', 'amino_acid'],
            how='outer',
            suffixes=('_hits', '_screen')
        ).fillna(0)
        
        # Calculate enrichment ratio with pseudocount to avoid division by zero
        pseudocount = 0.0001
        merged['enrichment_ratio'] = (merged['frequency_hits'] + pseudocount) / (merged['frequency_screen'] + pseudocount)
        merged['log2_enrichment'] = np.log2(merged['enrichment_ratio'])
        
        # Cap extreme values for visualization
        merged['log2_enrichment_capped'] = merged['log2_enrichment'].clip(-5, 5)
        
        self.enrichment_data = merged
        
        # Save results
        self.enrichment_data.to_csv(
            os.path.join(self.output_dir, 'position_enrichment_data.csv'),
            index=False
        )
        print(f"Saved enrichment data to {os.path.join(self.output_dir, 'position_enrichment_data.csv')}")
        
        return self.enrichment_data
    
    def calculate_overall_enrichment(self):
        """Calculate overall AA enrichment (not position-specific)"""
        print("\nCalculating overall AA enrichment...")
        
        # Get sequences starting from start_idx
        hits_trimmed = [''.join(seq[self.start_idx:]) for seq in self.hits_aa_sequences]
        screen_trimmed = [''.join(seq[self.start_idx:]) for seq in self.screen_aa_sequences]
        
        # Count all amino acids
        hits_counts = Counter(''.join(hits_trimmed))
        screen_counts = Counter(''.join(screen_trimmed))
        
        hits_total = sum(hits_counts.values())
        screen_total = sum(screen_counts.values())
        
        # Calculate enrichment for each AA
        enrichment_data = []
        all_aa = set(hits_counts.keys()) | set(screen_counts.keys())
        
        for aa in sorted(all_aa):
            hits_count = hits_counts.get(aa, 0)
            screen_count = screen_counts.get(aa, 0)
            
            hits_freq = hits_count / hits_total if hits_total > 0 else 0
            screen_freq = screen_count / screen_total if screen_total > 0 else 0
            
            # Calculate enrichment with pseudocount
            pseudocount = 0.0001
            enrichment_ratio = (hits_freq + pseudocount) / (screen_freq + pseudocount)
            log2_enrichment = np.log2(enrichment_ratio)
            
            # Fisher's exact test for significance
            contingency = [
                [hits_count, hits_total - hits_count],
                [screen_count, screen_total - screen_count]
            ]
            try:
                oddsratio, p_value = stats.fisher_exact(contingency)
            except:
                p_value = 1.0
                oddsratio = enrichment_ratio
            
            enrichment_data.append({
                'amino_acid': aa,
                'hits_count': hits_count,
                'screen_count': screen_count,
                'hits_frequency': hits_freq,
                'screen_frequency': screen_freq,
                'enrichment_ratio': enrichment_ratio,
                'log2_enrichment': log2_enrichment,
                'p_value': p_value,
                'odds_ratio': oddsratio
            })
        
        self.overall_enrichment = pd.DataFrame(enrichment_data)
        self.overall_enrichment = self.overall_enrichment.sort_values('enrichment_ratio', ascending=False)
        
        # Save results
        self.overall_enrichment.to_csv(
            os.path.join(self.output_dir, 'overall_aa_enrichment.csv'),
            index=False
        )
        
        return self.overall_enrichment
    
    def create_heatmaps(self):
        """Create comprehensive heatmaps for the enrichment analysis"""
        print("\nCreating heatmaps...")
        
        # Prepare data for heatmap - using enrichment ratio instead of log2
        pivot_enrichment = self.enrichment_data.pivot_table(
            values='enrichment_ratio',
            index='amino_acid',
            columns='position',
            fill_value=1  # Fill with 1 (no enrichment) instead of 0
        )
        
        pivot_hits_freq = self.hits_freq.pivot_table(
            values='frequency',
            index='amino_acid',
            columns='position',
            fill_value=0
        )
        
        pivot_screen_freq = self.screen_freq.pivot_table(
            values='frequency',
            index='amino_acid',
            columns='position',
            fill_value=0
        )
        
        # Limit positions if too many
        max_positions = 50
        if pivot_enrichment.shape[1] > max_positions:
            pivot_enrichment = pivot_enrichment.iloc[:, :max_positions]
            pivot_hits_freq = pivot_hits_freq.iloc[:, :max_positions]
            pivot_screen_freq = pivot_screen_freq.iloc[:, :max_positions]
        
        # Create custom colormap: beige/light orange for <1, red shades for >1
        from matplotlib.colors import LinearSegmentedColormap
        colors_list = [
            '#FFEAD3',  # Custom beige (depletion - low values ~0)
            '#FFEAD3',  # Light beige
            '#FFE8CC',  # Soft beige
            '#FFE5C5',  # Medium beige
            '#FFE2BE',  # Warm beige
            '#FFDEB3',  # Wheat beige (approaching neutral)
            '#FFD7A8',  # Light orange-beige (neutral - around 1.0)
            '#FFB6A3',  # Light coral-red (enrichment starts)
            '#FF9A8B',  # Soft salmon-red
            '#FF8075',  # Medium coral-red
            '#FF6B6B',  # Warm red
            '#E85D5D',  # Medium red
            '#D14545',  # Strong red (high enrichment ~2.5)
        ]
        custom_cmap = LinearSegmentedColormap.from_list('beige_red_diverging', colors_list, N=256)
        
        # Create figure with only the enrichment heatmap
        fig, ax1 = plt.subplots(1, 1, figsize=(20, 10))
        
        # Heatmap: Enrichment ratio (starting from 0)
        sns.heatmap(
            pivot_enrichment,
            cmap=custom_cmap,
            center=1,  # Center at 1 (no enrichment/depletion)
            vmin=0,    # Start from 0
            vmax=2.5,  # Max enrichment ratio to show
            cbar_kws={'label': 'Enrichment Ratio (Hits/Screen)'},
            ax=ax1,
            linewidths=0.3,
            linecolor='white'
        )
        ax1.set_title(f'UBR3 AA Enrichment Heatmap (Hits vs Screen) - Starting from Position {self.start_position}', 
                     fontsize=16, fontweight='bold', pad=20)
        ax1.set_xlabel('Position', fontsize=14)
        ax1.set_ylabel('Amino Acid', fontsize=14)
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'enrichment_heatmap_warm.png'), dpi=300, bbox_inches='tight')
        print(f"Saved enrichment heatmap to {os.path.join(self.output_dir, 'enrichment_heatmap_warm.png')}")
        plt.close()
        
        # Create a second figure for overall enrichment
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        
        # Plot 1: Enrichment ratio bar plot
        ax1 = axes[0]
        data = self.overall_enrichment.sort_values('enrichment_ratio', ascending=True)
        colors = ['red' if x > 1 else 'blue' for x in data['enrichment_ratio']]
        ax1.barh(data['amino_acid'], data['enrichment_ratio'], color=colors, alpha=0.7)
        ax1.axvline(x=1, color='black', linestyle='--', linewidth=2)
        ax1.set_xlabel('Enrichment Ratio (Hits/Screen)', fontsize=12)
        ax1.set_ylabel('Amino Acid', fontsize=12)
        ax1.set_title(f'Overall AA Enrichment (Position {self.start_position}+)', fontsize=14, fontweight='bold')
        ax1.set_xlim(left=0)  # Start scale from zero
        ax1.grid(axis='x', alpha=0.3)
        
        # Plot 2: Log2 enrichment
        ax2 = axes[1]
        data_log = self.overall_enrichment.sort_values('log2_enrichment', ascending=True)
        colors_log = ['red' if x > 0 else 'blue' for x in data_log['log2_enrichment']]
        ax2.barh(data_log['amino_acid'], data_log['log2_enrichment'], color=colors_log, alpha=0.7)
        ax2.axvline(x=0, color='black', linestyle='--', linewidth=2)
        ax2.set_xlabel('Log2 Enrichment', fontsize=12)
        ax2.set_ylabel('Amino Acid', fontsize=12)
        ax2.set_title(f'Overall AA Log2 Enrichment (Position {self.start_position}+)', fontsize=14, fontweight='bold')
        # Log2 enrichment scale keeps centered at 0 (can be negative or positive)
        ax2.grid(axis='x', alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'overall_enrichment_barplots.png'), dpi=300, bbox_inches='tight')
        print(f"Saved bar plots to {os.path.join(self.output_dir, 'overall_enrichment_barplots.png')}")
        plt.close()
    
    def generate_summary_report(self):
        """Generate a summary report"""
        print("\n" + "="*80)
        print("UBR3 COMPARATIVE ENRICHMENT ANALYSIS SUMMARY")
        print("="*80)
        print(f"\nHits file: {self.hits_file}")
        print(f"Screen file: {self.screen_file}")
        print(f"Analysis starting from position: {self.start_position}")
        print(f"\nTotal HITS sequences: {len(self.hits_aa_sequences)}")
        print(f"Total SCREEN sequences: {len(self.screen_aa_sequences)}")
        
        print("\n--- TOP 10 MOST ENRICHED AMINO ACIDS (in Hits vs Screen) ---")
        print(self.overall_enrichment[['amino_acid', 'enrichment_ratio', 'log2_enrichment', 
                                       'hits_frequency', 'screen_frequency', 'p_value']].head(10).to_string(index=False))
        
        print("\n--- TOP 10 MOST DEPLETED AMINO ACIDS (in Hits vs Screen) ---")
        print(self.overall_enrichment[['amino_acid', 'enrichment_ratio', 'log2_enrichment',
                                       'hits_frequency', 'screen_frequency', 'p_value']].tail(10).to_string(index=False))
        
        print("\n" + "="*80)
    
    def run_comparative_analysis(self):
        """Run the complete comparative analysis pipeline"""
        print("\n" + "="*80)
        print("STARTING UBR3 COMPARATIVE ENRICHMENT ANALYSIS")
        print("="*80)
        
        # Load both datasets
        self.hits_data, self.hits_aa_sequences = self.load_and_translate(self.hits_file, "HITS")
        self.screen_data, self.screen_aa_sequences = self.load_and_translate(self.screen_file, "SCREEN")
        
        # Calculate enrichment
        self.calculate_enrichment()
        
        # Calculate overall enrichment
        self.calculate_overall_enrichment()
        
        # Create heatmaps
        self.create_heatmaps()
        
        # Generate report
        self.generate_summary_report()
        
        print("\n" + "="*80)
        print("COMPARATIVE ANALYSIS COMPLETE!")
        print("="*80)
        print(f"Results saved to: {self.output_dir}/")


def main():
    """Main function to run the script"""
    parser = argparse.ArgumentParser(
        description='UBR3 Enrichment Analysis - Analyze amino acid enrichment in UBR3 NT screen'
    )
    parser.add_argument(
        'input_file',
        nargs='?',  # Make it optional
        default='ubr3_best (1).xlsx',
        help='Path to the UBR3 best file (FASTA, CSV, or TSV format) - default: ubr3_best (1).xlsx'
    )
    parser.add_argument(
        '-o', '--output',
        default='results',
        help='Output directory for results (default: results)'
    )
    parser.add_argument(
        '-w', '--window',
        type=int,
        default=3,
        help='Window size for motif analysis (default: 3)'
    )
    parser.add_argument(
        '-c', '--compare',
        action='store_true',
        help='Run comparative analysis (hits vs screen)'
    )
    parser.add_argument(
        '-s', '--screen',
        default='UBR3 Nt screen.xlsx',
        help='Path to the screen file for comparative analysis (default: UBR3 Nt screen.xlsx)'
    )
    parser.add_argument(
        '-p', '--start-position',
        type=int,
        default=2,
        help='Starting position for comparative analysis (1-based, default: 2)'
    )
    
    args = parser.parse_args()
    
    # Check if comparative analysis is requested
    if args.compare:
        # Run comparative analysis
        print("\n" + "="*80)
        print("COMPARATIVE MODE: Analyzing Hits vs Screen")
        print("="*80)
        
        # Check if files exist
        if not os.path.exists(args.input_file):
            print(f"Error: Hits file '{args.input_file}' not found!")
            return
        if not os.path.exists(args.screen):
            print(f"Error: Screen file '{args.screen}' not found!")
            return
        
        # Run comparative analysis
        comparative_output = args.output if args.output != 'results' else 'comparative_results'
        analyzer = UBR3ComparativeAnalyzer(
            hits_file=args.input_file,
            screen_file=args.screen,
            output_dir=comparative_output,
            start_position=args.start_position
        )
        analyzer.run_comparative_analysis()
    else:
        # Run standard single-file analysis
        if not os.path.exists(args.input_file):
            print(f"Error: Input file '{args.input_file}' not found!")
            return
        
        analyzer = UBR3EnrichmentAnalyzer(args.input_file, args.output)
        analyzer.run_full_analysis()


if __name__ == "__main__":
    main()
