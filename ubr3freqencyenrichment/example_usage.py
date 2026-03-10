#!/usr/bin/env python3
"""
Example usage of the UBR3 Enrichment Analyzer
This script demonstrates how to use the analyzer programmatically
"""

from ubr3_enrichment_analysis import UBR3EnrichmentAnalyzer

# Example 1: Run full analysis with default parameters
def example_full_analysis():
    """Run complete analysis pipeline"""
    analyzer = UBR3EnrichmentAnalyzer(
        input_file='ubr3_best_sequences.fasta',  # Replace with your file
        output_dir='results'
    )
    analyzer.run_full_analysis()


# Example 2: Step-by-step analysis with custom parameters
def example_custom_analysis():
    """Run analysis step by step with custom options"""
    
    # Initialize analyzer
    analyzer = UBR3EnrichmentAnalyzer(
        input_file='ubr3_best_sequences.csv',  # Replace with your file
        output_dir='custom_results'
    )
    
    # Load and translate
    analyzer.load_data()
    analyzer.translate_sequences()
    
    # Calculate frequencies and enrichment
    analyzer.calculate_aa_frequency()
    analyzer.calculate_overall_aa_enrichment()
    
    # Analyze motifs with custom window size
    analyzer.analyze_motifs(window_size=4)  # Look for 4-mers instead of 3-mers
    
    # Generate outputs
    analyzer.plot_aa_enrichment()
    analyzer.save_results()
    analyzer.generate_summary_report()
    
    # Access results programmatically
    print("\nTop 5 enriched amino acids:")
    print(analyzer.overall_enrichment.head())
    
    print("\nTop 5 motifs:")
    print(analyzer.motif_enrichment.head())


# Example 3: Analyze specific aspects only
def example_targeted_analysis():
    """Focus on specific analysis aspects"""
    
    analyzer = UBR3EnrichmentAnalyzer('ubr3_best_sequences.fasta', 'targeted_results')
    
    # Load and translate
    analyzer.load_data()
    analyzer.translate_sequences()
    
    # Only calculate overall enrichment
    enrichment_df = analyzer.calculate_overall_aa_enrichment()
    
    # Filter for significantly enriched AAs (p < 0.05 and enrichment > 1.5)
    significant = enrichment_df[
        (enrichment_df['p_value'] < 0.05) & 
        (enrichment_df['enrichment_ratio'] > 1.5)
    ]
    
    print("Significantly enriched amino acids:")
    print(significant)


if __name__ == "__main__":
    # Run the example you want
    # Uncomment the example you'd like to run:
    
    # example_full_analysis()
    # example_custom_analysis()
    # example_targeted_analysis()
    
    print("Please uncomment one of the examples and add your input file path!")
