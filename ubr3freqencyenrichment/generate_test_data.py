#!/usr/bin/env python3
"""
Generate test data for UBR3 enrichment analysis
Creates sample nucleotide sequences that can be used to test the analyzer
"""

import random
import argparse

def generate_random_nt_sequence(length):
    """Generate a random nucleotide sequence"""
    nucleotides = ['A', 'T', 'G', 'C']
    return ''.join(random.choices(nucleotides, k=length))

def generate_biased_sequence(length, enriched_codons):
    """
    Generate a sequence with enriched specific codons
    
    Args:
        length: Target sequence length (will be adjusted to nearest codon boundary)
        enriched_codons: List of codons to enrich
    """
    # Adjust length to codon boundary
    length = (length // 3) * 3
    
    # Standard codons for all amino acids
    all_codons = [
        'GCT', 'GCC', 'GCA', 'GCG',  # Alanine
        'CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG',  # Arginine
        'AAT', 'AAC',  # Asparagine
        'GAT', 'GAC',  # Aspartate
        'TGT', 'TGC',  # Cysteine
        'CAA', 'CAG',  # Glutamine
        'GAA', 'GAG',  # Glutamate
        'GGT', 'GGC', 'GGA', 'GGG',  # Glycine
        'CAT', 'CAC',  # Histidine
        'ATT', 'ATC', 'ATA',  # Isoleucine
        'TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG',  # Leucine
        'AAA', 'AAG',  # Lysine
        'ATG',  # Methionine (start)
        'TTT', 'TTC',  # Phenylalanine
        'CCT', 'CCC', 'CCA', 'CCG',  # Proline
        'TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC',  # Serine
        'ACT', 'ACC', 'ACA', 'ACG',  # Threonine
        'TGG',  # Tryptophan
        'TAT', 'TAC',  # Tyrosine
        'GTT', 'GTC', 'GTA', 'GTG',  # Valine
    ]
    
    sequence = 'ATG'  # Start with start codon
    remaining_length = length - 3
    
    # 70% enriched codons, 30% random
    while remaining_length >= 3:
        if random.random() < 0.7 and enriched_codons:
            codon = random.choice(enriched_codons)
        else:
            codon = random.choice(all_codons)
        sequence += codon
        remaining_length -= 3
    
    return sequence

def generate_fasta_file(filename, num_sequences, seq_length, enriched_codons=None):
    """Generate a FASTA file with test sequences"""
    with open(filename, 'w') as f:
        for i in range(num_sequences):
            f.write(f">sequence_{i+1}\n")
            if enriched_codons:
                seq = generate_biased_sequence(seq_length, enriched_codons)
            else:
                seq = generate_random_nt_sequence(seq_length)
            
            # Write sequence in 60-character lines
            for j in range(0, len(seq), 60):
                f.write(seq[j:j+60] + '\n')
    
    print(f"Generated {filename} with {num_sequences} sequences of length ~{seq_length}")

def generate_csv_file(filename, num_sequences, seq_length, enriched_codons=None):
    """Generate a CSV file with test sequences"""
    with open(filename, 'w') as f:
        f.write("id,sequence,score\n")
        for i in range(num_sequences):
            if enriched_codons:
                seq = generate_biased_sequence(seq_length, enriched_codons)
            else:
                seq = generate_random_nt_sequence(seq_length)
            score = random.randint(50, 100)
            f.write(f"seq_{i+1},{seq},{score}\n")
    
    print(f"Generated {filename} with {num_sequences} sequences of length ~{seq_length}")

def main():
    parser = argparse.ArgumentParser(
        description='Generate test data for UBR3 enrichment analysis'
    )
    parser.add_argument(
        '-n', '--num-sequences',
        type=int,
        default=100,
        help='Number of sequences to generate (default: 100)'
    )
    parser.add_argument(
        '-l', '--length',
        type=int,
        default=150,
        help='Length of each sequence in nucleotides (default: 150)'
    )
    parser.add_argument(
        '-f', '--format',
        choices=['fasta', 'csv', 'both'],
        default='fasta',
        help='Output format (default: fasta)'
    )
    parser.add_argument(
        '-e', '--enrich',
        action='store_true',
        help='Generate sequences with enriched specific amino acids'
    )
    
    args = parser.parse_args()
    
    # Define enriched codons if requested
    enriched_codons = None
    if args.enrich:
        # Enrich for specific amino acids (example: Lysine, Arginine, Leucine)
        enriched_codons = [
            'AAA', 'AAG',  # Lysine
            'CGT', 'CGC', 'AGA', 'AGG',  # Arginine
            'TTA', 'TTG', 'CTT', 'CTC',  # Leucine
            'GCT', 'GCC',  # Alanine
        ]
        print("Generating sequences enriched for: Lysine (K), Arginine (R), Leucine (L), Alanine (A)")
    
    # Generate files
    if args.format in ['fasta', 'both']:
        generate_fasta_file(
            'test_ubr3_sequences.fasta',
            args.num_sequences,
            args.length,
            enriched_codons
        )
    
    if args.format in ['csv', 'both']:
        generate_csv_file(
            'test_ubr3_sequences.csv',
            args.num_sequences,
            args.length,
            enriched_codons
        )
    
    print("\nTest data generated successfully!")
    print("You can now run:")
    if args.format in ['fasta', 'both']:
        print("  python ubr3_enrichment_analysis.py test_ubr3_sequences.fasta")
    if args.format in ['csv', 'both']:
        print("  python ubr3_enrichment_analysis.py test_ubr3_sequences.csv")

if __name__ == "__main__":
    main()
