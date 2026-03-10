# UBR3 Enrichment Analysis

A Python script for performing enrichment analysis on amino acid sequences from UBR3 nucleotide screen data.

## Features

- **Nucleotide to Amino Acid Translation**: Automatically translates NT sequences to AA sequences
- **Frequency Analysis**: Calculates amino acid frequencies at each position
- **Enrichment Calculation**: Computes enrichment ratios comparing observed vs expected frequencies
- **Statistical Testing**: Chi-square tests for significance
- **Motif Discovery**: Identifies enriched sequence motifs (3-mers by default)
- **Visualization**: Creates comprehensive plots including:
  - Enrichment ratio bar charts
  - Log2 enrichment plots
  - Observed vs expected frequency scatter plots
  - Position-specific frequency heatmaps
  - Top motif bar charts
- **Detailed Reports**: Generates CSV files and summary reports

## Installation

1. Install Python 3.8 or higher
2. Install required packages:

```bash
pip install -r requirements.txt
```

## Usage

### Basic Usage

```bash
python ubr3_enrichment_analysis.py <input_file>
```

### With Options

```bash
python ubr3_enrichment_analysis.py <input_file> -o output_folder -w 4
```

### Arguments

- `input_file`: Path to your UBR3 best file (required)
  - Supported formats: FASTA (.fasta, .fa), CSV (.csv), TSV (.tsv, .txt)
  - Should contain nucleotide sequences that will be translated to amino acids

- `-o, --output`: Output directory for results (default: `results`)
- `-w, --window`: Window size for motif analysis (default: 3)

## Input File Format

The script accepts multiple file formats:

### FASTA Format
```
>sequence1
ATGGCTAGCTAGC...
>sequence2
ATGCGATCGATCG...
```

### CSV/TSV Format
Must contain a column with nucleotide sequences. The script will auto-detect columns named:
- `sequence`, `seq`, `nt`, `nucleotide`, `nt_sequence`, etc.

Example CSV:
```
id,sequence,score
seq1,ATGGCTAGC...,100
seq2,ATGCGATCG...,95
```

## Output Files

The script creates several output files in the specified output directory:

1. **ubr3_overall_enrichment.csv** - Overall AA enrichment statistics
   - Enrichment ratios, log2 enrichment, p-values for each amino acid

2. **ubr3_position_frequencies.csv** - Position-specific AA frequencies
   - Frequency of each AA at each position in sequences

3. **ubr3_motif_enrichment.csv** - Top enriched motifs
   - Most common sequence patterns

4. **ubr3_translated_sequences.csv** - Original data with AA translations

5. **ubr3_enrichment_analysis.png** - Comprehensive visualization (4 subplots)

6. **ubr3_position_heatmap.png** - Position-specific frequency heatmap

## Example

```bash
# Run analysis on your UBR3 best file
python ubr3_enrichment_analysis.py ubr3_best_sequences.fasta

# Results will be saved in ./results/ directory
```

## Understanding the Results

### Enrichment Ratio
- **> 1**: Amino acid is enriched (more frequent than expected)
- **< 1**: Amino acid is depleted (less frequent than expected)
- **= 1**: Amino acid appears at expected frequency

### Log2 Enrichment
- **Positive values**: Enriched
- **Negative values**: Depleted
- **Zero**: No change from expected

### P-values
- **< 0.05**: Statistically significant enrichment/depletion
- **≥ 0.05**: Not significant

## Requirements

- Python 3.8+
- pandas
- numpy
- matplotlib
- seaborn
- scipy
- biopython

## Notes

- Expected amino acid frequencies are based on standard codon usage tables
- Chi-square tests are used for statistical significance
- Position-specific analysis shows which positions in your sequences have enriched AAs
- Motif analysis helps identify recurring patterns in your sequences
