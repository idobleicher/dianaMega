# Quick Start Guide - UBR3 Enrichment Analysis

## Installation (2 minutes)

1. **Install dependencies:**
```bash
pip install -r requirements.txt
```

## Usage

### Option 1: Test with Sample Data (Recommended First)

Generate test data and run analysis:

```bash
# Generate 100 test sequences with enriched amino acids
python generate_test_data.py -n 100 -e

# Run analysis on test data
python ubr3_enrichment_analysis.py test_ubr3_sequences.fasta
```

Results will be in the `results/` folder.

### Option 2: Analyze Your Own Data

```bash
# If you have a FASTA file
python ubr3_enrichment_analysis.py your_ubr3_best_file.fasta

# If you have a CSV or TSV file
python ubr3_enrichment_analysis.py your_ubr3_best_file.csv

# Specify custom output directory
python ubr3_enrichment_analysis.py your_file.fasta -o my_results
```

## What You'll Get

After running the analysis, you'll find these files in your output directory:

1. **Visualizations:**
   - `ubr3_enrichment_analysis.png` - 4-panel plot showing enrichment ratios, log2 enrichment, observed vs expected, and top motifs
   - `ubr3_position_heatmap.png` - Heatmap showing which amino acids are enriched at specific positions

2. **Data Files (CSV):**
   - `ubr3_overall_enrichment.csv` - Enrichment statistics for each amino acid
   - `ubr3_position_frequencies.csv` - Position-by-position AA frequencies
   - `ubr3_motif_enrichment.csv` - Top enriched sequence motifs
   - `ubr3_translated_sequences.csv` - Your sequences with translations

3. **Console Summary:**
   - Summary statistics printed to console during analysis

## Understanding Your Results

### Enrichment Ratio
- **> 1.5**: Strongly enriched (appears 50% more than expected)
- **0.7 - 1.5**: Normal range
- **< 0.7**: Depleted (appears less than expected)

### What to Look For

1. **Overall Enrichment:** Which amino acids are most enriched overall?
2. **Position-Specific:** Are certain amino acids enriched at specific positions?
3. **Motifs:** Are there recurring patterns (e.g., "KKK" or "RRR")?
4. **Statistical Significance:** Check p-values (< 0.05 is significant)

## Common Questions

**Q: What file format does my data need to be in?**  
A: FASTA (.fasta, .fa), CSV (.csv), or TSV (.tsv, .txt) are all supported.

**Q: What if my file has nucleotide sequences?**  
A: Perfect! The script automatically translates nucleotide sequences to amino acids.

**Q: What if my file already has amino acid sequences?**  
A: The script expects nucleotide sequences. If you have AA sequences, you can modify the script to skip translation.

**Q: Can I change the motif size?**  
A: Yes! Use the `-w` flag: `python ubr3_enrichment_analysis.py file.fasta -w 4` for 4-mers.

**Q: How do I interpret the results?**  
A: Check the console output summary first. Look for amino acids with high enrichment ratios and low p-values. These are significantly enriched in your screen.

## Troubleshooting

**Error: "Could not find nucleotide sequence column"**  
→ Make sure your CSV/TSV has a column named "sequence", "seq", "nt", or similar.

**Error: "ModuleNotFoundError"**  
→ Run `pip install -r requirements.txt`

**Warning: "Could not translate sequence"**  
→ Check that your sequences are valid nucleotides (A, T, G, C only). Sequences with ambiguous bases will be skipped.

## Advanced Usage

For programmatic access and custom analysis:

```python
from ubr3_enrichment_analysis import UBR3EnrichmentAnalyzer

analyzer = UBR3EnrichmentAnalyzer('your_file.fasta', 'results')
analyzer.run_full_analysis()

# Access results
print(analyzer.overall_enrichment)
print(analyzer.motif_enrichment)
```

See `example_usage.py` for more examples.

## Need Help?

Check the full `README.md` for detailed documentation.
