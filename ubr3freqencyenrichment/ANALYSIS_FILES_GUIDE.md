# UBR3 Screen Analysis - Files Guide

## 📊 Complete Analysis Results

### Quick Summary
- **Hits Consensus:** `MPDLSVLLLSLSGTGGGSSGGTRL`
- **Full Screen Consensus:** `MAALLLLLLLLLLLLLLLLLLLLL`
- **Key Finding:** Glycine-rich region (positions 15-21) is highly enriched in hits (up to 4.62x)

---

## 📁 Directory Structure

### `logo_results/` - Main Results (All Sequences)
**Key Files:**
- 🎨 **`ubr3_sequence_logo.png`** - Information content logo plot showing consensus
- 🎨 **`ubr3_sequence_heatmap.png`** - Heatmap of position-specific frequencies
- 🎨 **`comparison_hits_vs_screen.png`** ⭐ - Side-by-side comparison (MUST SEE!)
- 🎨 **`aa_enrichment_comparison.png`** - Amino acid enrichment analysis
- 🎨 **`ubr3_top_motifs_3mer.png`** - Top 3-mer motifs
- 🎨 **`ubr3_top_motifs_4mer.png`** - Top 4-mer motifs
- 🎨 **`ubr3_top_motifs_5mer.png`** - Top 5-mer motifs
- 📄 **`consensus_sequence.txt`** - Detailed consensus breakdown
- 📊 **`position_weight_matrix.csv`** - Position weight matrix
- 📊 **`information_content_matrix.csv`** - Information content values
- 📊 **`position_specific_enrichment.csv`** - Per-position enrichment data
- 📊 **`overall_aa_enrichment_comparison.csv`** - Hits vs screen AA enrichment

### `logo_results_hits/` - UBR3 Hits Only (91 sequences)
- All logo plots and analyses specific to enriched sequences
- Shows the pattern that UBR3 actually recognizes

### `logo_results_full_screen/` - Full Library (16,514 sequences)
- Background distribution from complete screen
- Shows what was NOT selected (mostly hydrophobic LLL repeats)

### `results/` - Original Enrichment Analysis
- Basic enrichment analysis from initial run
- Position heatmaps
- Overall AA enrichment

### `comparative_results/` - Comparative Analysis
- Hits vs Screen enrichment heatmaps
- Position-specific enrichment from comparative analysis

---

## 🎯 Most Important Files to View

### 1. **`logo_results/comparison_hits_vs_screen.png`** ⭐⭐⭐
   - **What:** Side-by-side logo comparison
   - **Shows:** Clear visual difference between hits and full screen
   - **Key insight:** Glycine-rich region is obvious in hits, absent in screen

### 2. **`logo_results_hits/ubr3_sequence_logo.png`** ⭐⭐
   - **What:** Sequence logo of enriched hits
   - **Shows:** Consensus motif with information content
   - **Key insight:** Positions 15-21 (GGGGSS G) stand out

### 3. **`logo_results/aa_enrichment_comparison.png`** ⭐⭐
   - **What:** Amino acid enrichment bar plot
   - **Shows:** Which AAs are enriched/depleted in hits
   - **Key insight:** G, H, R, T, V enriched; A, L, W depleted

### 4. **`LOGO_ANALYSIS_SUMMARY.md`** ⭐⭐⭐
   - **What:** Complete written analysis with interpretation
   - **Contains:** All findings, structural interpretation, recommendations
   - **Length:** ~250 lines of detailed analysis

### 5. **`logo_results_hits/ubr3_top_motifs_3mer.png`**
   - **What:** Common 3-amino acid patterns
   - **Shows:** MPT, MGE, GGG, GLL as top motifs
   - **Key insight:** Recurring patterns in hits

---

## 📈 Key Data Files

### For Further Analysis:

1. **`logo_results/position_weight_matrix.csv`**
   - Frequency of each amino acid at each position
   - Use for: Custom visualizations, statistical tests

2. **`logo_results/overall_aa_enrichment_comparison.csv`**
   - Amino acid enrichment ratios (hits/screen)
   - Columns: AA, hits_freq, screen_freq, enrichment_ratio, log2_enrichment
   - Use for: Identifying discriminating amino acids

3. **`logo_results_hits/position_specific_enrichment.csv`**
   - Per-position, per-amino acid enrichment
   - Use for: Finding critical positions

4. **`logo_results/consensus_sequence.txt`**
   - Position-by-position consensus with conservation scores
   - Use for: Understanding variability at each position

---

## 🔬 Analysis Scripts

### To Regenerate or Modify:

1. **`generate_logo_plots.py`**
   - Main script for logo plot generation
   - Usage: `python generate_logo_plots.py <input_file> -o <output_dir>`
   - Options: `-s` (start position), `-m` (max positions)

2. **`create_comparison_logos.py`**
   - Creates comparison visualizations
   - Usage: `python create_comparison_logos.py`
   - Requires: PWM files from logo_results_hits/ and logo_results_full_screen/

3. **`ubr3_enrichment_analysis.py`**
   - Original enrichment analysis script
   - Usage: `python ubr3_enrichment_analysis.py <input_file>`
   - Options: `-c` (comparative mode), `-s` (screen file)

---

## 📊 Analysis Workflow

```
Input Data
    ├── UBR3 Nt screen.xlsx (full library, 16,514 sequences)
    └── ubr3_best (1).xlsx (hits, 91 sequences)
           │
           ↓
    Enrichment Analysis (ubr3_enrichment_analysis.py)
           │
           ↓
    Logo Plot Generation (generate_logo_plots.py)
           │
           ├── logo_results_hits/
           ├── logo_results_full_screen/
           └── logo_results/
           │
           ↓
    Comparison Analysis (create_comparison_logos.py)
           │
           ↓
    Final Results
           ├── Sequence logos (PNG files)
           ├── Motif analyses (PNG files)
           ├── Data matrices (CSV files)
           └── Summary reports (MD/TXT files)
```

---

## 🎨 Visualization Types Generated

### 1. **Sequence Logos (Information Content)**
   - Height = information content (bits)
   - Shows conservation at each position
   - Color = amino acid properties (chemistry/hydrophobicity)

### 2. **Heatmaps**
   - Rows = amino acids, Columns = positions
   - Color intensity = frequency
   - Good for seeing overall patterns

### 3. **Motif Bar Plots**
   - Top recurring 3-mer, 4-mer, 5-mer sequences
   - Shows common patterns in hits

### 4. **Enrichment Plots**
   - Bar plots showing log2(hits/screen)
   - Identifies discriminating amino acids

### 5. **Comparison Plots**
   - Side-by-side logos
   - Difference heatmaps (hits - screen)
   - Position enrichment scores

---

## 📖 How to Interpret Results

### Reading Sequence Logos:
- **Tall letters** = high conservation (important for function)
- **Short letters** = variable (less important)
- **Stack height** = information content (0-4.3 bits for 20 amino acids)

### Reading Heatmaps:
- **Red/warm colors** = high frequency
- **Blue/cool colors** = low frequency
- **White** = amino acid not present

### Reading Enrichment Plots:
- **Positive log2** = enriched in hits (Red bars)
- **Negative log2** = depleted in hits (Blue bars)
- **|log2| > 1** = >2-fold change (significant)

---

## 🔑 Key Findings Summary

### 1. **Consensus Motif for UBR3 Recognition:**
```
M-P-[D/E]-[hydrophobic core]-[GLYCINE-RICH FLEXIBLE REGION]-[R/T/L]
│ │   │          │                      │                        │
1 2   3        4-10                  11-21                    22-24
```

### 2. **Most Discriminating Positions:**
- **Position 2:** P (proline) - 3.96x enriched
- **Position 15:** G (glycine) - 4.62x enriched ⭐ HIGHEST
- **Position 3:** D/E (acidic) - 2.64x enriched
- **Positions 15-21:** GGGG SSG pattern (glycine-rich)

### 3. **Enriched Amino Acids:**
- **G (Glycine):** 1.34x - flexibility
- **H (Histidine):** 1.25x - charge/coordination
- **R (Arginine):** 1.23x - positive charge
- **T (Threonine):** 1.21x - polar
- **V (Valine):** 1.14x - hydrophobic

### 4. **Depleted Amino Acids:**
- **W (Tryptophan):** 0.74x - bulky aromatic
- **A (Alanine):** 0.76x - small hydrophobic
- **L (Leucine):** 0.81x - hydrophobic (surprising!)

### 5. **Top Motifs in Hits:**
- **MPT, MGE** (N-terminal start motifs)
- **GGG, GGGG, GGGGG** (flexible linkers) ⭐
- **RRR** (charged regions)
- **LLL** (hydrophobic cores)

---

## 💡 Recommendations for Next Steps

### 1. **Validate Key Positions:**
   - Mutate P2 → A (remove kink)
   - Mutate G15-21 → A (remove flexibility)
   - Test binding to UBR3

### 2. **Search Proteome:**
   - Pattern: `M-P-[DE]-X{6,8}-G{3,5}`
   - Find endogenous UBR3 substrates

### 3. **Structural Studies:**
   - Co-crystallize UBR3 + consensus peptide
   - NMR to show G-rich region flexibility

### 4. **Biochemistry:**
   - Measure KD for different variants
   - Ubiquitination assays

---

## 📞 Questions?

If you need to:
- **Regenerate plots:** Run scripts in `logo_results*/`
- **Change analysis parameters:** Modify and rerun `generate_logo_plots.py`
- **Export different formats:** Modify matplotlib saving options in scripts
- **Analyze different positions:** Use `-s` and `-m` flags

---

**Analysis Completed:** January 17, 2026  
**Total Sequences Analyzed:** 16,605  
**Total Plots Generated:** 30+  
**Total Data Files:** 15+
