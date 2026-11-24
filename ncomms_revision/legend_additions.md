# NCOMMS-26-051439A: legend additions for the editor's comments

Text to **add to the end of the existing legends**. Keep the current legend text as it is.
Anything in **[brackets]** still has to be filled in.

Sources: Fig. 1d numbers from `ubr3enrichmentlogo/dipeptide_analysis/Source Data Fig 1d.xlsx`;
Fig. 1e numbers from `ubr3_PG_reviewer/data/significance_heatmap_cells.csv` (Fisher columns; figure
position = CSV position − 1, because the figure numbers residues from the first residue after Met).

---

## Fig. 1d: add to the end of the panel d legend

*n* = 104 (Library) and 8 (UBR3 substrates) for PD+PE+PT, and *n* = 109 (Library) and 12 (UBR3 substrates) for GE+GD, where *n* denotes independent N-terminal 24-mer peptide reporters. Library, all screen peptides carrying the indicated N-terminal dipeptides, excluding UBR3 substrates; UBR3 substrates, screen hits carrying the same dipeptides. PSI values are the mean of two control (AAVS) replicates **[biological / technical]**. Box plots show the median (centre line), the 25th and 75th percentiles (box bounds) and the minimum and maximum values (whiskers); no data points were excluded. *P* values, one-sided Brunner–Munzel test, not adjusted for multiple comparisons.

## Fig. 1e: add to the end of the panel e legend

Two-sided Fisher's exact test, one test per class and position. UBR3 substrates carrying an N-terminal P/G–D/E/T motif (*n* = 20 peptides) were compared with motif-carrying library peptides that were not stabilised upon UBR3 loss (*n* = 193 peptides). Colour shows −log10 *P*; numbers show the fold change (frequency in substrates / frequency in non-substrates) for cells with *P* < 0.05; \*, *P* < 0.05; \*\*, *P* < 0.01; \*\*\*, *P* < 0.001; not adjusted for multiple comparisons. Exact *P* values: basic residues at positions 4 (*P* = 0.0011), 6 (*P* = 0.0035) and 7 (*P* = 0.00039); aromatic residues at position 11 (*P* = 0.039). The test was run at every position from 3 to 23 (6 classes × 21 positions = 126 tests), and the panel shows positions 3–11. After Benjamini–Hochberg correction across all 126 tests, only basic residues at position 7 remain significant (*q* = 0.049).

**End of the Fig. 1e addition: choose one option** (the text above uses option A)

- **Option A (recommended), correct across all 126 tests:** "The test was run at every position from 3 to 23 (6 classes × 21 positions = 126 tests), and the panel shows positions 3–11. After Benjamini–Hochberg correction across all 126 tests, only basic residues at position 7 remain significant (*q* = 0.049)."
  Use this unless the Methods say that only positions 3–11 were analysed. The tests were run at all 21 positions, so correcting only over the positions shown can look like choosing the window after seeing the data.
- **Option B, correct across the 54 cells shown:** "After Benjamini–Hochberg correction across the 54 class–position combinations shown (6 classes × 9 positions, 3–11), basic residues at positions 4 (*q* = 0.029) and 7 (*q* = 0.021) remain significant."
  Stronger result (two cells survive instead of one; basic residues at position 6 still do not, *q* = 0.063), but only defensible if positions 3–11 were chosen before the analysis. If you use it, say so in the Methods.

If the legend does not already name the six classes, also add: "Classes: basic (K, H, R), acidic (D, E), aromatic (F, Y, W), aliphatic (G, P, A, M), hydrophobic (V, L, I), polar (S, T, C, N, Q)."

## Figs. 2a, 2b, 2c, 6a, 6c: add to the end of each of these panel legends

Representative of **[N]** independent experiments with similar results.

## Fig. 6c

Nothing to add to the legend: the missing "kDa" unit is fixed in the figure file itself.

## End of the Fig. 1, Fig. 2 and Fig. 6 legends: add

Source data are provided as a Source Data file.

---

## Recommended (not flagged by the editor, but the same rules apply)

- **Fig. 1b**, add: "Numbers show enrichment (frequency in the 54 UBR3 substrates / frequency in all 16,514 library peptides); \*, *P* < 0.05; \*\*, *P* < 0.01; \*\*\*, *P* < 0.001; two-sided Fisher's exact test, not adjusted for multiple comparisons; stars are shown for enriched residues only."
- **Figs. 1f and 6b**, add: "Representative of **[N]** independent experiments with similar results."
