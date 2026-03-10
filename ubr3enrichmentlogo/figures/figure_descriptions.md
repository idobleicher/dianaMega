# Figure Descriptions — UBR3 N-terminal Peptide Screen

**Dataset:** 16,514 peptides (24-mer, all starting with Met) screened for UBR3-dependent degradation. 54 best hits identified. All analyses compare the amino acid composition of the 54 hits to the full library at positions 2 and 3 (first two residues after Met).

---

### Fig 1 — Enrichment bar chart (positions 2 and 3)
Enrichment = frequency in hits / frequency in library, for each amino acid. Sorted by enrichment. Dashed line at 1.0 (no change).
- Pos 2: Pro 4.3x, Gly 3.4x, Thr 2.4x
- Pos 3: Asp 4.9x, Glu 2.2x, Pro 2.1x

### Fig 1b — Enrichment heatmap
Same values as Fig 1, displayed as a 20x2 heatmap (amino acids x positions). Sorted by position 2 enrichment.

### Fig 1c — Dipeptide enrichment heatmap (20x20)
Enrichment of all 400 possible position 2-3 dipeptide combinations. Top hits: PD 40x, PT 36x, GD 35x, GE 32x, TN 31x.

### Fig 2 — Frequency comparison (hits vs library)
Raw frequency (%) of each amino acid at positions 2 and 3. Gold = library, dark red = hits. Sorted by hit frequency.

### Fig 3 — Top dipeptide enrichment bars
Top 20 dipeptides at positions 2-3 ranked by enrichment. PD leads at ~40x.

---

### Fig 4 — Sequence logo (positions 2-3, best hits)
Information content logo (bits) from 54 hit sequences. Letter height = conservation. Colors by biochemical category: dark red = acidic (D,E), brick red = basic (R,K,H), orange = nonpolar (G,A,V,L,I,P,F,M,W), gold = polar (S,T,C,Y,N,Q).

### Fig 4b — Sequence logo with Met (positions 1-3)
Same as Fig 4, with Met at position 1 included. Met shows maximum information content (fully conserved).

### Fig 5 / 5b — Extended logos (positions 2-6 and 1-6)
Logo extended to 5 positions. Information content drops sharply after position 3, confirming specificity is restricted to the first two residues after Met.

### Logo color variants
Same logo in 8 alternative color schemes: canonical warm, classic (red/blue/black/green), red spectrum, sunset, rose, terracotta, coral/blush, charge-focused.

---

### Fig 6 — PSI by dipeptide position in the 24-mer
For PD, PE, and GE motifs: mean PSI (AAVS, average of two replicates) calculated for all library peptides containing each motif at each of the 23 positions in the 24-mer. Top panel = mean PSI per position. Bottom panel = sample size. All three motifs show lowest PSI (least stable) at positions 2-3.

### Fig 7 — PSI position trend (line plot)
Same data as Fig 6, as line plots. Filtered for n>=5. Positions 2-3 highlighted. Annotated with n and mean PSI at key positions.

### Fig 8 — PSI: positions 2-3 vs other positions
For each motif, peptides split into N-terminal (positions 2-3) vs downstream (positions 4-24). Mean PSI compared as paired bars with sample sizes. All three motifs show lower PSI at positions 2-3 (~0.1-0.3 units lower).

---

### Fig 9 — Enrichment by amino acid category
Amino acids grouped by property (acidic, basic, nonpolar, polar). Frequency compared between library and hits. Key result: acidic residues 3.3x enriched at position 3; basic residues depleted to 0.3x.

### Fig 10 — Frequency difference (iceLogo-style)
Delta frequency (hits minus library, in %). Red = enriched, gray = depleted. Pos 2: Gly +18%, Pro +17%, Ala -15%. Pos 3: Asp +18%.

### Fig 11 — Statistical enrichment (Fisher's exact test)
Same as Fig 1 with one-sided Fisher's exact test p-values. Significant (p<0.05) in dark red. Pos 2: Pro and Gly (p<0.001). Pos 3: Asp (p<0.001).

### Fig 12 — Position weight matrix (positions 2-8)
Probability of each amino acid at positions 2-8 in the 54 hits. Heatmap format. Positions 2-3 show clear preferences (G 26%, D 22%); positions 4+ approach uniform (~5%).

### Fig 13 — Consensus vs anti-consensus dipeptides
Top 10 most enriched and top 10 least enriched dipeptides at positions 2-3. Horizontal bar charts.

---

### Fig 14 — Enrichment heatmap (positions 2-12)
Enrichment for each amino acid across 11 positions. Strong signal confined to positions 2-3; downstream positions near 1.0.

### Fig 15 — Log2 enrichment heatmap (positions 2-12)
Same as Fig 14 on log2 scale. Diverging color map: red = enriched, blue = depleted, white = neutral.

### Fig 16 — Category enrichment heatmap (positions 2-12)
Enrichment by amino acid category across positions 2-12. Acidic 3.3x at position 3 is the dominant signal. Basic depleted at positions 2-3 but enriched at 7-8.

### Fig 17 — Frequency heatmap (library vs hits)
Side-by-side heatmaps of raw amino acid frequency (%) at positions 2-7. Library shows uniform distribution; hits show concentrated hot spots at positions 2-3.

### Fig 18 — Focused dipeptide heatmap
Dipeptide enrichment for the 6 most enriched position-2 residues (P, G, T, Q, E, D) against all 20 position-3 residues. Reveals sub-motifs: P/G prefer D/E; T prefers N; E prefers W.

---

## Methods

Enrichment = frequency in hits (n=54) / frequency in library (n=16,514), calculated per amino acid and per dipeptide at positions 2-3. Statistical significance by one-sided Fisher's exact test. Sequence logos generated with logomaker (Shannon information content, bits). PSI (AAVS) = mean of PSI-293a and PSI-293b. Position-dependent stability assessed by comparing mean PSI for peptides with PD/PE/GE motifs at positions 2-3 vs downstream positions. Python 3, pandas, numpy, matplotlib, logomaker, scipy.
