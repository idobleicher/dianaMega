# UBR3 N-Terminal Enrichment Analysis - Figure Guide
## Focus: First 3 Positions After Met (Positions 2, 3, 4)

**Dataset:** 53 hits (one per gene) vs 16,514 full screen sequences  
**Core N-Terminal Motif:** M - [P/G] - [D/E] - [Y/H]  
**Analysis Date:** January 2026

---

## Figure 1 - Comprehensive Enrichment Overview (Full-Length, Positions 2-24)

**File:** `Figure_1_Comprehensive_Enrichment_Overview.png`

**What it shows:** A 5-panel master overview of the entire enrichment analysis across all 24 positions. This is the "big picture" figure.

**Panel breakdown:**
- **Panel A (top):** Sequence logo of the 91 UBR3 **hits** (positions 2-24). Letter height = information content (bits). Shows strong P and G at position 2, D/E at position 3. Consensus: M-PDLSVLLLSLSGTGGGSSGGTRL.
- **Panel B:** Sequence logo of the **full screen library** (16,514 sequences). Dominated by Alanine at position 2 and Leucine throughout. Consensus: M-AALLLLLLLLLLLLLLLLLLLLL.
- **Panel C (KEY PLOT):** **Enrichment-weighted logo** where height = frequency x enrichment ratio. This highlights what is specifically selected by UBR3 vs background. Position 2 (P), position 3 (D), and position 4 (Y) stand out clearly at the N-terminus.
- **Panel D:** Log2 enrichment heatmap across all positions. Red = enriched in hits, Blue = depleted.
- **Panel E:** Position discrimination scores. Identifies positions 2, 3, 6, 12, and 15 as the top 5 most discriminating.

**Key finding for positions 2-4:** The enrichment-weighted logo (Panel C) clearly shows Proline (pos 2), Aspartate (pos 3), and Tyrosine/Serine (pos 4) as the dominant N-terminal signal, well above the background noise.

**Recommended use:** Overview figure for talks or as a supplementary figure showing the complete analysis.

---

## Figure 2 - Hits vs Library Statistical Comparison (With P-values)

**File:** `Figure_2_Hits_vs_Library_Statistical_Comparison.png`

**What it shows:** Side-by-side bar charts comparing the frequency of the top enriched amino acids in **hits (solid bars)** vs **full library (hatched bars)** at the first 5 N-terminal positions. Statistical significance (Fisher's exact test) is shown with asterisks above bars.

**Position +1 (= position 2 in sequence):**
- **P (Proline): ~19% in hits vs ~5% in library (p < 0.001, \*\*\*)**
- **G (Glycine): ~17% in hits vs ~7.5% in library (p < 0.01, \*\*)**
- T (Threonine): ~7.5% in both (not significant)

**Position +2 (= position 3):**
- **D (Aspartate): ~13% in hits vs ~4.5% in library (p < 0.001, \*\*\*)**
- **E (Glutamate): ~13% in hits vs ~6% in library (p < 0.05, \*)**
- T (Threonine): ~9% vs ~6% (not significant)

**Position +3 (= position 4):**
- **Y (Tyrosine): ~7% in hits vs ~2% in library (p < 0.01, \*\*)**
- I (Isoleucine): ~7% vs ~3% (trending)
- R (Arginine): ~12% vs ~7.5% (not significant)

**Position +4 (= position 5):**
- V (Valine): enriched (p < 0.05, *)
- I, C, D: moderate enrichment

**Key finding:** This is the **only figure with formal statistical testing (p-values)**. It confirms that P at position 2 and D at position 3 are the most statistically significant enrichments (p < 0.001). Y at position 4 is also significant (p < 0.01).

**Recommended use:** Main figure or supplementary for publication - provides statistical validation.

---

## Figure 3 - Annotated Enrichment Heatmap (First 5 Residues)

**File:** `Figure_3_Annotated_Enrichment_Heatmap.png`

**What it shows:** A heatmap of **Log2 enrichment ratios** for all 20 amino acids across positions 1-5. Color scale: **Red = enriched in hits, Blue = depleted in hits**. Numeric values are annotated only where |Log2| > 1.0 (i.e., >2-fold enrichment or <0.5-fold depletion).

**Key enrichments (red cells with positive values):**
- **Position 2, P (Proline): +2.1** (= 4.39x enrichment)
- **Position 3, D (Aspartate): +2.2** (= 4.57x enrichment) -- HIGHEST
- **Position 4, Y (Tyrosine): +2.1** (= 4.37x enrichment)
- Position 2, G (Glycine): +1.7
- Position 4, H (Histidine): +1.6
- Position 4, I (Isoleucine): +1.4
- Position 5, I (Isoleucine): +1.6
- Position 5, K (Lysine): +1.3
- Position 5, C (Cysteine): +1.3
- Position 3, E (Glutamate): +1.2
- Position 2, T (Threonine): +1.3

**Key depletions (blue cells with negative values):**
- Many amino acids are strongly depleted (values like -7 to -9), indicating they are essentially absent in the hits at these positions (e.g., F, M, N, C at positions 2-3).

**Key finding for positions 2-4:** The three strongest enrichments form a diagonal: P at pos 2 (+2.1), D at pos 3 (+2.2), Y at pos 4 (+2.1). This is the core **M-P-D-Y recognition motif**.

**Recommended use:** Excellent main figure for publication - clean, professional, shows both enrichment and depletion with exact values.

---

## Figure 4 - Frequency Heatmap: Hits vs Screen Side-by-Side

**File:** `Figure_4_Frequency_Heatmap_Hits_vs_Screen.png`

**What it shows:** Two stacked heatmaps showing **raw amino acid frequencies** (not enrichment ratios). Panel A shows the 53 hits, Panel B shows the 16,514 full screen library. Color intensity = frequency (0-40%).

**What stands out in Panel A (Hits):**
- Position 1: Met at 100% (dark, expected)
- **Position 2: P and G appear as bright orange/red spots (~20-25% each)** -- very different from the screen
- **Position 3: D and E appear as bright spots (~15-20% each)**
- Position 4: More spread out, but Y is visible

**What stands out in Panel B (Screen):**
- Position 1: Met at 100%
- Position 2: A (Alanine) dominates (~20%) -- typical library bias
- Positions 3-5: Much more uniform, no strong hotspots except L

**Key finding for positions 2-4:** The visual contrast between Panel A and B makes it immediately obvious that the hits have a very different amino acid composition at positions 2-3 compared to the background library. The "hot spots" for P, G (pos 2) and D, E (pos 3) in the hits panel are completely absent in the screen panel.

**Recommended use:** Great figure for showing the raw data difference between hits and library before any statistical processing. Very intuitive for audiences.

---

## Figure 5 - Publication Enrichment Profile (Top 8 AAs per Position)

**File:** `Figure_5_Publication_Enrichment_Profile.png`

**What it shows:** Four panels (positions 2-5), each showing the **top 8 amino acids ranked by enrichment ratio**. Bars are color-coded: dark red (>3x), red (2-3x), orange (1.5-2x), yellow (1-1.5x). Exact enrichment values are labeled above each bar. A text box summarizes the key enrichments.

**Position 2:**
- **P: 4.39x** (dark red, highest)
- **G: 3.25x** (dark red)
- T: 2.42x, Q: 2.10x (red)
- E: 1.33x (yellow)
- D, L, S below 1x (depleted)

**Position 3:**
- **D: 4.57x** (dark red, HIGHEST OVERALL)
- E: 2.27x, P: 2.17x (red)
- T: 1.66x (orange)
- N, L, W, S: around 1x

**Position 4:**
- **Y: 4.37x** (dark red, highest)
- H: 2.96x (red)
- I: 2.58x, D: 2.13x (red)
- R: 1.87x (orange)
- K: 1.48x, V: 1.31x (yellow/orange)

**Position 5:**
- I: 2.96x, K: 2.49x, C: 2.49x (red)
- D: 2.01x (red)
- Less dramatic than positions 2-4

**Key finding for positions 2-4:** The top enrichment values are remarkably consistent: P (4.39x), D (4.57x), and Y (4.37x) all show ~4.4-fold enrichment. This suggests they are equally important for recognition. The text box in the figure summarizes: Pos 3 D (4.57x) = HIGHEST, Pos 4 Y (4.37x) = KEY.

**Recommended use:** Best figure for publication main panel - clean, labeled with exact values, publication-ready.

---

## Figure 6 - Top 5 Enriched AAs Combined Bar Graph

**File:** `Figure_6_Top5_Enriched_Combined_Bars.png`

**What it shows:** A single combined bar chart showing the **top 5 most enriched amino acids at each position (2-5)**, grouped by position. Bars decrease in height within each group. Dashed lines mark 2x and 3x enrichment thresholds. Each bar is labeled with its amino acid letter.

**Position 2 group:** P (4.4x) > G (3.3x) > T (2.4x) > Q (2.1x) > E (1.3x)
**Position 3 group:** D (4.6x) > E (2.3x) > P (2.2x) > T (1.7x) > N (1.5x)
**Position 4 group:** Y (4.4x) > H (3.0x) > I (2.6x) > D (2.1x) > R (1.9x)
**Position 5 group:** I (3.0x) > K (2.5x) > C (2.5x) > D (2.0x) > (1.3x)

**Key finding for positions 2-4:** The three tallest bars across the entire chart are P (pos 2), D (pos 3), and Y (pos 4) -- all reaching ~4.4-4.6x. They clearly tower above all other amino acids, confirming the M-P-D-Y motif as the dominant signal.

**Recommended use:** Good overview figure for presentations - visually striking, easy to read at a glance.

---

## Figure 7 - Sequence Logo (Positions 2-10, Top AAs >5% Frequency)

**File:** `Figure_7_Sequence_Logo_Pos2_10.png`

**What it shows:** A traditional **sequence logo plot** for positions 2-10 from the 53 hit sequences. Letter height = frequency x information content (bits). Only amino acids with >5% frequency are shown for clarity. Letters are color-coded by chemical property: Red = acidic (D, E), Blue = basic (K, R), Light blue = His (H), Black = hydrophobic, Green = polar, Orange = G, P.

**Position 2:** G (orange, tallest) and P (orange) dominate, with T (green) and E (red) below. Total information content ~1.2 bits.
**Position 3:** D (red, tallest) and E (red) dominate, with P (orange) below. ~1.1 bits.
**Position 4:** Lower overall conservation (~0.6 bits). L (black), Y, and mixed residues.
**Positions 5-6:** Low conservation, mixed.
**Position 7:** G (orange) prominent with R (blue) and L (black).
**Positions 8-9:** R (blue) is tall -- Arginine enrichment in the hydrophobic core.
**Position 10:** Low conservation, mixed hydrophobic.

**Key finding for positions 2-4:** Positions 2 and 3 have the tallest letter stacks (highest information content), confirming they are the most conserved and important positions. The signal drops significantly at position 4 and beyond, suggesting positions 2-3 are the primary recognition determinants.

**Recommended use:** Classic sequence logo format -- immediately recognizable to the molecular biology community. Great for presentations and papers.

---

## Figure 8 - PSI Heatmap: Unstable vs Stable Peptide-GFP Fusions

**File:** `Figure_8_PSI_Unstable_vs_Stable_Heatmap.png`

**What it shows:** Side-by-side heatmaps comparing amino acid composition in **unstable (PSI <= 3.6)** vs **stable (PSI > 3.6)** peptide-GFP fusions across all 24 positions. Color = Log2(frequency ratio), where Red = enriched relative to full library, Blue = depleted.

**Panel A - Unstable peptides (UBR3 substrates, degraded):**
- **Position 2-3: D and E are DEPLETED (blue)** -- surprising!
- **Position 2: C is enriched (red)**
- Overall mild enrichment patterns, weaker signal
- F (Phenylalanine) shows enrichment at later positions

**Panel B - Stable peptides (NOT degraded by UBR3):**
- **Much stronger patterns overall** (more intense colors)
- **Position 2-4: D and E are STRONGLY ENRICHED (dark red)**
- **K (Lysine) strongly depleted (dark blue) across many positions**
- R, S, P show mixed patterns
- Clear distinction between N-terminal and C-terminal regions

**Key finding for positions 2-4:** The PSI analysis reveals that D/E at positions 2-3 are associated with **stable** (non-degraded) proteins, which seems counterintuitive given the enrichment analysis. This suggests the relationship between the N-terminal motif and protein stability is more complex than simple degradation -- the motif may play a role in protein folding or stability that is separate from UBR3-mediated degradation.

**Recommended use:** Important for discussing the relationship between sequence composition and protein stability. Good supplementary figure.

---

## Figure 9 - Motif Occurrence in Unstable vs Stable Peptides (Statistics)

**File:** `Figure_9_Motif_Occurrence_Statistics.png`

**What it shows:** A 4-panel analysis of how 5 key dipeptide motifs (PD, PE, GE, GD, PT) are distributed between unstable and stable peptide-GFP fusions.

**Panel A - Motif frequency in sequences:**
- All 5 motifs (PD, PE, GE, GD, PT) are **more frequent in STABLE peptides**
- PE: 14.5% stable vs 8.0% unstable
- GE: 14.4% stable vs 7.6% unstable
- GD: 9.9% stable vs 6.5% unstable
- PT: 8.2% stable vs 5.4% unstable
- PD: 6.0% stable vs 5.2% unstable (smallest difference)

**Panel B - Log2 enrichment (unstable/stable):**
- All bars are negative = all motifs are **depleted in unstable peptides**
- PE: 0.55x (p = 2.05e-17, \*\*\*)
- GE: 0.53x (p = 8.54e-20, \*\*\*)
- GD: 0.66x (p = 1.27e-06, \*\*\*)
- PT: 0.66x (p = 1.57e-05, \*\*\*)
- PD: 0.85x (p = 0.601, not significant)

**Panel C - Total motif occurrences:**
- Unstable peptides have far more total occurrences (larger group) but lower per-sequence frequency

**Panel D - Summary statistics table:**
- Clear statistical significance for PE, GE, GD, PT (p < 0.001)

**Key finding for positions 2-4:** The N-terminal dipeptide motifs that are enriched in the UBR3 hits (PD, PE, GE, GD) are paradoxically more common in **stable** peptides at the whole-library level. This strongly suggests that these motifs contribute to protein stability (possibly through favorable N-terminal interactions), and UBR3 may specifically recognize and counteract this stabilizing signal, or these motifs may have dual roles depending on context.

**Recommended use:** Critical for the biological interpretation -- supports a nuanced view of the N-terminal motif. Main figure or important supplementary.

---

## Summary Table

| Figure | Description | Key Data for Positions 2-4 | Best Use |
|--------|-------------|---------------------------|----------|
| **1** | Comprehensive 5-panel overview | Enrichment-weighted logo shows P, D, Y clearly | Overview / Supplementary |
| **2** | Hits vs Library with p-values | P (p<0.001), D (p<0.001), Y (p<0.01) | Main figure (statistics) |
| **3** | Annotated enrichment heatmap | P (+2.1), D (+2.2), Y (+2.1) Log2 values | Main figure (publication) |
| **4** | Frequency heatmap hits vs screen | Visual raw frequency difference | Supplementary |
| **5** | Publication enrichment bars | P (4.39x), D (4.57x), Y (4.37x) | Main figure (clean) |
| **6** | Combined top 5 bar graph | Same values, combined view | Presentations |
| **7** | Sequence logo (pos 2-10) | G/P at pos 2, D/E at pos 3 dominate | Main figure (classic) |
| **8** | PSI unstable vs stable heatmap | D/E enriched in stable, not unstable | Supplementary |
| **9** | Motif occurrence statistics | PE, GE depleted in unstable (p<0.001) | Main / Supplementary |

---

## Recommended Figures for Publication

### Main Figure (pick 2-3):
- **Figure 3** (Annotated heatmap) + **Figure 5** (Publication bars) + **Figure 2** (Statistics)

### Supplementary:
- **Figure 1** (Full overview) + **Figure 7** (Sequence logo) + **Figure 8** (PSI heatmap)

### For Presentations:
- **Figure 6** (Combined bars) + **Figure 7** (Sequence logo)

---

## The Core Motif (Positions 2-4)

```
Position:    2         3         4
             |         |         |
Enriched:    P(4.39x)  D(4.57x)  Y(4.37x)
             G(3.25x)  E(2.27x)  H(2.96x)
             T(2.42x)  P(2.17x)  I(2.58x)
             |         |         |
Function:    KINK      CHARGE    AROMATIC
             (turn)    (salt     (binding
             inducer)  bridge)   pocket)
```
