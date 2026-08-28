# UBR3 N24mer screen — Pro/Gly substrate tables and figures

Response package for the reviewer comment asking for (a) a clean, readable table of the
candidate UBR3 substrates with functional annotation, presented alongside the peptides that
carry the **[P/G]-[E/D]** feature but show **no** stabilisation on UBR3 loss, and (b) evidence
that only a subset of Pro-/Gly-initiated peptides is UBR3-regulated.

**Source data:** `Supplemental Data 1.xlsx` — sheet **(A)** N24mer library (16,514 peptides),
sheet **(B)** UBR3 substrates (54 peptides).

---

## The one-sentence answer to the reviewer

The [P/G]-[E/D] motif enriches for UBR3 substrates **27-fold** (Fisher exact
q = 3 × 10⁻¹⁸; odds ratio 50 vs non-P/G peptides), yet **163 of the 179 library peptides
carrying it (91%) are not stabilised when UBR3 is lost** — so the motif is necessary but
plainly not sufficient, exactly as the reviewer suspected.

| Group | Library peptides | UBR3-stabilised | Rate |
|---|---:|---:|---:|
| `[P/G]-[E/D]` | 179 | 16 | **8.94%** |
| `[P/G]`-other | 1,921 | 10 | 0.52% |
| non-P/G | 14,414 | 28 | 0.19% |

---

## Definitions used throughout

| Term | Definition |
|---|---|
| **pos1 / pos2 / pos3** | Residues 1–3 of the encoded 24-mer. pos1 is Met in **100%** of the library. |
| **Met excision** | Predicted when pos1 = Met and pos2 ∈ {A,C,G,P,S,T,V} (Sherman's rule). pos2 then becomes the true N-terminal residue. |
| **`[P/G]-[E/D]`** | pos2 ∈ {P,G} **and** pos3 ∈ {E,D} — a Pro/Gly N-terminus followed by an acidic residue. |
| **ΔPSI** | PSI(UBR3 KO) − PSI(control KO). Positive = stabilised when UBR3 is lost. |
| **Baseline strata** | Control PSI < 2.6 = "unstable at baseline" (8,112 peptides, 49.1%); ≥ 2.6 = "already stable" (8,402). The cut is the **antimode** of the control-PSI density (stable at 2.59–2.63 across kernel bandwidths 0.20–0.35), not a round number. |
| **Substrate** | The 54 peptides in sheet (B), used verbatim as the author-defined list. |
| **Amino-acid classes** | Acidic D,E · Basic K,R,H · Polar S,T,N,Q,C · Hydrophobic A,V,L,I,M · Aromatic F,W,Y · Special G,P |

**On the substrate definition.** The **54 peptides in sheet (B)** are the substrate set
throughout — in every table, every count and every figure. Any library peptide that is not one
of those 54 is treated as not stabilised. No secondary or intermediate hit list is used
anywhere in this package.

**Caveat on Pro.** Initiator-Met excision at Met-Pro is less efficient than at Met-Gly, so
`Nterm_after_MetAP` is a prediction for the Pro-starting peptides, not a measurement.

---

## The Excel workbook — `UBR3_PG_substrate_tables.xlsx`

19 sheets. Every sheet is frozen at the header row, auto-filtered, and colour-scaled on
ΔPSI / enrichment columns.

| Sheet | Rows | What it is |
|---|---:|---|
| `00_README` | 55 | Methods, definitions, headline numbers, sheet guide — all inside the file. |
| **`01_UBR3_substrates`** | 54 | **The clean readable substrate table the reviewer asked for.** Ranked by mean ΔPSI, with UniProt accession, protein name, one-line function, subcellular localisation, GO biological process, the N-terminal call (pos1/2/3, Met-excision prediction, exposed N-terminus, motif class), both replicate ΔPSI values, PSI per condition, and peptide-level properties (net charge, GRAVY hydropathy, acidic/basic counts). Every one of the 54 has a filled annotation. |
| **`02_PGED_motif_all`** | 179 | **The comparison table the reviewer asked for.** *All* library peptides carrying [P/G]-[E/D], stabilised and non-stabilised side by side, ranked by mean ΔPSI, with a `UBR3_substrate` yes/no flag. The 16 substrates and the 163 non-responders are directly comparable row by row. |
| `03_PG_library_all` | 2,100 | Every peptide with Pro or Gly at position 2, sorted by mean ΔPSI. |
| `04_substrates_PG_only` | 26 | The subset of the 54 substrates carrying P/G at position 2. |
| `05_counts_summary` | 27 | Every count behind Figures 1–2 as a flat table, plus the Mann-Whitney and Fisher test results. |
| `06_enrich_pos2_residue` | 20 | Fisher exact enrichment of each of the 20 residues at position 2, BH-FDR corrected. |
| `07_enrich_motif_class` | 31 | Same test on motif class, on pos3 identity, and on the predicted exposed N-terminus. |
| `08_AAclass_by_position` | 24 | Fraction of each of the six chemical classes at each of the 24 positions, substrates vs library, with log2 enrichment. |
| `09_posAA_freq_substrates` | 24 | 24 × 20 position-frequency matrix, substrates (logo input). |
| `10_posAA_freq_library` | 24 | Same for the library (logo background). |
| `11_posAA_log2_enrichment` | 24 | log2(substrate / library) per position per residue. |
| `12_posAA_fisher_all_cells` | 480 | **All 480 position × residue cells** tested by Fisher exact with BH-FDR across the whole matrix, sorted by q. Only 4 cells reach q < 0.05. |
| `13_posAA_qvalues` | 24 | The same q values as a 24 × 20 lookup matrix. |
| `14_PSI_strata_summary` | 13 | The library split at control PSI = 2.6 — size, substrate rate, mean ΔPSI and motif content of each stratum, plus how many peptides cross the line on UBR3 loss. |
| **`15_motif_x_PSI_strata`** | 8 | **The key control.** Substrate rate for each motif class *within* each baseline-stability stratum, with a Fisher test per stratum. |
| `16_substrates_PSI_transition` | 54 | Each substrate with its control PSI, UBR3-KO PSI, stratum on each side, and whether it crosses PSI 2.6 upward. |
| `17_dPSI_by_control_PSI_bin` | 6 | ΔPSI, substrate count and crossing rate in 0.5-wide bins of control PSI — the ceiling effect in numbers. |
| **`18_PSI_cutoff_robustness`** | 18 | **Proof the cut is not doing the work.** Motif odds ratio across cutoffs 2.0–3.4, a logistic regression using control PSI as a continuous covariate (no cutoff at all), the peaks and valleys of the control-PSI density, and the control-PSI vs ΔPSI correlation. |

---

## The figures — `figures/` (300-dpi PNG + vector PDF)

### Figure 1 — `Figure1_PG_landscape`
*Pro/Gly N-termini and the [P/G]-[E/D] motif in the screen.*

- **A** Position-2 residue composition, library vs substrates (grouped bars, Pro/Gly rows shaded).
  Pro and Gly together are 12.7% of the library but **48% of the substrates**.
- **B** Fisher exact log2 enrichment at position 2, BH-FDR. **Pro (4.3×, q = 3.4 × 10⁻⁴) and
  Gly (3.4×, q = 3.7 × 10⁻⁴) are the only enriched residues**; Ala is significantly depleted
  (0.27×, q = 0.025).
- **C** Nested funnel on a log scale: 16,514 library → 2,100 P/G → 179 [P/G]-[E/D] → 16 stabilised.
- **D** Stabilisation rate by motif class with Wilson 95% CIs. 8.94% vs 0.52% vs 0.19%;
  Fisher exact OR = 50.4, p = 4.1 × 10⁻²⁰.

### Figure 2 — `Figure2_motif_necessary_not_sufficient`
*The central claim, four ways.*

- **A** Every one of the 179 [P/G]-[E/D] peptides as one bar, ranked by mean ΔPSI, substrates in
  orange. The visual point: a short orange head and a long grey tail.
- **B** ECDF of mean ΔPSI by motif class. The whole [P/G]-[E/D] distribution shifts right
  (median 0.140 vs 0.066, Mann-Whitney p = 2.9 × 10⁻⁸) but the distributions still overlap almost
  completely — the motif biases stability without determining it.
- **C** Stacked composition per class with the non-stabilised counts spelled out.
- **D** Replicate-1 vs replicate-2 ΔPSI, showing the 16 substrates clustering together in the
  top-right corner (r = 0.84) — the stabilisation is reproducible, not replicate noise.

### Figure 3 — `Figure3_sequence_logos`
*Sequence logos coloured by amino-acid class.*

- **A** Information-content logo of the 54 substrates.
- **B** The same for the full library, on an identical y-scale — essentially flat.
- **C** **Statistical enrichment logo**: letter height is the signed −log10 Fisher p, and only
  residues with p < 0.05 are drawn, with p = 0.05 and Bonferroni reference lines. Pro and Gly at
  position 2 and Asp at position 3 tower above Bonferroni; everything else is at or near the
  p = 0.05 line.

Positions are plotted **2–24**: position 1 is Met in 100% of both sets, so it carries no
discriminating information and including it flattens every other column.

### Figure 4 — `Figure4_AA_class_analysis`
*Chemical-class architecture across all 24 positions.*

- **A** Stacked class composition at each of the 24 positions for the substrates. Position 2 is
  **48% Special (Pro/Gly)**; position 3 is **37% Acidic**.
- **B** The same for the library — flat from position 4 onward.
- **C** Class enrichment heatmap, **FDR-masked**: cells not reaching q < 0.05 are left beige.
  Four cells survive — Special↑ and Hydrophobic↓ at position 2, Acidic↑ at position 3, Basic↑ at
  position 7.
- **D** Whole-peptide composition averaged over all 24 positions: substrates are
  indistinguishable from the library, confirming the signal is positional, not compositional.

### Figure 5 — `Figure5_position_residue_heatmaps`
*Position × residue, all 24 positions × 20 residues, amino acids grouped by class.*

- **A** Residue frequency in the 54 substrates. Positions 2–3 boxed. Met at position 1 is 100%.
- **B** Enrichment over library, **FDR-masked**. **Only 4 of 480 cells survive q < 0.05**:
  Asp at position 3 (+2.4, q = 0.002), Pro at position 2 (+2.2, q = 0.004), Gly at position 2
  (+1.9, q = 0.006), and Arg at position 7 (+1.7, q = 0.042).

**Note on Arg-7.** It is reported honestly in every figure and sheet rather than smoothed away:
22% of substrates vs 7% of the library. With q = 0.042 it is the weakest of the four survivors
and sits far below the Bonferroni line, so it is best described as a borderline observation, not
a second motif.

### Figure 6 — `Figure6_PSI_baseline_stability`
*Peptides above and below PSI 2.6.*

- **A** Distribution of control PSI across the library, coloured by stratum, with the 54 substrates
  marked as ticks below the axis. The cut splits the library **almost evenly (49% / 51%)**, and
  **74% of the substrates** fall in the unstable half.
- **B** Control PSI vs UBR3-KO PSI for all 16,514 peptides, substrates highlighted, with the PSI 2.6
  lines and the no-change diagonal. **31 of the 40 substrates that start below the line rise above
  it** when UBR3 is lost.
- **C** Mean ΔPSI ± SEM in 0.5-wide bins of control PSI. Stabilisation **peaks in the 2.5–3.0 bin
  (+0.150) and turns negative above 3.5 (−0.067)** — a ceiling effect: peptides already at maximum
  stability have no headroom.
- **D** Substrate rate per stratum: **0.493% below vs 0.167% above** (Fisher OR 3.0,
  p = 3.0 × 10⁻⁴).

> **Caveat worth stating in the manuscript.** This particular contrast is the one result that *is*
> cutoff-sensitive: it is ~16-fold at PSI 3 but 3-fold at PSI 2.6, because the assay ceiling
> suppresses detectable stabilisation at high PSI and a higher cut therefore piles the substrates
> into the low stratum. The motif results below are not cutoff-sensitive.

### Figure 7 — `Figure7_motif_vs_baseline_stability`
*The control analysis: is the motif effect just baseline instability?*

Substrates start unstable, so the motif could in principle be enriched merely because Pro/Gly
N-termini are themselves destabilising. It is not:

- **A** Substrate rate for each motif class **within** each stratum. Fisher OR for [P/G]-[E/D] vs
  non-P/G is **52 among unstable peptides and 52 among stable ones** — essentially identical, which
  is exactly what "no interaction with baseline stability" looks like. (13.41% vs 0.295% in the
  unstable stratum, p = 2.5 × 10⁻¹⁴; 5.15% vs 0.105% in the stable one, p = 3.3 × 10⁻⁷.)
- **B** The control itself: [P/G]-[E/D] is carried by **1.01% of unstable and 1.15% of stable**
  peptides — evenly distributed, so it cannot be acting through baseline stability. Bare Pro/Gly
  *is* skewed (16.6% vs 9.0%), which is why the acidic residue at position 3, not Pro/Gly alone,
  is what carries UBR3 dependence.
- **C** Crossing rate above PSI 2.6, restricted to peptides starting below it. The effect is
  **graded**: substrates 77.5%, motif-bearing non-substrates 21.1%, background 8.9%. Motif-bearing
  peptides that never made the substrate list still respond at more than twice background — a real
  but sub-threshold effect.
- **D** ΔPSI distributions by motif class within the unstable stratum only
  (Mann-Whitney p = 1.7 × 10⁻⁶).

### Figure 8 — `Figure8_PSI_vs_dPSI_and_cutoff_robustness`
*Why stability is measured by control PSI and not by ΔPSI, and whether the cut matters.*

- **A** Control PSI vs ΔPSI for all 16,514 peptides. They are **near-orthogonal**
  (Spearman ρ = −0.20; 2.9% shared variance), so they measure genuinely different things:
  control PSI is how stable a peptide is *before* UBR3 is removed, ΔPSI is how much it *changes*
  when UBR3 goes.
- **B** Kernel density of control PSI. **Three robust modes** (≈1.52, ≈2.25, ≈3.49) with antimodes
  at 1.94 and **2.61** — the boundary between the two unstable modes and the stable one. This is
  where the cut is placed. (A fourth mode near 2.9 appears only at narrow bandwidths and is gone by
  bw 0.20, so it is not reproducible and is not claimed.)
- **C** Odds ratio for [P/G]-[E/D] vs non-P/G within the unstable stratum, across cutoffs from
  2.0 to 3.4. It stays between **42 and 88, with p ≤ 6 × 10⁻⁵ at every cutoff** — the conclusion
  is cutoff-independent.
- **D** **The cutoff-free analysis.** Logistic regression of substrate status on motif class with
  control PSI as a *continuous* covariate: [P/G]-[E/D] **OR 52.8 (p = 3.6 × 10⁻³⁴)**, [P/G]-other
  OR 2.33 (p = 0.023), and control PSI OR 0.61 per unit (p = 0.006). Both effects are real and
  independent; no threshold appears anywhere in the model.

**Why not split on ΔPSI.** ΔPSI is the variable the 54 substrates were *selected* on, so using it
to define the strata and then asking where the substrates fall is circular — every ΔPSI cut
contains all 54 by construction. Control PSI is the correct stratifying variable because it is
measured before UBR3 is perturbed and is nearly independent of ΔPSI.

---

## The arginine correction — read this before reusing any figure from before 2026-08-22

The arginine row of `data/AMINO_ACIDS_PG_motif.csv` is **byte-identical to the Basic-class row**
— it is not arginine at all. The other 19 residue rows match the authoritative file exactly, so
the error is confined to R. Every figure built from that file understated arginine:

| | drawn as | true value |
|---|---:|---:|
| R at position 7 | 2.89× | **5.63×** (7/20 vs 12/193) |
| R at position 8 | 4.14× | **6.43×** (6/20 vs 9/193) |
| R at position 17 | 1.86× | **5.36×** (5/20 vs 9/193) |

Since arginine is the one residue in this analysis with a real signal, the error suppressed
precisely the result worth reporting. **Figures 13, 13b, 14, 14b and 15 have been regenerated**
from the corrected data; any copy of them predating 2026-08-22, in a deck, a draft or a
submitted figure, is wrong and should be replaced.

To stop this recurring, all six logos and all five heatmaps now read one loader,
`pg_motif_data.py`, which is the only module that touches the source sheet.
`data/AMINO_ACIDS_PG_motif.csv` and `data/Categories_PG_motif.csv` are kept for provenance and
are **not read by anything** — the pre-aggregated `Categories` file is where the Basic row leaked
into R in the first place.

---

## Figures 13–15 — logos for positions 4–24

These plot the **P/G-D/E/T motif analysis**: n = 20 substrates against n = 193 non-substrate
controls *carrying the same motif*, i.e. the within-motif comparison, not substrates vs. the
whole library. All of them read `pg_motif_data.py`.

| Figure | Script | Letter height |
|---|---|---|
| `Figure13_foldchange_logo_pos4_24` | `make_foldchange_logo.py` | fold change, all enriched residues |
| `Figure13b_..._significant` | " | fold change, residues at p < 0.05 |
| `Figure14_logo_pos4_24_by_category` | `make_category_logo.py` | same heights, coloured by chemical class |
| `Figure14b_..._significant` | " | same, residues at p < 0.05 |
| `Figure15_significance_logo_pos4_24` | `make_significance_logo.py` | **−log₁₀ p** (chi-square), glyphs coloured by the heatmap ramp |
| `Figure15b_..._fisher` | " | **−log₁₀ p** (Fisher exact), same colouring |

Figure 13 and Figure 15 are the same data under two rulers, and they disagree on purpose.
Under fold change the tallest glyphs are W at 12/14/19 and M at 16 — every one of them 1
substrate vs 1 control, which is just what 1-vs-1 arithmetic gives at n = 20 vs 193. Under
significance the tall glyphs are **R at 7, 8 and 17**. Show Figure 15, not Figure 13, if the
question is "which positions matter".

The significance logos draw only residues reaching p < 0.05. Without that threshold all 20
residues are drawn at every position, most at p ≈ 0.3–0.9, and the stacks reach −log₁₀ p ≈ 8 out
of noise alone. Heights are exact −log₁₀ p (the first version of Figure 15 could only use the
`*/**/***` bin representatives 1.30 / 2.00 / 3.00, because the file it read had no p-values).

**Figure 15 is coloured to match Figures 16–17.** Its glyphs take their colour from the same
`CMAP_SIG` ramp, keyed to the same −log₁₀ p, and carry the same colourbar with the same p ticks —
so a letter that is tall is also dark, and the logo sits beside the heatmaps as one set. The ramp
is entered at its p = 0.05 point, since nothing below that is drawn, and there are no black
outlines on the glyphs: the hairline separating two stacked letters is the same white the heatmaps
use between cells. Chemical class is no longer carried by colour in Figure 15 — that is Figure
14's job on the same stacks — but `build(..., colour_by='class')` restores the class-coloured
version of the figure.

---

## Figures 16–17 — heatmaps coloured by significance (fold change nowhere in the figure)

`make_significance_heatmap.py` → `figures/Figure16*`, `figures/Figure17*`

The Colab heatmap embedded at the bottom of `AMINO ACIDS.csv` colours each position × residue
cell by **fold change** on a locked 0–5 scale. In this dataset that ranking is upside down: the
biggest fold changes come from the *fewest* observations, because with 20 substrates and 193
controls a single peptide in each group already gives 9.65×.

| | Loudest cells under **fold change** | Loudest cells under **significance** |
|---|---|---|
| 1 | W at 9 — 1 vs 0 peptides, FC = ∞, p = 0.82 | **R at 7** — 7/20 vs 12/193, 5.63×, p = 1.7 × 10⁻⁵ |
| 2 | W at 12 / 14 / 19 — 1 vs 1, 9.65×, p = 0.05 | **R at 8** — 6/20 vs 9/193, 6.43×, p = 2.5 × 10⁻⁵ |
| 3 | M at 16, C at 6 — 1–2 vs 1–2, 9.65× | **R at 17** — 5/20 vs 9/193, 5.36×, p = 4.8 × 10⁻⁴ |

Colouring by significance moves arginine to the front and drops the tryptophan cells into the
background where they belong.

**Source.** `data/AminoAcids_PG_motif_with_pvalues.csv` — the sheet of the same workbook that
carries the full 2 × 2 contingency table *and* a chi-square p-value for every cell, rather than
the `* / ** / ***` bins in `AMINO_ACIDS_PG_motif.csv`. n = 20 P/G-D/E/T substrates vs n = 193
non-substrate controls carrying the same motif, positions 4–24.

**Encoding — colour is the p value, the number is the effect size.** Colour is −log₁₀ p on a
single warm ramp — warm paper → straw → gold → orange → red → deep red, deeper red = more
significant, palest = not significant. The ramp lives in `pg_motif_data.py` as `CMAP_SIG`, next to
the loader every figure reads, and Figure 15 colours its glyphs from the same ramp on the same
scale, so the logo and the heatmaps are one visual set and cannot drift apart when one script is
edited alone.

**Every cell at p < 0.05 carries its fold change** (substrates / controls) and the star bin of its
uncorrected p (`*` p < 0.05, `**` p < 0.01, `***` p < 0.001) — and nothing else. The two
quantities are on different channels on purpose: p folds the effect size together with how many
peptides the cell rests on, so a reader given colour alone cannot tell 9.65× resting on one
peptide from 5.63× resting on seven. **The FDR survivors carry no mark inside the matrix** —
neither the old black outline nor a dagger; both put a third thing into a cell that already holds
a colour and a number. The caption names them in words instead. Grey = residue absent from both
groups, no test possible.

> Every cell's exact p and q, both counts, both percentages and the fold change are in
> `data/significance_heatmap_cells.csv` — read the numbers there, quote them in the text. Note that
> p is *not* a relabelled fold change: it folds the fold change together with how many peptides the
> cell rests on, which is exactly why the ranking changes: W at 12 is 9.65× at p = 0.05 (1 vs 1
> peptide) while R at 7 is 5.63× at p = 1.7 × 10⁻⁵ (7/20 vs 12/193). Quoting an effect size in the
> text therefore still needs the fold change from the CSV.

**No class colour panel.** The residue rows are still ordered by chemical class, but the class
colour strip and rotated class labels that used to run down the left edge are gone; the ordering
is stated in words in the caption instead. Nothing else moved.

**Type.** Helvetica throughout — labels, annotations and maths — falling through to Arial where
Helvetica is not installed (metrically identical, and the substitution journals accept). The SVGs
keep live text and declare `Helvetica, Arial, sans-serif`, so they open in real Helvetica on any
machine that has it.

**Why there is no cool arm for depletion.** **No depleted cell reaches p < 0.05 anywhere in this
dataset** — not under either test, not for residues, not for classes. The most significant
depleted cell is Glu at position 8 at p = 0.11 (Acidic at position 8, p = 0.062, for the classes).
A diverging scale would therefore spend half its range colouring sub-threshold noise while
competing for attention with the real signal. Every coloured cell in these figures is an
enrichment — which the caption states outright, so dropping the fold-change annotation costs no
information about direction.

**Workbook correction — chi-square columns I and T (positions 9 and 20).** Every one of the 26
cells in each of those two columns of the `CHI test` block (rows 246–271) disagreed with the
workbook's own contingency tables; the other 490 cells reproduce plain Pearson chi-square exactly.
The observed and expected blocks at those positions are both correct (each position's 20 residue
counts total exactly 20 substrates and 193 controls), so the fault was the CHITEST formula: at
those two columns its range covers only the second column of the 2 × 2 — the *does not carry the
residue* counts — instead of both, which reproduces 32 of the 52 wrong values exactly and inflates
all of them toward p = 1. The 52 cells were recomputed and written back into
`data/AminoAcids_PG_motif_with_pvalues.csv`; the untouched file is kept alongside it as
`…csv.ORIGINAL_before_chi_fix`. **If the source Google Sheet is still in use, fix the formula in
columns I and T there too — this repo's copy is an export.** Two calls changed:

| Cell | Counts | Workbook (wrong) | Corrected chi-square | Fisher |
|---|---|---|---|---|
| W at 9 | 1/20 vs 0/193 | 0.82 — n.s. | **0.0018** | 0.094 — n.s. |
| G at 20 | 5/20 vs 13/193, 3.71× | 0.39 — n.s. | **0.0052** | **0.017** |

After the correction chi-square calls 27 residue cells at p < 0.05 (was 25) against 21 expected by
chance; Fisher is unaffected throughout, since it is recomputed from the counts and never read
from the workbook.

**Two tests, both plotted.** Chi-square is the workbook's own test, but the workbook's expected-count
block shows most cells expect **fewer than 2** substrates, well under the ≥ 5 chi-square needs, so
Fisher exact is recomputed from the same tables and plotted alongside rather than assumed
equivalent. The difference is not cosmetic: chi-square calls 27 cells at p < 0.05, Fisher 18 —
against 21 expected by chance alone.

| Figure | Rows | p-value |
|---|---|---|
| `Figure16_significance_heatmap_residues` | 20 residues × 21 positions | chi-square, as reported in the workbook |
| `Figure16b_significance_heatmap_residues_fisher` | same | Fisher exact, recomputed |
| `Figure16c_significance_heatmap_residues_enriched_only` | same | Fisher, with depleted cells forced to the floor — the drop-in replacement for the Colab figure. Since no depleted cell is significant, it differs from `16b` only among cells that are already n.s. |
| `Figure17_significance_heatmap_categories` | 6 chemical classes × 21 positions | chi-square |
| `Figure17b_significance_heatmap_categories_fisher` | same | Fisher exact |

**What survives correction — state this plainly if the figure is used.** Across the 416 testable
residue cells, **27 reach p < 0.05 under chi-square (18 under Fisher) against 21 expected by
chance alone**. Under chi-square only
**R at positions 7 and 8** survive BH-FDR; under Fisher exact **nothing does**. At the class level
**Basic at position 8** survives both (9/20 vs 21/193, 4.14×, q = 0.049), with Basic at 5 and 7
surviving under chi-square only. Positions 4–24 are therefore best described as carrying **one
borderline arginine/basic signal and no other reproducible per-residue preference** — the motif
itself lives at positions 2–3, which this matrix does not cover.

`data/significance_heatmap_cells.csv` has all 546 cells (420 residue + 126 class) with counts,
percentages, fold change, both p-values and both q-values.

---

## Figures 18–19 — heatmaps coloured by frequency, alphabetical, residues and classes apart

`make_frequency_heatmap.py` → `figures/Figure18*`, `figures/Figure19*`

The third and plainest way to colour the same matrix: not fold change (the Colab original), not
−log₁₀ p (Figures 16–17), but the **percentage of a group that carries that residue — or that
chemical class — at that position**. Same loader, same 21 positions, same warm ramp, so the
whole set still reads as one system.

Two things differ from Figures 16–17:

- **Rows are alphabetical** — `A C D E F G H I K L M N P Q R S T V W Y` for the residues and
  `Acidic … Polar` for the classes — rather than grouped by chemical class. The class is still
  carried by the colour chip beside each row label (keyed to the logo palette, legend under the
  matrix), so the grouping stays recoverable without imposing the row order.
- **Residues and classes are fully separate figures with their own colour scales.** A class is
  the union of 2–5 residues and is necessarily commoner: the residue matrix tops out at 35% and
  the class matrix at 45%, so one shared scale would flatten the residue panel.

| Figure | Rows | Group |
|---|---|---|
| `Figure18_frequency_heatmap_residues` | 20 residues × 21 positions, alphabetical | the **20 substrates**, 0–35% scale |
| `Figure18b_frequency_heatmap_residues_controls` | same | the **193 controls**, same 0–35% scale |
| `Figure19_frequency_heatmap_categories` | 6 chemical classes × 21 positions, alphabetical | the **20 substrates**, 0–45% scale |
| `Figure19b_frequency_heatmap_categories_controls` | same | the **193 controls**, same 0–45% scale |

**Encoding.** Colour and cell text are the same number. Blank = 0%, no peptide in that group has
that residue there. With n = 20 substrates one peptide is 5%, so the substrate panels move in
5-point steps; the control panels are finer (1/193 = 0.52%, printed as `<1`). On the substrate
panels an asterisk marks cells that *also* differ between the two groups at Fisher exact
p < 0.05, so these figures never contradict Figures 16–17 — but nothing about the colour here is
a test.

**Read each substrate panel against its control panel.** A frequency heatmap on its own says
nothing about UBR3: its loudest cells are partly just the commonest amino acids. That is exactly
what `18b` / `19b` are for, and the comparison is not cosmetic — the top of each panel is a
different list.

| | Loudest cells among the **substrates** | Loudest cells among the **controls** |
|---|---|---|
| Residues | R at 7 (35%), R at 8 / G at 15 (30%), R at 17 · G at 20 · S at 18 · S at 23 (25%) | L at 21 (15%), L at 18 (14%), A at 11 · L at 23 (13%) — Leu is simply the commonest residue in this library |
| Classes | Basic at 7 and 8 · Polar at 22 (45%), Basic at 5 · Aliphatic at 15, 19, 20, 21 · Polar at 18 (40%) | Aliphatic at 8 and 11 (30%), Aliphatic at 15 (29%), Polar at 5 and 20 (28%) |

Of the 8 most frequent residue cells in the substrates, all 8 also reach Fisher p < 0.05 — but of
the 19 substrate cells at ≥ 20%, only 14 do, and **none of the residue cells survives FDR**
(1 class cell does: Basic at 8). The caveats recorded for Figures 16–17 apply here unchanged:
these panels describe composition, and the enrichment claim still rests on one borderline
arginine/basic signal.

`data/frequency_heatmap_cells.csv` has all 546 cells with both groups' counts and percentages
alongside fold change, both p-values and both q-values.

---

## Reproducing

```bash
python fetch_annotation.py      # UniProt REST -> data/annotation.json (cached; --refresh to redo)

# positions 4-24, P/G-D/E/T motif analysis -- all five read pg_motif_data.py
python pg_motif_data.py               # loader self-test: cell counts, what survives FDR
python make_foldchange_logo.py        # -> figures/Figure13*
python make_category_logo.py          # -> figures/Figure14*, data/significant_*.csv
python make_significance_logo.py      # -> figures/Figure15*
python make_significance_heatmap.py   # -> figures/Figure16*, Figure17*, data/significance_heatmap_cells.csv
python make_frequency_heatmap.py      # -> figures/Figure18*, Figure19*, data/frequency_heatmap_cells.csv

python build_workbook.py        # -> UBR3_PG_substrate_tables.xlsx
python make_figures.py          # -> figures/*.png and *.pdf   (8 composed figures)
python make_panels.py           # -> panels/*.png              (29 standalone panels)
python build_pptx.py            # -> UBR3_figures_explained.pptx (40 slides)
```

| File | Role |
|---|---|
| `ubr3_core.py` | Loading, N-terminal/motif annotation, class definitions, PSI strata, position matrices, the Fisher/FDR enrichment engine. Everything imports it, so tables, figures and slides cannot drift apart. |
| `ubr3_panels.py` | Every figure panel as a function that draws into a supplied Axes, plus the palette. Shared by the composed figures and the slide deck. |
| `fetch_annotation.py` | UniProt functional annotation for the 54 substrate genes. |
| `build_workbook.py` | Builds and styles the 19-sheet workbook. |
| `make_figures.py` | Composes the panels into the 8 multi-panel manuscript figures. |
| `make_panels.py` | Renders each panel on its own, large, for the deck. |
| `pptx_text.py` | All slide copy — definitions, per-figure blurbs, per-panel explanations. Edit here to reword the deck. |
| `build_pptx.py` | Assembles the deck; layout adapts to panel shape. |

## The slide deck — `UBR3_figures_explained.pptx`

40 slides, 16:9. Title → definitions and methods → then, for each of the 8 figures, a section
slide showing the composed figure followed by **one slide per panel**, each with the panel
rendered large and a written explanation covering *what is plotted*, *how it was computed*,
*what it shows*, and *what to watch out for*. Ends with a summary-of-conclusions slide.

Panels wider than 1.9:1 (the logos and the class-composition bars) get the image across the top
with the explanation in two columns beneath; all others get the image on the left and the
explanation on the right.

**Statistics.** Two-sided Fisher exact throughout; Benjamini-Hochberg FDR applied across the
residues/classes/cells tested within each analysis (20 for position 2, 480 for the full
position × residue matrix). Group ΔPSI comparisons use the two-sided Mann-Whitney U test.
Proportion confidence intervals are Wilson intervals.

Figures are light/print mode by design. The categorical palette was validated for
colour-vision-deficient separation; class identity is additionally carried by the letter glyph
in the logos and by direct labels elsewhere, never by colour alone.
