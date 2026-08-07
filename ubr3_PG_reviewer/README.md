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

## Reproducing

```bash
python fetch_annotation.py      # UniProt REST -> data/annotation.json (cached; --refresh to redo)
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
