#!/usr/bin/env python
"""Slide copy for the UBR3 panel deck.

FIGURES  - one entry per composed figure: title + what the figure argues.
PANELS   - one entry per panel: slide title + explanation sections.
           Each section is (heading, body). Keep bodies self-contained: a reader
           should be able to land on any single slide and understand the panel.
"""

DEFINITIONS = [
    ('The assay',
     'GPS (Global Protein Stability). Each of 16,514 human N-terminal 24-mers is fused to a '
     'fluorescent reporter, cells are FACS-sorted into 4 stability bins, and the bin distribution '
     'of each peptide is read out by Illumina sequencing.'),
    ('PSI',
     'Protein Stability Index - the read-weighted mean bin, running ~1 to 4. Higher = more stable. '
     '"Control PSI" is the mean of the two control-KO replicates.'),
    ('ΔPSI',
     'PSI(UBR3 KO) minus PSI(control KO). Positive means the peptide becomes more stable when UBR3 '
     'is lost, i.e. UBR3 normally destabilises it. Two independent replicates.'),
    ('The substrate set',
     'The 54 peptides in sheet (B) of Supplemental Data 1 are the candidate UBR3 substrates, used '
     'verbatim everywhere in this deck. Any library peptide not among those 54 is treated as not '
     'stabilised.'),
    ('Initiator-Met excision',
     'Position 1 is Met in 100% of the library. MetAP removes it when position 2 has a small side '
     'chain (A, C, G, P, S, T, V), so position 2 becomes the residue actually exposed at the '
     'N-terminus. Excision at Met-Pro is less efficient than at Met-Gly, so the call is a '
     'prediction for Pro-starting peptides.'),
    ('The motif under test',
     '[P/G]-[E/D] means Pro or Gly at position 2 AND Asp or Glu at position 3 - a Pro/Gly-initiated '
     'N-terminus immediately followed by an acidic residue. "[P/G]-other" carries Pro/Gly without '
     'the acidic residue; "non-P/G" is everything else.'),
    ('Baseline strata',
     'Control PSI < 2.6 = unstable at baseline (8,112 peptides); >= 2.6 = already stable (8,402). '
     'The cut is the antimode of the control-PSI density, not a round number. A sweep across '
     'cutoffs from 2.0 to 3.4 (workbook sheet 18) shows nothing depends on this choice.'),
    ('Amino-acid classes',
     'Acidic D,E | Basic K,R,H | Polar S,T,N,Q,C | Hydrophobic A,V,L,I,M | Aromatic F,W,Y | '
     'Special G,P.'),
    ('Statistics',
     'Two-sided Fisher exact throughout, with Benjamini-Hochberg FDR applied across the residues, '
     'classes or cells tested in each analysis. Group ΔPSI comparisons use the two-sided '
     'Mann-Whitney U test. Proportion error bars are Wilson 95% intervals.'),
]

FIGURES = {
    '2': ('The motif is necessary but not sufficient',
          'The reviewer\'s central point, argued four ways. The [P/G]-[E/D] motif enriches for '
          'substrates 27-fold, yet 163 of the 179 peptides carrying it are not stabilised at all. '
          'The motif biases stability without determining it.'),
    '3': ('Sequence logos of the substrate N-termini',
          'Where along the 24-mer the constraint actually lives. Panels A-B compare substrates '
          'against the library background on an identical scale; panel C converts the comparison '
          'into a statistical enrichment logo so only real signal is drawn.'),
    '4': ('Chemical-class architecture across all 24 positions',
          'The same question at the level of chemical classes rather than individual residues. '
          'Panels A-B show composition position by position, C tests every class x position cell '
          'for significance, and D checks whole-peptide composition as a control.'),
    '5': ('Position-by-residue heatmaps',
          'The full 24 x 20 picture. Panel A is raw substrate frequency; panel B is enrichment '
          'over the library with non-significant cells masked, which reduces the whole dataset to '
          'four surviving cells.'),
    '6': ('Baseline stability: peptides above and below PSI 2.6',
          'Splits the library by how stable each peptide is BEFORE UBR3 is removed. A degron '
          'substrate should start unstable and gain stability when its E3 ligase is lost - which '
          'is what the substrates do. Also exposes the assay\'s ceiling effect.'),
    '10': ('Positions 4-24: substrates vs non-substrates, both motif-bearing',
           'The comparison that controls for the motif. Every other analysis contrasts the 54 '
           'substrates with the whole library, so positions 2-3 dominate. Here both groups carry '
           'M-[P/G]-[E/D], making positions 1-3 identical by construction, so anything that '
           'separates them must lie downstream. 16 substrates vs 163 non-substrates.'),
    '11': ('Class composition of the Pro/Gly substrates against the library',
           'The stacked-composition view applied to the 26 substrates carrying Pro or Gly at '
           'position 2. Panel B supplies the P/G-matched background and C the whole library, so '
           'the position-2 selection can be separated from anything the substrates do differently.'),
    '9': ('Motif-bearing peptides classified by stability and substrate status',
          'Crosses the [P/G]-[E/D] motif with the baseline stability call to give a clean four-cell '
          'classification, and annotates the 16 substrates within it. Answers directly: how many '
          'motif-bearing peptides are stable, how many unstable, and which of each are UBR3 '
          'substrates.'),
}

PANELS = {
    # ---------------------------------------------------------------- FIG 10
    '10A': ('Chemical class by position, downstream of the motif', [
        ('What is plotted',
         'A heatmap of log2(substrate / non-substrate) for each of the six chemical classes at each '
         'position from 4 to 24. Red = over-represented in the 16 motif-bearing substrates, blue = '
         'over-represented in the 163 motif-bearing non-substrates. A boxed cell with an asterisk '
         'passed FDR at q < 0.05.'),
        ('Why this comparison and not the earlier ones',
         'Every other figure compares substrates against the whole library, where positions 2 and 3 '
         'carry almost all the signal and swamp everything else. Here BOTH groups carry '
         'M-[P/G]-[E/D], so positions 1-3 are identical between them and cannot contribute. Any '
         'difference that appears is genuinely downstream.'),
        ('What it shows',
         'One cell of 126 survives correction: Basic at position 5, present in 50% of substrates '
         '(8/16) versus 9.8% of non-substrates (16/163), q = 0.026. Basic also trends up at '
         'positions 7, 8 and 14 without reaching significance individually - a coherent direction '
         'rather than scattered noise.'),
        ('Watch out for',
         'With 16 substrates, one peptide moves a cell by 6.25%, so individual cells are volatile. '
         'Read the repeated Basic signal across nearby positions, not any one square. Panel C tests '
         'that aggregate properly.'),
    ]),
    '10B': ('No individual residue distinguishes them', [
        ('What is plotted',
         'A volcano over all 420 position x residue cells in the window: x is log2(substrate / '
         'non-substrate), y is -log10 of the Fisher exact p value. The dotted line is nominal '
         'p = 0.05. The four smallest p values are labelled by residue and position.'),
        ('What it shows',
         'Nothing survives FDR - not one of the 420 cells. Only 11 reach even nominal p < 0.05, '
         'which is FEWER than the ~21 you would expect by chance alone at that threshold. The '
         'strongest single cell is Arg at position 7 (37.5% vs 4.9%, p = 0.0004, q = 0.16).'),
        ('The correct interpretation',
         'At the resolution of "which residue sits at which position", motif-bearing substrates and '
         'motif-bearing non-substrates are indistinguishable. There is no second sequence motif '
         'downstream waiting to be found - at least not one this dataset can detect with 16 '
         'substrates.'),
        ('Why the asymmetry',
         'Almost all points sit right of zero because a residue absent from all 16 substrates gives '
         'a bounded negative log ratio, while one enriched in a few gives a large positive one. That '
         'is a small-sample artefact of the axis, not a biological signal.'),
    ]),
    '10C': ('Only basic residues differ, and only in aggregate', [
        ('What is plotted',
         'For each chemical class, the mean number of residues of that class per peptide summed '
         'across all 21 downstream positions. Orange = the 16 substrates, grey = the 163 '
         'non-substrates. Error bars are SEM; the p value above each pair is Mann-Whitney on the '
         'per-peptide counts.'),
        ('Why aggregate rather than per position',
         'Summing over 21 positions turns 21 underpowered tests into one well-powered one. If '
         'substrates carry more basic residues but at no fixed position, this is the only test that '
         'can see it - and panel B shows that is exactly the situation.'),
        ('What it shows',
         'Basic residues: 4.19 per substrate versus 2.91 per non-substrate, p = 0.014. Every other '
         'class is flat: Acidic p = 0.12, Polar p = 0.11, Aromatic p = 0.14, Hydrophobic p = 0.33, '
         'Special p = 0.80.'),
        ('Watch out for',
         'These six tests are not FDR-corrected against each other; at q < 0.05 across six, Basic '
         'would sit right at the boundary. The net-charge test in panel D is the stronger and more '
         'interpretable version of the same observation.'),
    ]),
    '10D': ('Substrates carry a net positive charge downstream', [
        ('What is plotted',
         'Net charge summed over positions 4-24 for every peptide (Lys and Arg +1, Asp and Glu -1, '
         'His +0.1), substrates versus non-substrates. Boxes are median and quartiles; every '
         'individual peptide is overlaid as a dot.'),
        ('What it shows',
         'Motif-bearing substrates average a net charge of +2.16 across the window; motif-bearing '
         'non-substrates average -0.22. Mann-Whitney p = 0.0021. This is the single strongest '
         'downstream difference in the whole comparison.'),
        ('Why it is the better statistic',
         'It combines the basic enrichment and the acidic depletion - which individually reach '
         'p = 0.014 and p = 0.12 - into one physically meaningful quantity, and it needs no '
         'position-by-position testing at all.'),
        ('How to read it biologically, and carefully',
         'It is consistent with a positively charged region downstream of the degron contributing to '
         'recognition, which would fit a ligase engaging an acidic surface. But it rests on 16 '
         'peptides, the distributions overlap heavily, and it is an observational correlation in a '
         'reporter assay. It should be presented as a hypothesis worth testing, not an established '
         'requirement.'),
    ]),

    # ---------------------------------------------------------------- FIG 11
    '11A': ('Class composition per position - the 26 Pro/Gly substrates', [
        ('What is plotted',
         'For each of the 24 positions, a 100% stacked bar giving the fraction of the 26 P/G '
         'substrates whose residue at that position belongs to each chemical class.'),
        ('What it shows',
         'Position 1 is entirely Hydrophobic (the invariant initiator Met) and position 2 entirely '
         'Special, because Pro/Gly at position 2 is the selection criterion - 14 Gly and 12 Pro. '
         'Position 3 is 42% Acidic. From position 4 the bars become an ordinary mixture.'),
        ('The right way to read it',
         'Positions 1 and 2 are constrained by definition and carry no information. The only place '
         'these substrates can reveal something is position 3 onward - compare against panel B, '
         'which has the identical position-2 constraint.'),
        ('Watch out for',
         'With 26 peptides one substrate is 3.8% of a bar, so bar-to-bar wobble after position 3 is '
         'sampling noise. Panel D applies the significance test.'),
    ]),
    '11B': ('Class composition per position - the 2,100 P/G library peptides', [
        ('What is plotted',
         'The same calculation over every library peptide carrying Pro or Gly at position 2, '
         'regardless of UBR3 dependence. This is the motif-matched background.'),
        ('Why this is the right background',
         'Comparing P/G substrates against the WHOLE library confounds two things: that they were '
         'selected for Pro/Gly, and whatever else makes them substrates. Holding position 2 fixed '
         'on both sides removes the first and leaves only the second.'),
        ('What it shows',
         'Positions 1 and 2 are constrained exactly as in panel A. From position 3 onward the '
         'composition is flat and typical - roughly 36% Hydrophobic, 22% Polar, 13% Special, '
         '12% Basic, 10% Acidic, 7% Aromatic.'),
        ('Watch out for',
         'Position 3 already carries a small Acidic bump here relative to the whole library, so '
         'even the matched background is not perfectly neutral at that position.'),
    ]),
    '11C': ('Class composition per position - the whole library', [
        ('What is plotted',
         'All 16,514 library peptides, the global background. Included so the figure is '
         'self-contained and the effect of the position-2 selection is visible directly.'),
        ('What it shows',
         'Position 1 is Met. Position 2 is unconstrained here and looks like ordinary protein '
         'sequence, which is precisely the difference from panels A and B. Everything from '
         'position 3 onward is flat.'),
        ('What comparing A with C would tell you',
         'It would show a huge Special signal at position 2 - but that is the selection criterion '
         'restated, not a finding. That is why panel D compares A with B instead.'),
        ('Watch out for',
         'With n = 16,514 these proportions are precise to well under a percentage point, so treat '
         'this panel as a fixed reference rather than an estimate.'),
    ]),
    '11D': ('P/G substrates vs the P/G-matched library, FDR-masked', [
        ('What is plotted',
         'log2(P/G substrates / P/G library) for every class at every position, with cells failing '
         'q < 0.05 left beige. Fisher exact per cell, Benjamini-Hochberg across all 144.'),
        ('What it shows',
         'Exactly one cell survives: Acidic at position 3, log2 = +2.8. Once the position-2 '
         'selection is held constant on both sides, the acidic residue at position 3 is the ONLY '
         'single-position compositional feature separating the substrates from their own '
         'background.'),
        ('Why this matters',
         'It is the cleanest statement of the central result. Pro/Gly alone is not what makes a '
         'peptide a UBR3 substrate; among peptides that all carry Pro/Gly, only the acidic residue '
         'at position 3 distinguishes the ones that are.'),
        ('The limit of the method',
         'This tests composition one position at a time. It cannot see a property distributed '
         'across many positions - and there is one: the downstream net-charge difference in the 10 '
         'series, which this panel is blind to by construction.'),
    ]),

    # ---------------------------------------------------------------- FIG 2
    '2A': ('163 of 179 motif-bearing peptides are not stabilised', [
        ('What is plotted',
         'Every one of the 179 library peptides carrying [P/G]-[E/D], one bar each, ranked left to '
         'right by mean ΔPSI. Orange bars are the 16 that are UBR3 substrates; grey bars are '
         'the 163 that are not.'),
        ('What it shows',
         'A short orange head and a very long grey tail. The motif-bearing peptides span the whole '
         'ΔPSI range, from +1.63 down to -0.57 - some are actually destabilised when UBR3 is '
         'lost. Carrying the motif tells you very little about what an individual peptide will do.'),
        ('Why it answers the reviewer',
         'This is the most direct possible answer to "only a subset of Pro-/Gly-initiated peptides '
         'are UBR3-regulated". 91% of motif-bearing peptides show no stabilisation whatsoever.'),
        ('Watch out for',
         'Orange and grey bars interleave near the top: some non-substrates have a higher mean '
         'ΔPSI than some substrates. Substrate status is the published call, not a pure rank '
         'on this axis, so the two orderings are not identical.'),
    ]),
    '2B': ('The whole motif class shifts, but only slightly', [
        ('What is plotted',
         'Empirical cumulative distribution functions of mean ΔPSI, one curve per motif class. '
         'For any x value the curve height is the fraction of that group at or below it. A curve '
         'sitting further right means a group that is more stabilised overall.'),
        ('What it shows',
         'The [P/G]-[E/D] curve (orange) sits to the right of the other two across most of the '
         'range: median ΔPSI 0.140 versus 0.066 for non-P/G, Mann-Whitney p = 2.9 x 10^-8. But '
         'the three curves overlap almost completely - the shift is real and highly significant, '
         'yet small.'),
        ('The nuance worth stating',
         'The motif is not a switch that some peptides flip and others do not. It shifts the entire '
         'distribution a little. Substrate status is the top tail of a continuum, not a separate '
         'population.'),
        ('Watch out for',
         'A very small p value here reflects the enormous n (14,414 non-P/G peptides), not a large '
         'effect. Read the separation between curves for effect size and the p value only for '
         'whether the shift is real.'),
    ]),
    '2C': ('Motif enrichment is 17-fold, but the motif is not sufficient', [
        ('What is plotted',
         'Each motif class as a 100% stacked bar: orange = UBR3-stabilised, pale grey = not '
         'stabilised. The percentage inside the orange segment is the substrate rate; the label to '
         'the right gives the absolute number of non-responders.'),
        ('What it shows',
         'The orange slivers are 8.94%, 0.52% and 0.19%. Even in the best group, 163 of 179 peptides '
         'are unaffected by UBR3 loss. In the other two groups the orange is barely visible: 1,911 '
         'and 14,386 non-responders respectively.'),
        ('Why this framing',
         'Panel 1D shows the same numbers as a rate comparison and makes the enrichment look '
         'dramatic. This panel shows the complement, and makes the insufficiency obvious. Both are '
         'true and the paper needs both.'),
        ('Watch out for',
         'Do not read the orange widths as comparable areas - the group sizes differ by two orders '
         'of magnitude, so these are proportions within each row, not counts.'),
    ]),
    '2D': ('Stabilisation is reproducible, not replicate noise', [
        ('What is plotted',
         'ΔPSI in replicate 1 against replicate 2. Pale dots are all other library peptides; '
         'dark grey dots are [P/G]-[E/D] peptides that are not substrates; large orange dots are the '
         '16 [P/G]-[E/D] substrates, with the top six labelled by gene.'),
        ('What it shows',
         'The 16 substrates cluster together in the upper-right, meaning both independent replicates '
         'agree that they gained stability. Replicate correlation across the [P/G]-[E/D] set is '
         'r = 0.84.'),
        ('Why it is here',
         'It rules out the most obvious alternative explanation - that the substrates are simply '
         'peptides that happened to fluctuate upward in one replicate. Agreement between two '
         'independently sorted replicates is not something noise produces.'),
        ('Watch out for',
         'The dense pale cloud is centred slightly above zero, so a small positive ΔPSI is the '
         'norm across the library. Judge the substrates against that cloud, not against zero.'),
    ]),

    # ---------------------------------------------------------------- FIG 3
    '3A': ('Sequence logo of the 54 candidate substrates', [
        ('What is plotted',
         'A sequence logo over positions 2 to 24. At each position the stack height is the '
         'information content in bits, and each letter\'s height is its share of that. Letters are '
         'coloured by chemical class (legend top right). The shaded block marks positions 2-3.'),
        ('Why position 1 is missing',
         'Position 1 is Met in 100% of both the substrates and the library, so it carries no '
         'information that distinguishes them. Including it produces a single tower that flattens '
         'every other column to invisibility.'),
        ('What it shows',
         'Almost all the height in the whole logo is at positions 2 and 3: Gly and Pro stacked at '
         'position 2, Asp and Glu at position 3. From position 4 rightward the stacks collapse to '
         'a low, even carpet.'),
        ('Watch out for',
         'Information content is computed against a uniform 1/20 background, so it measures how '
         'constrained a position is - not how different it is from the library. Panel B supplies '
         'the background and panel C the actual comparison.'),
    ]),
    '3B': ('Sequence logo of the full library (background)', [
        ('What is plotted',
         'The identical calculation applied to all 16,514 library peptides, drawn on exactly the '
         'same y-scale as panel A so the two can be compared directly by eye.'),
        ('What it shows',
         'Essentially flat. There is a small stack at position 2 - Ala, Ser, Gly are common there '
         'because of MetAP-driven codon usage - and nothing elsewhere. This is what "no constraint" '
         'looks like in this assay.'),
        ('Why it matters',
         'It establishes that the towers in panel A are a property of the substrates and not of the '
         'library construction. Without this panel a reviewer could argue the position-2 signal is '
         'just how the library was built.'),
        ('Watch out for',
         'The residual position-2 structure here is exactly why enrichment must be measured against '
         'this background rather than against a uniform expectation - which is what panel C does.'),
    ]),
    '3C': ('Statistical enrichment logo: what actually distinguishes substrates', [
        ('What is plotted',
         'Letter height is the signed -log10 of the Fisher exact p value for that residue at that '
         'position, substrates versus library. Letters above the line are over-represented in '
         'substrates; letters below are depleted. Only residues reaching p < 0.05 are drawn at all. '
         'The two red lines mark p = 0.05 and the Bonferroni threshold.'),
        ('Why not a simple ratio',
         'A plain log2 substrate/library logo stacks a dozen depleted residues at every position '
         'and becomes unreadable, with the y-axis running to -13. Plotting significance instead '
         'keeps only what the data can actually support.'),
        ('What it shows',
         'Pro and Gly at position 2 and Asp at position 3 tower far above the Bonferroni line. '
         'Everything else in the entire 23-position window sits at or just above the p = 0.05 line - '
         'the level you would expect from chance across this many tests.'),
        ('Watch out for',
         'Letters between the two red lines are nominally significant but do not survive correction '
         'for multiple testing. Treat them as noise; the FDR-corrected version is panel 5B.'),
    ]),

    # ---------------------------------------------------------------- FIG 4
    '4A': ('Class composition at each position - the 54 substrates', [
        ('What is plotted',
         'For each of the 24 positions, a 100% stacked bar giving the fraction of substrates whose '
         'residue at that position belongs to each of six chemical classes.'),
        ('What it shows',
         'Position 1 is entirely Hydrophobic - that is the invariant initiator Met. Position 2 is '
         '48% Special, i.e. Pro or Gly. Position 3 is 37% Acidic. From position 4 onward the bars '
         'settle into a mixture that looks like any random protein sequence.'),
        ('Why classes rather than residues',
         'Grouping 20 residues into six classes raises the counts per cell, which matters with only '
         '54 substrates, and it tests a chemically meaningful hypothesis - "an acidic residue" '
         'rather than "specifically Asp".'),
        ('Watch out for',
         'With n = 54, one peptide is 1.9% of a bar. Small differences between adjacent positions '
         'are not meaningful; panel C applies the significance test.'),
    ]),
    '4B': ('Class composition at each position - the full library', [
        ('What is plotted',
         'The same stacked-bar calculation applied to all 16,514 library peptides - the background '
         'against which panel A should be read.'),
        ('What it shows',
         'Position 1 is again all Met. Position 2 shows mild structure, then from position 4 onward '
         'the composition is essentially constant: roughly 36% Hydrophobic, 22% Polar, 13% Special, '
         '12% Basic, 10% Acidic, 7% Aromatic.'),
        ('Why it matters',
         'The flatness confirms there is no positional bias built into the library beyond the '
         'N-terminal region. Any structure in panel A after position 3 would therefore be real - '
         'and there is none.'),
        ('Watch out for',
         'With n = 16,514 these proportions are precise to well under a percentage point, so this '
         'panel can be treated as a fixed reference rather than an estimate.'),
    ]),
    '4C': ('Class enrichment by position, FDR-masked', [
        ('What is plotted',
         'A heatmap of log2(substrate / library) for each chemical class at each position. Red = '
         'enriched in substrates, blue = depleted. Beige cells are those that do NOT reach q < 0.05 '
         'and are deliberately left blank.'),
        ('How it was computed',
         'A Fisher exact test for every class x position cell, with Benjamini-Hochberg correction '
         'applied across the whole matrix. Ratios use Laplace smoothing so an empty cell cannot '
         'produce an infinite value.'),
        ('What it shows',
         'Four cells survive out of 144. At position 2: Special enriched (+1.9) and Hydrophobic '
         'depleted (-1.6). At position 3: Acidic enriched (+1.7). At position 7: Basic enriched '
         '(+1.2). Everything else is beige.'),
        ('Why the masking',
         'With 54 substrates, an unmasked version of this heatmap shows vivid red and blue cells '
         'scattered to position 24 that are pure sampling noise. Masking is what keeps this panel '
         'honest.'),
    ]),
    '4D': ('Whole-peptide composition is unchanged', [
        ('What is plotted',
         'Class composition averaged over all 24 positions of every peptide. Coloured bars are the '
         '54 substrates (coloured by class); grey bars are the library.'),
        ('What it shows',
         'The two are nearly identical: Acidic 10.1 vs 10.2, Basic 14.7 vs 12.2, Polar 22.5 vs 21.6, '
         'Hydrophobic 31.9 vs 35.8, Aromatic 6.2 vs 6.8, Special 14.6 vs 13.3.'),
        ('Why this is the control',
         'It rules out a bulk-composition explanation. If substrates were simply unusually acidic or '
         'unusually Pro/Gly-rich proteins overall, that would show here. It does not - the signal is '
         'strictly positional, confined to residues 2 and 3.'),
        ('Watch out for',
         'The Hydrophobic difference (31.9 vs 35.8) is the largest gap and goes in the direction of '
         'substrates being slightly less hydrophobic. With n = 54 it is not individually significant '
         'and no single position drives it.'),
    ]),

    # ---------------------------------------------------------------- FIG 5
    '5A': ('Residue frequency in the 54 substrates', [
        ('What is plotted',
         'A 24 x 20 heatmap: the percentage of substrates carrying each amino acid at each position. '
         'Rows are grouped by chemical class with horizontal rules between classes and row labels '
         'coloured to match. The black box outlines positions 2-3.'),
        ('What it shows',
         'Inside the box: Gly 26% and Pro 22% at position 2; Asp 22% and Glu 15% at position 3. The '
         'rest of the map is a pale, even wash, which is what an unconstrained sequence looks like.'),
        ('Reading the scale',
         'The colour ceiling is set at 28% so the motif cells are distinguishable. Met at position 1 '
         'is 100% and therefore saturates - it is labelled directly rather than left to the colour '
         'scale.'),
        ('Watch out for',
         'These are raw frequencies with no background correction. A residue can look dark simply '
         'because it is common in human proteins. Panel B is the version that corrects for that.'),
    ]),
    '5B': ('Enrichment over library background, FDR-masked', [
        ('What is plotted',
         'The same 24 x 20 grid, now showing log2(substrate / library) with every cell that fails '
         'q < 0.05 left beige. Red = enriched, blue = depleted.'),
        ('How it was computed',
         'All 480 position x residue cells tested by Fisher exact, with Benjamini-Hochberg '
         'correction across the entire matrix. The full table is sheet 12 of the workbook.'),
        ('What it shows',
         'Four cells survive out of 480. Asp at position 3 (+2.4, q = 0.002), Pro at position 2 '
         '(+2.2, q = 0.004), Gly at position 2 (+1.9, q = 0.006), and Arg at position 7 (+1.7, '
         'q = 0.042). The entire rest of the peptide is indistinguishable from background.'),
        ('The honest caveat',
         'Arg at position 7 is real but weak - 22% of substrates versus 7% of the library, at the '
         'edge of significance and far below Bonferroni. It is reported rather than hidden, but it '
         'should be described as a borderline observation, not a second motif.'),
    ]),

    # ---------------------------------------------------------------- FIG 6
    '6A': ('49% of the library starts below PSI 2.6', [
        ('What is plotted',
         'A histogram of control PSI across the library, coloured by stratum: purple below the cut, '
         'green above. Orange ticks below the axis mark where each of the 54 substrates falls.'),
        ('What it shows',
         'The distribution is clearly multimodal, with peaks near 1.5, 2.25 and 3.5. The cut at 2.6 '
         'splits it almost evenly (49% / 51%). 74% of the substrates sit in the unstable half.'),
        ('Why control PSI and not ΔPSI',
         'Control PSI is measured before UBR3 is touched, so it describes the peptide\'s intrinsic '
         'stability. ΔPSI is the variable the substrates were selected on, so '
         'stratifying on it would be circular.'),
        ('Watch out for',
         'The substrate ticks are not uniformly spread through the unstable half; they concentrate '
         'between roughly 2.0 and 2.9. Very unstable peptides (PSI < 1.5) are also poor substrates, '
         'presumably because they are degraded by other routes.'),
    ]),
    '6B': ('31 of the 40 substrates below the line cross it', [
        ('What is plotted',
         'Control PSI on the x-axis against UBR3-KO PSI on the y-axis for all 16,514 peptides. Pale '
         'dots are the library; orange dots are the 54 substrates. The dotted diagonal is "no '
         'change"; the dashed lines mark PSI 2.6. The shaded upper-left quadrant is where a peptide '
         'starts unstable and ends stable.'),
        ('What it shows',
         'Substrates sit well above the diagonal - they gain stability. 31 of the 40 substrates that '
         'start below PSI 2.6 rise above it in UBR3 KO, versus a library background of 8.9%.'),
        ('Why the diagonal matters',
         'Distance above the diagonal is ΔPSI. Plotting the two PSI values instead of ΔPSI '
         'alone shows both how much a peptide moved and where it started, which is what the '
         'stratified analysis depends on.'),
        ('Watch out for',
         'The pale cloud is thickest along the diagonal at both extremes because peptides pinned at '
         'PSI 1 or PSI 4 cannot move. That is the ceiling effect quantified in panel C.'),
    ]),
    '6C': ('Stabilisation peaks near PSI 2.75 and reverses above 3.5', [
        ('What is plotted',
         'Mean ΔPSI with standard-error bars, in 0.5-wide bins of control PSI. Bar colour '
         'follows the stratum of the bin midpoint; the dashed line marks the actual cut at 2.6. '
         'Bin sizes are printed along the axis.'),
        ('What it shows',
         'A clear inverted-U. Mean ΔPSI rises from +0.080 in the 1-1.5 bin to a maximum of '
         '+0.150 in the 2.5-3.0 bin, then falls to +0.085 and finally turns NEGATIVE at -0.067 in '
         'the 3.5-4.0 bin.'),
        ('Why it matters',
         'This is an assay ceiling, not biology. A peptide already at PSI 3.9 has almost no room to '
         'gain, so stabilisation becomes undetectable and measurement noise dominates. Any claim '
         'about substrates being concentrated at low PSI has to be read against this constraint.'),
        ('Watch out for',
         'Error bars are SEM across thousands of peptides and are therefore tiny; they show the mean '
         'is precisely estimated, not that individual peptides behave consistently.'),
    ]),
    '6D': ('Substrates are 3x more common below the line', [
        ('What is plotted',
         'The percentage of peptides in each baseline stratum that are UBR3 substrates, with Wilson '
         '95% confidence intervals and counts printed below each bar.'),
        ('What it shows',
         '0.493% of unstable peptides (40 of 8,112) are substrates versus 0.167% of stable ones '
         '(14 of 8,402) - a 3.0-fold difference, Fisher OR 3.0, p = 3.0 x 10^-4.'),
        ('The caveat that must travel with this panel',
         'This is the one result in the whole package that is genuinely sensitive to where the cut '
         'is placed: at PSI 3 the same comparison gives 16-fold. The reason is the ceiling effect '
         'in panel C - a higher cut mechanically herds substrates into the low stratum. The 3-fold '
         'value at the data-driven cut is the more honest number.'),
        ('What is NOT cutoff-sensitive',
         'The motif results are stable across every cut from 2.0 to 3.4 '
         '(workbook sheet 18). Only this baseline-stability contrast moves.'),
    ]),



    # ---------------------------------------------------------------- FIG 9
    '9A': ('The 179 motif-bearing peptides, classified two ways', [
        ('What is plotted',
         'A mosaic of the 179 [P/G]-[E/D] peptides. Column WIDTH is the share of peptides in that '
         'baseline-stability class; block HEIGHT within a column is the share that are UBR3 '
         'substrates (orange) versus not (grey). Every number is an actual peptide count.'),
        ('The four cells',
         'Unstable and a substrate: 11. Unstable and not: 71. Stable and a substrate: 5. Stable and '
         'not: 92. Totals: 82 unstable, 97 stable, 16 substrates, 163 non-substrates.'),
        ('What it shows',
         'The motif-bearing peptides split almost evenly between the two stability classes - 82 '
         'versus 97 - so carrying the motif does not by itself make a peptide unstable. The '
         'substrate share is 2.6x higher in the unstable column (13.4% vs 5.2%).'),
        ('Do not over-read the 2.6x',
         'Within these 179 peptides that difference is NOT statistically significant: Fisher '
         'OR 2.85, p = 0.067, because it rests on 11 substrates versus 5. Library-wide the same '
         'comparison IS significant (OR 2.97, p = 3 x 10^-4). Report the direction, not a '
         'significant effect within the motif class.'),
    ]),
    '9B': ('Substrate rate in each stability x motif cell', [
        ('What is plotted',
         'The same classification extended to all three motif groups. Each pair of bars is one motif '
         'group split by baseline stability; bar height is the percentage of peptides in that cell '
         'that are UBR3 substrates. Counts beneath each bar are substrates / peptides in the cell.'),
        ('What it shows',
         '[P/G]-[E/D]: 13.41% (11/82) unstable, 5.15% (5/97) stable. [P/G]-other: 0.71% (9/1,261) '
         'and 0.15% (1/660). non-P/G: 0.30% (20/6,769) and 0.10% (8/7,645).'),
        ('The two effects are separable',
         'Reading ACROSS groups, the motif multiplies the substrate rate by roughly 45-50x within '
         'either stability class - that is the large, highly significant effect. Reading WITHIN a '
         'pair, being unstable multiplies it by about 2.5-3x - a real but much smaller effect.'),
        ('Which comparisons are powered',
         'The motif comparison is significant in both strata (p = 2.5 x 10^-14 unstable, '
         'p = 3.3 x 10^-7 stable). The stability comparison is significant library-wide '
         '(OR 2.97, p = 3 x 10^-4) but NOT within any single motif group - within [P/G]-[E/D] it is '
         'OR 2.85, p = 0.067. Treat the within-pair differences as directional only.'),
    ]),
    '9C': ('Every motif-bearing substrate gains stability', [
        ('What is plotted',
         'One arrow per motif-bearing substrate, running from its control PSI (dot) to its UBR3-KO '
         'PSI (arrowhead). Genes are sorted by starting PSI. Purple = starts unstable (below the '
         'dashed line at PSI 2.6), green = starts stable.'),
        ('What it shows',
         'All 16 arrows point right - every motif-bearing substrate gains stability when UBR3 is '
         'lost, with no exceptions. Nine of the 11 purple arrows cross the dashed line; the two that '
         'do not (ERMAP, UVRAG) start so low that even a substantial gain leaves them below it.'),
        ('The five green ones',
         'MYLK, LPIN1, ATF7, THEG and SLC7A9 start ABOVE the line, so they are already relatively '
         'stable proteins - yet they still gain 0.57 to 0.93 PSI. They are genuinely UBR3-regulated; '
         'they simply had a higher starting point.'),
        ('Why this panel matters',
         'It shows the stability classification does not sort substrates from non-substrates. It '
         'sorts substrates by where they begin. UBR3 dependence is present in both classes.'),
    ]),
    '9D': ('Where the 179 motif-bearing peptides end up', [
        ('What is plotted',
         'The classification as a simple count cascade. The top grey bar is all 179 motif-bearing '
         'peptides; beneath it the split into 82 unstable and 97 stable, each followed by the '
         'orange count of UBR3 substrates within that class.'),
        ('The numbers to quote',
         '179 peptides carry [P/G]-[E/D]. 82 are unstable at baseline, of which 11 are UBR3 '
         'substrates. 97 are stable at baseline, of which 5 are substrates. 16 substrates in total, '
         '163 non-substrates.'),
        ('How to read it in the manuscript',
         'Roughly two thirds of the motif-bearing substrates (11 of 16) come from the unstable half '
         'of the library, but the stable half still contributes 5. Restricting attention to '
         'low-PSI peptides would therefore have missed nearly a third of them.'),
        ('Watch out for',
         'The two stability bars are siblings that sum to 179; each orange bar is a '
         'subset of the stability bar directly above it, not of the bar before it.'),
    ]),
}

CLOSING = [
    ('The reviewer asked for',
     'A clean table of the candidate UBR3 substrates with functional annotation, presented '
     'alongside the peptides that carry the [P/G]-[E/D] feature but show no stabilisation - and '
     'evidence that only a subset of Pro-/Gly-initiated peptides is UBR3-regulated.'),
    ('The answer, in one line',
     'The motif enriches for UBR3 substrates 27-fold (q = 3 x 10^-18, odds ratio 50 versus '
     'non-P/G), yet 163 of the 179 library peptides carrying it (91%) are not stabilised when UBR3 '
     'is lost. Necessary, but plainly not sufficient.'),
    ('What carries the specificity',
     'Not Pro/Gly on its own. Bare Pro/Gly raises the substrate rate only from 0.19% to 0.52%; '
     'adding the acidic residue at position 3 raises it to 8.94%. Asp at position 3 is the single '
     'strongest cell in the entire dataset (q = 0.002).'),
    ('How localised the signal is',
     'Of 480 position x residue cells tested across the whole 24-mer, four survive FDR correction: '
     'Pro and Gly at position 2, Asp at position 3, and Arg at position 7 (borderline, q = 0.04). '
     'Whole-peptide composition is indistinguishable from the library.'),
    ('The control that closes the loop',
     'The motif is evenly distributed across baseline-stability strata (1.01% vs 1.15%) yet gives '
     'the same odds ratio within each (52 and 52), and holds at OR 52.8 in a cutoff-free logistic '
     'regression adjusting for control PSI. The enrichment is not a by-product of Pro/Gly N-termini '
     'being intrinsically unstable.'),
    ('The nuance worth keeping',
     'Motif-bearing peptides that never made the substrate list still cross the stability line at '
     'twice background. Substrate status is the top of a continuum rather than a separate category, '
     'which is a more accurate framing than a hard necessary/sufficient dichotomy.'),
]
