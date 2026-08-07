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
     'The cut is the antimode of the control-PSI density, not a round number. Figure 8 shows '
     'nothing depends on this choice.'),
    ('Amino-acid classes',
     'Acidic D,E | Basic K,R,H | Polar S,T,N,Q,C | Hydrophobic A,V,L,I,M | Aromatic F,W,Y | '
     'Special G,P.'),
    ('Statistics',
     'Two-sided Fisher exact throughout, with Benjamini-Hochberg FDR applied across the residues, '
     'classes or cells tested in each analysis. Group ΔPSI comparisons use the two-sided '
     'Mann-Whitney U test. Proportion error bars are Wilson 95% intervals.'),
]

FIGURES = {
    '1': ('Figure 1 | Pro/Gly N-termini and the [P/G]-[E/D] motif in the screen',
          'Establishes the basic observation: Pro and Gly at position 2 are the only residues '
          'enriched among UBR3 substrates, and adding an acidic residue at position 3 sharpens '
          'that enrichment enormously. Panels A-B characterise position 2; C-D show what the '
          'motif buys you.'),
    '2': ('Figure 2 | The motif is necessary but not sufficient',
          'The reviewer\'s central point, argued four ways. The [P/G]-[E/D] motif enriches for '
          'substrates 27-fold, yet 163 of the 179 peptides carrying it are not stabilised at all. '
          'The motif biases stability without determining it.'),
    '3': ('Figure 3 | Sequence logos of the substrate N-termini',
          'Where along the 24-mer the constraint actually lives. Panels A-B compare substrates '
          'against the library background on an identical scale; panel C converts the comparison '
          'into a statistical enrichment logo so only real signal is drawn.'),
    '4': ('Figure 4 | Chemical-class architecture across all 24 positions',
          'The same question at the level of chemical classes rather than individual residues. '
          'Panels A-B show composition position by position, C tests every class x position cell '
          'for significance, and D checks whole-peptide composition as a control.'),
    '5': ('Figure 5 | Position-by-residue heatmaps',
          'The full 24 x 20 picture. Panel A is raw substrate frequency; panel B is enrichment '
          'over the library with non-significant cells masked, which reduces the whole dataset to '
          'four surviving cells.'),
    '6': ('Figure 6 | Baseline stability: peptides above and below PSI 2.6',
          'Splits the library by how stable each peptide is BEFORE UBR3 is removed. A degron '
          'substrate should start unstable and gain stability when its E3 ligase is lost - which '
          'is what the substrates do. Also exposes the assay\'s ceiling effect.'),
    '7': ('Figure 7 | The motif effect is not a by-product of baseline instability',
          'The key control. Pro/Gly N-termini are themselves destabilising, so the motif could in '
          'principle be enriched for that reason alone. Stratifying by control PSI shows it is '
          'not: the motif is evenly distributed across strata yet enriches for substrates within '
          'each one, with essentially the same odds ratio.'),
    '8': ('Figure 8 | Why control PSI, not ΔPSI - and why the cut does not matter',
          'Answers two methodological objections. First, ΔPSI cannot define the stability '
          'strata because the substrates were selected on it. Second, the PSI 2.6 cut is the '
          'empirical antimode, and every conclusion survives both a cutoff sweep and a fully '
          'cutoff-free model.'),
}

PANELS = {
    # ---------------------------------------------------------------- FIG 1
    '1A': ('Position 2 is Pro/Gly in half the substrates', [
        ('What is plotted',
         'Each of the 20 amino acids as a pair of horizontal bars: grey = its frequency at position '
         '2 across all 16,514 library peptides, blue = its frequency among the 54 UBR3 substrates. '
         'Rows are ordered by library frequency, so the grey bars descend. The two shaded rows are '
         'Pro and Gly.'),
        ('Why position 2',
         'Position 1 is Met in every peptide. When position 2 carries a small side chain, MetAP '
         'clips the initiator Met and position 2 becomes the residue actually exposed at the '
         'N-terminus - which is what an N-degron pathway would read.'),
        ('What it shows',
         'Pro and Gly together account for 12.7% of the library but 48% of the substrates (26% Gly '
         '+ 22% Pro). Alanine shows the opposite pattern: it is the single most common residue at '
         'position 2 in the library (21%) yet appears in only 6% of substrates.'),
        ('Watch out for',
         'These are raw frequencies, not a significance test. Because the library base rates differ '
         'so much between residues, judge enrichment from panel B, not from bar length here.'),
    ]),
    '1B': ('Pro and Gly are the only enriched residues at position 2', [
        ('What is plotted',
         'log2 fold enrichment of each residue at position 2 in substrates versus the library, '
         'sorted from most depleted to most enriched. Blue bars point right (enriched), red bars '
         'point left (depleted). The label on each bar gives the number of substrates carrying that '
         'residue, with asterisks when the FDR-corrected q value clears 0.05.'),
        ('How it was computed',
         'A two-sided Fisher exact test per residue - substrates carrying it versus not, against the '
         'library carrying it versus not - followed by Benjamini-Hochberg correction across all 20 '
         'tests.'),
        ('What it shows',
         'Only three residues survive correction. Pro is enriched 4.3-fold (q = 3.4 x 10^-4) and Gly '
         '3.4-fold (q = 3.7 x 10^-4). Alanine is significantly DEPLETED at 0.27-fold (q = 0.025). '
         'Everything else, including Thr and Gln which look elevated, is non-significant once the '
         '20 tests are accounted for.'),
        ('Watch out for',
         'Bars with n = 1 or 2 substrates carry almost no information; their apparent depletion is '
         'driven by small numbers and none of them reach significance.'),
    ]),
    '1C': ('From 16,514 peptides down to 16 motif-bearing substrates', [
        ('What is plotted',
         'A nested funnel on a log10 x-axis. Each bar is a strict subset of the bar above it, so '
         'the counts only ever shrink. The italic label between bars gives the percentage retained '
         'at that step.'),
        ('What it shows',
         'The full library is 16,514 peptides. Requiring Pro or Gly at position 2 keeps 2,100 '
         '(12.7%). Additionally requiring Asp or Glu at position 3 collapses that to 179 (8.5% of '
         'the P/G set, about 1% of the library). Of those 179, 16 are UBR3 substrates.'),
        ('Why it matters',
         'It shows how specific the motif is as a filter: roughly one library peptide in a hundred '
         'carries it. That rarity is what makes the enrichment in panel D so large.'),
        ('Watch out for',
         'The log scale compresses the differences - the top bar is 1,000x the bottom one. Read the '
         'printed numbers, not the bar lengths, and note that the final 8.9% is the substrate hit '
         'rate, not another sequence filter.'),
    ]),
    '1D': ('The acidic residue at position 3 is what matters', [
        ('What is plotted',
         'The percentage of peptides in each motif class that are UBR3 substrates. Error bars are '
         'Wilson 95% confidence intervals; the counts below each bar give substrates over group '
         'size.'),
        ('What it shows',
         'The gradient is steep and monotonic: 8.94% for [P/G]-[E/D] (16 of 179), 0.52% for '
         '[P/G]-other (10 of 1,921), and 0.19% for non-P/G (28 of 14,414). The bracket reports the '
         'Fisher test comparing the two ends: odds ratio 50.4, p = 4.1 x 10^-20.'),
        ('The key inference',
         'Pro/Gly alone raises the substrate rate only from 0.19% to 0.52%. Adding the acidic '
         'residue at position 3 raises it to 8.94% - a further 17-fold jump. The acidic residue, '
         'not the Pro/Gly by itself, is what carries UBR3 dependence.'),
        ('Watch out for',
         'The [P/G]-[E/D] confidence interval is wide (roughly 5.6% to 14%) because it rests on only '
         '16 substrates. The ordering of the three groups is solid; the exact 8.94% is not precise.'),
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
         'for multiple testing. Treat them as noise; the FDR-corrected version is Figure 5B.'),
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
         'scattered to position 24 that are pure sampling noise. Masking is what makes the figure '
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
         'stability. ΔPSI is the variable the substrates were selected on and cannot be used '
         'to stratify them - see Figure 8A.'),
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
         'stratified analysis in Figure 7 depends on.'),
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
         'figure at the data-driven cut is the more honest number.'),
        ('What is NOT cutoff-sensitive',
         'The motif results in Figure 7 are stable across every cut from 2.0 to 3.4 (Figure 8C). '
         'Only this baseline-stability contrast moves.'),
    ]),

    # ---------------------------------------------------------------- FIG 7
    '7A': ('The motif effect survives stratification', [
        ('What is plotted',
         'Substrate rate for each motif class, computed separately WITHIN each baseline-stability '
         'stratum. Purple bars are peptides that start unstable, green bars those that start stable. '
         'Counts are printed under each bar.'),
        ('What it shows',
         'Among unstable peptides: 13.41% for [P/G]-[E/D] versus 0.295% for non-P/G, Fisher OR 52.3, '
         'p = 2.5 x 10^-14. Among stable peptides: 5.15% versus 0.105%, OR 51.9, p = 3.3 x 10^-7.'),
        ('The point of the panel',
         'The two odds ratios are 52 and 52 - essentially identical. That is what "no interaction '
         'with baseline stability" looks like. The motif multiplies a peptide\'s chance of being a '
         'substrate by the same factor regardless of how stable it started.'),
        ('Watch out for',
         'The absolute rates differ between strata (13.4% vs 5.2%) even though the odds ratios do '
         'not. That is the ceiling effect again: stabilisation is simply harder to detect at high '
         'baseline PSI.'),
    ]),
    '7B': ('[P/G]-[E/D] is not over-represented among unstable peptides', [
        ('What is plotted',
         'How common each sequence feature is within each stratum. Left pair: the [P/G]-[E/D] motif. '
         'Right pair: bare Pro or Gly at position 2, with no requirement at position 3.'),
        ('What it shows',
         'The motif is carried by 1.01% of unstable and 1.15% of stable peptides - evenly split. '
         'Bare Pro/Gly is very unevenly split: 16.56% of unstable versus 9.01% of stable peptides.'),
        ('Why this is the crucial control',
         'If [P/G]-[E/D] were simply a marker of intrinsically unstable peptides, it would be '
         'concentrated in the purple bars. It is not. So its association with substrate status '
         'cannot be operating through baseline stability - the confound is empirically excluded, '
         'not just argued away.'),
        ('The mechanistic reading',
         'Pro/Gly N-termini ARE destabilising - that is the right-hand pair. But that instability is '
         'not UBR3-dependent. The acidic residue at position 3 is what converts a generically '
         'unstable N-terminus into a UBR3 substrate.'),
    ]),
    '7C': ('Crossing is graded, not all-or-nothing', [
        ('What is plotted',
         'Among peptides that start below PSI 2.6, the percentage that rise above it when UBR3 is '
         'lost. Three groups: the substrates, motif-bearing peptides that are NOT substrates, and '
         'all library peptides as background.'),
        ('What it shows',
         'A three-tier gradient: 77.5% of substrates cross (31 of 40), 21.1% of motif-bearing '
         'non-substrates (15 of 71), and 8.9% of library peptides overall (726 of 8,112).'),
        ('The finding worth flagging',
         'The 163 motif-bearing "non-responders" are not inert. They cross at more than twice the '
         'background rate. Carrying the motif produces a real effect that simply falls below the '
         'threshold for being called a substrate.'),
        ('How to phrase it',
         'This slightly softens "necessary but not sufficient" in a useful direction: substrate '
         'status is the top of a continuum rather than a distinct category. The motif confers a '
         'graded increase in UBR3 sensitivity.'),
    ]),
    '7D': ('Within unstable peptides, the motif still shifts ΔPSI', [
        ('What is plotted',
         'Box plots of mean ΔPSI by motif class, restricted to peptides that start below PSI '
         '2.6. Boxes give median and quartiles, whiskers the 1.5 IQR range. Individual points are '
         'overlaid where the group is small enough to show them.'),
        ('What it shows',
         'The [P/G]-[E/D] box sits visibly higher than the other two, with a longer upper whisker. '
         'Mann-Whitney against non-P/G gives p = 1.7 x 10^-6.'),
        ('Why restrict to the unstable stratum',
         'It removes baseline stability as an explanation by construction. Every peptide in this '
         'panel started in the same stability regime, so the remaining difference between boxes has '
         'to come from the sequence motif itself.'),
        ('Watch out for',
         'The boxes overlap substantially. As in Figure 2B, the effect is a distributional shift, '
         'not a separation - most motif-bearing peptides still behave like the background.'),
    ]),

    # ---------------------------------------------------------------- FIG 8
    '8A': ('PSI and ΔPSI measure different things', [
        ('What is plotted',
         'Control PSI against mean ΔPSI for all 16,514 peptides, with the 54 substrates in '
         'orange. The dashed vertical line is the stratification cut.'),
        ('What it shows',
         'The two axes are nearly independent: Spearman rho = -0.20, meaning baseline stability '
         'explains only about 3% of the variance in UBR3 dependence.'),
        ('Why ΔPSI cannot define the strata',
         'The 54 substrates were SELECTED on ΔPSI. Splitting the library on ΔPSI and then '
         'asking where the substrates fall is circular - at every threshold tested (0.1, 0.2, 0.3, '
         '0.5) all 54 land in the high group by construction. Control PSI is measured before UBR3 '
         'is perturbed and is therefore a legitimate stratifier.'),
        ('Watch out for',
         'The mild negative correlation is itself the ceiling effect: peptides starting at high PSI '
         'have less room to gain, which pushes their ΔPSI down.'),
    ]),
    '8B': ('The cut is data-driven, not a round number', [
        ('What is plotted',
         'A kernel density estimate of control PSI. Grey dotted lines mark the modes; red dots mark '
         'the density minima (antimodes). The dashed black line is the cut used throughout.'),
        ('What it shows',
         'Three robust modes at approximately 1.52, 2.25 and 3.49, with antimodes at 1.94 and 2.61. '
         'The cut at 2.6 sits on the second antimode - the natural boundary between the two unstable '
         'modes and the stable one.'),
        ('A modal structure that is NOT claimed',
         'At narrow kernel bandwidths a fourth mode appears near 2.9, which would make PSI 3 look '
         'like an antimode too. It disappears by bandwidth 0.20 and is therefore not reproducible. '
         'An earlier version of this analysis used PSI 3 on that basis; it was corrected.'),
        ('Why it matters less than it seems',
         'Panel C shows the conclusions are cutoff-independent anyway. The antimode is the '
         'principled choice, but nothing rests on it.'),
    ]),
    '8C': ('The conclusion does not depend on where the line goes', [
        ('What is plotted',
         'The Fisher odds ratio for [P/G]-[E/D] versus non-P/G within the unstable stratum, '
         'recomputed at every cutoff from 2.0 to 3.4. The y-axis is log-scaled. The orange point is '
         'the cut actually used.'),
        ('What it shows',
         'The odds ratio stays between 42 and 88 across the entire range, and the p value is at most '
         '6 x 10^-5 at every single cutoff. There is no threshold at which the motif result '
         'disappears or even weakens materially.'),
        ('Why include it',
         'It pre-empts the obvious methodological objection - that the cut was chosen to produce the '
         'result. Showing the full sweep is more convincing than defending any single choice.'),
        ('Watch out for',
         'The bounce at 2.2 (OR 88) is a small-sample artefact: at low cutoffs the unstable stratum '
         'contains few motif-bearing peptides, so the estimate is unstable. The trend across the '
         'range, not any individual point, is the message.'),
    ]),
    '8D': ('With no cutoff at all, the motif effect stands', [
        ('What is plotted',
         'A forest plot of adjusted odds ratios from a single logistic regression predicting '
         'substrate status from motif class plus control PSI as a CONTINUOUS covariate. Points are '
         'odds ratios, bars are 95% confidence intervals, and the dashed line at 1 is no effect. '
         'The x-axis is log-scaled.'),
        ('The model',
         'substrate ~ [P/G]-[E/D] + [P/G]-other + control PSI, fitted across all 16,514 peptides. No '
         'threshold appears anywhere in it, so the stratification question is bypassed entirely.'),
        ('What it shows',
         '[P/G]-[E/D] gives an adjusted OR of 52.8 (95% CI 27.9-99.9, p = 3.6 x 10^-34). '
         '[P/G]-other gives 2.33 (p = 0.023). Control PSI gives 0.61 per unit (p = 0.006) - more '
         'stable peptides are less likely to be substrates. Both effects are real and independent '
         'of each other.'),
        ('Recommended use',
         'Lead the manuscript with this analysis and present the stratified figures as the intuitive '
         'illustration. A cutoff-free model leaves nothing for a reviewer to attack on the '
         'threshold question.'),
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
