# ZZ- paper — ZYG11B / ZER1 N-terminome GPS screen

Started 2026-09-05. First task: read the Data S3 screen, and pull out the peptides most stabilised
on loss of ZYG11B, ZER1 and the double knockout — the equivalent of the `01_UBR3_substrates` list
in `../ubr3_PG_reviewer/`.

---

## The one-paragraph answer

**63 peptides are stabilised by a full stability bin in both the ZYG11B knockout and the
ZYG11B/ZER1 double knockout, and 55 of the 63 (87%) carry an N-terminal glycine** — 11× the 7.7%
glycine rate of the library, against an estimated false-discovery rate of 3%. That is the
headline list, `hits_high_confidence`. **ZER1 alone produces no such list:** once baseline
stability is controlled for, its glycine effect is not merely weak but absent (mean ΔPSI
difference −0.012, Mann-Whitney p = 0.45), and its top 127 peptides carry glycine at 8.7%, the
library rate. Everything ZYG11B-specific in this screen is carried by ZYG11B and by the double
knockout.

| List | n | % Gly | vs library | Est. FDR |
|---|---:|---:|---:|---:|
| **`hits_high_confidence`** — ZYG11B **and** double KO ≥ 1.0 | **63** | **87.3%** | **11.3×** | **3.2%** |
| `hits_double_ko` — double KO ≥ 1.0 | 194 | 62.4% | 8.1× | 1.0% |
| `hits_zyg11b_ko` — ZYG11B KO ≥ 1.0 | 214 | 46.7% | 6.0× | 9.8% |
| `hits_zer1_ko` — ZER1 KO ≥ 1.0 | 127 | 8.7% | 1.1× | 15.0% |

---

## Source data

`data/Data_S3_Nterminome_GPS_screen.xlsx` — "Data S3. N-terminome GPS screen data in different
genetic backgrounds", copied verbatim from the download. Read only through `zz_gps_data.py`.

| Sheet | Peptides | Control | Knockouts |
|---|---:|---|---|
| **Data S3A** | 20,720 | AAVS1 KO | UBR clones #1–3, **ZYG11B**, **ZER1** — one clone each |
| **Data S3B** | 18,099 | Wild-type | **ZYG11B/ZER1 double KO** |
| Data S3C | 23,356 | AAVS1 KO | NMT1/2 clones #1–3 — loaded, no hit list built |

17,901 transcripts appear in both S3A and S3B. Each row is one Ensembl transcript's N-terminal
24-mer; position 1 is Met in 100% of the library and every peptide is exactly 24 residues, both
checked on load. Read depth is already floored at 200 reads per peptide per sample.

The workbook's own ΔPSI columns are used as given. They reproduce PSI(KO) − PSI(control) to
1.5 × 10⁻¹⁵, i.e. exactly; the loader asserts it rather than assuming it.

## Definitions

| Term | Definition |
|---|---|
| **PSI** | Protein Stability Index: the read-weighted mean sorting bin, 1 (unstable) to 6 (stable). |
| **ΔPSI** | PSI(knockout) − PSI(control). **Positive = stabilised when the receptor is lost**, which is what makes a peptide that receptor's substrate. |
| **Position 2 / N-terminus** | Residue 2 of the 24-mer. Position 1 is Met throughout and is excised when position 2 is small (A C G P S T V), so position 2 is the residue actually exposed. 54.6% of the library qualifies. |
| **Gly N-degron** | Position 2 = Gly — the feature ZYG11B and ZER1 read. **7.7% of the library** (11.7% among peptides with headroom). |
| **Headroom** | Baseline PSI ≤ 4.0 in the experiment being read. |
| **Hit** | ΔPSI ≥ 1.0 — one full stability bin — with headroom. |

---

## Two things that decide every number here

**The ceiling.** A peptide already stable in the control cannot be stabilised further. Mean ΔPSI
falls monotonically with baseline PSI and turns negative above 4.5, so the left tail of every ΔPSI
distribution is heavier than the right for reasons that have nothing to do with any receptor.
Unfiltered, a mirror null returns a false-discovery rate **over 100%** for the single knockouts —
more peptides fall by 1.0 than rise by it. With the headroom filter the same estimate lands at 10%
for ZYG11B and 1% for the double KO. Nothing in this project is thresholded without it.

This is not a cosmetic filter, and it changes a conclusion. Across the whole library ZER1 looks
like it has a real glycine effect (+0.069 vs −0.016, p = 1.3 × 10⁻¹⁵). Glycine peptides are
unstable at baseline, and unstable peptides drift upward regardless of genotype; restrict to
peptides with headroom and **ZER1's glycine effect disappears entirely** (−0.012, p = 0.45) while
ZYG11B's and the double KO's survive untouched. The unfiltered ZER1 result is the confound, not a
finding.

| Genotype (headroom only) | Gly mean ΔPSI | Other | Difference | Mann-Whitney p |
|---|---:|---:|---:|---:|
| Double KO | +0.303 | +0.007 | **+0.296** | 4 × 10⁻⁶² |
| ZYG11B KO | +0.339 | +0.173 | **+0.166** | 2 × 10⁻³⁸ |
| ZER1 KO | +0.134 | +0.146 | −0.012 | 0.45 |
| UBR KO #1 | +0.065 | +0.128 | −0.063 | 8 × 10⁻⁹ |

UBR #1 is the specificity control and behaves correctly: it reads different degrons and its
glycine effect is slightly *negative*.

**The shared control.** Every ΔPSI in S3A is measured against the same AAVS1 sample, so noise in
that one control propagates into all five knockouts and correlates them whatever the biology does.
ZYG11B and UBR #1 correlate at **r = 0.46** despite reading different degrons — that number is
about the control, not about either receptor. ZYG11B and the double KO, which have different
controls, correlate at only **r = 0.077**. Agreement between S3A and S3B is therefore the only
cross-check in this workbook that is not inflated, and it is what the headline list is built on.

**There is no within-genotype replicate.** ZYG11B, ZER1 and the double KO are one clone each. Only
the UBR knockout has three. Every per-peptide ΔPSI here rests on a single measurement, which is
the second reason the headline list requires two experiments to agree.

---

## The hit lists

`make_hit_lists.py` → `ZZ_ZYG11B_ZER1_hit_lists.xlsx`, `data/hits_*.csv`, `data/hit_summary.csv`

A hit is **ΔPSI ≥ 1.0 with headroom**, one ruler for every genotype. `hits_high_confidence`
additionally requires both experiments to pass, ranked by the mean of the two ΔPSI.

False discovery is estimated with a **mirror null**: peptides at ΔPSI ≤ −1.0, counted over the
same headroom-filtered pool, standing in for how much of the +1.0 tail is noise. Strong
*destabilisation* on losing an E3 receptor has no mechanism behind it, so the negative tail is an
estimate of the noise. It is an approximation — the tails are not exactly symmetric even after
filtering — and it is quoted on every list rather than buried.

The glycine fraction is the built-in positive control. A correct ZYG11B list must be full of
glycine; a list sitting at the library's 7.7% is noise, whatever its ΔPSI.

**Top of the high-confidence list** (full 63 in the workbook, sheet `01_high_confidence`):

| Rank | Gene | Pos 2 | Mean ΔPSI | ZYG11B | Double KO | ZER1 | Baseline PSI |
|---:|---|:--:|---:|---:|---:|---:|---:|
| 1 | CCDC169 | P | 2.13 | 2.38 | 1.89 | 0.46 | 1.89 |
| 2 | C14orf178 | G | 1.88 | 2.03 | 1.74 | −0.05 | 3.03 |
| 3 | UCK2 | A | 1.85 | 1.76 | 1.94 | 0.31 | 3.34 |
| 4 | FXR2 | G | 1.79 | 1.52 | 2.05 | 0.28 | 2.60 |
| 5 | IFI27L1 | G | 1.76 | 1.84 | 1.69 | −0.14 | 2.29 |
| 6 | UBTD2 | G | 1.75 | 1.42 | 2.08 | 0.22 | 3.01 |
| 7 | HBG1 | G | 1.74 | 1.51 | 1.97 | 0.32 | 2.07 |
| 8 | ADARB2 | G | 1.74 | 1.17 | 2.31 | 0.34 | 2.10 |
| 9 | PPP2R2C | G | 1.73 | 1.57 | 1.89 | 0.53 | 2.52 |
| 10 | EPB42 | G | 1.70 | 1.43 | 1.97 | 0.14 | 2.64 |

Note how small every ZER1 column entry is beside the other two — that pattern holds down the
whole list.

**The 8 non-glycine members**, worth looking at precisely because they should not be there:
CCDC169 (Pro, and the single most stabilised peptide in the screen), KLHL32 (Pro), CTNNBL1 (Pro),
UCK2 (Ala), PPP2R2C (Arg), MAP3K4 (Arg), TRIM73 (Lys), KCTD20 (Asn). Three of the eight start
Met-Pro, where Met excision is inefficient — so the exposed N-terminus of those may not be what
position 2 predicts. They are either genuine non-glycine substrates or the list's few false
positives, and this screen cannot tell which.

PPP2R2C appears twice, as two transcripts with different N-termini (Gly and Arg). Every other gene
appears once; 62 genes across 63 peptides.

**Replication rate.** 63 of the 164 ZYG11B hits that were also measured in S3B pass in the double
KO too — 38%, against 2.2% of all headroom peptides passing there, a 17-fold enrichment. The
per-peptide correlation between the two experiments is near zero (r = 0.077) while the hits
replicate at 17× background: single-peptide ΔPSI is noisy across experiments, but the strong ones
reproduce.

---

## What this does and does not support

- **ZYG11B has a large, clean, glycine-specific substrate set in this screen.** Supported by two
  experiments with independent controls, an 11× glycine enrichment, and a 3% estimated FDR.
- **The double knockout is the single cleanest readout in the workbook** — 1% FDR at ΔPSI ≥ 1.0,
  62% glycine — and is the one to use where only one experiment is wanted.
- **ZER1 has no substrate list here.** Not "a weaker one": once baseline stability is controlled,
  no glycine effect and no glycine-enriched hits. Whether that means ZER1 is dispensable, is
  fully redundant with ZYG11B, or was poorly edited in this clone, the screen cannot say — but
  no ZER1-specific claim should be made from this data.
- **Nothing here is a within-genotype replicate.** Every ΔPSI is one clone, one measurement.
- **The mirror-null FDR is an estimate**, not a q value. There is no per-peptide p value in this
  design, and none is invented.

## Met-Pro peptides — is anything but the proline doing work?

`zz_pro_motif.py`, `make_pro_logo.py`, `make_pro_composition.py`
→ `figures/FigureZ1*`, `FigureZ2*`, `FigureZ3*`

Both sides of this comparison start Met-Pro, so the proline cannot separate them and positions
3–24 are the question. 36 stabilised Met-Pro peptides (ΔPSI ≥ 1.0 with headroom in any of ZYG11B
KO, ZER1 KO, double KO) against 424 that are not stabilised anywhere (ΔPSI < 0.5 everywhere
measured). The 662 Met-Pro peptides in between are in neither group — with 36 substrates, letting
peptides at ΔPSI 0.5–1.0 into the background would seed it with probable weak substrates.

**No, there is no second motif position.** Across the 440 position × residue cells, **25 reach
p < 0.05 against 22 expected by chance alone, and not one survives BH-FDR.** The two versions of
the logo — all three genotypes, and ZYG11B + double KO with the ZER1-only calls dropped — do not
even agree on which cells are tall: Gly at position 3 is the second-tallest glyph in the first and
absent from the second, Lys at 14 the reverse. That disagreement is what noise looks like.
Figures Z1 and Z2 are published as a negative result, because "there is nothing here" is a claim
a reader should be able to check.

**But a per-cell test is the wrong instrument for the question.** It spends its power on 440
hypotheses, each about one residue at one position, when the thing that differs is spread across
the whole sequence. Asking about **composition** instead — 27 tests, each pooling all 22 free
positions — finds it at once, and the effects are not subtle:

| Positions 3–24 | Stabilised | Not stabilised | q |
|---|---:|---:|---:|
| **D+E content** | **13.3%** | **7.3%** | 3 × 10⁻⁴ |
| Glu content | 8.2% | 4.2% | 0.003 |
| Hydropathy (Kyte–Doolittle) | −0.20 | +0.24 | 0.005 |
| Net charge | −0.75 | +0.50 | 0.006 |
| Asp content | 5.1% | 3.0% | 0.010 |
| Leu content | 10.9% | 16.7% | 0.016 |
| Arg content | 3.5% | 6.4% | 0.016 |

Seven of 27 tests survive FDR against 1.4 expected at p < 0.05 by chance. The acidic excess is
spread evenly along the peptide — 11.6%, 12.5%, 15.6%, 13.9% across positions 3–8, 9–14, 15–19 and
20–24 against roughly 7–8% in the controls — so it is bulk character, not a position.

**The control that makes it believable.** An acidic peptide could look stabilised merely by
sitting lower at baseline with more room to rise. It does not: among the controls, acidic content
correlates **positively** with baseline PSI (Spearman r = +0.22, p = 8 × 10⁻⁶), so the ceiling
works against this result rather than producing it. Compared band by band inside matched baseline
PSI the difference holds anyway — combined p = 0.0018. Panels D and E of Figure Z3 are that
control.

**It is not proline-specific, and that is the real conclusion.** The same composition separates
stabilised from unstabilised **Met-Gly** peptides, and far more strongly: 13.1% vs 6.2% acidic,
16 of 27 tests surviving FDR, combined p = 7 × 10⁻¹⁹ after baseline matching. Met-Gly peptides are
the canonical ZYG11B substrates. So an acidic, hydrophilic, negatively charged N-terminal region
is a property of what gets degraded in this screen generally — plausibly because such regions are
disordered and accessible — rather than a second determinant of how proline is read.

One more caution specific to these peptides: Met excision is inefficient at Met-Pro, so most of
this group probably keeps its initiator methionine and presents Met, not Pro, at the N-terminus.
Whatever is happening to them may not be an N-degron mechanism at all.

**Palette.** This project's figures are blue / purple / pink and deliberately do **not** share the
warm ramp of `../ubr3_PG_reviewer/`; they are a different paper. The three values live in
`zz_pro_motif.py` and every figure reads them from there.

Colour in the logos groups residues by **charge** — acidic D E in blue, basic K H R in pink,
everything uncharged in purple — rather than by the six chemical classes the sibling project uses.
That is a consequence of the palette, and it is worth stating why: six hues cannot be told apart
inside blue-purple-pink. The family is roughly a third of the hue wheel, and the best six-colour
set that fits in it still leaves a pair below the ΔE 15 floor for **normal** colour vision, which
no amount of legend or labelling excuses. The rule for that case is to cut the series count, so
three it is — and charge happens to be the axis this paper's result is about anyway.

The three (`#2E6FD0`, `#E0407A`, `#8E44AD`) pass chroma and the normal-vision floor at worst
ΔE 15.9. Colour-vision separation sits at ΔE 6.2 for blue against purple, which is in the band
that is legal only where the mark carries its own label — in a sequence logo it does, since the
glyph is the residue's name. Figure Z3's pair (`#D6408F` stabilised, `#2E6FD0` not) passes every
check outright. The −log₁₀ p ramp runs pale lavender → periwinkle → violet → purple → deep plum,
verified monotonic in luminance.

| Figure | What it shows |
|---|---|
| `FigureZ1_pro_foldchange_logo` | Fold change, stabilised vs not, enriched residues only |
| `FigureZ1b_..._significant` | The same, restricted to p < 0.05 |
| `FigureZ2_pro_significance_logo` | Height and colour = −log₁₀ p (Fisher) — **the one to read** |
| `FigureZ2b_..._strict` | The same without the ZER1-only calls; it disagrees with Z2 |
| `FigureZ3_composition_stabilised_vs_not` | Composition, with the baseline-matched control |

---

## Reproducing

```bash
python zz_gps_data.py          # loader self-test; writes data/gps_tidy.csv.gz
python make_hit_lists.py       # -> ZZ_ZYG11B_ZER1_hit_lists.xlsx, data/hits_*.csv
python zz_pro_motif.py         # Met-Pro group sizes and the top cells
python make_pro_logo.py        # -> figures/FigureZ1*, FigureZ2*, data/pro_motif_cells.csv
python make_pro_composition.py # -> figures/FigureZ3*, data/pro_composition_tests.csv
```

| File | Role |
|---|---|
| `zz_gps_data.py` | The single loader. Reads the workbook, checks its ΔPSI columns against its own PSI columns, joins S3A to S3B, and adds position 2, Met-excision and headroom. Everything else imports it. |
| `make_hit_lists.py` | Hit definition, the four lists, the mirror-null FDR, the glycine control, and the Excel workbook. |
| `zz_pro_motif.py` | The Met-Pro within-motif comparison: the two groups, the 440 position × residue cells, Fisher and BH-FDR. |
| `make_pro_logo.py` | Fold-change and significance logos for that comparison. |
| `make_pro_composition.py` | Composition and peptide-property tests, plus the baseline-matched control. |
| `data/gps_tidy.csv.gz` | Cached tidy table, one row per transcript. Delete it to force a re-read of the xlsx. |

## Still open

- Annotation. The UBR3 list carries UniProt function, localisation and GO process per protein;
  these lists carry gene symbol and transcript only. `../ubr3_PG_reviewer/fetch_annotation.py`
  does that job from UniProt and could be pointed at these 63.
- Figures. The Met-Pro question has Z1–Z3; there is still no overview panel of ΔPSI by
  position-2 residue for the three genotypes side by side, which is the whole screen in one image.
- Data S3C (NMT1/2) is loaded but unused. N-myristoyltransferase competes for the same
  N-terminal glycine, so it is the natural next comparison.
