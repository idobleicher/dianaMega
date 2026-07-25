# Is the UBR3 P–[D/E] motif enriched after the 2A ribosomal skip?

Enrichment analysis of the residue at **position 2** of 2A downstream products — the
residue UBR3 reads — against three backgrounds built from viral proteins.

## Reproduce

```
python enrichment.py       # -> enrichment_results.csv, _composition.csv, _summary.json
python plot_enrichment.py  # -> enrichment_figure.png/pdf/svg
python build_pptx.py       # -> 2A_UBR3_enrichment.pptx  (9 slides)
```

`viral_swissprot.fasta.gz` (17,517 reviewed viral proteins, UniProt taxonomy 10239) is
committed so the analysis runs offline.

## Setup

The 2A skip occurs between the terminal G and P of `…D(V/I)ExNPG↓P`, so every downstream
product starts with that proline — **position 1 is P by construction**. Position 2 is the
only variable position, and it comes from the protein downstream of the skip site, not
from the 2A peptide. The question is whether D/E appears there more often than chance.

**Foreground:** `2a_combined_UBR3/combined_all_instances.csv`, 35,573 instances with a
resolved position 2.

**Backgrounds** (all from the 17,517 reviewed viral proteins, 7,989,533 residues):

| # | Background | D/E |
|---|---|---|
| 1 | All amino acids, all positions | 11.05% |
| 2 | Position 2 of viral proteins (after initiator Met; 17,136 proteins) | 21.24% |
| 3 | Residue immediately after any proline (P+1) | 11.15% |

Background 3 is the best-matched control — position 2 sits directly after a proline, and
proline constrains what follows it.

## The counting unit decides the answer — read this first

84% of instances come from **UniParc**, which archives every sequence variant separately.
The 35,573 instances collapse to **2,693 distinct downstream contexts**, a 12.3-fold
duplication. Counting raw instances counts the same context over and over.

| Counting unit | n | D/E at position 2 |
|---|---|---|
| All instances (raw) | 35,573 | 2.65% |
| Unique source accession | 25,326 | 3.46% |
| Unique 2A core peptide | 1,430 | 11.05% |
| **Unique downstream context (used here)** | **2,693** | **12.14%** |

Position 2 is a property of the downstream protein, so each distinct downstream product
should count once. That is the primary unit for every result below.

## Result

**Observed: 12.14% D/E at position 2** (D 9.43%, E 2.71%; 95% CI 10.93–13.44).

| Comparison | Expected | Fold | Fisher p | Verdict |
|---|---|---|---|---|
| vs all amino acids, all positions | 11.05% | 1.10× | 0.075 | no difference |
| vs residue after any proline | 11.15% | 1.09× | 0.104 | no difference |
| vs position 2 of viral proteins | 21.24% | 0.57× | 1.1e-30 | depleted |

The viral-only subset behaves identically (12.45%, n = 522).

**D/E follows the skip proline at essentially the rate chance predicts.** There is no
selection for or against the UBR3 motif at 2A skip sites.

The depletion against viral N-termini is real but means less than it appears. Position 2
of real N-termini is acidic-rich (21%) because Met-aminopeptidase removes the initiator
methionine only when residue 2 has a small side chain (A, C, G, P, S, T, V); D/E at
position 2 blocks cleavage and therefore accumulates there. That comparison measures
N-terminal biology, not 2A biology.

## This corrects an earlier conclusion

`2a_combined_UBR3/README.md` reported that D/E at position 2 occurs in "only 2.4% of all
resolved 2A sequences — so it is genuinely uncommon, not a default feature." That figure
came from raw instance counts dominated by UniParc duplication. On distinct downstream
contexts it is 12.1%, five times higher and statistically at background. **The motif is
not rare; it is ordinary.**

This does not change the practical output of the screen: the 942 P–[D/E] instances /
280 distinct peptides are still the candidate UBR3 substrate list. It changes the claim
you can make about them — they are not an unusual class, they are what chance produces.

## Caveats

- 6.8% of instances still have no resolved downstream residue, almost all MGnify (MGYP).
  `2a_rnabioco_compare/resolve_mgnify.py` is written but needs the ~270 GB MGnify dump.
- The background is reviewed viral proteins; the 2A set also contains host and
  unclassified sequences. The viral-only subset gives the same answer.
- De-duplication is on the 20-aa downstream product, so homologous proteins with slightly
  different downstream sequence still count separately. Residual redundancy makes the
  12.1% estimate conservative, not inflated.
- Pooled composition treats related viruses as independent observations. A
  phylogeny-aware background would be stronger.

## Files

| file | contents |
|---|---|
| `2A_UBR3_enrichment.pptx` | 9-slide deck: question, backgrounds, counting problem, result, figure, interpretation, caveats, methods |
| `enrichment_figure.png/pdf/svg` | 3-panel figure (counting unit · full composition · fold change) |
| `enrichment_results.csv` | every foreground × background comparison with fold, CI, Fisher and chi-square p |
| `enrichment_composition.csv` | full 20-amino-acid composition for every foreground and background |
| `enrichment_summary.json` | headline numbers, redundancy stats, parameters |
| `viral_swissprot.fasta.gz` | the background protein set |
| `enrichment.py`, `plot_enrichment.py`, `build_pptx.py` | the pipeline |
