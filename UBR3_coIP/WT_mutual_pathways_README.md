# Mutual pathways among the UBR3-WT best hits

**Question.** Do the best hits enriched in wild-type UBR3 share pathways with each
other, or does every enriched term come from a different protein?

**Answer.** They do share. Four mutual modules account for the specific enrichment,
three of them robust. But the convergence is partial: half of the 50 best hits sit in
no specific term at all, and a large block of the Reactome signalling terms in
`GO_WT_besthits.csv` turns out to be carried by the histones rather than by a real
multi-protein module.

## How to reproduce

```
python wt_mutual_pathways.py      # enrichment + modules + cross-checks -> CSV/JSON
python plot_wt_mutual_pathways.py # -> WT_mutual_pathways.png/pdf/svg
```

`go_analysis_separate.py` ran g:Profiler with `no_evidences=True`, which discards the
genes supporting each term — so the existing tables cannot say which hits share a
pathway. These scripts re-run the identical query with `no_evidences=False`
(reproducing the same 59 WT terms) and keep the membership.

## Method

1. **Best hits** — the same definition as `go_analysis_separate.py`: the 50 significant
   WT interactors (`hit_UBR3wt`, q<0.05 vs 293T), ranked by `WT_293T_log2diff`.
2. **Enrichment** — g:Profiler g:GOSt, hsapiens, GO:BP/CC/MF + KEGG + REAC + CORUM,
   g:SCS correction at 0.05, whole annotated genome background.
3. **Specific terms** — analysis of shared pathways uses only terms of size ≤ 500
   (29 of 59). The broad terms ("nucleus", "protein binding", "regulation of
   macromolecule metabolic process") contain nearly every hit and bridge everything
   into one blob; they carry no information about who shares what.
4. **Modules** — terms are linked when their supporting hits coincide (overlap
   coefficient ≥ 0.5) and clustered by greedy modularity.
5. **Histone control** — 7 histone genes are among the WT hits, and histones are
   annotated to a very large number of Reactome "gene expression" pathways. The
   enrichment is re-run without them; a term that disappears was carried by the
   histones, not by a genuine module.
6. **STRING cross-check** — modules are also scored against `STRING_edges.csv` /
   `STRING_clusters.csv`, which is independent of the annotation databases.

## Result

| Module | Theme | Hits | Best p | Survives histone removal | STRING edges / pairs |
|---|---|---|---|---|---|
| M1 | mRNA silencing & RNP granules | 14 | 1.2e-05 | 8/15 terms | 17 / 91 |
| M2 | ATP-dependent chromatin remodeling | 10 | 5.5e-03 | 6/9 terms | 16 / 45 |
| M3 | Nucleosome / histone core | 6 | 1.5e-02 | 0/4 terms | 5 / 15 |
| M4 | BAF–LAP2β nuclear-envelope complex | 2 | 5.0e-02 | 1/1 term | 1 / 1 |

- **M1** is the strongest and cleanest: DCP1A, MCRIP2, LARP1B, TNRC6B, TNRC6C, GIGYF2,
  UBAP2, UBAP2L in cytoplasmic RNP granule, P-body, stress granule, regulation of mRNA
  stability, post-transcriptional silencing by small RNAs. It gets *more* significant
  without the histones (p 1.2e-05 → 5.7e-06).
- **M2** is anchored by SMARCA5 + SMARCA1, which together form the CERF and NURF
  complexes (both terms are 2-gene CORUM/GO complexes, so their significance rests
  entirely on that pair), extended by CHD7 and, through "catalytic activity acting on
  DNA", to BLM, TET1 and POLE.
- **M3 is the histone block itself** — nucleosome, structural constituent of chromatin.
  Its "survives = 0" is tautological (removing the histones removes the module), but the
  KEGG *Systemic lupus erythematosus* and *beta-catenin:TCF* terms in it are histone
  annotation spillover, not chromatin biology.
- **M4** is BANF1 + TMPO (LAP2β), a 2-protein CORUM complex that is also a direct STRING
  edge and a shared STRING cluster.

**Histone caveat.** 21 of the 59 WT terms vanish when the 7 histone genes are dropped,
including *Estrogen-dependent gene expression*, *ESR-mediated signaling*, *Signaling by
Nuclear Receptors*, *Signaling by WNT*, both *Pre-NOTCH* terms, *RUNX1 regulates…* and
*DNA Repair*. These should not be reported as UBR3-WT pathways. They are marked in
`WT_mutual_terms.csv` (`survives_without_histones`) and drawn as open circles / italic
labels in the figure.

**Coverage caveat.** 25 of the 50 best hits — including UBR3 itself and the 3rd- and
8th-ranked hits TDRD3 and COL14A1 — fall in no specific term, so no pathway statement
covers them. 22 hits share at least one specific term with another hit.

**Genotype specificity.** Only 5 of the 59 WT terms are also enriched in the RING
mutant, and all 5 are generic nuclear compartments (nucleoplasm, nuclear lumen,
organelle lumen ×3). None of the four modules is enriched in RING. Consistently, 24 of
the 25 module members are classed *WT-preferring (RING-dependent binding)* in
`WTvsRING_differential.csv`; the single exception is POLE.

## Outputs

| File | Contents |
|---|---|
| `WT_mutual_terms.csv` | 59 WT terms + supporting genes + module + histone robustness |
| `WT_mutual_modules.csv` | module → core genes, all genes, member terms, STRING cross-check |
| `WT_mutual_genes.csv` | per hit: #terms, modules, mutual partners, rank, log2diff, role |
| `WT_mutual_gene_matrix.csv` | gene × term binary membership matrix |
| `WT_mutual_vs_RING.csv` | each WT term: mutual with RING or WT-exclusive |
| `WT_mutual_summary.json` | counts + all parameters |
| `WT_mutual_pathways.png/pdf/svg` | 2-panel figure (membership matrix; pairwise sharing) |
