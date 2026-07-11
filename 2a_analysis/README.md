# 2A downstream products — UBR3 P–[D/E] motif screen

**Source paper:** Rao et al. (2025) *Cell Reports* — "Systematic identification and
characterization of eukaryotic and viral 2A peptide-bond-skipping sequences."
PMID 40536874 / PMC12338488. DOI 10.1016/j.celrep.2025.115822.
Sequences taken from **Supplementary Table S2** (`NIHMS2099490-supplement-3.xlsx`),
the paper's full HMMER-identified 2A set (36,008 entries).

## Biology / definitions
- 2A peptides end in the conserved `...PG↓P` motif. Ribosomal skipping occurs between
  the terminal **G** and **P**.
- The **downstream product therefore begins with that proline (position 1 = P) — always.**
- **UBR3 motif screened:** downstream position 1 = **P**, position 2 = **D or E**
  (i.e. downstream products of the form `PD...` or `PE...`).

## Method
1. Parsed all 36,008 rows of Table S2 → source database + accession + coordinates.
2. Fetched the **full source protein** for each accession:
   - UniProtKB (2,531 unique acc) via `rest.uniprot.org/uniprotkb/accessions`
   - UniParc (20,935 unique acc) via `rest.uniprot.org/uniparc/stream`
   - MGnify (`MGYP…`, 2,175 unique acc) — **not retrievable** (no working per-ID API;
     ESMAtlas is a different MGnify release). These 2,240 rows are unresolved.
3. Located the skip proline in each source via the `PG|P` motif (not by trusting the
   trimmed table sequence, which for ~588 entries stops one residue short at `...NPG`).
4. Downstream product = 20 aa starting **at** the skip proline. Position 2 = next residue.

## Coverage
| | rows |
|---|---|
| Total 2A entries | 36,008 |
| Downstream resolved (position 2 known) | 31,644 |
| 2A at protein C-terminus (no residue after skip-P) | 498 |
| Unresolved (MGnify + obsolete accessions + no PGP found) | 3,866 |

## Result — UBR3 P–[D/E] motif
- **893** entries have downstream `P` + `D/E` at position 2.
- **327 unique** downstream sequences (**254** `PD…`, **73** `PE…`).
- By origin (unique): Host 90, Virus 52, Unclassified 159, mixed 26, class B 1.

## Files
| file | contents |
|---|---|
| `downstream_all.csv` | every row: 2A seq, source id, resolved downstream product (20 aa), position-2 residue, status |
| `downstream_all.fasta` | all resolved downstream products (FASTA) |
| `motif_PDE_hits.csv` / `.fasta` | the 893 P–[D/E] hits, sorted by amino-acid sequence |
| `motif_PDE_unique.csv` / `.fasta` | the 327 **unique** P–[D/E] downstream sequences (with occurrence counts), sorted |
| `parsed_rows.json`, `fasta_cache.json` | intermediate data (re-runnable) |
| `fetch_sources.py`, `reconstruct.py` | the pipeline |

## Caveats
- Downstream products are truncated to 20 aa; re-run `reconstruct.py` with a larger
  `DOWN_LEN` for longer context.
- 2,240 MGnify (`MGYP`) rows could not be resolved (6.2%); if needed these require the
  MGnify protein-database FTP dump.
- "position 2 = D/E" reflects the *natural* residue following the skip-P in each source
  protein. Classic engineered peptides (P2A/T2A, Table S1) have no fixed natural
  downstream and are not included here.
