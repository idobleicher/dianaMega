# rnabioco pipeline dataset vs. our prior Table S2 UBR3 screen

Compares the 2A hits shipped by **rnabioco/2a-peptide-search** (their HMMER
`combined-class-1/2.sto.gz` database searches) against our earlier
**Rao et al. Table S2** UBR3 P–[D/E] downstream screen (`../2a_analysis`).

## Datasets

| | rnabioco pipeline | our Table S2 screen |
|---|---|---|
| Total hits | **7,291** | 36,008 |
| Unique source accessions | 6,119 | 25,641 |
| Unique 2A peptide sequences | 1,824 | 1,989 |
| UniProtKB | 3,208 | 2,644 |
| MGnify (MGYP) | 2,383 | 2,240 |
| IMG/VR (viral metagenome) | **1,700** | 0 |
| UniParc | **0** | 31,124 |

The two sets are **complementary, not redundant**: Table S2 is UniParc-dominated
(31k, absent from rnabioco), while rnabioco adds 1,700 IMG/VR viral-metagenome
hits that our prior screen never saw.

## Overlap (identity)

- **656** hits shared by exact accession **and** coordinates (same subseq).
- **3,984** source accessions shared.
- **355** unique 2A peptide sequences shared (of 1,824 rnabioco / 1,989 ours) → only ~18%.
- **2,135** rnabioco accessions and **1,469** unique 2A sequences are **not** in our
  S2 set — overwhelmingly IMG/VR + MGnify. Genuinely new candidates.
- **21,657** S2 accessions are absent from rnabioco (the UniParc bulk).

## Source databases — what we could fetch

The rnabioco pipeline searched five databases. To resolve downstream position-2 for
hits that truncate at the skip-proline, the full source protein must be re-fetched:

| Database | In rnabioco set | Fetched | Notes |
|---|---|---|---|
| UniProt (curated) + Reference Proteomes | 3,208 | **yes** | REST API; 526 unresolved accs → 312 live, 214 obsolete |
| UniParc | 0 | n/a | not present in this set |
| MGnify (MGYP) | 2,383 | **not yet** | no per-ID API; needs the ~270 GB FASTA dump (`resolve_mgnify.py` ready) |
| IMG/VR | 1,700 | **yes** | JGI login; resolved from `IMG_VR_2022-12-19_7` full-proteins FASTA (31 GB) |

- UniProt: fetched 526 unresolved accessions → resolved **375** hits.
- IMG/VR: streamed the v7 `IMGVR_all_proteins.faa.gz` → **all 895** unresolved hits matched
  exactly by protein id (confirming v7 is the release the search used); 892 resolved,
  3 had a subseq mismatch.

## UBR3 P–[D/E] downstream angle

Downstream position-2 (the residue after the skip-proline, `…PG↓P·X`) comes from the
alignment column past `NPGP`, our prior S2 table, or a freshly-fetched full source protein.

- pos2 resolved for **5,991 / 7,291 (82%)** — 4,680 from the alignment column,
  162 from S2, 375 from fetched UniProt, 892 from fetched IMG/VR.
- **1,300 (18%) still unresolved** — MGnify 969 (needs the DB dump), obsolete/C-terminal
  UniProt 328, IMG/VR 3.
- Of resolved hits, **231 are P–[D/E]** (170 class-1, 61 class-2; **88 unique seqs**;
  by source: UniProt 120, MGnify 58, IMG/VR 53).
- **186** of those are **not** in our prior 893-hit set → new UBR3 candidates
  (**66 unique sequences**; UniProt 87, IMG/VR 52, MGnify 47). The 52 IMG/VR ones are
  viral-metagenome downstream products our Table S2 screen never had access to.

For reference our prior screen found **893** P–[D/E] hits (327 unique downstream
products / 209 unique 2A peptides) out of 31,644 resolved.

### Caveat
Position-2 is dominated by **F (33% of resolved)**, then S/T — an artifact of
redundant metagenomic hits (related viral genomes sharing the same downstream
protein). Judge by **unique-sequence** counts, not raw hit counts.

## How to resolve the remaining MGnify + IMG/VR hits

Both resolver scripts are written and tested; they stream the DB files (never load them
into RAM) and update `rnabioco_hits_annotated.csv` in place. Each run only targets rows
still lacking a pos2, so you can process **one split at a time and delete it** — you
never need all 270 GB on disk at once.

### MGnify (969 hits) — no login needed
```bash
# base url from the pipeline's own config:
URL=https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/current_release
mkdir -p mgy && cd mgy
# option A – all 25 splits (~270 GB), then one pass:
for i in $(seq 1 25); do wget -c "$URL/mgy_proteins_$i.fa.gz"; done
cd .. && python resolve_mgnify.py . mgy

# option B – disk-friendly, one split at a time:
for i in $(seq 1 25); do
  wget -c "$URL/mgy_proteins_$i.fa.gz" -O mgy/mgy_proteins_$i.fa.gz
  python resolve_mgnify.py . mgy         # resolves what this split contains, saves progress
  rm mgy/mgy_proteins_$i.fa.gz
done
```

### IMG/VR (895 hits) — free JGI account required
1. Register / log in at **https://img.jgi.doe.gov/vr/** (DOE JGI single sign-on).
2. Download the IMG/VR v4 protein FASTA (e.g. `IMGVR_all_proteins-high_confidence.faa.gz`).
3. Run:
   ```bash
   python resolve_imgvr.py . /path/to/IMGVR_all_proteins-high_confidence.faa.gz
   ```
   The script prints how many hits it matched. IMG/VR header formats vary by release —
   if the match count is low, tell me the header format (`zcat file | head -1`) and I'll
   adjust the key matching in `keys_for()`.

After either run, `rnabioco_hits_annotated.csv`, the P–[D/E] counts, and
`rnabioco_PDE_new.csv` all update automatically; re-run `compare.py` only if you want
the JSON summary refreshed.

## Files
| file | contents |
|---|---|
| `rnabioco_hits_annotated.csv` | all 7,291 hits: class, db, acc, coords, 2A seq, aln pos2, S2 pos2, PDE flag, in-S2 flags, S2 downstream product |
| `rnabioco_PDE_new.csv` | the 172 P–[D/E] hits new vs. our prior set |
| `compare_summary.json` | all counts above, machine-readable |
| `compare.py` | the comparison script (re-runnable) |
| `fetch_uniprot.py` | resolves UniProt hits via REST (already run) |
| `resolve_mgnify.py` | resolves MGnify hits from local `mgy_proteins_*.fa.gz` |
| `resolve_imgvr.py` | resolves IMG/VR hits from local IMG/VR protein FASTA |
