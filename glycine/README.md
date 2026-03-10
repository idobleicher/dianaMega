## Find polyprotein cleavage products that start with Glycine (G)

This repo contains a small script that:
- reads viral **polyprotein** sequences from a FASTA file
- applies protease cleavage motifs (3C/3CL examples included)
- "digests" the polyprotein into fragments
- reports the fragments whose **new N-terminus matches what you require** (default: Glycine `G`) after cleavage
- writes results as a **CSV table**

### Files
- `find_glycine_ntermini.py`: main script
- `protease_motifs.json`: default motif definitions (edit to add other proteases)

### Run

```bash
python find_glycine_ntermini.py --fasta your_polyproteins.fasta --motifs protease_motifs.json --out glycine_fragments.csv
```

By default, the script only reports fragments whose N-terminus match was **created by a cleavage**.
If you also want to include a polyprotein that already starts with `G` (no cleavage):

```bash
python find_glycine_ntermini.py --fasta your_polyproteins.fasta --include_polyprotein_nterm
```

If you need **position 1 = Gly (G)** and **position 2 = Glu (E)** (i.e. fragments start with **`GE`**):

```bash
python find_glycine_ntermini.py --fasta your_polyproteins.fasta --nterm GE --out glycine_fragments_GE.csv
```

### Motif config format (`protease_motifs.json`)

Each motif is a regex match window plus where the cut occurs within that window:

- `pattern`: regex (matched against the amino-acid sequence)
- `cut_after`: cut position = `match_start + cut_after`

Example: `"pattern": "QG", "cut_after": 1` means cut between `Q|G`, so the downstream fragment starts with `G`.

### More proteases that can cut before Gly (P1′ = G)

If you mean “after cleavage, the **first residue of the new fragment is Glycine**”, that’s **P1′ = G**.
Common examples (often used in cloning / annotation) include:
- **3C/3CL-like**: `...Q|G...` (and related `...Q|S/A...`)
- **TEV / NIaPro**: `ENLYFQ|G`
- **PreScission (HRV 3C)**: `LEVLFQ|GP`
- **Thrombin**: `LVPR|G...` (often written `LVPR|GS`)

To add another protease, add a new entry in `protease_motifs.json` with an appropriate `pattern` and `cut_after` so the cut lands immediately before `G`.

### Output columns

The CSV includes:
- polyprotein id
- protease motif name
- fragment start/end (1-based)
- fragment length
- N-terminus (and first 5 aa)
- cleavage-site context with a `|` marking the cut


TO RUN:
python .\find_glycine_ntermini.py --fasta .\sample_polyprotein.fasta --motifs .\protease_motifs.json --nterm GE --out .\fragments_GE.csv