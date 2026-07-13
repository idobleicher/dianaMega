#!/usr/bin/env python3
"""Resolve downstream pos2 for MGnify (MGYP) rnabioco hits by streaming the local
MGnify protein FASTA splits. Never loads the DB into memory.

USAGE:
    python resolve_mgnify.py <compare_dir> <dir_with_mgy_proteins_*.fa.gz>

Files expected (from the pipeline's config):
    mgy_proteins_1.fa.gz ... mgy_proteins_25.fa.gz
    downloaded from
    https://ftp.ebi.ac.uk/pub/databases/metagenomics/peptide_database/current_release/
"""
import csv, os, sys, glob, gzip

COMPARE = sys.argv[1]
DBDIR   = sys.argv[2]

csv_path = os.path.join(COMPARE, 'rnabioco_hits_annotated.csv')
rows = list(csv.DictReader(open(csv_path, encoding='utf-8')))

# needed: acc -> list of hit rows (an accession may host >1 skip site)
need = {}
for r in rows:
    if r['db'] == 'MGnify' and not r['pos2_final']:
        need.setdefault(r['acc'], []).append(r)
print(f"{sum(len(v) for v in need.values())} unresolved MGnify hits / {len(need)} unique MGYP accessions")

def pos2_from_source(seq2a, source):
    idx = source.find(seq2a)
    if idx == -1:
        tail = seq2a[-6:]                 # …NPGP tail
        idx = source.find(tail)
        if idx == -1: return None
        endP = idx + len(tail) - 1
    else:
        endP = idx + len(seq2a) - 1       # last residue of subseq = skip-P
    p2 = endP + 1
    return source[p2] if p2 < len(source) else ''   # '' = 2A at C-terminus

def stream_fasta(fh):
    name, buf = None, []
    for line in fh:
        if line.startswith('>'):
            if name: yield name, ''.join(buf)
            name = line[1:].split()[0].strip()
            buf = []
        else:
            buf.append(line.strip())
    if name: yield name, ''.join(buf)

files = sorted(glob.glob(os.path.join(DBDIR, 'mgy_proteins_*.fa.gz'))) \
        or sorted(glob.glob(os.path.join(DBDIR, '*.fa.gz'))) \
        or sorted(glob.glob(os.path.join(DBDIR, '*.fasta*')))
if not files:
    sys.exit(f"No MGnify FASTA found in {DBDIR}")
print("scanning:", *[os.path.basename(f) for f in files])

remaining = set(need)
found = 0
for fp in files:
    if not remaining: break
    op = gzip.open if fp.endswith('.gz') else open
    with op(fp, 'rt') as fh:
        for name, seq in stream_fasta(fh):
            if name in remaining:
                for r in need[name]:
                    p2 = pos2_from_source(r['twoA_seq'], seq)
                    if p2:
                        r['pos2_final'] = p2
                        r['is_PDE'] = 'True' if p2 in ('D', 'E') else 'False'
                found += 1
                remaining.discard(name)
                if not remaining: break
    print(f"  after {os.path.basename(fp)}: matched {found}/{len(need)}")

print(f"unmatched MGYP accessions: {len(remaining)}")

# rewrite CSV + refresh derived files
with open(csv_path, 'w', newline='', encoding='utf-8') as fh:
    w = csv.DictWriter(fh, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)

pde = [r for r in rows if r['is_PDE'] == 'True']
resolved = sum(1 for r in rows if r['pos2_final'])
print(f"pos2 resolved now: {resolved}/{len(rows)} ({100*resolved//len(rows)}%)")
print(f"P-[D/E] total now: {len(pde)} (unique seq {len(set(r['twoA_seq'] for r in pde))})")
