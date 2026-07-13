#!/usr/bin/env python3
"""Resolve downstream pos2 for IMG/VR rnabioco hits by streaming the local IMG/VR
protein FASTA. IMG/VR requires a JGI account: download from https://img.jgi.doe.gov/vr/
(file is typically IMGVR_all_proteins.faa.gz).

USAGE:
    python resolve_imgvr.py <compare_dir> <path/to/IMGVR_all_proteins.faa.gz>

IMG/VR hit names in the rnabioco set look like:
    IMGVR_UViG_2700989091_000001|2700989091|2701390664/692-706
The protein-FASTA header format varies by IMG/VR release, so this script matches a
record if ANY pipe-delimited field of the hit name appears in the header's first
token. It prints how many matched so you can adjust `keys_for` if the release differs.
"""
import csv, os, sys, gzip, zipfile, io

COMPARE = sys.argv[1]
DBFILE  = sys.argv[2]
# optional 3rd arg: member path inside a .zip (defaults to the full v7 proteins)
MEMBER  = sys.argv[3] if len(sys.argv) > 3 else \
          "Cus_VR/IMG_VR_2022-12-19_7/IMGVR_all_proteins.faa.gz"

csv_path = os.path.join(COMPARE, 'rnabioco_hits_annotated.csv')
rows = list(csv.DictReader(open(csv_path, encoding='utf-8')))

# The hit `acc` is the FULL protein id (UViG|taxon|geneid) — the IMG/VR protein-FASTA
# header's first token has exactly this form, so match on the full id, exactly.
key_index = {}
n_hits = 0
for r in rows:
    if r['db'] == 'IMGVR' and not r['pos2_final']:
        n_hits += 1
        key_index.setdefault(r['acc'] or r['name'], []).append(r)
print(f"{n_hits} unresolved IMG/VR hits, {len(key_index)} unique protein ids")

def pos2_from_source(seq2a, source):
    idx = source.find(seq2a)
    if idx == -1:
        tail = seq2a[-6:]
        idx = source.find(tail)
        if idx == -1: return None
        endP = idx + len(tail) - 1
    else:
        endP = idx + len(seq2a) - 1
    p2 = endP + 1
    return source[p2] if p2 < len(source) else ''

def stream_fasta(fh):
    name, buf = None, []
    for line in fh:
        if line.startswith('>'):
            if name: yield name, ''.join(buf)
            name = line[1:].rstrip()
            buf = []
        else:
            buf.append(line.strip())
    if name: yield name, ''.join(buf)

if not os.path.exists(DBFILE):
    sys.exit(f"IMG/VR file not found: {DBFILE}")

def open_source(path, member):
    """Yield a text FASTA stream from a .faa/.faa.gz OR a member inside a .zip."""
    if path.endswith('.zip'):
        z = zipfile.ZipFile(path)
        names = z.namelist()
        if member not in names:
            cand = [n for n in names if n.endswith('.faa.gz') or n.endswith('.faa')]
            sys.exit(f"member '{member}' not in zip.\nAvailable protein files:\n  " +
                     "\n  ".join(cand))
        raw = z.open(member)                       # deflate-decompressed .faa.gz bytes
        inner = gzip.GzipFile(fileobj=raw) if member.endswith('.gz') else raw
        return io.TextIOWrapper(inner, encoding='utf-8', errors='replace')
    if path.endswith('.gz'):
        return gzip.open(path, 'rt', encoding='utf-8', errors='replace')
    return open(path, 'rt', encoding='utf-8', errors='replace')

print(f"reading: {DBFILE}" + (f" :: {MEMBER}" if DBFILE.endswith('.zip') else ""))
matched = 0
with open_source(DBFILE, MEMBER) as fh:
    for header, seq in stream_fasta(fh):
        hitrows = key_index.get(header.split()[0])
        if not hitrows:
            continue
        for r in hitrows:
            p2 = pos2_from_source(r['twoA_seq'], seq)
            if p2 is not None:
                r['pos2_final'] = p2
                r['is_PDE'] = 'True' if p2 in ('D', 'E') else 'False'
                matched += 1
print(f"resolved {matched}/{n_hits} IMG/VR hits")

with open(csv_path, 'w', newline='', encoding='utf-8') as fh:
    w = csv.DictWriter(fh, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)

pde = [r for r in rows if r['is_PDE'] == 'True']
resolved = sum(1 for r in rows if r['pos2_final'])
print(f"pos2 resolved now: {resolved}/{len(rows)} ({100*resolved//len(rows)}%)")
print(f"P-[D/E] total now: {len(pde)} (unique seq {len(set(r['twoA_seq'] for r in pde))})")
