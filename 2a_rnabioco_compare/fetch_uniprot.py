#!/usr/bin/env python3
"""Resolve downstream pos2 for UniProt rnabioco hits whose alignment truncated at
the skip-proline, by fetching the full source protein and reading the residue
after the C-terminal PG|P."""
import csv, json, os, sys, time, urllib.request, urllib.parse

OUT = sys.argv[1]           # 2a_rnabioco_compare dir
CACHE = sys.argv[2]         # scratchpad fasta cache path
rows = list(csv.DictReader(open(os.path.join(OUT, 'rnabioco_hits_annotated.csv'), encoding='utf-8')))

need = [r for r in rows if r['db'] == 'UniProtKB' and not r['pos2_final']]
accs = sorted(set(r['acc'] for r in need))
print(f"{len(need)} unresolved UniProt hits / {len(accs)} unique accessions")

cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}

def fetch_chunk(chunk):
    url = "https://rest.uniprot.org/uniprotkb/accessions?" + urllib.parse.urlencode(
        {"accessions": ",".join(chunk), "format": "fasta"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read().decode()
        except Exception as e:
            print("  retry", attempt, e); time.sleep(3)
    return ""

def parse_fasta(txt):
    d, name, buf = {}, None, []
    for line in txt.splitlines():
        if line.startswith('>'):
            if name: d[name] = ''.join(buf)
            name = line[1:].split('|')[1] if '|' in line else line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name: d[name] = ''.join(buf)
    return d

todo = [a for a in accs if a not in cache]
print(f"fetching {len(todo)} not cached")
for i in range(0, len(todo), 200):
    chunk = todo[i:i+200]
    d = parse_fasta(fetch_chunk(chunk))
    for a in chunk:
        cache[a] = d.get(a, "")     # "" = obsolete/deleted
    print(f"  {i+len(chunk)}/{len(todo)}  (+{sum(1 for a in chunk if cache[a])} found)")
    json.dump(cache, open(CACHE, 'w'))
json.dump(cache, open(CACHE, 'w'))

# resolve pos2: locate 2A subseq in source, read residue after terminal skip-P
def resolve(seq2a, source):
    if not source: return '', 'obsolete'
    idx = source.find(seq2a)
    if idx == -1:
        # try locating just the …NPGP tail
        tail = seq2a[-6:]
        idx = source.find(tail)
        if idx == -1: return '', 'not_found'
        endP = idx + len(tail) - 1
    else:
        endP = idx + len(seq2a) - 1        # last residue of subseq = skip-P
    p2 = endP + 1
    return (source[p2] if p2 < len(source) else ''), ('resolved' if p2 < len(source) else 'C-terminus')

nres = 0
stat = {}
for r in need:
    p2, st = resolve(r['twoA_seq'], cache.get(r['acc'], ''))
    stat[st] = stat.get(st, 0) + 1
    if p2:
        r['pos2_final'] = p2
        r['is_PDE'] = 'True' if p2 in ('D', 'E') else 'False'
        r['S2_pos2'] = r['S2_pos2'] or ''
        nres += 1
print("resolve outcomes:", stat, "| newly resolved:", nres)

# rewrite annotated csv
with open(os.path.join(OUT, 'rnabioco_hits_annotated.csv'), 'w', newline='', encoding='utf-8') as fh:
    w = csv.DictWriter(fh, fieldnames=rows[0].keys())
    w.writeheader(); w.writerows(rows)

# refresh new-PDE csv + a small summary delta
pde = [r for r in rows if r['is_PDE'] == 'True']
new_pde = [r for r in pde if r['in_S2_aid'] == 'False' and r['in_S2_seq'] == 'False']
with open(os.path.join(OUT, 'rnabioco_PDE_new.csv'), 'w', newline='', encoding='utf-8') as fh:
    w = csv.writer(fh)
    w.writerow(['class','name','db','acc','twoA_seq','pos2_final','source'])
    for r in sorted(new_pde, key=lambda x: x['twoA_seq']):
        w.writerow([r['class'],r['name'],r['db'],r['acc'],r['twoA_seq'],r['pos2_final'],
                    'S2' if r['S2_pos2'] else 'alignment/source'])

resolved_total = sum(1 for r in rows if r['pos2_final'])
print(f"pos2 resolved now: {resolved_total}/{len(rows)} ({100*resolved_total//len(rows)}%)")
print(f"P-[D/E] total now: {len(pde)}  (unique seq {len(set(r['twoA_seq'] for r in pde))})")
print(f"new-vs-prior P-[D/E]: {len(new_pde)}  (unique seq {len(set(r['twoA_seq'] for r in new_pde))})")
