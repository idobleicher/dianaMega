#!/usr/bin/env python3
"""Compare rnabioco/2a-peptide-search pipeline hits against our prior Rao et al.
Table S2 UBR3 P-[D/E] downstream screen."""
import gzip, json, csv, os, re, collections, sys

RNA = sys.argv[1]      # scratchpad/2a-peptide-search
PRIOR = sys.argv[2]    # 2a_analysis dir (relative, cwd handles unicode)
OUT = sys.argv[3]      # output dir (relative)
os.makedirs(OUT, exist_ok=True)

# ---------- 1. load rnabioco Stockholm hits ----------
def db_of(name):
    if name.startswith(('tr|', 'sp|')): return 'UniProtKB'
    if name.startswith('MGYP'):         return 'MGnify'
    if name.startswith('IMGVR'):        return 'IMGVR'
    return 'other'

def acc_of(name):
    # tr|ACC|ENTRY/coords  -> ACC ; else strip trailing /coords
    if name.startswith(('tr|', 'sp|')):
        return name.split('|')[1]
    base = name.rsplit('/', 1)[0]
    return base

def coords_of(name):
    m = re.search(r'/(\d+)-(\d+)$', name)
    return (int(m.group(1)), int(m.group(2))) if m else (None, None)

def pos2_from_seq(s):
    """downstream position-2 residue = char after terminal skip-proline (…PG|P)."""
    i = s.rfind('PGP')
    if i == -1:
        return None, None            # no terminal PGP captured
    skipP = i + 2
    p2 = i + 3
    return (s[p2] if p2 < len(s) else None), skipP

def load_sto(path, cls):
    hits = []
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            if line.startswith('#') or line.startswith('//') or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            raw = parts[0]
            aln = parts[1]
            name = re.sub(r'^\d+\|', '', raw)          # strip uniquifying NNNN|
            seq = aln.replace('-', '').replace('.', '').upper()
            p2, _ = pos2_from_seq(seq)
            hits.append(dict(cls=cls, name=name, db=db_of(name),
                             acc=acc_of(name), coords=coords_of(name),
                             seq=seq, aln_pos2=p2))
    return hits

rna = load_sto(os.path.join(RNA, 'legacy/db-searches-historical/combined-class-1.sto.gz'), '1')
rna += load_sto(os.path.join(RNA, 'legacy/db-searches-historical/combined-class-2.sto.gz'), '2')

# ---------- 2. load prior Table S2 sets ----------
prior = json.load(open(os.path.join(PRIOR, 'parsed_rows.json'), encoding='utf-8'))
# downstream_all: adds resolved pos2 / downstream_product / found
downstream = {}
with open(os.path.join(PRIOR, 'downstream_all.csv'), encoding='utf-8') as fh:
    for r in csv.DictReader(fh):
        downstream[r['aid']] = r
# the 893 P-[D/E] hits
pde_hits = []
with open(os.path.join(PRIOR, 'motif_PDE_hits.csv'), encoding='utf-8') as fh:
    pde_hits = list(csv.DictReader(fh))

prior_acc = set(r['acc'] for r in prior if r.get('acc'))
prior_aid = set(r['aid'] for r in prior)
prior_seq = set(r['seq'].upper() for r in prior if r.get('seq'))
pde_acc = set(r['aid'].split('|')[1] if r['aid'].startswith(('tr|','sp|')) else r['aid'].rsplit('/',1)[0]
              for r in pde_hits)
pde_seq = set(r['twoA_seq'].upper() for r in pde_hits)

# ---------- 3. compare ----------
rna_acc = set(h['acc'] for h in rna)
rna_seq = set(h['seq'] for h in rna)
rna_aid = set(h['name'] for h in rna)

def bd(hits):  # db breakdown
    return dict(collections.Counter(h['db'] for h in hits))

# per-hit annotation
for h in rna:
    h['in_S2_acc'] = h['acc'] in prior_acc
    h['in_S2_aid'] = h['name'] in prior_aid
    h['in_S2_seq'] = h['seq'] in prior_seq
    # resolved downstream from prior table if aid matches
    d = downstream.get(h['name'])
    h['S2_pos2'] = d['pos2'] if d else ''
    h['S2_downstream'] = d['downstream_product'] if d else ''
    # PDE call: prefer prior resolved, else alignment-derived
    p2 = (d['pos2'] if d and d['pos2'] else h['aln_pos2']) or ''
    h['pos2_final'] = p2
    h['is_PDE'] = p2 in ('D', 'E')

# write per-hit csv
with open(os.path.join(OUT, 'rnabioco_hits_annotated.csv'), 'w', newline='', encoding='utf-8') as fh:
    w = csv.writer(fh)
    w.writerow(['class','name','db','acc','start','end','twoA_seq','aln_pos2',
                'S2_pos2','pos2_final','is_PDE','in_S2_acc','in_S2_aid','in_S2_seq','S2_downstream'])
    for h in rna:
        s,e = h['coords']
        w.writerow([h['cls'],h['name'],h['db'],h['acc'],s,e,h['seq'],h['aln_pos2'],
                    h['S2_pos2'],h['pos2_final'],h['is_PDE'],h['in_S2_acc'],h['in_S2_aid'],
                    h['in_S2_seq'],h['S2_downstream']])

# unique-to-rnabioco (by accession, restricted to fetchable/comparable dbs)
rna_pde = [h for h in rna if h['is_PDE']]
rna_pde_new = [h for h in rna_pde if not h['in_S2_aid'] and not h['in_S2_seq']]

R = {}
R['rna_total'] = len(rna)
R['rna_db'] = bd(rna)
R['rna_unique_acc'] = len(rna_acc)
R['rna_unique_seq'] = len(rna_seq)
R['prior_total'] = len(prior)
R['prior_db'] = dict(collections.Counter(r['db'] for r in prior))
R['prior_unique_acc'] = len(prior_acc)
R['prior_unique_seq'] = len(prior_seq)
# overlap
R['overlap_aid'] = len(rna_aid & prior_aid)
R['overlap_acc'] = len(rna_acc & prior_acc)
R['overlap_seq'] = len(rna_seq & prior_seq)
R['rna_hits_in_S2_by_aid'] = sum(h['in_S2_aid'] for h in rna)
R['rna_hits_in_S2_by_seq'] = sum(h['in_S2_seq'] for h in rna)
R['rna_acc_not_in_S2'] = len(rna_acc - prior_acc)
R['rna_seq_not_in_S2'] = len(rna_seq - prior_seq)
R['S2_acc_not_in_rna'] = len(prior_acc - rna_acc)
# UBR3 P-[D/E]
R['rna_PDE_hits'] = len(rna_pde)
R['rna_PDE_db'] = bd(rna_pde)
R['rna_PDE_unique_seq'] = len(set(h['seq'] for h in rna_pde))
R['rna_PDE_in_prior_893_by_seq'] = sum(h['seq'] in pde_seq for h in rna_pde)
R['rna_PDE_new_vs_prior'] = len(rna_pde_new)
R['rna_PDE_new_unique_seq'] = len(set(h['seq'] for h in rna_pde_new))
R['prior_893_PDE'] = len(pde_hits)
R['prior_327_PDE_unique'] = len(pde_seq)

json.dump(R, open(os.path.join(OUT, 'compare_summary.json'), 'w'), indent=2)

# write the NEW P-[D/E] hits unique to rnabioco
with open(os.path.join(OUT, 'rnabioco_PDE_new.csv'), 'w', newline='', encoding='utf-8') as fh:
    w = csv.writer(fh)
    w.writerow(['class','name','db','acc','twoA_seq','pos2_final','source'])
    for h in sorted(rna_pde_new, key=lambda x: x['seq']):
        w.writerow([h['cls'],h['name'],h['db'],h['acc'],h['seq'],h['pos2_final'],
                    'alignment' if not h['S2_pos2'] else 'S2'])

for k,v in R.items():
    print(f"{k}: {v}")
