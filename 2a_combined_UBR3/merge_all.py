#!/usr/bin/env python3
"""Build the MAXIMAL combined 2A dataset from both sources of the same project
(Rao et al. 2025 = rnabioco/2a-peptide-search) and screen it for the UBR3 P-[D/E]
downstream motif. Reports at two levels: INSTANCES and UNIQUE SEQUENCES."""
import csv, os, sys, json, collections

S2DIR    = sys.argv[1]   # 2a_analysis
RNADIR   = sys.argv[2]   # 2a_rnabioco_compare
OUT      = sys.argv[3]   # 2a_combined_UBR3
os.makedirs(OUT, exist_ok=True)

def core(seq):
    """Normalize a 2A peptide to a comparable form: the <=13 residues ending at the
    G immediately before the skip-proline (…P-G↓P). Absorbs the two pipelines'
    different trim windows so the same biological peptide collapses to one key."""
    seq = seq.upper().replace('-', '')
    i = seq.rfind('PG')
    end = (i + 2) if i != -1 else len(seq)   # index just past the pre-skip G
    return seq[:end][-13:]

def parse_coords(aid):
    import re
    m = re.search(r'/(\d+)-(\d+)$', aid)
    return (m.group(1), m.group(2)) if m else ('', '')

recs = []   # unified instance records

# ---- Source 1: paper's full Table S2 (downstream_all.csv, 36,008) ----
with open(os.path.join(S2DIR, 'downstream_all.csv'), encoding='utf-8') as fh:
    for r in csv.DictReader(fh):
        db = r['db']
        acc = (r['aid'].split('|')[1] if r['aid'].startswith(('tr|','sp|'))
               else r['aid'].rsplit('/',1)[0])
        s,e = parse_coords(r['aid'])
        recs.append(dict(src='TableS2', db=db, acc=acc, start=s, end=e,
                         seq=r['twoA_seq'], cls=r['class'], origin=r.get('origin',''),
                         pos2=(r['pos2'] if r['found']=='True' else ''),
                         downstream=r.get('downstream_product',''),
                         organism='', lineage=r.get('lineage','')))

# ---- Source 2: repo snapshot (rnabioco_hits_annotated.csv, 7,291) ----
# organism names we fetched for the PDE candidates
orgmap = {}
opath = os.path.join(RNADIR, 'rnabioco_PDE_candidates_with_organism.csv')
if os.path.exists(opath):
    for r in csv.DictReader(open(opath, encoding='utf-8')):
        orgmap[r['accession']] = r['organism_species']

with open(os.path.join(RNADIR, 'rnabioco_hits_annotated.csv'), encoding='utf-8') as fh:
    for r in csv.DictReader(fh):
        db = 'MGnify' if r['db']=='MGnify' else ('IMGVR' if r['db']=='IMGVR' else 'UniProtKB')
        recs.append(dict(src='repo', db=db, acc=r['acc'], start=r['start'], end=r['end'],
                         seq=r['twoA_seq'], cls=('A' if r['class']=='1' else 'B'),
                         origin='', pos2=r['pos2_final'], downstream=r.get('S2_downstream',''),
                         organism=orgmap.get(r['acc'],''), lineage=''))

# ---- merge into instances: key = (db, acc, core) ----
inst = {}
for r in recs:
    r['core'] = core(r['seq'])
    k = (r['db'], r['acc'], r['core'])
    cur = inst.get(k)
    if cur is None:
        inst[k] = r
    else:
        # merge: keep a resolved pos2 / organism / longer seq if the other has it
        if not cur['pos2'] and r['pos2']:            cur['pos2'] = r['pos2']
        if not cur['downstream'] and r['downstream']: cur['downstream'] = r['downstream']
        if not cur['organism'] and r['organism']:    cur['organism'] = r['organism']
        if not cur['lineage'] and r['lineage']:      cur['lineage'] = r['lineage']
        if len(r['seq']) > len(cur['seq']):           cur['seq'] = r['seq']
        cur['src'] = 'both' if cur['src'] != r['src'] and cur['src'] != 'both' else cur['src']

instances = list(inst.values())

def is_pde(r): return r['pos2'] in ('D','E')

# unique-sequence level: a UNIQUE PEPTIDE = a distinct exact amino-acid string.
# (core-13 is only used above to merge instances of the *same site*; it is far too
#  coarse to count peptides — it collapses whole conserved families into one.)
def ukey(r): return r['seq'].upper()
uniq = {}
for r in instances:
    u = uniq.setdefault(ukey(r), dict(n=0, dbs=set(), srcs=set(), pos2set=set(),
                                      pde=False, rep=r))
    u['n'] += 1; u['dbs'].add(r['db']); u['srcs'].add(r['src'])
    if r['pos2']: u['pos2set'].add(r['pos2'])
    if is_pde(r): u['pde'] = True

pde_inst = [r for r in instances if is_pde(r)]
# a unique P-[D/E] peptide = exact seq that is P-[D/E] in at least one instance
pde_uniq = {}
for r in pde_inst:
    u = pde_uniq.setdefault(ukey(r), dict(n=0, dbs=set(), pos2set=set(), rep=r))
    u['n'] += 1; u['dbs'].add(r['db']); u['pos2set'].add(r['pos2'])

# ---------- write outputs ----------
# all instances
with open(os.path.join(OUT,'combined_all_instances.csv'),'w',newline='',encoding='utf-8') as fh:
    w=csv.writer(fh); w.writerow(['source','class','db','accession','start','end',
        'twoA_seq','core','pos2','is_UBR3_PDE','downstream_product','organism','lineage'])
    for r in instances:
        w.writerow([r['src'],r['cls'],r['db'],r['acc'],r['start'],r['end'],r['seq'],
                    r['core'],r['pos2'],is_pde(r),r['downstream'],r['organism'],r['lineage']])

# UBR3 P-[D/E] instances
with open(os.path.join(OUT,'combined_UBR3_PDE_instances.csv'),'w',newline='',encoding='utf-8') as fh:
    w=csv.writer(fh); w.writerow(['source','class','db','accession','organism',
        'twoA_seq','downstream_P+pos2','pos2','downstream_product','lineage'])
    for r in sorted(pde_inst,key=lambda x:(x['db'],x['seq'])):
        w.writerow([r['src'],r['cls'],r['db'],r['acc'],r['organism'],r['seq'],
                    'P'+r['pos2'],r['pos2'],r['downstream'],r['lineage']])

# UBR3 P-[D/E] unique sequences (distinct exact peptides)
with open(os.path.join(OUT,'combined_UBR3_PDE_unique.csv'),'w',newline='',encoding='utf-8') as fh:
    w=csv.writer(fh); w.writerow(['twoA_seq','n_instances','databases','pos2_observed',
        'example_organism'])
    for c,u in sorted(pde_uniq.items(), key=lambda kv:-kv[1]['n']):
        w.writerow([u['rep']['seq'],u['n'],'|'.join(sorted(u['dbs'])),
                    ''.join(sorted(u['pos2set'])),u['rep']['organism']])

# ---------- summary ----------
resolved = sum(1 for r in instances if r['pos2'])
S = {}
S['TOTAL_instances']           = len(instances)
S['TOTAL_unique_sequences']    = len(uniq)
S['instances_by_db']           = dict(collections.Counter(r['db'] for r in instances))
S['instances_by_source']       = dict(collections.Counter(r['src'] for r in instances))
S['pos2_resolved_instances']   = resolved
S['pos2_resolved_pct']         = round(100*resolved/len(instances),1)
S['UBR3_PDE_instances']        = len(pde_inst)
S['UBR3_PDE_unique_sequences'] = len(pde_uniq)
S['UBR3_PDE_instances_by_db']  = dict(collections.Counter(r['db'] for r in pde_inst))
S['UBR3_PDE_by_pos2']          = dict(collections.Counter(r['pos2'] for r in pde_inst))
json.dump(S, open(os.path.join(OUT,'combined_summary.json'),'w'), indent=2)
for k,v in S.items(): print(f"{k}: {v}")
