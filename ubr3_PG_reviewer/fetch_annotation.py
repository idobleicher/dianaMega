#!/usr/bin/env python
"""Fetch UniProt functional annotation for the 54 candidate UBR3 substrates.

Writes ubr3_PG_reviewer/data/annotation.json  {gene: {...}}
Re-run is cheap: existing entries are kept unless --refresh is passed.
"""
import json
import os
import sys
import time

import pandas as pd
import requests

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
CACHE = os.path.join(DATA, 'annotation.json')
XLSX = r'C:\Users\User\Downloads\Supplemental Data 1.xlsx'
FIELDS = 'accession,gene_primary,protein_name,cc_function,cc_subcellular_location,go_p,length'


def short(text, n=260):
    """Trim a UniProt FUNCTION blurb to a one-line annotation."""
    text = ' '.join(text.split())
    for tag in (' (By similarity)', ' (Probable)', ' (PubMed', ' (By similarity'):
        if tag in text:
            text = text.split(tag)[0]
    if len(text) > n:
        cut = text[:n]
        text = cut[:cut.rfind(' ')] + '...'
    return text.strip().rstrip('.') + '.' if text else ''


def query(gene):
    for q in (f'gene_exact:{gene} AND organism_id:9606 AND reviewed:true',
              f'gene:{gene} AND organism_id:9606 AND reviewed:true',
              f'{gene} AND organism_id:9606'):
        try:
            r = requests.get('https://rest.uniprot.org/uniprotkb/search',
                             params={'query': q, 'fields': FIELDS, 'format': 'json', 'size': 1},
                             timeout=30)
            if r.status_code != 200:
                continue
            res = r.json().get('results') or []
            if res:
                return res[0]
        except Exception as e:              # network hiccup -> try the next query form
            print('   retry after', type(e).__name__, file=sys.stderr)
            time.sleep(2)
    return None


def parse(j):
    desc = j.get('proteinDescription', {})
    name = ((desc.get('recommendedName') or {}).get('fullName') or {}).get('value')
    if not name:
        subs = desc.get('submissionNames') or []
        name = subs[0]['fullName']['value'] if subs else ''
    func, loc = '', ''
    for c in j.get('comments', []):
        if c['commentType'] == 'FUNCTION' and not func and c.get('texts'):
            func = short(c['texts'][0]['value'])
        if c['commentType'] == 'SUBCELLULAR LOCATION' and not loc:
            locs = [x['location']['value'] for x in c.get('subcellularLocations', [])]
            loc = '; '.join(dict.fromkeys(locs))[:120]
    gos = []
    for x in j.get('uniProtKBCrossReferences', []):
        if x.get('database') != 'GO':
            continue
        for p in x.get('properties', []):
            if p.get('key') == 'GoTerm' and p.get('value', '').startswith('P:'):
                gos.append(p['value'][2:])
    return {'uniprot': j['primaryAccession'],
            'protein_name': name,
            'function': func,
            'localization': loc,
            'go_biological_process': '; '.join(gos[:6]),
            'protein_length': j.get('sequence', {}).get('length')}


def main():
    os.makedirs(DATA, exist_ok=True)
    cache = {}
    if os.path.exists(CACHE) and '--refresh' not in sys.argv:
        cache = json.load(open(CACHE, encoding='utf-8'))

    B = pd.read_excel(XLSX, sheet_name=' (B) UBR3 substrates', header=[0, 1])
    genes = [str(g).strip() for g in B[B.columns[1]].tolist()]

    for i, g in enumerate(genes, 1):
        if g in cache:
            continue
        print(f'[{i}/{len(genes)}] {g}', flush=True)
        j = query(g)
        cache[g] = parse(j) if j else {'uniprot': '', 'protein_name': '', 'function': '',
                                       'localization': '', 'go_biological_process': '',
                                       'protein_length': None}
        time.sleep(0.15)

    json.dump(cache, open(CACHE, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    missing = [g for g in genes if not cache[g]['uniprot']]
    print(f'wrote {CACHE}  ({len(genes)} genes, {len(missing)} unresolved: {missing})')


if __name__ == '__main__':
    main()
