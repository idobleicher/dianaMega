"""Fetch UniProt sequences for UBR3 substrate candidates (cached to seqs.json)."""
import urllib.request, urllib.error, json, time, os
import pandas as pd

PATH = r'C:\Users\User\Downloads\95473-81 Protein groups statistics (1).xlsx'
df = pd.read_excel(PATH, header=8)

pA = "Student's T-test p-value UBR3 KO_AAVS"; dA = "Student's T-test Difference UBR3 KO_AAVS"
pW = "Student's T-test p-value UBR3 KO_UBR3 WT"; dW = "Student's T-test Difference UBR3 KO_UBR3 WT"
N = 'N.Sequences'
multi = df[N] > 1
listB = (df[dA] > 0) & (df[pA] < 0.05) & multi

# first accession of each protein group (representative)
accs = (df.loc[listB, 'Protein.Group'].astype(str)
          .str.split(';').str[0].str.strip().unique().tolist())
print('accessions to fetch:', len(accs))

CACHE = 'seqs.json'
seqs = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
todo = [a for a in accs if a not in seqs]
print('already cached:', len(seqs), '| to fetch:', len(todo))

def fetch_batch(batch):
    q = ','.join(batch)
    url = ('https://rest.uniprot.org/uniprotkb/accessions?accessions=' + q +
           '&fields=accession,sequence&format=tsv')
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'ubr3-analysis'})
            txt = urllib.request.urlopen(req, timeout=60).read().decode()
            out = {}
            for line in txt.strip().split('\n')[1:]:
                parts = line.split('\t')
                if len(parts) >= 2 and parts[1]:
                    out[parts[0]] = parts[1]
            return out
        except Exception as e:
            print('  retry', attempt, repr(e)[:80]); time.sleep(2 + attempt * 2)
    return {}

B = 100
for i in range(0, len(todo), B):
    batch = todo[i:i + B]
    got = fetch_batch(batch)
    seqs.update(got)
    json.dump(seqs, open(CACHE, 'w'))
    print(f'  fetched {i + len(batch)}/{len(todo)}  (+{len(got)})')

# report any unresolved (obsolete/secondary accessions)
missing = [a for a in accs if a not in seqs]
print('DONE. total cached:', len(seqs), '| unresolved:', len(missing))
if missing:
    print('unresolved sample:', missing[:20])
    json.dump(missing, open('unresolved_accs.json', 'w'))
