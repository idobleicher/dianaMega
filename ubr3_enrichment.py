"""Enrichment of 2nd N-terminal residue among rescue-confirmed UBR3 substrates."""
import json
import pandas as pd, numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

PATH = r'C:\Users\User\Downloads\95473-81 Protein groups statistics (1).xlsx'
df = pd.read_excel(PATH, header=8)
seqs = json.load(open('seqs.json'))

pA = "Student's T-test p-value UBR3 KO_AAVS"; dA = "Student's T-test Difference UBR3 KO_AAVS"
pW = "Student's T-test p-value UBR3 KO_UBR3 WT"; dW = "Student's T-test Difference UBR3 KO_UBR3 WT"
N = 'N.Sequences'
df['acc1'] = df['Protein.Group'].astype(str).str.split(';').str[0].str.strip()

MET_EXCISED = set('ACGPSTV')
def info(acc):
    s = seqs.get(acc)
    if not s or len(s) < 2: return (None, None)
    p1, p2 = s[0], s[1]
    eff = p2 if (p1 == 'M' and p2 in MET_EXCISED) else p1
    return (p2, eff)
df['p2'], df['eff'] = zip(*df['acc1'].map(info))

multi   = df[N] > 1
listB   = (df[dA] > 0) & (df[pA] < 0.05) & multi
rescue  = (df[dW] > 0) & (df[pW] < 0.05)
listA   = listB & rescue                                  # rescue-confirmed foreground

detected = df['p2'].notna()                               # all proteins with a resolved sequence
order = list('ACDEFGHIKLMNPQRSTVWY')

def enrich(col, fg_mask, label):
    fg = df[fg_mask & df[col].notna()]
    bg = df[detected & ~fg_mask]                          # rest of detected proteome
    rows = []
    nfg, nbg = len(fg), len(bg)
    for r in order:
        a = int((fg[col] == r).sum()); b = nfg - a
        c = int((bg[col] == r).sum()); d = nbg - c
        orr, p = fisher_exact([[a, b], [c, d]], alternative='two-sided')
        rows.append([r, a, nfg, 100*a/nfg, c, nbg, 100*c/nbg, orr, p])
    t = pd.DataFrame(rows, columns=['residue', 'fg_n', 'fg_tot', 'fg_pct',
                                    'bg_n', 'bg_tot', 'bg_pct', 'odds_ratio', 'p_value'])
    t['FDR_q'] = multipletests(t['p_value'], method='fdr_bh')[1]
    t['enriched'] = np.where((t['odds_ratio'] > 1) & (t['FDR_q'] < 0.05), 'ENRICHED',
                     np.where((t['odds_ratio'] < 1) & (t['FDR_q'] < 0.05), 'depleted', ''))
    t = t.sort_values('odds_ratio', ascending=False)
    print(f'\n===== {label}  (foreground n={nfg}, background n={nbg}) =====')
    print(t.round(3).to_string(index=False))
    return t

t_p2  = enrich('p2',  listA, '2nd residue (literal position 2)')
t_eff = enrich('eff', listA, 'Effective N-terminal residue (after Met-excision)')

# N-degron category grouping on the EFFECTIVE N-terminus
cat = {**{r: 'Type1_basic' for r in 'RKH'},
       **{r: 'Type2_bulkyhydrophobic' for r in 'FLWYI'},
       **{r: 'Tertiary/pre(N,Q,C,D,E)' for r in 'NQCDE'}}
def categ(r): return cat.get(r, 'Stabilizing/other')
for label, mask in [('rescue-confirmed (List A)', listA), ('all detected (background)', detected)]:
    sub = df[mask & df['eff'].notna()]
    vc = sub['eff'].map(categ).value_counts()
    tot = len(sub)
    print(f'\nN-degron category by EFFECTIVE N-term [{label}] n={tot}:')
    for k in ['Type1_basic', 'Type2_bulkyhydrophobic', 'Tertiary/pre(N,Q,C,D,E)', 'Stabilizing/other']:
        n = int(vc.get(k, 0)); print(f'  {k:28s} {n:5d}  ({100*n/tot:4.1f}%)')

with pd.ExcelWriter('UBR3_2ndAA_enrichment.xlsx') as xw:
    t_p2.to_excel(xw, 'enrich_2nd_residue', index=False)
    t_eff.to_excel(xw, 'enrich_effective_Nterm', index=False)
print('\nWrote UBR3_2ndAA_enrichment.xlsx')
