"""Plot 2nd-AA enrichment and write grouped candidate view."""
import json
import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PATH = r'C:\Users\User\Downloads\95473-81 Protein groups statistics (1).xlsx'
df = pd.read_excel(PATH, header=8)
seqs = json.load(open('seqs.json'))
en = pd.read_excel('UBR3_2ndAA_enrichment.xlsx', sheet_name='enrich_2nd_residue')

# ---- enrichment bar plot (log2 odds ratio, colored by nominal p) ----
en = en.sort_values('odds_ratio', ascending=False)
en['log2OR'] = np.log2(en['odds_ratio'])
colors = ['#c0392b' if (o > 1 and p < 0.05) else '#2980b9' if (o < 1 and p < 0.05) else '#bdc3c7'
          for o, p in zip(en['odds_ratio'], en['p_value'])]
fig, ax = plt.subplots(figsize=(11, 5))
ax.bar(en['residue'], en['log2OR'], color=colors, edgecolor='black', linewidth=0.4)
for i, (r, o, p, q) in enumerate(zip(en['residue'], en['log2OR'], en['p_value'], en['FDR_q'])):
    if p < 0.06:
        ax.text(i, o + (0.03 if o >= 0 else -0.08), f'p={p:.3f}', ha='center', fontsize=8)
ax.axhline(0, color='black', lw=0.8)
ax.set_ylabel('log2 odds ratio  (substrate vs proteome)')
ax.set_xlabel('2nd amino acid (position 2 of canonical sequence)')
ax.set_title('UBR3 rescue-confirmed substrates (n=501): 2nd-residue enrichment\n'
             'red=nominally up, blue=nominally down (p<0.05); NONE survive FDR q<0.05', fontsize=11)
plt.tight_layout()
plt.savefig('UBR3_2ndAA_enrichment.png', dpi=150)
print('wrote UBR3_2ndAA_enrichment.png')

# ---- grouped candidate view: List A genes grouped by 2nd AA ----
pA = "Student's T-test p-value UBR3 KO_AAVS"; dA = "Student's T-test Difference UBR3 KO_AAVS"
pW = "Student's T-test p-value UBR3 KO_UBR3 WT"; dW = "Student's T-test Difference UBR3 KO_UBR3 WT"
N = 'N.Sequences'
df['acc1'] = df['Protein.Group'].astype(str).str.split(';').str[0].str.strip()
df['p2'] = df['acc1'].map(lambda a: (seqs.get(a) or '..')[1] if (seqs.get(a) and len(seqs.get(a)) > 1) else None)
listA = (df[dA] > 0) & (df[pA] < 0.05) & (df[N] > 1) & (df[dW] > 0) & (df[pW] < 0.05)
A = df[listA].copy()
A['minD'] = A[[dA, dW]].min(1)
print('\nRescue-confirmed substrates grouped by 2nd AA (top 8 per residue by min KO-effect):')
for r in list('ALDEGKSVRLPNTFMWQCHIY'):
    g = A[A['p2'] == r].sort_values('minD', ascending=False)
    if len(g) == 0: continue
    names = ', '.join(g['Genes'].head(8).astype(str))
    print(f'  [{r}] n={len(g):3d}: {names}{" ..." if len(g) > 8 else ""}')
