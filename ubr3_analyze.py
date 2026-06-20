"""UBR3 substrate analysis: candidate lists + classification by 2nd N-terminal residue."""
import json
import pandas as pd, numpy as np

PATH = r'C:\Users\User\Downloads\95473-81 Protein groups statistics (1).xlsx'
df = pd.read_excel(PATH, header=8)
seqs = json.load(open('seqs.json'))

aavs = ['Seq95473 AAVS 1 ', 'Seq95474 AAVS 2', 'Seq95475 AAVS 3']
wt   = ['Seq95476 UBR3 WT 1 ', 'Seq95477 UBR3 WT 2', 'Seq95478 UBR3 WT 3']
ko   = ['Seq95479 UBR3 KO 1 ', 'Seq95480 UBR3 KO 2', 'Seq95481 UBR3 KO 3']
pA = "Student's T-test p-value UBR3 KO_AAVS"; qA = "Student's T-test q-value UBR3 KO_AAVS"; dA = "Student's T-test Difference UBR3 KO_AAVS"
pW = "Student's T-test p-value UBR3 KO_UBR3 WT"; qW = "Student's T-test q-value UBR3 KO_UBR3 WT"; dW = "Student's T-test Difference UBR3 KO_UBR3 WT"
N = 'N.Sequences'

df['acc1'] = df['Protein.Group'].astype(str).str.split(';').str[0].str.strip()
df['mean_AAVS_log2'] = df[aavs].mean(1).round(3)
df['mean_WT_log2']   = df[wt].mean(1).round(3)
df['mean_KO_log2']   = df[ko].mean(1).round(3)

MET_EXCISED = set('ACGPSTV')  # 2nd residue triggers initiator-Met removal by MetAP

def nterm_info(acc):
    s = seqs.get(acc)
    if not s or len(s) < 2:
        return pd.Series([None, None, None, None])
    p1, p2 = s[0], s[1]
    excised = (p1 == 'M' and p2 in MET_EXCISED)
    eff_nterm = p2 if excised else p1          # residue actually exposed at the N-terminus
    return pd.Series([p1, p2, 'yes' if excised else 'no', eff_nterm])

df[['Nterm_pos1', 'SecondAA_pos2', 'Met_excised_pred', 'Effective_Nterm']] = df['acc1'].apply(nterm_info)

multi   = df[N] > 1
listB_m = (df[dA] > 0) & (df[pA] < 0.05) & multi               # broad: KO>AAVS
rescue  = (df[dW] > 0) & (df[pW] < 0.05)                       # KO>WT (rescue agrees)
listA_m = listB_m & rescue                                    # rescue-confirmed
df['rescue_confirmed'] = np.where(listA_m, 'yes', np.where(listB_m, 'no', ''))

cols = ['Protein.Group', 'Genes', 'Protein.Names', 'First.Protein.Description', 'acc1',
        'SecondAA_pos2', 'Nterm_pos1', 'Met_excised_pred', 'Effective_Nterm',
        'mean_AAVS_log2', 'mean_WT_log2', 'mean_KO_log2',
        dA, pA, qA, dW, pW, qW, N, 'rescue_confirmed']

A = df[listA_m][cols].copy().sort_values(dW, ascending=False)
B = df[listB_m][cols].copy().sort_values(dA, ascending=False)

print('=== List A: rescue-confirmed (KO>AAVS AND KO>WT, p<0.05, >=2 pep) ===', len(A))
print('=== List B: broad (KO>AAVS, p<0.05, >=2 pep) ===', len(B))
print('  of which rescue-confirmed:', (B['rescue_confirmed'] == 'yes').sum())
print('  sequence resolved for 2nd-AA in A:', A['SecondAA_pos2'].notna().sum(), '/', len(A))

order = list('ACDEFGHIKLMNPQRSTVWY')
def tab(d, label):
    vc = d['SecondAA_pos2'].value_counts()
    s = pd.Series({r: int(vc.get(r, 0)) for r in order})
    print(f'\n2nd-AA distribution [{label}] (n={int(s.sum())}):')
    for r in order:
        bar = '#' * s[r]
        print(f'  {r}: {s[r]:3d}  {bar}')
    return s

sA = tab(A, 'List A rescue-confirmed')
sB = tab(B, 'List B broad')

# Met-excision summary for rescue-confirmed
print('\nMet-excision (rescue-confirmed): excised =',
      int((A['Met_excised_pred'] == 'yes').sum()), '| retained =', int((A['Met_excised_pred'] == 'no').sum()))

# write Excel workbook
summary = pd.DataFrame({'second_AA': order,
                        'ListA_rescue_confirmed': [int(sA[r]) for r in order],
                        'ListB_broad': [int(sB[r]) for r in order]})
summary['Met_excised_if_pos1M'] = ['yes' if r in MET_EXCISED else 'no' for r in order]

with pd.ExcelWriter('UBR3_substrates_classified.xlsx') as xw:
    A.to_excel(xw, 'A_rescue_confirmed', index=False)
    B.to_excel(xw, 'B_broad_KOvsAAVS', index=False)
    summary.to_excel(xw, 'Summary_by_2ndAA', index=False)
print('\nWrote UBR3_substrates_classified.xlsx (sheets: A_rescue_confirmed, B_broad_KOvsAAVS, Summary_by_2ndAA)')
