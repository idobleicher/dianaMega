"""
Export raw PSI values for PD+PE+PT and GE+GD groups (Top Hits vs All Screen)
so the statistical tests can be independently verified.
"""
import pandas as pd
import numpy as np
from scipy import stats
import os, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

OUTPUT_DIR = r'c:\Users\User\Desktop\תינוקת\dianaMega\ubr3enrichmentlogo\dipeptide_analysis'
SCREEN_EXCEL = r'c:\Users\User\Downloads\UBR3 Nt screen (1).xlsx'

COMBINED_GROUPS = {
    'PD+PE+PT': ['PD', 'PE', 'PT'],
    'GE+GD':    ['GE', 'GD'],
}

df_all = pd.read_excel(SCREEN_EXCEL, sheet_name='Nprot_5_analyzed')
df_hits = pd.read_excel(SCREEN_EXCEL, sheet_name='sub_high')

df_all['PSI_AAVS'] = (df_all['PSI-293a'] + df_all['PSI-293b']) / 2
df_hits['PSI_AAVS'] = (df_hits['PSI-293a'] + df_hits['PSI-293b']) / 2
df_all['dipeptide']  = df_all['AA2'].astype(str)  + df_all['AA3'].astype(str)
df_hits['dipeptide'] = df_hits['AA2'].astype(str) + df_hits['AA3'].astype(str)

cols_keep = ['Gene_ID', 'AA_seq', 'AA2', 'AA3', 'dipeptide',
             'PSI-293a', 'PSI-293b', 'PSI_AAVS']
cols_keep = [c for c in cols_keep if c in df_all.columns]

out_xlsx = os.path.join(OUTPUT_DIR, 'RAW_DATA_for_verification.xlsx')

with pd.ExcelWriter(out_xlsx, engine='xlsxwriter') as writer:
    readme = pd.DataFrame({
        'Column / Concept': [
            'Source file',
            'All Screen sheet',
            'Top Hits sheet',
            'PSI_AAVS',
            'dipeptide',
            'Group: PD+PE+PT',
            'Group: GE+GD',
            '',
            'Statistical test (original)',
            'Statistical test (new)',
            'Alternative hypothesis (H1)',
            'Cohen\'s d (effect size)',
        ],
        'Description': [
            SCREEN_EXCEL,
            "'Nprot_5_analyzed' (full library, 16,514 peptides)",
            "'sub_high' (best hits, 54 peptides)",
            'Mean PSI from control replicates: (PSI-293a + PSI-293b) / 2  -- AAVS = control cell line',
            'String = AA2 + AA3 (positions 2 and 3, after the N-terminal Met)',
            'Peptides whose P2-P3 is PD, PE, or PT',
            'Peptides whose P2-P3 is GE or GD',
            '',
            'One-tailed Mann-Whitney U  (scipy.stats.mannwhitneyu, alternative="less")',
            'One-tailed Brunner-Munzel  (scipy.stats.brunnermunzel, alternative="less")',
            'Top Hits have LOWER PSI than All Screen (i.e. more destabilized)',
            'd = (mean_all - mean_hits) / pooled_SD,  pooled_SD = sqrt((var_hits + var_all)/2)',
        ],
    })
    readme.to_excel(writer, sheet_name='README', index=False)

    for grp_name, members in COMBINED_GROUPS.items():
        safe = grp_name.replace('+', '_')
        hits_rows = pd.concat(
            [df_hits[df_hits['dipeptide'] == m] for m in members], ignore_index=True
        )[cols_keep].copy()
        hits_rows.insert(0, 'Group_label', 'Top Hits')
        hits_rows.insert(0, 'Combined_group', grp_name)

        all_rows = pd.concat(
            [df_all[df_all['dipeptide'] == m] for m in members], ignore_index=True
        )[cols_keep].copy()
        all_rows.insert(0, 'Group_label', 'All Screen')
        all_rows.insert(0, 'Combined_group', grp_name)

        combined = pd.concat([hits_rows, all_rows], ignore_index=True)
        combined.to_excel(writer, sheet_name=f'RAW_{safe}', index=False)

    stat_rows = []
    for grp_name, members in COMBINED_GROUPS.items():
        hits_psi = pd.concat(
            [df_hits[df_hits['dipeptide'] == m] for m in members]
        )['PSI_AAVS'].dropna().values
        all_psi = pd.concat(
            [df_all[df_all['dipeptide'] == m] for m in members]
        )['PSI_AAVS'].dropna().values

        mean_h, mean_a = np.mean(hits_psi), np.mean(all_psi)
        sd_h, sd_a = np.std(hits_psi, ddof=1), np.std(all_psi, ddof=1)
        pooled = np.sqrt((np.var(hits_psi, ddof=1) + np.var(all_psi, ddof=1)) / 2)
        d = (mean_a - mean_h) / pooled if pooled > 0 else np.nan

        mw_U, mw_p1 = stats.mannwhitneyu(hits_psi, all_psi, alternative='less')
        _,   mw_p2 = stats.mannwhitneyu(hits_psi, all_psi, alternative='two-sided')
        bm_W, bm_p1 = stats.brunnermunzel(hits_psi, all_psi, alternative='less')
        _,   bm_p2 = stats.brunnermunzel(hits_psi, all_psi, alternative='two-sided')
        t_stat, t_p2 = stats.ttest_ind(hits_psi, all_psi, equal_var=False)
        t_p1 = t_p2 / 2 if t_stat < 0 else 1 - t_p2 / 2

        stat_rows.append({
            'Group': grp_name,
            'n_hits': len(hits_psi), 'n_all': len(all_psi),
            'mean_hits': mean_h, 'mean_all': mean_a,
            'SD_hits': sd_h, 'SD_all': sd_a,
            'mean_diff (hits-all)': mean_h - mean_a,
            'Cohen_d': d,
            'MannWhitney_U': mw_U, 'MW_p_1tail': mw_p1, 'MW_p_2tail': mw_p2,
            'BrunnerMunzel_W': bm_W, 'BM_p_1tail': bm_p1, 'BM_p_2tail': bm_p2,
            'Welch_t': t_stat, 'Welch_p_1tail': t_p1, 'Welch_p_2tail': t_p2,
        })
    pd.DataFrame(stat_rows).to_excel(writer, sheet_name='Statistics_reproduced', index=False)

print(f"Wrote: {out_xlsx}")

print("\n" + "=" * 80)
print("RAW PSI_AAVS VALUES (for hand verification)")
print("=" * 80)
for grp_name, members in COMBINED_GROUPS.items():
    hits_rows = pd.concat(
        [df_hits[df_hits['dipeptide'] == m] for m in members], ignore_index=True
    )
    all_rows = pd.concat(
        [df_all[df_all['dipeptide'] == m] for m in members], ignore_index=True
    )
    print(f"\n--- {grp_name} : TOP HITS  (n={len(hits_rows)}) ---")
    display_cols = [c for c in ['Gene_ID', 'dipeptide', 'PSI-293a', 'PSI-293b', 'PSI_AAVS']
                    if c in hits_rows.columns]
    print(hits_rows[display_cols].to_string(index=False, float_format='{:.4f}'.format))
    vals = hits_rows['PSI_AAVS'].dropna().values
    print(f"   values: {np.round(vals, 4).tolist()}")
    print(f"   mean = {vals.mean():.4f}, SD = {vals.std(ddof=1):.4f}")

    print(f"\n--- {grp_name} : ALL SCREEN  (n={len(all_rows)}) ---")
    vals_all = all_rows['PSI_AAVS'].dropna().values
    print(f"   values: {np.round(vals_all, 4).tolist()}")
    print(f"   mean = {vals_all.mean():.4f}, SD = {vals_all.std(ddof=1):.4f}")
