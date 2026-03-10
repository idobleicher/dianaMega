import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os, sys, io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

OUTPUT_DIR = r'c:\Users\User\Desktop\תינוקת\ubr3paper\analysis_output'
os.makedirs(OUTPUT_DIR, exist_ok=True)

EXCEL_PATH = r'c:\Users\User\Downloads\Data S3. N-terminome GPS screen data in different genetic backgrounds.xlsx'

plt.rcParams.update({
    'font.size': 12,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'figure.dpi': 150,
    'savefig.dpi': 150,
})

# ── 1. LOAD AND PARSE DATA ──────────────────────────────────────────────────

print("=" * 70)
print("LOADING DATA FROM SHEET: Data S3A")
print("=" * 70)

raw = pd.read_excel(EXCEL_PATH, sheet_name='Data S3A', header=None)

PSI_COLS = {
    'AAVS1 KO': 9,
    'UBR KO #1': 17,
    'UBR KO #2': 25,
    'UBR KO #3': 33,
    'ZYG11B KO': 41,
    'ZER1 KO': 49,
}
AA_COL = 55

data_rows = raw.iloc[2:].copy().reset_index(drop=True)

records = []
for _, row in data_rows.iterrows():
    seq = str(row[AA_COL]) if pd.notna(row[AA_COL]) else ''
    if len(seq) < 3:
        continue
    rec = {
        'Transcript': row[0],
        'Gene': row[1],
        'Sequence': seq,
        'P1': seq[0],
        'P2': seq[1],
        'P3': seq[2],
    }
    for name, col_idx in PSI_COLS.items():
        val = row[col_idx]
        try:
            rec[f'PSI_{name}'] = float(val)
        except (ValueError, TypeError):
            rec[f'PSI_{name}'] = np.nan
    records.append(rec)

df = pd.DataFrame(records)
print(f"Total peptides parsed: {len(df)}")
print(f"PSI columns: {[c for c in df.columns if c.startswith('PSI_')]}")

psi_cols = [c for c in df.columns if c.startswith('PSI_')]
df['PSI_mean_all_conditions'] = df[psi_cols].mean(axis=1)

# ── 2. PRO AT POSITION 2 ────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("SECTION 1: PEPTIDES WITH PROLINE (P) AT POSITION 2")
print("=" * 70)

pro_mask = df['P2'] == 'P'
df_pro = df[pro_mask].copy()
df_non_pro = df[~pro_mask].copy()

print(f"Peptides with Pro at P2: {len(df_pro)}")
print(f"Peptides without Pro at P2: {len(df_non_pro)}")
print(f"Percentage with Pro at P2: {len(df_pro)/len(df)*100:.1f}%")

print("\n--- Mean PSI values for Pro-at-P2 peptides ---")
for col in psi_cols:
    mean_val = df_pro[col].dropna().mean()
    std_val = df_pro[col].dropna().std()
    n = df_pro[col].dropna().shape[0]
    print(f"  {col}: mean={mean_val:.4f}, std={std_val:.4f}, n={n}")

mean_all = df_pro['PSI_mean_all_conditions'].dropna().mean()
std_all = df_pro['PSI_mean_all_conditions'].dropna().std()
print(f"  Overall mean PSI (across conditions): {mean_all:.4f} ± {std_all:.4f}")

print("\n--- Mean PSI values for ALL library peptides (for reference) ---")
for col in psi_cols:
    mean_val = df[col].dropna().mean()
    std_val = df[col].dropna().std()
    print(f"  {col}: mean={mean_val:.4f}, std={std_val:.4f}")

# ── 3. PRO-P2 PEPTIDES GROUPED BY P3 RESIDUE ────────────────────────────────

print("\n" + "=" * 70)
print("SECTION 2: PRO-AT-P2 PEPTIDES GROUPED BY 3RD RESIDUE")
print("=" * 70)

p3_groups = df_pro.groupby('P3')

print(f"\nUnique P3 residues among Pro-P2 peptides: {sorted(df_pro['P3'].unique())}")
print(f"Number of unique P3 residues: {df_pro['P3'].nunique()}")

summary_rows = []
for p3_res, group in sorted(p3_groups, key=lambda x: -len(x[1])):
    row = {'P3_residue': p3_res, 'Count': len(group)}
    for col in psi_cols:
        row[f'{col}_mean'] = group[col].dropna().mean()
        row[f'{col}_std'] = group[col].dropna().std()
    row['PSI_overall_mean'] = group['PSI_mean_all_conditions'].dropna().mean()
    row['PSI_overall_std'] = group['PSI_mean_all_conditions'].dropna().std()
    summary_rows.append(row)
    print(f"\n  P3={p3_res} (n={len(group)}):")
    for col in psi_cols:
        m = group[col].dropna().mean()
        s = group[col].dropna().std()
        print(f"    {col}: {m:.4f} ± {s:.4f}")

df_p3_summary = pd.DataFrame(summary_rows)

# ── 4. GRAPHS ────────────────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("SECTION 3: GENERATING GRAPHS")
print("=" * 70)

# ── GRAPH 1: Bar chart – Mean PSI Pro vs All Library per condition ────────

fig, ax = plt.subplots(figsize=(12, 6))
conditions = [c.replace('PSI_', '') for c in psi_cols]
pro_means = [df_pro[c].dropna().mean() for c in psi_cols]
all_means = [df[c].dropna().mean() for c in psi_cols]

x = np.arange(len(conditions))
width = 0.35
bars1 = ax.bar(x - width/2, all_means, width, label='All Library', color='#4C72B0', edgecolor='black', linewidth=0.5)
bars2 = ax.bar(x + width/2, pro_means, width, label='Pro at P2', color='#DD8452', edgecolor='black', linewidth=0.5)

ax.set_ylabel('Mean PSI')
ax.set_title('Mean PSI: Pro-at-P2 Peptides vs. Entire Library')
ax.set_xticks(x)
ax.set_xticklabels(conditions, rotation=30, ha='right')
ax.legend()
ax.set_ylim(0, max(max(pro_means), max(all_means)) * 1.15)

for bar_group in [bars1, bars2]:
    for bar in bar_group:
        height = bar.get_height()
        ax.annotate(f'{height:.2f}', xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=9)

plt.tight_layout()
path1 = os.path.join(OUTPUT_DIR, '01_bar_pro_vs_all_by_condition.png')
plt.savefig(path1)
plt.close()
print(f"  Saved: {path1}")

# ── GRAPH 2: Grouped bar chart – Pro-P2 by P3 residue (AAVS1 KO PSI) ────

fig, ax = plt.subplots(figsize=(14, 6))
p3_order = df_p3_summary.sort_values('PSI_overall_mean')['P3_residue'].tolist()
means = [df_p3_summary[df_p3_summary['P3_residue'] == r]['PSI_AAVS1 KO_mean'].values[0] for r in p3_order]
stds = [df_p3_summary[df_p3_summary['P3_residue'] == r]['PSI_AAVS1 KO_std'].values[0] for r in p3_order]
counts = [df_p3_summary[df_p3_summary['P3_residue'] == r]['Count'].values[0] for r in p3_order]

colors = plt.cm.RdYlGn(np.linspace(0, 1, len(p3_order)))
bars = ax.bar(range(len(p3_order)), means, yerr=stds, capsize=3,
              color=colors, edgecolor='black', linewidth=0.5)

all_lib_mean = df['PSI_AAVS1 KO'].dropna().mean()
ax.axhline(y=all_lib_mean, color='red', linestyle='--', linewidth=1.5, label=f'Library mean ({all_lib_mean:.2f})')

ax.set_xticks(range(len(p3_order)))
ax.set_xticklabels([f"{r}\n(n={c})" for r, c in zip(p3_order, counts)], fontsize=9)
ax.set_ylabel('Mean PSI (AAVS1 KO)')
ax.set_title('Pro-at-P2 Peptides: Mean PSI by 3rd Residue (AAVS1 KO)')
ax.legend()
plt.tight_layout()
path2 = os.path.join(OUTPUT_DIR, '02_bar_pro_p2_by_p3_residue.png')
plt.savefig(path2)
plt.close()
print(f"  Saved: {path2}")

# ── GRAPH 3: Heatmap – Mean PSI for Pro-P2 by P3 residue × condition ─────

heatmap_data = []
for p3_res in sorted(df_pro['P3'].unique()):
    grp = df_pro[df_pro['P3'] == p3_res]
    row = {'P3': p3_res}
    for col in psi_cols:
        row[col.replace('PSI_', '')] = grp[col].dropna().mean()
    heatmap_data.append(row)

df_heat = pd.DataFrame(heatmap_data).set_index('P3')

fig, ax = plt.subplots(figsize=(12, 10))
sns.heatmap(df_heat, annot=True, fmt='.2f', cmap='RdYlGn_r', linewidths=0.5,
            ax=ax, cbar_kws={'label': 'Mean PSI'})
ax.set_title('Heatmap: Mean PSI of Pro-at-P2 Peptides\nGrouped by P3 Residue × Condition')
ax.set_ylabel('3rd Residue (P3)')
ax.set_xlabel('Condition')
plt.tight_layout()
path3 = os.path.join(OUTPUT_DIR, '03_heatmap_pro_p2_by_p3_condition.png')
plt.savefig(path3)
plt.close()
print(f"  Saved: {path3}")

# ── GRAPH 4: Violin + Strip plot – PSI distribution Pro vs All ────────────

fig, axes = plt.subplots(2, 3, figsize=(18, 12))
axes = axes.flatten()

for idx, col in enumerate(psi_cols):
    ax = axes[idx]
    condition = col.replace('PSI_', '')
    plot_df = pd.DataFrame({
        'PSI': pd.concat([df[col].dropna(), df_pro[col].dropna()]),
        'Group': ['All Library'] * df[col].dropna().shape[0] + ['Pro at P2'] * df_pro[col].dropna().shape[0]
    })
    sns.violinplot(data=plot_df, x='Group', y='PSI', ax=ax, inner='quartile',
                   palette=['#4C72B0', '#DD8452'], cut=0)
    ax.set_title(condition)
    ax.set_xlabel('')

plt.suptitle('PSI Distribution: All Library vs. Pro-at-P2 Peptides', fontsize=16, y=1.01)
plt.tight_layout()
path4 = os.path.join(OUTPUT_DIR, '04_violin_pro_vs_all.png')
plt.savefig(path4)
plt.close()
print(f"  Saved: {path4}")

# ── GRAPH 5: Box plot – Pro-P2 by P3 residue (AAVS1 KO) ─────────────────

fig, ax = plt.subplots(figsize=(16, 7))
order = df_pro.groupby('P3')['PSI_AAVS1 KO'].median().sort_values().index.tolist()
sns.boxplot(data=df_pro, x='P3', y='PSI_AAVS1 KO', order=order, ax=ax,
            palette='viridis', fliersize=2)
ax.axhline(y=all_lib_mean, color='red', linestyle='--', linewidth=1.5, label=f'Library mean ({all_lib_mean:.2f})')
ax.set_title('Boxplot: PSI Distribution of Pro-at-P2 Peptides by 3rd Residue (AAVS1 KO)')
ax.set_xlabel('3rd Residue (P3)')
ax.set_ylabel('PSI (AAVS1 KO)')
ax.legend()
plt.tight_layout()
path5 = os.path.join(OUTPUT_DIR, '05_boxplot_pro_p2_by_p3.png')
plt.savefig(path5)
plt.close()
print(f"  Saved: {path5}")

# ── GRAPH 6: Histogram – PSI distribution overlay ────────────────────────

fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(df['PSI_AAVS1 KO'].dropna(), bins=50, alpha=0.5, color='#4C72B0',
        label=f'All Library (n={df["PSI_AAVS1 KO"].dropna().shape[0]})', density=True, edgecolor='black', linewidth=0.3)
ax.hist(df_pro['PSI_AAVS1 KO'].dropna(), bins=50, alpha=0.6, color='#DD8452',
        label=f'Pro at P2 (n={df_pro["PSI_AAVS1 KO"].dropna().shape[0]})', density=True, edgecolor='black', linewidth=0.3)
ax.axvline(df['PSI_AAVS1 KO'].dropna().mean(), color='#4C72B0', linestyle='--', linewidth=2, label=f'Library mean ({df["PSI_AAVS1 KO"].dropna().mean():.2f})')
ax.axvline(df_pro['PSI_AAVS1 KO'].dropna().mean(), color='#DD8452', linestyle='--', linewidth=2, label=f'Pro-P2 mean ({df_pro["PSI_AAVS1 KO"].dropna().mean():.2f})')
ax.set_xlabel('PSI (AAVS1 KO)')
ax.set_ylabel('Density')
ax.set_title('PSI Distribution: All Library vs. Pro-at-P2 Peptides (AAVS1 KO)')
ax.legend()
plt.tight_layout()
path6 = os.path.join(OUTPUT_DIR, '06_histogram_psi_overlay.png')
plt.savefig(path6)
plt.close()
print(f"  Saved: {path6}")

# ── GRAPH 7: Heatmap – Mean PSI by P2 residue × condition (all residues) ─

p2_groups = df.groupby('P2')
p2_heat_data = []
for p2_res in sorted(df['P2'].unique()):
    grp = df[df['P2'] == p2_res]
    row = {'P2': p2_res}
    for col in psi_cols:
        row[col.replace('PSI_', '')] = grp[col].dropna().mean()
    p2_heat_data.append(row)

df_p2_heat = pd.DataFrame(p2_heat_data).set_index('P2')

fig, ax = plt.subplots(figsize=(12, 10))
sns.heatmap(df_p2_heat, annot=True, fmt='.2f', cmap='RdYlGn_r', linewidths=0.5,
            ax=ax, cbar_kws={'label': 'Mean PSI'})
ax.set_title('Heatmap: Mean PSI by P2 Residue × Condition\n(Proline = P row)')
ax.set_ylabel('2nd Residue (P2)')
ax.set_xlabel('Condition')
plt.tight_layout()
path7 = os.path.join(OUTPUT_DIR, '07_heatmap_all_p2_residues.png')
plt.savefig(path7)
plt.close()
print(f"  Saved: {path7}")

# ── GRAPH 8: Bar chart – Mean PSI by P2 residue (AAVS1 KO) ─────────────

fig, ax = plt.subplots(figsize=(14, 6))
p2_stats = df.groupby('P2')['PSI_AAVS1 KO'].agg(['mean', 'std', 'count']).sort_values('mean')
colors_p2 = ['#DD8452' if idx == 'P' else '#4C72B0' for idx in p2_stats.index]
bars = ax.bar(range(len(p2_stats)), p2_stats['mean'], yerr=p2_stats['std'], capsize=3,
              color=colors_p2, edgecolor='black', linewidth=0.5)
ax.set_xticks(range(len(p2_stats)))
ax.set_xticklabels([f"{r}\n(n={int(c)})" for r, c in zip(p2_stats.index, p2_stats['count'])], fontsize=9)
ax.axhline(y=all_lib_mean, color='red', linestyle='--', linewidth=1.5, label=f'Library mean ({all_lib_mean:.2f})')
ax.set_ylabel('Mean PSI (AAVS1 KO)')
ax.set_title('Mean PSI by 2nd Residue (P2)\nProline highlighted in orange')
ax.legend()
plt.tight_layout()
path8 = os.path.join(OUTPUT_DIR, '08_bar_all_p2_residues.png')
plt.savefig(path8)
plt.close()
print(f"  Saved: {path8}")

# ── GRAPH 9: Cumulative distribution (ECDF) – Pro vs All ────────────────

fig, ax = plt.subplots(figsize=(10, 6))
all_psi_sorted = np.sort(df['PSI_AAVS1 KO'].dropna().values)
pro_psi_sorted = np.sort(df_pro['PSI_AAVS1 KO'].dropna().values)
ax.plot(all_psi_sorted, np.linspace(0, 1, len(all_psi_sorted)), color='#4C72B0', linewidth=2, label='All Library')
ax.plot(pro_psi_sorted, np.linspace(0, 1, len(pro_psi_sorted)), color='#DD8452', linewidth=2, label='Pro at P2')
ax.set_xlabel('PSI (AAVS1 KO)')
ax.set_ylabel('Cumulative Fraction')
ax.set_title('Empirical CDF: PSI Distribution (AAVS1 KO)')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
path9 = os.path.join(OUTPUT_DIR, '09_ecdf_pro_vs_all.png')
plt.savefig(path9)
plt.close()
print(f"  Saved: {path9}")

# ── 5. STATISTICAL TESTS ────────────────────────────────────────────────────

print("\n" + "=" * 70)
print("SECTION 4: STATISTICAL COMPARISONS")
print("=" * 70)

for col in psi_cols:
    condition = col.replace('PSI_', '')
    pro_vals = df_pro[col].dropna().values
    all_vals = df[col].dropna().values
    non_pro_vals = df_non_pro[col].dropna().values

    t_stat, t_pval = stats.ttest_ind(pro_vals, non_pro_vals, equal_var=False)
    ks_stat, ks_pval = stats.ks_2samp(pro_vals, non_pro_vals)
    mw_stat, mw_pval = stats.mannwhitneyu(pro_vals, non_pro_vals, alternative='two-sided')

    print(f"\n  {condition}:")
    print(f"    Pro mean={np.mean(pro_vals):.4f}, Non-Pro mean={np.mean(non_pro_vals):.4f}, Diff={np.mean(pro_vals)-np.mean(non_pro_vals):.4f}")
    print(f"    Welch's t-test:         t={t_stat:.4f}, p={t_pval:.2e}")
    print(f"    KS test:                D={ks_stat:.4f}, p={ks_pval:.2e}")
    print(f"    Mann-Whitney U test:    U={mw_stat:.0f}, p={mw_pval:.2e}")

# ── 6. EXPORT DATA TO EXCEL ─────────────────────────────────────────────────

print("\n" + "=" * 70)
print("SECTION 5: EXPORTING RESULTS")
print("=" * 70)

output_excel = os.path.join(OUTPUT_DIR, 'Pro_P2_analysis_results.xlsx')

with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
    # Sheet 1: Pro-P2 peptides with PSI values
    export_pro = df_pro[['Transcript', 'Gene', 'Sequence', 'P2', 'P3'] + psi_cols + ['PSI_mean_all_conditions']].copy()
    export_pro.to_excel(writer, sheet_name='Pro_P2_peptides', index=False)

    # Sheet 2: P3 residue summary
    df_p3_summary.to_excel(writer, sheet_name='P3_residue_summary', index=False)

    # Sheet 3: Comparison statistics
    stat_rows = []
    for col in psi_cols:
        condition = col.replace('PSI_', '')
        pro_vals = df_pro[col].dropna().values
        non_pro_vals = df_non_pro[col].dropna().values
        t_stat, t_pval = stats.ttest_ind(pro_vals, non_pro_vals, equal_var=False)
        ks_stat, ks_pval = stats.ks_2samp(pro_vals, non_pro_vals)
        mw_stat, mw_pval = stats.mannwhitneyu(pro_vals, non_pro_vals, alternative='two-sided')
        stat_rows.append({
            'Condition': condition,
            'Pro_P2_mean': np.mean(pro_vals),
            'Pro_P2_std': np.std(pro_vals),
            'Pro_P2_n': len(pro_vals),
            'NonPro_mean': np.mean(non_pro_vals),
            'NonPro_std': np.std(non_pro_vals),
            'NonPro_n': len(non_pro_vals),
            'All_Library_mean': df[col].dropna().mean(),
            'All_Library_std': df[col].dropna().std(),
            'Difference': np.mean(pro_vals) - np.mean(non_pro_vals),
            'Welch_t_stat': t_stat,
            'Welch_t_pval': t_pval,
            'KS_D_stat': ks_stat,
            'KS_pval': ks_pval,
            'MannWhitney_U': mw_stat,
            'MannWhitney_pval': mw_pval,
        })
    pd.DataFrame(stat_rows).to_excel(writer, sheet_name='Statistical_tests', index=False)

    # Sheet 4: All P2 residue summary
    p2_summary_rows = []
    for p2_res in sorted(df['P2'].unique()):
        grp = df[df['P2'] == p2_res]
        row_data = {'P2_residue': p2_res, 'Count': len(grp)}
        for col in psi_cols:
            row_data[f'{col}_mean'] = grp[col].dropna().mean()
        p2_summary_rows.append(row_data)
    pd.DataFrame(p2_summary_rows).to_excel(writer, sheet_name='All_P2_residues_summary', index=False)

print(f"  Results saved to: {output_excel}")

# ── 7. FINAL SUMMARY REPORT ─────────────────────────────────────────────────

print("\n" + "=" * 70)
print("FINAL SUMMARY REPORT")
print("=" * 70)

print(f"""
DATA SOURCE: Data S3. N-terminome GPS screen data in different genetic backgrounds.xlsx
SHEET USED: Data S3A

ANALYSIS OVERVIEW:
1. Parsed {len(df)} peptides with amino acid sequences and PSI values
   across 6 conditions: {', '.join(conditions)}.

2. Identified {len(df_pro)} peptides with Proline (P) at the 2nd position
   ({len(df_pro)/len(df)*100:.1f}% of all peptides).

3. Mean PSI for Pro-at-P2 peptides (AAVS1 KO): {df_pro['PSI_AAVS1 KO'].dropna().mean():.4f}
   Mean PSI for all library peptides (AAVS1 KO): {df['PSI_AAVS1 KO'].dropna().mean():.4f}
   Difference: {df_pro['PSI_AAVS1 KO'].dropna().mean() - df['PSI_AAVS1 KO'].dropna().mean():.4f}

4. Pro-at-P2 peptides were grouped by their 3rd residue (P3).
   Found {df_pro['P3'].nunique()} unique P3 residues.
   Highest mean PSI (AAVS1 KO) at P3: {df_p3_summary.loc[df_p3_summary['PSI_AAVS1 KO_mean'].idxmax(), 'P3_residue']} ({df_p3_summary['PSI_AAVS1 KO_mean'].max():.4f})
   Lowest mean PSI (AAVS1 KO) at P3:  {df_p3_summary.loc[df_p3_summary['PSI_AAVS1 KO_mean'].idxmin(), 'P3_residue']} ({df_p3_summary['PSI_AAVS1 KO_mean'].min():.4f})

5. Statistical tests comparing Pro-P2 vs non-Pro-P2 peptides were performed
   using Welch's t-test, Kolmogorov-Smirnov test, and Mann-Whitney U test
   for each condition.

GRAPHS GENERATED:
  01 - Bar chart: Mean PSI Pro-at-P2 vs All Library (per condition)
  02 - Bar chart: Pro-P2 peptides grouped by P3 residue (AAVS1 KO)
  03 - Heatmap: Pro-P2 mean PSI by P3 residue x Condition
  04 - Violin plots: PSI distribution All Library vs Pro-P2 (all conditions)
  05 - Boxplot: Pro-P2 by P3 residue (AAVS1 KO)
  06 - Histogram overlay: PSI distribution Pro-P2 vs All Library
  07 - Heatmap: Mean PSI by ALL P2 residues x Condition (context)
  08 - Bar chart: Mean PSI by ALL P2 residues (Proline highlighted)
  09 - ECDF: Cumulative distribution Pro-P2 vs All Library

OUTPUT FILES:
  Graphs:  {OUTPUT_DIR}
  Excel:   {output_excel}
""")

print("Analysis complete!")
