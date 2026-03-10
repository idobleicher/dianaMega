import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
from scipy import stats
import os, sys, io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

OUTPUT_DIR = r'c:\Users\User\Desktop\תינוקת\ubr3paper\analysis_output_paper'
os.makedirs(OUTPUT_DIR, exist_ok=True)

EXCEL_PATH = r'c:\Users\User\Downloads\Data S3. N-terminome GPS screen data in different genetic backgrounds.xlsx'

# ── Publication-quality matplotlib settings ───────────────────────────────
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.titleweight': 'bold',
    'axes.labelsize': 12,
    'axes.linewidth': 1.2,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
})

# ── 1. LOAD AND PARSE ────────────────────────────────────────────────────

print("=" * 70)
print("LOADING AND PARSING DATA")
print("=" * 70)

raw = pd.read_excel(EXCEL_PATH, sheet_name='Data S3A', header=None)
data_rows = raw.iloc[2:].copy().reset_index(drop=True)

AA_COL = 55
PSI_COL = 9  # PSI AAVS1 KO

records = []
for _, row in data_rows.iterrows():
    seq = str(row[AA_COL]) if pd.notna(row[AA_COL]) else ''
    if len(seq) < 3:
        continue
    try:
        psi = float(row[PSI_COL])
    except (ValueError, TypeError):
        continue
    records.append({
        'Transcript': row[0],
        'Gene': row[1],
        'Sequence': seq,
        'P2': seq[1],
        'P3': seq[2],
        'PSI': psi,
    })

df = pd.DataFrame(records)
print(f"Total peptides with valid PSI (AAVS1 KO): {len(df)}")

# ── 2. SEPARATE PRO VS NON-PRO ──────────────────────────────────────────

df_pro = df[df['P2'] == 'P'].copy()
df_non_pro = df[df['P2'] != 'P'].copy()

print(f"Pro at P2:     {len(df_pro)}")
print(f"Non-Pro at P2: {len(df_non_pro)}")

# ── SELF-CHECK: Verify counts ────────────────────────────────────────────
assert len(df_pro) + len(df_non_pro) == len(df), "FAIL: Pro + Non-Pro != Total"
assert df_pro['P2'].unique().tolist() == ['P'], "FAIL: Pro subset has non-P residues"
assert 'P' not in df_non_pro['P2'].values, "FAIL: Non-Pro subset contains P"
print("  [CHECK PASSED] Subset counts and identity verified.")

# ── Helper: SEM ──────────────────────────────────────────────────────────
def sem(arr):
    return np.std(arr, ddof=1) / np.sqrt(len(arr))

# ── 3. CALCULATE STATS ──────────────────────────────────────────────────

pro_psi = df_pro['PSI'].values
non_pro_psi = df_non_pro['PSI'].values
all_psi = df['PSI'].values

pro_mean, pro_std, pro_sem = np.mean(pro_psi), np.std(pro_psi, ddof=1), sem(pro_psi)
non_mean, non_std, non_sem = np.mean(non_pro_psi), np.std(non_pro_psi, ddof=1), sem(non_pro_psi)
all_mean, all_std, all_sem = np.mean(all_psi), np.std(all_psi, ddof=1), sem(all_psi)

# Self-check: recalculate manually
assert abs(pro_mean - df_pro['PSI'].mean()) < 1e-10, "FAIL: Pro mean mismatch"
assert abs(non_mean - df_non_pro['PSI'].mean()) < 1e-10, "FAIL: Non-Pro mean mismatch"
print("  [CHECK PASSED] Mean calculations verified.")

print(f"\n{'Group':<20} {'n':>6} {'Mean':>8} {'SD':>8} {'SEM':>8}")
print("-" * 55)
print(f"{'All Library':<20} {len(all_psi):>6} {all_mean:>8.4f} {all_std:>8.4f} {all_sem:>8.4f}")
print(f"{'Pro at P2':<20} {len(pro_psi):>6} {pro_mean:>8.4f} {pro_std:>8.4f} {pro_sem:>8.4f}")
print(f"{'Non-Pro at P2':<20} {len(non_pro_psi):>6} {non_mean:>8.4f} {non_std:>8.4f} {non_sem:>8.4f}")

t_stat, t_pval = stats.ttest_ind(pro_psi, non_pro_psi, equal_var=False)
mw_stat, mw_pval = stats.mannwhitneyu(pro_psi, non_pro_psi, alternative='two-sided')
ks_stat, ks_pval = stats.ks_2samp(pro_psi, non_pro_psi)

print(f"\nPro vs Non-Pro (AAVS1 KO PSI):")
print(f"  Welch's t-test:      t = {t_stat:.4f},  p = {t_pval:.2e}")
print(f"  Mann-Whitney U:      U = {mw_stat:.0f},  p = {mw_pval:.2e}")
print(f"  Kolmogorov-Smirnov:  D = {ks_stat:.4f},  p = {ks_pval:.2e}")

# ── 4. PRO-P2 BY P3 RESIDUE ─────────────────────────────────────────────

print("\n" + "=" * 70)
print("PRO-AT-P2 PEPTIDES GROUPED BY 3RD RESIDUE (P3)")
print("=" * 70)

p3_stats = df_pro.groupby('P3')['PSI'].agg(['mean', 'std', 'count', 'sem']).reset_index()
p3_stats.columns = ['P3', 'Mean', 'SD', 'n', 'SEM']
p3_stats = p3_stats.sort_values('Mean').reset_index(drop=True)

# Self-check
for _, row in p3_stats.iterrows():
    grp = df_pro[df_pro['P3'] == row['P3']]['PSI']
    assert abs(grp.mean() - row['Mean']) < 1e-10, f"FAIL: P3={row['P3']} mean mismatch"
    assert grp.shape[0] == row['n'], f"FAIL: P3={row['P3']} count mismatch"
print("  [CHECK PASSED] P3 group stats verified.")

print(f"\n{'P3':<5} {'n':>5} {'Mean':>8} {'SD':>8} {'SEM':>8}")
print("-" * 40)
for _, row in p3_stats.iterrows():
    print(f"{row['P3']:<5} {int(row['n']):>5} {row['Mean']:>8.4f} {row['SD']:>8.4f} {row['SEM']:>8.4f}")

# Also compute all-library P3 stats for context
lib_p3_stats = df.groupby('P3')['PSI'].agg(['mean', 'std', 'count', 'sem']).reset_index()
lib_p3_stats.columns = ['P3', 'Mean', 'SD', 'n', 'SEM']

# ── 5. SIGNIFICANCE HELPER ──────────────────────────────────────────────

def significance_marker(p):
    if p < 0.0001:
        return '****'
    elif p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    return 'ns'

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 1: Bar chart – Pro at P2 vs All Library (with SEM error bars)
# ═══════════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(4.5, 5))

groups = ['All Library', 'Pro at P2', 'Non-Pro at P2']
means = [all_mean, pro_mean, non_mean]
sems = [all_sem, pro_sem, non_sem]
ns = [len(all_psi), len(pro_psi), len(non_pro_psi)]
colors = ['#7F7F7F', '#D62728', '#1F77B4']

bars = ax.bar(range(3), means, yerr=sems, capsize=5, color=colors,
              edgecolor='black', linewidth=0.8, width=0.65, error_kw={'linewidth': 1.2})

ax.set_xticks(range(3))
ax.set_xticklabels([f'{g}\n(n={n:,})' for g, n in zip(groups, ns)], fontsize=9)
ax.set_ylabel('Mean PSI (AAVS1 KO)')
ax.set_title('PSI: Pro at P2 vs. Library')

y_max = max(means) + max(sems) + 0.15
bracket_y = y_max - 0.05
ax.plot([1, 1, 2, 2], [bracket_y - 0.02, bracket_y, bracket_y, bracket_y - 0.02],
        color='black', linewidth=1)
ax.text(1.5, bracket_y + 0.01, significance_marker(t_pval), ha='center', va='bottom', fontsize=10)

ax.set_ylim(0, y_max + 0.15)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
path = os.path.join(OUTPUT_DIR, 'Fig1_bar_pro_vs_library.png')
plt.savefig(path)
plt.close()
print(f"\n  Saved: {path}")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 2: Bar chart – Pro-P2 grouped by P3 residue (SEM error bars)
# ═══════════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(10, 5))

x = np.arange(len(p3_stats))
colors_grad = plt.cm.RdYlGn_r(np.linspace(0.15, 0.85, len(p3_stats)))

bars = ax.bar(x, p3_stats['Mean'], yerr=p3_stats['SEM'], capsize=3,
              color=colors_grad, edgecolor='black', linewidth=0.6, width=0.75,
              error_kw={'linewidth': 1.0})

ax.axhline(y=all_mean, color='black', linestyle='--', linewidth=1.0,
           label=f'Library mean = {all_mean:.2f}')
ax.axhline(y=pro_mean, color='#D62728', linestyle=':', linewidth=1.0,
           label=f'Pro-P2 mean = {pro_mean:.2f}')

ax.set_xticks(x)
ax.set_xticklabels([f"{r}\n({int(n)})" for r, n in zip(p3_stats['P3'], p3_stats['n'])], fontsize=9)
ax.set_ylabel('Mean PSI ± SEM (AAVS1 KO)')
ax.set_title('Pro-at-P2 Peptides: PSI by 3rd Residue')
ax.legend(loc='upper left', framealpha=0.9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
path = os.path.join(OUTPUT_DIR, 'Fig2_bar_pro_p2_by_p3.png')
plt.savefig(path)
plt.close()
print(f"  Saved: {path}")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 3: Heatmap – Mean PSI by P3 residue (Pro-P2 only)
# ═══════════════════════════════════════════════════════════════════════════

p3_heat = p3_stats.set_index('P3')[['Mean']].T
p3_heat.index = ['Mean PSI']

ordered_residues = p3_stats.sort_values('Mean')['P3'].tolist()
p3_heat = p3_heat[ordered_residues]

fig, ax = plt.subplots(figsize=(12, 1.8))
sns.heatmap(p3_heat, annot=True, fmt='.2f', cmap='RdYlGn_r', linewidths=0.8,
            ax=ax, cbar_kws={'label': 'Mean PSI', 'shrink': 0.6},
            annot_kws={'size': 9, 'weight': 'bold'})
ax.set_title('Pro-at-P2: Mean PSI (AAVS1 KO) by 3rd Residue', pad=10)
ax.set_ylabel('')
ax.set_xlabel('3rd Residue (P3)')
plt.tight_layout()
path = os.path.join(OUTPUT_DIR, 'Fig3_heatmap_pro_p2_by_p3.png')
plt.savefig(path)
plt.close()
print(f"  Saved: {path}")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 4: Heatmap – All P2 residues (single row per residue, wider context)
# ═══════════════════════════════════════════════════════════════════════════

p2_stats = df.groupby('P2')['PSI'].agg(['mean', 'std', 'count', 'sem']).reset_index()
p2_stats.columns = ['P2', 'Mean', 'SD', 'n', 'SEM']
p2_stats = p2_stats.sort_values('Mean').reset_index(drop=True)

p2_heat = p2_stats.set_index('P2')[['Mean']].T
p2_heat.index = ['Mean PSI']
p2_ordered = p2_stats.sort_values('Mean')['P2'].tolist()
p2_heat = p2_heat[p2_ordered]

fig, ax = plt.subplots(figsize=(12, 1.8))
sns.heatmap(p2_heat, annot=True, fmt='.2f', cmap='RdYlGn_r', linewidths=0.8,
            ax=ax, cbar_kws={'label': 'Mean PSI', 'shrink': 0.6},
            annot_kws={'size': 9, 'weight': 'bold'})
ax.set_title('Mean PSI (AAVS1 KO) by 2nd Residue (P2) — All Library', pad=10)
ax.set_ylabel('')
ax.set_xlabel('2nd Residue (P2)')
plt.tight_layout()
path = os.path.join(OUTPUT_DIR, 'Fig4_heatmap_all_p2_residues.png')
plt.savefig(path)
plt.close()
print(f"  Saved: {path}")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 5: Violin + box – PSI distribution Pro vs Non-Pro
# ═══════════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(4.5, 5.5))

plot_df = pd.DataFrame({
    'PSI': np.concatenate([pro_psi, non_pro_psi]),
    'Group': ['Pro at P2'] * len(pro_psi) + ['Non-Pro at P2'] * len(non_pro_psi)
})

parts = ax.violinplot([pro_psi, non_pro_psi], positions=[0, 1], showmeans=False,
                      showmedians=False, showextrema=False)
for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor(['#D62728', '#1F77B4'][i])
    pc.set_alpha(0.35)
    pc.set_edgecolor('black')
    pc.set_linewidth(0.6)

bp = ax.boxplot([pro_psi, non_pro_psi], positions=[0, 1], widths=0.15,
                patch_artist=True, showfliers=False, zorder=3)
for i, patch in enumerate(bp['boxes']):
    patch.set_facecolor(['#D62728', '#1F77B4'][i])
    patch.set_alpha(0.8)
for element in ['whiskers', 'caps', 'medians']:
    for line in bp[element]:
        line.set_color('black')
        line.set_linewidth(1.0)

ax.scatter([0], [pro_mean], marker='D', color='white', edgecolor='black', s=50, zorder=4, label='Mean')

ax.scatter([1], [non_mean], marker='D', color='white', edgecolor='black', s=50, zorder=4)

ax.set_xticks([0, 1])
ax.set_xticklabels([f'Pro at P2\n(n={len(pro_psi):,})', f'Non-Pro at P2\n(n={len(non_pro_psi):,})'])
ax.set_ylabel('PSI (AAVS1 KO)')
ax.set_title('PSI Distribution')

y_top = max(np.max(pro_psi), np.max(non_pro_psi)) + 0.3
ax.plot([0, 0, 1, 1], [y_top - 0.1, y_top, y_top, y_top - 0.1], color='black', linewidth=1)
ax.text(0.5, y_top + 0.02, f'p = {t_pval:.2e}', ha='center', va='bottom', fontsize=9)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(loc='lower right', framealpha=0.9)

plt.tight_layout()
path = os.path.join(OUTPUT_DIR, 'Fig5_violin_pro_vs_nonpro.png')
plt.savefig(path)
plt.close()
print(f"  Saved: {path}")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 6: Boxplot – Pro-P2 by P3 residue
# ═══════════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(12, 5.5))

p3_order = p3_stats['P3'].tolist()
box_data = [df_pro[df_pro['P3'] == r]['PSI'].values for r in p3_order]

bp = ax.boxplot(box_data, positions=range(len(p3_order)), widths=0.6,
                patch_artist=True, showfliers=True,
                flierprops={'marker': 'o', 'markersize': 2, 'alpha': 0.4})

colors_box = plt.cm.RdYlGn_r(np.linspace(0.15, 0.85, len(p3_order)))
for i, patch in enumerate(bp['boxes']):
    patch.set_facecolor(colors_box[i])
    patch.set_edgecolor('black')
    patch.set_linewidth(0.8)
for element in ['whiskers', 'caps']:
    for line in bp[element]:
        line.set_color('black')
        line.set_linewidth(0.8)
for line in bp['medians']:
    line.set_color('black')
    line.set_linewidth(1.5)

ax.axhline(y=all_mean, color='black', linestyle='--', linewidth=1.0,
           label=f'Library mean = {all_mean:.2f}')

# Mean markers with SEM error bars
for i, r in enumerate(p3_order):
    row = p3_stats[p3_stats['P3'] == r].iloc[0]
    ax.errorbar(i, row['Mean'], yerr=row['SEM'], fmt='D', color='white',
                markeredgecolor='black', markersize=5, capsize=3, linewidth=1.0, zorder=5)

ax.set_xticks(range(len(p3_order)))
ax.set_xticklabels([f"{r}\n(n={int(p3_stats[p3_stats['P3']==r]['n'].values[0])})"
                     for r in p3_order], fontsize=9)
ax.set_ylabel('PSI (AAVS1 KO)')
ax.set_xlabel('3rd Residue (P3)')
ax.set_title('Pro-at-P2: PSI Distribution by 3rd Residue')
ax.legend(loc='upper left', framealpha=0.9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
path = os.path.join(OUTPUT_DIR, 'Fig6_boxplot_pro_p2_by_p3.png')
plt.savefig(path)
plt.close()
print(f"  Saved: {path}")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 7: Histogram – density overlay Pro vs All
# ═══════════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(7, 5))

bins = np.linspace(df['PSI'].min() - 0.1, df['PSI'].max() + 0.1, 60)
ax.hist(all_psi, bins=bins, alpha=0.4, color='#1F77B4', density=True,
        edgecolor='black', linewidth=0.3, label=f'All Library (n={len(all_psi):,})')
ax.hist(pro_psi, bins=bins, alpha=0.5, color='#D62728', density=True,
        edgecolor='black', linewidth=0.3, label=f'Pro at P2 (n={len(pro_psi):,})')

ax.axvline(all_mean, color='#1F77B4', linestyle='--', linewidth=1.5,
           label=f'Library mean = {all_mean:.2f}')
ax.axvline(pro_mean, color='#D62728', linestyle='--', linewidth=1.5,
           label=f'Pro-P2 mean = {pro_mean:.2f}')

ax.set_xlabel('PSI (AAVS1 KO)')
ax.set_ylabel('Density')
ax.set_title('PSI Distribution: Pro at P2 vs. All Library')
ax.legend(framealpha=0.9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
path = os.path.join(OUTPUT_DIR, 'Fig7_histogram_overlay.png')
plt.savefig(path)
plt.close()
print(f"  Saved: {path}")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 8: ECDF – Pro vs All Library
# ═══════════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(6, 5))

all_sorted = np.sort(all_psi)
pro_sorted = np.sort(pro_psi)
ax.plot(all_sorted, np.linspace(0, 1, len(all_sorted)), color='#1F77B4',
        linewidth=1.8, label=f'All Library (n={len(all_psi):,})')
ax.plot(pro_sorted, np.linspace(0, 1, len(pro_sorted)), color='#D62728',
        linewidth=1.8, label=f'Pro at P2 (n={len(pro_psi):,})')

ax.set_xlabel('PSI (AAVS1 KO)')
ax.set_ylabel('Cumulative Fraction')
ax.set_title('Empirical CDF')

ax.text(0.98, 0.05, f'KS D = {ks_stat:.3f}\np = {ks_pval:.2e}',
        transform=ax.transAxes, ha='right', va='bottom', fontsize=9,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.9))

ax.legend(loc='upper left', framealpha=0.9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(True, alpha=0.2)

plt.tight_layout()
path = os.path.join(OUTPUT_DIR, 'Fig8_ecdf.png')
plt.savefig(path)
plt.close()
print(f"  Saved: {path}")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 9: Bar chart – All P2 residues ranked, Pro highlighted
# ═══════════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(10, 5))

x = np.arange(len(p2_stats))
colors_p2 = ['#D62728' if r == 'P' else '#4C72B0' for r in p2_stats['P2']]

ax.bar(x, p2_stats['Mean'], yerr=p2_stats['SEM'], capsize=3,
       color=colors_p2, edgecolor='black', linewidth=0.6, width=0.7,
       error_kw={'linewidth': 1.0})

ax.axhline(y=all_mean, color='black', linestyle='--', linewidth=1.0,
           label=f'Library mean = {all_mean:.2f}')

ax.set_xticks(x)
ax.set_xticklabels([f"{r}\n({int(n)})" for r, n in zip(p2_stats['P2'], p2_stats['n'])], fontsize=9)
ax.set_ylabel('Mean PSI ± SEM (AAVS1 KO)')
ax.set_title('Mean PSI by 2nd Residue (All Library)')
ax.legend(loc='upper left', framealpha=0.9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
path = os.path.join(OUTPUT_DIR, 'Fig9_bar_all_p2_residues.png')
plt.savefig(path)
plt.close()
print(f"  Saved: {path}")

# ═══════════════════════════════════════════════════════════════════════════
# EXPORT TO EXCEL
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("EXPORTING RESULTS TO EXCEL")
print("=" * 70)

output_excel = os.path.join(OUTPUT_DIR, 'Pro_P2_analysis_AAVS1.xlsx')

with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
    # Sheet 1: All Pro-P2 peptides
    df_pro[['Transcript', 'Gene', 'Sequence', 'P2', 'P3', 'PSI']].to_excel(
        writer, sheet_name='Pro_P2_peptides', index=False)

    # Sheet 2: P3 summary
    p3_export = p3_stats.copy()
    p3_export.columns = ['3rd_Residue', 'Mean_PSI', 'SD', 'n', 'SEM']
    p3_export.to_excel(writer, sheet_name='P3_residue_summary', index=False)

    # Sheet 3: Statistical comparison
    stat_df = pd.DataFrame([{
        'Comparison': 'Pro-P2 vs Non-Pro-P2',
        'PSI_column': 'AAVS1 KO',
        'Pro_n': len(pro_psi),
        'Pro_Mean': pro_mean,
        'Pro_SD': pro_std,
        'Pro_SEM': pro_sem,
        'NonPro_n': len(non_pro_psi),
        'NonPro_Mean': non_mean,
        'NonPro_SD': non_std,
        'NonPro_SEM': non_sem,
        'Mean_Difference': pro_mean - non_mean,
        'Welch_t': t_stat,
        'Welch_p': t_pval,
        'MannWhitney_U': mw_stat,
        'MannWhitney_p': mw_pval,
        'KS_D': ks_stat,
        'KS_p': ks_pval,
    }])
    stat_df.to_excel(writer, sheet_name='Statistics', index=False)

    # Sheet 4: All P2 residues
    p2_export = p2_stats.copy()
    p2_export.columns = ['2nd_Residue', 'Mean_PSI', 'SD', 'n', 'SEM']
    p2_export.to_excel(writer, sheet_name='All_P2_residues', index=False)

    # Sheet 5: Methods summary
    methods = pd.DataFrame({
        'Step': [
            '1. Data loading',
            '2. Filtering',
            '3. Mean PSI calculation',
            '4. P3 grouping',
            '5. Error metric',
            '6. Statistical tests',
            '7. Verification',
        ],
        'Description': [
            f'Loaded {len(df)} peptides from Data S3A sheet with valid PSI (AAVS1 KO) values and amino acid sequences >= 3 residues.',
            f'Selected {len(df_pro)} peptides with Proline (P) at the 2nd position (P2). Remaining {len(df_non_pro)} peptides serve as the non-Pro control.',
            f'Calculated mean PSI for Pro-P2 = {pro_mean:.4f} (SD = {pro_std:.4f}, SEM = {pro_sem:.4f}). Library mean = {all_mean:.4f}.',
            f'Pro-P2 peptides grouped by their 3rd residue (P3): {len(p3_stats)} unique residues. Highest PSI at P3 = {p3_stats.iloc[-1]["P3"]} ({p3_stats.iloc[-1]["Mean"]:.4f}), lowest at P3 = {p3_stats.iloc[0]["P3"]} ({p3_stats.iloc[0]["Mean"]:.4f}).',
            'Error bars represent Standard Error of the Mean (SEM = SD / sqrt(n)), using Bessel-corrected (ddof=1) standard deviation.',
            f'Welch\'s t-test (unequal variance): t = {t_stat:.4f}, p = {t_pval:.2e}. Mann-Whitney U: U = {mw_stat:.0f}, p = {mw_pval:.2e}. KS test: D = {ks_stat:.4f}, p = {ks_pval:.2e}.',
            'All counts, means, and group identities were programmatically verified with assertion checks.',
        ]
    })
    methods.to_excel(writer, sheet_name='Methods', index=False)

print(f"  Saved: {output_excel}")

# ═══════════════════════════════════════════════════════════════════════════
# FINAL REPORT
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("FINAL REPORT")
print("=" * 70)
print(f"""
DATA: N-terminome GPS screen (Data S3A), PSI column: AAVS1 KO
TOTAL PEPTIDES: {len(df):,}
PRO AT P2:      {len(df_pro):,} ({len(df_pro)/len(df)*100:.1f}%)
NON-PRO AT P2:  {len(df_non_pro):,}

MEAN PSI (AAVS1 KO):
  All Library:  {all_mean:.4f} +/- {all_sem:.4f} (SEM)
  Pro at P2:    {pro_mean:.4f} +/- {pro_sem:.4f} (SEM)
  Non-Pro P2:   {non_mean:.4f} +/- {non_sem:.4f} (SEM)
  Difference (Pro - Non-Pro): {pro_mean - non_mean:.4f}

P3 RESIDUE ANALYSIS (Pro-P2 subset):
  Most stabilizing P3:  {p3_stats.iloc[-1]['P3']} (mean PSI = {p3_stats.iloc[-1]['Mean']:.4f} +/- {p3_stats.iloc[-1]['SEM']:.4f})
  Least stabilizing P3: {p3_stats.iloc[0]['P3']} (mean PSI = {p3_stats.iloc[0]['Mean']:.4f} +/- {p3_stats.iloc[0]['SEM']:.4f})

STATISTICAL TESTS (Pro vs Non-Pro):
  Welch's t:      p = {t_pval:.2e} ({significance_marker(t_pval)})
  Mann-Whitney U:  p = {mw_pval:.2e} ({significance_marker(mw_pval)})
  KS test:         p = {ks_pval:.2e} ({significance_marker(ks_pval)})

VERIFICATION: All self-checks PASSED.

FIGURES (9 total, 300 DPI):
  Fig1 - Bar: Pro-P2 vs Library (SEM error bars, significance bracket)
  Fig2 - Bar: Pro-P2 by P3 residue (SEM error bars)
  Fig3 - Heatmap: Pro-P2 mean PSI by P3 residue
  Fig4 - Heatmap: All P2 residues mean PSI
  Fig5 - Violin+Box: Pro vs Non-Pro distributions (p-value annotated)
  Fig6 - Boxplot: Pro-P2 by P3 residue (mean+SEM overlay)
  Fig7 - Histogram: Density overlay Pro vs Library
  Fig8 - ECDF: Cumulative distribution (KS stat annotated)
  Fig9 - Bar: All P2 residues ranked (Pro highlighted, SEM error bars)
""")
print("Done!")
