import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os, sys, io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

OUTPUT_DIR = r'c:\Users\User\Desktop\תינוקת\ubr3paper\analysis_output_unstable'
os.makedirs(OUTPUT_DIR, exist_ok=True)

EXCEL_PATH = r'c:\Users\User\Downloads\Data S3. N-terminome GPS screen data in different genetic backgrounds.xlsx'

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

def sem(arr):
    arr = np.array(arr)
    return np.std(arr, ddof=1) / np.sqrt(len(arr)) if len(arr) > 1 else 0

def significance_marker(p):
    if p < 0.0001: return '****'
    elif p < 0.001: return '***'
    elif p < 0.01: return '**'
    elif p < 0.05: return '*'
    return 'ns'

# ══════════════════════════════════════════════════════════════════════════
# LOAD DATA
# ══════════════════════════════════════════════════════════════════════════

raw = pd.read_excel(EXCEL_PATH, sheet_name='Data S3A', header=None)
data_rows = raw.iloc[2:].copy().reset_index(drop=True)

records = []
for _, row in data_rows.iterrows():
    seq = str(row[55]) if pd.notna(row[55]) else ''
    if len(seq) < 3:
        continue
    try:
        psi = float(row[9])
    except (ValueError, TypeError):
        continue
    records.append({
        'Transcript': row[0], 'Gene': row[1], 'Sequence': seq,
        'P2': seq[1], 'P3': seq[2], 'PSI': psi,
    })

df = pd.DataFrame(records)
df_pro = df[df['P2'] == 'P'].copy()
df_non_pro = df[df['P2'] != 'P'].copy()

all_mean = df['PSI'].mean()
Q1 = df['PSI'].quantile(0.25)
Q1_pro = df_pro['PSI'].quantile(0.25)

print("=" * 70)
print("DATA LOADED")
print("=" * 70)
print(f"Total: {len(df)}, Pro-P2: {len(df_pro)}, Non-Pro-P2: {len(df_non_pro)}")
print(f"Library Q1 (25th percentile): {Q1:.4f}")
print(f"Pro-P2 Q1: {Q1_pro:.4f}")
print(f"Library mean: {all_mean:.4f}")

# ══════════════════════════════════════════════════════════════════════════
# DEFINE UNSTABLE: PSI below library Q1 (bottom 25%)
# ══════════════════════════════════════════════════════════════════════════

THRESHOLD = Q1
df_pro_unstable = df_pro[df_pro['PSI'] <= THRESHOLD].copy()
df_pro_stable = df_pro[df_pro['PSI'] > THRESHOLD].copy()

print(f"\nUnstable threshold (Library Q1): PSI <= {THRESHOLD:.4f}")
print(f"Unstable Pro-P2: {len(df_pro_unstable)} ({len(df_pro_unstable)/len(df_pro)*100:.1f}%)")
print(f"Stable Pro-P2:   {len(df_pro_stable)} ({len(df_pro_stable)/len(df_pro)*100:.1f}%)")

lib_frac_unstable = (df['PSI'] <= THRESHOLD).mean()
pro_frac_unstable = (df_pro['PSI'] <= THRESHOLD).mean()
print(f"Fraction unstable in library: {lib_frac_unstable:.3f}")
print(f"Fraction unstable in Pro-P2:  {pro_frac_unstable:.3f}")

# ══════════════════════════════════════════════════════════════════════════
# FIG 1: Unstable Pro-P2 — P3 residue composition vs library
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("FIG 1: P3 COMPOSITION OF UNSTABLE PRO-P2 PEPTIDES")
print("=" * 70)

unstable_p3 = df_pro_unstable.groupby('P3')['PSI'].agg(['mean', 'count', 'std']).reset_index()
unstable_p3['SEM'] = unstable_p3.apply(lambda r: r['std']/np.sqrt(r['count']) if r['count'] > 1 else 0, axis=1)
unstable_p3 = unstable_p3.sort_values('mean').reset_index(drop=True)

fig, ax = plt.subplots(figsize=(11, 5.5))
x = np.arange(len(unstable_p3))
colors = plt.cm.RdYlGn_r(np.linspace(0.1, 0.6, len(unstable_p3)))

ax.bar(x, unstable_p3['mean'], yerr=unstable_p3['SEM'], capsize=3,
       color=colors, edgecolor='black', linewidth=0.6, width=0.7,
       error_kw={'linewidth': 1.0})

ax.axhline(THRESHOLD, color='red', linestyle='--', linewidth=1.2,
           label=f'Unstable threshold (Q1 = {THRESHOLD:.2f})')

ax.set_xticks(x)
ax.set_xticklabels([f"{r}\n(n={int(n)})" for r, n in
                     zip(unstable_p3['P3'], unstable_p3['count'])], fontsize=9)
ax.set_ylabel('Mean PSI ± SEM (AAVS1 KO)')
ax.set_xlabel('3rd Residue (P3)')
ax.set_title(f'Unstable Pro-P2 Peptides (PSI ≤ {THRESHOLD:.2f}): Mean PSI by P3')
ax.legend(framealpha=0.9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'Fig01_unstable_pro_p2_by_p3.png'))
plt.close()
print("  Saved Fig01")

# ══════════════════════════════════════════════════════════════════════════
# FIG 2: Fraction of Pro-P2 peptides that are unstable, by P3 residue
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("FIG 2: FRACTION UNSTABLE BY P3 RESIDUE")
print("=" * 70)

frac_data = []
for p3 in sorted(df_pro['P3'].unique()):
    grp = df_pro[df_pro['P3'] == p3]
    n_total = len(grp)
    n_unstable = (grp['PSI'] <= THRESHOLD).sum()
    frac = n_unstable / n_total if n_total > 0 else 0
    ci = 1.96 * np.sqrt(frac * (1 - frac) / n_total) if n_total > 0 else 0
    frac_data.append({'P3': p3, 'n': n_total, 'n_unstable': n_unstable,
                      'frac_unstable': frac, 'CI_95': ci})

frac_df = pd.DataFrame(frac_data).sort_values('frac_unstable', ascending=False).reset_index(drop=True)

print(f"{'P3':<5} {'n':>5} {'n_unst':>7} {'frac':>8} {'CI95':>8}")
for _, r in frac_df.iterrows():
    print(f"{r['P3']:<5} {int(r['n']):>5} {int(r['n_unstable']):>7} {r['frac_unstable']:>8.3f} {r['CI_95']:>8.3f}")

fig, ax = plt.subplots(figsize=(11, 5.5))
x = np.arange(len(frac_df))
colors_frac = plt.cm.Reds(np.linspace(0.3, 0.9, len(frac_df)))

ax.bar(x, frac_df['frac_unstable'], yerr=frac_df['CI_95'], capsize=3,
       color=colors_frac, edgecolor='black', linewidth=0.6, width=0.7,
       error_kw={'linewidth': 1.0})

ax.axhline(lib_frac_unstable, color='black', linestyle='--', linewidth=1.2,
           label=f'Library baseline ({lib_frac_unstable:.2f})')
ax.axhline(pro_frac_unstable, color='#D62728', linestyle=':', linewidth=1.2,
           label=f'Pro-P2 overall ({pro_frac_unstable:.2f})')

ax.set_xticks(x)
ax.set_xticklabels([f"{r}\n({int(n)})" for r, n in zip(frac_df['P3'], frac_df['n'])], fontsize=9)
ax.set_ylabel('Fraction Unstable (PSI ≤ Q1)')
ax.set_xlabel('3rd Residue (P3)')
ax.set_title('Fraction of Pro-P2 Peptides that are Unstable, by P3 Residue')
ax.legend(framealpha=0.9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'Fig02_fraction_unstable_by_p3.png'))
plt.close()
print("  Saved Fig02")

# ══════════════════════════════════════════════════════════════════════════
# FIG 3: PAIRED COMPARISON — Pro-P2 vs Non-Pro-P2, same P3 residue
#   This isolates the EFFECT OF PROLINE at P2, controlling for P3 identity
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("FIG 3: PAIRED P3 ANALYSIS (Pro-P2 vs Non-Pro-P2, same P3)")
print("=" * 70)

paired_data = []
for p3 in sorted(df['P3'].unique()):
    pro_grp = df_pro[df_pro['P3'] == p3]['PSI']
    non_grp = df_non_pro[df_non_pro['P3'] == p3]['PSI']
    if len(pro_grp) < 5 or len(non_grp) < 5:
        continue
    t, p = stats.ttest_ind(pro_grp, non_grp, equal_var=False)
    paired_data.append({
        'P3': p3,
        'Pro_mean': pro_grp.mean(), 'Pro_SEM': sem(pro_grp), 'Pro_n': len(pro_grp),
        'NonPro_mean': non_grp.mean(), 'NonPro_SEM': sem(non_grp), 'NonPro_n': len(non_grp),
        'Delta': pro_grp.mean() - non_grp.mean(),
        'p_value': p, 'sig': significance_marker(p),
    })

paired_df = pd.DataFrame(paired_data).sort_values('Delta').reset_index(drop=True)

print(f"{'P3':<4} {'Pro mean':>9} {'NonPro':>9} {'Delta':>8} {'p':>10} {'sig':>5}")
for _, r in paired_df.iterrows():
    print(f"{r['P3']:<4} {r['Pro_mean']:>9.4f} {r['NonPro_mean']:>9.4f} {r['Delta']:>8.4f} {r['p_value']:>10.2e} {r['sig']:>5}")

fig, ax = plt.subplots(figsize=(11, 5.5))
x = np.arange(len(paired_df))
width = 0.35

bars1 = ax.bar(x - width/2, paired_df['Pro_mean'], width, yerr=paired_df['Pro_SEM'],
               capsize=3, color='#D62728', edgecolor='black', linewidth=0.6,
               label='Pro at P2', error_kw={'linewidth': 1.0})
bars2 = ax.bar(x + width/2, paired_df['NonPro_mean'], width, yerr=paired_df['NonPro_SEM'],
               capsize=3, color='#1F77B4', edgecolor='black', linewidth=0.6,
               label='Non-Pro at P2', error_kw={'linewidth': 1.0})

for i, row in paired_df.iterrows():
    y_max = max(row['Pro_mean'] + row['Pro_SEM'], row['NonPro_mean'] + row['NonPro_SEM'])
    if row['p_value'] < 0.05:
        ax.text(i, y_max + 0.08, row['sig'], ha='center', va='bottom', fontsize=8, fontweight='bold')

ax.set_xticks(x)
ax.set_xticklabels(paired_df['P3'], fontsize=10)
ax.set_ylabel('Mean PSI ± SEM (AAVS1 KO)')
ax.set_xlabel('3rd Residue (P3)')
ax.set_title('Pro-P2 vs. Non-Pro-P2: Same 3rd Residue\n(isolates the effect of Proline at P2)')
ax.legend(framealpha=0.9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'Fig03_paired_pro_vs_nonpro_by_p3.png'))
plt.close()
print("  Saved Fig03")

# ══════════════════════════════════════════════════════════════════════════
# FIG 4: DELTA PSI (Pro-P2 minus Non-Pro-P2) for each P3 residue
#   Positive = Pro is MORE stable, Negative = Pro is LESS stable
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("FIG 4: DELTA PSI (Pro effect, controlling for P3)")
print("=" * 70)

fig, ax = plt.subplots(figsize=(11, 5))
x = np.arange(len(paired_df))
colors_delta = ['#D62728' if d < 0 else '#2CA02C' for d in paired_df['Delta']]

prop_sem = [np.sqrt(r['Pro_SEM']**2 + r['NonPro_SEM']**2) for _, r in paired_df.iterrows()]

ax.bar(x, paired_df['Delta'], yerr=prop_sem, capsize=3,
       color=colors_delta, edgecolor='black', linewidth=0.6, width=0.7,
       error_kw={'linewidth': 1.0})

ax.axhline(0, color='black', linewidth=1.0)

for i, (_, row) in enumerate(paired_df.iterrows()):
    if row['p_value'] < 0.05:
        y_pos = row['Delta'] + prop_sem[i] + 0.02 if row['Delta'] > 0 else row['Delta'] - prop_sem[i] - 0.06
        ax.text(i, y_pos, row['sig'], ha='center', va='bottom' if row['Delta'] > 0 else 'top',
                fontsize=8, fontweight='bold')

ax.set_xticks(x)
ax.set_xticklabels(paired_df['P3'], fontsize=10)
ax.set_ylabel('ΔPSI (Pro-P2 minus Non-Pro-P2) ± SEM')
ax.set_xlabel('3rd Residue (P3)')
ax.set_title('Effect of Proline at P2 on Stability, by P3 Residue\n(red = Pro destabilizes, green = Pro stabilizes)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'Fig04_delta_psi_pro_effect_by_p3.png'))
plt.close()
print("  Saved Fig04")

# ══════════════════════════════════════════════════════════════════════════
# FIG 5: Stacked bar — Stable vs Unstable composition by P3
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("FIG 5: STACKED COMPOSITION (STABLE VS UNSTABLE)")
print("=" * 70)

stack_data = []
for p3 in sorted(df_pro['P3'].unique()):
    grp = df_pro[df_pro['P3'] == p3]
    n_total = len(grp)
    n_unstable = (grp['PSI'] <= THRESHOLD).sum()
    n_stable = n_total - n_unstable
    stack_data.append({'P3': p3, 'Unstable': n_unstable/n_total, 'Stable': n_stable/n_total, 'n': n_total})

stack_df = pd.DataFrame(stack_data).sort_values('Unstable', ascending=False).reset_index(drop=True)

fig, ax = plt.subplots(figsize=(11, 5))
x = np.arange(len(stack_df))

ax.bar(x, stack_df['Unstable'], color='#D62728', edgecolor='black', linewidth=0.5,
       label=f'Unstable (PSI ≤ {THRESHOLD:.2f})', width=0.75)
ax.bar(x, stack_df['Stable'], bottom=stack_df['Unstable'], color='#2CA02C',
       edgecolor='black', linewidth=0.5, label=f'Stable (PSI > {THRESHOLD:.2f})', width=0.75)

ax.axhline(lib_frac_unstable, color='black', linestyle='--', linewidth=1.2,
           label=f'Library fraction unstable ({lib_frac_unstable:.2f})')

ax.set_xticks(x)
ax.set_xticklabels([f"{r}\n({int(n)})" for r, n in zip(stack_df['P3'], stack_df['n'])], fontsize=9)
ax.set_ylabel('Fraction')
ax.set_xlabel('3rd Residue (P3)')
ax.set_title('Pro-P2 Peptides: Fraction Stable vs. Unstable by P3 Residue')
ax.legend(loc='upper right', framealpha=0.9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'Fig05_stacked_stable_unstable.png'))
plt.close()
print("  Saved Fig05")

# ══════════════════════════════════════════════════════════════════════════
# FIG 6: Enrichment/depletion of P3 residues in unstable Pro-P2
#   Log2(observed/expected) where expected = library P3 frequency
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("FIG 6: ENRICHMENT OF P3 RESIDUES IN UNSTABLE PRO-P2")
print("=" * 70)

lib_p3_freq = df['P3'].value_counts(normalize=True)
unstable_p3_freq = df_pro_unstable['P3'].value_counts(normalize=True)

enrich_data = []
for p3 in sorted(df_pro_unstable['P3'].unique()):
    obs = unstable_p3_freq.get(p3, 0)
    exp = lib_p3_freq.get(p3, 0.001)
    log2fc = np.log2(obs / exp) if obs > 0 and exp > 0 else 0
    enrich_data.append({'P3': p3, 'Observed_freq': obs, 'Expected_freq': exp, 'Log2FC': log2fc})

enrich_df = pd.DataFrame(enrich_data).sort_values('Log2FC').reset_index(drop=True)

print(f"{'P3':<5} {'Obs_freq':>9} {'Exp_freq':>9} {'Log2FC':>8}")
for _, r in enrich_df.iterrows():
    print(f"{r['P3']:<5} {r['Observed_freq']:>9.4f} {r['Expected_freq']:>9.4f} {r['Log2FC']:>8.3f}")

fig, ax = plt.subplots(figsize=(11, 5))
x = np.arange(len(enrich_df))
colors_enrich = ['#D62728' if v > 0 else '#2CA02C' for v in enrich_df['Log2FC']]

ax.bar(x, enrich_df['Log2FC'], color=colors_enrich, edgecolor='black',
       linewidth=0.6, width=0.7)
ax.axhline(0, color='black', linewidth=1.0)

ax.set_xticks(x)
ax.set_xticklabels(enrich_df['P3'], fontsize=10)
ax.set_ylabel('Log₂(Enrichment)')
ax.set_xlabel('3rd Residue (P3)')
ax.set_title('P3 Residue Enrichment in Unstable Pro-P2 Peptides\n(red = enriched among unstable, green = depleted)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'Fig06_enrichment_p3_in_unstable.png'))
plt.close()
print("  Saved Fig06")

# ══════════════════════════════════════════════════════════════════════════
# FIG 7: Heatmap — P2×P3 combination mean PSI (full library)
#   Shows where Pro at P2 sits relative to all other P2 residues
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("FIG 7: P2 x P3 HEATMAP (FULL LIBRARY)")
print("=" * 70)

top_p2 = df.groupby('P2').size().nlargest(10).index.tolist()
top_p3 = df.groupby('P3').size().nlargest(10).index.tolist()
if 'P' not in top_p2:
    top_p2.append('P')

pivot = df[df['P2'].isin(top_p2) & df['P3'].isin(top_p3)].groupby(['P2', 'P3'])['PSI'].mean().unstack(fill_value=np.nan)

p2_order = pivot.mean(axis=1).sort_values().index.tolist()
p3_order_heat = pivot.mean(axis=0).sort_values().index.tolist()
pivot = pivot.loc[p2_order, p3_order_heat]

fig, ax = plt.subplots(figsize=(9, 7))
sns.heatmap(pivot, annot=True, fmt='.2f', cmap='RdYlGn_r', linewidths=0.5,
            ax=ax, cbar_kws={'label': 'Mean PSI (AAVS1 KO)'},
            annot_kws={'size': 9})

for i, p2 in enumerate(pivot.index):
    if p2 == 'P':
        ax.add_patch(plt.Rectangle((0, i), len(pivot.columns), 1,
                                    fill=False, edgecolor='red', linewidth=2.5))
        break

ax.set_title('Mean PSI by P2 × P3 Combination\n(Proline row outlined in red)')
ax.set_ylabel('2nd Residue (P2)')
ax.set_xlabel('3rd Residue (P3)')
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'Fig07_heatmap_p2_x_p3.png'))
plt.close()
print("  Saved Fig07")

# ══════════════════════════════════════════════════════════════════════════
# FIG 8: Scatter — Pro-P2 mean PSI vs Non-Pro-P2 mean PSI per P3 residue
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("FIG 8: SCATTER — Pro vs Non-Pro PSI PER P3")
print("=" * 70)

fig, ax = plt.subplots(figsize=(6.5, 6))

lims = [2.5, 5.0]
ax.plot(lims, lims, '--', color='gray', linewidth=1, zorder=1, label='y = x (no Pro effect)')

for _, row in paired_df.iterrows():
    color = '#D62728' if row['p_value'] < 0.05 else '#7F7F7F'
    size = 60 if row['p_value'] < 0.05 else 40
    ax.scatter(row['NonPro_mean'], row['Pro_mean'], s=size, c=color,
               edgecolor='black', linewidth=0.6, zorder=3)
    offset_x = 0.03
    offset_y = 0.03
    ax.annotate(row['P3'], (row['NonPro_mean'], row['Pro_mean']),
                textcoords="offset points", xytext=(6, 4), fontsize=9)

ax.errorbar(paired_df['NonPro_mean'], paired_df['Pro_mean'],
            xerr=paired_df['NonPro_SEM'], yerr=paired_df['Pro_SEM'],
            fmt='none', ecolor='gray', alpha=0.4, linewidth=0.8, zorder=2)

ax.set_xlabel('Non-Pro-P2 Mean PSI ± SEM')
ax.set_ylabel('Pro-P2 Mean PSI ± SEM')
ax.set_title('Pro vs. Non-Pro Stability per P3 Residue')
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_aspect('equal')
ax.legend(loc='upper left', framealpha=0.9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.text(0.95, 0.05, 'Below diagonal:\nPro destabilizes',
        transform=ax.transAxes, ha='right', va='bottom', fontsize=8,
        style='italic', color='#D62728')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'Fig08_scatter_pro_vs_nonpro_per_p3.png'))
plt.close()
print("  Saved Fig08")

# ══════════════════════════════════════════════════════════════════════════
# FIG 9: Violin by P3 for unstable Pro-P2 only
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("FIG 9: VIOLIN — UNSTABLE PRO-P2 BY P3")
print("=" * 70)

top_p3_unstable = df_pro_unstable['P3'].value_counts().head(10).index.tolist()
df_plot = df_pro_unstable[df_pro_unstable['P3'].isin(top_p3_unstable)].copy()
p3_med_order = df_plot.groupby('P3')['PSI'].median().sort_values().index.tolist()

fig, ax = plt.subplots(figsize=(10, 5.5))
sns.boxplot(data=df_plot, x='P3', y='PSI', order=p3_med_order, ax=ax,
            hue='P3', palette='Reds_r', legend=False, fliersize=2, linewidth=0.8)

for i, p3 in enumerate(p3_med_order):
    grp = df_pro_unstable[df_pro_unstable['P3'] == p3]['PSI']
    ax.errorbar(i, grp.mean(), yerr=sem(grp), fmt='D', color='white',
                markeredgecolor='black', markersize=5, capsize=3, linewidth=1.0, zorder=5)

counts = df_pro_unstable.groupby('P3').size()
ax.set_xticklabels([f"{p3}\n(n={counts.get(p3, 0)})" for p3 in p3_med_order], fontsize=9)
ax.set_ylabel('PSI (AAVS1 KO)')
ax.set_xlabel('3rd Residue (P3)')
ax.set_title(f'Unstable Pro-P2 Peptides (PSI ≤ {THRESHOLD:.2f}): Distribution by P3')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'Fig09_boxplot_unstable_by_p3.png'))
plt.close()
print("  Saved Fig09")

# ══════════════════════════════════════════════════════════════════════════
# FIG 10: Cumulative fraction unstable — Pro-P2 vs Library by P3
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("FIG 10: ODDS RATIO — PRO-P2 INSTABILITY VS LIBRARY BY P3")
print("=" * 70)

or_data = []
for p3 in sorted(df_pro['P3'].unique()):
    pro_grp = df_pro[df_pro['P3'] == p3]
    lib_grp = df_non_pro[df_non_pro['P3'] == p3]
    if len(pro_grp) < 5 or len(lib_grp) < 5:
        continue

    a = (pro_grp['PSI'] <= THRESHOLD).sum()
    b = (pro_grp['PSI'] > THRESHOLD).sum()
    c = (lib_grp['PSI'] <= THRESHOLD).sum()
    d = (lib_grp['PSI'] > THRESHOLD).sum()

    if a == 0 or c == 0 or b == 0 or d == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5

    OR = (a * d) / (b * c)
    log_or = np.log2(OR)
    se_log_or = np.sqrt(1/a + 1/b + 1/c + 1/d)
    ci_low = log_or - 1.96 * se_log_or
    ci_high = log_or + 1.96 * se_log_or

    or_data.append({'P3': p3, 'OR': OR, 'Log2OR': log_or,
                    'CI_low': ci_low, 'CI_high': ci_high, 'n_pro': len(pro_grp)})

or_df = pd.DataFrame(or_data).sort_values('Log2OR').reset_index(drop=True)

fig, ax = plt.subplots(figsize=(11, 5.5))
x = np.arange(len(or_df))
colors_or = ['#D62728' if v > 0 else '#2CA02C' for v in or_df['Log2OR']]
yerr_low = or_df['Log2OR'] - or_df['CI_low']
yerr_high = or_df['CI_high'] - or_df['Log2OR']

ax.bar(x, or_df['Log2OR'], yerr=[yerr_low, yerr_high], capsize=3,
       color=colors_or, edgecolor='black', linewidth=0.6, width=0.7,
       error_kw={'linewidth': 1.0})
ax.axhline(0, color='black', linewidth=1.0)

ax.set_xticks(x)
ax.set_xticklabels([f"{r}\n({int(n)})" for r, n in zip(or_df['P3'], or_df['n_pro'])], fontsize=9)
ax.set_ylabel('Log₂(Odds Ratio) ± 95% CI')
ax.set_xlabel('3rd Residue (P3)')
ax.set_title('Odds of Instability: Pro-P2 vs. Non-Pro-P2, by P3 Residue\n(red = Pro more likely unstable, green = Pro less likely unstable)')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'Fig10_odds_ratio_by_p3.png'))
plt.close()
print("  Saved Fig10")

# ══════════════════════════════════════════════════════════════════════════
# EXPORT TO EXCEL
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("EXPORTING")
print("=" * 70)

output_excel = os.path.join(OUTPUT_DIR, 'Unstable_Pro_P2_analysis.xlsx')
with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
    df_pro_unstable[['Transcript', 'Gene', 'Sequence', 'P3', 'PSI']].sort_values('PSI').to_excel(
        writer, sheet_name='Unstable_Pro_P2_peptides', index=False)

    unstable_p3.to_excel(writer, sheet_name='Unstable_P3_summary', index=False)
    frac_df.to_excel(writer, sheet_name='Fraction_unstable_by_P3', index=False)
    paired_df.to_excel(writer, sheet_name='Pro_vs_NonPro_by_P3', index=False)
    enrich_df.to_excel(writer, sheet_name='P3_enrichment', index=False)
    or_df.to_excel(writer, sheet_name='Odds_ratio_by_P3', index=False)

    methods = pd.DataFrame({
        'Step': [
            '1. Threshold definition',
            '2. Unstable Pro-P2 by P3',
            '3. Fraction unstable',
            '4. Paired P3 comparison',
            '5. Delta PSI',
            '6. Enrichment analysis',
            '7. P2 x P3 heatmap',
            '8. Scatter plot',
            '9. Odds ratio',
        ],
        'Description': [
            f'Defined unstable peptides as PSI <= library Q1 ({THRESHOLD:.4f}). {len(df_pro_unstable)} of {len(df_pro)} Pro-P2 peptides ({len(df_pro_unstable)/len(df_pro)*100:.1f}%) fall below this threshold.',
            'Grouped unstable Pro-P2 peptides by their 3rd residue and calculated mean PSI ± SEM.',
            'For each P3 residue, calculated the fraction of Pro-P2 peptides classified as unstable, with 95% binomial CI.',
            'For each P3 residue, compared Pro-P2 vs Non-Pro-P2 mean PSI using Welch\'s t-test. This controls for P3 identity and isolates the Pro effect.',
            'Calculated Delta PSI = Pro-P2 mean - Non-Pro-P2 mean for each P3, with propagated SEM. Negative = Pro destabilizes.',
            'Calculated log2(observed/expected) frequency of each P3 residue among unstable Pro-P2 peptides vs library baseline.',
            'Created a P2 x P3 mean PSI heatmap showing combination effects across the full library.',
            'Scatter plot of Pro-P2 vs Non-Pro-P2 mean PSI per P3 residue. Points below diagonal = Pro is destabilizing.',
            'Calculated odds ratio of instability for Pro-P2 vs Non-Pro-P2 at each P3, with 95% CI on log2 scale.',
        ]
    })
    methods.to_excel(writer, sheet_name='Methods', index=False)

print(f"  Saved: {output_excel}")

# ══════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

n_sig = sum(1 for _, r in paired_df.iterrows() if r['p_value'] < 0.05)
sig_residues = [r['P3'] for _, r in paired_df.iterrows() if r['p_value'] < 0.05]

print(f"""
UNSTABLE THRESHOLD: PSI <= {THRESHOLD:.4f} (Library Q1)
UNSTABLE PRO-P2: {len(df_pro_unstable)}/{len(df_pro)} ({len(df_pro_unstable)/len(df_pro)*100:.1f}%)

P3 RESIDUES MOST ENRICHED AMONG UNSTABLE PRO-P2:
  {enrich_df.iloc[-1]['P3']} (Log2FC = {enrich_df.iloc[-1]['Log2FC']:+.3f})
  {enrich_df.iloc[-2]['P3']} (Log2FC = {enrich_df.iloc[-2]['Log2FC']:+.3f})

P3 RESIDUES MOST DEPLETED AMONG UNSTABLE PRO-P2:
  {enrich_df.iloc[0]['P3']} (Log2FC = {enrich_df.iloc[0]['Log2FC']:+.3f})
  {enrich_df.iloc[1]['P3']} (Log2FC = {enrich_df.iloc[1]['Log2FC']:+.3f})

PAIRED COMPARISON (Pro vs Non-Pro, same P3):
  Significant at p<0.05 for {n_sig}/20 P3 residues: {sig_residues}

HIGHEST ODDS OF INSTABILITY (Pro vs Non-Pro):
  {or_df.iloc[-1]['P3']} (OR = {or_df.iloc[-1]['OR']:.2f}, Log2OR = {or_df.iloc[-1]['Log2OR']:+.3f})

LOWEST ODDS OF INSTABILITY (Pro vs Non-Pro):
  {or_df.iloc[0]['P3']} (OR = {or_df.iloc[0]['OR']:.2f}, Log2OR = {or_df.iloc[0]['Log2OR']:+.3f})
""")
print("Done!")
