import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import os, sys, io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

OUTPUT_DIR = r'c:\Users\User\Desktop\תינוקת\ubr3paper\pro_gly_analysis'
os.makedirs(OUTPUT_DIR, exist_ok=True)

SCREEN_EXCEL = r'c:\Users\User\Downloads\UBR3 Nt screen (1).xlsx'

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

# ══════════════════════════════════════════════════════════════════════════
# LOAD DATA
# ══════════════════════════════════════════════════════════════════════════

print("=" * 70)
print("LOADING DATA")
print("=" * 70)

df_all = pd.read_excel(SCREEN_EXCEL, sheet_name='Nprot_5_analyzed')
df_hits = pd.read_excel(SCREEN_EXCEL, sheet_name='sub_high')

df_all['PSI_AAVS'] = (df_all['PSI-293a'] + df_all['PSI-293b']) / 2
df_hits['PSI_AAVS'] = (df_hits['PSI-293a'] + df_hits['PSI-293b']) / 2

print(f"All screen peptides: {len(df_all)}")
print(f"Top hits:            {len(df_hits)}")

# ══════════════════════════════════════════════════════════════════════════
# SEPARATE BY P2 RESIDUE
# ══════════════════════════════════════════════════════════════════════════

hits_pro = df_hits[df_hits['AA2'] == 'P'].copy()
hits_gly = df_hits[df_hits['AA2'] == 'G'].copy()
all_pro  = df_all[df_all['AA2'] == 'P'].copy()
all_gly  = df_all[df_all['AA2'] == 'G'].copy()

print(f"\nTop hits  — Pro at P2: {len(hits_pro)},  Gly at P2: {len(hits_gly)}")
print(f"All screen — Pro at P2: {len(all_pro)},  Gly at P2: {len(all_gly)}")

# ══════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════

def significance_marker(p):
    if p < 0.0001: return '****'
    if p < 0.001:  return '***'
    if p < 0.01:   return '**'
    if p < 0.05:   return '*'
    return 'ns'

def draw_bracket(ax, x1, x2, y, h, text, fontsize=11):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.2, color='black')
    ax.text((x1 + x2) / 2, y + h, text, ha='center', va='bottom',
            fontsize=fontsize, fontweight='bold')


def make_boxplot(ax, data_list, labels, colors, title, ylabel='Protein Stability Index'):
    bp = ax.boxplot(data_list, positions=range(len(data_list)), widths=0.55,
                    patch_artist=True, showfliers=True,
                    flierprops={'marker': 'o', 'markersize': 2.5, 'alpha': 0.4,
                                'markerfacecolor': '#999', 'markeredgecolor': '#999'},
                    whiskerprops={'linewidth': 1.2, 'color': 'black'},
                    capprops={'linewidth': 1.2, 'color': 'black'},
                    medianprops={'linewidth': 1.8, 'color': 'black'})

    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor('black')
        patch.set_linewidth(1.0)
        patch.set_alpha(0.85)

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=35, ha='right', fontsize=10)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=13, fontweight='bold', pad=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    return bp

# ══════════════════════════════════════════════════════════════════════════
# STATISTICS
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("STATISTICAL ANALYSIS")
print("=" * 70)

def run_stats(group1, group2, label1, label2):
    g1 = group1.dropna().values
    g2 = group2.dropna().values
    t_stat, t_p = stats.ttest_ind(g1, g2, equal_var=False)
    mw_stat, mw_p = stats.mannwhitneyu(g1, g2, alternative='two-sided')
    print(f"  {label1} (n={len(g1)}, mean={np.mean(g1):.4f}) vs "
          f"{label2} (n={len(g2)}, mean={np.mean(g2):.4f})")
    print(f"    Welch t: t={t_stat:.3f}, p={t_p:.2e} ({significance_marker(t_p)})")
    print(f"    Mann-Whitney: U={mw_stat:.0f}, p={mw_p:.2e} ({significance_marker(mw_p)})")
    return mw_p

print("\n--- PROLINE (P) at Position 2 ---")
p_hits_vs_all = run_stats(hits_pro['PSI_AAVS'], all_pro['PSI_AAVS'],
                          'Top Hits Pro-P2', 'All Screen Pro-P2')

print("\n--- GLYCINE (G) at Position 2 ---")
p_gly_hits_vs_all = run_stats(hits_gly['PSI_AAVS'], all_gly['PSI_AAVS'],
                              'Top Hits Gly-P2', 'All Screen Gly-P2')

print("\n--- CROSS COMPARISONS ---")
p_all_pro_vs_gly = run_stats(all_pro['PSI_AAVS'], all_gly['PSI_AAVS'],
                             'All Screen Pro-P2', 'All Screen Gly-P2')
p_hits_pro_vs_gly = run_stats(hits_pro['PSI_AAVS'], hits_gly['PSI_AAVS'],
                              'Top Hits Pro-P2', 'Top Hits Gly-P2')

# ══════════════════════════════════════════════════════════════════════════
# FIGURE 1: PROLINE ANALYSIS — Box Plot (Top Hits vs All Screen)
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("GENERATING FIGURES")
print("=" * 70)

fig, ax = plt.subplots(figsize=(5, 5.5))

pro_data = [
    all_pro['PSI_AAVS'].dropna().values,
    hits_pro['PSI_AAVS'].dropna().values,
]
pro_labels = [
    f'All Screen\nPro-P2\n(n={len(pro_data[0])})',
    f'Top Hits\nPro-P2\n(n={len(pro_data[1])})',
]
pro_colors = ['#B0B0B0', '#E74C3C']

bp = make_boxplot(ax, pro_data, pro_labels, pro_colors,
                  'Proline at Position 2: Top Hits vs All Screen')

y_max = max(np.percentile(d, 99) for d in pro_data)
bracket_y = y_max + 0.15
draw_bracket(ax, 0, 1, bracket_y, 0.08, significance_marker(p_hits_vs_all))
ax.set_ylim(top=bracket_y + 0.35)

plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, 'Fig1_Pro_P2_boxplot.png'))
fig.savefig(os.path.join(OUTPUT_DIR, 'Fig1_Pro_P2_boxplot.pdf'))
plt.close()
print("  Saved Fig1_Pro_P2_boxplot")

# ══════════════════════════════════════════════════════════════════════════
# FIGURE 2: GLYCINE ANALYSIS — Box Plot (Top Hits vs All Screen)
# ══════════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(5, 5.5))

gly_data = [
    all_gly['PSI_AAVS'].dropna().values,
    hits_gly['PSI_AAVS'].dropna().values,
]
gly_labels = [
    f'All Screen\nGly-P2\n(n={len(gly_data[0])})',
    f'Top Hits\nGly-P2\n(n={len(gly_data[1])})',
]
gly_colors = ['#B0B0B0', '#2E86C1']

bp = make_boxplot(ax, gly_data, gly_labels, gly_colors,
                  'Glycine at Position 2: Top Hits vs All Screen')

y_max = max(np.percentile(d, 99) for d in gly_data)
bracket_y = y_max + 0.15
draw_bracket(ax, 0, 1, bracket_y, 0.08, significance_marker(p_gly_hits_vs_all))
ax.set_ylim(top=bracket_y + 0.35)

plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, 'Fig2_Gly_P2_boxplot.png'))
fig.savefig(os.path.join(OUTPUT_DIR, 'Fig2_Gly_P2_boxplot.pdf'))
plt.close()
print("  Saved Fig2_Gly_P2_boxplot")

# ══════════════════════════════════════════════════════════════════════════
# FIGURE 3: COMBINED — All 4 groups in one plot (reference image style)
# ══════════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(7.5, 5.5))

combined_data = [
    all_pro['PSI_AAVS'].dropna().values,
    hits_pro['PSI_AAVS'].dropna().values,
    all_gly['PSI_AAVS'].dropna().values,
    hits_gly['PSI_AAVS'].dropna().values,
]

combined_labels = [
    f'All Screen\n(n={len(combined_data[0])})',
    f'Top Hits\n(n={len(combined_data[1])})',
    f'All Screen\n(n={len(combined_data[2])})',
    f'Top Hits\n(n={len(combined_data[3])})',
]
combined_colors = ['#B0B0B0', '#E74C3C', '#B0B0B0', '#E67E22']

bp = make_boxplot(ax, combined_data, combined_labels, combined_colors,
                  'Protein Stability Index: Pro vs Gly at Position 2',
                  ylabel='Protein Stability Index')
ax.set_xticklabels(combined_labels, rotation=0, ha='center', fontsize=9)

# Group labels below x-axis
ax.text(0.5, -0.18, 'Pro at P2', transform=ax.get_xaxis_transform(),
        ha='center', fontsize=12, fontweight='bold')
ax.text(2.5, -0.18, 'Gly at P2', transform=ax.get_xaxis_transform(),
        ha='center', fontsize=12, fontweight='bold')

# Dividing lines under group labels
ax.annotate('', xy=(0, -0.14), xytext=(1, -0.14),
            xycoords=('data', 'axes fraction'), textcoords=('data', 'axes fraction'),
            arrowprops=dict(arrowstyle='-', lw=1.5, color='black'))
ax.annotate('', xy=(2, -0.14), xytext=(3, -0.14),
            xycoords=('data', 'axes fraction'), textcoords=('data', 'axes fraction'),
            arrowprops=dict(arrowstyle='-', lw=1.5, color='black'))

y_max_all = max(np.percentile(d, 99) for d in combined_data)

# Within-group brackets
bracket_base = y_max_all + 0.12
draw_bracket(ax, 0, 1, bracket_base, 0.06, significance_marker(p_hits_vs_all), fontsize=10)
draw_bracket(ax, 2, 3, bracket_base, 0.06, significance_marker(p_gly_hits_vs_all), fontsize=10)

# Cross-group bracket (All Screen Pro vs All Screen Gly)
cross_y = bracket_base + 0.25
draw_bracket(ax, 0, 2, cross_y, 0.06, significance_marker(p_all_pro_vs_gly), fontsize=10)

# Cross-group bracket (Top Hits Pro vs Top Hits Gly)
cross_y2 = cross_y + 0.22
draw_bracket(ax, 1, 3, cross_y2, 0.06, significance_marker(p_hits_pro_vs_gly), fontsize=10)

ax.set_ylim(top=cross_y2 + 0.4)

plt.tight_layout()
fig.subplots_adjust(bottom=0.22)
fig.savefig(os.path.join(OUTPUT_DIR, 'Fig3_Combined_Pro_Gly_boxplot.png'))
fig.savefig(os.path.join(OUTPUT_DIR, 'Fig3_Combined_Pro_Gly_boxplot.pdf'))
plt.close()
print("  Saved Fig3_Combined_Pro_Gly_boxplot")

# ══════════════════════════════════════════════════════════════════════════
# FIGURE 4: COMBINED — Cleaner version matching reference style exactly
# ══════════════════════════════════════════════════════════════════════════

fig, ax = plt.subplots(figsize=(6.5, 5.5))

positions = [0, 1, 2.5, 3.5]
bp = ax.boxplot(combined_data, positions=positions, widths=0.6,
                patch_artist=True, showfliers=True,
                flierprops={'marker': 'o', 'markersize': 2, 'alpha': 0.3,
                            'markerfacecolor': '#999', 'markeredgecolor': '#999'},
                whiskerprops={'linewidth': 1.2, 'color': 'black'},
                capprops={'linewidth': 1.2, 'color': 'black'},
                medianprops={'linewidth': 1.8, 'color': 'black'})

for patch, color in zip(bp['boxes'], combined_colors):
    patch.set_facecolor(color)
    patch.set_edgecolor('black')
    patch.set_linewidth(1.0)
    patch.set_alpha(0.85)

ax.set_xticks(positions)
ax.set_xticklabels([
    f'All Screen\n(n={len(combined_data[0])})',
    f'Top Hits\n(n={len(combined_data[1])})',
    f'All Screen\n(n={len(combined_data[2])})',
    f'Top Hits\n(n={len(combined_data[3])})',
], fontsize=9)

ax.set_ylabel('Protein Stability Index', fontsize=12)
ax.set_title('PSI: Pro vs Gly at Position 2\nTop Hits vs All Screen', fontsize=13, fontweight='bold', pad=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Group labels
ax.text(0.5, -0.17, 'Pro at P2', transform=ax.get_xaxis_transform(),
        ha='center', fontsize=11, fontweight='bold')
ax.text(3.0, -0.17, 'Gly at P2', transform=ax.get_xaxis_transform(),
        ha='center', fontsize=11, fontweight='bold')

# Underlines for group labels
ax.plot([0, 1], [-0.13, -0.13], transform=ax.get_xaxis_transform(),
        color='black', linewidth=1.5, clip_on=False)
ax.plot([2.5, 3.5], [-0.13, -0.13], transform=ax.get_xaxis_transform(),
        color='black', linewidth=1.5, clip_on=False)

y_max_all = max(np.percentile(d, 99) for d in combined_data)

# Within-group significance brackets
bracket_h = 0.06
b1 = y_max_all + 0.15
draw_bracket(ax, positions[0], positions[1], b1, bracket_h,
             significance_marker(p_hits_vs_all), fontsize=10)
draw_bracket(ax, positions[2], positions[3], b1, bracket_h,
             significance_marker(p_gly_hits_vs_all), fontsize=10)

# Across-group significance bracket (top hits vs top hits)
b2 = b1 + 0.28
draw_bracket(ax, positions[1], positions[2], b2, bracket_h,
             significance_marker(p_hits_pro_vs_gly), fontsize=10)

# Across-group bracket (all screen vs all screen)
b3 = b2 + 0.25
draw_bracket(ax, positions[0], positions[3], b3, bracket_h,
             significance_marker(p_all_pro_vs_gly), fontsize=9)

ax.set_ylim(top=b3 + 0.40)

plt.tight_layout()
fig.subplots_adjust(bottom=0.22)
fig.savefig(os.path.join(OUTPUT_DIR, 'Fig4_Combined_Pro_Gly_publication.png'))
fig.savefig(os.path.join(OUTPUT_DIR, 'Fig4_Combined_Pro_Gly_publication.pdf'))
plt.close()
print("  Saved Fig4_Combined_Pro_Gly_publication")

# ══════════════════════════════════════════════════════════════════════════
# SUMMARY TABLE
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

summary_rows = []
for name, data in [('All Screen Pro-P2', all_pro['PSI_AAVS']),
                    ('Top Hits Pro-P2', hits_pro['PSI_AAVS']),
                    ('All Screen Gly-P2', all_gly['PSI_AAVS']),
                    ('Top Hits Gly-P2', hits_gly['PSI_AAVS'])]:
    d = data.dropna()
    summary_rows.append({
        'Group': name, 'n': len(d),
        'Mean': d.mean(), 'Median': d.median(),
        'SD': d.std(), 'SEM': d.std() / np.sqrt(len(d)),
        'Q1': d.quantile(0.25), 'Q3': d.quantile(0.75),
    })

summary_df = pd.DataFrame(summary_rows)
print(summary_df.to_string(index=False, float_format='{:.4f}'.format))

output_excel = os.path.join(OUTPUT_DIR, 'Pro_Gly_P2_analysis.xlsx')
with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
    summary_df.to_excel(writer, sheet_name='Summary', index=False)

    hits_pro[['Gene_ID', 'AA_seq', 'AA2', 'AA3', 'PSI_AAVS']].to_excel(
        writer, sheet_name='TopHits_Pro_P2', index=False)
    hits_gly[['Gene_ID', 'AA_seq', 'AA2', 'AA3', 'PSI_AAVS']].to_excel(
        writer, sheet_name='TopHits_Gly_P2', index=False)

    stat_df = pd.DataFrame([
        {'Comparison': 'Top Hits Pro-P2 vs All Screen Pro-P2',
         'p_value': p_hits_vs_all, 'Significance': significance_marker(p_hits_vs_all)},
        {'Comparison': 'Top Hits Gly-P2 vs All Screen Gly-P2',
         'p_value': p_gly_hits_vs_all, 'Significance': significance_marker(p_gly_hits_vs_all)},
        {'Comparison': 'All Screen Pro-P2 vs All Screen Gly-P2',
         'p_value': p_all_pro_vs_gly, 'Significance': significance_marker(p_all_pro_vs_gly)},
        {'Comparison': 'Top Hits Pro-P2 vs Top Hits Gly-P2',
         'p_value': p_hits_pro_vs_gly, 'Significance': significance_marker(p_hits_pro_vs_gly)},
    ])
    stat_df.to_excel(writer, sheet_name='Statistics', index=False)

print(f"\n  Excel saved: {output_excel}")
print(f"  Figures saved to: {OUTPUT_DIR}")
print("\nDone!")
