import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import os, sys, io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

OUTPUT_DIR = r'c:\Users\User\Desktop\תינוקת\ubr3paper\dipeptide_analysis'
os.makedirs(OUTPUT_DIR, exist_ok=True)

SCREEN_EXCEL = r'c:\Users\User\Downloads\UBR3 Nt screen (1).xlsx'

DIPEPTIDES = ['PD', 'PE', 'PT', 'GE', 'GD']

DIPEPTIDE_COLORS = {
    'PD': '#E74C3C',
    'PE': '#E67E22',
    'PT': '#F1C40F',
    'GE': '#2ECC71',
    'GD': '#2E86C1',
}

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

df_all['dipeptide'] = df_all['AA2'].astype(str) + df_all['AA3'].astype(str)
df_hits['dipeptide'] = df_hits['AA2'].astype(str) + df_hits['AA3'].astype(str)

print(f"All screen peptides: {len(df_all)}")
print(f"Top hits:            {len(df_hits)}")

# ══════════════════════════════════════════════════════════════════════════
# SEPARATE BY DIPEPTIDE
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("DIPEPTIDE COUNTS")
print("=" * 70)

dp_data = {}
for dp in DIPEPTIDES:
    h = df_hits[df_hits['dipeptide'] == dp].copy()
    a = df_all[df_all['dipeptide'] == dp].copy()
    dp_data[dp] = {'hits': h, 'all': a}
    print(f"  {dp}: Top Hits = {len(h)},  All Screen = {len(a)}")

# ══════════════════════════════════════════════════════════════════════════
# HELPERS
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

# ══════════════════════════════════════════════════════════════════════════
# STATISTICS
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("STATISTICAL ANALYSIS")
print("=" * 70)

pvals = {}
for dp in DIPEPTIDES:
    hits_psi = dp_data[dp]['hits']['PSI_AAVS'].dropna().values
    all_psi  = dp_data[dp]['all']['PSI_AAVS'].dropna().values
    print(f"\n--- {dp} ---")
    print(f"  Top Hits:   n={len(hits_psi)}, mean={np.mean(hits_psi):.4f}" if len(hits_psi) > 0 else f"  Top Hits:   n=0")
    print(f"  All Screen: n={len(all_psi)},  mean={np.mean(all_psi):.4f}")

    if len(hits_psi) >= 2 and len(all_psi) >= 2:
        mw_stat, mw_p = stats.mannwhitneyu(hits_psi, all_psi, alternative='two-sided')
        t_stat, t_p = stats.ttest_ind(hits_psi, all_psi, equal_var=False)
        print(f"  Welch t: t={t_stat:.3f}, p={t_p:.2e} ({significance_marker(t_p)})")
        print(f"  Mann-Whitney: U={mw_stat:.0f}, p={mw_p:.2e} ({significance_marker(mw_p)})")
        pvals[dp] = mw_p
    elif len(hits_psi) == 1:
        print(f"  Only 1 hit — no test possible (single value: {hits_psi[0]:.4f})")
        pvals[dp] = np.nan
    else:
        print(f"  No hits — skipping")
        pvals[dp] = np.nan

# ══════════════════════════════════════════════════════════════════════════
# INDIVIDUAL FIGURES — one per dipeptide
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("GENERATING INDIVIDUAL FIGURES")
print("=" * 70)

for dp in DIPEPTIDES:
    hits_psi = dp_data[dp]['hits']['PSI_AAVS'].dropna().values
    all_psi  = dp_data[dp]['all']['PSI_AAVS'].dropna().values

    if len(hits_psi) == 0:
        print(f"  Skipping {dp} — no top hits")
        continue

    fig, ax = plt.subplots(figsize=(5, 5.5))

    data = [all_psi, hits_psi]
    colors = ['#B0B0B0', DIPEPTIDE_COLORS[dp]]
    labels = [
        f'All Screen\n{dp}\n(n={len(all_psi)})',
        f'Top Hits\n{dp}\n(n={len(hits_psi)})',
    ]

    bp = ax.boxplot(data, positions=[0, 1], widths=0.55,
                    patch_artist=True, showfliers=True,
                    flierprops={'marker': 'o', 'markersize': 2.5, 'alpha': 0.4,
                                'markerfacecolor': '#999', 'markeredgecolor': '#999'},
                    whiskerprops={'linewidth': 1.2, 'color': 'black'},
                    capprops={'linewidth': 1.2, 'color': 'black'},
                    medianprops={'linewidth': 1.8, 'color': 'black'})

    for patch, c in zip(bp['boxes'], colors):
        patch.set_facecolor(c)
        patch.set_edgecolor('black')
        patch.set_linewidth(1.0)
        patch.set_alpha(0.85)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel('Protein Stability Index', fontsize=12)
    aa2, aa3 = dp[0], dp[1]
    aa2_name = 'Pro' if aa2 == 'P' else 'Gly'
    ax.set_title(f'{dp} ({aa2_name}-{aa3}): Top Hits vs All Screen',
                 fontsize=13, fontweight='bold', pad=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if not np.isnan(pvals.get(dp, np.nan)):
        y_max = max(np.percentile(all_psi, 99), np.max(hits_psi) if len(hits_psi) > 0 else 0)
        bracket_y = y_max + 0.15
        draw_bracket(ax, 0, 1, bracket_y, 0.08, significance_marker(pvals[dp]))
        ax.set_ylim(top=bracket_y + 0.40)

    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, f'Fig_{dp}_boxplot.png'))
    fig.savefig(os.path.join(OUTPUT_DIR, f'Fig_{dp}_boxplot.pdf'))
    plt.close()
    print(f"  Saved Fig_{dp}_boxplot")

# ══════════════════════════════════════════════════════════════════════════
# COMBINED FIGURE — all dipeptides side by side (publication style)
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("GENERATING COMBINED FIGURE")
print("=" * 70)

valid_dp = [dp for dp in DIPEPTIDES if len(dp_data[dp]['hits']) > 0]
n_dp = len(valid_dp)

fig, ax = plt.subplots(figsize=(3.2 * n_dp, 5.5))

positions = []
all_data = []
all_colors = []
tick_positions = []
tick_labels = []
group_centers = []
group_labels = []
bracket_pairs = []

pos = 0
for i, dp in enumerate(valid_dp):
    hits_psi = dp_data[dp]['hits']['PSI_AAVS'].dropna().values
    all_psi  = dp_data[dp]['all']['PSI_AAVS'].dropna().values

    p_all = pos
    p_hits = pos + 1

    positions.extend([p_all, p_hits])
    all_data.extend([all_psi, hits_psi])
    all_colors.extend(['#B0B0B0', DIPEPTIDE_COLORS[dp]])

    tick_positions.extend([p_all, p_hits])
    tick_labels.extend([
        f'All Screen\n(n={len(all_psi)})',
        f'Top Hits\n(n={len(hits_psi)})',
    ])

    group_centers.append((p_all + p_hits) / 2)
    group_labels.append(dp)
    bracket_pairs.append((p_all, p_hits, dp))

    pos += 3

bp = ax.boxplot(all_data, positions=positions, widths=0.6,
                patch_artist=True, showfliers=True,
                flierprops={'marker': 'o', 'markersize': 2, 'alpha': 0.3,
                            'markerfacecolor': '#999', 'markeredgecolor': '#999'},
                whiskerprops={'linewidth': 1.2, 'color': 'black'},
                capprops={'linewidth': 1.2, 'color': 'black'},
                medianprops={'linewidth': 1.8, 'color': 'black'})

for patch, c in zip(bp['boxes'], all_colors):
    patch.set_facecolor(c)
    patch.set_edgecolor('black')
    patch.set_linewidth(1.0)
    patch.set_alpha(0.85)

ax.set_xticks(tick_positions)
ax.set_xticklabels(tick_labels, fontsize=8)
ax.set_ylabel('Protein Stability Index', fontsize=12)
ax.set_title('Dipeptide Motifs (P2-P3): Top Hits vs All Screen\nPSI Comparison',
             fontsize=13, fontweight='bold', pad=10)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Group labels under each pair
for center, label in zip(group_centers, group_labels):
    ax.text(center, -0.16, label, transform=ax.get_xaxis_transform(),
            ha='center', fontsize=12, fontweight='bold',
            color=DIPEPTIDE_COLORS[label])

# Underlines
for p_all, p_hits, dp in bracket_pairs:
    ax.plot([p_all - 0.3, p_hits + 0.3], [-0.12, -0.12],
            transform=ax.get_xaxis_transform(),
            color=DIPEPTIDE_COLORS[dp], linewidth=2.0, clip_on=False)

# Significance brackets
global_max = max(np.percentile(d, 99) for d in all_data if len(d) > 0)
bracket_y = global_max + 0.15

for p_all, p_hits, dp in bracket_pairs:
    p = pvals.get(dp, np.nan)
    if not np.isnan(p):
        draw_bracket(ax, p_all, p_hits, bracket_y, 0.06,
                     significance_marker(p), fontsize=9)

ax.set_ylim(top=bracket_y + 0.45)

plt.tight_layout()
fig.subplots_adjust(bottom=0.22)
fig.savefig(os.path.join(OUTPUT_DIR, 'Fig_Combined_dipeptides_boxplot.png'))
fig.savefig(os.path.join(OUTPUT_DIR, 'Fig_Combined_dipeptides_boxplot.pdf'))
plt.close()
print("  Saved Fig_Combined_dipeptides_boxplot")

# ══════════════════════════════════════════════════════════════════════════
# SUMMARY TABLE & EXCEL EXPORT
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

summary_rows = []
for dp in DIPEPTIDES:
    for label, subset in [('All Screen', dp_data[dp]['all']),
                           ('Top Hits',  dp_data[dp]['hits'])]:
        d = subset['PSI_AAVS'].dropna()
        summary_rows.append({
            'Dipeptide': dp, 'Group': label, 'n': len(d),
            'Mean': d.mean() if len(d) > 0 else np.nan,
            'Median': d.median() if len(d) > 0 else np.nan,
            'SD': d.std() if len(d) > 1 else np.nan,
            'SEM': d.std() / np.sqrt(len(d)) if len(d) > 1 else np.nan,
            'Q1': d.quantile(0.25) if len(d) > 0 else np.nan,
            'Q3': d.quantile(0.75) if len(d) > 0 else np.nan,
        })

summary_df = pd.DataFrame(summary_rows)
print(summary_df.to_string(index=False, float_format='{:.4f}'.format))

stat_rows = []
for dp in DIPEPTIDES:
    p = pvals.get(dp, np.nan)
    stat_rows.append({
        'Dipeptide': dp,
        'n_hits': len(dp_data[dp]['hits']),
        'n_all': len(dp_data[dp]['all']),
        'Hits_mean_PSI': dp_data[dp]['hits']['PSI_AAVS'].mean() if len(dp_data[dp]['hits']) > 0 else np.nan,
        'All_mean_PSI': dp_data[dp]['all']['PSI_AAVS'].mean(),
        'Delta_PSI': (dp_data[dp]['hits']['PSI_AAVS'].mean() - dp_data[dp]['all']['PSI_AAVS'].mean()) if len(dp_data[dp]['hits']) > 0 else np.nan,
        'MannWhitney_p': p,
        'Significance': significance_marker(p) if not np.isnan(p) else 'N/A',
    })

stat_df = pd.DataFrame(stat_rows)
print("\n" + stat_df.to_string(index=False, float_format='{:.4f}'.format))

output_excel = os.path.join(OUTPUT_DIR, 'Dipeptide_P2P3_analysis.xlsx')
with pd.ExcelWriter(output_excel, engine='xlsxwriter') as writer:
    summary_df.to_excel(writer, sheet_name='Summary', index=False)
    stat_df.to_excel(writer, sheet_name='Statistics', index=False)

    for dp in DIPEPTIDES:
        h = dp_data[dp]['hits']
        if len(h) > 0:
            h[['Gene_ID', 'AA_seq', 'AA2', 'AA3', 'dipeptide', 'PSI_AAVS']].to_excel(
                writer, sheet_name=f'TopHits_{dp}', index=False)

print(f"\n  Excel saved: {output_excel}")
print(f"  Figures saved to: {OUTPUT_DIR}")
print("\nDone!")
