import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import os, sys, io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

OUTPUT_DIR = r'c:\Users\User\Desktop\תינוקת\dianaMega\ubr3enrichmentlogo\dipeptide_analysis'
os.makedirs(OUTPUT_DIR, exist_ok=True)

SCREEN_EXCEL = r'c:\Users\User\Downloads\UBR3 Nt screen (1).xlsx'

DIPEPTIDES = ['PD', 'PE', 'PT', 'GE', 'GD']
COMBINED_GROUPS = {
    'PD+PE+PT': ['PD', 'PE', 'PT'],
    'GE+GD':    ['GE', 'GD'],
}
ALL_KEYS = DIPEPTIDES + list(COMBINED_GROUPS.keys())

DIPEPTIDE_COLORS = {
    'PD': '#E74C3C',
    'PE': '#E67E22',
    'PT': '#F1C40F',
    'GE': '#2ECC71',
    'GD': '#2E86C1',
    'PD+PE+PT': '#E67E22',
    'GE+GD':    '#DC143C',
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

for grp_name, members in COMBINED_GROUPS.items():
    h = pd.concat([df_hits[df_hits['dipeptide'] == m] for m in members], ignore_index=True)
    a = pd.concat([df_all[df_all['dipeptide'] == m] for m in members], ignore_index=True)
    dp_data[grp_name] = {'hits': h, 'all': a}
    print(f"  {grp_name}: Top Hits = {len(h)},  All Screen = {len(a)}")

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
for dp in ALL_KEYS:
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

for dp in ALL_KEYS:
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
    ax.set_title(f'{dp}: Top Hits vs All Screen',
                 fontsize=13, fontweight='bold', pad=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if not np.isnan(pvals.get(dp, np.nan)):
        y_max = max(np.percentile(all_psi, 99), np.max(hits_psi) if len(hits_psi) > 0 else 0)
        bracket_y = y_max + 0.15
        draw_bracket(ax, 0, 1, bracket_y, 0.08, significance_marker(pvals[dp]))
        ax.set_ylim(top=bracket_y + 0.40)

    plt.tight_layout()
    safe_name = dp.replace('+', '_')
    fig.savefig(os.path.join(OUTPUT_DIR, f'Fig_{safe_name}_boxplot.png'))
    fig.savefig(os.path.join(OUTPUT_DIR, f'Fig_{safe_name}_boxplot.pdf'))
    plt.close()
    print(f"  Saved Fig_{safe_name}_boxplot")

# ══════════════════════════════════════════════════════════════════════════
# COMBINED FIGURE — all dipeptides side by side (publication style)
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("GENERATING COMBINED FIGURE")
print("=" * 70)

valid_dp = [dp for dp in ALL_KEYS if len(dp_data[dp]['hits']) > 0]
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

# ── Combined-only figure (just PD+PE+PT and GE+GD) ──
# Uses one-tailed Mann-Whitney (hypothesis: top hits have LOWER PSI)
print("\n  Generating combined-only figure...")
combined_only = [k for k in COMBINED_GROUPS.keys() if len(dp_data[k]['hits']) > 0]
n_co = len(combined_only)

pvals_1t = {}
cohen_d_co = {}
for dp in combined_only:
    hits_psi = dp_data[dp]['hits']['PSI_AAVS'].dropna().values
    all_psi  = dp_data[dp]['all']['PSI_AAVS'].dropna().values
    if len(hits_psi) >= 2 and len(all_psi) >= 2:
        _, p_1t = stats.mannwhitneyu(hits_psi, all_psi, alternative='less')
        pvals_1t[dp] = p_1t
        pooled = np.sqrt((np.var(hits_psi, ddof=1) + np.var(all_psi, ddof=1)) / 2)
        cohen_d_co[dp] = (np.mean(all_psi) - np.mean(hits_psi)) / pooled if pooled > 0 else 0
    else:
        pvals_1t[dp] = np.nan
        cohen_d_co[dp] = np.nan
    print(f"    {dp}: one-tailed p={pvals_1t[dp]:.4f} ({significance_marker(pvals_1t[dp])}), "
          f"Cohen's d={cohen_d_co[dp]:.2f}")

fig, ax = plt.subplots(figsize=(4.0 * n_co, 6))

positions_co = []
data_co = []
colors_co = []
ticks_co = []
tlabels_co = []
centers_co = []
brackets_co = []

pos = 0
for dp in combined_only:
    hits_psi = dp_data[dp]['hits']['PSI_AAVS'].dropna().values
    all_psi  = dp_data[dp]['all']['PSI_AAVS'].dropna().values

    p_all = pos
    p_hits = pos + 1

    positions_co.extend([p_all, p_hits])
    data_co.extend([all_psi, hits_psi])
    colors_co.extend(['#D5D8DC', DIPEPTIDE_COLORS[dp]])

    ticks_co.extend([p_all, p_hits])
    tlabels_co.extend([
        f'All Screen\n(n={len(all_psi)})',
        f'Top Hits\n(n={len(hits_psi)})',
    ])

    centers_co.append((p_all + p_hits) / 2)
    brackets_co.append((p_all, p_hits, dp))
    pos += 3.5

bp = ax.boxplot(data_co, positions=positions_co, widths=0.55,
                patch_artist=True, showfliers=False,
                whiskerprops={'linewidth': 1.3, 'color': 'black'},
                capprops={'linewidth': 1.3, 'color': 'black'},
                medianprops={'linewidth': 2.0, 'color': 'black'})

for patch, c in zip(bp['boxes'], colors_co):
    patch.set_facecolor(c)
    patch.set_edgecolor('black')
    patch.set_linewidth(1.0)
    patch.set_alpha(0.85)

ax.set_xticks(ticks_co)
ax.set_xticklabels(tlabels_co, fontsize=10)
ax.set_ylabel('Protein Stability Index (PSI)', fontsize=12)
ax.set_title('Combined Dipeptide Groups (P2-P3): Top Hits vs All Screen',
             fontsize=13, fontweight='bold', pad=18)
ax.text(0.5, 1.01, r'One-tailed Mann-Whitney U ($H_1$: Top Hits < All Screen)',
        transform=ax.transAxes, ha='center', va='bottom', fontsize=10, color='#555')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

for center, (_, _, dp) in zip(centers_co, brackets_co):
    ax.text(center, -0.18, dp, transform=ax.get_xaxis_transform(),
            ha='center', fontsize=13, fontweight='bold',
            color=DIPEPTIDE_COLORS[dp])
    members_str = ' + '.join(COMBINED_GROUPS[dp])
    ax.text(center, -0.25, f'({members_str})', transform=ax.get_xaxis_transform(),
            ha='center', fontsize=9, color='#666')

for p_all, p_hits, dp in brackets_co:
    ax.plot([p_all - 0.3, p_hits + 0.3], [-0.14, -0.14],
            transform=ax.get_xaxis_transform(),
            color=DIPEPTIDE_COLORS[dp], linewidth=2.5, clip_on=False)

global_max_co = max(np.percentile(d, 99) for d in data_co if len(d) > 0)
bracket_y_co = global_max_co + 0.15

for p_all, p_hits, dp in brackets_co:
    p = pvals_1t.get(dp, np.nan)
    d_val = cohen_d_co.get(dp, np.nan)
    if not np.isnan(p):
        sig = significance_marker(p)
        draw_bracket(ax, p_all, p_hits, bracket_y_co, 0.06, sig, fontsize=12)
        ax.text((p_all + p_hits) / 2, bracket_y_co + 0.18,
                f'p = {p:.3f}, d = {d_val:.2f}',
                ha='center', va='bottom', fontsize=8, color='#555')

ax.set_ylim(top=bracket_y_co + 0.55)

plt.tight_layout()
fig.subplots_adjust(bottom=0.28)
fig.savefig(os.path.join(OUTPUT_DIR, 'Fig_Combined_only_boxplot.png'))
fig.savefig(os.path.join(OUTPUT_DIR, 'Fig_Combined_only_boxplot.pdf'))
plt.close()
print("  Saved Fig_Combined_only_boxplot")

# ── Cohen's d version ──
print("\n  Generating Cohen's d version...")

def cohen_d_label(d):
    if abs(d) >= 0.8: return 'large'
    if abs(d) >= 0.5: return 'medium'
    if abs(d) >= 0.2: return 'small'
    return 'negligible'

fig, ax = plt.subplots(figsize=(4.0 * n_co, 6))

positions_cd = []
data_cd = []
colors_cd = []
ticks_cd = []
tlabels_cd = []
centers_cd = []
brackets_cd = []

pos = 0
for dp in combined_only:
    hits_psi = dp_data[dp]['hits']['PSI_AAVS'].dropna().values
    all_psi  = dp_data[dp]['all']['PSI_AAVS'].dropna().values

    p_all = pos
    p_hits = pos + 1

    positions_cd.extend([p_all, p_hits])
    data_cd.extend([all_psi, hits_psi])
    colors_cd.extend(['#D5D8DC', DIPEPTIDE_COLORS[dp]])

    ticks_cd.extend([p_all, p_hits])
    tlabels_cd.extend([
        f'All Screen\n(n={len(all_psi)})',
        f'Top Hits\n(n={len(hits_psi)})',
    ])

    centers_cd.append((p_all + p_hits) / 2)
    brackets_cd.append((p_all, p_hits, dp))
    pos += 3.5

bp = ax.boxplot(data_cd, positions=positions_cd, widths=0.55,
                patch_artist=True, showfliers=False,
                whiskerprops={'linewidth': 1.3, 'color': 'black'},
                capprops={'linewidth': 1.3, 'color': 'black'},
                medianprops={'linewidth': 2.0, 'color': 'black'})

for patch, c in zip(bp['boxes'], colors_cd):
    patch.set_facecolor(c)
    patch.set_edgecolor('black')
    patch.set_linewidth(1.0)
    patch.set_alpha(0.85)

ax.set_xticks(ticks_cd)
ax.set_xticklabels(tlabels_cd, fontsize=10)
ax.set_ylabel('Protein Stability Index (PSI)', fontsize=12)
ax.set_title('Combined Dipeptide Groups (P2-P3): Top Hits vs All Screen',
             fontsize=13, fontweight='bold', pad=18)
ax.text(0.5, 1.01, "Cohen's d effect size (All Screen vs Top Hits)",
        transform=ax.transAxes, ha='center', va='bottom', fontsize=10, color='#555')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

for center, (_, _, dp) in zip(centers_cd, brackets_cd):
    ax.text(center, -0.18, dp, transform=ax.get_xaxis_transform(),
            ha='center', fontsize=13, fontweight='bold',
            color=DIPEPTIDE_COLORS[dp])
    members_str = ' + '.join(COMBINED_GROUPS[dp])
    ax.text(center, -0.25, f'({members_str})', transform=ax.get_xaxis_transform(),
            ha='center', fontsize=9, color='#666')

for p_all, p_hits, dp in brackets_cd:
    ax.plot([p_all - 0.3, p_hits + 0.3], [-0.14, -0.14],
            transform=ax.get_xaxis_transform(),
            color=DIPEPTIDE_COLORS[dp], linewidth=2.5, clip_on=False)

global_max_cd = max(np.percentile(d, 99) for d in data_cd if len(d) > 0)
bracket_y_cd = global_max_cd + 0.15

for p_all, p_hits, dp in brackets_cd:
    d_val = cohen_d_co.get(dp, np.nan)
    p_val = pvals_1t.get(dp, np.nan)
    if not np.isnan(d_val):
        label = cohen_d_label(d_val)
        sig = significance_marker(p_val) if not np.isnan(p_val) else ''
        draw_bracket(ax, p_all, p_hits, bracket_y_cd, 0.06, sig, fontsize=12)
        ax.text((p_all + p_hits) / 2, bracket_y_cd + 0.16,
                f'd = {d_val:.2f} ({label} effect)',
                ha='center', va='bottom', fontsize=8.5, color='#555',
                style='italic')
        ax.text((p_all + p_hits) / 2, bracket_y_cd + 0.30,
                f'p = {p_val:.3f}',
                ha='center', va='bottom', fontsize=8, color='#555')

ax.set_ylim(top=bracket_y_cd + 0.65)

plt.tight_layout()
fig.subplots_adjust(bottom=0.28)
fig.savefig(os.path.join(OUTPUT_DIR, 'Fig_Combined_only_cohen_d.png'))
fig.savefig(os.path.join(OUTPUT_DIR, 'Fig_Combined_only_cohen_d.pdf'))
plt.close()
print("  Saved Fig_Combined_only_cohen_d")

# ══════════════════════════════════════════════════════════════════════════
# SUMMARY TABLE & EXCEL EXPORT
# ══════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

summary_rows = []
for dp in ALL_KEYS:
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
for dp in ALL_KEYS:
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

    for dp in ALL_KEYS:
        h = dp_data[dp]['hits']
        if len(h) > 0:
            safe = dp.replace('+', '_')
            cols = [c for c in ['Gene_ID', 'AA_seq', 'AA2', 'AA3', 'dipeptide', 'PSI_AAVS'] if c in h.columns]
            h[cols].to_excel(writer, sheet_name=f'TopHits_{safe}', index=False)

print(f"\n  Excel saved: {output_excel}")
print(f"  Figures saved to: {OUTPUT_DIR}")
print("\nDone!")
