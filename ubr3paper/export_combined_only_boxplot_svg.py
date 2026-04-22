"""
Regenerate Fig_Combined_only_boxplot as a true vector SVG
(extracted from dipeptide_boxplot_analysis.py, same data and styling).
"""
import os, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats

OUTPUT_DIR = r'c:\Users\User\Desktop\תינוקת\dianaMega\ubr3enrichmentlogo\dipeptide_analysis'
SCREEN_EXCEL = r'c:\Users\User\Downloads\UBR3 Nt screen (1).xlsx'

DIPEPTIDES = ['PD', 'PE', 'PT', 'GE', 'GD']
COMBINED_GROUPS = {
    'PD+PE+PT': ['PD', 'PE', 'PT'],
    'GE+GD':    ['GE', 'GD'],
}
ALL_KEYS = DIPEPTIDES + list(COMBINED_GROUPS.keys())

DIPEPTIDE_COLORS = {
    'PD': '#E74C3C', 'PE': '#E67E22', 'PT': '#F1C40F',
    'GE': '#2ECC71', 'GD': '#2E86C1',
    'PD+PE+PT': '#E67E22', 'GE+GD': '#DC143C',
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
    'svg.fonttype': 'none',
})

df_all  = pd.read_excel(SCREEN_EXCEL, sheet_name='Nprot_5_analyzed')
df_hits = pd.read_excel(SCREEN_EXCEL, sheet_name='sub_high')

df_all['PSI_AAVS']  = (df_all['PSI-293a']  + df_all['PSI-293b'])  / 2
df_hits['PSI_AAVS'] = (df_hits['PSI-293a'] + df_hits['PSI-293b']) / 2
df_all['dipeptide']  = df_all['AA2'].astype(str)  + df_all['AA3'].astype(str)
df_hits['dipeptide'] = df_hits['AA2'].astype(str) + df_hits['AA3'].astype(str)

dp_data = {}
for grp_name, members in COMBINED_GROUPS.items():
    h = pd.concat([df_hits[df_hits['dipeptide'] == m] for m in members], ignore_index=True)
    a = pd.concat([df_all[df_all['dipeptide'] == m]  for m in members], ignore_index=True)
    dp_data[grp_name] = {'hits': h, 'all': a}

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

fig, ax = plt.subplots(figsize=(4.0 * n_co, 6))

positions_co, data_co, colors_co = [], [], []
ticks_co, tlabels_co = [], []
centers_co, brackets_co = [], []

pos = 0
for dp in combined_only:
    hits_psi = dp_data[dp]['hits']['PSI_AAVS'].dropna().values
    all_psi  = dp_data[dp]['all']['PSI_AAVS'].dropna().values
    p_all, p_hits = pos, pos + 1
    positions_co.extend([p_all, p_hits])
    data_co.extend([all_psi, hits_psi])
    colors_co.extend(['#D5D8DC', DIPEPTIDE_COLORS[dp]])
    ticks_co.extend([p_all, p_hits])
    tlabels_co.extend([f'All Screen\n(n={len(all_psi)})',
                       f'Top Hits\n(n={len(hits_psi)})'])
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

svg_path = os.path.join(OUTPUT_DIR, 'Fig_Combined_only_boxplot.svg')
fig.savefig(svg_path, format='svg', bbox_inches='tight')
plt.close(fig)
print(f"Saved: {svg_path}")
