"""
Regenerate Fig_Combined_only_cohen_d as a true vector SVG
(same data + styling as dipeptide_boxplot_analysis.py, Brunner-Munzel version).
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

COMBINED_GROUPS = {
    'PD+PE+PT': ['PD', 'PE', 'PT'],
    'GE+GD':    ['GE', 'GD'],
}

DIPEPTIDE_COLORS = {
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

def cohen_d_label(d):
    if abs(d) >= 0.8: return 'large'
    if abs(d) >= 0.5: return 'medium'
    if abs(d) >= 0.2: return 'small'
    return 'negligible'

def draw_bracket(ax, x1, x2, y, h, text, fontsize=11):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.2, color='black')
    ax.text((x1 + x2) / 2, y + h, text, ha='center', va='bottom',
            fontsize=fontsize, fontweight='bold')

combined_only = list(COMBINED_GROUPS.keys())
n_co = len(combined_only)

pvals_bm = {}
cohen_d_co = {}
for dp in combined_only:
    hits_psi = dp_data[dp]['hits']['PSI_AAVS'].dropna().values
    all_psi  = dp_data[dp]['all']['PSI_AAVS'].dropna().values
    _, p_bm = stats.brunnermunzel(hits_psi, all_psi, alternative='less')
    pvals_bm[dp] = p_bm
    pooled = np.sqrt((np.var(hits_psi, ddof=1) + np.var(all_psi, ddof=1)) / 2)
    cohen_d_co[dp] = (np.mean(all_psi) - np.mean(hits_psi)) / pooled

fig, ax = plt.subplots(figsize=(4.0 * n_co, 6))

positions_cd, data_cd, colors_cd = [], [], []
ticks_cd, tlabels_cd = [], []
centers_cd, brackets_cd = [], []

pos = 0
for dp in combined_only:
    hits_psi = dp_data[dp]['hits']['PSI_AAVS'].dropna().values
    all_psi  = dp_data[dp]['all']['PSI_AAVS'].dropna().values
    p_all, p_hits = pos, pos + 1
    positions_cd.extend([p_all, p_hits])
    data_cd.extend([all_psi, hits_psi])
    colors_cd.extend(['#D5D8DC', DIPEPTIDE_COLORS[dp]])
    ticks_cd.extend([p_all, p_hits])
    tlabels_cd.extend([f'All Screen\n(n={len(all_psi)})',
                       f'Top Hits\n(n={len(hits_psi)})'])
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
ax.text(0.5, 1.01, r"Brunner-Munzel test ($H_1$: Top Hits < All Screen) + Cohen's d",
        transform=ax.transAxes, ha='center', va='bottom', fontsize=10, color='#555')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

for center, (_, _, dp) in zip(centers_cd, brackets_cd):
    ax.text(center, -0.26, dp, transform=ax.get_xaxis_transform(),
            ha='center', fontsize=13, fontweight='bold',
            color=DIPEPTIDE_COLORS[dp])
    members_str = ' + '.join(COMBINED_GROUPS[dp])
    ax.text(center, -0.33, f'({members_str})', transform=ax.get_xaxis_transform(),
            ha='center', fontsize=9, color='#666')

for p_all, p_hits, dp in brackets_cd:
    ax.plot([p_all - 0.3, p_hits + 0.3], [-0.22, -0.22],
            transform=ax.get_xaxis_transform(),
            color=DIPEPTIDE_COLORS[dp], linewidth=2.5, clip_on=False)

global_max_cd = max(np.percentile(d, 99) for d in data_cd if len(d) > 0)
bracket_y_cd = global_max_cd + 0.15

for p_all, p_hits, dp in brackets_cd:
    d_val = cohen_d_co.get(dp, np.nan)
    p_val = pvals_bm.get(dp, np.nan)
    label = cohen_d_label(d_val)
    sig = significance_marker(p_val)
    draw_bracket(ax, p_all, p_hits, bracket_y_cd, 0.06, sig, fontsize=12)
    ax.text((p_all + p_hits) / 2, bracket_y_cd + 0.30,
            f'd = {d_val:.2f} ({label} effect)',
            ha='center', va='bottom', fontsize=8.5, color='#555',
            style='italic')
    ax.text((p_all + p_hits) / 2, bracket_y_cd + 0.46,
            f'p = {p_val:.4f}',
            ha='center', va='bottom', fontsize=8, color='#555')

ax.set_ylim(top=bracket_y_cd + 0.80)

plt.tight_layout()
fig.subplots_adjust(bottom=0.34)

svg_path = os.path.join(OUTPUT_DIR, 'Fig_Combined_only_cohen_d.svg')
fig.savefig(svg_path, format='svg', bbox_inches='tight')
plt.close(fig)
print(f"Saved: {svg_path}")
