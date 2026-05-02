"""
Regenerate fig1b_enrichment_heatmap as a true vector SVG
(extracted from peptide_analysis.py, same data and styling).
"""
import os, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
from collections import Counter
from scipy import stats

EXCEL_PATH = r'c:\Users\User\Downloads\UBR3 Nt screen (1).xlsx'
OUTPUT_DIR = r'c:\Users\User\Desktop\תינוקת\dianaMega\ubr3enrichmentlogo\figures'
os.makedirs(OUTPUT_DIR, exist_ok=True)

AMINO_ACIDS = list('ACDEFGHIKLMNPQRSTVWY')

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'Liberation Sans', 'DejaVu Sans'],
    'font.size': 8,
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'axes.linewidth': 0.8,
    'svg.fonttype': 'none',
    'figure.dpi': 600,
    'savefig.dpi': 600,
    'savefig.bbox': 'tight',
    'axes.spines.top': False,
    'axes.spines.right': False,
})

df_all  = pd.read_excel(EXCEL_PATH, sheet_name='Nprot_5_analyzed')
df_hits = pd.read_excel(EXCEL_PATH, sheet_name='sub_high')

def get_freq(series):
    c = Counter(series.dropna())
    t = sum(c.values())
    return {aa: c.get(aa, 0) / t for aa in AMINO_ACIDS}

def get_counts(series):
    c = Counter(series.dropna())
    return {aa: c.get(aa, 0) for aa in AMINO_ACIDS}

def enrichment(hit_freq, lib_freq):
    return {aa: (hit_freq.get(aa, 0) / lib_freq[aa]) if lib_freq.get(aa, 0) > 0 else 0
            for aa in AMINO_ACIDS}

lib_aa2_freq = get_freq(df_all['AA2'])
lib_aa3_freq = get_freq(df_all['AA3'])
hit_aa2_freq = get_freq(df_hits['AA2'])
hit_aa3_freq = get_freq(df_hits['AA3'])
enrich_aa2 = enrichment(hit_aa2_freq, lib_aa2_freq)
enrich_aa3 = enrichment(hit_aa3_freq, lib_aa3_freq)

enrich_matrix = pd.DataFrame({'Position 2': enrich_aa2, 'Position 3': enrich_aa3})
sorted_by_pos2 = enrich_matrix.sort_values('Position 2', ascending=False)

lib_counts_2 = get_counts(df_all['AA2'])
lib_counts_3 = get_counts(df_all['AA3'])
hit_counts_2 = get_counts(df_hits['AA2'])
hit_counts_3 = get_counts(df_hits['AA3'])
n_lib = len(df_all)
n_hit = len(df_hits)

pval_dict = {}
for aa in AMINO_ACIDS:
    for pos_label, hc, lc in [('Position 2', hit_counts_2, lib_counts_2),
                              ('Position 3', hit_counts_3, lib_counts_3)]:
        table = [[hc[aa], n_hit - hc[aa]],
                 [lc[aa], n_lib - lc[aa]]]
        _, p = stats.fisher_exact(table, alternative='two-sided')
        pval_dict[(aa, pos_label)] = p

fig, ax = plt.subplots(figsize=(2.75, 4.7))
cmap_warm = LinearSegmentedColormap.from_list(
    'warm', ['#FEF9E7', '#F5B041', '#E67E22', '#C0392B', '#78281F'])
im = ax.imshow(sorted_by_pos2.values, aspect='auto', cmap=cmap_warm, vmin=0,
               vmax=max(sorted_by_pos2.values.max(), 5))
ax.set_xticks([0, 1])
ax.set_xticklabels(['Position 2', 'Position 3'], fontsize=7.5)
ax.set_yticks(range(len(sorted_by_pos2)))
ax.set_yticklabels(sorted_by_pos2.index, fontsize=7.5)

col_labels = ['Position 2', 'Position 3']
for i, aa in enumerate(sorted_by_pos2.index):
    for j in range(2):
        v = sorted_by_pos2.values[i, j]
        p = pval_dict[(aa, col_labels[j])]
        stars = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
        color = 'white' if v > 2.5 else 'black'
        if stars:
            ax.text(j + 0.03, i, f'{v:.1f}', ha='right', va='center',
                    fontsize=5.8, color=color, fontweight='normal')
            ax.text(j + 0.03, i, stars, ha='left', va='center',
                    fontsize=6.8, color=color, fontweight='bold')
        else:
            ax.text(j, i, f'{v:.1f}', ha='center', va='center',
                    fontsize=5.8, color=color, fontweight='normal')

cbar = plt.colorbar(im, ax=ax, shrink=0.62, pad=0.06, label='Enrichment')
cbar.ax.tick_params(labelsize=6, length=2, width=0.6)
cbar.outline.set_linewidth(0.6)
ax.set_title(f'Enrichment Heatmap (Best Hits / Library)\n'
             f'(n={n_hit} hits, n={n_lib:,} library)', fontweight='bold')
fig.text(0.5, 0.005, '* p<0.05  ** p<0.01  *** p<0.001 (Fisher exact test)',
         ha='center', fontsize=5.4, color='#555')
ax.tick_params(length=2, width=0.7, pad=1.5)
plt.tight_layout(rect=[0, 0.025, 1, 1])

base = 'fig1b_enrichment_heatmap'
for ext in ('png', 'pdf', 'svg'):
    path = os.path.join(OUTPUT_DIR, f'{base}.{ext}')
    fig.savefig(path, format=ext)
    print(f"Saved: {path}")
plt.close(fig)
