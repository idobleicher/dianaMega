"""
Regenerate fig4_logo_pos2_pos3_besthits without the horizontal gridlines.
Saves PNG + PDF + SVG (vector).
"""
import os, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
import logomaker
from collections import Counter

EXCEL_PATH = r'c:\Users\User\Downloads\UBR3 Nt screen (1).xlsx'
OUTPUT_DIR = r'c:\Users\User\Desktop\תינוקת\dianaMega\ubr3enrichmentlogo\figures'
os.makedirs(OUTPUT_DIR, exist_ok=True)

AMINO_ACIDS = list('ACDEFGHIKLMNPQRSTVWY')

AA_CATEGORIES = {
    'D': 'Acidic',   'E': 'Acidic',
    'R': 'Basic',    'K': 'Basic',    'H': 'Basic',
    'G': 'Nonpolar', 'A': 'Nonpolar', 'V': 'Nonpolar', 'L': 'Nonpolar',
    'I': 'Nonpolar', 'P': 'Nonpolar', 'F': 'Nonpolar', 'M': 'Nonpolar', 'W': 'Nonpolar',
    'S': 'Polar',    'T': 'Polar',    'C': 'Polar',
    'Y': 'Polar',    'N': 'Polar',    'Q': 'Polar',
}
CAT_COLORS = {
    'Acidic':   '#922B21',
    'Basic':    '#C0392B',
    'Nonpolar': '#D35400',
    'Polar':    '#F39C12',
}
AA_COLOR_SCHEME = {aa: CAT_COLORS[cat] for aa, cat in AA_CATEGORIES.items()}

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'Liberation Sans', 'DejaVu Sans'],
    'font.size': 8,
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

def enrichment(hit_freq, lib_freq):
    return {aa: (hit_freq.get(aa, 0) / lib_freq[aa]) if lib_freq.get(aa, 0) > 0 else 0
            for aa in AMINO_ACIDS}

enrich_aa2 = enrichment(get_freq(df_hits['AA2']), get_freq(df_all['AA2']))
enrich_aa3 = enrichment(get_freq(df_hits['AA3']), get_freq(df_all['AA3']))

fc_matrix = pd.DataFrame({0: enrich_aa2, 1: enrich_aa3}).T
fc_display = fc_matrix.copy()
fc_display[fc_display <= 1] = 0

fig, ax = plt.subplots(figsize=(3.15, 2.55))
logo = logomaker.Logo(fc_display, ax=ax, color_scheme=AA_COLOR_SCHEME,
                      font_name='Helvetica', vpad=0.0, width=0.95,
                      stack_order='big_on_top')
for patch in ax.patches:
    patch.set_edgecolor('black')
    patch.set_linewidth(0.35)

ax.axhline(y=0, color='black', linewidth=0.5)
ax.yaxis.grid(False)
ax.xaxis.grid(False)
ax.set_axisbelow(True)

ax.set_ylabel('Fold Change (hits / library)', fontsize=8)
ax.set_xticks([0, 1])
ax.set_xticklabels(['Pos 2', 'Pos 3'], fontsize=8)
ax.tick_params(axis='both', labelsize=6.5, length=2.5, width=0.7, pad=2)
ax.set_title('Fold-Change Logo — Best Hits vs Library',
             fontweight='bold', fontsize=8.5, pad=4)

max_stack = max(sum(v for v in enr.values() if v > 1)
                for enr in [enrich_aa2, enrich_aa3])
ax.set_ylim(-0.3, max_stack + 2.0)

ax.text(0.5, -0.10, f'n = {len(df_hits)} best hits  vs  n = {len(df_all):,} library',
        transform=ax.transAxes, ha='center', fontsize=5.8, color='#666', style='italic')

legend_patches = [mpatches.Patch(color=CAT_COLORS[c], label=c) for c in
                  ['Acidic', 'Basic', 'Nonpolar', 'Polar']]
ax.legend(handles=legend_patches, fontsize=5.5, frameon=False,
          loc='upper left', bbox_to_anchor=(1.02, 1.0),
          handlelength=0.8, handletextpad=0.3, borderpad=0.2)
plt.tight_layout(pad=0.3)

base = 'fig4_logo_pos2_pos3_besthits'
for ext in ('png', 'pdf', 'svg'):
    path = os.path.join(OUTPUT_DIR, f'{base}.{ext}')
    fig.savefig(path, format=ext)
    print(f"  saved: {path}  ({os.path.getsize(path)/1024:.1f} KB)")

plt.close(fig)
