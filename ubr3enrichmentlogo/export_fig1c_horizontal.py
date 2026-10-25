"""
Fig. 1c of Fig_V3 revision.pptx -- the Best Hits / Library enrichment heatmap
(the vertical fig1b_enrichment_heatmap) laid out HORIZONTALLY: residues run
along the x axis, Position 2 and Position 3 are the two rows. Same data, order
(by Position 2 enrichment), colour ramp and 0-5 scale as export_fig1b_svg.py;
the cells are a little smaller and every label reads left to right.

The significance stars are drawn in the figure, above each number, since the
star text boxes laid over the vertical panel in the slide cannot follow the
cells into the new layout.
"""
import os, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import font_manager, patheffects
from matplotlib.colors import LinearSegmentedColormap
from collections import Counter
from scipy import stats

EXCEL_PATH = r'c:\Users\User\Desktop\תינוקת\dianaMega\ubr3freqencyenrichment\UBR3 Nt screen.xlsx'
OUTPUT_DIR = r'c:\Users\User\Desktop\תינוקת\dianaMega\ubr3enrichmentlogo\figures'
os.makedirs(OUTPUT_DIR, exist_ok=True)

AMINO_ACIDS = list('ACDEFGHIKLMNPQRSTVWY')

# Helvetica is not installed on this machine, so matplotlib was falling back to
# Arial. TeX Gyre Heros (fonts/, GUST font licence) is a free Helvetica clone.
FONT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fonts')
for f in ('texgyreheros-regular.otf', 'texgyreheros-bold.otf'):
    font_manager.fontManager.addfont(os.path.join(FONT_DIR, f))

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['TeX Gyre Heros', 'Helvetica', 'Arial'],
    'font.size': 8,
    'axes.titlesize': 9,
    'axes.labelsize': 8,
    'axes.linewidth': 0.8,
    'svg.fonttype': 'path',   # outlines, so the SVG shows Helvetica without the font installed
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
M = sorted_by_pos2.T          # rows = positions, columns = residues

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

# 20 columns x 2 rows; 5.0 x 1.965 in matches the box the panel was resized to on the slide
# Look of Fig. 1g (ubr3_PG_reviewer Figure 17b, make_significance_heatmap.py):
# its warm ramp (pg_motif_data.CMAP_SIG), white cell separators, no frame or
# tick marks, near-black text and a hairline colour-bar outline.
CMAP_SIG = LinearSegmentedColormap.from_list('sig_warm', [
    '#F7F4EC', '#F4EAD0', '#F1DCA6', '#EDC76E', '#E8A33D',
    '#DC7C31', '#C0392B', '#9B2A20', '#71201A'])
GRID, INK, RULE = '#FFFFFF', '#1A1A1A', '#C9C4BA'

fig, ax = plt.subplots(figsize=(5.0, 1.965))
nrow, ncol = M.shape
im = ax.imshow(M.values, aspect='auto', cmap=CMAP_SIG, vmin=0,
               vmax=max(M.values.max(), 5),
               extent=(-0.5, ncol - 0.5, nrow - 0.5, -0.5))
for x in np.arange(-0.5, ncol, 1):
    ax.axvline(x, color=GRID, linewidth=0.8)
for y in np.arange(-0.5, nrow, 1):
    ax.axhline(y, color=GRID, linewidth=0.8)
for sp in ax.spines.values():
    sp.set_visible(False)
ax.set_xticks(range(ncol))
ax.set_xticklabels(M.columns, fontsize=6, color=INK)
ax.set_yticks([0, 1])
# The slide numbers residues after the initiator Met: sheet AA2 is its
# "N-terminal Position 1", AA3 its "Position 2" (as in panel e).
ax.set_yticklabels(['N-terminal\nPosition 1', 'N-terminal\nPosition 2'],
                   fontsize=5.5, color=INK)

# values from ubr3_PG_reviewer/pg_motif_data.py (CELL_TEXT, CELL_HALO, ...)
CELL_TEXT, CELL_HALO = '#FFFFFF', '#4A1410'
HALO_LW, STAR_HALO_LW = 0.9, 0.6   # 17b's 1.15 / 0.75, scaled to the smaller text
NUM_SIZE, STAR_BUMP = 5.4, 1.8

def halo(lw):
    return [patheffects.withStroke(linewidth=lw, foreground=CELL_HALO, alpha=0.9)]

for i, pos in enumerate(M.index):
    for j, aa in enumerate(M.columns):
        v = M.values[i, j]
        p = pval_dict[(aa, pos)]
        stars = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
        if v <= 1:
            stars = ''    # the slide stars enrichments only (A pos 2, G pos 3 are depleted)
        # Cell text styled as in ubr3_PG_reviewer Figures 16c/17b
        # (pg_motif_data.annotate_cell): white bold on every cell with a thin
        # dark halo, stars stacked above and set larger. Every number sits on
        # the same line, starred or not, so the row reads level.
        if stars:
            ax.text(j, i - 0.12, stars, ha='center', va='center',
                    fontsize=NUM_SIZE + STAR_BUMP, fontweight='bold',
                    color=CELL_TEXT, path_effects=halo(STAR_HALO_LW))
        ax.text(j, i + 0.08, f'{v:.1f}', ha='center', va='center',
                fontsize=NUM_SIZE, fontweight='bold', color=CELL_TEXT,
                path_effects=halo(HALO_LW))

cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.015, aspect=8)
cbar.set_ticks(range(0, 6))
cbar.set_label('Enrichment', fontsize=6.5, color=INK)
cbar.ax.tick_params(labelsize=6, length=2, width=0.6, colors=INK)
cbar.outline.set_linewidth(0.5)
cbar.outline.set_edgecolor(RULE)
# no title: the slide carries its own ("Amino acid enrichment at the N-terminus")
ax.tick_params(length=0, pad=2)
plt.tight_layout()

print(f'n = {n_hit} hits, n = {n_lib:,} library')
for (aa, pos), p in sorted(pval_dict.items(), key=lambda kv: kv[1]):
    if p < 0.05:
        print(f'  {pos}  {aa}  enrichment {enrich_matrix.loc[aa, pos]:.2f}  p = {p:.2g}')

base = 'fig1c_enrichment_heatmap_horizontal'
for ext in ('png', 'pdf', 'svg'):
    path = os.path.join(OUTPUT_DIR, f'{base}.{ext}')
    fig.savefig(path, format=ext)
    print(f"Saved: {path}")
plt.close(fig)
