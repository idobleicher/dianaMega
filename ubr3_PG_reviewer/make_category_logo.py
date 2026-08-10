"""
Enrichment logo for positions 4-24 with amino acids coloured by residue
CATEGORY, in the house style of export_fig4_logo_no_grid.py: black glyph
outlines at 0.35, no gridlines, same rcParams, font sizes, legend placement
and italic n-caption.

Letters are the real amino acids and keep their own per-residue fold change
(as in Figure13); only the colouring changes, so nothing is aggregated away.
Category membership is the user's grouping -- a complete partition of all 20
residues, no overlaps:

    Basic        K R H
    Acidic       D E
    Aromatic     F W Y
    Aliphatic    G P A M
    Hydrophobic  L I V
    Polar        S T C N Q

Note: the pre-aggregated Categories.csv is NOT used. Its Basic row is
byte-identical to the arginine row of AMINO ACIDS.csv across all 21
positions, and its category totals cannot be reproduced from the per-residue
values by this (or any tested) grouping. Rebuilding from the amino-acid file
avoids that.

Source: data/AMINO_ACIDS_PG_motif.csv -- the P/G-D/E/T motif analysis.
Substrates (n = 20) vs. non-substrate controls with the same motif (n = 193).
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import logomaker

HERE       = os.path.dirname(os.path.abspath(__file__))
SRC        = os.path.join(HERE, 'data', 'AMINO_ACIDS_PG_motif.csv')
OUTPUT_DIR = os.path.join(HERE, 'figures')
os.makedirs(OUTPUT_DIR, exist_ok=True)

AA  = ["K", "H", "R", "D", "E", "F", "Y", "W", "G", "P",
       "A", "M", "V", "L", "I", "S", "T", "C", "N", "Q"]
POS = list(range(4, 25))
N_SUB, N_CTRL = 20, 193

# ------------------------------------------------------------------ grouping
CATEGORY_MEMBERS = {
    'Basic':       list('KRH'),
    'Acidic':      list('DE'),
    'Aromatic':    list('FWY'),
    'Aliphatic':   list('GPAM'),
    'Hydrophobic': list('LIV'),
    'Polar':       list('STCNQ'),
}
assert sorted(sum(CATEGORY_MEMBERS.values(), [])) == sorted(AA), \
    "grouping must be a partition of all 20 residues"

# Warm house family, ordered by lightness so the six stay as separable as this
# palette allows. Six warm hues cannot be fully separated (fig4's own four
# already sit at dE 7.5 for Nonpolar vs Basic); the letters and the member
# lists in the legend are what actually identify a category.
CAT_COLORS = {
    'Acidic':      '#7B241C',
    'Aromatic':    '#A0522D',
    'Basic':       '#C0392B',
    'Hydrophobic': '#D35400',
    'Aliphatic':   '#E8A33D',
    'Polar':       '#F1C40F',
}
LEGEND_ORDER = ['Basic', 'Acidic', 'Aromatic', 'Aliphatic', 'Hydrophobic', 'Polar']
AA_COLOR_SCHEME = {a: CAT_COLORS[c] for c, mem in CATEGORY_MEMBERS.items() for a in mem}

# ------------------------------------------------------------------ load
raw = pd.read_csv(SRC, header=None, dtype=str, keep_default_na=False)


def block(row0, conv):
    b = raw.iloc[row0:row0 + 20, 1:22].map(conv)
    b.index, b.columns = raw.iloc[row0:row0 + 20, 0].tolist(), POS
    return b


def to_num(x):
    x = str(x).strip()
    if x in ("N/A", "", "-", "#NUM!", "nan"):
        return np.nan
    try:
        return float(x)
    except ValueError:
        return np.nan


V = block(1, to_num)                                       # log2 fold change
S = block(23, lambda x: str(x).strip().count("*"))         # significance stars
assert list(V.index) == AA and list(S.index) == AA, "row order mismatch"

fc_matrix = (2.0 ** V).T                                   # linear fold change
fc_matrix.index = POS

# ------------------------------------------------------------------ style
mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'Liberation Sans', 'DejaVu Sans'],
    'font.size': 8,
    'axes.linewidth': 0.8,
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'figure.dpi': 600,
    'savefig.dpi': 600,
    'savefig.bbox': 'tight',
    'axes.spines.top': False,
    'axes.spines.right': False,
})


def build(fc_display, title, stem, figsize):
    fig, ax = plt.subplots(figsize=figsize)
    logomaker.Logo(fc_display, ax=ax, color_scheme=AA_COLOR_SCHEME,
                   font_name='Helvetica', vpad=0.0, width=0.95,
                   stack_order='big_on_top')
    for patch in ax.patches:
        patch.set_edgecolor('black')
        patch.set_linewidth(0.35)

    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.yaxis.grid(False)
    ax.xaxis.grid(False)
    ax.set_axisbelow(True)

    ax.set_ylabel('Fold Change (substrates / controls)', fontsize=8)
    ax.set_xlabel('Position', fontsize=8, labelpad=2)
    ax.set_xticks(POS)
    ax.set_xticklabels(POS, fontsize=6.5)
    ax.tick_params(axis='both', labelsize=6.5, length=2.5, width=0.7, pad=2)
    ax.set_title(title, fontweight='bold', fontsize=8.5, pad=4)

    max_stack = float(fc_display.sum(axis=1).max())
    ax.set_ylim(-0.3, max_stack + 2.0)
    ax.set_xlim(POS[0] - 0.6, POS[-1] + 0.6)

    ax.text(0.5, -0.15,
            f'n = {N_SUB} P/G-D/E/T substrates  vs  n = {N_CTRL} non-substrate controls',
            transform=ax.transAxes, ha='center', fontsize=5.8,
            color='#666', style='italic')

    legend_patches = [
        mpatches.Patch(color=CAT_COLORS[c],
                       label=f'{c}  ({" ".join(CATEGORY_MEMBERS[c])})')
        for c in LEGEND_ORDER]
    ax.legend(handles=legend_patches, fontsize=5.5, frameon=False,
              loc='upper left', bbox_to_anchor=(1.02, 1.0),
              handlelength=0.8, handletextpad=0.3, borderpad=0.2)
    plt.tight_layout(pad=0.3)

    for ext in ('png', 'pdf', 'svg'):
        path = os.path.join(OUTPUT_DIR, f'{stem}.{ext}')
        fig.savefig(path, format=ext)
        print(f"  saved: {os.path.basename(path)}  ({os.path.getsize(path)/1024:.1f} KB)")
    plt.close(fig)


disp_all = fc_matrix.copy()
disp_all[~(disp_all > 1)] = 0
build(disp_all,
      'Enrichment Logo by Residue Category — P/G-D/E/T Substrates vs Non-Substrate Controls',
      'Figure14_logo_pos4_24_by_category', (10.5, 2.9))

disp_sig = disp_all.where(S.T.set_axis(POS).ge(1), 0)
build(disp_sig,
      'Enrichment Logo by Residue Category — Significantly Enriched (Positions 4–24)',
      'Figure14b_logo_pos4_24_by_category_significant', (10.5, 2.9))

# ------------------------------------------------------------------ table
cat_of = {a: c for c, mem in CATEGORY_MEMBERS.items() for a in mem}
rows = [{'Position': p, 'Residue': a, 'Category': cat_of[a],
         'log2FC': round(float(V.loc[a, p]), 3),
         'Fold change': round(float(2 ** V.loc[a, p]), 2),
         'Significance': '*' * int(S.loc[a, p])}
        for p in POS for a in AA if S.loc[a, p] >= 1 and V.loc[a, p] > 0]
tab = (pd.DataFrame(rows).sort_values(['Position', 'log2FC'], ascending=[True, False])
         .reset_index(drop=True))
tab.to_csv(os.path.join(HERE, 'data', 'significant_residues_by_category.csv'), index=False)

print()
print(tab.groupby('Category').size().sort_values(ascending=False).to_string())
