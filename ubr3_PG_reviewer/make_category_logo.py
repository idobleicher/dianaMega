"""
Enrichment logo for positions 4-24 with amino acids coloured by residue
CATEGORY, in the house style of export_fig4_logo_no_grid.py: black glyph
outlines at 0.35, no gridlines, same rcParams, font sizes, legend placement
and italic n-caption.

Letters are the real amino acids and keep their own per-residue fold change
(as in Figure13); only the colouring changes, so nothing is aggregated away.
Category membership is a complete partition of all 20 residues, no overlaps:

    Basic        K H R
    Acidic       D E
    Aromatic     F Y W
    Aliphatic    G P A M
    Hydrophobic  V L I
    Polar        S T C N Q

Note: the pre-aggregated Categories.csv is NOT used. Its Basic row is
byte-identical to the arginine row of the older AMINO_ACIDS_PG_motif.csv --
because that file's arginine row is in fact the Basic-class row. Everything
here is rebuilt from the per-residue contingency tables via pg_motif_data.py,
which avoids that.

Data: pg_motif_data.py. Substrates (n = 20) vs. non-substrate controls
carrying the same motif (n = 193).
"""
import os

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import logomaker

import pg_motif_data as pg

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(HERE, 'figures')
os.makedirs(OUTPUT_DIR, exist_ok=True)

POS, AA = pg.POS, pg.AA
N_SUB, N_CTRL = pg.N_SUB, pg.N_CTRL
CATEGORY_MEMBERS, CAT_COLORS = pg.CATEGORY_MEMBERS, pg.CAT_COLORS
LEGEND_ORDER = ['Basic', 'Acidic', 'Aromatic', 'Aliphatic', 'Hydrophobic', 'Polar']
AA_COLOR_SCHEME = {a: CAT_COLORS[c] for c, mem in CATEGORY_MEMBERS.items() for a in mem}

cells = pg.load()
P = pg.matrix('p_chi2', cells=cells)

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


build(pg.logo_frame('fold_change', cells=cells),
      'Enrichment Logo by Residue Category — P/G-D/E/T Substrates vs Non-Substrate Controls',
      'Figure14_logo_pos4_24_by_category', (10.5, 2.9))

build(pg.logo_frame('fold_change', cells=cells, only=P.lt(0.05)),
      'Enrichment Logo by Residue Category — Significantly Enriched (Positions 4–24)',
      'Figure14b_logo_pos4_24_by_category_significant', (10.5, 2.9))

# ------------------------------------------------------------------ tables
res = cells[(cells.kind == 'residue') & cells.defined].copy()
sig = res[(res.p_chi2 < 0.05) & (res.direction > 0)].copy()
sig['Class'] = sig.label.map(pg.CAT_OF_AA)
sig['Significance'] = sig.p_chi2.map(pg.stars)
tab = (sig.rename(columns={'position': 'Position', 'label': 'Residue',
                           'fold_change': 'Fold change', 'p_chi2': 'p'})
          .assign(log2FC=lambda d: d.log2FC.round(3),
                  **{'Fold change': lambda d: d['Fold change'].round(2)})
          [['Position', 'Residue', 'Class', 'log2FC', 'Fold change', 'p',
            'q_chi2', 'Significance', 'n_sub', 'n_ctrl']]
          .sort_values(['Position', 'log2FC'], ascending=[True, False])
          .reset_index(drop=True))

tab.to_csv(os.path.join(HERE, 'data', 'significant_residues_by_category.csv'), index=False)
(tab.rename(columns={'Class': 'Category'})
    .to_csv(os.path.join(HERE, 'data', 'significant_enriched_residues_pos4_24.csv'),
            index=False))

print()
print(tab.groupby('Class').size().sort_values(ascending=False).to_string())
print(f"\n{len(tab)} residue cells enriched at p < 0.05 "
      f"(of {len(res)} testable; {0.05 * len(res):.0f} expected by chance)")
