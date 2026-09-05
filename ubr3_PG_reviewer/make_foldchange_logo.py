"""
Fold-change logo for positions 4-24, drawn exactly in the style of
export_fig4_logo_no_grid.py (fig4_logo_pos2_pos3_besthits).

Same palette, same black glyph outlines, same absence of gridlines, same
legend placement, same italic n-caption, same rcParams and font sizes.
Only the data and the number of positions differ.

Data comes from pg_motif_data.py (see that module for the source and for the
arginine correction). Substrates (n = 20) vs. non-substrate controls carrying
the same motif (n = 193). Letter height is linear fold change; as in fig4,
residues with fold change <= 1 are dropped.

Figure13b keeps only the residues significant at p < 0.05. Significance is
taken from the workbook's chi-square p values and binned here, rather than
read from its asterisk matrix -- that matrix is not internally consistent with
its own p values (it marks Basic at position 5 with one star at p = 1.1e-4).
"""
import os
import textwrap

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import logomaker

import pg_motif_data as pg

HERE = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(HERE, 'figures')
os.makedirs(OUTPUT_DIR, exist_ok=True)

POS, AA = pg.POS, pg.AA
N_SUB, N_CTRL = pg.N_SUB, pg.N_CTRL

cells = pg.load()
fc_matrix = pg.matrix('fold_change', cells=cells)          # AA x POS, linear
P = pg.matrix('p_chi2', cells=cells)

def fig_width(ncol):
    """Constant glyph column width across windows: the 21-position logo is
    10.5 in wide, so the 9-position one is narrower by the columns it drops
    rather than the same width with stretched letters."""
    return 3.2 + 0.35 * ncol


def wrap(text, figw, per_inch=16.5):
    return textwrap.fill(text, width=max(40, int(figw * per_inch)))


# Every logo is drawn twice: over the full 4-24 range, and over the N-terminal
# window out to the 12th amino acid. Heights, colours and p values are
# identical in both -- the window is a display crop, not a re-analysis.
WINDOWS = [(POS, "pos4_24", ""),
           (pg.POS_N12, "pos4_12", " (positions 4–12)")]


# ------------------------------------------------------------------ style (verbatim from fig4)
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
    'pdf.fonttype': 42,
    'figure.dpi': 600,
    'savefig.dpi': 600,
    'savefig.bbox': 'tight',
    'axes.spines.top': False,
    'axes.spines.right': False,
})


def build(fc_display, title, stem, positions):
    figsize = (fig_width(len(positions)), 2.9)
    fig, ax = plt.subplots(figsize=figsize)
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

    ax.set_ylabel('Fold Change (substrates / controls)', fontsize=8)
    ax.set_xlabel('Position', fontsize=8, labelpad=2)
    ax.set_xticks(positions)
    ax.set_xticklabels(positions, fontsize=6.5)
    ax.tick_params(axis='both', labelsize=6.5, length=2.5, width=0.7, pad=2)
    ax.set_title(wrap(title, figsize[0]), fontweight='bold', fontsize=8.5,
                 pad=4)

    max_stack = float(fc_display.sum(axis=1).max())
    ax.set_ylim(-0.3, max_stack + 2.0)
    ax.set_xlim(positions[0] - 0.6, positions[-1] + 0.6)

    ax.text(0.5, -0.15,
            f'n = {N_SUB} P/G-D/E/T substrates  vs  n = {N_CTRL} non-substrate controls',
            transform=ax.transAxes, ha='center', fontsize=5.8,
            color='#666', style='italic')

    legend_patches = [mpatches.Patch(color=CAT_COLORS[c], label=c) for c in
                      ['Acidic', 'Basic', 'Nonpolar', 'Polar']]
    ax.legend(handles=legend_patches, fontsize=5.5, frameon=False,
              loc='upper left', bbox_to_anchor=(1.02, 1.0),
              handlelength=0.8, handletextpad=0.3, borderpad=0.2)
    plt.tight_layout(pad=0.3)

    for ext in ('png', 'pdf', 'svg'):
        path = os.path.join(OUTPUT_DIR, f'{stem}.{ext}')
        fig.savefig(path, format=ext)
        print(f"  saved: {os.path.basename(path)}  ({os.path.getsize(path)/1024:.1f} KB)")
    plt.close(fig)


for positions, tag, wt in WINDOWS:
    # all enriched residues -- the direct fig4 analogue
    build(pg.logo_frame('fold_change', cells=cells, positions=positions),
          'Fold-Change Logo — P/G-D/E/T Substrates vs Non-Substrate '
          'Controls' + wt,
          f'Figure13_foldchange_logo_{tag}', positions)

    # same style, restricted to the residues significant at p < 0.05
    build(pg.logo_frame('fold_change', cells=cells, only=P.lt(0.05),
                        positions=positions),
          f'Fold-Change Logo — Significantly Enriched Residues '
          f'(Positions {positions[0]}–{positions[-1]})',
          f'Figure13b_foldchange_logo_{tag}_significant', positions)
