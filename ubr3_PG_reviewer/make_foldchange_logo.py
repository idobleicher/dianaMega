"""
Fold-change logo for positions 4-24, drawn exactly in the style of
export_fig4_logo_no_grid.py (fig4_logo_pos2_pos3_besthits).

Same palette, same black glyph outlines, same absence of gridlines, same
legend placement, same italic n-caption, same rcParams and font sizes.
Only the data and the number of positions differ.

Source: data/AMINO_ACIDS_PG_motif.csv -- the P/G-D/E/T motif analysis.
Substrates (n = 20) vs. non-substrate controls with the same motif (n = 193).
The CSV stores log2 fold change, so it is converted back to linear fold
change here; as in fig4, residues with fold change <= 1 are dropped.
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

# ------------------------------------------------------------------ load
raw = pd.read_csv(SRC, header=None, dtype=str, keep_default_na=False)


def block(row0):
    b = raw.iloc[row0:row0 + 20, 1:22].copy()
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


V = block(1).map(to_num)                                  # log2 fold change
S = block(23).map(lambda x: str(x).strip().count("*"))    # significance stars
assert list(V.index) == AA and list(S.index) == AA, "row order mismatch"

fc_matrix = (2.0 ** V).T          # positions as rows, linear fold change
fc_matrix.index = POS

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


def build(fc_display, title, stem, figsize):
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

    legend_patches = [mpatches.Patch(color=CAT_COLORS[c], label=c) for c in
                      ['Acidic', 'Basic', 'Nonpolar', 'Polar']]
    ax.legend(handles=legend_patches, fontsize=5.5, frameon=False,
              loc='upper left', bbox_to_anchor=(1.02, 1.0),
              handlelength=0.8, handletextpad=0.3, borderpad=0.2)
    plt.tight_layout(pad=0.3)

    for ext in ('png', 'pdf', 'svg'):
        path = os.path.join(OUTPUT_DIR, f'{stem}.{ext}')
        fig.savefig(path, format=ext)
        print(f"  saved: {path}  ({os.path.getsize(path)/1024:.1f} KB)")
    plt.close(fig)


# all enriched residues -- the direct fig4 analogue
disp_all = fc_matrix.copy()
disp_all[~(disp_all > 1)] = 0
build(disp_all,
      'Fold-Change Logo \u2014 P/G-D/E/T Substrates vs Non-Substrate Controls',
      'Figure13_foldchange_logo_pos4_24', (10.5, 2.9))

# same style, restricted to the residues carrying significance stars
disp_sig = disp_all.where(S.T.set_axis(POS).ge(1), 0)
build(disp_sig,
      'Fold-Change Logo \u2014 Significantly Enriched Residues (Positions 4\u201324)',
      'Figure13b_foldchange_logo_pos4_24_significant', (10.5, 2.9))
