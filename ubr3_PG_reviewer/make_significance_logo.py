"""
Significance logo for positions 4-24: letter height is -log10(p) instead of
fold change, so a residue that is MORE SIGNIFICANT is TALLER regardless of how
large its fold change happens to be. This is the pLogo / kpLogo convention.

Why this exists: fold change and significance are different quantities. W at
position 12 has one of the largest fold changes in the data (9.65x) but only
reaches p = 0.048, because 9.65x is exactly what you get when a residue has
equal raw counts in both groups -- here 1 peptide vs 1 peptide. R at positions
7, 8 and 17 has a smaller fold change but a far smaller p, because it rests on
many more observations. A fold-change logo makes W look like the main result;
a significance logo makes R the main result.

Heights are EXACT -log10(p), not the * / ** / *** bin representatives the first
version of this figure had to use: pg_motif_data.py reads the workbook sheet
that carries a p value per cell. That module also documents the arginine
correction -- the earlier source file's R row was the Basic-class row, so R was
drawn far too short here.

COLOUR -- the glyphs are coloured by the same warm significance ramp as the
heatmaps (pg.CMAP_SIG), keyed to the same -log10 p, with the same colourbar and
the same p ticks. Height and colour therefore carry one quantity twice, which
is the point: a letter that is tall is also dark, so the figure reads at a
glance and sits beside Figures 16-17 as one visual set rather than as a second
palette. Chemical class is no longer carried by colour here -- that is Figure
14's job, which colours the same stacks by class. Pass colour_by="class" to
build() to get the old class-coloured version of this figure back.

There are no black outlines on the glyphs. The hairline that separates two
stacked letters is white, the same separator the heatmaps use between cells.

Two versions are produced, each over two position windows (pos4_24, the full
range, and pos4_12, out to the 12th amino acid -- a display crop on the same
scale, with q still corrected across the full 4-24 family):
  Figure15   chi-square, the workbook's own test.
  Figure15b  Fisher exact, recomputed from the same 2 x 2 tables. Chi-square
             needs expected counts >= 5 and most cells here expect fewer than
             2 substrates, so this is the version to trust.

Substrates (n = 20) vs. non-substrate controls carrying the same motif
(n = 193).
"""
import os
import textwrap

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap, Normalize
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
DIR = pg.matrix('direction', cells=cells)

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


# Only cells at p < 0.05 are ever drawn, so the ramp is entered at its p = 0.05
# point rather than at its pale end: the faintest glyph in the figure is the
# gold the heatmaps use for a one-star cell, not the near-white they use for
# n.s. VMAX matches the heatmaps' colour ceiling exactly.
VMIN, VMAX = 1.30103, 5.0
RAMP = LinearSegmentedColormap.from_list(
    'sig_glyph', pg.CMAP_SIG(np.linspace(0.34, 1.0, 256)))
NORM = Normalize(vmin=VMIN, vmax=VMAX, clip=True)

# ------------------------------------------------------------------ style
mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'Liberation Sans', 'DejaVu Sans'],
    'font.size': 8, 'axes.linewidth': 0.8,
    'svg.fonttype': 'none', 'pdf.fonttype': 42,
    'figure.dpi': 600, 'savefig.dpi': 600, 'savefig.bbox': 'tight',
    'figure.facecolor': 'white', 'savefig.facecolor': 'white',
    'axes.spines.top': False, 'axes.spines.right': False,
})


def build(pcol, test_label, stem, q_note, positions=None,
          colour_by='significance'):
    positions = list(positions or POS)
    figsize = (fig_width(len(positions)), 2.9)
    P = pg.matrix(pcol, cells=cells, positions=positions)
    Q = pg.matrix(pcol.replace('p_', 'q_'), cells=cells, positions=positions)
    direction = DIR[positions]

    # Height = -log10(p), for residues ENRICHED in the substrates and reaching
    # p < 0.05. The threshold matters: without it all 20 residues are drawn at
    # every position, most of them at p ~ 0.3-0.9, and the stacks pile up to
    # -log10 p ~ 8 out of noise alone. Drawing only what passes is the pLogo
    # convention and is what the first version of this figure did (via its
    # star bins); the heights are now exact rather than bin representatives.
    with np.errstate(divide='ignore', invalid='ignore'):
        H = -np.log10(P.astype(float))
    H = H.where(np.isfinite(H) & direction.gt(0) & P.lt(0.05),
                0.0).clip(lower=0.0)
    H = H.T.set_axis(positions)

    fig, ax = plt.subplots(figsize=figsize)
    base = AA_COLOR_SCHEME if colour_by == 'class' else '#CFCAC1'
    logo = logomaker.Logo(H, ax=ax, color_scheme=base,
                          font_name='Helvetica', vpad=0.0, width=0.95,
                          stack_order='big_on_top')

    # One colour per glyph, from that glyph's own p -- logomaker's colour_scheme
    # is per-residue and cannot do this, so the glyphs are restyled after the
    # stacks are laid out.
    if colour_by != 'class':
        for p_ in H.index:
            for a in H.columns:
                h = float(H.at[p_, a])
                if h > 0:
                    logo.style_single_glyph(p_, a, color=RAMP(NORM(h)))

    # No black outlines. A white hairline is the same separator the heatmaps
    # draw between cells, and it only shows where two glyphs actually touch.
    for patch in ax.patches:
        patch.set_edgecolor(pg.GRID)
        patch.set_linewidth(0.35)

    for y, lab in [(1.30103, 'p = 0.05'), (2.0, 'p = 0.01'), (3.0, 'p = 0.001')]:
        ax.axhline(y=y, color=pg.RULE, linewidth=0.5, linestyle=(0, (3, 2)),
                   zorder=1)
        ax.text(positions[-1] + 0.55, y, lab, fontsize=5, color=pg.MUTED,
                va='center', ha='left')

    ax.axhline(y=0, color=pg.RULE, linewidth=0.6)
    ax.yaxis.grid(False)
    ax.xaxis.grid(False)
    ax.set_axisbelow(True)
    ax.set_ylabel('$-$log$_{10}$ $p$   (stacked)', fontsize=8, color=pg.INK)
    ax.set_xlabel('Position', fontsize=8, labelpad=2, color=pg.INK)
    ax.set_xticks(positions)
    ax.set_xticklabels(positions, fontsize=6.5)
    ax.tick_params(axis='both', labelsize=6.5, length=2.5, width=0.7, pad=2,
                   colors=pg.INK)
    for sp in ax.spines.values():
        sp.set_edgecolor(pg.RULE)
    wt = '' if positions == list(POS) else (
        f'  —  positions {positions[0]}–{positions[-1]}')
    ax.set_title(wrap(f'Significance Logo — height and colour are '
                      f'$-$log$_{{10}}$ $p$, not fold change '
                      f'({test_label}){wt}', figsize[0]),
                 fontweight='bold', fontsize=8.5, pad=4, color=pg.INK)
    ax.set_ylim(-0.05, float(H.sum(axis=1).max()) * 1.12)
    ax.set_xlim(positions[0] - 0.6, positions[-1] + 0.6)
    # The caption wraps on the narrow window figure, so it is pushed down by
    # a line's worth of height per extra line: at -0.15 a three-line caption
    # would sit on top of the position numbers and the x label.
    cap = wrap(f'n = {N_SUB} P/G-D/E/T substrates  vs  n = {N_CTRL} '
               f'non-substrate controls.   {q_note}', figsize[0], 24)
    ax.text(0.5, -0.15 - 0.085 * cap.count('\n'), cap,
            transform=ax.transAxes, ha='center', fontsize=5.8, color=pg.MUTED,
            style='italic', linespacing=1.5)

    if colour_by == 'class':
        ax.legend(handles=[mpatches.Patch(color=CAT_COLORS[c],
                                          label=f'{c}  ({" ".join(CATEGORY_MEMBERS[c])})')
                           for c in LEGEND_ORDER],
                  fontsize=5.5, frameon=False, loc='upper left',
                  bbox_to_anchor=(1.10, 1.0), handlelength=0.8,
                  handletextpad=0.3, borderpad=0.2)
    else:
        # the same key the heatmaps carry, on the same ramp and the same ticks
        sm = mpl.cm.ScalarMappable(cmap=RAMP, norm=NORM)
        cbar = fig.colorbar(sm, ax=ax, fraction=0.020, pad=0.055)
        cbar.ax.tick_params(labelsize=6, length=2, width=0.6, colors=pg.INK)
        cbar.outline.set_linewidth(0.5)
        cbar.outline.set_edgecolor(pg.RULE)
        cbar.set_ticks([VMIN, 2, 3, 4, VMAX])
        cbar.set_ticklabels(['0.05', '0.01', '0.001', '1e-4', '1e-5'])
        cbar.set_label('$p$ value', fontsize=6.5, labelpad=5, color=pg.INK)

    plt.tight_layout(pad=0.3)

    for ext in ('png', 'pdf', 'svg'):
        p = os.path.join(OUTPUT_DIR, f'{stem}.{ext}')
        fig.savefig(p, format=ext)
        print(f"  saved: {os.path.basename(p)}  ({os.path.getsize(p)/1024:.1f} KB)")
    plt.close(fig)
    return P, Q


res = cells[(cells.kind == 'residue') & cells.defined]
n_test = len(res)


def note(pcol, positions):
    """The caption's own arithmetic, counted over the positions actually
    drawn. q stays the BH-FDR of the full 4-24 family: cropping the figure
    must not make a cell look more significant than it is."""
    q = pcol.replace('p_', 'q_')
    win = res[res.position.isin(positions)]
    n_win = len(win)
    n_hit, n_fdr = int((win[pcol] < 0.05).sum()), int((win[q] < 0.05).sum())
    survivors = ', '.join(f'{r.label}{r.position}' for _, r in
                          win[win[q] < 0.05].sort_values(q).iterrows())
    crop = ('' if len(positions) == len(POS) else
            f' This is the N-terminal window only — positions '
            f'{positions[0]}–{positions[-1]}, out to the '
            f'{positions[-1]}th amino acid; $q$ is still BH-FDR across all '
            f'{n_test} cells of the full {POS[0]}–{POS[-1]} matrix.')
    return (f'Heights are uncorrected $p$: {n_hit} of the {n_win} cells shown '
            f'reach $p$ < 0.05 against {0.05 * n_win:.0f} expected by chance, '
            + (f'and only {survivors} survive'
               f'{"s" if n_fdr == 1 else ""} BH-FDR.' if n_fdr else
               'and none survives BH-FDR.') + crop)


for positions, tag, _wt in WINDOWS:
    build('p_chi2', 'chi-square, as reported in the workbook',
          f'Figure15_significance_logo_{tag}', note('p_chi2', positions),
          positions)
    build('p_fisher', 'Fisher exact, recomputed',
          f'Figure15b_significance_logo_{tag}_fisher',
          note('p_fisher', positions), positions)

# ---------------------------------------------- fold change vs significance
enr = res[res.direction > 0].copy()
enr = enr[enr.p_chi2 < 0.05]
enr['Stars'] = enr.p_chi2.map(pg.stars)
cols = ['position', 'label', 'n_sub', 'n_ctrl', 'fold_change', 'p_chi2',
        'p_fisher', 'Stars']
fmt = {'fold_change': '{:.2f}'.format, 'p_chi2': '{:.2e}'.format,
       'p_fisher': '{:.2e}'.format}

print("\nfold change and significance rank residues differently:")
print("\n  top 5 by FOLD CHANGE:")
print(enr.nlargest(5, 'fold_change')[cols].to_string(index=False, formatters=fmt))
print("\n  top 5 by SIGNIFICANCE:")
print(enr.nsmallest(5, 'p_chi2')[cols].to_string(index=False, formatters=fmt))
