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

Heights are now EXACT -log10(p), not the * / ** / *** bin representatives the
first version of this figure had to use: pg_motif_data.py reads the workbook
sheet that carries a p value per cell. That module also documents the arginine
correction -- the earlier source file's R row was the Basic-class row, so R was
drawn far too short here.

Two versions are produced:
  Figure15   chi-square, the workbook's own test.
  Figure15b  Fisher exact, recomputed from the same 2 x 2 tables. Chi-square
             needs expected counts >= 5 and most cells here expect fewer than
             2 substrates, so this is the version to trust.

Substrates (n = 20) vs. non-substrate controls carrying the same motif
(n = 193).
"""
import os

import numpy as np
import pandas as pd
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
CATEGORY_MEMBERS, CAT_COLORS = pg.CATEGORY_MEMBERS, pg.CAT_COLORS
LEGEND_ORDER = ['Basic', 'Acidic', 'Aromatic', 'Aliphatic', 'Hydrophobic', 'Polar']
AA_COLOR_SCHEME = {a: CAT_COLORS[c] for c, mem in CATEGORY_MEMBERS.items() for a in mem}

cells = pg.load()
DIR = pg.matrix('direction', cells=cells)

# ------------------------------------------------------------------ style
mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'Liberation Sans', 'DejaVu Sans'],
    'font.size': 8, 'axes.linewidth': 0.8,
    'svg.fonttype': 'none', 'pdf.fonttype': 42,
    'figure.dpi': 600, 'savefig.dpi': 600, 'savefig.bbox': 'tight',
    'axes.spines.top': False, 'axes.spines.right': False,
})


def build(pcol, test_label, stem, q_note):
    P = pg.matrix(pcol, cells=cells)
    Q = pg.matrix(pcol.replace('p_', 'q_'), cells=cells)

    # Height = -log10(p), for residues ENRICHED in the substrates and reaching
    # p < 0.05. The threshold matters: without it all 20 residues are drawn at
    # every position, most of them at p ~ 0.3-0.9, and the stacks pile up to
    # -log10 p ~ 8 out of noise alone. Drawing only what passes is the pLogo
    # convention and is what the first version of this figure did (via its
    # star bins); the heights are now exact rather than bin representatives.
    with np.errstate(divide='ignore', invalid='ignore'):
        H = -np.log10(P.astype(float))
    H = H.where(np.isfinite(H) & DIR.gt(0) & P.lt(0.05), 0.0).clip(lower=0.0)
    H = H.T.set_axis(POS)

    fig, ax = plt.subplots(figsize=(10.5, 2.9))
    logomaker.Logo(H, ax=ax, color_scheme=AA_COLOR_SCHEME,
                   font_name='Helvetica', vpad=0.0, width=0.95,
                   stack_order='big_on_top')
    for patch in ax.patches:
        patch.set_edgecolor('black')
        patch.set_linewidth(0.35)

    for y, lab in [(1.30103, 'p = 0.05'), (2.0, 'p = 0.01'), (3.0, 'p = 0.001')]:
        ax.axhline(y=y, color='#999', linewidth=0.5, linestyle=(0, (3, 2)), zorder=1)
        ax.text(POS[-1] + 0.55, y, lab, fontsize=5, color='#666', va='center', ha='left')

    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.yaxis.grid(False)
    ax.xaxis.grid(False)
    ax.set_axisbelow(True)
    ax.set_ylabel('$-$log$_{10}$ $p$   (stacked)', fontsize=8)
    ax.set_xlabel('Position', fontsize=8, labelpad=2)
    ax.set_xticks(POS)
    ax.set_xticklabels(POS, fontsize=6.5)
    ax.tick_params(axis='both', labelsize=6.5, length=2.5, width=0.7, pad=2)
    ax.set_title(f'Significance Logo — height is $-$log$_{{10}}$ $p$, '
                 f'not fold change ({test_label})',
                 fontweight='bold', fontsize=8.5, pad=4)
    ax.set_ylim(-0.05, float(H.sum(axis=1).max()) * 1.12)
    ax.set_xlim(POS[0] - 0.6, POS[-1] + 0.6)
    ax.text(0.5, -0.15,
            f'n = {N_SUB} P/G-D/E/T substrates  vs  n = {N_CTRL} non-substrate '
            f'controls.   {q_note}',
            transform=ax.transAxes, ha='center', fontsize=5.8, color='#666',
            style='italic')

    ax.legend(handles=[mpatches.Patch(color=CAT_COLORS[c],
                                      label=f'{c}  ({" ".join(CATEGORY_MEMBERS[c])})')
                       for c in LEGEND_ORDER],
              fontsize=5.5, frameon=False, loc='upper left', bbox_to_anchor=(1.10, 1.0),
              handlelength=0.8, handletextpad=0.3, borderpad=0.2)
    plt.tight_layout(pad=0.3)

    for ext in ('png', 'pdf', 'svg'):
        p = os.path.join(OUTPUT_DIR, f'{stem}.{ext}')
        fig.savefig(p, format=ext)
        print(f"  saved: {os.path.basename(p)}  ({os.path.getsize(p)/1024:.1f} KB)")
    plt.close(fig)
    return P, Q


res = cells[(cells.kind == 'residue') & cells.defined]
n_test = len(res)


def note(pcol):
    q = pcol.replace('p_', 'q_')
    n_hit, n_fdr = int((res[pcol] < 0.05).sum()), int((res[q] < 0.05).sum())
    survivors = ', '.join(f'{r.label}{r.position}' for _, r in
                          res[res[q] < 0.05].sort_values(q).iterrows())
    return (f'Heights are uncorrected $p$: {n_hit} of {n_test} cells reach '
            f'$p$ < 0.05 against {0.05 * n_test:.0f} expected by chance, and '
            + (f'only {survivors} survive'
               f'{"s" if n_fdr == 1 else ""} BH-FDR.' if n_fdr else
               'none survives BH-FDR.'))


build('p_chi2', 'chi-square, as reported in the workbook',
      'Figure15_significance_logo_pos4_24', note('p_chi2'))
build('p_fisher', 'Fisher exact, recomputed',
      'Figure15b_significance_logo_pos4_24_fisher', note('p_fisher'))

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
