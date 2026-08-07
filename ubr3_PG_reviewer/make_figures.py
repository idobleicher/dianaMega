#!/usr/bin/env python
"""Publication figures for the UBR3 N24mer Pro/Gly analysis.

Renders 5 multi-panel figures as 300-dpi PNG + vector PDF into figures/.
Light/print mode only - these are manuscript figures, not screen artifacts.
"""
import os
import warnings

import logomaker
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from scipy import stats

import ubr3_core as U

warnings.filterwarnings('ignore')

# ---------------------------------------------------------------- palette
SURFACE = '#fcfcfb'
INK = '#0b0b0b'
INK2 = '#52514e'
MUTED = '#898781'
GRID = '#e1e0d9'
AXIS = '#c3c2b7'

BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN = (
    '#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300')
RED, VIOLET = '#e34948', '#4a3aa7'

# validated adjacent-pair order, mapped onto CLASS_ORDER
CLASS_COLOR = dict(zip(U.CLASS_ORDER, [BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN]))
AA_COLOR = {a: CLASS_COLOR[U.AA_CLASS[a]] for a in U.AAS}

SEQ_BLUE = LinearSegmentedColormap.from_list(
    'seqblue', ['#fcfcfb', '#cde2fb', '#9ec5f4', '#5598e7', '#2a78d6', '#256abf', '#0d366b'])
DIVERGE = LinearSegmentedColormap.from_list(
    'bluered', ['#0d366b', '#2a78d6', '#9ec5f4', '#f0efec', '#f6b3b3', '#e34948', '#8f1f1f'])

mpl.rcParams.update({
    'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE,
    'font.family': 'sans-serif',
    'font.sans-serif': ['Segoe UI', 'DejaVu Sans', 'Arial'],
    'font.size': 9, 'axes.titlesize': 11, 'axes.labelsize': 9.5,
    'axes.edgecolor': AXIS, 'axes.labelcolor': INK, 'axes.linewidth': 0.9,
    'xtick.color': MUTED, 'ytick.color': MUTED, 'text.color': INK,
    'xtick.labelsize': 8.5, 'ytick.labelsize': 8.5,
    'legend.frameon': False, 'legend.fontsize': 8.5,
    'grid.color': GRID, 'grid.linewidth': 0.7,
    'axes.spines.top': False, 'axes.spines.right': False,
    'figure.dpi': 110, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})


def panel_tag(ax, letter, xoff=-46, yoff=30):
    """Panel letter, placed in point offsets so it never collides with the title."""
    ax.annotate(letter, xy=(0, 1), xycoords='axes fraction',
                xytext=(xoff, yoff), textcoords='offset points',
                fontsize=14, fontweight='bold', va='bottom', ha='left', color=INK)


def title(ax, text, sub=None):
    """Bold title above a smaller subtitle; both in point offsets from the axes top."""
    ax.set_title(text, loc='left', fontweight='bold', color=INK,
                 pad=26 if sub else 10)
    if sub:
        ax.annotate(sub, xy=(0, 1), xycoords='axes fraction',
                    xytext=(0, 7), textcoords='offset points',
                    fontsize=8.3, color=INK2, va='bottom', ha='left')


def save(fig, name):
    for ext in ('png', 'pdf'):
        fig.savefig(os.path.join(U.FIGS, f'{name}.{ext}'))
    plt.close(fig)
    print('  wrote', name + '.png / .pdf')


def stars(q):
    return '***' if q < 1e-3 else '**' if q < 1e-2 else '*' if q < 0.05 else 'ns'


# ================================================================ FIGURE 1
def figure1(lib, hit, enr):
    fig = plt.figure(figsize=(13.6, 10.6))
    gs = fig.add_gridspec(2, 2, hspace=0.46, wspace=0.26,
                          left=0.07, right=0.975, top=0.855, bottom=0.065)

    # ---- A: position-2 residue composition -------------------------------
    ax = fig.add_subplot(gs[0, 0])
    lf = 100 * lib.pos2.value_counts(normalize=True).reindex(U.AAS).fillna(0)
    hf = 100 * hit.pos2.value_counts(normalize=True).reindex(U.AAS).fillna(0)
    order = lf.sort_values(ascending=False).index.tolist()
    y = np.arange(len(order))
    ax.barh(y + 0.21, lf[order], height=0.4, color=MUTED, label=f'Library (n={len(lib):,})',
            zorder=3, edgecolor=SURFACE, linewidth=0.8)
    ax.barh(y - 0.21, hf[order], height=0.4, color=BLUE, label=f'UBR3 substrates (n={len(hit)})',
            zorder=3, edgecolor=SURFACE, linewidth=0.8)
    for i, a in enumerate(order):
        if a in 'PG':
            ax.add_patch(Rectangle((0, i - 0.5), 100, 1, color=YELLOW, alpha=0.16, zorder=0))
        if hf[a] > 0:
            ax.text(hf[a] + 0.5, i - 0.21, f'{hf[a]:.0f}%', va='center', fontsize=7.4, color=BLUE)
    ax.set_yticks(y)
    ax.set_yticklabels(order, fontsize=8.5)
    for t in ax.get_yticklabels():
        if t.get_text() in 'PG':
            t.set_fontweight('bold')
            t.set_color(INK)
    ax.invert_yaxis()
    ax.set_xlabel('Peptides with this residue at position 2 (%)')
    ax.set_ylabel('Residue at position 2')
    ax.grid(axis='x', zorder=0)
    ax.legend(loc='lower right')
    ax.set_xlim(0, max(lf.max(), hf.max()) * 1.14)
    title(ax, 'Position 2 is Pro/Gly in half the substrates',
          'Shaded rows = Pro and Gly, the residues exposed after initiator-Met excision')
    panel_tag(ax, 'A')

    # ---- B: Fisher enrichment at pos2 ------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    e = enr.dropna(subset=['log2_fold']).sort_values('log2_fold')
    cols = [BLUE if v > 0 else RED for v in e.log2_fold]
    yy = np.arange(len(e))
    ax.barh(yy, e.log2_fold, color=cols, height=0.66, zorder=3,
            edgecolor=SURFACE, linewidth=0.8)
    ax.axvline(0, color=AXIS, lw=1)
    for i, (_, r) in enumerate(e.iterrows()):
        off = 0.09 if r.log2_fold > 0 else -0.09
        ha = 'left' if r.log2_fold > 0 else 'right'
        lab = f'{stars(r["q_value_BH"])}  n={r.n_hits}' if r['q_value_BH'] < 0.05 else f'n={r.n_hits}'
        ax.text(r.log2_fold + off, i, lab, va='center', ha=ha, fontsize=7.4,
                color=INK if r['q_value_BH'] < 0.05 else MUTED,
                fontweight='bold' if r['q_value_BH'] < 0.05 else 'normal')
    ax.set_yticks(yy)
    ax.set_yticklabels(e.residue_at_pos2, fontsize=8.5)
    ax.set_xlabel('log$_2$ fold enrichment (substrates / library)')
    ax.set_ylabel('Residue at position 2')
    ax.grid(axis='x', zorder=0)
    ax.set_xlim(-2.9, 3.2)
    title(ax, 'Pro and Gly are the only enriched residues',
          'Two-sided Fisher exact, Benjamini-Hochberg FDR;  *** q<0.001, * q<0.05')
    panel_tag(ax, 'B')

    # ---- C: funnel --------------------------------------------------------
    ax = fig.add_subplot(gs[1, 0])
    steps = [('All library peptides', len(lib), MUTED),
             ('Pro or Gly at position 2', int(lib.is_PG.sum()), AQUA),
             ('[P/G]-[E/D] motif', int(lib.is_PG_ED.sum()), YELLOW),
             ('[P/G]-[E/D] and UBR3-stabilised',
              int(lib[lib.is_PG_ED].is_UBR3_substrate.sum()), ORANGE)]
    ys = np.arange(len(steps))[::-1]
    for yv, (lab, n, c) in zip(ys, steps):
        ax.barh(yv, np.log10(n), color=c, height=0.62, zorder=3,
                edgecolor=SURFACE, linewidth=0.9)
        ax.text(np.log10(n) + 0.06, yv, f'{n:,}', va='center', fontweight='bold',
                fontsize=10.5, color=INK)
        ax.text(0.06, yv + 0.40, lab, va='bottom', ha='left', fontsize=8.6, color=INK2)
    for i in range(len(steps) - 1):
        a, b = steps[i][1], steps[i + 1][1]
        ax.text(np.log10(b) + 0.62, ys[i] - 0.5, f'{100 * b / a:.1f}% retained',
                va='center', fontsize=7.6, color=MUTED, style='italic')
    ax.set_yticks([])
    ax.set_xlabel('Number of peptides (log$_{10}$ scale)')
    ax.set_xlim(0, np.log10(len(lib)) * 1.22)
    ax.grid(axis='x', zorder=0)
    ax.spines['left'].set_visible(False)
    title(ax, 'From 16,514 peptides to 16 motif-bearing substrates',
          'Each bar is a nested subset of the one above it')
    panel_tag(ax, 'C')

    # ---- D: stabilisation rate by motif class ----------------------------
    ax = fig.add_subplot(gs[1, 1])
    groups = ['[P/G]-[E/D]', '[P/G]-other', 'non-P/G']
    cols3 = [ORANGE, AQUA, MUTED]
    rate, ns, ks = [], [], []
    for g in groups:
        s = lib[lib.motif_class == g]
        rate.append(100 * s.is_UBR3_substrate.mean())
        ns.append(len(s))
        ks.append(int(s.is_UBR3_substrate.sum()))
    x = np.arange(3)
    ax.bar(x, rate, color=cols3, width=0.56, zorder=3, edgecolor=SURFACE, linewidth=1.2)
    tops = []
    for i, (r, n, k) in enumerate(zip(rate, ns, ks)):
        # Wilson 95% CI
        z, p, nn = 1.96, k / n, n
        c = (p + z * z / (2 * nn)) / (1 + z * z / nn)
        h = z * np.sqrt(p * (1 - p) / nn + z * z / (4 * nn * nn)) / (1 + z * z / nn)
        lo, hi = 100 * (c - h), 100 * (c + h)
        tops.append(hi)
        ax.plot([i, i], [lo, hi], color=INK2, lw=1.3, zorder=4)
        ax.text(i + 0.34, r, f'{r:.2f}%', ha='left', va='center',
                fontweight='bold', fontsize=10.5)
        ax.text(i, -0.75, f'{k} / {n:,}', ha='center', fontsize=8.2, color=INK2)
    ax.set_xticks(x)
    ax.set_xticklabels(groups, fontsize=9.5)
    ax.set_ylabel('UBR3-stabilised peptides (% of group)')
    top = max(tops)
    ax.set_ylim(-1.4, top * 1.30)
    ax.set_xlim(-0.6, 2.95)
    ax.grid(axis='y', zorder=0)
    yb = top * 1.10
    ax.plot([0, 0, 2, 2], [yb - 0.3, yb, yb, yb - 0.3], color=INK2, lw=1.1, zorder=5)
    ax.text(1, yb + 0.18, 'Fisher exact  OR = 50.4,  p = 4.1 $\\times$ 10$^{-20}$',
            ha='center', fontsize=8.4, color=INK)
    title(ax, 'The acidic residue at position 3 is what matters',
          'Error bars are Wilson 95% confidence intervals;  counts below each bar')
    panel_tag(ax, 'D')

    fig.suptitle('Figure 1 | Pro/Gly N-termini and the [P/G]-[E/D] motif in the N24mer stability screen',
                 x=0.07, y=0.968, ha='left', fontsize=14.5, fontweight='bold', color=INK)
    fig.text(0.07, 0.938,
             'GPS screen of 16,514 human N-terminal 24-mers in control-KO vs UBR3-KO cells. '
             'Position 2 becomes the N-terminal residue after initiator-Met excision.',
             ha='left', fontsize=9.2, color=INK2)
    save(fig, 'Figure1_PG_landscape')


# ================================================================ FIGURE 2
def figure2(lib, hit):
    pged = lib[lib.motif_class == '[P/G]-[E/D]'].sort_values('mean_dPSI', ascending=False)
    fig = plt.figure(figsize=(13.6, 10.6))
    gs = fig.add_gridspec(2, 2, hspace=0.50, wspace=0.26,
                          left=0.07, right=0.975, top=0.855, bottom=0.065)

    # ---- A: ranked dPSI of every [P/G]-[E/D] peptide ---------------------
    ax = fig.add_subplot(gs[0, 0])
    xr = np.arange(1, len(pged) + 1)
    st = pged.is_UBR3_substrate.values
    ax.bar(xr[~st], pged.mean_dPSI.values[~st], color=MUTED, width=1.0, zorder=3,
           label=f'Not stabilised (n={int((~st).sum())})')
    ax.bar(xr[st], pged.mean_dPSI.values[st], color=ORANGE, width=1.0, zorder=4,
           label=f'UBR3 substrate (n={int(st.sum())})')
    ax.axhline(0, color=AXIS, lw=1)
    ax.set_xlabel('[P/G]-[E/D] peptides, ranked by mean $\\Delta$PSI')
    ax.set_ylabel('Mean $\\Delta$PSI  (UBR3 KO $-$ control)')
    ax.grid(axis='y', zorder=0)
    ax.legend(loc='upper right')
    ax.set_xlim(0, len(pged) + 1)
    title(ax, '163 of 179 motif-bearing peptides are not stabilised',
          'One bar per peptide; the 16 substrates are orange')
    panel_tag(ax, 'A')

    # ---- B: ECDF by motif class ------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    for g, c in [('[P/G]-[E/D]', ORANGE), ('[P/G]-other', AQUA), ('non-P/G', MUTED)]:
        v = np.sort(lib[lib.motif_class == g].mean_dPSI.values)
        ax.plot(v, np.arange(1, len(v) + 1) / len(v), color=c, lw=2,
                label=f'{g}  (n={len(v):,})', zorder=3)
    ax.set_xlabel('Mean $\\Delta$PSI')
    ax.set_ylabel('Cumulative fraction of peptides')
    ax.set_xlim(-0.8, 1.2)
    ax.grid(zorder=0)
    ax.legend(loc='upper left')
    p = stats.mannwhitneyu(lib[lib.motif_class == '[P/G]-[E/D]'].mean_dPSI,
                           lib[lib.motif_class == 'non-P/G'].mean_dPSI)[1]
    ax.text(0.98, 0.33, f'[P/G]-[E/D] vs non-P/G\nMann-Whitney p = {p:.1e}',
            transform=ax.transAxes, ha='right', fontsize=8.4, color=INK)
    title(ax, 'The whole motif class shifts, but only slightly',
          'The [P/G]-[E/D] curve sits right of the other two, but the overlap is large')
    panel_tag(ax, 'B')

    # ---- C: stacked composition ------------------------------------------
    ax = fig.add_subplot(gs[1, 0])
    groups = ['[P/G]-[E/D]', '[P/G]-other', 'non-P/G']
    y = np.arange(3)[::-1]
    for yv, g in zip(y, groups):
        s = lib[lib.motif_class == g]
        k, n = int(s.is_UBR3_substrate.sum()), len(s)
        pk = 100 * k / n
        ax.barh(yv, pk, color=ORANGE, height=0.5, zorder=3, edgecolor=SURFACE, linewidth=1.4)
        ax.barh(yv, 100 - pk, left=pk, color='#dedcd4', height=0.5, zorder=3,
                edgecolor=SURFACE, linewidth=1.4)
        ax.text(101, yv, f'{n - k:,} not stabilised', va='center', fontsize=8.6, color=INK2)
        ax.text(pk / 2 if pk > 14 else pk + 1.2, yv, f'{pk:.2f}%',
                va='center', ha='center' if pk > 14 else 'left',
                fontsize=9.4, fontweight='bold',
                color='white' if pk > 14 else ORANGE, zorder=5)
    ax.set_yticks(y)
    ax.set_yticklabels(groups, fontsize=10)
    ax.set_xlim(0, 100)
    ax.set_xlabel('Percent of peptides in the group')
    ax.grid(axis='x', zorder=0)
    ax.spines['left'].set_visible(False)
    ax.set_ylim(-0.9, 2.6)
    ax.legend(handles=[Patch(color=ORANGE, label='UBR3-stabilised'),
                       Patch(color='#dedcd4', label='Not stabilised')],
              loc='lower center', ncol=2, columnspacing=1.6)
    title(ax, 'Motif enrichment is 17-fold, but the motif is not sufficient',
          'Even in the best group, 91% of motif-bearing peptides are unaffected by UBR3 loss')
    panel_tag(ax, 'C')

    # ---- D: replicate concordance ----------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    bg = lib[~lib.is_PG_ED]
    ax.scatter(bg[U.D1], bg[U.D2], s=3, color='#e6e4dc', linewidths=0, zorder=1,
               label=f'Other library peptides (n={len(bg):,})')
    nn = pged[~pged.is_UBR3_substrate]
    ss = pged[pged.is_UBR3_substrate]
    ax.scatter(nn[U.D1], nn[U.D2], s=26, color=MUTED, linewidths=0.6, edgecolors=SURFACE,
               zorder=3, label=f'[P/G]-[E/D], not stabilised (n={len(nn)})')
    ax.scatter(ss[U.D1], ss[U.D2], s=52, color=ORANGE, linewidths=0.9, edgecolors=SURFACE,
               zorder=4, label=f'[P/G]-[E/D], UBR3 substrate (n={len(ss)})')
    # label the top hits, nudging each away from the ones already placed
    placed = []
    for _, r in ss.nlargest(6, 'mean_dPSI').iterrows():
        x0, y0 = r[U.D1], r[U.D2]
        for dx, dy in [(8, 4), (8, -11), (-8, 8), (-8, -13), (8, 13), (-8, 18)]:
            if all(abs(x0 + dx / 90 - px) > 0.16 or abs(y0 + dy / 90 - py) > 0.075
                   for px, py in placed):
                break
        placed.append((x0 + dx / 90, y0 + dy / 90))
        ax.annotate(r.Gene_ID, (x0, y0), textcoords='offset points', xytext=(dx, dy),
                    fontsize=8, color=INK, fontweight='bold',
                    ha='left' if dx > 0 else 'right')
    ax.axhline(0, color=AXIS, lw=0.8)
    ax.axvline(0, color=AXIS, lw=0.8)
    ax.set_xlabel('$\\Delta$PSI replicate 1')
    ax.set_ylabel('$\\Delta$PSI replicate 2')
    ax.set_xlim(-1.0, 1.85)
    ax.set_ylim(-1.0, 1.85)
    ax.grid(zorder=0)
    ax.legend(loc='lower right', markerscale=1.4)
    r = stats.pearsonr(pged[U.D1], pged[U.D2])
    ax.text(0.03, 0.96, f'[P/G]-[E/D] replicate r = {r[0]:.2f}', transform=ax.transAxes,
            fontsize=8.6, va='top', color=INK)
    title(ax, 'Stabilisation is reproducible, not replicate noise',
          'Both replicates plotted; the 16 substrates cluster together in the top-right corner')
    panel_tag(ax, 'D')

    fig.suptitle('Figure 2 | The [P/G]-[E/D] motif is necessary but not sufficient for UBR3 regulation',
                 x=0.07, y=0.968, ha='left', fontsize=14.5, fontweight='bold', color=INK)
    fig.text(0.07, 0.938,
             'Only 16 of the 179 library peptides carrying a Pro/Gly-initiated acidic motif are '
             'stabilised when UBR3 is lost - a 27-fold enrichment over background, but still a small minority.',
             ha='left', fontsize=9.2, color=INK2)
    save(fig, 'Figure2_motif_necessary_not_sufficient')


# ================================================================ FIGURE 3
def figure3(lib, hit, sig):
    """Logos over positions 2-24.  Position 1 is Met in 100% of both sets, so it
    carries no information and is dropped - keeping it flattens every other column."""
    hs, ls = list(hit.peptide_24mer), list(lib.peptide_24mer)
    lr, qv, signed = sig
    POS = list(range(2, 25))

    def trim(m):
        m = m.loc[POS].copy()
        m.index = range(len(POS))          # logomaker needs 0-based contiguous rows
        return m

    info = trim(logomaker.transform_matrix(U.position_matrix(hs), from_type='counts',
                                           to_type='information'))
    linfo = trim(logomaker.transform_matrix(U.position_matrix(ls), from_type='counts',
                                            to_type='information'))
    # statistical enrichment logo: height = signed -log10(p), only p<0.05 drawn
    praw = trim(signed)
    praw = praw.where(praw.abs() >= -np.log10(0.05), 0.0)

    fig = plt.figure(figsize=(15.4, 11.9))
    gs = fig.add_gridspec(3, 1, hspace=0.62, left=0.065, right=0.845, top=0.858, bottom=0.055)

    def dress(ax, ylab, sub, tag, ttl, ymax=None):
        ax.set_xticks(range(len(POS)))
        ax.set_xticklabels(POS, fontsize=8.5)
        ax.set_xlabel('Position in the encoded 24-mer   '
                      '(position 1 is the invariant initiator Met and is omitted)')
        ax.set_ylabel(ylab)
        ax.grid(axis='y', zorder=0, alpha=0.55)
        ax.set_axisbelow(True)
        ax.axvspan(-0.5, 1.5, color=YELLOW, alpha=0.16, zorder=0)
        title(ax, ttl, sub)
        panel_tag(ax, tag, xoff=-52)

    ax = fig.add_subplot(gs[0])
    logomaker.Logo(info, ax=ax, color_scheme=AA_COLOR, show_spines=False,
                   stack_order='big_on_top')
    top = max(info.sum(1).max(), linfo.sum(1).max()) * 1.18
    ax.set_ylim(0, top)
    dress(ax, 'Information content (bits)',
          f'Letter height = information content; the shaded block is positions 2-3 (n={len(hs)} substrates)',
          'A', 'Sequence logo of the 54 candidate UBR3 substrates')
    ax.annotate('Pro/Gly at position 2,\nAsp/Glu at position 3',
                xy=(0.9, info.sum(1).max() * 0.99), xytext=(3.2, top * 0.80),
                fontsize=9, fontweight='bold', color=INK,
                arrowprops=dict(arrowstyle='->', color=INK2, lw=1.1))

    ax = fig.add_subplot(gs[1])
    logomaker.Logo(linfo, ax=ax, color_scheme=AA_COLOR, show_spines=False,
                   stack_order='big_on_top')
    ax.set_ylim(0, top)
    dress(ax, 'Information content (bits)',
          f'Background from all {len(ls):,} library peptides, on the same y-scale as panel A - essentially flat',
          'B', 'Sequence logo of the full N24mer library (background)')

    ax = fig.add_subplot(gs[2])
    logomaker.Logo(praw, ax=ax, color_scheme=AA_COLOR, show_spines=False,
                   stack_order='big_on_top', flip_below=True)
    ax.axhline(0, color=AXIS, lw=1)
    for lv, lab, st in [(-np.log10(0.05), 'p = 0.05', ':'),
                        (-np.log10(0.05 / (23 * 20)), 'Bonferroni', '--')]:
        for s in (1, -1):
            ax.axhline(s * lv, color=RED, lw=0.9, ls=st, zorder=5)
        ax.text(len(POS) - 0.4, lv + 0.12, lab, fontsize=7.6, color=RED, ha='right')
    dress(ax, 'signed $-$log$_{10}$ $p$  (Fisher exact)',
          'Only residues with p < 0.05 are drawn. Above the line = enriched in substrates, below = depleted.',
          'C', 'Statistical enrichment logo: what actually distinguishes substrates from the library')

    handles = [Patch(color=CLASS_COLOR[c], label=f'{c}  ({", ".join(U.CLASS_MEMBERS[c])})')
               for c in U.CLASS_ORDER]
    fig.legend(handles=handles, loc='center left', bbox_to_anchor=(0.855, 0.5),
               title='Amino-acid class', title_fontsize=9.5, labelspacing=0.95,
               alignment='left')

    fig.suptitle('Figure 3 | Sequence logos of the UBR3 substrate N-termini, coloured by amino-acid class',
                 x=0.065, y=0.968, ha='left', fontsize=14.5, fontweight='bold', color=INK)
    fig.text(0.065, 0.938,
             'The constraint is concentrated at positions 2 and 3 - Pro/Gly followed by Asp/Glu. '
             'Only one cell further along the peptide (Arg at position 7) survives FDR correction.',
             ha='left', fontsize=9.2, color=INK2)
    save(fig, 'Figure3_sequence_logos')


# ================================================================ FIGURE 4
def figure4(lib, hit):
    hs, ls = list(hit.peptide_24mer), list(lib.peptide_24mer)
    ch, cl = U.class_matrix(hs), U.class_matrix(ls)
    clr, clq, _ = U.enrichment_test(hs, ls, groups=U.CLASS_MEMBERS)

    fig = plt.figure(figsize=(14.2, 11.7))
    gs = fig.add_gridspec(3, 2, hspace=0.62, wspace=0.24, height_ratios=[1, 1, 1.05],
                          left=0.07, right=0.98, top=0.872, bottom=0.062)

    # ---- A/B stacked class composition ------------------------------------
    for j, (mat, lab, tag) in enumerate([(ch, f'54 UBR3 substrates', 'A'),
                                         (cl, f'{len(ls):,} library peptides', 'B')]):
        ax = fig.add_subplot(gs[j, :])
        bottom = np.zeros(24)
        xs = np.arange(1, 25)
        for c in U.CLASS_ORDER:
            v = 100 * mat[c].values
            ax.bar(xs, v, bottom=bottom, color=CLASS_COLOR[c], width=0.76,
                   label=f'{c} ({",".join(U.CLASS_MEMBERS[c])})',
                   zorder=3, edgecolor=SURFACE, linewidth=1.1)
            bottom += v
        ax.set_xticks(xs)
        ax.set_xticklabels(xs, fontsize=8.3)
        ax.set_ylim(0, 100)
        ax.set_xlim(0.4, 24.6)
        ax.set_ylabel('Composition (%)')
        ax.set_xlabel('Position in the 24-mer')
        ax.set_axisbelow(True)
        if j == 0:
            ax.legend(ncol=6, loc='lower center', bbox_to_anchor=(0.5, 1.20),
                      columnspacing=1.4, handlelength=1.4, fontsize=8.2)
            title(ax, f'Amino-acid class composition at each position - {lab}',
                  'Position 2 of the substrates is 48% Special (Pro/Gly); position 3 is 37% Acidic')
        else:
            title(ax, f'Amino-acid class composition at each position - {lab}',
                  'The background is flat from position 4 onward')
        panel_tag(ax, tag, xoff=-52)

    # ---- C heatmap of class enrichment (FDR-masked) ------------------------
    ax = fig.add_subplot(gs[2, 0])
    m = clr[U.CLASS_ORDER].T.values
    q = clq[U.CLASS_ORDER].T.values
    v = np.nanmax(np.abs(np.where(q < 0.05, m, 0))) or np.nanmax(np.abs(m))
    shown = np.where(q < 0.05, m, np.nan)
    ax.imshow(np.zeros_like(m), cmap=mpl.colors.ListedColormap(['#f2f1ec']), aspect='auto')
    im = ax.imshow(shown, cmap=DIVERGE, norm=TwoSlopeNorm(0, -v, v), aspect='auto')
    ax.set_yticks(range(len(U.CLASS_ORDER)))
    ax.set_yticklabels(U.CLASS_ORDER, fontsize=9)
    for t, c in zip(ax.get_yticklabels(), U.CLASS_ORDER):
        t.set_color(CLASS_COLOR[c])
        t.set_fontweight('bold')
    ax.set_xticks(range(24))
    ax.set_xticklabels(range(1, 25), fontsize=7.6)
    ax.set_xlabel('Position in the 24-mer')
    for i in range(len(U.CLASS_ORDER)):
        for k in range(24):
            if q[i, k] < 0.05:
                ax.text(k, i, f'{m[i, k]:+.1f}', ha='center', va='center', fontsize=6.4,
                        color='white' if abs(m[i, k]) > v * 0.55 else INK, fontweight='bold')
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label('log$_2$ (substrate / library)', fontsize=8.5)
    cb.ax.tick_params(labelsize=7.5)
    ax.set_xticks(np.arange(-0.5, 24, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 6, 1), minor=True)
    ax.grid(which='minor', color=SURFACE, lw=1.1)
    ax.tick_params(which='minor', length=0)
    title(ax, 'Class enrichment by position (FDR-masked)',
          'Beige = not significant at q < 0.05. Positions 2 and 3 dominate; Basic at position 7 '
          'is the only other survivor.')
    panel_tag(ax, 'C', xoff=-70)

    # ---- D overall class composition --------------------------------------
    ax = fig.add_subplot(gs[2, 1])
    ph = np.array([np.mean([1 for p in hs for a in p if U.AA_CLASS[a] == c]) * 0 +
                   sum(U.AA_CLASS[a] == c for p in hs for a in p) / (len(hs) * 24)
                   for c in U.CLASS_ORDER]) * 100
    pl = np.array([sum(U.AA_CLASS[a] == c for p in ls for a in p) / (len(ls) * 24)
                   for c in U.CLASS_ORDER]) * 100
    x = np.arange(len(U.CLASS_ORDER))
    ax.bar(x - 0.2, ph, width=0.38, color=[CLASS_COLOR[c] for c in U.CLASS_ORDER],
           zorder=3, edgecolor=SURFACE, linewidth=1.1)
    ax.bar(x + 0.2, pl, width=0.38, color=MUTED, zorder=3, edgecolor=SURFACE, linewidth=1.1)
    for i in range(len(x)):
        ax.text(x[i] - 0.2, ph[i] + 0.5, f'{ph[i]:.1f}', ha='center', fontsize=7.6, color=INK)
        ax.text(x[i] + 0.2, pl[i] + 0.5, f'{pl[i]:.1f}', ha='center', fontsize=7.6, color=MUTED)
    ax.set_xticks(x)
    ax.set_xticklabels(U.CLASS_ORDER, rotation=22, ha='right', fontsize=8.8)
    ax.set_ylabel('Residues of this class (% of all 24 positions)')
    ax.grid(axis='y', zorder=0)
    ax.set_ylim(0, max(ph.max(), pl.max()) * 1.24)
    ax.legend(handles=[Patch(facecolor='white', edgecolor=INK,
                             label='UBR3 substrates (bar coloured by class)'),
                       Patch(facecolor=MUTED, label='Library')],
              loc='upper left')
    title(ax, 'Whole-peptide composition is unchanged',
          'Averaged over all 24 positions the substrates look like the library')
    panel_tag(ax, 'D', xoff=-62)

    fig.suptitle('Figure 4 | Chemical-class architecture of the substrate N-termini across all 24 positions',
                 x=0.07, y=0.982, ha='left', fontsize=14.5, fontweight='bold', color=INK)
    fig.text(0.07, 0.955,
             'Every residue of every peptide is assigned to one of six chemical classes. '
             'Only positions 2, 3 and 7 differ from the library after FDR correction.',
             ha='left', fontsize=9.2, color=INK2)
    save(fig, 'Figure4_AA_class_analysis')


# ================================================================ FIGURE 5
def figure5(lib, hit, sig):
    hs, ls = list(hit.peptide_24mer), list(lib.peptide_24mer)
    fh = U.freq_matrix(hs)
    lr, qv, _ = sig
    aa_by_class = [a for c in U.CLASS_ORDER for a in U.CLASS_MEMBERS[c]]

    fig = plt.figure(figsize=(15.0, 10.1))
    gs = fig.add_gridspec(1, 2, wspace=0.24, left=0.075, right=0.965, top=0.825, bottom=0.09)

    def frame(ax, tag, ttl, sub):
        ax.set_yticks(range(len(aa_by_class)))
        ax.set_yticklabels(aa_by_class, fontsize=8.5, fontweight='bold')
        for t, a in zip(ax.get_yticklabels(), aa_by_class):
            t.set_color(CLASS_COLOR[U.AA_CLASS[a]])
        ax.set_xticks(range(24))
        ax.set_xticklabels(range(1, 25), fontsize=7.6)
        ax.set_xlabel('Position in the 24-mer')
        ax.set_ylabel('Amino acid (grouped by chemical class)')
        acc = 0
        for c in U.CLASS_ORDER[:-1]:
            acc += len(U.CLASS_MEMBERS[c])
            ax.axhline(acc - 0.5, color=INK, lw=1.4)
        # box the motif region (positions 2-3)
        ax.add_patch(Rectangle((0.5, -0.5), 2, len(aa_by_class), fill=False,
                               edgecolor=INK, lw=1.8, zorder=6))
        title(ax, ttl, sub)
        panel_tag(ax, tag, xoff=-64, yoff=34)

    # ---- A: raw frequency --------------------------------------------------
    ax = fig.add_subplot(gs[0])
    m = (100 * fh)[aa_by_class].T.values
    im = ax.imshow(m, cmap=SEQ_BLUE, vmin=0, vmax=28, aspect='auto')
    ax.text(0, aa_by_class.index('M'), '100', ha='center', va='center', fontsize=6,
            color='white', fontweight='bold')
    for lab, pos, aa in [('P', 1, 'P'), ('G', 1, 'G'), ('D', 2, 'D'), ('E', 2, 'E')]:
        i = aa_by_class.index(aa)
        ax.text(pos, i, f'{m[i, pos]:.0f}', ha='center', va='center', fontsize=6.6,
                color='white' if m[i, pos] > 24 else INK, fontweight='bold')
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label('Substrates carrying this residue (%)', fontsize=8.5)
    cb.ax.tick_params(labelsize=7.5)
    frame(ax, 'A', 'Residue frequency in the 54 substrates',
          'Boxed = positions 2-3. Met at position 1 is 100% (above the colour ceiling).')

    # ---- B: FDR-masked enrichment -----------------------------------------
    ax = fig.add_subplot(gs[1])
    L = lr[aa_by_class].T.values
    Q = qv[aa_by_class].T.values
    shown = np.where(Q < 0.05, L, np.nan)
    v = np.nanmax(np.abs(shown))
    ax.imshow(np.zeros_like(L), cmap=mpl.colors.ListedColormap(['#f2f1ec']), aspect='auto')
    im = ax.imshow(shown, cmap=DIVERGE, norm=TwoSlopeNorm(0, -v, v), aspect='auto')
    for i in range(L.shape[0]):
        for k in range(L.shape[1]):
            if Q[i, k] < 0.05:
                ax.text(k, i, f'{L[i, k]:+.1f}', ha='center', va='center', fontsize=6.2,
                        color='white' if abs(L[i, k]) > v * 0.55 else INK, fontweight='bold')
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label('log$_2$ (substrate / library)', fontsize=8.5)
    cb.ax.tick_params(labelsize=7.5)
    nsig = int((Q < 0.05).sum())
    n_motif = int((Q[:, :3] < 0.05).sum())
    frame(ax, 'B', 'Enrichment over the library background (FDR-masked)',
          f'Beige = not significant. {nsig} of {Q.size} position x residue cells survive q < 0.05: '
          f'{n_motif} at positions 2-3, plus Arg at position 7 (q = 0.04, borderline).')

    fig.suptitle('Figure 5 | Position-by-residue heatmaps of the 54 candidate UBR3 substrates',
                 x=0.075, y=0.962, ha='left', fontsize=14.5, fontweight='bold', color=INK)
    fig.text(0.075, 0.918,
             'The dominant signal is a single block: Pro and Gly at position 2 and Asp at position 3. '
             'The only other cell surviving FDR is Arg at position 7 (22% of substrates vs 7% of the library).',
             ha='left', fontsize=9.2, color=INK2)
    save(fig, 'Figure5_position_residue_heatmaps')


# ================================================================ FIGURE 6
def figure6(lib, hit):
    """Baseline stability: the PSI < 3 / PSI >= 3 split."""
    CUT = U.PSI_CUT
    lo = lib.mean_PSI_control < CUT
    fig = plt.figure(figsize=(13.6, 10.6))
    gs = fig.add_gridspec(2, 2, hspace=0.50, wspace=0.26,
                          left=0.07, right=0.975, top=0.855, bottom=0.065)

    # ---- A: distribution of control PSI ----------------------------------
    ax = fig.add_subplot(gs[0, 0])
    bins = np.linspace(1, 4, 61)
    ax.hist(lib.mean_PSI_control[lo], bins=bins, color=VIOLET, zorder=3,
            label=f'Unstable, PSI < {CUT:.0f}  (n={int(lo.sum()):,})')
    ax.hist(lib.mean_PSI_control[~lo], bins=bins, color=AQUA, zorder=3,
            label=f'Stable, PSI $\\geq$ {CUT:.0f}  (n={int((~lo).sum()):,})')
    ax.axvline(CUT, color=INK, lw=1.4, ls='--', zorder=5)
    top = ax.get_ylim()[1]
    for x, s in zip(hit.mean_PSI_control, [1] * len(hit)):
        ax.plot([x, x], [-top * 0.055, -top * 0.012], color=ORANGE, lw=1.1, zorder=4)
    ax.text(1.05, -top * 0.085, f'the {len(hit)} substrates', fontsize=8.2, color=ORANGE,
            fontweight='bold', va='top')
    ax.set_ylim(-top * 0.14, top * 1.02)
    ax.set_xlabel('Control PSI (mean of both control replicates)')
    ax.set_ylabel('Number of library peptides')
    ax.grid(axis='y', zorder=0)
    ax.legend(loc='upper left')
    title(ax, f'{100 * lo.mean():.0f}% of the library starts below PSI {CUT:.0f}',
          f'Ticks below the axis mark the {len(hit)} substrates - almost all sit in the unstable half')
    panel_tag(ax, 'A')

    # ---- B: control vs KO PSI --------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    ax.scatter(lib.mean_PSI_control, lib.mean_PSI_UBR3KO, s=2.5, color='#e6e4dc',
               linewidths=0, zorder=1, label=f'Library (n={len(lib):,})')
    ax.scatter(hit.mean_PSI_control, hit.mean_PSI_UBR3KO, s=42, color=ORANGE,
               edgecolors=SURFACE, linewidths=0.8, zorder=4,
               label=f'UBR3 substrates (n={len(hit)})')
    ax.plot([1, 4], [1, 4], color=MUTED, lw=1, ls=':', zorder=2)
    ax.axvline(CUT, color=INK, lw=1.2, ls='--', zorder=3)
    ax.axhline(CUT, color=INK, lw=1.2, ls='--', zorder=3)
    ncross = int((hit.crosses_PSI3_up == 'yes').sum())
    ax.add_patch(Rectangle((1, CUT), CUT - 1, 4 - CUT, facecolor=ORANGE, alpha=0.09, zorder=0))
    ax.text(1.08, 3.88, f'crosses upward\n{ncross} of {int((hit.mean_PSI_control < CUT).sum())} substrates',
            fontsize=8.4, color=ORANGE, fontweight='bold', va='top')
    ax.set_xlim(1, 4)
    ax.set_ylim(1, 4)
    ax.set_xlabel('Control PSI')
    ax.set_ylabel('UBR3-KO PSI')
    ax.grid(zorder=0)
    ax.legend(loc='lower right', markerscale=1.3)
    title(ax, 'Over half the substrates cross the PSI 3 line',
          'Dotted diagonal = no change. Shaded quadrant = unstable in control, stable in UBR3 KO.')
    panel_tag(ax, 'B')

    # ---- C: ceiling effect ------------------------------------------------
    ax = fig.add_subplot(gs[1, 0])
    edges = [1, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    cats = pd.cut(lib.mean_PSI_control, edges)
    g = lib.groupby(cats, observed=True)
    means = g.mean_dPSI.mean()
    sems = g.mean_dPSI.sem()
    ns = g.size()
    x = np.arange(len(means))
    cols = [VIOLET if i.right <= CUT else AQUA for i in means.index]
    ax.bar(x, means.values, yerr=sems.values, color=cols, width=0.62, zorder=3,
           edgecolor=SURFACE, linewidth=1.1, error_kw=dict(ecolor=INK2, lw=1.2))
    ax.axhline(0, color=AXIS, lw=1)
    span = means.max() - means.min()
    ax.set_ylim(means.min() - span * 0.34, means.max() + span * 0.20)
    for i, (v, n) in enumerate(zip(means.values, ns.values)):
        off = span * 0.045 if v >= 0 else -span * 0.045
        ax.text(i, v + off, f'{v:+.3f}', ha='center', fontsize=8,
                va='bottom' if v >= 0 else 'top', fontweight='bold')
        # n= labels pinned just above the axis line, independent of the data range
        ax.text(i, 0.015, f'n={n:,}', ha='center', va='bottom', fontsize=7.4, color=MUTED,
                transform=ax.get_xaxis_transform())
    ax.set_xticks(x)
    ax.set_xticklabels([f'{i.left:g}-{i.right:g}' for i in means.index], fontsize=8.4)
    ax.set_xlabel('Control PSI bin')
    ax.set_ylabel('Mean $\\Delta$PSI $\\pm$ SEM')
    ax.grid(axis='y', zorder=0)
    title(ax, 'Stabilisation saturates - and reverses - above PSI 3',
          'Peptides that are already maximally stable have no headroom left to gain')
    panel_tag(ax, 'C')

    # ---- D: substrate rate by stratum ------------------------------------
    ax = fig.add_subplot(gs[1, 1])
    labs = [f'PSI < {CUT:.0f}\n(unstable)', f'PSI $\\geq$ {CUT:.0f}\n(stable)']
    subs = [lib[lo], lib[~lo]]
    rates = [100 * s.is_UBR3_substrate.mean() for s in subs]
    ax.bar([0, 1], rates, color=[VIOLET, AQUA], width=0.5, zorder=3,
           edgecolor=SURFACE, linewidth=1.2)
    tops = []
    for i, s in enumerate(subs):
        k, n = int(s.is_UBR3_substrate.sum()), len(s)
        z, p, nn = 1.96, k / n, n
        c = (p + z * z / (2 * nn)) / (1 + z * z / nn)
        h = z * np.sqrt(p * (1 - p) / nn + z * z / (4 * nn * nn)) / (1 + z * z / nn)
        tops.append(100 * (c + h))
        ax.plot([i, i], [100 * (c - h), 100 * (c + h)], color=INK2, lw=1.3, zorder=4)
        ax.text(i + 0.31, rates[i], f'{rates[i]:.3f}%', ha='left', va='center',
                fontweight='bold', fontsize=10.5)
        ax.text(i, 0.015, f'{k} / {n:,}', ha='center', va='bottom', fontsize=8.4, color=INK2,
                transform=ax.get_xaxis_transform())
    ax.set_xticks([0, 1])
    ax.set_xticklabels(labs, fontsize=9.5)
    ax.set_xlim(-0.6, 1.9)
    top = max(tops)
    ax.set_ylim(-top * 0.14, top * 1.30)
    ax.set_ylabel('UBR3-stabilised peptides (% of stratum)')
    ax.grid(axis='y', zorder=0)
    orv, pv = stats.fisher_exact([[int(lib[lo].is_UBR3_substrate.sum()),
                                   int((~lib[lo].is_UBR3_substrate).sum())],
                                  [int(lib[~lo].is_UBR3_substrate.sum()),
                                   int((~lib[~lo].is_UBR3_substrate).sum())]])
    yb = top * 1.12
    ax.plot([0, 0, 1, 1], [yb - top * 0.025, yb, yb, yb - top * 0.025], color=INK2, lw=1.1, zorder=5)
    ax.text(0.5, yb + top * 0.015, f'Fisher exact  OR = {orv:.1f},  p = {pv:.1e}',
            ha='center', fontsize=8.4, color=INK)
    title(ax, 'Substrates are 16x more common below the line',
          'Error bars are Wilson 95% confidence intervals; counts below each bar')
    panel_tag(ax, 'D')

    fig.suptitle('Figure 6 | Baseline stability: peptides above and below PSI 3',
                 x=0.07, y=0.968, ha='left', fontsize=14.5, fontweight='bold', color=INK)
    fig.text(0.07, 0.938,
             'PSI is the read-weighted mean FACS bin (~1-4); PSI = 3 splits the library into unstable and '
             'already-stable halves. A degron substrate should start unstable and gain stability when its '
             'E3 ligase is lost - which is what the 54 substrates do.',
             ha='left', fontsize=9.2, color=INK2)
    save(fig, 'Figure6_PSI_baseline_stability')


# ================================================================ FIGURE 7
def figure7(lib, hit):
    """The motif effect is not a by-product of baseline instability."""
    CUT = U.PSI_CUT
    lo = lib.mean_PSI_control < CUT
    groups = ['[P/G]-[E/D]', '[P/G]-other', 'non-P/G']
    strata = [(f'PSI < {CUT:.0f} (unstable)', lo, VIOLET),
              (f'PSI $\\geq$ {CUT:.0f} (stable)', ~lo, AQUA)]

    fig = plt.figure(figsize=(13.6, 10.6))
    gs = fig.add_gridspec(2, 2, hspace=0.50, wspace=0.28,
                          left=0.07, right=0.975, top=0.855, bottom=0.065)

    # ---- A: stratified substrate rate ------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    w = 0.36
    x = np.arange(len(groups))
    for j, (slab, m, col) in enumerate(strata):
        sub = lib[m]
        rates, ns, ks = [], [], []
        for gname in groups:
            g = sub[sub.motif_class == gname]
            rates.append(100 * g.is_UBR3_substrate.mean())
            ns.append(len(g))
            ks.append(int(g.is_UBR3_substrate.sum()))
        pos = x + (j - 0.5) * w
        ax.bar(pos, rates, width=w * 0.92, color=col, zorder=3,
               edgecolor=SURFACE, linewidth=1.1, label=slab)
        for xi, r, k, n in zip(pos, rates, ks, ns):
            ax.text(xi, r + 0.28, f'{r:.2f}%', ha='center', fontsize=7.8, fontweight='bold')
            ax.text(xi, 0.015, f'{k}/{n:,}', ha='center', va='bottom', fontsize=7,
                    color=MUTED, transform=ax.get_xaxis_transform())
    ax.set_xticks(x)
    ax.set_xticklabels(groups, fontsize=9.4)
    ax.set_ylabel('UBR3-stabilised peptides (% of that group)')
    ax.set_ylim(-1.9, 15.6)
    ax.grid(axis='y', zorder=0)
    ax.legend(loc='upper right', title='Baseline stratum', title_fontsize=8.5)
    title(ax, 'The motif effect survives stratification',
          'Within unstable peptides alone: 13.1% vs 0.33%, Fisher OR 45.9, p = 3 $\\times$ 10$^{-17}$')
    panel_tag(ax, 'A')

    # ---- B: motif frequency per stratum (the control) ---------------------
    ax = fig.add_subplot(gs[0, 1])
    labs = ['[P/G]-[E/D]', 'P or G at pos 2']
    vals = [[100 * lib[m].is_PG_ED.mean() for _, m, _ in strata],
            [100 * lib[m].is_PG.mean() for _, m, _ in strata]]
    x2 = np.arange(2)
    for j, (slab, m, col) in enumerate(strata):
        ax.bar(x2 + (j - 0.5) * w, [vals[0][j], vals[1][j]], width=w * 0.92, color=col,
               zorder=3, edgecolor=SURFACE, linewidth=1.1, label=slab)
        for xi, v in zip(x2 + (j - 0.5) * w, [vals[0][j], vals[1][j]]):
            ax.text(xi, v + 0.22, f'{v:.2f}%', ha='center', fontsize=8, fontweight='bold')
    ax.set_xticks(x2)
    ax.set_xticklabels(labs, fontsize=9.4)
    ax.set_ylabel('Peptides carrying the feature (% of stratum)')
    ax.set_ylim(0, 19)
    ax.grid(axis='y', zorder=0)
    ax.legend(loc='upper left', title='Baseline stratum', title_fontsize=8.5)
    title(ax, '[P/G]-[E/D] is not over-represented among unstable peptides',
          'Bare Pro/Gly is skewed (16.0% vs 7.5%), so the acidic residue carries the UBR3 dependence')
    panel_tag(ax, 'B')

    # ---- C: crossing rate --------------------------------------------------
    ax = fig.add_subplot(gs[1, 0])
    below_lib = lib[lo]
    sets = [('UBR3 substrates', hit[hit.mean_PSI_control < CUT], ORANGE),
            ('[P/G]-[E/D], not stabilised',
             below_lib[(below_lib.motif_class == '[P/G]-[E/D]') & ~below_lib.is_UBR3_substrate], YELLOW),
            ('All library peptides', below_lib, MUTED)]
    ys = np.arange(len(sets))[::-1]
    for yv, (lab, s, col) in zip(ys, sets):
        pct = 100 * (s.crosses_PSI3_up == 'yes').mean()
        ax.barh(yv, pct, color=col, height=0.5, zorder=3, edgecolor=SURFACE, linewidth=1.2)
        ax.text(pct + 1.2, yv, f'{pct:.1f}%   ({int((s.crosses_PSI3_up == "yes").sum())} of {len(s):,})',
                va='center', fontsize=9, fontweight='bold', color=INK)
        ax.text(0.6, yv + 0.36, lab, va='bottom', fontsize=8.8, color=INK2)
    ax.set_yticks([])
    ax.set_xlim(0, 82)
    ax.set_xlabel(f'Peptides rising above PSI {CUT:.0f} on UBR3 loss (%)')
    ax.grid(axis='x', zorder=0)
    ax.spines['left'].set_visible(False)
    title(ax, 'Crossing is graded, not all-or-nothing',
          f'Peptides starting below PSI {CUT:.0f}; non-substrates with the motif still cross at '
          'twice background')
    panel_tag(ax, 'C')

    # ---- D: dPSI distribution within the unstable stratum -----------------
    ax = fig.add_subplot(gs[1, 1])
    data, labels, cols = [], [], []
    for gname, col in zip(groups, [ORANGE, AQUA, MUTED]):
        v = lib[lo & (lib.motif_class == gname)].mean_dPSI.values
        data.append(v)
        labels.append(f'{gname}\n(n={len(v):,})')
        cols.append(col)
    bp = ax.boxplot(data, patch_artist=True, widths=0.5, showfliers=False,
                    medianprops=dict(color=INK, lw=1.6),
                    whiskerprops=dict(color=INK2, lw=1),
                    capprops=dict(color=INK2, lw=1),
                    boxprops=dict(edgecolor=SURFACE, lw=1.2))
    for patch, c in zip(bp['boxes'], cols):
        patch.set_facecolor(c)
    for i, v in enumerate(data, 1):
        jitter = (np.random.default_rng(0).random(len(v)) - 0.5) * 0.26
        if len(v) < 400:
            ax.scatter(i + jitter, v, s=5, color=INK2, alpha=0.35, linewidths=0, zorder=4)
    ax.axhline(0, color=AXIS, lw=1)
    ax.set_xticks(range(1, 4))
    ax.set_xticklabels(labels, fontsize=8.8)
    ax.set_ylabel('Mean $\\Delta$PSI')
    ax.grid(axis='y', zorder=0)
    p = stats.mannwhitneyu(data[0], data[2])[1]
    ax.text(0.5, 0.96, f'[P/G]-[E/D] vs non-P/G\nMann-Whitney p = {p:.1e}',
            transform=ax.transAxes, ha='center', va='top', fontsize=8.4, color=INK)
    title(ax, 'Within unstable peptides, the motif still shifts $\\Delta$PSI',
          f'Only peptides starting below PSI {CUT:.0f}; boxes are median and quartiles')
    panel_tag(ax, 'D')

    fig.suptitle('Figure 7 | The [P/G]-[E/D] effect is not a by-product of baseline instability',
                 x=0.07, y=0.968, ha='left', fontsize=14.5, fontweight='bold', color=INK)
    fig.text(0.07, 0.938,
             'Substrates start unstable, so the motif could in principle be enriched simply because '
             'Pro/Gly N-termini are unstable. Stratifying by control PSI shows it is not: the motif is '
             'evenly distributed across strata yet enriches for substrates within each one.',
             ha='left', fontsize=9.2, color=INK2)
    save(fig, 'Figure7_motif_vs_baseline_stability')


# ================================================================ FIGURE 8
def figure8(lib, hit):
    """Why control PSI (not dPSI) defines the strata, and why the cut is not arbitrary."""
    import statsmodels.api as sm
    from scipy.signal import argrelextrema
    from scipy.stats import gaussian_kde
    CUT = U.PSI_CUT

    fig = plt.figure(figsize=(13.6, 10.6))
    gs = fig.add_gridspec(2, 2, hspace=0.50, wspace=0.28,
                          left=0.07, right=0.975, top=0.855, bottom=0.065)

    # ---- A: the two axes are nearly independent --------------------------
    ax = fig.add_subplot(gs[0, 0])
    ax.scatter(lib.mean_PSI_control, lib.mean_dPSI, s=2.5, color='#e6e4dc',
               linewidths=0, zorder=1, label=f'Library (n={len(lib):,})')
    ax.scatter(hit.mean_PSI_control, hit.mean_dPSI, s=40, color=ORANGE,
               edgecolors=SURFACE, linewidths=0.8, zorder=4,
               label=f'UBR3 substrates (n={len(hit)})')
    ax.axvline(CUT, color=INK, lw=1.2, ls='--', zorder=3)
    ax.axhline(0, color=AXIS, lw=1)
    rho = stats.spearmanr(lib.mean_PSI_control, lib.mean_dPSI)
    ax.text(0.03, 0.97, f'Spearman $\\rho$ = {rho.statistic:.2f}\nshared variance '
                        f'{100 * stats.pearsonr(lib.mean_PSI_control, lib.mean_dPSI)[0] ** 2:.1f}%',
            transform=ax.transAxes, va='top', fontsize=8.6, color=INK)
    ax.set_xlabel('Control PSI  (baseline stability)')
    ax.set_ylabel('Mean $\\Delta$PSI  (UBR3 dependence)')
    ax.set_xlim(1, 4)
    ax.set_ylim(-1.0, 1.8)
    ax.grid(zorder=0)
    ax.legend(loc='lower right', markerscale=1.3)
    title(ax, 'PSI and $\\Delta$PSI measure different things',
          'Baseline stability explains ~3% of the variance in UBR3 dependence - they are near-orthogonal')
    panel_tag(ax, 'A')

    # ---- B: control PSI density and its antimodes -------------------------
    ax = fig.add_subplot(gs[0, 1])
    gx = np.linspace(1, 4, 900)
    # bw 0.25: the modal structure is stable here. Narrower bandwidths add a
    # shoulder near 2.9 that disappears again by 0.20, so it is not a real mode.
    gy = gaussian_kde(lib.mean_PSI_control, bw_method=0.25)(gx)
    ax.fill_between(gx, gy, color=BLUE, alpha=0.16, zorder=2)
    ax.plot(gx, gy, color=BLUE, lw=2, zorder=3)
    vall, peaks = gx[argrelextrema(gy, np.less)[0]], gx[argrelextrema(gy, np.greater)[0]]
    for v in peaks:
        h = gy[np.argmin(abs(gx - v))]
        ax.plot([v, v], [0, h], color=MUTED, lw=0.9, ls=':', zorder=4)
        ax.text(v, h + 0.008, f'{v:.2f}', ha='center', fontsize=7.6, color=MUTED)
    for v in vall:
        ax.scatter([v], [gy[np.argmin(abs(gx - v))]], s=46, color=RED, zorder=6,
                   edgecolors=SURFACE, linewidths=0.8)
    anti = vall[-1]
    ax.axvline(anti, color=RED, lw=1.3, ls=':', zorder=5)
    ax.axvline(CUT, color=INK, lw=1.4, ls='--', zorder=5)
    ax.annotate(f'nearest real antimode\nis {anti:.2f}, not {CUT:.0f}',
                xy=(anti, gy[np.argmin(abs(gx - anti))]),
                xytext=(1.06, max(gy) * 0.72), fontsize=8.6, color=INK, fontweight='bold',
                arrowprops=dict(arrowstyle='->', color=INK2, lw=1.1))
    ax.text(CUT + 0.04, max(gy) * 0.98, f'cut used\nin Fig 6-7', fontsize=8.2, color=INK,
            va='top', fontweight='bold')
    ax.set_xlabel('Control PSI')
    ax.set_ylabel('Kernel density')
    ax.set_xlim(1, 4)
    ax.set_ylim(0, max(gy) * 1.16)
    ax.grid(zorder=0)
    ax.legend(handles=[Line2D([], [], marker='o', ls='none', color=RED, markersize=7,
                              label='density minimum (antimode)')], loc='upper left')
    title(ax, 'PSI 3 is a convention, not a natural boundary',
          'Three robust modes. The data-driven split sits at 2.61; panel C shows the choice '
          'changes nothing.')
    panel_tag(ax, 'B')

    # ---- C: sensitivity of the conclusion to the cutoff -------------------
    ax = fig.add_subplot(gs[1, 0])
    cuts = [2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4]
    ors, ps = [], []
    for c in cuts:
        s = lib[lib.mean_PSI_control < c]
        a, b = s[s.motif_class == '[P/G]-[E/D]'], s[s.motif_class == 'non-P/G']
        o, p = stats.fisher_exact([[int(a.is_UBR3_substrate.sum()), int((~a.is_UBR3_substrate).sum())],
                                   [int(b.is_UBR3_substrate.sum()), int((~b.is_UBR3_substrate).sum())]])
        ors.append(o)
        ps.append(p)
    cols = [ORANGE if c == CUT else BLUE for c in cuts]
    ax.plot(cuts, ors, color=BLUE, lw=1.8, zorder=3)
    ax.scatter(cuts, ors, s=[95 if c == CUT else 52 for c in cuts], color=cols,
               edgecolors=SURFACE, linewidths=1.1, zorder=4)
    for c, o, p in zip(cuts, ors, ps):
        ax.text(c, o * 1.10, f'{o:.0f}', ha='center', fontsize=7.8, fontweight='bold',
                color=ORANGE if c == CUT else INK)
    ax.annotate('the cut used\nin Figures 6-7', xy=(CUT, ors[cuts.index(CUT)]),
                xytext=(3.05, 26), fontsize=8.4, color=ORANGE, fontweight='bold',
                arrowprops=dict(arrowstyle='->', color=ORANGE, lw=1.1))
    ax.set_yscale('log')
    ax.set_ylim(18, 190)
    ax.set_xlabel('PSI cutoff defining "unstable at baseline"')
    ax.set_ylabel('Odds ratio, [P/G]-[E/D] vs non-P/G\n(within the unstable stratum)')
    ax.grid(zorder=0, which='both')
    ax.text(0.97, 0.06, f'every cutoff: p $\\leq$ {max(ps):.0e}', transform=ax.transAxes,
            ha='right', fontsize=8.4, color=INK)
    title(ax, 'The conclusion does not depend on where the line goes',
          'Odds ratio stays between 42 and 88 for every cutoff from 2.0 to 3.4')
    panel_tag(ax, 'C')

    # ---- D: cutoff-free logistic regression -------------------------------
    ax = fig.add_subplot(gs[1, 1])
    d = lib.copy()
    d['y'] = d.is_UBR3_substrate.astype(int)
    d['PGED'] = (d.motif_class == '[P/G]-[E/D]').astype(int)
    d['PGother'] = (d.motif_class == '[P/G]-other').astype(int)
    fit = sm.Logit(d.y, sm.add_constant(d[['PGED', 'PGother', 'mean_PSI_control']])).fit(disp=0)
    ci = fit.conf_int()
    terms = ['PGED', 'PGother', 'mean_PSI_control']
    labs = ['[P/G]-[E/D]\nvs non-P/G', '[P/G]-other\nvs non-P/G',
            'control PSI\n(per +1 unit)']
    y = np.arange(len(terms))[::-1]
    for yv, t, c in zip(y, terms, [ORANGE, AQUA, VIOLET]):
        o, lo_, hi_ = np.exp(fit.params[t]), np.exp(ci.loc[t, 0]), np.exp(ci.loc[t, 1])
        ax.plot([lo_, hi_], [yv, yv], color=c, lw=2.4, zorder=3, solid_capstyle='round')
        ax.scatter([o], [yv], s=95, color=c, edgecolors=SURFACE, linewidths=1.2, zorder=4)
        ax.text(hi_ * 1.25, yv, f'OR {o:.2f}   p = {fit.pvalues[t]:.1e}', va='center',
                fontsize=8.6, fontweight='bold', color=INK)
    ax.axvline(1, color=AXIS, lw=1.2, ls='--', zorder=2)
    ax.set_xscale('log')
    ax.set_xlim(0.28, 700)
    ax.set_yticks(y)
    ax.set_yticklabels(labs, fontsize=9)
    ax.set_xlabel('Adjusted odds ratio for being a UBR3 substrate (log scale)')
    ax.grid(axis='x', zorder=0, which='both')
    ax.spines['left'].set_visible(False)
    ax.set_ylim(-0.7, len(terms) - 0.3)
    title(ax, 'With no cutoff at all, the motif effect stands',
          'Logistic regression with control PSI as a continuous covariate; bars are 95% CIs')
    panel_tag(ax, 'D')

    fig.suptitle('Figure 8 | Why baseline stability is measured by control PSI, not by $\\Delta$PSI',
                 x=0.07, y=0.968, ha='left', fontsize=14.5, fontweight='bold', color=INK)
    fig.text(0.07, 0.938,
             'Control PSI is how stable a peptide is before UBR3 is removed; $\\Delta$PSI is how much it '
             'changes when UBR3 goes. The 54 substrates were selected on $\\Delta$PSI, so splitting on it '
             'would be circular - every $\\Delta$PSI cut contains all 54 by construction.',
             ha='left', fontsize=9.2, color=INK2)
    save(fig, 'Figure8_PSI_vs_dPSI_and_cutoff_robustness')


def main():
    os.makedirs(U.FIGS, exist_ok=True)
    lib, hit = U.load()
    enr = pd.read_excel(os.path.join(U.HERE, 'UBR3_PG_substrate_tables.xlsx'),
                        '06_enrich_pos2_residue')
    print('rendering figures ...')
    print('  running 480 position x residue Fisher tests ...')
    sig = U.enrichment_test(list(hit.peptide_24mer), list(lib.peptide_24mer))
    figure1(lib, hit, enr)
    figure2(lib, hit)
    figure3(lib, hit, sig)
    figure4(lib, hit)
    figure5(lib, hit, sig)
    figure6(lib, hit)
    figure7(lib, hit)
    figure8(lib, hit)
    print('done ->', U.FIGS)


if __name__ == '__main__':
    main()
