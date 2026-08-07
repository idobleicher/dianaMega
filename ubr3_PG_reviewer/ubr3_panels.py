#!/usr/bin/env python
"""Every figure panel as a standalone function that draws into a given Axes.

`make_figures.py` composes these into the multi-panel manuscript figures and
`make_panels.py` renders each one on its own for the slide deck, so the two
outputs are guaranteed to show the same thing.

Each panel is `p<fig><letter>(ax, C)` where C is the context from `context()`.
"""
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
SURFACE, INK, INK2, MUTED, GRID, AXIS = (
    '#fcfcfb', '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7')
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN = (
    '#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300')
RED, VIOLET = '#e34948', '#4a3aa7'

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

AA_BY_CLASS = [a for c in U.CLASS_ORDER for a in U.CLASS_MEMBERS[c]]


def stars(q):
    return '***' if q < 1e-3 else '**' if q < 1e-2 else '*' if q < 0.05 else 'ns'


def panel_tag(ax, letter, C, xoff=-46, yoff=30):
    """Panel letter - suppressed when the panel is rendered on its own slide."""
    if not C.get('tags', True):
        return
    ax.annotate(letter, xy=(0, 1), xycoords='axes fraction',
                xytext=(xoff, yoff), textcoords='offset points',
                fontsize=14, fontweight='bold', va='bottom', ha='left', color=INK)


def title(ax, text, sub=None):
    ax.set_title(text, loc='left', fontweight='bold', color=INK, pad=26 if sub else 10)
    if sub:
        ax.annotate(sub, xy=(0, 1), xycoords='axes fraction', xytext=(0, 7),
                    textcoords='offset points', fontsize=8.3, color=INK2,
                    va='bottom', ha='left')


def context(lib, hit, sig, tags=True, standalone=False):
    """Everything the panels need, computed once."""
    hs, ls = list(hit.peptide_24mer), list(lib.peptide_24mer)
    CUT = U.PSI_CUT
    lo = lib.mean_PSI_control < CUT
    C = dict(lib=lib, hit=hit, sig=sig, tags=tags, standalone=standalone,
             hs=hs, ls=ls, CUT=CUT, lo=lo,
             fh=U.freq_matrix(hs), fl=U.freq_matrix(ls),
             ch=U.class_matrix(hs), cl=U.class_matrix(ls),
             pged=lib[lib.motif_class == '[P/G]-[E/D]'].sort_values('mean_dPSI', ascending=False),
             groups=['[P/G]-[E/D]', '[P/G]-other', 'non-P/G'],
             strata=[(f'PSI < {CUT:g} (unstable)', lo, VIOLET),
                     (f'PSI $\\geq$ {CUT:g} (stable)', ~lo, AQUA)])
    C['clr'], C['clq'], _ = U.enrichment_test(hs, ls, groups=U.CLASS_MEMBERS)
    return C


# ============================================================ FIGURE 2
def p2a(ax, C):
    pged = C['pged']
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
    panel_tag(ax, 'A', C)


def p2b(ax, C):
    lib = C['lib']
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
    panel_tag(ax, 'B', C)


def p2c(ax, C):
    lib, groups = C['lib'], C['groups']
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
                va='center', ha='center' if pk > 14 else 'left', fontsize=9.4,
                fontweight='bold', color='white' if pk > 14 else ORANGE, zorder=5)
    ax.set_yticks(y)
    ax.set_yticklabels(groups, fontsize=10)
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.9, 2.6)
    ax.set_xlabel('Percent of peptides in the group')
    ax.grid(axis='x', zorder=0)
    ax.spines['left'].set_visible(False)
    ax.legend(handles=[Patch(color=ORANGE, label='UBR3-stabilised'),
                       Patch(color='#dedcd4', label='Not stabilised')],
              loc='lower center', ncol=2, columnspacing=1.6)
    title(ax, 'Motif enrichment is 17-fold, but the motif is not sufficient',
          'Even in the best group, 91% of motif-bearing peptides are unaffected by UBR3 loss')
    panel_tag(ax, 'C', C)


def p2d(ax, C):
    lib, pged = C['lib'], C['pged']
    bg = lib[~lib.is_PG_ED]
    ax.scatter(bg[U.D1], bg[U.D2], s=3, color='#e6e4dc', linewidths=0, zorder=1,
               label=f'Other library peptides (n={len(bg):,})')
    nn = pged[~pged.is_UBR3_substrate]
    ss = pged[pged.is_UBR3_substrate]
    ax.scatter(nn[U.D1], nn[U.D2], s=26, color=MUTED, linewidths=0.6, edgecolors=SURFACE,
               zorder=3, label=f'[P/G]-[E/D], not stabilised (n={len(nn)})')
    ax.scatter(ss[U.D1], ss[U.D2], s=52, color=ORANGE, linewidths=0.9, edgecolors=SURFACE,
               zorder=4, label=f'[P/G]-[E/D], UBR3 substrate (n={len(ss)})')
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
    panel_tag(ax, 'D', C)


# ============================================================ FIGURE 3
POS3 = list(range(2, 25))


def _trim(m):
    m = m.loc[POS3].copy()
    m.index = range(len(POS3))
    return m


def _dress3(ax, C, ylab, sub, tag, ttl):
    ax.set_xticks(range(len(POS3)))
    ax.set_xticklabels(POS3, fontsize=8.5)
    ax.set_xlabel('Position in the encoded 24-mer   '
                  '(position 1 is the invariant initiator Met and is omitted)')
    ax.set_ylabel(ylab)
    ax.grid(axis='y', zorder=0, alpha=0.55)
    ax.set_axisbelow(True)
    ax.axvspan(-0.5, 1.5, color=YELLOW, alpha=0.16, zorder=0)
    title(ax, ttl, sub)
    panel_tag(ax, tag, C, xoff=-52)
    if C.get('standalone'):
        ax.legend(handles=[Patch(color=CLASS_COLOR[c], label=c) for c in U.CLASS_ORDER],
                  loc='upper right', ncol=3, fontsize=8, columnspacing=1.1)


def _info3(C):
    info = _trim(logomaker.transform_matrix(U.position_matrix(C['hs']),
                                            from_type='counts', to_type='information'))
    linfo = _trim(logomaker.transform_matrix(U.position_matrix(C['ls']),
                                             from_type='counts', to_type='information'))
    return info, linfo, max(info.sum(1).max(), linfo.sum(1).max()) * 1.18


def p3a(ax, C):
    info, _, top = _info3(C)
    logomaker.Logo(info, ax=ax, color_scheme=AA_COLOR, show_spines=False,
                   stack_order='big_on_top')
    ax.set_ylim(0, top)
    _dress3(ax, C, 'Information content (bits)',
            f'Letter height = information content; the shaded block is positions 2-3 '
            f'(n={len(C["hs"])} substrates)',
            'A', 'Sequence logo of the 54 candidate UBR3 substrates')
    ax.annotate('Pro/Gly at position 2,\nAsp/Glu at position 3',
                xy=(0.9, info.sum(1).max() * 0.99), xytext=(3.2, top * 0.80),
                fontsize=9, fontweight='bold', color=INK,
                arrowprops=dict(arrowstyle='->', color=INK2, lw=1.1))


def p3b(ax, C):
    _, linfo, top = _info3(C)
    logomaker.Logo(linfo, ax=ax, color_scheme=AA_COLOR, show_spines=False,
                   stack_order='big_on_top')
    ax.set_ylim(0, top)
    _dress3(ax, C, 'Information content (bits)',
            f'Background from all {len(C["ls"]):,} library peptides, on the same y-scale as '
            'panel A - essentially flat',
            'B', 'Sequence logo of the full N24mer library (background)')


def p3c(ax, C):
    _, _, signed = C['sig']
    praw = _trim(signed)
    praw = praw.where(praw.abs() >= -np.log10(0.05), 0.0)
    logomaker.Logo(praw, ax=ax, color_scheme=AA_COLOR, show_spines=False,
                   stack_order='big_on_top', flip_below=True)
    ax.axhline(0, color=AXIS, lw=1)
    for lv, lab, st in [(-np.log10(0.05), 'p = 0.05', ':'),
                        (-np.log10(0.05 / (23 * 20)), 'Bonferroni', '--')]:
        for s in (1, -1):
            ax.axhline(s * lv, color=RED, lw=0.9, ls=st, zorder=5)
        ax.text(len(POS3) - 0.4, lv + 0.12, lab, fontsize=7.6, color=RED, ha='right')
    _dress3(ax, C, 'signed $-$log$_{10}$ $p$  (Fisher exact)',
            'Only residues with p < 0.05 are drawn. Above the line = enriched in substrates, '
            'below = depleted.',
            'C', 'Statistical enrichment logo: what actually distinguishes substrates from the library')


# ============================================================ FIGURE 4
def _stack4(ax, C, mat, lab, tag, sub, legend):
    bottom = np.zeros(24)
    xs = np.arange(1, 25)
    for c in U.CLASS_ORDER:
        v = 100 * mat[c].values
        ax.bar(xs, v, bottom=bottom, color=CLASS_COLOR[c], width=0.76,
               label=f'{c} ({",".join(U.CLASS_MEMBERS[c])})', zorder=3,
               edgecolor=SURFACE, linewidth=1.1)
        bottom += v
    ax.set_xticks(xs)
    ax.set_xticklabels(xs, fontsize=8.3)
    ax.set_ylim(0, 100)
    ax.set_xlim(0.4, 24.6)
    ax.set_ylabel('Composition (%)')
    ax.set_xlabel('Position in the 24-mer')
    ax.set_axisbelow(True)
    if legend:
        if C.get('standalone'):
            ax.legend(ncol=3, loc='lower center', bbox_to_anchor=(0.5, 1.02),
                      columnspacing=1.3, handlelength=1.4, fontsize=8)
        else:
            ax.legend(ncol=6, loc='lower center', bbox_to_anchor=(0.5, 1.20),
                      columnspacing=1.4, handlelength=1.4, fontsize=8.2)
    title(ax, f'Amino-acid class composition at each position - {lab}', sub)
    panel_tag(ax, tag, C, xoff=-52)


def p4a(ax, C):
    _stack4(ax, C, C['ch'], '54 UBR3 substrates', 'A',
            'Position 2 of the substrates is 48% Special (Pro/Gly); position 3 is 37% Acidic',
            legend=True)


def p4b(ax, C):
    _stack4(ax, C, C['cl'], f'{len(C["ls"]):,} library peptides', 'B',
            'The background is flat from position 4 onward', legend=C.get('standalone', False))


def p4c(ax, C):
    m = C['clr'][U.CLASS_ORDER].T.values
    q = C['clq'][U.CLASS_ORDER].T.values
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
    cb = ax.get_figure().colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label('log$_2$ (substrate / library)', fontsize=8.5)
    cb.ax.tick_params(labelsize=7.5)
    ax.set_xticks(np.arange(-0.5, 24, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, 6, 1), minor=True)
    ax.grid(which='minor', color=SURFACE, lw=1.1)
    ax.tick_params(which='minor', length=0)
    title(ax, 'Class enrichment by position (FDR-masked)',
          'Beige = not significant at q < 0.05. Positions 2 and 3 dominate; Basic at position 7 '
          'is the only other survivor.')
    panel_tag(ax, 'C', C, xoff=-70)


def p4d(ax, C):
    hs, ls = C['hs'], C['ls']
    ph = np.array([sum(U.AA_CLASS[a] == c for p in hs for a in p) / (len(hs) * 24)
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
                       Patch(facecolor=MUTED, label='Library')], loc='upper left')
    title(ax, 'Whole-peptide composition is unchanged',
          'Averaged over all 24 positions the substrates look like the library')
    panel_tag(ax, 'D', C, xoff=-62)


# ============================================================ FIGURE 5
def _frame5(ax, C, tag, ttl, sub):
    ax.set_yticks(range(len(AA_BY_CLASS)))
    ax.set_yticklabels(AA_BY_CLASS, fontsize=8.5, fontweight='bold')
    for t, a in zip(ax.get_yticklabels(), AA_BY_CLASS):
        t.set_color(CLASS_COLOR[U.AA_CLASS[a]])
    ax.set_xticks(range(24))
    ax.set_xticklabels(range(1, 25), fontsize=7.6)
    ax.set_xlabel('Position in the 24-mer')
    ax.set_ylabel('Amino acid (grouped by chemical class)')
    acc = 0
    for c in U.CLASS_ORDER[:-1]:
        acc += len(U.CLASS_MEMBERS[c])
        ax.axhline(acc - 0.5, color=INK, lw=1.4)
    ax.add_patch(Rectangle((0.5, -0.5), 2, len(AA_BY_CLASS), fill=False,
                           edgecolor=INK, lw=1.8, zorder=6))
    title(ax, ttl, sub)
    panel_tag(ax, tag, C, xoff=-64, yoff=34)


def p5a(ax, C):
    m = (100 * C['fh'])[AA_BY_CLASS].T.values
    im = ax.imshow(m, cmap=SEQ_BLUE, vmin=0, vmax=28, aspect='auto')
    ax.text(0, AA_BY_CLASS.index('M'), '100', ha='center', va='center', fontsize=6,
            color='white', fontweight='bold')
    for pos, aa in [(1, 'P'), (1, 'G'), (2, 'D'), (2, 'E')]:
        i = AA_BY_CLASS.index(aa)
        ax.text(pos, i, f'{m[i, pos]:.0f}', ha='center', va='center', fontsize=6.6,
                color='white' if m[i, pos] > 24 else INK, fontweight='bold')
    cb = ax.get_figure().colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label('Substrates carrying this residue (%)', fontsize=8.5)
    cb.ax.tick_params(labelsize=7.5)
    _frame5(ax, C, 'A', 'Residue frequency in the 54 substrates',
            'Boxed = positions 2-3. Met at position 1 is 100% (above the colour ceiling).')


def p5b(ax, C):
    lr, qv, _ = C['sig']
    L = lr[AA_BY_CLASS].T.values
    Q = qv[AA_BY_CLASS].T.values
    shown = np.where(Q < 0.05, L, np.nan)
    v = np.nanmax(np.abs(shown))
    ax.imshow(np.zeros_like(L), cmap=mpl.colors.ListedColormap(['#f2f1ec']), aspect='auto')
    im = ax.imshow(shown, cmap=DIVERGE, norm=TwoSlopeNorm(0, -v, v), aspect='auto')
    for i in range(L.shape[0]):
        for k in range(L.shape[1]):
            if Q[i, k] < 0.05:
                ax.text(k, i, f'{L[i, k]:+.1f}', ha='center', va='center', fontsize=6.2,
                        color='white' if abs(L[i, k]) > v * 0.55 else INK, fontweight='bold')
    cb = ax.get_figure().colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label('log$_2$ (substrate / library)', fontsize=8.5)
    cb.ax.tick_params(labelsize=7.5)
    nsig, nmot = int((Q < 0.05).sum()), int((Q[:, :3] < 0.05).sum())
    _frame5(ax, C, 'B', 'Enrichment over the library background (FDR-masked)',
            f'Beige = not significant. {nsig} of {Q.size} position x residue cells survive '
            f'q < 0.05: {nmot} at positions 2-3, plus Arg at position 7 (q = 0.04, borderline).')


# ============================================================ FIGURE 6
def p6a(ax, C):
    lib, hit, CUT, lo = C['lib'], C['hit'], C['CUT'], C['lo']
    bins = np.linspace(1, 4, 61)
    ax.hist(lib.mean_PSI_control[lo], bins=bins, color=VIOLET, zorder=3,
            label=f'Unstable, PSI < {CUT:g}  (n={int(lo.sum()):,})')
    ax.hist(lib.mean_PSI_control[~lo], bins=bins, color=AQUA, zorder=3,
            label=f'Stable, PSI $\\geq$ {CUT:g}  (n={int((~lo).sum()):,})')
    ax.axvline(CUT, color=INK, lw=1.4, ls='--', zorder=5)
    top = ax.get_ylim()[1]
    for x in hit.mean_PSI_control:
        ax.plot([x, x], [-top * 0.055, -top * 0.012], color=ORANGE, lw=1.1, zorder=4)
    ax.text(1.05, -top * 0.085, f'the {len(hit)} substrates', fontsize=8.2, color=ORANGE,
            fontweight='bold', va='top')
    ax.set_ylim(-top * 0.14, top * 1.02)
    ax.set_xlabel('Control PSI (mean of both control replicates)')
    ax.set_ylabel('Number of library peptides')
    ax.grid(axis='y', zorder=0)
    ax.legend(loc='upper left')
    pct_hit_lo = 100 * (hit.mean_PSI_control < CUT).mean()
    title(ax, f'{100 * lo.mean():.0f}% of the library starts below PSI {CUT:g}',
          f'Ticks below the axis mark the {len(hit)} substrates; {pct_hit_lo:.0f}% of them sit in '
          'the unstable half')
    panel_tag(ax, 'A', C)


def p6b(ax, C):
    lib, hit, CUT = C['lib'], C['hit'], C['CUT']
    ax.scatter(lib.mean_PSI_control, lib.mean_PSI_UBR3KO, s=2.5, color='#e6e4dc',
               linewidths=0, zorder=1, label=f'Library (n={len(lib):,})')
    ax.scatter(hit.mean_PSI_control, hit.mean_PSI_UBR3KO, s=42, color=ORANGE,
               edgecolors=SURFACE, linewidths=0.8, zorder=4,
               label=f'UBR3 substrates (n={len(hit)})')
    ax.plot([1, 4], [1, 4], color=MUTED, lw=1, ls=':', zorder=2)
    ax.axvline(CUT, color=INK, lw=1.2, ls='--', zorder=3)
    ax.axhline(CUT, color=INK, lw=1.2, ls='--', zorder=3)
    ncross = int((hit.crosses_PSI3_up == 'yes').sum())
    nbelow = int((hit.mean_PSI_control < CUT).sum())
    ax.add_patch(Rectangle((1, CUT), CUT - 1, 4 - CUT, facecolor=ORANGE, alpha=0.09, zorder=0))
    ax.text(1.08, 3.88, f'crosses upward\n{ncross} of {nbelow} substrates',
            fontsize=8.4, color=ORANGE, fontweight='bold', va='top')
    ax.set_xlim(1, 4)
    ax.set_ylim(1, 4)
    ax.set_xlabel('Control PSI')
    ax.set_ylabel('UBR3-KO PSI')
    ax.grid(zorder=0)
    ax.legend(loc='lower right', markerscale=1.3)
    title(ax, f'{ncross} of the {nbelow} substrates below the line cross it',
          'Dotted diagonal = no change. Shaded quadrant = unstable in control, stable in UBR3 KO.')
    panel_tag(ax, 'B', C)


def p6c(ax, C):
    lib, CUT = C['lib'], C['CUT']
    edges = [1, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    cats = pd.cut(lib.mean_PSI_control, edges)
    g = lib.groupby(cats, observed=True)
    means, sems, ns = g.mean_dPSI.mean(), g.mean_dPSI.sem(), g.size()
    x = np.arange(len(means))
    cols = [VIOLET if (i.left + i.right) / 2 < CUT else AQUA for i in means.index]
    ax.bar(x, means.values, yerr=sems.values, color=cols, width=0.62, zorder=3,
           edgecolor=SURFACE, linewidth=1.1, error_kw=dict(ecolor=INK2, lw=1.2))
    ax.axhline(0, color=AXIS, lw=1)
    cut_x = next(i - 0.5 + (CUT - iv.left) / (iv.right - iv.left)
                 for i, iv in enumerate(means.index) if iv.left <= CUT < iv.right)
    ax.axvline(cut_x, color=INK, lw=1.4, ls='--', zorder=5)
    ax.text(cut_x, 0.985, f'cut, PSI {CUT:g}', fontsize=8.2, color=INK, fontweight='bold',
            ha='center', va='top', transform=ax.get_xaxis_transform())
    span = means.max() - means.min()
    ax.set_ylim(means.min() - span * 0.34, means.max() + span * 0.20)
    for i, (v, n) in enumerate(zip(means.values, ns.values)):
        off = span * 0.045 if v >= 0 else -span * 0.045
        ax.text(i, v + off, f'{v:+.3f}', ha='center', fontsize=8,
                va='bottom' if v >= 0 else 'top', fontweight='bold')
        ax.text(i, 0.015, f'n={n:,}', ha='center', va='bottom', fontsize=7.4, color=MUTED,
                transform=ax.get_xaxis_transform())
    ax.set_xticks(x)
    ax.set_xticklabels([f'{i.left:g}-{i.right:g}' for i in means.index], fontsize=8.4)
    ax.set_xlabel('Control PSI bin')
    ax.set_ylabel('Mean $\\Delta$PSI $\\pm$ SEM')
    ax.grid(axis='y', zorder=0)
    title(ax, 'Stabilisation peaks near PSI 2.75 and reverses above 3.5',
          'Peptides that are already maximally stable have no headroom left to gain')
    panel_tag(ax, 'C', C)


def p6d(ax, C):
    lib, CUT, lo = C['lib'], C['CUT'], C['lo']
    labs = [f'PSI < {CUT:g}\n(unstable)', f'PSI $\\geq$ {CUT:g}\n(stable)']
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
    title(ax, f'Substrates are {rates[0] / rates[1]:.1f}x more common below the line',
          'Error bars are Wilson 95% confidence intervals; counts below each bar')
    panel_tag(ax, 'D', C)


# ============================================================ FIGURE 9
def p9a(ax, C):
    """2x2 mosaic of the motif-bearing peptides."""
    lib, CUT = C['lib'], C['CUT']
    m = lib[lib.motif_class == '[P/G]-[E/D]']
    cells = [('unstable', 'yes'), ('unstable', 'no'), ('stable', 'yes'), ('stable', 'no')]
    counts = {(s, u): int(((m.stability_class == s) & (m.UBR3_substrate == u)).sum())
              for s, u in cells}
    n_un = counts[('unstable', 'yes')] + counts[('unstable', 'no')]
    n_st = counts[('stable', 'yes')] + counts[('stable', 'no')]
    total = n_un + n_st
    wu, ws_ = n_un / total, n_st / total
    x0 = 0.0
    for stab, wcol in [('unstable', wu), ('stable', ws_)]:
        k = counts[(stab, 'yes')]
        n = counts[(stab, 'yes')] + counts[(stab, 'no')]
        frac = k / n
        ax.add_patch(Rectangle((x0, 0), wcol - 0.012, frac, facecolor=ORANGE,
                               edgecolor=SURFACE, lw=2, zorder=3))
        ax.add_patch(Rectangle((x0, frac), wcol - 0.012, 1 - frac, facecolor='#dedcd4',
                               edgecolor=SURFACE, lw=2, zorder=3))
        cx = x0 + (wcol - 0.012) / 2
        # a thin substrate band cannot hold its own label - park it just above
        if frac > 0.11:
            ax.text(cx, frac / 2, f'{k}', ha='center', va='center', fontsize=17,
                    fontweight='bold', color='white', zorder=5)
        else:
            ax.text(cx, frac + 0.018, f'{k}', ha='center', va='bottom', fontsize=17,
                    fontweight='bold', color=ORANGE, zorder=5)
        ax.text(cx, frac + (1 - frac) / 2, f'{n - k}', ha='center', va='center',
                fontsize=17, fontweight='bold', color=INK2, zorder=5)
        ax.text(cx, -0.045,
                f'{stab}\ncontrol PSI {"<" if stab == "unstable" else "$\\geq$"} {CUT:g}\n'
                f'n = {n}   ({100 * frac:.1f}% substrates)',
                ha='center', va='top', fontsize=9.2, color=INK)
        x0 += wcol
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.46, 1.02)
    ax.axis('off')
    ax.legend(handles=[Patch(color=ORANGE, label='UBR3 substrate'),
                       Patch(color='#dedcd4', label='Not a substrate')],
              loc='lower center', bbox_to_anchor=(0.5, 0.0), ncol=2, columnspacing=1.8)
    title(ax, f'The {total} [P/G]-[E/D] peptides, classified two ways',
          'Column width is the share of peptides in that stability class; height is substrate share')
    panel_tag(ax, 'A', C)


def p9b(ax, C):
    """Substrate rate per stratum for each motif group."""
    lib = C['lib']
    groups = [('[P/G]-[E/D]', lib[lib.motif_class == '[P/G]-[E/D]']),
              ('[P/G]-other', lib[lib.motif_class == '[P/G]-other']),
              ('non-P/G', lib[lib.motif_class == 'non-P/G'])]
    w, x = 0.36, np.arange(len(groups))
    for j, (stab, col) in enumerate([('unstable', VIOLET), ('stable', AQUA)]):
        rates, labs = [], []
        for _, g in groups:
            s = g[g.stability_class == stab]
            rates.append(100 * s.is_UBR3_substrate.mean())
            labs.append(f'{int(s.is_UBR3_substrate.sum())}/{len(s):,}')
        pos = x + (j - 0.5) * w
        ax.bar(pos, rates, width=w * 0.92, color=col, zorder=3, edgecolor=SURFACE,
               linewidth=1.1, label=stab)
        for xi, r, lb in zip(pos, rates, labs):
            ax.text(xi, r + 0.3, f'{r:.2f}%', ha='center', fontsize=8.2, fontweight='bold')
            ax.text(xi, 0.015, lb, ha='center', va='bottom', fontsize=7.2, color=MUTED,
                    transform=ax.get_xaxis_transform())
    ax.set_xticks(x)
    ax.set_xticklabels([g for g, _ in groups], fontsize=9.4)
    ax.set_ylabel('UBR3 substrates (% of that cell)')
    ax.set_ylim(-2.0, 16.5)
    ax.grid(axis='y', zorder=0)
    ax.legend(title='Baseline stability', title_fontsize=8.5, loc='upper right')
    title(ax, 'Substrate rate in each stability x motif cell',
          'Counts under each bar are substrates / peptides in that cell')
    panel_tag(ax, 'B', C)


def p9c(ax, C):
    """The 16 motif substrates: where they start and where they end."""
    lib, CUT = C['lib'], C['CUT']
    m = lib[(lib.motif_class == '[P/G]-[E/D]') & lib.is_UBR3_substrate].copy()
    m = m.sort_values('mean_PSI_control')
    y = np.arange(len(m))
    for yi, (_, r) in zip(y, m.iterrows()):
        col = VIOLET if r.stability_class == 'unstable' else AQUA
        ax.annotate('', xy=(r.mean_PSI_UBR3KO, yi), xytext=(r.mean_PSI_control, yi),
                    arrowprops=dict(arrowstyle='-|>', color=col, lw=2.0,
                                    shrinkA=0, shrinkB=0), zorder=3)
        ax.scatter([r.mean_PSI_control], [yi], s=34, color=col, zorder=4,
                   edgecolors=SURFACE, linewidths=0.8)
    ax.axvline(CUT, color=INK, lw=1.4, ls='--', zorder=5)
    ax.text(CUT + 0.03, len(m) - 0.4, f'PSI {CUT:g}', fontsize=8.4, color=INK, fontweight='bold')
    ax.set_yticks(y)
    ax.set_yticklabels(m.Gene_ID, fontsize=8.6)
    for t, s in zip(ax.get_yticklabels(), m.stability_class):
        t.set_color(VIOLET if s == 'unstable' else AQUA)
        t.set_fontweight('bold')
    ax.set_xlabel('PSI   (arrow runs from control to UBR3 KO)')
    ax.set_xlim(1.2, 4.0)
    ax.set_ylim(-0.8, len(m) - 0.2)
    ax.grid(axis='x', zorder=0)
    ax.legend(handles=[Line2D([], [], color=VIOLET, lw=2.4, label='starts unstable'),
                       Line2D([], [], color=AQUA, lw=2.4, label='starts stable')],
              loc='lower right')
    title(ax, 'Every motif-bearing substrate gains stability',
          'Each arrow runs from control PSI to UBR3-KO PSI; all 16 point right')
    panel_tag(ax, 'C', C, xoff=-72)


def p9d(ax, C):
    """How the motif-bearing peptides split, as a flow of counts."""
    lib = C['lib']
    m = lib[lib.motif_class == '[P/G]-[E/D]']
    steps = [
        ('[P/G]-[E/D] peptides', len(m), MUTED),
        ('unstable at baseline', int((m.stability_class == 'unstable').sum()), VIOLET),
        ('  ... and a UBR3 substrate',
         int(((m.stability_class == 'unstable') & m.is_UBR3_substrate).sum()), ORANGE),
        ('stable at baseline', int((m.stability_class == 'stable').sum()), AQUA),
        ('  ... and a UBR3 substrate',
         int(((m.stability_class == 'stable') & m.is_UBR3_substrate).sum()), ORANGE),
    ]
    ys = np.arange(len(steps))[::-1]
    for yv, (lab, n, c) in zip(ys, steps):
        ax.barh(yv, n, color=c, height=0.6, zorder=3, edgecolor=SURFACE, linewidth=1.0)
        ax.text(n + 3, yv, f'{n}', va='center', fontweight='bold', fontsize=11.5, color=INK)
        ax.text(2, yv + 0.40, lab, va='bottom', ha='left', fontsize=8.8, color=INK2)
    ax.set_yticks([])
    ax.set_xlim(0, len(m) * 1.13)
    ax.set_ylim(-0.7, len(steps) - 0.15)
    ax.set_xlabel('Number of peptides')
    ax.grid(axis='x', zorder=0)
    ax.spines['left'].set_visible(False)
    title(ax, 'Where the 179 motif-bearing peptides end up',
          '11 of the 16 substrates start unstable; the other 5 were already stable and '
          'gained more')
    panel_tag(ax, 'D', C)


# ============================================================ FIGURE 10
# Downstream of the motif: what separates motif-bearing substrates from
# motif-bearing non-substrates at positions 4-24, where both sets are identical
# by construction at positions 1-3.
def _down(C):
    if 'down' not in C:
        m = C['lib'][C['lib'].motif_class == '[P/G]-[E/D]']
        A = list(m[m.is_UBR3_substrate].peptide_24mer)
        B = list(m[~m.is_UBR3_substrate].peptide_24mer)
        C['down'] = U.downstream_compare(A, B) + (len(A), len(B))
    return C['down']


def p10a(ax, C):
    per_res, per_cls, _, na, nb = _down(C)
    piv = per_cls.pivot(index='position', columns='group', values='log2_enrichment')[U.CLASS_ORDER]
    qv = per_cls.pivot(index='position', columns='group', values='q_value_BH')[U.CLASS_ORDER]
    m, q = piv.T.values, qv.T.values
    v = np.nanmax(np.abs(m))
    im = ax.imshow(m, cmap=DIVERGE, norm=TwoSlopeNorm(0, -v, v), aspect='auto')
    ax.set_yticks(range(len(U.CLASS_ORDER)))
    ax.set_yticklabels(U.CLASS_ORDER, fontsize=9)
    for t, c in zip(ax.get_yticklabels(), U.CLASS_ORDER):
        t.set_color(CLASS_COLOR[c])
        t.set_fontweight('bold')
    pos = list(piv.index)
    ax.set_xticks(range(len(pos)))
    ax.set_xticklabels(pos, fontsize=7.8)
    ax.set_xlabel('Position in the 24-mer')
    for i in range(m.shape[0]):
        for k in range(m.shape[1]):
            if q[i, k] < 0.05:
                ax.add_patch(Rectangle((k - 0.5, i - 0.5), 1, 1, fill=False,
                                       edgecolor=INK, lw=2.2, zorder=6))
                ax.text(k, i, '*', ha='center', va='center', fontsize=15,
                        fontweight='bold', color=INK, zorder=7)
    cb = ax.get_figure().colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label('log$_2$ (substrate / non-substrate)', fontsize=8.5)
    cb.ax.tick_params(labelsize=7.5)
    nsig = int((q < 0.05).sum())
    title(ax, 'Chemical class by position, downstream of the motif',
          f'{na} motif-bearing substrates vs {nb} motif-bearing non-substrates. Boxed * = q < 0.05 '
          f'({nsig} of {q.size} cells): Basic at position 5.')
    panel_tag(ax, 'A', C, xoff=-70)


def p10b(ax, C):
    """Volcano over all position x residue cells - the null result, shown honestly."""
    per_res, _, _, na, nb = _down(C)
    x = per_res.log2_enrichment.values
    y = -np.log10(per_res.p_value.values)
    sig = per_res.q_value_BH.values < 0.05
    ax.scatter(x[~sig], y[~sig], s=16, color=MUTED, alpha=0.55, linewidths=0, zorder=3,
               label=f'not significant (n={int((~sig).sum())})')
    if sig.any():
        ax.scatter(x[sig], y[sig], s=44, color=ORANGE, edgecolors=SURFACE, linewidths=0.8,
                   zorder=4, label=f'q < 0.05 (n={int(sig.sum())})')
    ax.axhline(-np.log10(0.05), color=RED, lw=1.1, ls=':', zorder=2)
    ax.text(ax.get_xlim()[1], -np.log10(0.05) + 0.06, 'p = 0.05', ha='right',
            fontsize=7.8, color=RED)
    ax.axvline(0, color=AXIS, lw=1)
    placed = []
    for _, r in per_res.nsmallest(4, 'p_value').iterrows():
        x0, y0 = r.log2_enrichment, -np.log10(r.p_value)
        for dx, dy in [(7, 3), (7, -11), (-7, 5), (-7, -12)]:
            if all(abs(x0 + dx / 60 - px) > 0.45 or abs(y0 + dy / 60 - py) > 0.22
                   for px, py in placed):
                break
        placed.append((x0 + dx / 60, y0 + dy / 60))
        ax.annotate(f'{r.group}{int(r.position)}', (x0, y0), textcoords='offset points',
                    xytext=(dx, dy), fontsize=8.2, fontweight='bold', color=INK,
                    ha='left' if dx > 0 else 'right')
    n_nom = int((per_res.p_value < 0.05).sum())
    exp = 0.05 * len(per_res)
    ax.set_xlabel('log$_2$ (substrate / non-substrate)')
    ax.set_ylabel('$-$log$_{10}$ $p$   (Fisher exact)')
    ax.grid(zorder=0)
    ax.legend(loc='upper left')
    title(ax, 'No individual residue distinguishes them',
          f'{n_nom} of {len(per_res)} cells reach p < 0.05 - fewer than the {exp:.0f} '
          'expected by chance; none survive FDR')
    panel_tag(ax, 'B', C)


def p10c(ax, C):
    """Aggregate class composition over the whole downstream window."""
    from scipy import stats as st
    _, _, pep, na, nb = _down(C)
    A = pep[pep.group == 'UBR3 substrate']
    B = pep[pep.group == 'not a substrate']
    x = np.arange(len(U.CLASS_ORDER))
    w = 0.36
    for j, (df, lab, col) in enumerate([(A, f'substrates (n={na})', ORANGE),
                                        (B, f'non-substrates (n={nb})', MUTED)]):
        vals = [df[f'n_{c}'].mean() for c in U.CLASS_ORDER]
        errs = [df[f'n_{c}'].sem() for c in U.CLASS_ORDER]
        ax.bar(x + (j - 0.5) * w, vals, yerr=errs, width=w * 0.92, color=col, zorder=3,
               edgecolor=SURFACE, linewidth=1.1, label=lab,
               error_kw=dict(ecolor=INK2, lw=1.1))
    for i, c in enumerate(U.CLASS_ORDER):
        p = st.mannwhitneyu(A[f'n_{c}'], B[f'n_{c}'])[1]
        hi = max(A[f'n_{c}'].mean(), B[f'n_{c}'].mean())
        ax.text(i, hi + 0.55, f'p = {p:.3f}' if p >= 0.001 else f'p = {p:.0e}',
                ha='center', fontsize=7.8,
                fontweight='bold' if p < 0.05 else 'normal',
                color=INK if p < 0.05 else MUTED)
    ax.set_xticks(x)
    ax.set_xticklabels(U.CLASS_ORDER, rotation=20, ha='right', fontsize=8.8)
    ax.set_ylabel('Residues per peptide, positions 4-24')
    ax.set_ylim(0, 8.6)
    ax.grid(axis='y', zorder=0)
    ax.legend(loc='upper right')
    title(ax, 'Only basic residues differ, and only in aggregate',
          'Mean per peptide over positions 4-24; bars are SEM, p from Mann-Whitney')
    panel_tag(ax, 'C', C)


def p10d(ax, C):
    """Net charge downstream - the headline separation."""
    from scipy import stats as st
    _, _, pep, na, nb = _down(C)
    A = pep[pep.group == 'UBR3 substrate'].net_charge.values
    B = pep[pep.group == 'not a substrate'].net_charge.values
    bp = ax.boxplot([A, B], patch_artist=True, widths=0.45, showfliers=False,
                    medianprops=dict(color=INK, lw=1.8),
                    whiskerprops=dict(color=INK2, lw=1), capprops=dict(color=INK2, lw=1),
                    boxprops=dict(edgecolor=SURFACE, lw=1.2))
    for patch, c in zip(bp['boxes'], [ORANGE, MUTED]):
        patch.set_facecolor(c)
    rng = np.random.default_rng(0)
    for i, v in enumerate([A, B], 1):
        ax.scatter(i + (rng.random(len(v)) - 0.5) * 0.28, v, s=13, color=INK2,
                   alpha=0.4, linewidths=0, zorder=4)
    ax.axhline(0, color=AXIS, lw=1)
    p = st.mannwhitneyu(A, B)[1]
    ax.set_xticks([1, 2])
    ax.set_xticklabels([f'UBR3 substrate\n(n={na})', f'not a substrate\n(n={nb})'], fontsize=9.2)
    ax.set_ylabel('Net charge over positions 4-24')
    ax.grid(axis='y', zorder=0)
    top = max(A.max(), B.max())
    ax.plot([1, 1, 2, 2], [top + 0.7, top + 1.3, top + 1.3, top + 0.7], color=INK2, lw=1.1)
    ax.text(1.5, top + 1.5, f'Mann-Whitney p = {p:.4f}', ha='center', fontsize=9,
            fontweight='bold', color=INK)
    ax.text(0.03, 0.03, f'means:  {A.mean():+.2f}   vs   {B.mean():+.2f}',
            transform=ax.transAxes, fontsize=9, color=INK)
    ax.set_ylim(min(A.min(), B.min()) - 1, top + 2.6)
    title(ax, 'Substrates carry a net positive charge downstream',
          'Both groups carry the motif - only the downstream window differs')
    panel_tag(ax, 'D', C)


# ============================================================ FIGURE 11
# The 4A/4B stacked-composition view, but for the P/G substrates against both
# the whole library and the P/G-matched background.
def _pg(C):
    if 'pg' not in C:
        lib, hit = C['lib'], C['hit']
        sub = list(hit[hit.is_PG].peptide_24mer)
        pglib = list(lib[lib.is_PG].peptide_24mer)
        C['pg'] = (sub, pglib, U.class_matrix(sub), U.class_matrix(pglib))
    return C['pg']


def p11a(ax, C):
    sub, _, cs, _ = _pg(C)
    _stack4(ax, C, cs, f'{len(sub)} UBR3 substrates with Pro/Gly at position 2', 'A',
            'Position 2 is Pro or Gly by definition (14 Gly, 12 Pro); position 3 is 42% Acidic',
            legend=True)


def p11b(ax, C):
    _, pglib, _, cl = _pg(C)
    _stack4(ax, C, cl, f'{len(pglib):,} library peptides with Pro/Gly at position 2', 'B',
            'The motif-matched background: same position-2 constraint, no selection for '
            'UBR3 dependence', legend=C.get('standalone', False))


def p11c(ax, C):
    _stack4(ax, C, C['cl'], f'{len(C["ls"]):,} library peptides (whole library)', 'C',
            'The global background, for reference - position 2 here is unconstrained',
            legend=C.get('standalone', False))


def p11d(ax, C):
    """P/G substrates vs the P/G-matched library, FDR-masked."""
    sub, pglib, _, _ = _pg(C)
    lr, qv, _ = U.enrichment_test(sub, pglib, groups=U.CLASS_MEMBERS)
    m, q = lr[U.CLASS_ORDER].T.values, qv[U.CLASS_ORDER].T.values
    shown = np.where(q < 0.05, m, np.nan)
    v = np.nanmax(np.abs(shown)) if np.isfinite(shown).any() else np.nanmax(np.abs(m))
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
    for i in range(m.shape[0]):
        for k in range(m.shape[1]):
            if q[i, k] < 0.05:
                ax.text(k, i, f'{m[i, k]:+.1f}', ha='center', va='center', fontsize=6.4,
                        color='white' if abs(m[i, k]) > v * 0.55 else INK, fontweight='bold')
    cb = ax.get_figure().colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label('log$_2$ (P/G substrates / P/G library)', fontsize=8.5)
    cb.ax.tick_params(labelsize=7.5)
    nsig = int((q < 0.05).sum())
    title(ax, 'P/G substrates vs the P/G-matched library (FDR-masked)',
          f'Beige = not significant. {nsig} of {q.size} cells survive q < 0.05 - only Acidic at '
          'position 3. Beyond position 3 the substrates match their own background.')
    panel_tag(ax, 'D', C, xoff=-70)


# ---------------------------------------------------------------- registry
PANELS = {
    '2A': p2a, '2B': p2b, '2C': p2c, '2D': p2d,
    '3A': p3a, '3B': p3b, '3C': p3c,
    '4A': p4a, '4B': p4b, '4C': p4c, '4D': p4d,
    '5A': p5a, '5B': p5b,
    '6A': p6a, '6B': p6b, '6C': p6c, '6D': p6d,
    '9A': p9a, '9B': p9b, '9C': p9c, '9D': p9d,
    '10A': p10a, '10B': p10b, '10C': p10c, '10D': p10d,
    '11A': p11a, '11B': p11b, '11C': p11c, '11D': p11d,
}

# aspect ratio hint per panel for standalone rendering (width, height) in inches
PANEL_SIZE = {
    '2A': (9.8, 6.4), '2B': (9.2, 6.6), '2C': (9.8, 6.0), '2D': (8.8, 7.2),
    '3A': (14.0, 5.2), '3B': (14.0, 5.2), '3C': (14.0, 5.6),
    '4A': (14.0, 5.4), '4B': (14.0, 5.0), '4C': (11.0, 5.2), '4D': (9.0, 6.4),
    '5A': (9.4, 8.6), '5B': (9.4, 8.6),
    '6A': (9.6, 6.4), '6B': (8.8, 7.2), '6C': (9.4, 6.4), '6D': (8.4, 6.6),
    '9A': (9.0, 6.8), '9B': (9.4, 6.4), '9C': (8.8, 7.2), '9D': (9.8, 5.8),
    '10A': (12.4, 5.2), '10B': (9.4, 6.6), '10C': (9.6, 6.4), '10D': (8.6, 6.8),
    '11A': (14.0, 5.4), '11B': (14.0, 5.0), '11C': (14.0, 5.0), '11D': (12.6, 5.2),
}
