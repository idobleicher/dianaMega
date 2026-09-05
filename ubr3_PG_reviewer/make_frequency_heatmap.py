"""
Position x residue / position x class heatmaps coloured by FREQUENCY.

Companion to make_significance_heatmap.py (Figures 16-17), which colours by
-log10 p, and to the original Colab figure, which colours by fold change. This
one colours by the plainest quantity in the table: the percentage of peptides in
a group that carry that residue (or that chemical class) at that position.

Two differences from Figures 16-17, both requested:
  * rows are ALPHABETICAL -- A C D E F G H I K L M N P Q R S T V W Y for the
    residues, Acidic ... Polar for the classes -- not grouped by class. The
    class of each row is still shown as a colour chip beside its label, keyed
    to the same palette the logos use, so the grouping is recoverable by eye.
  * residues and classes are in SEPARATE figures with their own colour scales.
    They were separate figures before too, but they shared one -log10 p scale;
    a class is the union of 2-5 residues and is necessarily commoner, so a
    shared frequency scale would flatten the residue panel.

Data: pg_motif_data.py -- the shared loader. n = 20 P/G-D/E/T substrates vs
n = 193 non-substrate controls carrying the same motif, positions 4-24
(positions 1-3 are the motif itself).

READ THE PAIRED PANELS TOGETHER. A frequency heatmap on its own says nothing
about UBR3: the loudest cells are partly just the commonest amino acids. Gly at
position 15 is 30% of substrates, but Gly is common among these controls too
(10.9%); Arg at 7 is 35% against 6.2%. The control panel (Figure18b / 19b) is
the same matrix for the 193 controls on the same colour scale, and is what
separates those two cases. Enrichment is Figures 13-17's job, not this one's.

ANNOTATION matches Figures 16-17 exactly, from the same two constants: cell
text is WHITE on every cell with a thin dark halo (pg.CELL_TEXT /
pg.cell_text_effects), and the row labels are one step larger. The halo carries
more weight here than there -- most cells of this figure sit at the pale end of
the ramp, and a plain white number on a 5% cell would not be visible at all.

WINDOWS -- every panel is written twice, over the full POS (4-24) and over
pg.POS_N12 (4-12, out to the 12th amino acid), on the SAME colour scale, so the
crop changes what is shown and nothing about how it is read.

Outputs (figures/, 600 dpi PNG + vector PDF + SVG), each also as ..._pos4_12:
  Figure18_frequency_heatmap_residues            % of the 20 substrates
  Figure18b_frequency_heatmap_residues_controls  % of the 193 controls, same scale
  Figure19_frequency_heatmap_categories          % of the 20 substrates, 6 classes
  Figure19b_frequency_heatmap_categories_controls
  data/frequency_heatmap_cells.csv               every cell, both groups
"""
import os
import textwrap

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import pg_motif_data as pg

HERE = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(HERE, "figures")
DATA_DIR = os.path.join(HERE, "data")
os.makedirs(FIG_DIR, exist_ok=True)

POS = pg.POS
N_SUB, N_CTRL = pg.N_SUB, pg.N_CTRL
CAT_COLORS, CAT_OF_AA = pg.CAT_COLORS, pg.CAT_OF_AA


def fig_width(ncol):
    """Constant cell width across windows -- the same rule Figures 16-17 use,
    so a cropped panel is the full one with columns removed, not squeezed."""
    return 3.0 + 0.495 * ncol


def wrap(text, figw, fontsize=5.5):
    return textwrap.fill(text, width=max(60, int(figw * 25)))

# ---- alphabetical row orders (the whole point of this variant) --------------
AA_ALPHA = sorted(pg.AA)          # A C D E F G H I K L M N P Q R S T V W Y
CAT_ALPHA = sorted(pg.CAT_ORDER)  # Acidic Aliphatic Aromatic Basic Hydrophobic Polar
assert AA_ALPHA == list("ACDEFGHIKLMNPQRSTVWY")

cells = pg.load()
cells[["kind", "label", "position", "n_sub", "pct_sub", "n_ctrl", "pct_ctrl",
       "fold_change", "p_chi2", "p_fisher", "q_chi2", "q_fisher"]].to_csv(
    os.path.join(DATA_DIR, "frequency_heatmap_cells.csv"), index=False)

# ---------------------------------------------------------------- style
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8, "axes.linewidth": 0.8,
    "svg.fonttype": "none", "pdf.fonttype": 42,
    "figure.dpi": 600, "savefig.dpi": 600, "savefig.bbox": "tight",
})

# Same warm ramp as Figures 16-17 (which reuses the logo palette #E8A33D /
# #C0392B / #7B241C), so every figure in the set reads as one system. Here it
# is a plain sequential magnitude scale: near-white = 0% of the group carries
# this residue here, dark red = the top of the scale.
CMAP = LinearSegmentedColormap.from_list("freq_warm", [
    "#FFFDF8", "#FDF3D4", "#FAE3A0", "#F5C463", "#E8A33D",
    "#DD7B2F", "#C0392B", "#96271F", "#7B241C"])


def pct_text(v):
    """Cell annotation. 0 is left blank so the shape of the matrix stays
    legible; a control cell can sit below 1% (1/193 = 0.52%) and must not be
    rounded down to a printed 0."""
    if not np.isfinite(v) or v <= 0:
        return ""
    if v < 1:
        return "<1"
    return f"{v:.0f}"


def heatmap(kind, group, order, vmax, title, subtitle, stem, height,
            positions=None, annot_size=6.0, chips=True, legend=False):
    positions = list(positions or POS)
    window = positions != list(POS)
    figsize = (fig_width(len(positions)), height)
    pcol = "pct_sub" if group == "sub" else "pct_ctrl"
    ncell = "n_sub" if group == "sub" else "n_ctrl"
    n_group = N_SUB if group == "sub" else N_CTRL

    F = pg.matrix(pcol, kind=kind, order=order, cells=cells,
                  positions=positions).values.astype(float)
    PF = pg.matrix("p_fisher", kind=kind, order=order, cells=cells,
                   positions=positions).values.astype(float)

    C = np.clip(F, 0.0, vmax)
    nrow, ncols = C.shape

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(C, cmap=CMAP, aspect="auto", vmin=0.0, vmax=vmax,
                   extent=(-0.5, ncols - 0.5, nrow - 0.5, -0.5))

    for x in np.arange(-0.5, ncols, 1):
        ax.axvline(x, color="white", linewidth=0.7)
    for y in np.arange(-0.5, nrow, 1):
        ax.axhline(y, color="white", linewidth=0.7)

    # Every non-zero cell carries its percentage. On the substrate panel an
    # asterisk flags the cells that also differ between the two groups (Fisher
    # exact p < 0.05), so this figure never contradicts Figures 16-17 -- but
    # colour here is frequency alone and most dark cells are NOT marked.
    for r in range(nrow):
        for c in range(ncols):
            txt = pct_text(F[r, c])
            if not txt:
                continue
            if group == "sub" and np.isfinite(PF[r, c]) and PF[r, c] < 0.05:
                txt += "*"
            # White on every cell, as in Figures 16-17. Most cells of this
            # figure sit at the pale end of the ramp -- a 5% cell is barely
            # tinted -- so the white glyph carries the thin dark halo defined
            # in pg_motif_data. Without it a white number on a 5% cell would
            # not be there at all.
            ax.text(c, r, txt, ha="center", va="center", fontsize=annot_size,
                    fontweight="bold", color=pg.CELL_TEXT,
                    path_effects=pg.cell_text_effects())

    ax.set_xticks(range(ncols))
    ax.set_xticklabels(positions, fontsize=7)
    ax.set_yticks(range(nrow))
    # One step up, matching Figures 16-17: the residue letter is the key to
    # its row and should not be the quietest thing on the axis.
    ax.set_yticklabels(order, fontsize=12.0 if kind == "residue" else 10.5,
                       fontweight="bold" if kind == "residue" else "normal")
    ax.tick_params(length=0, pad=3)
    ax.set_xlabel("Position in the 24-mer", fontsize=8.5, labelpad=4)
    for sp in ax.spines.values():
        sp.set_visible(False)

    # Class colour chip outside each row label. The rows are alphabetical, so
    # this is the only thing left carrying the chemical grouping. The gap has
    # to clear the tick label itself, which is one character wide for the
    # residues and up to "Hydrophobic" wide for the classes.
    if chips:
        gap = 0.75 if kind == "residue" else 1.70
        for r, lab in enumerate(order):
            cat = CAT_OF_AA[lab] if kind == "residue" else lab
            ax.add_patch(Rectangle((-0.5 - gap, r - 0.34), 0.30, 0.68,
                                   facecolor=CAT_COLORS[cat], clip_on=False,
                                   edgecolor="white", linewidth=0.5, zorder=6))
        ax.set_xlim(-0.5, ncols - 0.5)

    # The class key goes UNDER the matrix, where it competes with nothing; the
    # top of the figure stays title -> subtitle -> matrix on all four panels.
    if legend:
        handles = [Rectangle((0, 0), 1, 1, facecolor=CAT_COLORS[c],
                             edgecolor="white", linewidth=0.5) for c in CAT_ALPHA]
        labels = [f"{c}  {''.join(sorted(pg.CATEGORY_MEMBERS[c]))}" for c in CAT_ALPHA]
        ax.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.075),
                  ncol=6, frameon=False, fontsize=6.2, handlelength=0.9,
                  handleheight=0.9, handletextpad=0.35, columnspacing=1.6)

    cbar = fig.colorbar(im, ax=ax, fraction=0.022, pad=0.015)
    cbar.ax.tick_params(labelsize=6.5, length=2, width=0.6)
    cbar.outline.set_linewidth(0.6)
    ticks = [t for t in (0, 5, 10, 15, 20, 25, 30, 35, 40, 45) if t <= vmax]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{t}%" for t in ticks])
    # The label runs along a colourbar as tall as the axes. On the 3-in-tall
    # class panels that bar is shorter than the label, so the text is wrapped
    # to two lines there instead of overrunning the title.
    cbar_label = ("% of the group carrying this "
                  + ("residue" if kind == "residue" else "class")
                  + " at this position")
    if height < 5:
        cbar_label = textwrap.fill(cbar_label, width=30)
    cbar.set_label(cbar_label, fontsize=7, labelpad=6, linespacing=1.4)

    ax.set_title(title, fontweight="bold", fontsize=9.5, pad=16)
    ax.text(0.5, 1.012, subtitle, transform=ax.transAxes,
            ha="center", va="bottom", fontsize=6.4, color="#666", style="italic")

    step = 100.0 / n_group
    who = ("P/G-D/E/T substrates" if group == "sub"
           else "non-substrate controls carrying the same motif")
    other = "controls" if group == "sub" else "substrates"
    # counts for THIS family, not the other one
    fam = cells[(cells.kind == kind) & cells.defined]
    n_test, n_hit = len(fam), int((fam.p_fisher < 0.05).sum())
    n_fdr = int((fam.q_fisher < 0.05).sum())
    tail = (f"* = the cell also differs between the two groups at Fisher exact "
            f"p < 0.05 ({n_hit} of {n_test} testable "
            f"{'residue' if kind == 'residue' else 'class'} cells, against "
            f"{0.05 * n_test:.0f} expected by chance alone"
            + ("; none survives FDR correction" if not n_fdr else
               f"; {n_fdr} survives FDR correction" if n_fdr == 1 else
               f"; {n_fdr} survive FDR correction") + ").  "
            if group == "sub" else
            "This is the background panel: the same matrix, on the same colour "
            "scale, for the controls.  ")
    companion = "18b / 19b" if group == "sub" else "18 / 19"
    if window:
        companion += " for positions 4–12"
    lines = [
        f"Colour and cell text are the same number: the percentage of the "
        f"n = {n_group} {who} with that "
        f"{'residue' if kind == 'residue' else 'class'} at that position "
        f"(one peptide = {step:.2g}%).   Blank = 0%, no peptide in the group.   "
        f"Rows are alphabetical; the chip beside each label is its chemical class.",

        f"This is composition, NOT enrichment — a residue frequent here may be "
        f"just as frequent among the {other}.  " + tail
        + f"Compare with Figure {companion} "
          f"for the other group, and with Figures 16–17 for significance.",
    ]
    if window:
        lines.append(
            f"This panel is the N-terminal window only — positions "
            f"{positions[0]}–{positions[-1]}, out to the {positions[-1]}th "
            f"amino acid of the 24-mer — on the same 0–{vmax:.0f}% scale as "
            f"the full {POS[0]}–{POS[-1]} version, so the two are read the "
            f"same way.")
    ax.text(0.5, -0.145 if kind == "residue" else -0.36,
            "\n".join(wrap(l, figsize[0]) for l in lines),
            transform=ax.transAxes, ha="center", va="top", fontsize=5.5,
            color="#666")

    plt.tight_layout(pad=0.4)
    for ext in ("png", "pdf", "svg"):
        path = os.path.join(FIG_DIR, f"{stem}.{ext}")
        fig.savefig(path, format=ext)
        print(f"  saved: {os.path.basename(path)}  ({os.path.getsize(path)/1024:.1f} KB)")
    plt.close(fig)


# ---------------------------------------------------------------- figures
# One scale per family, shared by that family's substrate and control panel so
# the two are directly comparable. Residues top out at 35% (R at 7), classes at
# 45% (Polar at 12): a class is a union of residues and is necessarily commoner.
VMAX_RES, VMAX_CAT = 35.0, 45.0
for _kind, _vmax in (("residue", VMAX_RES), ("category", VMAX_CAT)):
    _m = cells[cells.kind == _kind][["pct_sub", "pct_ctrl"]].max().max()
    assert _m <= _vmax + 1e-9, f"{_kind}: {_m:.1f}% exceeds the {_vmax:.0f}% scale top"

print("building frequency heatmaps\n")

# Each panel twice: the full 4-24 matrix and the N-terminal window out to the
# 12th amino acid, on the SAME colour scale, so the crop changes what is shown
# and nothing about how it is read.
WINDOWS = [(POS, "", ""),
           (pg.POS_N12, "_pos4_12", " (positions 4–12)")]

for positions, suffix, wt in WINDOWS:
    heatmap("residue", "sub", AA_ALPHA, VMAX_RES,
            "Residue frequency by position — the 20 P/G-D/E/T substrates" + wt,
            "colour = % of substrates carrying that residue at that position · "
            "rows alphabetical · this is not an enrichment scale",
            "Figure18_frequency_heatmap_residues" + suffix, 8.0, positions,
            legend=True)

    heatmap("residue", "ctrl", AA_ALPHA, VMAX_RES,
            "Residue frequency by position — the 193 non-substrate controls "
            "(background panel)" + wt,
            "same matrix, same 0–35% scale, for the controls · what the "
            "substrate panel has to be read against",
            "Figure18b_frequency_heatmap_residues_controls" + suffix,
            8.0, positions, legend=True)

    heatmap("category", "sub", CAT_ALPHA, VMAX_CAT,
            "Chemical class frequency by position — the 20 P/G-D/E/T "
            "substrates" + wt,
            "colour = % of substrates with a residue of that class at that "
            "position · rows alphabetical · the six classes partition the 20 "
            "residues, so every column sums to 100%",
            "Figure19_frequency_heatmap_categories" + suffix, 3.0, positions,
            annot_size=6.6)

    heatmap("category", "ctrl", CAT_ALPHA, VMAX_CAT,
            "Chemical class frequency by position — the 193 non-substrate "
            "controls (background panel)" + wt,
            "same matrix, same 0–45% scale, for the controls",
            "Figure19b_frequency_heatmap_categories_controls" + suffix,
            3.0, positions, annot_size=6.6)

# ---------------------------------------------------------------- readout
print("\n" + "=" * 78)
print("frequency does not rank the same cells as significance")
print("=" * 78)
cols = ["position", "label", "n_sub", "pct_sub", "n_ctrl", "pct_ctrl",
        "fold_change", "p_fisher"]
fmt = {"pct_sub": "{:.1f}".format, "pct_ctrl": "{:.1f}".format,
       "fold_change": "{:.2f}".format, "p_fisher": "{:.2e}".format}
res = cells[cells.kind == "residue"]
top = res.nlargest(8, "pct_sub")
print("\ntop 8 residue cells by FREQUENCY IN SUBSTRATES (what this figure makes loudest):")
print(top[cols].to_string(index=False, formatters=fmt))
print(f"\n...of which p < 0.05 (Fisher): {int((top.p_fisher < 0.05).sum())} of 8")
print("\ntop 8 residue cells by FREQUENCY IN CONTROLS (the background it sits on):")
print(res.nlargest(8, "pct_ctrl")[cols].to_string(index=False, formatters=fmt))

cat = cells[cells.kind == "category"]
print("\ntop 6 class cells by frequency in substrates:")
print(cat.nlargest(6, "pct_sub")[cols].to_string(index=False, formatters=fmt))

print("\ncolumn sums (sanity: every position must be 100% in both groups)")
for kind in ("residue", "category"):
    k = cells[cells.kind == kind]
    for col, lab in (("pct_sub", "substrates"), ("pct_ctrl", "controls")):
        s = k.groupby("position")[col].sum()
        print(f"   {kind:<9} {lab:<11} min {s.min():.1f}%  max {s.max():.1f}%")

print(f"\nwrote data/frequency_heatmap_cells.csv ({len(cells)} cells)")
