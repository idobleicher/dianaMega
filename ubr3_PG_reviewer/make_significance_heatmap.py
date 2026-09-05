"""
Position x residue heatmaps coloured by SIGNIFICANCE instead of fold change.

The user's own Colab heatmap (embedded at the bottom of AMINO ACIDS.csv) colours
each cell by fold change on a locked 0-5 scale. That makes the largest fold
changes the loudest cells -- but in this dataset the largest fold changes are
the ones resting on the fewest observations. W at position 12 is 9.65x because
it is 1 substrate vs 1 non-substrate; R at position 7 is only 5.6x but it is
7/20 vs 12/193, which is where the actual signal is. Colouring by significance
inverts that ranking and puts the well-supported cells in front.

Data: pg_motif_data.py -- the shared loader. See that module for the source
sheet, the two tests, and the arginine correction.

n = 20 P/G-D/E/T substrates  vs  n = 193 non-substrate controls carrying the
same motif. Positions 4-24 (positions 1-3 are the motif itself).

Two p-values are carried for every cell:
  chi2   -- the author's, read verbatim from the workbook.
  Fisher -- recomputed here from the same 2x2 table. Chi-square needs expected
            counts >= 5; the workbook's own expected-count block shows most
            cells expect < 2 substrates, so Fisher exact is the valid test.
            Both are plotted so the difference is visible rather than assumed.
BH-FDR is applied across all 420 residue cells (126 category cells).

ANNOTATION -- colour is the p value, cell text is the EFFECT SIZE. Every cell
at p < 0.05 carries its fold change (substrates / controls) under the star bin
of its uncorrected p -- stars on their own line above the number, set a little
larger -- and nothing else. The two quantities sit on different
channels on purpose: p folds the effect size together with how many peptides
the cell rests on, so a reader given only colour cannot tell 9.65x resting on
one peptide from 5.63x resting on seven. The FDR survivors carry NO mark of
their own inside the matrix -- neither the old black box nor a dagger, both of
which put a third thing in a cell that already holds a colour and a number.
The caption names them instead, and every cell's exact p, q, both counts and
both percentages are in data/significance_heatmap_cells.csv.

That cell text is WHITE on every cell, pale ground or dark, with a thin dark
halo (pg.CELL_TEXT / pg.cell_text_effects) so the white keeps an edge at the
pale end of the ramp. It used to flip to near-black below 62% of the scale,
which put two text colours in one figure and made the flip read as a threshold
in the data. Row labels are one step larger for the same reason the text is one
colour: the residue letter is the key to its row.

Type is Helvetica throughout (Arial where Helvetica is absent; see HELV), and
no left-hand class colour panel is drawn -- the residue rows are still ordered
by chemical class, which the caption states in words.

WINDOWS -- every figure is written twice, over the full POS (4-24) and over
pg.POS_N12 (4-12, out to the 12th amino acid). The crop is a display crop: the
p values are per cell and unchanged, the colour scale is identical, the cells
are the same width in inches, and q stays the BH-FDR of the FULL family rather
than being recomputed inside the window -- recomputing would be choosing the
region after seeing the data. The cropped captions say both.

Outputs (figures/, 600 dpi PNG + vector PDF + SVG), each also as ..._pos4_12:
  Figure16_significance_heatmap_residues        -log10 chi2 p, warm scale
  Figure16b_significance_heatmap_residues_fisher  -log10 Fisher p
  Figure16c_significance_heatmap_residues_enriched_only  depleted cells forced to
                                                the floor; drop-in for the Colab
                                                figure's look
  Figure17_significance_heatmap_categories      the six chemical classes
  Figure17b_significance_heatmap_categories_fisher  same, Fisher exact
  data/significance_heatmap_cells.csv           every cell, both tests, q values
"""
import os
import re
import textwrap

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pg_motif_data as pg

HERE = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(HERE, "figures")
DATA_DIR = os.path.join(HERE, "data")
os.makedirs(FIG_DIR, exist_ok=True)

POS, AA = pg.POS, pg.AA
N_SUB, N_CTRL = pg.N_SUB, pg.N_CTRL
CATEGORY_MEMBERS, CAT_ORDER = pg.CATEGORY_MEMBERS, pg.CAT_ORDER
stars = pg.stars

cells = pg.load()
cells.to_csv(os.path.join(DATA_DIR, "significance_heatmap_cells.csv"), index=False)


def matrix(kind, col, order, positions=None):
    return pg.matrix(col, kind=kind, order=order, cells=cells,
                     positions=positions)


def fig_width(ncol):
    """Constant cell width across windows. The 21-position figure is 13.4 in
    wide; a 9-position one is narrower by the columns it drops rather than the
    same width with fatter cells, so the two are directly comparable."""
    return 3.0 + 0.495 * ncol


def wrap(text, figw, fontsize=5.6):
    """Wrap a caption line to the figure width. At 5.6 pt a character averages
    about 0.039 in, so roughly 25 fit per inch of figure."""
    return textwrap.fill(text, width=max(60, int(figw * 25)))


# ---------------------------------------------------------------- style
# Helvetica everywhere, including inside the maths. Helvetica itself is not
# installed on every machine, so the stack falls through to Arial, which is
# metrically identical to Helvetica and is the substitution every journal
# accepts. mathtext is pointed at the same stack so a "p" set in maths cannot
# come out in a different face -- and the strings below avoid mathtext anyway.
HELV = ["Helvetica", "Helvetica Neue", "Arial", "Nimbus Sans",
        "Liberation Sans", "DejaVu Sans"]
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": HELV,
    "mathtext.fontset": "custom", "mathtext.default": "regular",
    "font.size": 8, "axes.linewidth": 0.8,
    "svg.fonttype": "none", "pdf.fonttype": 42,
    "figure.dpi": 600, "savefig.dpi": 600, "savefig.bbox": "tight",
    "figure.facecolor": "white", "savefig.facecolor": "white",
})

# the first member of the stack that is actually installed -- mathtext resolves
# a family by name, so it has to be given a real one or it silently drops to
# DejaVu and the maths comes out in a different face from the labels
_INSTALLED = {f.name for f in mpl.font_manager.fontManager.ttflist}
FONT = next((f for f in HELV if f in _INSTALLED), "DejaVu Sans")
mpl.rcParams.update({"mathtext.rm": FONT, "mathtext.it": f"{FONT}:italic",
                     "mathtext.bf": f"{FONT}:bold"})

# The palette lives in pg_motif_data.py, next to the loader every figure reads,
# so the heatmaps and the significance logo cannot drift apart: they are one
# visual set keyed to one ramp.
CMAP_MAIN = pg.CMAP_SIG
CMAP_ENRICH = CMAP_MAIN
GREY = pg.NO_TEST        # residue absent from both groups -- no test possible
GRID = pg.GRID           # cell separator
INK = pg.INK             # body text
MUTED = pg.MUTED         # captions


def fmt_p(p):
    """p as it is written in a paper: 0.012, 4.8e-4."""
    if not np.isfinite(p):
        return "--"
    if p >= 1e-3:
        return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e")
    return f"{m}e{int(e)}"


def save(fig, stem):
    """PNG + PDF + SVG. The SVG keeps live text (svg.fonttype = none), so its
    font-family is rewritten to the full Helvetica stack: it then renders in
    real Helvetica wherever Helvetica exists and in Arial everywhere else."""
    for ext in ("png", "pdf", "svg"):
        path = os.path.join(FIG_DIR, f"{stem}.{ext}")
        fig.savefig(path, format=ext)
        if ext == "svg":
            with open(path, encoding="utf-8") as fh:
                svg = fh.read()
            # only where a single resolved family was written; the stack
            # matplotlib emits from font.sans-serif is already Helvetica-first
            svg = re.sub(r"font-family: ?'?" + FONT + r"'?(?=[;\"<])",
                         "font-family: Helvetica, " + FONT + ", sans-serif", svg)
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(svg)
        print(f"  saved: {os.path.basename(path)}  "
              f"({os.path.getsize(path)/1024:.1f} KB)")
    plt.close(fig)


CLASS_NOTE = ("Rows are ordered by chemical class — Basic K H R · "
              "Acidic D E · Aromatic F Y W · Aliphatic G P A M · "
              "Hydrophobic V L I · Polar S T C N Q.   ")


def heatmap(kind, pcol, qcol, order, signed, title, subtitle, stem,
            height, positions=None, vmax=5.0, annot_size=6.0,
            star_size=None):
    positions = list(positions or POS)
    star_size = annot_size + 1.8 if star_size is None else star_size
    window = positions != list(POS)
    figsize = (fig_width(len(positions)), height)
    P = matrix(kind, pcol, order, positions).values.astype(float)
    Q = matrix(kind, qcol, order, positions).values.astype(float)
    FC = matrix(kind, "fold_change", order, positions).values.astype(float)
    D = matrix(kind, "direction", order, positions).values.astype(float)
    OK = matrix(kind, "defined", order, positions).values.astype(float) == 1

    with np.errstate(divide="ignore", invalid="ignore"):
        L = -np.log10(P)
    L = np.where(np.isfinite(L), L, np.nan)
    # `signed` now means "colour every testable cell by its p"; the alternative
    # forces depleted cells to the floor so only enrichment carries colour.
    # With no depleted cell significant anywhere, the two differ only among
    # cells that are already n.s.
    C = L if signed else np.where(D > 0, L, 0.0)
    C = np.clip(C, 0.0, vmax)
    C = np.ma.masked_where(~OK | np.isnan(C), C)

    nrow, ncol = C.shape
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_facecolor(GREY)
    cmap = (CMAP_MAIN if signed else CMAP_ENRICH).copy()
    cmap.set_bad(GREY)
    im = ax.imshow(C, cmap=cmap, aspect="auto", vmin=0.0, vmax=vmax,
                   extent=(-0.5, ncol - 0.5, nrow - 0.5, -0.5))

    # Cell separators, as in the Colab figure. They stay white while the pale
    # end of the ramp moved to a warm off-white, so an n.s. cell now reads as a
    # quiet ground with a visible rule around it rather than as a hole in the
    # figure with an invisible one.
    for x in np.arange(-0.5, ncol, 1):
        ax.axvline(x, color=GRID, linewidth=0.8)
    for y in np.arange(-0.5, nrow, 1):
        ax.axhline(y, color=GRID, linewidth=0.8)

    # Cell text = FOLD CHANGE with its STAR BIN STACKED ABOVE IT, on every
    # cell at p < 0.05. Colour is the p value, so the number carries the effect
    # size the colour cannot.
    #
    # The stars sit on their own line above the number rather than trailing
    # it, and the layout lives in pg.annotate_cell -- shared with Figures
    # 18-19, so the two heatmap families annotate identically.
    # The FDR survivors get no mark of their own here -- the caption names
    # them, so nothing extra is laid over the colour.
    #
    # The text is WHITE on every cell, dark ground or pale. It used to flip to
    # near-black below 62% of the scale, which put two text colours in one
    # figure and made the flip itself read as a threshold in the data. A white
    # number needs an edge on a straw cell, so every glyph carries the thin
    # halo defined in pg_motif_data: the letterform stays white, and the same
    # white, everywhere in the matrix.
    for r in range(nrow):
        for c_ in range(ncol):
            p_, q, fc = P[r, c_], Q[r, c_], FC[r, c_]
            if not (np.isfinite(p_) and p_ < 0.05):
                continue
            fc_txt = "∞" if np.isinf(fc) else (
                "--" if not np.isfinite(fc) else f"{fc:.2f}")
            pg.annotate_cell(ax, c_, r, fc_txt, stars(p_),
                             size=annot_size, star_size=star_size)

    ax.set_xticks(range(ncol))
    ax.set_xticklabels(positions, fontsize=7, color=INK)
    ax.set_yticks(range(nrow))
    # The residue letter is the key to its row and was set smaller than the
    # position numbers carry visually; it is one step up from that here.
    ax.set_yticklabels(order, fontsize=12.0 if kind == "residue" else 10.5,
                       fontweight="bold" if kind == "residue" else "normal",
                       color=INK)
    ax.tick_params(length=0, pad=3)
    ax.set_xlabel("Position in the 24-mer", fontsize=8.5, labelpad=4, color=INK)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xlim(-0.5, ncol - 0.5)

    cbar = fig.colorbar(im, ax=ax, fraction=0.022, pad=0.015)
    cbar.ax.tick_params(labelsize=6.5, length=2, width=0.6, colors=INK)
    cbar.outline.set_linewidth(0.5)
    cbar.outline.set_edgecolor(pg.RULE)
    cbar.set_ticks([0, 1.30103, 2, 3, 4, vmax])
    cbar.set_ticklabels(["n.s.", "0.05", "0.01", "0.001", "1e-4", "1e-5"])
    cbar.set_label("$p$ value  (colour = $-$log$_{10}$ $p$)" if signed
                   else "$p$ value  (enriched cells only)",
                   fontsize=7, labelpad=6, color=INK)
    cbar.ax.axhline(1.30103, color=INK, linewidth=0.7, linestyle=(0, (2, 1.6)))

    full = cells[(cells.kind == kind) & cells.defined]
    k = full[full.position.isin(positions)]
    n_test, n_full = len(k), len(full)
    n_hit = int((k[pcol] < 0.05).sum())
    n_fdr = int((k[qcol] < 0.05).sum())
    # the FDR survivors are named in the caption rather than marked in the
    # matrix, so the cells keep to one colour and one number
    fdr_names = ", ".join(
        f"{r.label} at {r.position}"
        for _, r in k[k[qcol] < 0.05].sort_values(qcol).iterrows())

    ax.set_title(title, fontweight="bold", fontsize=9.5, pad=16, color=INK)
    ax.text(0.5, 1.012, subtitle, transform=ax.transAxes, ha="center",
            va="bottom", fontsize=6.4, color=MUTED, style="italic")
    foot = [
        f"n = {N_SUB} P/G-D/E/T substrates  vs  n = {N_CTRL} non-substrate "
        f"controls carrying the same motif."
        + ("   " + CLASS_NOTE if kind == "residue" else ""),

        "Cell text is the fold change (substrates / controls) of every cell at "
        "$p$ < 0.05, under the star bin of its uncorrected $p$:  "
        "* $p$ < 0.05, ** $p$ < 0.01, *** $p$ < 0.001.   "
        + ("Grey = residue absent from both groups, no test possible."
           if kind == "residue" else ""),

        "No depleted cell reaches $p$ < 0.05 anywhere in this dataset, so "
        "every coloured cell above is an enrichment.   "
        f"{n_hit} of the {n_test} testable cells "
        + ("shown here " if window else "")
        + "reach $p$ < 0.05 uncorrected, against "
        f"{0.05 * n_test:.0f} expected by chance alone — "
        + (f"and only {fdr_names} survive"
            f"{'s' if n_fdr == 1 else ''} BH-FDR ($q$ < 0.05) — the rest "
            f"is not separable from that expectation."
            if n_fdr else
            "and none survives FDR correction, so no single cell here is "
            "more than a multiple-testing artefact."),
    ]
    if window:
        foot.append(
            f"This panel is the N-terminal window only — positions "
            f"{positions[0]}–{positions[-1]}, out to the {positions[-1]}th "
            f"amino acid of the 24-mer. The crop is a display crop: every "
            f"cell's $p$ is unchanged, and $q$ is still BH-FDR across all "
            f"{n_full} cells of the full {POS[0]}–{POS[-1]} matrix rather "
            f"than recomputed inside the window — recomputing it would be "
            f"choosing the region after seeing the data.")
    ax.text(0.5, -0.085 if kind == "residue" else -0.26,
            "\n".join(wrap(f, figsize[0]) for f in foot),
            transform=ax.transAxes, ha="center",
            va="top", fontsize=5.6, color=MUTED, linespacing=1.5)

    plt.tight_layout(pad=0.4)
    save(fig, stem)


# ---------------------------------------------------------------- figures
print("building significance heatmaps\n")

# Every figure is drawn twice: over the full 4-24 matrix, and over the
# N-terminal window out to the 12th amino acid. Same data, same colour scale,
# same encoding - the second is the first with its right-hand columns cropped.
WINDOWS = [(POS, "", ""),
           (pg.POS_N12, "_pos4_12", " (positions 4–12)")]

for positions, suffix, wt in WINDOWS:
    heatmap("residue", "p_chi2", "q_chi2", AA, True,
            "Position × residue enrichment, coloured by significance" + wt,
            "colour = $-$log$_{10}$ $p$, chi-square as reported in the "
            "workbook · cell text = fold change · "
            "deeper red = more significant",
            "Figure16_significance_heatmap_residues" + suffix, 8.0, positions)

    heatmap("residue", "p_fisher", "q_fisher", AA, True,
            "Position × residue enrichment, coloured by significance "
            "(Fisher exact)" + wt,
            "colour = $-$log$_{10}$ $p$, Fisher exact recomputed from the "
            "same 2×2 tables — valid at the low cell counts here",
            "Figure16b_significance_heatmap_residues_fisher" + suffix,
            8.0, positions)

    heatmap("residue", "p_fisher", "q_fisher", AA, False,
            "Position × residue enrichment, coloured by significance "
            "(enriched cells only)" + wt,
            "same layout as the original Colab heatmap, but the colour scale "
            "is $-$log$_{10}$ $p$ (Fisher exact) and depleted cells are held "
            "at the floor",
            "Figure16c_significance_heatmap_residues_enriched_only" + suffix,
            8.0, positions)

    heatmap("category", "p_chi2", "q_chi2", CAT_ORDER, True,
            "Chemical class × position, coloured by significance" + wt,
            "colour = $-$log$_{10}$ $p$, chi-square as reported in the "
            "workbook · Basic K H R · Acidic D E · "
            "Aromatic F Y W · Aliphatic G P A M · "
            "Hydrophobic V L I · Polar S T C N Q",
            "Figure17_significance_heatmap_categories" + suffix,
            3.0, positions, annot_size=6.6)

    heatmap("category", "p_fisher", "q_fisher", CAT_ORDER, True,
            "Chemical class × position, coloured by significance "
            "(Fisher exact)" + wt,
            "colour = $-$log$_{10}$ $p$, Fisher exact recomputed from the "
            "same 2×2 tables",
            "Figure17b_significance_heatmap_categories_fisher" + suffix,
            3.0, positions, annot_size=6.6)

# ---------------------------------------------------------------- readout
res = cells[(cells.kind == "residue") & cells.defined].copy()
res["stars_chi2"] = res.p_chi2.map(stars)

print("\n" + "=" * 78)
print("fold change and significance do NOT rank the same cells")
print("=" * 78)
enr = res[res.direction > 0].copy()
cols = ["position", "label", "n_sub", "n_ctrl", "fold_change", "p_chi2",
        "p_fisher", "q_fisher"]
fmt = {"fold_change": "{:.2f}".format, "p_chi2": "{:.2e}".format,
       "p_fisher": "{:.2e}".format, "q_fisher": "{:.3f}".format}
print("\ntop 6 by FOLD CHANGE (what the old heatmap makes loudest):")
print(enr.nlargest(6, "fold_change")[cols].to_string(index=False, formatters=fmt))
print("\ntop 6 by SIGNIFICANCE (what this heatmap makes loudest):")
print(enr.nsmallest(6, "p_chi2")[cols].to_string(index=False, formatters=fmt))

print("\n" + "=" * 78)
print("what survives multiple-testing correction")
print("=" * 78)
for kind, n in (("residue", 420), ("category", 126)):
    k = cells[(cells.kind == kind) & cells.defined]
    print(f"\n{kind}s ({len(k)} testable cells of {n}):")
    for col, lab in (("p_chi2", "chi-square"), ("p_fisher", "Fisher exact")):
        q = col.replace("p_", "q_")
        print(f"   {lab:<13} p < 0.05: {int((k[col] < 0.05).sum()):>3} cells"
              f"    q < 0.05 (BH): {int((k[q] < 0.05).sum()):>3} cells"
              f"    expected by chance at p<0.05: {0.05*len(k):.0f}")
    surv = k[k.q_fisher < 0.05].sort_values("q_fisher")
    if len(surv):
        print("   FDR survivors:")
        print(surv[cols + ["q_chi2"]].to_string(index=False, formatters=fmt))
    else:
        print("   FDR survivors: none")

print("\nwrote data/significance_heatmap_cells.csv "
      f"({len(cells)} cells: {int((cells.kind=='residue').sum())} residue + "
      f"{int((cells.kind=='category').sum())} category)")
