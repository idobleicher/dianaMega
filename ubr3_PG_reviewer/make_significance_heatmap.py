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
BH-FDR is applied across all 420 residue cells (126 category cells) and the
survivors are outlined in black.

Nothing is written inside a cell except the star bin for cells at p < 0.05:
the colour is the p value, so a number in the cell would only repeat it, and
fold change is deliberately absent -- this figure is about significance alone.
Every cell's exact p, q, both counts, both percentages and the fold change are
in data/significance_heatmap_cells.csv for anything that needs quoting in text
(p is not a substitute for effect size: it folds the fold change together with
how many peptides the cell rests on).

Type is Helvetica throughout (Arial where Helvetica is absent; see HELV), and
no left-hand class colour panel is drawn -- the residue rows are still ordered
by chemical class, which the caption states in words.

Outputs (figures/, 600 dpi PNG + vector PDF + SVG):
  Figure16_significance_heatmap_residues        -log10 chi2 p, warm scale
  Figure16b_significance_heatmap_residues_fisher  -log10 Fisher p
  Figure16c_significance_heatmap_residues_enriched_only  depleted cells forced to
                                                the floor; drop-in for the Colab
                                                figure's look
  Figure17_significance_heatmap_categories      the six chemical classes
  data/significance_heatmap_cells.csv           every cell, both tests, q values
"""
import os
import re

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
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


def matrix(kind, col, order):
    return pg.matrix(col, kind=kind, order=order, cells=cells)


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
})

# the first member of the stack that is actually installed -- mathtext resolves
# a family by name, so it has to be given a real one or it silently drops to
# DejaVu and the maths comes out in a different face from the labels
_INSTALLED = {f.name for f in mpl.font_manager.fontManager.ttflist}
FONT = next((f for f in HELV if f in _INSTALLED), "DejaVu Sans")
mpl.rcParams.update({"mathtext.rm": FONT, "mathtext.it": f"{FONT}:italic",
                     "mathtext.bf": f"{FONT}:bold"})


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


# One warm ramp, cream -> yellow -> orange -> red -> dark red, in the house
# palette (#E8A33D / #C0392B / #7B241C are the Aliphatic / Basic / Acidic
# colours used by the logos), so the heatmaps and the logos read as one set.
#
# Why not a diverging scale with a cool arm for depletion: NO depleted cell
# reaches p < 0.05 anywhere in this dataset -- not in either test, not for
# residues, not for classes. The most significant depleted cell is E at
# position 8 at p = 0.11. A cool arm would therefore colour nothing but
# sub-threshold noise while competing for attention with the real signal.
# Direction is still recoverable from the annotation: the fold change printed
# in a cell is > 1 for enrichment and < 1 for depletion.
CMAP_MAIN = LinearSegmentedColormap.from_list("sig_warm", [
    "#FFFDF8", "#FDF3D4", "#FAE3A0", "#F5C463", "#E8A33D",
    "#DD7B2F", "#C0392B", "#96271F", "#7B241C"])
CMAP_ENRICH = CMAP_MAIN
GREY = "#E9E9E7"          # residue absent from both groups -- no test possible


CLASS_NOTE = ("Rows are ordered by chemical class — Basic K H R · Acidic D E · "
              "Aromatic F Y W · Aliphatic G P A M · Hydrophobic V L I · "
              "Polar S T C N Q.   ")


def heatmap(kind, pcol, qcol, order, signed, title, subtitle, stem,
            figsize, vmax=5.0, annot_size=8.5):
    P = matrix(kind, pcol, order).values.astype(float)
    Q = matrix(kind, qcol, order).values.astype(float)
    D = matrix(kind, "direction", order).values.astype(float)
    OK = matrix(kind, "defined", order).values.astype(float) == 1

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

    # white gridlines, as in the Colab figure
    for x in np.arange(-0.5, ncol, 1):
        ax.axvline(x, color="white", linewidth=0.7)
    for y in np.arange(-0.5, nrow, 1):
        ax.axhline(y, color="white", linewidth=0.7)

    # No numbers inside the cells: colour carries the p value, and the only
    # marks are the conventional star bins on cells at p < 0.05 plus a black
    # outline on the FDR survivors. Every cell's exact p, q, counts,
    # percentages and fold change are in data/significance_heatmap_cells.csv.
    for r in range(nrow):
        for c_ in range(ncol):
            p_, q = P[r, c_], Q[r, c_]
            if not (np.isfinite(p_) and p_ < 0.05):
                continue
            shade = C[r, c_]
            dark = np.ma.is_masked(shade) or abs(shade) > 0.62 * vmax
            # asterisks hang from the cap line, so nudge them down to sit
            # optically centred in the cell
            ax.text(c_, r + 0.16, stars(p_), ha="center", va="center",
                    fontsize=annot_size, fontweight="bold",
                    color="white" if dark else "#1A1A1A")
            if np.isfinite(q) and q < 0.05:
                ax.add_patch(Rectangle((c_ - 0.5, r - 0.5), 1, 1, fill=False,
                                       edgecolor="black", linewidth=1.4, zorder=5))

    ax.set_xticks(range(ncol))
    ax.set_xticklabels(POS, fontsize=7)
    ax.set_yticks(range(nrow))
    ax.set_yticklabels(order, fontsize=9.5 if kind == "residue" else 8.5,
                       fontweight="bold" if kind == "residue" else "normal")
    ax.tick_params(length=0, pad=3)
    ax.set_xlabel("Position in the 24-mer", fontsize=8.5, labelpad=4)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xlim(-0.5, ncol - 0.5)

    cbar = fig.colorbar(im, ax=ax, fraction=0.022, pad=0.015)
    cbar.ax.tick_params(labelsize=6.5, length=2, width=0.6)
    cbar.outline.set_linewidth(0.6)
    cbar.set_ticks([0, 1.30103, 2, 3, 4, vmax])
    cbar.set_ticklabels(["n.s.", "0.05", "0.01", "0.001", "1e-4", "1e-5"])
    cbar.set_label("$p$ value  (colour = $-$log$_{10}$ $p$)" if signed
                   else "$p$ value  (enriched cells only)",
                   fontsize=7, labelpad=6)
    cbar.ax.axhline(1.30103, color="#333", linewidth=0.7, linestyle=(0, (2, 1.6)))

    k = cells[(cells.kind == kind) & cells.defined]
    n_test = len(k)
    n_hit = int((k[pcol] < 0.05).sum())
    n_fdr = int((k[qcol] < 0.05).sum())

    ax.set_title(title, fontweight="bold", fontsize=9.5, pad=16)
    ax.text(0.5, 1.012, subtitle, transform=ax.transAxes, ha="center",
            va="bottom", fontsize=6.4, color="#666", style="italic")
    foot = [
        f"n = {N_SUB} P/G-D/E/T substrates  vs  n = {N_CTRL} non-substrate "
        f"controls carrying the same motif."
        + ("   " + CLASS_NOTE if kind == "residue" else ""),

        "Stars are the uncorrected $p$ value of the cell:  * $p$ < 0.05, "
        "** $p$ < 0.01, *** $p$ < 0.001.   "
        f"Black outline = survives BH-FDR across all {n_test} cells "
        "($q$ < 0.05).   "
        + ("Grey = residue absent from both groups, no test possible."
           if kind == "residue" else ""),

        "No depleted cell reaches $p$ < 0.05 anywhere in this dataset, so "
        "every coloured cell above is an enrichment.   "
        f"{n_hit} of {n_test} testable cells reach $p$ < 0.05 uncorrected, "
        f"against {0.05 * n_test:.0f} expected by chance alone — "
        + ("only the outlined cells are more than a multiple-testing "
           "artefact." if n_fdr else
           "and none survives FDR correction, so no single cell here is "
           "more than a multiple-testing artefact."),
    ]
    ax.text(0.5, -0.085 if kind == "residue" else -0.26,
            "\n".join(foot), transform=ax.transAxes, ha="center",
            va="top", fontsize=5.6, color="#666", linespacing=1.5)

    plt.tight_layout(pad=0.4)
    save(fig, stem)


# ---------------------------------------------------------------- figures
print("building significance heatmaps\n")

heatmap("residue", "p_chi2", "q_chi2", AA, True,
        "Position × residue enrichment, coloured by significance",
        "colour = $-$log$_{10}$ $p$, chi-square as reported in the "
        "workbook · deeper red = more significant",
        "Figure16_significance_heatmap_residues", (13.4, 8.0))

heatmap("residue", "p_fisher", "q_fisher", AA, True,
        "Position × residue enrichment, coloured by significance "
        "(Fisher exact)",
        "colour = $-$log$_{10}$ $p$, Fisher exact recomputed from the "
        "same 2×2 tables — valid at the low cell counts here",
        "Figure16b_significance_heatmap_residues_fisher", (13.4, 8.0))

heatmap("residue", "p_fisher", "q_fisher", AA, False,
        "Position × residue enrichment, coloured by significance "
        "(enriched cells only)",
        "same layout as the original Colab heatmap, but the colour scale is "
        "$-$log$_{10}$ $p$ (Fisher exact) and depleted cells are held "
        "at the floor",
        "Figure16c_significance_heatmap_residues_enriched_only", (13.4, 8.0))

heatmap("category", "p_chi2", "q_chi2", CAT_ORDER, True,
        "Chemical class × position, coloured by significance",
        "colour = $-$log$_{10}$ $p$, chi-square as reported in the "
        "workbook · Basic K H R · Acidic D E · "
        "Aromatic F Y W · Aliphatic G P A M · "
        "Hydrophobic V L I · Polar S T C N Q",
        "Figure17_significance_heatmap_categories", (13.4, 3.0),
        annot_size=9.0)

heatmap("category", "p_fisher", "q_fisher", CAT_ORDER, True,
        "Chemical class × position, coloured by significance "
        "(Fisher exact)",
        "colour = $-$log$_{10}$ $p$, Fisher exact recomputed from the "
        "same 2×2 tables",
        "Figure17b_significance_heatmap_categories_fisher", (13.4, 3.0),
        annot_size=9.0)

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
