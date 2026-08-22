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
CAT_COLORS, CAT_OF_AA = pg.CAT_COLORS, pg.CAT_OF_AA
stars = pg.stars

cells = pg.load()
cells.to_csv(os.path.join(DATA_DIR, "significance_heatmap_cells.csv"), index=False)


def matrix(kind, col, order):
    return pg.matrix(col, kind=kind, order=order, cells=cells)


# ---------------------------------------------------------------- style
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8, "axes.linewidth": 0.8,
    "svg.fonttype": "none", "pdf.fonttype": 42,
    "figure.dpi": 600, "savefig.dpi": 600, "savefig.bbox": "tight",
})

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


def heatmap(kind, pcol, qcol, order, signed, title, subtitle, stem,
            figsize, vmax=5.0, annot_size=5.4, group_rows=True):
    P = matrix(kind, pcol, order).values.astype(float)
    Q = matrix(kind, qcol, order).values.astype(float)
    FC = matrix(kind, "fold_change", order).values.astype(float)
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

    # annotate significant cells with their fold change; outline FDR survivors
    for r in range(nrow):
        for c_ in range(ncol):
            p, q, fc = P[r, c_], Q[r, c_], FC[r, c_]
            if not (np.isfinite(p) and p < 0.05):
                continue
            s = stars(p)
            fc_txt = "∞" if np.isinf(fc) else ("--" if not np.isfinite(fc)
                                                    else f"{fc:.2f}")
            shade = C[r, c_]
            dark = np.ma.is_masked(shade) or abs(shade) > 0.62 * vmax
            ax.text(c_, r, f"{fc_txt}{s}", ha="center", va="center",
                    fontsize=annot_size, fontweight="bold",
                    color="white" if dark else "#1A1A1A", linespacing=0.9)
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

    # class strip + separators down the left edge of the residue heatmap
    if group_rows and kind == "residue":
        edges, seen = [], None
        for r, a in enumerate(order):
            cat = CAT_OF_AA[a]
            if cat != seen:
                if seen is not None:
                    ax.axhline(r - 0.5, color="#4D4D4D", linewidth=1.1)
                edges.append([cat, r, r])
                seen = cat
            else:
                edges[-1][2] = r
        for cat, r0, r1 in edges:
            ax.add_patch(Rectangle((-0.5 - 1.55, r0 - 0.5), 0.42, r1 - r0 + 1,
                                   facecolor=CAT_COLORS[cat], clip_on=False,
                                   edgecolor="white", linewidth=0.6, zorder=6))
            ax.text(-0.5 - 1.72, (r0 + r1) / 2, cat, rotation=90, ha="center",
                    va="center", fontsize=6.2, color="#4D4D4D", clip_on=False)
        ax.set_xlim(-0.5, ncol - 0.5)

    cbar = fig.colorbar(im, ax=ax, fraction=0.022, pad=0.015)
    cbar.ax.tick_params(labelsize=6.5, length=2, width=0.6)
    cbar.outline.set_linewidth(0.6)
    cbar.set_ticks([0, 1.30103, 2, 3, 4, vmax])
    cbar.set_ticklabels(["n.s.", "0.05", "0.01", "0.001", "1e-4", "1e-5"])
    cbar.set_label("$p$ value  (colour = $-$log$_{10}\\,p$)"
                   if signed else "$p$ value  (enriched cells only)",
                   fontsize=7, labelpad=6)
    cbar.ax.axhline(1.30103, color="#333", linewidth=0.7, linestyle=(0, (2, 1.6)))

    k = cells[(cells.kind == kind) & cells.defined]
    n_test = len(k)
    n_hit = int((k[pcol] < 0.05).sum())
    n_fdr = int((k[qcol] < 0.05).sum())

    ax.set_title(title, fontweight="bold", fontsize=9.5, pad=16)
    ax.text(0.5, 1.012, subtitle, transform=ax.transAxes, ha="center",
            va="bottom", fontsize=6.4, color="#666", style="italic")
    ax.text(0.5, -0.115 if kind == "residue" else -0.30,
            f"n = {N_SUB} P/G-D/E/T substrates  vs  n = {N_CTRL} non-substrate "
            f"controls carrying the same motif.   Cell text is the fold change "
            f"(substrates / controls) for cells at $p$ < 0.05; "
            f"* $p$ < 0.05, ** $p$ < 0.01, *** $p$ < 0.001.   "
            f"Black outline = survives BH-FDR across all "
            f"{n_test} cells ($q$ < 0.05).   "
            f"Grey = residue absent from both groups, no test possible.   "
            f"No depleted cell reaches $p$ < 0.05 anywhere in this dataset, so "
            f"every coloured cell above is an enrichment.\n"
            f"{n_hit} of {n_test} testable cells reach $p$ < 0.05 uncorrected, "
            f"against {0.05 * n_test:.0f} expected by chance alone — "
            + ("only the outlined cells are more than a multiple-testing "
               "artefact." if n_fdr else
               "and none survives FDR correction, so no single cell here is "
               "more than a multiple-testing artefact."),
            transform=ax.transAxes, ha="center", va="top", fontsize=5.5,
            color="#666")

    plt.tight_layout(pad=0.4)
    for ext in ("png", "pdf", "svg"):
        path = os.path.join(FIG_DIR, f"{stem}.{ext}")
        fig.savefig(path, format=ext)
        print(f"  saved: {os.path.basename(path)}  ({os.path.getsize(path)/1024:.1f} KB)")
    plt.close(fig)


# ---------------------------------------------------------------- figures
print("building significance heatmaps\n")

heatmap("residue", "p_chi2", "q_chi2", AA, True,
        "Position × residue enrichment, coloured by SIGNIFICANCE "
        "(not fold change)",
        "colour = $-$log$_{10}$ $p$, chi-square as reported in the workbook · "
        "deeper red = more significant",
        "Figure16_significance_heatmap_residues", (13.4, 8.0))

heatmap("residue", "p_fisher", "q_fisher", AA, True,
        "Position × residue enrichment, coloured by SIGNIFICANCE "
        "(Fisher exact)",
        "colour = $-$log$_{10}$ $p$, Fisher exact recomputed from the same "
        "2×2 tables — valid at the low cell counts here",
        "Figure16b_significance_heatmap_residues_fisher", (13.4, 8.0))

heatmap("residue", "p_fisher", "q_fisher", AA, False,
        "Position × residue enrichment, coloured by SIGNIFICANCE "
        "(enriched cells only)",
        "same layout and palette as the original fold-change heatmap, but the "
        "colour scale is $-$log$_{10}$ $p$ (Fisher exact)",
        "Figure16c_significance_heatmap_residues_enriched_only", (13.4, 8.0))

heatmap("category", "p_chi2", "q_chi2", CAT_ORDER, True,
        "Chemical class × position, coloured by SIGNIFICANCE "
        "(not fold change)",
        "colour = $-$log$_{10}$ $p$, chi-square as reported in the workbook · "
        "Basic K H R · Acidic D E · Aromatic F Y W · "
        "Aliphatic G P A M · Hydrophobic V L I · Polar S T C N Q",
        "Figure17_significance_heatmap_categories", (13.4, 3.0),
        annot_size=6.0, group_rows=False)

heatmap("category", "p_fisher", "q_fisher", CAT_ORDER, True,
        "Chemical class × position, coloured by SIGNIFICANCE "
        "(Fisher exact)",
        "colour = $-$log$_{10}$ $p$, Fisher exact recomputed from the same "
        "2×2 tables",
        "Figure17b_significance_heatmap_categories_fisher", (13.4, 3.0),
        annot_size=6.0, group_rows=False)

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
