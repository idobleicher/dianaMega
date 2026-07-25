"""
Figures for the mutual-pathway analysis of the UBR3-WT best hits.

Panel A  membership matrix: which best hit supports which enriched term. Terms are
         grouped into the modules found by wt_mutual_pathways.py, so a module reads
         as a filled block = hits that share pathways with each other. The right
         strip is each term's significance; open circles mark terms that disappear
         when the 7 histone genes are dropped (histone-carried annotation, not a
         genuine multi-protein module).
Panel B  hit x hit co-membership: how many specific terms each pair of hits shares.
         Block structure on the diagonal = the mutual modules; off-diagonal signal =
         hits that bridge two modules.

Outputs: WT_mutual_pathways.(png|pdf|svg)
"""
import pathlib

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

HERE = pathlib.Path(__file__).parent

# validated categorical slots (CVD-safe on all pairs, light surface)
MODULE_COLOR = {1: "#2a78d6", 2: "#eb6834", 3: "#1baf7a", 4: "#4a3aa7"}
MODULE_NAME = {
    1: "mRNA silencing & RNP granules",
    2: "ATP-dependent chromatin remodeling",
    3: "Nucleosome / histone core",
    4: "BAF–LAP2β nuclear-envelope complex",
}
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"
GRID = "#e2e1dc"
SEQ = LinearSegmentedColormap.from_list(
    "seq_blue", ["#eef4fc", "#c3d9f4", "#8fb8e8", "#5892db", "#2a78d6", "#17457d"])

plt.rcParams.update({
    "font.family": "DejaVu Sans", "axes.edgecolor": GRID, "axes.linewidth": 0.8,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "xtick.major.size": 0, "ytick.major.size": 0, "figure.facecolor": "white",
})

terms = pd.read_csv(HERE / "WT_mutual_terms.csv")
genes_df = pd.read_csv(HERE / "WT_mutual_genes.csv")
spec = terms[terms["specific"] & (terms["module"] > 0)].copy()
spec["genes_list"] = spec["genes"].fillna("").apply(lambda s: s.split(";") if s else [])
spec = spec.sort_values(["module", "p_value"]).reset_index(drop=True)

# ---- each hit's dominant module (module holding most of its specific terms) --
hit_mod_counts = {}
for r in spec.itertuples():
    for g in r.genes_list:
        hit_mod_counts.setdefault(g, {}).setdefault(r.module, 0)
        hit_mod_counts[g][r.module] += 1
dominant = {g: max(c.items(), key=lambda kv: (kv[1], -kv[0]))[0]
            for g, c in hit_mod_counts.items()}
n_spec_of = dict(zip(genes_df["gene"], genes_df["n_specific_terms"]))

hits = sorted(dominant, key=lambda g: (dominant[g], -n_spec_of[g], g))
xi = {g: i for i, g in enumerate(hits)}
gene_breaks = [i for i in range(1, len(hits)) if dominant[hits[i]] != dominant[hits[i - 1]]]

# ---- pairwise co-membership over the specific terms -------------------------
n = len(hits)
C = np.zeros((n, n))
for r in spec.itertuples():
    mem = [g for g in r.genes_list if g in xi]
    for i in range(len(mem)):
        for j in range(len(mem)):
            if i != j:
                C[xi[mem[i]], xi[mem[j]]] += 1
np.fill_diagonal(C, np.nan)

# ============================================================ figure
fig = plt.figure(figsize=(13.6, 17.4))
gs = GridSpec(2, 3, height_ratios=[1.18, 1.0], width_ratios=[1, 0.030, 0.145],
              hspace=0.46, wspace=0.015, figure=fig)
axA = fig.add_subplot(gs[0, 0])
axM = fig.add_subplot(gs[0, 1], sharey=axA)
axP = fig.add_subplot(gs[0, 2], sharey=axA)
axB = fig.add_subplot(gs[1, 0])

# ---------------------------------------------------------------- panel A
for r in spec.itertuples():
    xs = [xi[g] for g in r.genes_list if g in xi]
    axA.scatter(xs, [r.Index] * len(xs), s=84, marker="s",
                color=MODULE_COLOR[r.module], edgecolor="white", linewidth=1.1, zorder=3)

bounds = {m: (d.index.min(), d.index.max()) for m, d in spec.groupby("module")}
for m, (lo, hi) in bounds.items():
    axA.axhspan(lo - 0.5, hi + 0.5, color=MODULE_COLOR[m], alpha=0.055, zorder=0)
    if hi + 0.5 < len(spec) - 0.5:
        axA.axhline(hi + 0.5, color=INK3, lw=0.8, ls=(0, (4, 3)), zorder=2)
for b in gene_breaks:
    axA.axvline(b - 0.5, color=INK3, lw=0.8, ls=(0, (4, 3)), zorder=2)

axA.set_yticks(range(len(spec)))
axA.set_yticklabels(
    [f"{r.term_name[:54]}{'…' if len(r.term_name) > 54 else ''}  ({r.source.replace('GO:', '')})"
     for r in spec.itertuples()], fontsize=8.3)
for tick, r in zip(axA.get_yticklabels(), spec.itertuples()):
    if not r.survives_without_histones:
        tick.set_color(INK3)
        tick.set_style("italic")
axA.set_xticks(range(n))
axA.set_xticklabels(hits, rotation=90, fontsize=8.7)
for tick, g in zip(axA.get_xticklabels(), hits):
    tick.set_color(MODULE_COLOR[dominant[g]])
    tick.set_fontweight("bold")
axA.set_xlim(-0.6, n - 0.4)
axA.set_ylim(len(spec) - 0.5, -0.5)
axA.grid(color=GRID, lw=0.5, zorder=1)
axA.set_axisbelow(True)
axA.set_title("A   Which best hit supports which enriched term",
              fontsize=12.5, fontweight="bold", loc="left", pad=10, color=INK)
axA.set_xlabel("UBR3-WT best hit   (colour = module holding most of its terms)",
               fontsize=9.3, color=INK2, labelpad=7)

# module strip between matrix and significance
axM.set_xlim(0, 1)
for m, (lo, hi) in bounds.items():
    axM.add_patch(Rectangle((0.08, lo - 0.5), 0.84, hi - lo + 1,
                            color=MODULE_COLOR[m], zorder=3))
    axM.text(0.5, (lo + hi) / 2, f"M{m}", ha="center", va="center", rotation=90,
             fontsize=9.5, fontweight="bold", color="white", zorder=4)
axM.set_axis_off()

# significance strip
axP.barh(range(len(spec)), spec["neg_log10_p"], height=0.55,
         color=[MODULE_COLOR[m] for m in spec["module"]], zorder=3)
axP.axvline(-np.log10(0.05), color=INK3, lw=0.9, ls=(0, (3, 2)), zorder=2)
xmax = spec["neg_log10_p"].max()
for r in spec.itertuples():
    axP.scatter(xmax * 1.16, r.Index, s=32, marker="o", zorder=4,
                facecolor=MODULE_COLOR[r.module] if r.survives_without_histones else "white",
                edgecolor=MODULE_COLOR[r.module], linewidth=1.2)
axP.set_xlim(0, xmax * 1.30)
axP.set_xticks([0, 2, 4])
plt.setp(axP.get_yticklabels(), visible=False)
axP.tick_params(axis="y", left=False)
axP.set_xlabel("$-\\log_{10}p$", fontsize=9, color=INK2)
axP.grid(axis="x", color=GRID, lw=0.5, zorder=1)
axP.set_axisbelow(True)

axA.legend(handles=[
    Line2D([], [], marker="o", ls="", markersize=7, color=INK2,
           label="term survives removal of the 7 histone hits"),
    Line2D([], [], marker="o", ls="", markersize=7, markerfacecolor="white",
           markeredgecolor=INK2, color="none",
           label="term lost without histones — histone-carried annotation"),
], loc="upper left", bbox_to_anchor=(0.0, -0.235), frameon=False, fontsize=8.7,
    ncol=2, labelcolor=INK2, handletextpad=0.6, columnspacing=1.8)

# ---------------------------------------------------------------- panel B
im = axB.imshow(C, cmap=SEQ, vmin=0, vmax=np.nanmax(C), zorder=2)
axB.set_facecolor("#f7f7f5")
for b in gene_breaks:
    axB.axvline(b - 0.5, color="white", lw=2.2, zorder=3)
    axB.axhline(b - 0.5, color="white", lw=2.2, zorder=3)
    axB.axvline(b - 0.5, color=INK3, lw=0.8, ls=(0, (4, 3)), zorder=4)
    axB.axhline(b - 0.5, color=INK3, lw=0.8, ls=(0, (4, 3)), zorder=4)

axB.set_xticks(range(n))
axB.set_xticklabels(hits, rotation=90, fontsize=8.7)
axB.set_yticks(range(n))
axB.set_yticklabels(hits, fontsize=8.7)
for ticks in (axB.get_xticklabels(), axB.get_yticklabels()):
    for tick, g in zip(ticks, hits):
        tick.set_color(MODULE_COLOR[dominant[g]])
        tick.set_fontweight("bold")
axB.set_xlim(-0.5, n - 0.5)
axB.set_ylim(n - 0.5, -0.5)
axB.set_title("B   Specific terms shared by each pair of best hits",
              fontsize=12.5, fontweight="bold", loc="left", pad=12, color=INK)

cb = fig.colorbar(im, ax=axB, fraction=0.032, pad=0.02, shrink=0.62)
cb.set_label("shared specific terms", fontsize=9, color=INK2)
cb.outline.set_edgecolor(GRID)
cb.ax.tick_params(labelsize=8.5, color=INK2)

axB.legend(handles=[
    Line2D([], [], marker="s", ls="", markersize=10, color=MODULE_COLOR[m],
           label=f"M{m}  {MODULE_NAME[m]}") for m in sorted(MODULE_COLOR)],
    loc="upper left", bbox_to_anchor=(0.0, -0.185), frameon=False,
    fontsize=10, ncol=2, labelcolor=INK, handletextpad=0.7, columnspacing=2.6)

fig.suptitle("UBR3-WT co-IP best hits converge on four mutual pathway modules",
             fontsize=14.5, fontweight="bold", x=0.055, ha="left", y=0.972, color=INK)
fig.text(0.055, 0.955,
         "50 significant WT interactors (q<0.05 vs 293T) · g:Profiler g:GOSt, g:SCS-corrected · "
         "29 specific terms (size ≤ 500) shown; broad terms omitted · "
         "25 of the 50 hits fall in no specific term and are not plotted",
         fontsize=9.2, color=INK2, ha="left")

fig.subplots_adjust(left=0.30, right=0.83, top=0.928, bottom=0.115)
for ext in ("png", "pdf", "svg"):
    fig.savefig(HERE / f"WT_mutual_pathways.{ext}", dpi=200, bbox_inches="tight",
                facecolor="white")
print("figure -> WT_mutual_pathways.png/pdf/svg")
print(f"panel A: {len(spec)} specific terms x {n} hits")
print(f"panel B: {n} x {n} co-membership, max shared = {int(np.nanmax(C))}")
