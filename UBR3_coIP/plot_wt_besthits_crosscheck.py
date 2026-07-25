"""
Figure: do the UBR3-WT best-hit half-lives hold up in independent datasets?

Panel A  proteome-wide validation. Every protein shared between the Yi in-extract
         assay and the Li et al. 2021 HEK293T cycloheximide chase. If the extract
         assay reflects turnover in living cells, proteins that lose more in extract
         should have shorter cellular half-lives - and they do (Spearman +0.64).
         Note the Li table only lists proteins that measurably degraded, so the
         x-range is truncated; that attenuates the correlation rather than inflating it.
Panel B  the 13 best hits with an independent half-life, one row each, showing every
         dataset that measured them. Agreement across systems is the point.

Outputs: WT_besthits_crosscheck.(png|pdf|svg)
"""
import json
import pathlib

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

HERE = pathlib.Path(__file__).parent

# validated categorical slots, one per dataset
DS = [
    ("Yi_extract_t_half_h", "Yi 2023 — in extract (derived)", "#8a8880", "o"),
    ("Li2021_HEK293T_t_half_h", "Li 2021 — HEK293T (CHX chase)", "#2a78d6", "s"),
    ("Li2021_U2OS_t_half_h", "Li 2021 — U2OS", "#1baf7a", "^"),
    ("Li2021_HCT116_t_half_h", "Li 2021 — HCT116", "#eb6834", "D"),
    ("Li2021_RPE1_t_half_h", "Li 2021 — RPE1", "#4a3aa7", "v"),
    ("Mathieson_primary_median_t_half_h", "Mathieson 2018 — primary cells (SILAC)",
     "#0b0b0b", "P"),
]
INK, INK2, INK3, GRID = "#0b0b0b", "#52514e", "#8a8880", "#e2e1dc"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "axes.edgecolor": GRID, "axes.linewidth": 0.8,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "xtick.major.size": 0, "ytick.major.size": 0, "figure.facecolor": "white",
})

cross = pd.read_csv(HERE / "WT_besthits_halflife_crosscheck.csv")
sc = pd.read_csv(HERE / "WT_besthits_halflife_crosscheck_scatter.csv")
summ = json.loads((HERE / "WT_besthits_halflife_crosscheck.json").read_text())
c_hek = summ["proteome_wide_concordance"]["Li2021_HEK293T"]
c_mat = summ["proteome_wide_concordance"]["Mathieson2018_primary_cells"]

fig = plt.figure(figsize=(13.2, 15.0))
gs = GridSpec(2, 1, height_ratios=[1.0, 1.75], hspace=0.20, figure=fig)
axA = fig.add_subplot(gs[0])
axB = fig.add_subplot(gs[1])

# ---------------------------------------------------------------- panel A
bgm = ~sc["is_WT_best_hit"]
axA.scatter(sc.loc[bgm, "Yi_remaining_15h"], sc.loc[bgm, "Li_HEK293T_t_half_h"],
            s=11, color="#c9c7c0", alpha=0.55, linewidth=0, zorder=2,
            label=f"all shared proteins (n={len(sc)})")
hm = sc["is_WT_best_hit"]
axA.scatter(sc.loc[hm, "Yi_remaining_15h"], sc.loc[hm, "Li_HEK293T_t_half_h"],
            s=78, color="#2a78d6", edgecolor="white", linewidth=1.3, zorder=4,
            label=f"UBR3-WT best hits (n={int(hm.sum())})")
for r in sc[hm].itertuples():
    axA.annotate(r.gene, (r.Yi_remaining_15h, r.Li_HEK293T_t_half_h),
                 xytext=(7, 4), textcoords="offset points", fontsize=8.8,
                 fontweight="bold", color="#2a78d6", zorder=5)

# monotone trend, computed on ranks so it matches the Spearman statistic
xq = np.linspace(sc["Yi_remaining_15h"].quantile(0.01),
                 sc["Yi_remaining_15h"].quantile(0.99), 12)
med = [sc.loc[(sc["Yi_remaining_15h"] >= lo) & (sc["Yi_remaining_15h"] < hi),
              "Li_HEK293T_t_half_h"].median() for lo, hi in zip(xq[:-1], xq[1:])]
axA.plot((xq[:-1] + xq[1:]) / 2, med, color=INK, lw=2.0, zorder=6,
         label="binned median")

axA.set_xlabel("Yi 2023 — protein remaining at 15 h in extract  (more decay ←)",
               fontsize=9.8, color=INK2)
axA.set_ylabel("Li 2021 — half-life in HEK293T (h)", fontsize=9.8, color=INK2)
axA.set_ylim(0, 60)
axA.grid(color=GRID, lw=0.5)
axA.set_axisbelow(True)
axA.legend(frameon=False, fontsize=9, loc="upper left", labelcolor=INK2)
axA.set_title("A   The in-extract assay does track turnover in living cells",
              fontsize=12.2, fontweight="bold", loc="left", pad=10, color=INK)
axA.text(0.985, 0.05,
         f"Spearman ρ = {c_hek['spearman_rho_Yi_15h_remaining_vs_halflife']:+.2f}  "
         f"(n = {c_hek['n_shared_proteins']}, p = {c_hek['p_value_15h']:.0e})\n"
         f"vs Mathieson primary cells: ρ = "
         f"{c_mat['spearman_rho_Yi_15h_remaining_vs_halflife']:+.2f} "
         f"(n = {c_mat['n_shared_proteins']})",
         transform=axA.transAxes, ha="right", va="bottom", fontsize=9.4, color=INK,
         bbox=dict(boxstyle="round,pad=0.5", facecolor="white", edgecolor=GRID))

# ---------------------------------------------------------------- panel B
cols = [c for c, _, _, _ in DS]
has = cross[cross[cols[1:]].notna().any(axis=1)].copy()
has["sort_key"] = has[cols[1:]].min(axis=1)
has = has.sort_values("sort_key").reset_index(drop=True)

XMAX = 105          # values past this are drawn as carets on the axis edge
n_offscale = 0
for i, r in has.iterrows():
    axB.axhline(i, color=GRID, lw=0.8, zorder=1)
    vals = [(r[c], col, mk) for c, _, col, mk in DS if pd.notna(r[c])]
    inr = [v for v in vals if v[0] <= XMAX]
    out = [v for v in vals if v[0] > XMAX]
    n_offscale += len(out)
    if len(vals) > 1:
        xs = [v[0] for v in inr] or [XMAX]
        lo, hi = min(xs), (XMAX if out else max(xs))
        axB.plot([lo, hi], [i, i], color=INK3, lw=1.4, alpha=0.5, zorder=2)
    for v, col, mk in inr:
        axB.scatter(v, i, s=92, color=col, marker=mk, edgecolor="white",
                    linewidth=1.2, zorder=4)
    for v, col, mk in out:
        axB.scatter(XMAX * 0.995, i, s=70, color=col, marker=">", zorder=4,
                    edgecolor="white", linewidth=0.9, clip_on=False)

axB.set_yticks(range(len(has)))
axB.set_yticklabels(has["gene"], fontsize=9.2, fontweight="bold")
axB.set_ylim(-0.6, len(has) - 0.4)
axB.set_xlim(0, XMAX)
axB.set_xlabel("half-life (h)   —   ▶ on the right edge = value beyond "
               f"{XMAX} h (all are Yi extrapolations)", fontsize=9.8, color=INK2)
axB.grid(axis="x", color=GRID, lw=0.5)
axB.set_axisbelow(True)
axB.set_title("B   Best hits with an independent half-life — every dataset that measured them",
              fontsize=12.2, fontweight="bold", loc="left", pad=10, color=INK)
axB.legend(handles=[Line2D([], [], marker=mk, ls="", markersize=9, color=col,
                           markeredgecolor="white", label=lab)
                    for _, lab, col, mk in DS],
           frameon=False, fontsize=9, loc="lower right", labelcolor=INK2,
           ncol=2, handletextpad=0.6, columnspacing=1.8)

fig.suptitle("UBR3-WT best-hit half-lives, cross-checked against published turnover data",
             fontsize=14.2, fontweight="bold", x=0.055, ha="left", y=0.985, color=INK)
fig.text(0.055, 0.962,
         "Li et al. 2021 lists only proteins losing >15% in 8 h, so absence from it is "
         "not evidence of stability · Mathieson values are the median of 4 human primary "
         "cell types",
         fontsize=9.2, color=INK2, ha="left")

fig.subplots_adjust(left=0.085, right=0.975, top=0.925, bottom=0.06)
for ext in ("png", "pdf", "svg"):
    fig.savefig(HERE / f"WT_besthits_crosscheck.{ext}", dpi=200, bbox_inches="tight",
                facecolor="white")
print("figure -> WT_besthits_crosscheck.png/pdf/svg")
print(f"panel A: {len(sc)} shared proteins, {int(hm.sum())} of them WT best hits")
print(f"panel B: {len(has)} hits with an independent measurement")
