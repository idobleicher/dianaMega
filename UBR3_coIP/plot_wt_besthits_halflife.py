"""
Figure for the stability analysis of the UBR3-WT best hits in the Yi et al. (2023)
degradomics time course.

Panel A  cumulative distribution of protein remaining at 15 h — the 45 quantified
         best hits against the 8098-protein background, with the M2 module drawn
         separately. Answers "are the hits unusually short-lived?".
Panel B  every quantified best hit ranked by protein remaining at 15 h, coloured by
         pathway module, against the background median and IQR.
Panel C  measured decay curves (0/5/10/15 h) for the least stable hits, next to the
         background median trajectory — the data the half-lives are fitted to.

Outputs: WT_besthits_halflife.(png|pdf|svg)
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

MODULE_COLOR = {1: "#2a78d6", 2: "#eb6834", 3: "#1baf7a", 4: "#4a3aa7"}
MODULE_NAME = {1: "M1 mRNA silencing & RNP granules",
               2: "M2 ATP-dependent chromatin remodeling",
               3: "M3 Nucleosome / histone core",
               4: "M4 BAF–LAP2β nuclear envelope"}
NOMOD = "#9a9891"
INK, INK2, INK3 = "#0b0b0b", "#52514e", "#8a8880"
GRID, BGFILL = "#e2e1dc", "#c9c7c0"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "axes.edgecolor": GRID, "axes.linewidth": 0.8,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "xtick.major.size": 0, "ytick.major.size": 0, "figure.facecolor": "white",
})

res = pd.read_csv(HERE / "WT_besthits_halflife.csv")
bg = pd.read_csv(HERE / "WT_besthits_halflife_background.csv")
curves = pd.read_csv(HERE / "WT_besthits_halflife_curves.csv")
summ = json.loads((HERE / "WT_besthits_halflife_summary.json").read_text())

q = res[res["quantified"]].copy()


def first_module(s):
    if pd.isna(s) or str(s).strip() == "":
        return 0
    return int(str(s).split(";")[0])


q["mod"] = q["pathway_modules"].apply(first_module)
color_of = lambda m: MODULE_COLOR.get(m, NOMOD)

TIMES = [0.0, 0.0, 5.0, 5.0, 10.0, 10.0, 15.0, 15.0]
RCOLS = ["WT_UT1_ratio", "WT_UT2_ratio", "WT_UT5h1_ratio", "WT_UT5h2_ratio",
         "WT_UT10h1_ratio", "WT_UT10h2_ratio", "WT_UT15h1_ratio", "WT_UT15h2_ratio"]
TPTS = [0.0, 5.0, 10.0, 15.0]

fig = plt.figure(figsize=(13.4, 14.0))
gs = GridSpec(3, 2, height_ratios=[0.85, 1.55, 0.95], width_ratios=[1, 1],
              hspace=0.42, wspace=0.20, figure=fig)
axA = fig.add_subplot(gs[0, 0])
axL = fig.add_subplot(gs[0, 1])
axB = fig.add_subplot(gs[1, :])
axC = fig.add_subplot(gs[2, :])


def ecdf(v):
    v = np.sort(np.asarray(v, dtype=float))
    return v, np.arange(1, len(v) + 1) / len(v)


# ---------------------------------------------------------------- panel A
x, y = ecdf(bg["frac_remaining_15h"].dropna())
axA.plot(x, y, color=BGFILL, lw=2.6, label=f"all quantified proteins (n={len(x)})", zorder=2)
x, y = ecdf(q["frac_remaining_15h"])
axA.plot(x, y, color=INK, lw=2.0, label=f"UBR3-WT best hits (n={len(x)})", zorder=4)
m2 = q[q["pathway_modules"].fillna("").astype(str).str.split(";").apply(lambda L: "2" in L)]
x, y = ecdf(m2["frac_remaining_15h"])
axA.plot(x, y, color=MODULE_COLOR[2], lw=2.0, label=f"M2 module (n={len(x)})", zorder=5)
axA.set_xlim(0, 1.4)
axA.set_ylim(0, 1)
axA.set_xlabel("protein remaining at 15 h  (ratio to 0 h)", fontsize=9.5, color=INK2)
axA.set_ylabel("cumulative fraction", fontsize=9.5, color=INK2)
axA.grid(color=GRID, lw=0.5)
axA.set_axisbelow(True)
axA.legend(frameon=False, fontsize=8.6, loc="lower right", labelcolor=INK2)
axA.set_title("A   Best hits are no less stable than the proteome",
              fontsize=11.5, fontweight="bold", loc="left", pad=8, color=INK)

# ---------------------------------------------------------------- stats box
axL.set_axis_off()
c15, ck = summ["fraction_remaining_15h"], summ["decay_rate_k"]
hl = summ["derived_halflife"]
lines = [
    ("protein remaining at 15 h", f"hits {c15['median_hits']:.3f}   background "
                                  f"{c15['median_background']:.3f}   p = {c15['p_value']:.2f}"),
    ("decay rate k (1/h)", f"hits {ck['median_hits']:.4f}   background "
                           f"{ck['median_background']:.4f}   p = {ck['p_value']:.2f}"),
    ("derived median half-life", f"hits {hl['hits_median_h']:.0f} h   background "
                                 f"{hl['background_median_h']:.0f} h"),
    ("M2 module vs background", f"0.576 vs {c15['median_background']:.3f}   "
                               f"p = {summ['modules'][1]['p_vs_background']:.3f} "
                               f"(Bonferroni {summ['modules'][1]['p_bonferroni']:.2f})"),
]
axL.text(0, 1.0, "Medians, Mann–Whitney U (two-sided)", fontsize=10.5,
         fontweight="bold", color=INK, va="top")
for i, (k, v) in enumerate(lines):
    axL.text(0, 0.78 - i * 0.175, k, fontsize=9.3, color=INK2, va="top")
    axL.text(0, 0.705 - i * 0.175, v, fontsize=9.8, color=INK, va="top")
axL.text(0, 0.0,
         f"{hl['hits_beyond_15h_window']} of {summ['quantified']} hits decay too slowly to reach "
         f"50% within the\n15 h assay, so their half-lives are extrapolations.",
         fontsize=8.6, color=INK3, va="bottom", style="italic")

# ---------------------------------------------------------------- panel B
o = q.sort_values("frac_remaining_15h").reset_index(drop=True)
bmed = bg["frac_remaining_15h"].median()
b25, b75 = bg["frac_remaining_15h"].quantile([0.25, 0.75])
axB.axvspan(b25, b75, color=BGFILL, alpha=0.35, zorder=0,
            label=f"background IQR ({b25:.2f}–{b75:.2f})")
axB.axvline(bmed, color=INK3, lw=1.2, ls=(0, (4, 3)), zorder=2,
            label=f"background median ({bmed:.2f})")
for r in o.itertuples():
    c = color_of(r.mod)
    axB.plot([0, r.frac_remaining_15h], [r.Index, r.Index], color=c, lw=1.4,
             alpha=0.55, zorder=3, solid_capstyle="round")
    axB.scatter(r.frac_remaining_15h, r.Index, s=68, color=c, edgecolor="white",
                linewidth=1.2, zorder=4)
axB.set_yticks(range(len(o)))
axB.set_yticklabels(o["gene"], fontsize=8.3)
for tick, m in zip(axB.get_yticklabels(), o["mod"]):
    tick.set_color(color_of(m))
    tick.set_fontweight("bold")
axB.set_ylim(-0.7, len(o) - 0.3)
axB.set_xlim(0, 1.32)
axB.set_xlabel("protein remaining at 15 h  (ratio to 0 h)", fontsize=9.5, color=INK2)
axB.grid(axis="x", color=GRID, lw=0.5)
axB.set_axisbelow(True)
axB.set_title("B   Every quantified best hit, ranked by stability",
              fontsize=11.5, fontweight="bold", loc="left", pad=8, color=INK)
handles = [Line2D([], [], marker="o", ls="", markersize=8, color=MODULE_COLOR[m],
                  label=MODULE_NAME[m]) for m in sorted(MODULE_COLOR)]
handles.append(Line2D([], [], marker="o", ls="", markersize=8, color=NOMOD,
                      label="no specific pathway module"))
axB.legend(handles=handles + list(axB.get_legend_handles_labels()[0][:0]),
           frameon=False, fontsize=8.8, loc="lower right", labelcolor=INK2,
           ncol=1, handletextpad=0.6)

# ---------------------------------------------------------------- panel C
bgcurve = [bg[[RCOLS[i], RCOLS[i + 1]]].stack().median() for i in (0, 2, 4, 6)]
axC.plot(TPTS, bgcurve, color=BGFILL, lw=3.0, marker="o", markersize=7,
         markeredgecolor="white", markeredgewidth=1.4, zorder=2,
         label="background median")
worst = o.head(8)["gene"].tolist()
cmap = dict(zip(q["gene"], q["mod"]))
curve_handles = []
for g in worst:
    row = curves[curves["gene"] == g]
    if row.empty:
        continue
    row = row.iloc[0]
    yv = [np.mean([row[RCOLS[i]], row[RCOLS[i + 1]]]) for i in (0, 2, 4, 6)]
    c = color_of(cmap.get(g, 0))
    ln, = axC.plot(TPTS, yv, color=c, lw=1.7, marker="o", markersize=5.5,
                   markeredgecolor="white", markeredgewidth=1.0, zorder=3,
                   alpha=0.9, label=g)
    curve_handles.append(ln)
axC.axhline(0.5, color=INK3, lw=1.0, ls=(0, (3, 2)), zorder=1)
axC.text(0.15, 0.515, "50% remaining", fontsize=8.2, color=INK3, va="bottom")
axC.set_xticks(TPTS)
axC.set_xlim(-0.4, 15.4)
axC.set_ylim(0, 1.32)
axC.set_xlabel("time in extract (h)", fontsize=9.5, color=INK2)
axC.set_ylabel("ratio to 0 h", fontsize=9.5, color=INK2)
axC.grid(color=GRID, lw=0.5)
axC.set_axisbelow(True)
bg_handle = Line2D([], [], color=BGFILL, lw=3.0, marker="o", markersize=7,
                   markeredgecolor="white", label="background median")
axC.legend(handles=[bg_handle] + curve_handles, frameon=False, fontsize=8.6,
           loc="upper right", ncol=3, labelcolor="linecolor",
           handletextpad=0.6, columnspacing=1.6)
axC.set_title("C   Measured decay curves of the eight least stable best hits",
              fontsize=11.5, fontweight="bold", loc="left", pad=8, color=INK)

fig.suptitle("Stability of the UBR3-WT co-IP best hits (Yi et al. 2023 degradomics)",
             fontsize=14.5, fontweight="bold", x=0.055, ha="left", y=0.990, color=INK)
fig.text(0.055, 0.9665,
         "TMTpro-MS3 time course in WT extract, untreated · 45 of the 50 WT best hits quantified · "
         "no half-life column exists in the source table; half-lives fitted here",
         fontsize=9.2, color=INK2, ha="left")

fig.subplots_adjust(left=0.085, right=0.965, top=0.912, bottom=0.05)
for ext in ("png", "pdf", "svg"):
    fig.savefig(HERE / f"WT_besthits_halflife.{ext}", dpi=200, bbox_inches="tight",
                facecolor="white")
print("figure -> WT_besthits_halflife.png/pdf/svg")
