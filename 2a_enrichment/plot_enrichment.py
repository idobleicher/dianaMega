"""
Figure for the UBR3 P-[D/E] enrichment analysis.

A  The counting unit decides the answer. D/E at P2 rises from 2.6% to 12.1% as
   database redundancy is removed, crossing the background lines on the way.
B  Full amino-acid composition at P2 (deduplicated) against the three backgrounds.
C  Fold change with 95% CI for the primary unit against each background.

Outputs: enrichment_figure.(png|pdf|svg)
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

FG = "#2a78d6"
BG_COLOR = {
    "All amino acids, all positions": "#8a8880",
    "Position 2 of viral proteins": "#eb6834",
    "Residue after any proline (P+1)": "#1baf7a",
}
BG_SHORT = {
    "All amino acids, all positions": "all AA, all positions",
    "Position 2 of viral proteins": "viral protein position 2",
    "Residue after any proline (P+1)": "residue after any P",
}
INK, INK2, INK3, GRID = "#0b0b0b", "#52514e", "#8a8880", "#e2e1dc"
AA = list("ACDEFGHIKLMNPQRSTVWY")

plt.rcParams.update({
    "font.family": "DejaVu Sans", "axes.edgecolor": GRID, "axes.linewidth": 0.8,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
    "xtick.major.size": 0, "ytick.major.size": 0, "figure.facecolor": "white",
})

res = pd.read_csv(HERE / "enrichment_results.csv")
comp = pd.read_csv(HERE / "enrichment_composition.csv")
summ = json.loads((HERE / "enrichment_summary.json").read_text())
PRIMARY = summ["primary_foreground"]

order = ["All 2A instances (raw)", "Unique source accession", "Unique 2A core peptide",
         "Unique downstream context"]
uni = res.drop_duplicates("foreground").set_index("foreground")

fig = plt.figure(figsize=(14.0, 11.4))
gs = GridSpec(2, 2, height_ratios=[1.0, 1.0], width_ratios=[1.05, 1],
              hspace=0.40, wspace=0.22, figure=fig)
axA = fig.add_subplot(gs[0, 0])
axC = fig.add_subplot(gs[0, 1])
axB = fig.add_subplot(gs[1, :])

# ---------------------------------------------------------------- panel A
x = np.arange(len(order))
vals = [uni.loc[o, "observed_DE_pct"] for o in order]
lo = [uni.loc[o, "observed_CI95_low_pct"] for o in order]
hi = [uni.loc[o, "observed_CI95_high_pct"] for o in order]
axA.bar(x, vals, color=FG, width=0.6, zorder=3)
axA.errorbar(x, vals, yerr=[np.array(vals) - lo, np.array(hi) - np.array(vals)],
             fmt="none", ecolor=INK, elinewidth=1.4, capsize=5, zorder=4)
for b, v, h, n in zip(x, vals, hi, [uni.loc[o, "foreground_n"] for o in order]):
    axA.text(b, h + 0.6, f"{v:.1f}%", ha="center", fontsize=10, fontweight="bold",
             color=INK, zorder=5)
    axA.text(b, 0.5, f"n={n:,}", ha="center", fontsize=8.2, color="white", zorder=5)

bg_handles = []
for bname, c in BG_COLOR.items():
    y = res[res["background"] == bname]["background_DE_pct"].iloc[0]
    axA.axhline(y, color=c, lw=1.8, ls=(0, (5, 3)), zorder=2)
    bg_handles.append(Line2D([], [], color=c, lw=1.8, ls=(0, (5, 3)),
                             label=f"{BG_SHORT[bname]} — {y:.1f}%"))
axA.legend(handles=bg_handles, frameon=False, fontsize=8.6, loc="upper left",
           bbox_to_anchor=(0.0, 0.82), labelcolor=INK2, title="background",
           title_fontsize=8.6, handlelength=2.4, borderpad=0.2)

axA.set_xticks(x)
axA.set_xticklabels(["all\ninstances", "unique\naccession", "unique 2A\ncore",
                     "unique downstream\ncontext"], fontsize=8.8)
axA.set_ylabel("D or E at position 2  (%)", fontsize=10, color=INK2)
axA.set_ylim(0, 24)
axA.grid(axis="y", color=GRID, lw=0.5)
axA.set_axisbelow(True)
axA.set_title("A   Removing database redundancy changes the answer",
              fontsize=12, fontweight="bold", loc="left", pad=10, color=INK)

# ---------------------------------------------------------------- panel C
sub = res[res["foreground"] == PRIMARY].copy()
o_lo = uni.loc[PRIMARY, "observed_CI95_low_pct"]
o_hi = uni.loc[PRIMARY, "observed_CI95_high_pct"]
y = np.arange(len(sub))
for i, r in enumerate(sub.itertuples()):
    c = BG_COLOR[r.background]
    f = r.fold_change
    flo, fhi = o_lo / r.background_DE_pct, o_hi / r.background_DE_pct
    axC.plot([flo, fhi], [i, i], color=c, lw=2.6, solid_capstyle="round", zorder=3)
    axC.scatter(f, i, s=150, color=c, edgecolor="white", linewidth=1.6, zorder=4)
    sig = "n.s." if r.fisher_p > 0.05 else (f"p = {r.fisher_p:.0e}")
    axC.text(fhi + 0.06, i, f"{f:.2f}×   {sig}", va="center", fontsize=9.4,
             color=INK, fontweight="bold")
axC.axvline(1.0, color=INK, lw=1.4, zorder=2)
axC.text(1.02, len(sub) - 0.35, "no difference", fontsize=8.6, color=INK2)
axC.set_yticks(y)
axC.set_yticklabels([BG_SHORT[b] for b in sub["background"]], fontsize=9.4)
axC.set_xlim(0.3, 1.9)
axC.set_ylim(-0.6, len(sub) - 0.3)
axC.set_xlabel("fold change vs background   (observed ÷ expected)",
               fontsize=10, color=INK2)
axC.grid(axis="x", color=GRID, lw=0.5)
axC.set_axisbelow(True)
axC.set_title(f"C   {PRIMARY} vs each background",
              fontsize=12, fontweight="bold", loc="left", pad=10, color=INK)

# ---------------------------------------------------------------- panel B
xb = np.arange(len(AA))
w = 0.21
axB.bar(xb - 1.5 * w, comp[PRIMARY], width=w, color=FG, zorder=3,
        label=f"2A position 2 — {PRIMARY.lower()} (n={uni.loc[PRIMARY,'foreground_n']:,})")
for j, (bname, c) in enumerate(BG_COLOR.items()):
    axB.bar(xb + (j - 0.5) * w, comp[bname], width=w, color=c, zorder=3,
            label=f"background: {BG_SHORT[bname]}")
for i, a in enumerate(AA):
    if a in ("D", "E"):
        axB.axvspan(i - 0.5, i + 0.5, color="#eda100", alpha=0.13, zorder=0)
axB.set_xticks(xb)
axB.set_xticklabels(AA, fontsize=10.5, fontweight="bold")
for t in axB.get_xticklabels():
    if t.get_text() in ("D", "E"):
        t.set_color("#b06d00")
axB.set_ylabel("frequency (%)", fontsize=10, color=INK2)
axB.set_xlabel("amino acid at position 2 (the residue UBR3 reads)",
               fontsize=10, color=INK2)
axB.set_ylim(0, 21)
axB.grid(axis="y", color=GRID, lw=0.5)
axB.set_axisbelow(True)
axB.legend(frameon=False, fontsize=9.2, loc="upper right", labelcolor=INK2, ncol=2)
axB.set_title("B   Full composition at position 2 — D and E highlighted",
              fontsize=12, fontweight="bold", loc="left", pad=10, color=INK)

fig.suptitle("Is the UBR3 P–[D/E] motif enriched after the 2A ribosomal skip?",
             fontsize=15, fontweight="bold", x=0.045, ha="left", y=0.985, color=INK)
fig.text(0.045, 0.960,
         "Foreground: 2A downstream products (Rao et al. 2025 combined set) · "
         "Background: 17,517 reviewed viral proteins, 8.0 M residues (UniProt Swiss-Prot)",
         fontsize=9.6, color=INK2, ha="left")

fig.subplots_adjust(left=0.07, right=0.975, top=0.925, bottom=0.065)
for ext in ("png", "pdf", "svg"):
    fig.savefig(HERE / f"enrichment_figure.{ext}", dpi=200, bbox_inches="tight",
                facecolor="white")
print("figure -> enrichment_figure.png/pdf/svg")
