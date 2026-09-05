"""
Per-peptide sort-bin profiles, in the layout of the UBR3 manuscript's
Supplementary Fig. 1 (Fig_V3.pptx, slide 6).

That figure is a small-multiple grid: one panel per reporter peptide, titled
with the gene and its N-terminal triplet, plotting the PROPORTION OF READS in
each sorting bin, control in grey and knockout in colour. A peptide that is
degraded in the control sits in the low bins; when the E3 that degrades it is
removed, the whole curve slides to the right. It is the same measurement PSI
summarises -- PSI is exactly the read-weighted mean of these bins, asserted on
load -- but shown at full resolution instead of collapsed to one number, so a
reader can see whether a dPSI comes from a clean shift of the population or
from noise in one bin.

DIFFERENCES FROM THE ORIGINAL, all forced by this dataset:
  * six bins, not four -- this screen sorts into six.
  * five lines, not four. The original has two controls and two knockout
    replicates. Here there are no within-genotype replicates at all, so the
    lines are the two experiments' controls and the three knockouts:
        AAVS1 control      dark grey    Data S3A
        wild-type control  light grey   Data S3B
        ZYG11B KO          blue         S3A, against the AAVS1 control
        ZER1 KO            blue dashed  S3A, same control, same complex
        double KO          pink         S3B, against the wild-type control
    ZER1 is dashed rather than a third hue: it is ZYG11B's paralogue in the
    same CRL2 complex, the dash says so, and the project's palette is two
    hues on purpose (see zz_pro_motif.py).
  * the grey pair is not a replicate pair. Read each knockout against the
    control of ITS OWN experiment -- blue against dark grey, pink against
    light grey. Comparing blue to light grey compares two experiments.

WHICH PEPTIDES. Twelve, in three rows of four:
  row 1  the two peptides confirmed in the lab, then the two strongest other
         Met-Pro-Gly candidates
  row 2  four more Met-Pro-Gly
  row 3  four canonical Met-Gly substrates, for comparison -- the shape a
         known Gly/N-degron substrate makes in this assay
Met-Pro-Gly is the group the lab-confirmed peptides both belong to; see the
README for why it behaves like Met-Gly.

Outputs:
  figures/FigureZ4_bin_profiles          the twelve panels
  data/bin_profiles_plotted.csv          the numbers behind every line
"""
import os

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import zz_gps_data as zz

HERE = os.path.dirname(os.path.abspath(__file__))
FIG, DATA = os.path.join(HERE, "figures"), os.path.join(HERE, "data")
os.makedirs(FIG, exist_ok=True)

BLUE, PINK = "#2E6FD0", "#D6408F"
GREY_A, GREY_B = "#6E6E76", "#B6B6BE"
INK, MUTED, RULE = "#1A1A1A", "#6B6B6B", "#C6C6CE"

# name, colour, linestyle, which control it belongs with
SERIES = [
    ("AAVS1 control", GREY_A, "-", 1.3),
    ("ZYG11B KO", BLUE, "-", 1.5),
    ("ZER1 KO", BLUE, (0, (3, 1.6)), 1.5),
    ("wild-type control", GREY_B, "-", 1.3),
    ("double KO", PINK, "-", 1.5),
]

# transcript given where a specific isoform is meant; otherwise the gene's
# best-scoring transcript in the double KO is used
PANELS = [
    ("VWA5B1", "ENST00000485375.6"), ("SEPTIN4", "ENST00000412945.7"),
    ("FOSB", None), ("SSBP3", None),
    ("ZNF254", None), ("PLD4", None), ("ZNF729", None), ("NDUFA8", None),
    ("C14orf178", None), ("FXR2", None), ("UBTD2", None), ("HBG1", None),
]
NCOL = 4

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8, "axes.linewidth": 0.8,
    "svg.fonttype": "none", "pdf.fonttype": 42,
    "figure.dpi": 600, "savefig.dpi": 600, "savefig.bbox": "tight",
    "figure.facecolor": "white", "savefig.facecolor": "white",
})

d = zz.load()
b = zz.bins()
merged = b.merge(d[["transcript", "dpsi_zyg11b", "dpsi_zer1", "dpsi_dko",
                    "psi_ctrl_a", "psi_ctrl_b"]], on="transcript", how="left")


def resolve(gene, transcript):
    rows = merged[merged.gene == gene]
    assert len(rows), f"{gene} is not in the screen"
    if transcript:
        r = rows[rows.transcript == transcript]
        assert len(r), f"{gene}: transcript {transcript} not found"
        return r.iloc[0]
    return rows.loc[rows.dpsi_dko.idxmax()]


rows = [resolve(g, t) for g, t in PANELS]

# every number that reaches the figure, written out beside it
out = []
for r in rows:
    for name, *_ in SERIES:
        for i in zz.BINS:
            out.append({"gene": r.gene, "transcript": r.transcript,
                        "nterm3": r.nterm3, "condition": name, "bin": i,
                        "proportion": r[f"{name} bin{i}"],
                        "reads": r[f"{name} reads"]})
pd.DataFrame(out).to_csv(os.path.join(DATA, "bin_profiles_plotted.csv"),
                         index=False, float_format="%.5g")

nrow = int(np.ceil(len(rows) / NCOL))
fig, axes = plt.subplots(nrow, NCOL, figsize=(10.6, 2.05 * nrow))
axes = np.atleast_2d(axes)

for k, r in enumerate(rows):
    ax = axes[k // NCOL][k % NCOL]
    for name, colour, style, lw in SERIES:
        y = [r[f"{name} bin{i}"] for i in zz.BINS]
        if not np.isfinite(y).all():
            continue                      # peptide absent from that experiment
        ax.plot(zz.BINS, y, color=colour, linestyle=style, linewidth=lw,
                marker="o", markersize=2.6, markeredgewidth=0, zorder=3,
                clip_on=False)

    ax.set_title(f"{r.gene} ({r.nterm3}-)", fontsize=8, fontweight="bold",
                 color=INK, pad=9)
    dp = "   ".join(
        f"{lab} {r[c]:+.2f}" for lab, c in
        (("ZYG", "dpsi_zyg11b"), ("ZER", "dpsi_zer1"), ("DKO", "dpsi_dko"))
        if np.isfinite(r[c]))
    ax.text(0.5, 1.015, f"ΔPSI   {dp}", transform=ax.transAxes, ha="center",
            va="bottom", fontsize=5.2, color=MUTED)

    ax.set_xticks(zz.BINS)
    ax.set_xticklabels([f"Bin{i}" for i in zz.BINS], fontsize=5.6, color=MUTED)
    ax.tick_params(axis="y", labelsize=6, length=2.5, width=0.7, colors=INK)
    ax.tick_params(axis="x", length=0, pad=3)
    ax.set_ylabel("Proportion of reads", fontsize=6.4, color=INK, labelpad=3)
    ax.set_ylim(0, None)
    ax.margins(x=0.06)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(RULE)

for k in range(len(rows), nrow * NCOL):
    axes[k // NCOL][k % NCOL].set_visible(False)

handles = [Line2D([0], [0], color=c, linestyle=s, linewidth=lw, marker="o",
                  markersize=3, markeredgewidth=0, label=n)
           for n, c, s, lw in SERIES]
fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
           fontsize=6.6, handlelength=2.2, handletextpad=0.5,
           columnspacing=2.0, bbox_to_anchor=(0.5, -0.005))

fig.suptitle("Sort-bin profiles of individual peptides", fontsize=10,
             fontweight="bold", color=INK, y=1.005)
fig.text(0.5, 0.978,
         "proportion of reads in each of the six stability bins · low bins = degraded, "
         "high bins = stable · a knockout that stabilises a peptide moves its curve to the right",
         ha="center", fontsize=6.4, color=MUTED, style="italic")
fig.text(0.5, -0.035,
         "Read each knockout against the control of its own experiment: blue and blue-dashed "
         "against dark grey (Data S3A), pink against light grey (Data S3B).  The two greys are "
         "different experiments, not replicates, and there is no within-genotype replicate "
         "anywhere in this workbook.\nPSI is the read-weighted mean of these six bins, so these "
         "curves and every ΔPSI quoted in this project are one measurement at two resolutions.  "
         "Rows 1–2 are Met-Pro-Gly peptides, row 3 canonical Met-Gly substrates for comparison.",
         ha="center", va="top", fontsize=5.6, color=MUTED, linespacing=1.6)

plt.tight_layout(rect=(0, 0.03, 1, 0.965), h_pad=2.4, w_pad=1.6)
for ext in ("png", "pdf", "svg"):
    path = os.path.join(FIG, f"FigureZ4_bin_profiles.{ext}")
    fig.savefig(path, format=ext)
    print(f"  saved: {os.path.basename(path)}  ({os.path.getsize(path)/1024:.1f} KB)")
plt.close(fig)

print("\npanels drawn:")
for r in rows:
    print(f"   {r.gene:<10} {r.nterm3}-  {r.transcript:<20} "
          f"ZYG {r.dpsi_zyg11b:+.2f}  ZER {r.dpsi_zer1:+.2f}  DKO {r.dpsi_dko:+.2f}")
