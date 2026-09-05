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

TWO FIGURES COME OUT OF THIS
  FigureZ4  twelve curated panels: the two peptides confirmed in the lab, the
            strongest other Met-Pro-Gly candidates, and four canonical Met-Gly
            substrates on the bottom row for comparison.
  FigureZ5  EVERY Met-Pro substrate -- all 180 of tiers A, B and C from
            make_mp_substrate_list.py -- paginated 20 to a page, in tier order
            and ranked by best dPSI within each tier. One multi-page PDF plus
            a PNG per page. Each panel carries its tier, so a page can be read
            without the tier table beside it.

Outputs:
  figures/FigureZ4_bin_profiles                    the twelve curated panels
  figures/FigureZ5_bin_profiles_all_MP_substrates.pdf   all 180, multi-page
  figures/FigureZ5_bin_profiles_all_MP_substrates_p01..pNN.png
  data/bin_profiles_plotted.csv                    numbers behind FigureZ4
"""
import os

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

import zz_gps_data as zz

HERE = os.path.dirname(os.path.abspath(__file__))
FIG, DATA = os.path.join(HERE, "figures"), os.path.join(HERE, "data")
os.makedirs(FIG, exist_ok=True)

BLUE, PINK = "#2E6FD0", "#D6408F"
GREY_A, GREY_B = "#6E6E76", "#B6B6BE"
INK, MUTED, RULE = "#1A1A1A", "#6B6B6B", "#C6C6CE"

SERIES = [
    ("AAVS1 control", GREY_A, "-", 1.3),
    ("ZYG11B KO", BLUE, "-", 1.5),
    ("ZER1 KO", BLUE, (0, (3, 1.6)), 1.5),
    ("wild-type control", GREY_B, "-", 1.3),
    ("double KO", PINK, "-", 1.5),
]

CURATED = [
    ("VWA5B1", "ENST00000485375.6"), ("SEPTIN4", "ENST00000412945.7"),
    ("FOSB", None), ("SSBP3", None),
    ("ZNF254", None), ("PLD4", None), ("ZNF729", None), ("NDUFA8", None),
    ("C14orf178", None), ("FXR2", None), ("UBTD2", None), ("HBG1", None),
]

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

LEGEND = [Line2D([0], [0], color=c, linestyle=s, linewidth=lw, marker="o",
                 markersize=3, markeredgewidth=0, label=n)
          for n, c, s, lw in SERIES]
FOOT = ("Read each knockout against the control of its own experiment: blue and blue-dashed "
        "against dark grey (Data S3A), pink against light grey (Data S3B).  The two greys are "
        "different experiments, not replicates, and there is no within-genotype replicate "
        "anywhere in this workbook.\nPSI is the read-weighted mean of these six bins, so these "
        "curves and every ΔPSI quoted in this project are one measurement at two resolutions.")


def draw_panel(ax, r, tier=None):
    for name, colour, style, lw in SERIES:
        y = [r[f"{name} bin{i}"] for i in zz.BINS]
        if not np.isfinite(y).all():
            continue                      # peptide absent from that experiment
        ax.plot(zz.BINS, y, color=colour, linestyle=style, linewidth=lw,
                marker="o", markersize=2.6, markeredgewidth=0, zorder=3,
                clip_on=False)

    title = f"{r.gene} ({r.nterm3}-)"
    ax.set_title(title, fontsize=8, fontweight="bold", color=INK, pad=9)
    dp = "   ".join(
        f"{lab} {r[c]:+.2f}" for lab, c in
        (("ZYG", "dpsi_zyg11b"), ("ZER", "dpsi_zer1"), ("DKO", "dpsi_dko"))
        if np.isfinite(r[c]))
    ax.text(0.5, 1.015, f"ΔPSI   {dp}", transform=ax.transAxes, ha="center",
            va="bottom", fontsize=5.2, color=MUTED)
    if tier:
        ax.text(0.0, 1.015, f"tier {tier}", transform=ax.transAxes, ha="left",
                va="bottom", fontsize=5.2, color=MUTED, fontweight="bold")

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


def grid(rows, tiers, title, subtitle, foot, ncol=4, panel_h=2.05, panel_w=2.65):
    nrow = int(np.ceil(len(rows) / ncol))
    fig, axes = plt.subplots(nrow, ncol,
                             figsize=(panel_w * ncol, panel_h * nrow))
    axes = np.atleast_2d(axes)
    for k, r in enumerate(rows):
        draw_panel(axes[k // ncol][k % ncol], r, tiers[k] if tiers else None)
    for k in range(len(rows), nrow * ncol):
        axes[k // ncol][k % ncol].set_visible(False)

    fig.legend(handles=LEGEND, loc="lower center", ncol=5, frameon=False,
               fontsize=6.6, handlelength=2.2, handletextpad=0.5,
               columnspacing=2.0, bbox_to_anchor=(0.5, -0.004))
    fig.suptitle(title, fontsize=10, fontweight="bold", color=INK, y=1.004)
    top = 1.0 - 0.026 * (2.05 * 3) / (panel_h * nrow)
    fig.text(0.5, top, subtitle, ha="center", fontsize=6.4, color=MUTED,
             style="italic")
    fig.text(0.5, -0.028, foot, ha="center", va="top", fontsize=5.6,
             color=MUTED, linespacing=1.6)
    plt.tight_layout(rect=(0, 0.03, 1, 1 - 0.035 * (2.05 * 3) / (panel_h * nrow)),
                     h_pad=2.4, w_pad=1.6)
    return fig


def save(fig, stem, exts=("png", "pdf", "svg")):
    for ext in exts:
        p = os.path.join(FIG, f"{stem}.{ext}")
        fig.savefig(p, format=ext)
        print(f"  saved: {os.path.basename(p)}  ({os.path.getsize(p)/1024:.1f} KB)")


# ----------------------------------------------------------- Figure Z4
def resolve(gene, transcript):
    rows = merged[merged.gene == gene]
    assert len(rows), f"{gene} is not in the screen"
    if transcript:
        r = rows[rows.transcript == transcript]
        assert len(r), f"{gene}: transcript {transcript} not found"
        return r.iloc[0]
    return rows.loc[rows.dpsi_dko.idxmax()]


rows = [resolve(g, t) for g, t in CURATED]
out = [{"gene": r.gene, "transcript": r.transcript, "nterm3": r.nterm3,
        "condition": name, "bin": i, "proportion": r[f"{name} bin{i}"],
        "reads": r[f"{name} reads"]}
       for r in rows for name, *_ in SERIES for i in zz.BINS]
pd.DataFrame(out).to_csv(os.path.join(DATA, "bin_profiles_plotted.csv"),
                         index=False, float_format="%.5g")

fig = grid(rows, None, "Sort-bin profiles of individual peptides",
           "proportion of reads in each of the six stability bins · low bins = degraded, "
           "high bins = stable · a knockout that stabilises a peptide moves its curve to the right",
           FOOT + "  Rows 1–2 are Met-Pro-Gly peptides, row 3 canonical Met-Gly substrates "
                  "for comparison.")
save(fig, "FigureZ4_bin_profiles")
plt.close(fig)

# ----------------------------------------------------------- Figure Z5
subs = pd.read_csv(os.path.join(DATA, "mp_substrates.csv"))
subs = subs.sort_values(["tier", "best_dpsi"], ascending=[True, False])
panels = merged.set_index("transcript").loc[subs.transcript].reset_index()
tiers = list(subs.tier)

PER_PAGE, NCOL = 20, 4
pages = int(np.ceil(len(panels) / PER_PAGE))
stem = "FigureZ5_bin_profiles_all_MP_substrates"
print(f"\nFigureZ5: {len(panels)} Met-Pro substrates over {pages} pages")

with PdfPages(os.path.join(FIG, f"{stem}.pdf")) as pdf:
    for pg in range(pages):
        sl = slice(pg * PER_PAGE, (pg + 1) * PER_PAGE)
        chunk = [panels.iloc[i] for i in range(len(panels))][sl]
        t_chunk = tiers[sl]
        counts = pd.Series(t_chunk).value_counts().reindex(list("ABC")).fillna(0)
        fig = grid(
            chunk, t_chunk,
            f"Sort-bin profiles — every Met-Pro substrate  ({pg + 1} of {pages})",
            f"tiers on this page: "
            + ", ".join(f"{int(counts[t])} × {t}" for t in "ABC" if counts[t])
            + "  ·  A = ΔPSI ≥ 1.0 in two or more genotypes, B = in exactly one, "
              "C = best 0.5–1.0  ·  ranked by best ΔPSI within tier",
            FOOT + "  Tier is how much evidence there is, not how good a substrate is: with no "
                   "within-genotype replicate, 0.6 in one clone is not evidence of absence.",
            ncol=NCOL, panel_h=2.05)
        pdf.savefig(fig, bbox_inches="tight")
        # the PDF is the archival artefact and is vector; the page PNGs exist
        # to be looked at, so they are written at screen resolution rather than
        # 600 dpi -- nine 600-dpi pages would add 25 MB to the repository for
        # nothing the PDF does not already hold.
        png = os.path.join(FIG, f"{stem}_p{pg + 1:02d}.png")
        fig.savefig(png, format="png", dpi=200)
        plt.close(fig)
        print(f"  page {pg + 1:>2}: {len(chunk):>2} panels  "
              f"({os.path.getsize(png)/1024:.0f} KB png)")

size = os.path.getsize(os.path.join(FIG, f"{stem}.pdf")) / 1024
print(f"  saved: {stem}.pdf  ({size:.1f} KB, {pages} pages)")
print(f"\ntier counts across the whole figure: "
      f"{ {t: int((subs.tier == t).sum()) for t in 'ABC'} }")
