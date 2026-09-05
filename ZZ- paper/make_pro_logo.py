"""
Logos comparing stabilised Met-Pro peptides with unstabilised Met-Pro peptides.

Both groups start Met-Pro, so the proline cannot separate them; positions 3-24
are the question. Drawn in the layout of ../ubr3_PG_reviewer/ -- black glyph
outlines on the fold-change logos, white hairlines on the significance logos,
no gridlines -- but in THIS project's blue/purple/pink palette rather than that
project's warm one. The two are different papers and are painted apart. The
palette, the three charge groups and the reasoning behind both live in
zz_pro_motif.py.

  FigureZ1_pro_foldchange_logo        letter height = fold change (substrates /
                                      controls), enriched residues only
  FigureZ1b_..._significant           the same, restricted to p < 0.05
  FigureZ2_pro_significance_logo      letter height AND colour = -log10 p
                                      (Fisher exact), the version to read
  FigureZ2b_..._strict                the same for the ZYG11B + double-KO
                                      subset, without the ZER1-only calls
  data/pro_motif_cells.csv            every cell, both groups, p and q

READ THE CAPTION BEFORE THE GLYPHS. 25 of the 440 cells reach p < 0.05 against
22 expected by chance alone, and not one survives BH-FDR. These figures are
published here as a negative result: after conditioning on Met-Pro, nothing in
positions 3-24 separates the stabilised peptides from the rest. The letters are
drawn because "there is nothing here" is a claim a reader should be able to
check, not because any of them is a finding.
"""
import os

import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap, Normalize
import logomaker

import zz_pro_motif as zp

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
DATA = os.path.join(HERE, "data")
os.makedirs(FIG, exist_ok=True)

POS = zp.POS
AA_COLOR = {a: zp.CAT_COLORS[c] for c, mem in zp.CATEGORY_MEMBERS.items() for a in mem}
CMAP_SIG = LinearSegmentedColormap.from_list("sig_warm", zp.CMAP_STOPS)
RAMP = LinearSegmentedColormap.from_list("sig_glyph", CMAP_SIG(np.linspace(0.34, 1.0, 256)))
VMIN, VMAX = 1.30103, 4.0
NORM = Normalize(vmin=VMIN, vmax=VMAX, clip=True)

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8, "axes.linewidth": 0.8,
    "svg.fonttype": "none", "pdf.fonttype": 42,
    "figure.dpi": 600, "savefig.dpi": 600, "savefig.bbox": "tight",
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "axes.spines.top": False, "axes.spines.right": False,
})


def save(fig, stem):
    for ext in ("png", "pdf", "svg"):
        p = os.path.join(FIG, f"{stem}.{ext}")
        fig.savefig(p, format=ext)
        print(f"  saved: {os.path.basename(p)}  ({os.path.getsize(p) / 1024:.1f} KB)")
    plt.close(fig)


def counts(table):
    k = table[table.defined]
    return (len(k), int((k.p_fisher < 0.05).sum()), int((k.q_fisher < 0.05).sum()),
            int(table.n_group_sub.iloc[0]), int(table.n_group_ctrl.iloc[0]))


def caption(table):
    n_cell, n_hit, n_fdr, n_s, n_c = counts(table)
    return (f"n = {n_s} stabilised Met-Pro peptides  vs  n = {n_c} Met-Pro peptides that are "
            f"not stabilised.   {n_hit} of {n_cell} cells reach $p$ < 0.05 against "
            f"{0.05 * n_cell:.0f} expected by chance alone, and "
            + ("none survives BH-FDR — nothing here is separable from noise."
               if not n_fdr else f"{n_fdr} survive BH-FDR."))


# --------------------------------------------------------------- fold change
def foldchange_logo(table, title, stem, only_sig=False):
    only = table.pivot(index="label", columns="position", values="p_fisher").reindex(
        index=zp.AA, columns=POS).lt(0.05) if only_sig else None
    H = zp.logo_frame("fold_change", table=table, only=only)

    fig, ax = plt.subplots(figsize=(10.5, 2.9))
    logomaker.Logo(H, ax=ax, color_scheme=AA_COLOR, font_name="Helvetica",
                   vpad=0.0, width=0.95, stack_order="big_on_top")
    for patch in ax.patches:
        patch.set_edgecolor("black"); patch.set_linewidth(0.35)

    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_ylabel("Fold change (stabilised / not)", fontsize=8)
    ax.set_xlabel("Position in the 24-mer", fontsize=8, labelpad=2)
    ax.set_xticks(POS); ax.set_xticklabels(POS, fontsize=6.5)
    ax.tick_params(labelsize=6.5, length=2.5, width=0.7, pad=2)
    ax.set_title(title, fontweight="bold", fontsize=8.5, pad=4)
    ax.set_ylim(-0.3, float(H.sum(axis=1).max()) + 2.0)
    ax.set_xlim(POS[0] - 0.6, POS[-1] + 0.6)
    ax.text(0.5, -0.155, caption(table), transform=ax.transAxes, ha="center",
            fontsize=5.6, color=zp.MUTED, style="italic")
    ax.legend(handles=[mpatches.Patch(color=zp.CAT_COLORS[c],
                                      label=f'{c}  ({" ".join(zp.CATEGORY_MEMBERS[c])})')
                       for c in zp.CAT_ORDER],
              fontsize=5.5, frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0),
              handlelength=0.8, handletextpad=0.3, borderpad=0.2)
    plt.tight_layout(pad=0.3)
    save(fig, stem)


# -------------------------------------------------------------- significance
def significance_logo(table, title, stem):
    P = table.pivot(index="label", columns="position", values="p_fisher").reindex(
        index=zp.AA, columns=POS)
    D = table.pivot(index="label", columns="position", values="direction").reindex(
        index=zp.AA, columns=POS)
    with np.errstate(divide="ignore", invalid="ignore"):
        H = -np.log10(P.astype(float))
    H = H.where(np.isfinite(H) & D.gt(0) & P.lt(0.05), 0.0).clip(lower=0.0).T.set_axis(POS)

    fig, ax = plt.subplots(figsize=(10.5, 2.9))
    logo = logomaker.Logo(H, ax=ax, color_scheme="#CFCAC1", font_name="Helvetica",
                          vpad=0.0, width=0.95, stack_order="big_on_top")
    for p_ in H.index:
        for a in H.columns:
            h = float(H.at[p_, a])
            if h > 0:
                logo.style_single_glyph(p_, a, color=RAMP(NORM(h)))
    for patch in ax.patches:
        patch.set_edgecolor(zp.GRID); patch.set_linewidth(0.35)

    for y, lab in [(1.30103, "p = 0.05"), (2.0, "p = 0.01"), (3.0, "p = 0.001")]:
        ax.axhline(y, color=zp.RULE, linewidth=0.5, linestyle=(0, (3, 2)), zorder=1)
        ax.text(POS[-1] + 0.55, y, lab, fontsize=5, color=zp.MUTED, va="center", ha="left")

    ax.axhline(0, color=zp.RULE, linewidth=0.6)
    ax.set_ylabel("$-$log$_{10}$ $p$   (stacked)", fontsize=8, color=zp.INK)
    ax.set_xlabel("Position in the 24-mer", fontsize=8, labelpad=2, color=zp.INK)
    ax.set_xticks(POS); ax.set_xticklabels(POS, fontsize=6.5)
    ax.tick_params(labelsize=6.5, length=2.5, width=0.7, pad=2, colors=zp.INK)
    for sp in ax.spines.values():
        sp.set_edgecolor(zp.RULE)
    ax.set_title(title, fontweight="bold", fontsize=8.5, pad=4, color=zp.INK)
    ax.set_ylim(-0.05, max(float(H.sum(axis=1).max()) * 1.12, 3.6))
    ax.set_xlim(POS[0] - 0.6, POS[-1] + 0.6)
    ax.text(0.5, -0.155, caption(table), transform=ax.transAxes, ha="center",
            fontsize=5.6, color=zp.MUTED, style="italic")

    sm = mpl.cm.ScalarMappable(cmap=RAMP, norm=NORM)
    cbar = fig.colorbar(sm, ax=ax, fraction=0.020, pad=0.055)
    cbar.ax.tick_params(labelsize=6, length=2, width=0.6, colors=zp.INK)
    cbar.outline.set_linewidth(0.5); cbar.outline.set_edgecolor(zp.RULE)
    cbar.set_ticks([VMIN, 2, 3, VMAX]); cbar.set_ticklabels(["0.05", "0.01", "0.001", "1e-4"])
    cbar.set_label("$p$ value", fontsize=6.5, labelpad=5, color=zp.INK)
    plt.tight_layout(pad=0.3)
    save(fig, stem)


print("building Met-Pro logos\n")
union, strict = zp.cells(False), zp.cells(True)
out = union.assign(set="ZYG11B/ZER1/DKO")._append(strict.assign(set="ZYG11B/DKO only")) \
    if hasattr(union, "_append") else None
import pandas as pd
pd.concat([union.assign(hit_set="ZYG11B_ZER1_DKO"),
           strict.assign(hit_set="ZYG11B_DKO_only")]).to_csv(
    os.path.join(DATA, "pro_motif_cells.csv"), index=False, float_format="%.6g")

foldchange_logo(union,
                "Met-Pro peptides — fold change, stabilised vs not stabilised "
                "(ZYG11B / ZER1 / double KO)",
                "FigureZ1_pro_foldchange_logo")
foldchange_logo(union,
                "Met-Pro peptides — fold change, residues at $p$ < 0.05 only",
                "FigureZ1b_pro_foldchange_logo_significant", only_sig=True)
significance_logo(union,
                  "Met-Pro peptides — height and colour are $-$log$_{10}$ $p$ "
                  "(Fisher exact), stabilised vs not",
                  "FigureZ2_pro_significance_logo")
significance_logo(strict,
                  "Met-Pro peptides — $-$log$_{10}$ $p$, ZYG11B and double KO only "
                  "(ZER1-only calls dropped)",
                  "FigureZ2b_pro_significance_logo_strict")

print("\n" + "=" * 78)
print("is there anything beyond the proline?")
print("=" * 78)
for tag, t in (("ZYG11B / ZER1 / DKO", union), ("ZYG11B + DKO only", strict)):
    n_cell, n_hit, n_fdr, n_s, n_c = counts(t)
    print(f"\n{tag}:  {n_s} vs {n_c} peptides")
    print(f"   {n_hit} of {n_cell} cells at p < 0.05, {0.05 * n_cell:.0f} expected by chance, "
          f"{n_fdr} survive BH-FDR")
    k = t[t.defined & (t.direction > 0)].nsmallest(5, "p_fisher")
    print(k[["position", "label", "n_sub", "n_ctrl", "fold_change", "p_fisher", "q_fisher"]]
          .to_string(index=False, formatters={"fold_change": "{:.2f}".format,
                                              "p_fisher": "{:.2e}".format,
                                              "q_fisher": "{:.2f}".format}))
