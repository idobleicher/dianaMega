"""
Combined figure: Conservation plot (A) on top, domain architecture (B) below,
both sharing the exact same amino acid x-axis.
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.image as mpimg
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from collections import Counter
from Bio import AlignIO

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
ALIGNED_FASTA = os.path.join(OUTPUT_DIR, "ubr3_aligned.fasta")
ICONS_DIR = os.path.join(OUTPUT_DIR, "icons")

SPECIES_ORDER = [
    "H_sapiens", "M_musculus", "X_laevis", "P_marinus",
    "C_intestinalis", "B_floridae", "S_purpuratus", "O_sinensis",
    "D_melanogaster", "A_mellifera", "P_tepidariorum",
    "H_vulgaris", "N_vectensis",
]
ICON_FILES = {
    "H_sapiens": "icon_human.png", "M_musculus": "icon_mouse.png",
    "X_laevis": "icon_frog.png", "P_marinus": "icon_lamprey.png",
    "C_intestinalis": "icon_seasquirt.png", "B_floridae": "icon_amphioxus.png",
    "S_purpuratus": "icon_seaurchin.png", "O_sinensis": "icon_octopus.png",
    "D_melanogaster": "icon_fruitfly.png", "A_mellifera": "icon_honeybee.png",
    "P_tepidariorum": "icon_spider.png", "H_vulgaris": "icon_hydra.png",
    "N_vectensis": "icon_seaanemone.png",
}
DISPLAY_NAMES = {
    "H_sapiens": "H. sapiens", "M_musculus": "M. musculus",
    "X_laevis": "X. laevis", "P_marinus": "P. marinus",
    "C_intestinalis": "C. intestinalis", "B_floridae": "B. floridae",
    "S_purpuratus": "S. purpuratus", "O_sinensis": "O. sinensis",
    "D_melanogaster": "D. melanogaster", "A_mellifera": "A. mellifera",
    "P_tepidariorum": "P. tepidariorum", "H_vulgaris": "H. vulgaris",
    "N_vectensis": "N. vectensis",
}
COMMON_NAMES = {
    "H_sapiens": "Human", "M_musculus": "Mouse", "X_laevis": "Frog",
    "P_marinus": "Lamprey", "C_intestinalis": "Sea squirt",
    "B_floridae": "Amphioxus", "S_purpuratus": "Sea urchin",
    "O_sinensis": "Octopus", "D_melanogaster": "Fruit fly",
    "A_mellifera": "Honeybee", "P_tepidariorum": "Spider",
    "H_vulgaris": "Hydra", "N_vectensis": "Sea anemone",
}
PHYLUM_COLORS = {
    "H_sapiens": "#B71C1C", "M_musculus": "#B71C1C",
    "X_laevis": "#BF360C", "P_marinus": "#E65100",
    "C_intestinalis": "#F57F17", "B_floridae": "#827717",
    "S_purpuratus": "#33691E", "O_sinensis": "#00695C",
    "D_melanogaster": "#01579B", "A_mellifera": "#01579B",
    "P_tepidariorum": "#0D47A1", "H_vulgaris": "#4A148C",
    "N_vectensis": "#4A148C",
}

DOMAINS = [
    ("UBR box",       118,  189,  "#7E57C2"),
    ("Region 400-600", 400, 600,  "#1E88E5"),
    ("RING domain",   1232, 1376, "#EF6C00"),
]


def auto_crop(img, bg_thresh=245):
    gray = np.mean(img[:, :, :3], axis=2) if img.ndim == 3 else img
    mask = gray < bg_thresh
    rows = np.any(mask, axis=1)
    cols = np.any(mask, axis=0)
    if not np.any(rows) or not np.any(cols):
        return img
    r0, r1 = np.where(rows)[0][[0, -1]]
    c0, c1 = np.where(cols)[0][[0, -1]]
    pad = 6
    r0 = max(0, r0 - pad)
    r1 = min(img.shape[0], r1 + pad + 1)
    c0 = max(0, c0 - pad)
    c1 = min(img.shape[1], c1 + pad + 1)
    return img[r0:r1, c0:c1]


def load_icon(sp_id):
    fname = ICON_FILES.get(sp_id)
    if not fname:
        return None
    path = os.path.join(ICONS_DIR, fname)
    if not os.path.exists(path):
        return None
    img = mpimg.imread(path)
    if img.max() <= 1.0:
        img = (img * 255).astype(np.uint8)
    return auto_crop(img)


def main():
    print("=" * 60)
    print("Combined Conservation + Domain Architecture Figure")
    print("=" * 60)

    alignment = AlignIO.read(ALIGNED_FASTA, "fasta")
    n_cols = alignment.get_alignment_length()
    n_seqs = len(alignment)

    human_idx = 0
    for i, rec in enumerate(alignment):
        if "sapiens" in rec.id.lower():
            human_idx = i
            break

    col_to_human = {}
    human_to_col = {}
    res_num = 0
    for col in range(n_cols):
        if alignment[human_idx, col] != "-":
            res_num += 1
            col_to_human[col] = res_num
            human_to_col[res_num] = col
        else:
            col_to_human[col] = None
    human_len = res_num

    scores = np.zeros(n_cols)
    for col in range(n_cols):
        residues = [alignment[row, col] for row in range(n_seqs)]
        non_gap = [r for r in residues if r != "-"]
        if non_gap:
            counts = Counter(non_gap)
            scores[col] = counts.most_common(1)[0][1] / len(non_gap)

    human_x = []
    last_h = 0
    for col in range(n_cols):
        h = col_to_human.get(col)
        if h is not None:
            last_h = h
        human_x.append(last_h)
    human_x = np.array(human_x, dtype=float)

    window = 20
    kernel = np.ones(window) / window
    smoothed = np.convolve(scores, kernel, mode="same")

    seq_map = {rec.id: rec for rec in alignment}
    ordered = [seq_map[sp] for sp in SPECIES_ORDER if sp in seq_map]
    icons = {sp: load_icon(sp) for sp in SPECIES_ORDER}

    species_lengths = []
    for rec in ordered:
        sl = sum(1 for c in range(n_cols) if rec.seq[c] != "-")
        species_lengths.append(sl)

    n = len(ordered)
    x_max = max(human_len, max(species_lengths)) * 1.03

    # ── Figure with shared x-axis ──
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(18, 15), sharex=True,
        gridspec_kw={"height_ratios": [1, 2.2], "hspace": 0.05}
    )

    # ── Panel A: Conservation ──
    ax_top.fill_between(human_x, smoothed, alpha=0.3, color="#2196F3")
    ax_top.plot(human_x, smoothed, linewidth=0.8, color="#0D47A1")
    ax_top.axhline(0.9, ls="--", color="red", alpha=0.5, lw=1)
    ax_top.axhline(0.7, ls="--", color="orange", alpha=0.5, lw=1)
    ax_top.set_ylim(0, 1.05)
    ax_top.set_ylabel("Conservation", fontsize=12)
    ax_top.set_title("UBR3 Conservation and Domain Architecture",
                     fontsize=16, fontweight="bold", pad=12)

    for dname, h_start, h_end, color in DOMAINS:
        ax_top.axvspan(h_start, h_end, alpha=0.15, color=color)
        ax_bot.axvspan(h_start, h_end, alpha=0.08, color=color)

    domain_patches = [
        mpatches.Patch(facecolor=c, edgecolor="none", label=nm, alpha=0.5)
        for nm, _, _, c in DOMAINS
    ]
    ax_top.legend(handles=domain_patches, fontsize=8, loc="lower right",
                  framealpha=0.9)
    ax_top.spines["top"].set_visible(False)
    ax_top.spines["right"].set_visible(False)
    ax_top.text(0.01, 0.95, "A", transform=ax_top.transAxes,
                fontsize=22, fontweight="bold", va="top")

    # ── Panel B: Domain architecture ──
    row_h = 1.0
    bar_h = 0.4

    for i, rec in enumerate(ordered):
        sp_id = rec.id
        y = (n - 1 - i) * row_h
        seq_len = species_lengths[i]
        ph_color = PHYLUM_COLORS.get(sp_id, "#888")

        ax_bot.barh(y, seq_len, height=bar_h,
                    color="#F5F5F5", edgecolor="#CCCCCC", linewidth=0.4)
        ax_bot.barh(y, seq_len, height=bar_h,
                    color=ph_color, alpha=0.08, edgecolor="none")

        for dname, h_start, h_end, color in DOMAINS:
            col_s = human_to_col.get(h_start)
            col_e = human_to_col.get(h_end)
            if col_s is None or col_e is None:
                continue
            sp_start = sp_end = None
            res_count = 0
            for c in range(n_cols):
                if rec.seq[c] != "-":
                    res_count += 1
                    if c >= col_s and sp_start is None:
                        sp_start = res_count
                    if c <= col_e:
                        sp_end = res_count
            if sp_start and sp_end and sp_end > sp_start:
                w = sp_end - sp_start
                ax_bot.barh(y, w, left=sp_start, height=bar_h,
                            color=color, edgecolor="white", linewidth=0.4,
                            alpha=0.88)

        display = DISPLAY_NAMES.get(sp_id, sp_id)
        common = COMMON_NAMES.get(sp_id, "")
        is_human = "sapiens" in sp_id.lower()
        label = f"  {display}  ({common})"

        ax_bot.annotate(
            label, xy=(0, y), xytext=(-15, 0),
            textcoords="offset points",
            va="center", ha="right",
            fontsize=9.5, fontstyle="italic",
            fontweight="bold" if is_human else "normal",
            color=ph_color, annotation_clip=False,
        )

        icon = icons.get(sp_id)
        if icon is not None:
            im = OffsetImage(icon, zoom=0.042)
            im.image.axes = ax_bot
            ab = AnnotationBbox(
                im, (0, y),
                xybox=(-135, 0), boxcoords="offset points",
                frameon=False, annotation_clip=False,
            )
            ax_bot.add_artist(ab)

        ax_bot.text(seq_len + 15, y, f"{seq_len} aa",
                    va="center", ha="left", fontsize=7, color="#999")

    ax_bot.set_yticks([])
    ax_bot.set_ylim(-0.8, (n - 1) * row_h + 0.8)
    ax_bot.set_xlim(0, x_max)
    ax_bot.set_xlabel("Amino acid position (human numbering)", fontsize=12,
                      labelpad=10)
    ax_bot.spines["top"].set_visible(False)
    ax_bot.spines["right"].set_visible(False)
    ax_bot.spines["left"].set_visible(False)
    ax_bot.text(0.01, 0.99, "B", transform=ax_bot.transAxes,
                fontsize=22, fontweight="bold", va="top")

    fig.subplots_adjust(left=0.18, right=0.95)

    out_path = os.path.join(OUTPUT_DIR,
                            "ubr3_combined_conservation_architecture.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
