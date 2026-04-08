"""
Generate domain-architecture figures with animal icons for RPL22 and RPL22L1.
Reuses icons from the UBR3 analysis.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.image as mpimg
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from Bio import AlignIO

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
ICONS_DIR = os.path.join(os.path.dirname(OUTPUT_DIR), "ubr3_alignment", "icons")

SPECIES_ORDER = [
    "H_sapiens", "M_musculus", "X_laevis", "G_gallus", "D_rerio",
    "P_marinus", "C_intestinalis", "B_floridae", "S_purpuratus",
    "O_sinensis", "D_melanogaster", "A_mellifera", "P_tepidariorum",
    "H_vulgaris", "N_vectensis",
]

ICON_FILES = {
    "H_sapiens": "icon_human.png", "M_musculus": "icon_mouse.png",
    "X_laevis": "icon_frog.png", "G_gallus": "icon_chicken.png",
    "D_rerio": "icon_zebrafish.png", "P_marinus": "icon_lamprey.png",
    "C_intestinalis": "icon_seasquirt.png", "B_floridae": "icon_amphioxus.png",
    "S_purpuratus": "icon_seaurchin.png", "O_sinensis": "icon_octopus.png",
    "D_melanogaster": "icon_fruitfly.png", "A_mellifera": "icon_honeybee.png",
    "P_tepidariorum": "icon_spider.png", "H_vulgaris": "icon_hydra.png",
    "N_vectensis": "icon_seaanemone.png",
}

DISPLAY_NAMES = {
    "H_sapiens": "H. sapiens", "M_musculus": "M. musculus",
    "X_laevis": "X. laevis", "G_gallus": "G. gallus",
    "D_rerio": "D. rerio", "P_marinus": "P. marinus",
    "C_intestinalis": "C. intestinalis", "B_floridae": "B. floridae",
    "S_purpuratus": "S. purpuratus", "O_sinensis": "O. sinensis",
    "D_melanogaster": "D. melanogaster", "A_mellifera": "A. mellifera",
    "P_tepidariorum": "P. tepidariorum", "H_vulgaris": "H. vulgaris",
    "N_vectensis": "N. vectensis",
}
COMMON_NAMES = {
    "H_sapiens": "Human", "M_musculus": "Mouse", "X_laevis": "Frog",
    "G_gallus": "Chicken", "D_rerio": "Zebrafish",
    "P_marinus": "Lamprey", "C_intestinalis": "Sea squirt",
    "B_floridae": "Amphioxus", "S_purpuratus": "Sea urchin",
    "O_sinensis": "Octopus", "D_melanogaster": "Fruit fly",
    "A_mellifera": "Honeybee", "P_tepidariorum": "Spider",
    "H_vulgaris": "Hydra", "N_vectensis": "Sea anemone",
}
PHYLUM_LABELS = {
    "H_sapiens": "Mammalia", "M_musculus": "Mammalia",
    "X_laevis": "Amphibia", "G_gallus": "Aves",
    "D_rerio": "Actinopterygii", "P_marinus": "Cyclostomata",
    "C_intestinalis": "Tunicata", "B_floridae": "Cephalochordata",
    "S_purpuratus": "Echinodermata", "O_sinensis": "Mollusca",
    "D_melanogaster": "Insecta", "A_mellifera": "Insecta",
    "P_tepidariorum": "Arachnida", "H_vulgaris": "Cnidaria",
    "N_vectensis": "Cnidaria",
}
PHYLUM_COLORS = {
    "H_sapiens": "#B71C1C", "M_musculus": "#B71C1C",
    "X_laevis": "#BF360C", "G_gallus": "#D84315",
    "D_rerio": "#E65100", "P_marinus": "#E65100",
    "C_intestinalis": "#F57F17", "B_floridae": "#827717",
    "S_purpuratus": "#33691E", "O_sinensis": "#00695C",
    "D_melanogaster": "#01579B", "A_mellifera": "#01579B",
    "P_tepidariorum": "#0D47A1", "H_vulgaris": "#4A148C",
    "N_vectensis": "#4A148C",
}


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
        print(f"  Warning: icon not found: {path}")
        return None
    img = mpimg.imread(path)
    if img.max() <= 1.0:
        img = (img * 255).astype(np.uint8)
    return auto_crop(img)


def alignment_col_to_human_pos(alignment):
    human_idx = 0
    for i, rec in enumerate(alignment):
        if "sapiens" in rec.id.lower():
            human_idx = i
            break
    col_to_human = {}
    res_num = 0
    for col in range(alignment.get_alignment_length()):
        if alignment[human_idx, col] != "-":
            res_num += 1
            col_to_human[col] = res_num
        else:
            col_to_human[col] = None
    human_to_col = {v: k for k, v in col_to_human.items() if v is not None}
    return col_to_human, human_to_col, human_idx


def plot_animal_architecture(fasta_path, protein_name, output_filename,
                              l22e_start, l22e_end):
    alignment = AlignIO.read(fasta_path, "fasta")
    col_to_human, human_to_col, human_idx = alignment_col_to_human_pos(alignment)

    seq_map = {rec.id: rec for rec in alignment}
    ordered = [seq_map[sp] for sp in SPECIES_ORDER if sp in seq_map]
    icons = {sp: load_icon(sp) for sp in SPECIES_ORDER}

    domains = [
        (f"L22e domain ({l22e_start}-{l22e_end})", l22e_start, l22e_end, "#7E57C2"),
    ]

    species_lengths = []
    for rec in ordered:
        seq_len = sum(1 for c in range(alignment.get_alignment_length())
                      if rec.seq[c] != "-")
        species_lengths.append(seq_len)
    max_len = max(species_lengths)

    n = len(ordered)
    row_h = 1.15
    fig_h = max(12, n * row_h + 3.5)
    fig, ax = plt.subplots(figsize=(16, fig_h))

    bar_h = 0.32
    icon_zoom = 0.085
    label_zone = 4.2
    icon_center = 0.6
    label_left = 1.2
    bar_left = label_zone
    bar_right_max = 16.0 - 0.3
    bar_width_avail = bar_right_max - bar_left
    px_per_aa = bar_width_avail / max_len

    for i, rec in enumerate(ordered):
        sp_id = rec.id
        y = (n - 1 - i) * row_h
        seq_len = species_lengths[i]
        ph_color = PHYLUM_COLORS.get(sp_id, "#888")

        bar_len = seq_len * px_per_aa
        ax.barh(y, bar_len, left=bar_left, height=bar_h,
                color="#F5F5F5", edgecolor="#CCCCCC", linewidth=0.5)
        ax.barh(y, bar_len, left=bar_left, height=bar_h,
                color=ph_color, alpha=0.10, edgecolor="none")

        for dname, h_start, h_end, color in domains:
            col_s = human_to_col.get(h_start)
            col_e = human_to_col.get(h_end)
            if col_s is None or col_e is None:
                continue
            sp_start = sp_end = None
            res_count = 0
            for c in range(alignment.get_alignment_length()):
                if rec.seq[c] != "-":
                    res_count += 1
                    if c >= col_s and sp_start is None:
                        sp_start = res_count
                    if c <= col_e:
                        sp_end = res_count
            if sp_start and sp_end and sp_end > sp_start:
                w = (sp_end - sp_start) * px_per_aa
                ax.barh(y, w, left=bar_left + sp_start * px_per_aa,
                        height=bar_h, color=color, edgecolor="white",
                        linewidth=0.5, alpha=0.88)

        icon = icons.get(sp_id)
        if icon is not None:
            im = OffsetImage(icon, zoom=icon_zoom)
            im.image.axes = ax
            ab = AnnotationBbox(im, (icon_center, y), frameon=False,
                                box_alignment=(0.5, 0.5))
            ax.add_artist(ab)

        display = DISPLAY_NAMES.get(sp_id, sp_id)
        common = COMMON_NAMES.get(sp_id, "")
        phylum = PHYLUM_LABELS.get(sp_id, "")
        is_human = "sapiens" in sp_id.lower()

        ax.text(label_left, y + 0.18, f"{display}  ({common})",
                va="center", ha="left", fontsize=12, fontstyle="italic",
                fontweight="bold" if is_human else "normal",
                color=ph_color, clip_on=False)
        ax.text(label_left, y - 0.20, phylum,
                va="center", ha="left", fontsize=9, color="#999",
                clip_on=False)
        ax.text(bar_left + bar_len + 0.06, y, f"{seq_len} aa",
                va="center", ha="left", fontsize=8, color="#999")

    yticks = [(n - 1 - i) * row_h for i in range(n)]
    ax.set_yticks(yticks)
    ax.set_yticklabels([""] * n)
    ax.tick_params(axis="y", length=0)

    aa_ticks = list(range(0, max_len + 1, 20))
    if aa_ticks[-1] < max_len:
        aa_ticks.append(max_len)
    ax.set_xticks([bar_left + t * px_per_aa for t in aa_ticks])
    ax.set_xticklabels([str(t) for t in aa_ticks], fontsize=10)

    ax.set_xlim(-0.1, bar_right_max + 0.5)
    ax.set_ylim(-0.8, (n - 1) * row_h + 1.0)
    ax.set_xlabel("Amino acid position", fontsize=12, labelpad=10)
    ax.set_title(f"{protein_name} Protein Domain Architecture Across the Animal Kingdom",
                 fontsize=15, fontweight="bold", pad=18)

    legend_patches = [
        mpatches.Patch(facecolor=c, edgecolor="white", label=nm, alpha=0.88)
        for nm, _, _, c in domains
    ]
    ax.legend(handles=legend_patches, fontsize=10, loc="upper right",
              framealpha=0.95, edgecolor="#DDD", borderpad=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_bounds(bar_left, bar_left + max_len * px_per_aa)

    fig.tight_layout()
    out_path = os.path.join(OUTPUT_DIR, output_filename)
    fig.savefig(out_path, dpi=300, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"Saved: {out_path}")
    return out_path


def main():
    print("=" * 60)
    print("RPL22 / RPL22L1 Domain Architecture with Animal Icons")
    print("=" * 60)

    rpl22_fasta = os.path.join(OUTPUT_DIR, "rpl22_aligned.fasta")
    rpl22l1_fasta = os.path.join(OUTPUT_DIR, "rpl22l1_aligned.fasta")

    if os.path.exists(rpl22_fasta):
        plot_animal_architecture(rpl22_fasta, "RPL22",
                                  "rpl22_animal_domain_architecture.png",
                                  l22e_start=15, l22e_end=128)
    else:
        print(f"SKIP RPL22: {rpl22_fasta} not found")

    if os.path.exists(rpl22l1_fasta):
        plot_animal_architecture(rpl22l1_fasta, "RPL22L1",
                                  "rpl22l1_animal_domain_architecture.png",
                                  l22e_start=16, l22e_end=122)
    else:
        print(f"SKIP RPL22L1: {rpl22l1_fasta} not found")

    print(f"\n{'='*60}")
    print("DONE")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
