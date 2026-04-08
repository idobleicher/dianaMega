"""
Publication-style RPL22 / RPL22L1 alignment figures:
  - ClustalW-style text alignment for the full protein (L22e domain)
  - Domain architecture diagram across species
  - Combined multi-panel figure
"""
import os
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

from Bio import AlignIO

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

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

STRONG_GROUPS = [
    set("STA"), set("NEQK"), set("NHQK"), set("NDEQ"), set("QHRK"),
    set("MILV"), set("MILF"), set("HY"), set("FYW"),
]
WEAK_GROUPS = [
    set("CSA"), set("ATV"), set("SAG"), set("STNK"), set("STPA"),
    set("SGND"), set("SNDEQK"), set("NDEQHK"), set("NEQHRK"),
    set("FVLIM"), set("HFY"),
]

AA_BG = {
    "D": "#FFCDD2", "E": "#FFCDD2",
    "R": "#BBDEFB", "K": "#BBDEFB", "H": "#C5CAE9",
    "S": "#C8E6C9", "T": "#C8E6C9", "N": "#DCEDC8", "Q": "#DCEDC8",
    "A": "#FFE0B2", "V": "#FFE0B2", "I": "#FFE0B2", "L": "#FFE0B2",
    "M": "#FFCCBC", "F": "#FFCCBC", "W": "#FFCCBC", "Y": "#FFCCBC",
    "P": "#E1BEE7", "G": "#CFD8DC", "C": "#FCE4EC",
}
AA_FG = {
    "D": "#B71C1C", "E": "#B71C1C",
    "R": "#0D47A1", "K": "#0D47A1", "H": "#1A237E",
    "S": "#1B5E20", "T": "#1B5E20", "N": "#33691E", "Q": "#33691E",
    "A": "#E65100", "V": "#E65100", "I": "#E65100", "L": "#E65100",
    "M": "#BF360C", "F": "#BF360C", "W": "#BF360C", "Y": "#BF360C",
    "P": "#4A148C", "G": "#37474F", "C": "#880E4F",
}


def get_label(rec_id):
    dn = DISPLAY_NAMES.get(rec_id, rec_id)
    cn = COMMON_NAMES.get(rec_id, "")
    return f"{dn} ({cn})" if cn else dn


def alignment_col_to_human_pos(alignment):
    human_idx = 0
    for i, rec in enumerate(alignment):
        if "sapiens" in rec.id.lower() or "human" in rec.id.lower():
            human_idx = i
            break
    col_to_human = {}
    res_num = 0
    for col in range(alignment.get_alignment_length()):
        aa = alignment[human_idx, col]
        if aa != "-":
            res_num += 1
            col_to_human[col] = res_num
        else:
            col_to_human[col] = None
    human_to_col = {v: k for k, v in col_to_human.items() if v is not None}
    return col_to_human, human_to_col, human_idx


def conservation_symbol(alignment, col):
    residues = [alignment[row, col] for row in range(len(alignment))]
    non_gap = [r for r in residues if r != "-"]
    if len(non_gap) < 2:
        return " "
    unique = set(non_gap)
    if len(unique) == 1:
        return "*"
    for g in STRONG_GROUPS:
        if unique.issubset(g):
            return ":"
    for g in WEAK_GROUPS:
        if unique.issubset(g):
            return "."
    return " "


def get_alignment_slice(human_to_col, human_start, human_end):
    col_start = human_to_col.get(human_start)
    col_end = human_to_col.get(human_end)
    if col_start is None or col_end is None:
        for h in range(human_start, human_end + 1):
            if h in human_to_col:
                if col_start is None:
                    col_start = human_to_col[h]
                col_end = human_to_col[h]
    return col_start, col_end


def render_clustal_figure(alignment, col_start, col_end, col_to_human,
                           title, filename, block_width=60,
                           highlight_positions=None):
    n_seqs = len(alignment)
    labels = [get_label(rec.id) for rec in alignment]
    max_label = max(len(l) for l in labels)
    cols = list(range(col_start, col_end + 1))
    n_blocks = (len(cols) + block_width - 1) // block_width
    extra_per_block = 2.0 if highlight_positions else 0
    lines_per_block = n_seqs + 2 + extra_per_block

    total_lines = n_blocks * lines_per_block + 3
    char_w_pts = 6.2
    label_space = max_label * 0.52
    fig_w = max(10, (block_width + label_space + 2) * char_w_pts / 72 + 0.5)
    fig_h = max(3, total_lines * 0.155 + 0.8)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, fig_w)
    ax.set_ylim(0, total_lines + 3)
    ax.invert_yaxis()
    ax.axis("off")

    font = {"family": "monospace", "fontsize": 6}
    bold = {"family": "monospace", "fontsize": 6, "fontweight": "bold"}

    char_w = fig_w / (block_width + label_space + 2)
    x_seq = (label_space + 0.1) * char_w
    y = 0.0

    ax.text(fig_w / 2, y, title, ha="center", va="top",
            fontsize=10, fontweight="bold", fontstyle="italic",
            family="sans-serif")
    y += 1.8 if not highlight_positions else 2.8

    for blk in range(n_blocks):
        chunk_cols = cols[blk * block_width:(blk + 1) * block_width]

        first_h = last_h = None
        for c in chunk_cols:
            h = col_to_human.get(c)
            if h is not None:
                if first_h is None:
                    first_h = h
                last_h = h

        if highlight_positions:
            hits_top = []
            for ci, c in enumerate(chunk_cols):
                h = col_to_human.get(c)
                if h is not None and h in highlight_positions:
                    hits_top.append((ci, h))
            if hits_top:
                for ci, h in hits_top:
                    x = x_seq + ci * char_w
                    ax.plot(x, y + 0.55, "v", color="#C62828",
                            markersize=2.5, zorder=10, clip_on=False)
                    ax.text(x, y + 0.35, highlight_positions[h],
                            fontsize=4.2, ha="left", va="bottom",
                            color="#B71C1C", fontweight="bold",
                            rotation=45, rotation_mode="anchor",
                            clip_on=False)
                y += 1.2

        num_y = y + 0.0
        if first_h:
            ax.text(x_seq - char_w * 0.8, num_y, str(first_h),
                    fontsize=5, family="monospace", color="#999",
                    va="center", ha="right")
        if last_h:
            ax.text(x_seq + len(chunk_cols) * char_w + char_w * 0.15,
                    num_y, str(last_h),
                    fontsize=5, family="monospace", color="#999",
                    va="center", ha="left")

        block_top = y - 0.5
        block_bot = y + n_seqs * 1 + 0.5
        if highlight_positions:
            for ci, c in enumerate(chunk_cols):
                h = col_to_human.get(c)
                if h is not None and h in highlight_positions:
                    x = x_seq + ci * char_w
                    ax.add_patch(plt.Rectangle(
                        (x - char_w * 0.45, block_top),
                        char_w * 0.9, block_bot - block_top,
                        fc="#FFCDD2", ec="none", alpha=0.35, zorder=-1))

        for seq_i in range(n_seqs):
            rec_id = alignment[seq_i].id
            label = labels[seq_i]
            is_human = "sapiens" in rec_id.lower()
            lbl_style = bold if is_human else font
            phylum_c = PHYLUM_COLORS.get(rec_id, "#333")
            ax.text(0.1, y, label, va="center", color=phylum_c, **lbl_style)

            for ci, c in enumerate(chunk_cols):
                aa = alignment[seq_i, c]
                x = x_seq + ci * char_w
                if aa == "-":
                    ax.text(x, y, "-", va="center", ha="center",
                            color="#DDD", **font)
                    continue

                residues = [alignment[r, c] for r in range(n_seqs)]
                non_gap = [r for r in residues if r != "-"]
                is_identical = len(set(non_gap)) == 1 and len(non_gap) > 1
                sym = conservation_symbol(alignment, c)

                if is_identical:
                    bg = AA_BG.get(aa, "#E0E0E0")
                    ax.add_patch(plt.Rectangle(
                        (x - char_w * 0.45, y - 0.42),
                        char_w * 0.9, 0.84,
                        fc=bg, ec="none", zorder=0))
                    ax.text(x, y, aa, va="center", ha="center",
                            color=AA_FG.get(aa, "#333"), **bold)
                elif sym == ":":
                    bg = AA_BG.get(aa, "#F5F5F5")
                    ax.add_patch(plt.Rectangle(
                        (x - char_w * 0.45, y - 0.42),
                        char_w * 0.9, 0.84,
                        fc=bg, ec="none", alpha=0.4, zorder=0))
                    ax.text(x, y, aa, va="center", ha="center",
                            color=AA_FG.get(aa, "#555"), **font)
                else:
                    ax.text(x, y, aa, va="center", ha="center",
                            color="#666", **font)
            y += 1

        for ci, c in enumerate(chunk_cols):
            sym = conservation_symbol(alignment, c)
            x = x_seq + ci * char_w
            if sym == "*":
                ax.text(x, y, sym, va="center", ha="center",
                        color="#000", fontweight="bold", **font)
            elif sym == ":":
                ax.text(x, y, sym, va="center", ha="center",
                        color="#555", **font)
            elif sym == ".":
                ax.text(x, y, sym, va="center", ha="center",
                        color="#AAA", **font)
        y += 0.7
        y += 1.3

    fig.tight_layout(pad=0.3)
    path = os.path.join(OUTPUT_DIR, filename)
    fig.savefig(path, dpi=300, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"  Saved {filename}")
    return path


def plot_domain_architecture(alignment, col_to_human, human_to_col,
                              human_idx, protein_name, prefix,
                              l22e_start, l22e_end):
    n_seqs = len(alignment)
    labels = [get_label(rec.id) for rec in alignment]

    domains = [
        (f"L22e domain ({l22e_start}-{l22e_end})", l22e_start, l22e_end, "#7E57C2"),
    ]

    species_lengths = []
    for seq_i in range(n_seqs):
        seq_len = sum(1 for c in range(alignment.get_alignment_length())
                      if alignment[seq_i, c] != "-")
        species_lengths.append(seq_len)
    max_len = max(species_lengths)

    fig_h = max(5, n_seqs * 0.55 + 2.5)
    fig_w = max(10, (60 + max(len(l) for l in labels) + 6) * 6.2 / 72 + 0.5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    bar_h = 0.35
    y_spacing = 0.7

    for seq_i in range(n_seqs):
        y = (n_seqs - 1 - seq_i) * y_spacing
        seq_len = species_lengths[seq_i]
        rec_id = alignment[seq_i].id
        ph_color = PHYLUM_COLORS.get(rec_id, "#888")

        ax.barh(y, seq_len, height=bar_h,
                color="#EEEEEE", edgecolor="#CCCCCC", linewidth=0.4)
        ax.barh(y, seq_len, height=bar_h,
                color=ph_color, alpha=0.08, edgecolor="none")

        for dname, h_start, h_end, color in domains:
            col_s = human_to_col.get(h_start)
            col_e = human_to_col.get(h_end)
            if col_s is None or col_e is None:
                continue
            sp_start = sp_end = None
            res_count = 0
            for c in range(alignment.get_alignment_length()):
                if alignment[seq_i, c] != "-":
                    res_count += 1
                    if c >= col_s and sp_start is None:
                        sp_start = res_count
                    if c <= col_e:
                        sp_end = res_count
            if sp_start and sp_end and sp_end > sp_start:
                width = sp_end - sp_start
                ax.barh(y, width, left=sp_start, height=bar_h,
                        color=color, edgecolor="white", linewidth=0.4,
                        alpha=0.85)

    yticks = [(n_seqs - 1 - i) * y_spacing for i in range(n_seqs)]
    ax.set_yticks(yticks)
    ax.set_yticklabels(labels, fontsize=7.5)

    for i, rec in enumerate(alignment):
        rec_id = rec.id
        ph_color = PHYLUM_COLORS.get(rec_id, "#333")
        is_human = "sapiens" in rec_id.lower()
        ax.get_yticklabels()[i].set_color(ph_color)
        if is_human:
            ax.get_yticklabels()[i].set_fontweight("bold")

    ax.set_xlabel("Amino acid position", fontsize=10)
    ax.set_xlim(0, max_len + 10)
    ax.set_title(f"{protein_name} Domain Architecture across Species",
                 fontsize=12, fontweight="bold", pad=12)

    legend_patches = [
        mpatches.Patch(facecolor=c, edgecolor="white", label=n, alpha=0.85)
        for n, _, _, c in domains
    ]
    ax.legend(handles=legend_patches, fontsize=8, loc="upper right",
              framealpha=0.95, edgecolor="#DDD")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, f"{prefix}_domain_architecture.png")
    fig.savefig(path, dpi=300, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"  Saved {os.path.basename(path)}")
    return path


def create_combined_figure(paths_a, path_b, prefix):
    from PIL import Image, ImageDraw, ImageFont

    imgs_a = [Image.open(p) for p in paths_a]
    img_b = Image.open(path_b)

    target_w = imgs_a[0].width

    def resize_w(img, tw):
        if img.width == tw:
            return img
        r = tw / img.width
        return img.resize((tw, int(img.height * r)), Image.LANCZOS)

    imgs_a = [resize_w(img, target_w) for img in imgs_a]
    img_b = resize_w(img_b, target_w)

    label_h = 50
    gap = 15
    total_h = (sum(img.height for img in imgs_a) +
               img_b.height + label_h * 2 + gap * 3)

    combined = Image.new("RGB", (target_w, total_h), "white")
    draw = ImageDraw.Draw(combined)
    try:
        font = ImageFont.truetype("arialbd.ttf", 32)
    except:
        try:
            font = ImageFont.truetype("arial.ttf", 32)
        except:
            font = ImageFont.load_default()

    y = gap
    draw.text((12, y), "A", fill="black", font=font)
    y += label_h
    for img in imgs_a:
        combined.paste(img, (0, y))
        y += img.height

    y += gap
    draw.text((12, y), "B", fill="black", font=font)
    y += label_h
    combined.paste(img_b, (0, y))

    path = os.path.join(OUTPUT_DIR, f"{prefix}_combined_panels.png")
    combined.save(path, dpi=(300, 300))
    print(f"  Saved {os.path.basename(path)}")
    return path


def process_protein(fasta_path, protein_name, prefix, l22e_start, l22e_end):
    print(f"\n{'='*60}")
    print(f"{protein_name} Publication Figures")
    print(f"{'='*60}")

    alignment = AlignIO.read(fasta_path, "fasta")
    col_to_human, human_to_col, human_idx = alignment_col_to_human_pos(alignment)
    human_len = max(v for v in col_to_human.values() if v is not None)
    print(f"  {len(alignment)} sequences x {alignment.get_alignment_length()} columns")
    print(f"  Human length: {human_len} aa")

    panel_a_paths = []

    print(f"\n--- Full Protein Alignment ---")
    col_start, col_end = get_alignment_slice(human_to_col, 1, human_len)
    if col_start is not None and col_end is not None:
        path = render_clustal_figure(
            alignment, col_start, col_end, col_to_human,
            f"{protein_name} Full Protein Alignment (aa 1-{human_len})",
            f"{prefix}_full_alignment.png", block_width=60)
        panel_a_paths.append(path)

    print(f"\n--- L22e Domain Alignment ---")
    col_start, col_end = get_alignment_slice(human_to_col, l22e_start, l22e_end)
    if col_start is not None and col_end is not None:
        path = render_clustal_figure(
            alignment, col_start, col_end, col_to_human,
            f"{protein_name} L22e Domain (aa {l22e_start}-{l22e_end})",
            f"{prefix}_l22e_domain.png", block_width=60)
        panel_a_paths.append(path)

    print(f"\n--- Domain Architecture ---")
    path_b = plot_domain_architecture(
        alignment, col_to_human, human_to_col, human_idx,
        protein_name, prefix, l22e_start, l22e_end)

    if panel_a_paths and path_b:
        print(f"\n--- Combined Figure ---")
        try:
            create_combined_figure(panel_a_paths, path_b, prefix)
        except Exception as exc:
            print(f"  Combined figure failed: {exc}")

    return alignment, col_to_human, human_to_col, human_idx


def main():
    rpl22_fasta = os.path.join(OUTPUT_DIR, "rpl22_aligned.fasta")
    rpl22l1_fasta = os.path.join(OUTPUT_DIR, "rpl22l1_aligned.fasta")

    if os.path.exists(rpl22_fasta):
        process_protein(rpl22_fasta, "RPL22", "rpl22_pub",
                         l22e_start=15, l22e_end=128)
    else:
        print(f"SKIP RPL22: {rpl22_fasta} not found")

    if os.path.exists(rpl22l1_fasta):
        process_protein(rpl22l1_fasta, "RPL22L1", "rpl22l1_pub",
                         l22e_start=16, l22e_end=122)
    else:
        print(f"SKIP RPL22L1: {rpl22l1_fasta} not found")

    print(f"\n{'='*60}")
    print("DONE")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
