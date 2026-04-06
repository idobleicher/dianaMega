"""
Publication-style UBR3 alignment figures with broad phylogenetic panel.
  Panel A: ClustalW-style text alignment for UBR box, 400-600 region, RING domain
  Panel B: Domain architecture diagram across species
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
ALIGNED_FASTA = os.path.join(OUTPUT_DIR, "ubr3_aligned.fasta")

DISPLAY_NAMES = {
    "H_sapiens":       "H. sapiens",
    "M_musculus":      "M. musculus",
    "X_laevis":        "X. laevis",
    "P_marinus":       "P. marinus",
    "C_intestinalis":  "C. intestinalis",
    "B_floridae":      "B. floridae",
    "S_purpuratus":    "S. purpuratus",
    "O_sinensis":      "O. sinensis",
    "D_melanogaster":  "D. melanogaster",
    "A_mellifera":     "A. mellifera",
    "P_tepidariorum":  "P. tepidariorum",
    "H_vulgaris":      "H. vulgaris",
    "N_vectensis":     "N. vectensis",
}

COMMON_NAMES = {
    "H_sapiens":       "Human",
    "M_musculus":      "Mouse",
    "X_laevis":        "Frog",
    "P_marinus":       "Lamprey",
    "C_intestinalis":  "Sea squirt",
    "B_floridae":      "Amphioxus",
    "S_purpuratus":    "Sea urchin",
    "O_sinensis":      "Octopus",
    "D_melanogaster":  "Fruit fly",
    "A_mellifera":     "Honeybee",
    "P_tepidariorum":  "Spider",
    "H_vulgaris":      "Hydra",
    "N_vectensis":     "Sea anemone",
}

PHYLUM_COLORS = {
    "H_sapiens":       "#B71C1C",
    "M_musculus":      "#B71C1C",
    "X_laevis":        "#BF360C",
    "P_marinus":       "#E65100",
    "C_intestinalis":  "#F57F17",
    "B_floridae":      "#827717",
    "S_purpuratus":    "#33691E",
    "O_sinensis":      "#00695C",
    "D_melanogaster":  "#01579B",
    "A_mellifera":     "#01579B",
    "P_tepidariorum":  "#0D47A1",
    "H_vulgaris":      "#4A148C",
    "N_vectensis":     "#4A148C",
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


def load_alignment():
    return AlignIO.read(ALIGNED_FASTA, "fasta")


def get_display_name(rec_id):
    return DISPLAY_NAMES.get(rec_id, rec_id)


def get_label(rec_id):
    dn = DISPLAY_NAMES.get(rec_id, rec_id)
    cn = COMMON_NAMES.get(rec_id, "")
    if cn:
        return f"{dn} ({cn})"
    return dn


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


# ── Panel A: ClustalW-style alignment figure ─────────────────────────────
def render_clustal_figure(alignment, col_start, col_end, col_to_human,
                           title, filename, block_width=60):
    n_seqs = len(alignment)
    labels = [get_label(rec.id) for rec in alignment]
    max_label = max(len(l) for l in labels) + 1
    cols = list(range(col_start, col_end + 1))
    n_blocks = (len(cols) + block_width - 1) // block_width
    lines_per_block = n_seqs + 2

    total_lines = n_blocks * lines_per_block + 4
    char_w_pts = 6.2
    fig_w = max(10, (block_width + max_label + 6) * char_w_pts / 72 + 0.5)
    fig_h = max(3, total_lines * 0.155 + 0.8)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, fig_w)
    ax.set_ylim(0, total_lines + 1)
    ax.invert_yaxis()
    ax.axis("off")

    font = {"family": "monospace", "fontsize": 6}
    bold = {"family": "monospace", "fontsize": 6, "fontweight": "bold"}

    char_w = fig_w / (block_width + max_label + 6)
    x_seq = (max_label + 1) * char_w
    y = 0.5

    ax.text(fig_w / 2, y, title, ha="center", va="top",
            fontsize=10, fontweight="bold", fontstyle="italic",
            family="sans-serif")
    y += 1.8

    for blk in range(n_blocks):
        chunk_cols = cols[blk * block_width:(blk + 1) * block_width]

        first_h = last_h = None
        for c in chunk_cols:
            h = col_to_human.get(c)
            if h is not None:
                if first_h is None:
                    first_h = h
                last_h = h
        if first_h:
            ax.text(x_seq - char_w * 0.3, y - 0.25, str(first_h),
                    fontsize=5, family="monospace", color="#999", va="bottom")
        if last_h:
            ax.text(x_seq + len(chunk_cols) * char_w + char_w * 0.5,
                    y - 0.25, str(last_h),
                    fontsize=5, family="monospace", color="#999",
                    va="bottom", ha="left")

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
        y += 2

    fig.tight_layout(pad=0.3)
    path = os.path.join(OUTPUT_DIR, filename)
    fig.savefig(path, dpi=300, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"  Saved {filename}")
    return path


# ── Panel B: Domain architecture ─────────────────────────────────────────
def plot_domain_architecture(alignment, col_to_human, human_to_col,
                              human_idx, prefix):
    n_seqs = len(alignment)
    labels = [get_label(rec.id) for rec in alignment]

    domains = [
        ("UBR box (118-189)",     118,  189, "#7E57C2"),
        ("Region 400-600",        400,  600, "#1E88E5"),
        ("RING domain (1306-1364)", 1306, 1364, "#EF6C00"),
    ]

    species_lengths = []
    for seq_i in range(n_seqs):
        seq_len = sum(1 for c in range(alignment.get_alignment_length())
                      if alignment[seq_i, c] != "-")
        species_lengths.append(seq_len)
    max_len = max(species_lengths)

    fig_h = max(5, n_seqs * 0.55 + 2.5)
    fig, ax = plt.subplots(figsize=(13, fig_h))

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
    ax.set_xlim(0, max_len + 30)
    ax.set_title("UBR3 Domain Architecture across Species",
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


# ── Combined figure ──────────────────────────────────────────────────────
def create_combined_figure(paths_a, path_b, prefix):
    from PIL import Image, ImageDraw, ImageFont

    imgs_a = [Image.open(p) for p in paths_a]
    img_b = Image.open(path_b)

    max_w = max(img.width for img in imgs_a + [img_b])

    def resize_w(img, tw):
        r = tw / img.width
        return img.resize((tw, int(img.height * r)), Image.LANCZOS)

    imgs_a = [resize_w(img, max_w) for img in imgs_a]
    img_b = resize_w(img_b, max_w)

    label_h = 50
    gap = 15
    total_h = (sum(img.height for img in imgs_a) +
               img_b.height + label_h * 2 + gap * 3)

    combined = Image.new("RGB", (max_w, total_h), "white")
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


def main():
    print("=" * 60)
    print("UBR3 Publication-Style Figures (Broad Phylogeny)")
    print("=" * 60)

    alignment = load_alignment()
    col_to_human, human_to_col, human_idx = alignment_col_to_human_pos(alignment)
    print(f"  {len(alignment)} sequences x {alignment.get_alignment_length()} columns")
    print(f"  Species: {', '.join(get_display_name(r.id) for r in alignment)}")

    prefix = "ubr3_pub"

    regions = [
        ("UBR Box Domain (aa 118-189)", 118, 189, "ubr_box"),
        ("Conserved Region (aa 400-600)", 400, 600, "region_400_600"),
        ("RING Domain (aa 1306-1364)", 1306, 1364, "ring_domain"),
    ]

    panel_a_paths = []
    for title, h_start, h_end, suffix in regions:
        print(f"\n--- {title} ---")
        col_start, col_end = get_alignment_slice(human_to_col, h_start, h_end)
        if col_start is None or col_end is None:
            print(f"  Could not find alignment columns for {h_start}-{h_end}")
            continue
        print(f"  Alignment columns: {col_start}-{col_end}")
        path = render_clustal_figure(
            alignment, col_start, col_end, col_to_human,
            title, f"{prefix}_{suffix}.png", block_width=60)
        panel_a_paths.append(path)

    print(f"\n--- Domain Architecture ---")
    path_b = plot_domain_architecture(
        alignment, col_to_human, human_to_col, human_idx, prefix)

    if panel_a_paths and path_b:
        print(f"\n--- Combined Figure ---")
        try:
            create_combined_figure(panel_a_paths, path_b, prefix)
        except Exception as exc:
            print(f"  Combined figure failed: {exc}")

    print(f"\n{'='*60}")
    print("DONE")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
