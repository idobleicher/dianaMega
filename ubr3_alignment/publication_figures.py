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


SUBSTRATE_RESIDUES = {
    132: "Arg132", 140: "Pro140", 141: "Cys141",
    169: "Ala169", 171: "Asp171", 177: "Val177",
    444: "Val444", 448: "Val448", 449: "Gln449",
    511: "Trp511", 577: "Glu577",
}

RING_RESIDUES = {
    1239: "Cys1239", 1334: "Cys1334",
    1360: "Cys1360", 1363: "Cys1363",
}


# ── Panel A: ClustalW-style alignment figure ─────────────────────────────
def render_clustal_figure(alignment, col_start, col_end, col_to_human,
                           title, filename, block_width=60,
                           highlight_positions=None,
                           highlight_positions_2=None,
                           for_presentation=False):
    n_seqs = len(alignment)
    labels = [get_label(rec.id) for rec in alignment]
    max_label = max(len(l) for l in labels)
    cols = list(range(col_start, col_end + 1))
    len_cols = len(cols)

    if for_presentation:
        # Prefer wide rows (fewer stacked blocks); cap width for very long spans.
        block_width = max(72, min(96, len_cols))

    n_blocks = (len(cols) + block_width - 1) // block_width
    any_highlights = highlight_positions or highlight_positions_2
    n_highlight_rows = (1 if highlight_positions else 0) + (1 if highlight_positions_2 else 0)
    extra_per_block = 2.0 * n_highlight_rows if any_highlights else 0
    lines_per_block = n_seqs + 2 + extra_per_block

    total_lines = n_blocks * lines_per_block + 3

    if for_presentation:
        line_h = 1.32
        char_w_pts = 10.8
        fs_aa = 10.0
        fs_title = 13.5
        fs_coord = 8.2
        fs_ann = 7.8
        fs_cons_star = 12.5
        fs_cons_colon = 9.8
        fs_cons_dot = 9.0
        dpi = 450
        tri_ms = 4.8
        y_after_title = 2.05 if not any_highlights else 3.2
        y_hi = 1.55
        y_cons_gap = 0.92
        y_block_gap = 1.62
        box_half_h = 0.50 * line_h
    else:
        line_h = 1.0
        char_w_pts = 6.2
        fs_aa = 6
        fs_title = 10
        fs_coord = 5
        fs_ann = 4.2
        fs_cons_star = 6
        fs_cons_colon = 6
        fs_cons_dot = 6
        dpi = 300
        tri_ms = 2.5
        y_after_title = 1.8 if not any_highlights else 2.8
        y_hi = 1.2
        y_cons_gap = 0.7
        y_block_gap = 1.3
        box_half_h = 0.42

    label_space = max_label * 0.52
    fig_w = max(10, (block_width + label_space + 2) * char_w_pts / 72 + 0.5)
    line_unit = 0.155 * line_h
    fig_h = max(3, total_lines * line_unit + 0.8 + n_blocks * line_h)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, fig_w)
    ax.set_ylim(0, total_lines * line_h + 8)
    ax.invert_yaxis()
    ax.axis("off")

    font = {"family": "monospace", "fontsize": fs_aa}
    bold = {"family": "monospace", "fontsize": fs_aa, "fontweight": "bold"}

    char_w = fig_w / (block_width + label_space + 2)
    x_seq = (label_space + 0.1) * char_w
    y = 0.0

    ax.text(fig_w / 2, y, title, ha="center", va="top",
            fontsize=fs_title, fontweight="bold", fontstyle="italic",
            family="sans-serif")
    y += y_after_title

    box_h_aa = box_half_h * 2

    for blk in range(n_blocks):
        chunk_cols = cols[blk * block_width:(blk + 1) * block_width]

        first_h = last_h = None
        for c in chunk_cols:
            h = col_to_human.get(c)
            if h is not None:
                if first_h is None:
                    first_h = h
                last_h = h
        for hp, tri_c, txt_c in [
            (highlight_positions, "#C62828", "#B71C1C"),
            (highlight_positions_2, "#1565C0", "#0D47A1"),
        ]:
            if hp:
                hits_top = []
                for ci, c in enumerate(chunk_cols):
                    h = col_to_human.get(c)
                    if h is not None and h in hp:
                        hits_top.append((ci, h))
                if hits_top:
                    for ci, h in hits_top:
                        x = x_seq + ci * char_w
                        ax.plot(x, y + 0.55 * line_h, "v", color=tri_c,
                                markersize=tri_ms, zorder=10, clip_on=False)
                        ax.text(x, y + 0.35 * line_h, hp[h],
                                fontsize=fs_ann, ha="left", va="bottom",
                                color=txt_c, fontweight="bold",
                                rotation=45, rotation_mode="anchor",
                                clip_on=False)
                    y += y_hi * line_h

        num_y = y + 0.0
        if first_h:
            ax.text(x_seq - char_w * 0.8, num_y, str(first_h),
                    fontsize=fs_coord, family="monospace", color="#999",
                    va="center", ha="right")
        if last_h:
            ax.text(x_seq + len(chunk_cols) * char_w + char_w * 0.15,
                    num_y, str(last_h),
                    fontsize=fs_coord, family="monospace", color="#999",
                    va="center", ha="left")

        block_pad = 0.5 * line_h
        block_top = y - block_pad
        block_bot = y + n_seqs * line_h + block_pad
        for hp, col_bg in [
            (highlight_positions, "#FFCDD2"),
            (highlight_positions_2, "#BBDEFB"),
        ]:
            if hp:
                for ci, c in enumerate(chunk_cols):
                    h = col_to_human.get(c)
                    if h is not None and h in hp:
                        x = x_seq + ci * char_w
                        ax.add_patch(plt.Rectangle(
                            (x - char_w * 0.45, block_top),
                            char_w * 0.9, block_bot - block_top,
                            fc=col_bg, ec="none", alpha=0.35, zorder=-1))

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
                        (x - char_w * 0.45, y - box_half_h),
                        char_w * 0.9, box_h_aa,
                        fc=bg, ec="none", zorder=0))
                    ax.text(x, y, aa, va="center", ha="center",
                            color=AA_FG.get(aa, "#333"), **bold)
                elif sym == ":":
                    bg = AA_BG.get(aa, "#F5F5F5")
                    ax.add_patch(plt.Rectangle(
                        (x - char_w * 0.45, y - box_half_h),
                        char_w * 0.9, box_h_aa,
                        fc=bg, ec="none", alpha=0.4, zorder=0))
                    ax.text(x, y, aa, va="center", ha="center",
                            color=AA_FG.get(aa, "#555"), **font)
                else:
                    ax.text(x, y, aa, va="center", ha="center",
                            color="#666", **font)
            y += line_h

        sym_font_star = dict(
            family="monospace", fontsize=fs_cons_star, fontweight="bold")
        sym_font_colon = dict(family="monospace", fontsize=fs_cons_colon,
                              fontweight="bold")
        sym_font_dot = dict(family="monospace", fontsize=fs_cons_dot)

        for ci, c in enumerate(chunk_cols):
            sym = conservation_symbol(alignment, c)
            x = x_seq + ci * char_w
            if sym == "*":
                ax.text(x, y, sym, va="center", ha="center",
                        color="#000", **sym_font_star)
            elif sym == ":":
                ax.text(x, y, sym, va="center", ha="center",
                        color="#555", **sym_font_colon)
            elif sym == ".":
                ax.text(x, y, sym, va="center", ha="center",
                        color="#AAA", **sym_font_dot)
        y += y_cons_gap * line_h

        y += y_block_gap * line_h

    # Keep the alignment coordinate system top-to-bottom. Calling set_ylim with
    # ascending values after invert_yaxis() flips the rendered panel, moving the
    # title to the bottom and reversing species order.
    ax.set_ylim(y + 1.5 * line_h, -0.5)

    fig.tight_layout(pad=0.3)
    path = os.path.join(OUTPUT_DIR, filename)
    fig.savefig(path, dpi=dpi, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    suffix = " (presentation)" if for_presentation else ""
    print(f"  Saved {filename}{suffix}")
    return path


# ── Panel B: Domain architecture ─────────────────────────────────────────
def plot_domain_architecture(alignment, col_to_human, human_to_col,
                              human_idx, prefix):
    n_seqs = len(alignment)
    labels = [get_label(rec.id) for rec in alignment]

    domains = [
        ("UBR box (118-189)",     118,  189, "#7E57C2"),
        ("Region 400-600",        400,  600, "#1E88E5"),
        ("RING domain (1232-1376)", 1232, 1376, "#EF6C00"),
    ]

    species_lengths = []
    for seq_i in range(n_seqs):
        seq_len = sum(1 for c in range(alignment.get_alignment_length())
                      if alignment[seq_i, c] != "-")
        species_lengths.append(seq_len)
    max_len = max(species_lengths)

    fig_h = max(5, n_seqs * 0.55 + 2.5)
    fig_w_arch = max(10, (60 + max(len(l) for l in labels) + 6) * 6.2 / 72 + 0.5)
    fig, ax = plt.subplots(figsize=(fig_w_arch, fig_h))

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


def _slide_label_font():
    from PIL import ImageFont
    try:
        return ImageFont.truetype("arialbd.ttf", 34)
    except OSError:
        try:
            return ImageFont.truetype("C:/Windows/Fonts/arialbd.ttf", 34)
        except OSError:
            try:
                return ImageFont.truetype("C:/Windows/Fonts/calibrib.ttf", 34)
            except OSError:
                return ImageFont.load_default()


def _prepare_slide_montage_arrays(paths_a, path_b):
    """Resize panels to a capped width suitable for 16×9 slides."""
    from PIL import Image

    imgs_a = [Image.open(p) for p in paths_a]
    img_arch = Image.open(path_b)
    if len(imgs_a) < 3:
        return None

    def resize_to_width(im, tw):
        if im.width == tw:
            return im
        scale = tw / im.width
        return im.resize((tw, max(1, int(im.height * scale))), Image.LANCZOS)

    base_w = max(im.width for im in imgs_a)
    target_inner = min(4000, max(3000, int(1.15 * base_w)))
    margin = 28
    col_gap = 24
    row_gap = 22
    label_h = 44
    target_w = target_inner + 2 * margin
    half_w = (target_inner - col_gap) // 2

    img_tl = resize_to_width(imgs_a[0], half_w)
    img_tr = resize_to_width(imgs_a[1], half_w)
    row_ab_h = max(img_tl.height, img_tr.height)
    img_c = resize_to_width(imgs_a[2], target_inner)
    img_d = resize_to_width(img_arch, target_inner)

    return dict(
        img_tl=img_tl, img_tr=img_tr, img_c=img_c, img_d=img_d,
        row_ab_h=row_ab_h,
        target_inner=target_inner, half_w=half_w, margin=margin,
        col_gap=col_gap, row_gap=row_gap, label_h=label_h,
        target_w=target_w,
    )


def create_combined_figure_slide(paths_a, path_b, prefix):
    """Widescreen composite for PowerPoint: (a)|(b), then (c), then (d).

    One file for convenience; tall — prefer split slide PNGs for full-screen use."""
    from PIL import Image, ImageDraw

    pack = _prepare_slide_montage_arrays(paths_a, path_b)
    if pack is None:
        print("  Slide layout skipped (need three alignment PNGs)")
        return None

    margin = pack["margin"]
    col_gap = pack["col_gap"]
    row_gap = pack["row_gap"]
    label_h = pack["label_h"]
    target_w = pack["target_w"]
    half_w = pack["half_w"]
    img_tl = pack["img_tl"]
    img_tr = pack["img_tr"]
    img_c = pack["img_c"]
    img_d = pack["img_d"]
    row_ab_h = pack["row_ab_h"]

    total_h = (
        margin +
        row_gap + label_h + row_ab_h +
        row_gap + label_h + img_c.height +
        row_gap + label_h + img_d.height +
        margin
    )

    combined = Image.new("RGB", (target_w, total_h), "white")
    draw = ImageDraw.Draw(combined)
    lab_font = _slide_label_font()

    def paste_row_label(letter: str, x_off: int, y: int):
        draw.text((x_off, y), f"{letter}.", fill="black", font=lab_font)

    y = margin

    paste_row_label("a", margin, y)
    paste_row_label("b", margin + half_w + col_gap, y)
    y += label_h
    paste_y = y + (row_ab_h - img_tl.height) // 2
    combined.paste(img_tl, (margin, paste_y))
    paste_y_tr = y + (row_ab_h - img_tr.height) // 2
    combined.paste(img_tr, (margin + half_w + col_gap, paste_y_tr))
    y += row_ab_h + row_gap

    paste_row_label("c", margin, y)
    y += label_h
    combined.paste(img_c, (margin, y))
    y += img_c.height + row_gap

    paste_row_label("d", margin, y)
    y += label_h
    combined.paste(img_d, (margin, y))

    path = os.path.join(OUTPUT_DIR, f"{prefix}_combined_panels_slide.png")
    combined.save(path, dpi=(300, 300))
    print(f"  Saved {os.path.basename(path)}")
    print(f"  Slide layout size: {target_w} x {total_h} px (use for PPTX)")
    return path


def create_slide_figure_splits(paths_a, path_b, prefix):
    """Separate PNGs (a+b | c | d) so each fits ~one 16×9 slide."""
    from PIL import Image, ImageDraw

    pack = _prepare_slide_montage_arrays(paths_a, path_b)
    if pack is None:
        return []

    margin = pack["margin"]
    col_gap = pack["col_gap"]
    label_h = pack["label_h"]
    target_w = pack["target_w"]
    half_w = pack["half_w"]
    img_tl = pack["img_tl"]
    img_tr = pack["img_tr"]
    img_c = pack["img_c"]
    img_d = pack["img_d"]
    row_ab_h = pack["row_ab_h"]
    lab_font = _slide_label_font()

    out = []

    h_ab = margin + label_h + row_ab_h + margin
    cab = Image.new("RGB", (target_w, h_ab), "white")
    dab = ImageDraw.Draw(cab)
    dab.text((margin, margin), "a.", fill="black", font=lab_font)
    dab.text((margin + half_w + col_gap, margin), "b.", fill="black",
             font=lab_font)
    y0 = margin + label_h
    cab.paste(img_tl, (margin, y0 + (row_ab_h - img_tl.height) // 2))
    cab.paste(img_tr, (margin + half_w + col_gap,
                       y0 + (row_ab_h - img_tr.height) // 2))
    p_ab = os.path.join(OUTPUT_DIR, f"{prefix}_slide_panel_ab.png")
    cab.save(p_ab, dpi=(300, 300))
    print(f"  Saved {os.path.basename(p_ab)}")
    out.append(p_ab)

    h_c = margin + label_h + img_c.height + margin
    cc = Image.new("RGB", (target_w, h_c), "white")
    dc = ImageDraw.Draw(cc)
    dc.text((margin, margin), "c.", fill="black", font=lab_font)
    cc.paste(img_c, (margin, margin + label_h))
    p_c = os.path.join(OUTPUT_DIR, f"{prefix}_slide_panel_c.png")
    cc.save(p_c, dpi=(300, 300))
    print(f"  Saved {os.path.basename(p_c)}")
    out.append(p_c)

    h_d = margin + label_h + img_d.height + margin
    cd = Image.new("RGB", (target_w, h_d), "white")
    dd = ImageDraw.Draw(cd)
    dd.text((margin, margin), "d.", fill="black", font=lab_font)
    cd.paste(img_d, (margin, margin + label_h))
    p_d = os.path.join(OUTPUT_DIR, f"{prefix}_slide_panel_d.png")
    cd.save(p_d, dpi=(300, 300))
    print(f"  Saved {os.path.basename(p_d)}")
    out.append(p_d)

    return out


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
        ("RING Domain (aa 1232-1376)", 1232, 1376, "ring_domain"),
    ]

    panel_a_paths = []
    panel_a_paths_ppt = []
    for title, h_start, h_end, suffix in regions:
        print(f"\n--- {title} ---")
        col_start, col_end = get_alignment_slice(human_to_col, h_start, h_end)
        if col_start is None or col_end is None:
            print(f"  Could not find alignment columns for {h_start}-{h_end}")
            continue
        print(f"  Alignment columns: {col_start}-{col_end}")
        region_subs = {p: n for p, n in SUBSTRATE_RESIDUES.items()
                       if h_start <= p <= h_end}
        region_ring = {p: n for p, n in RING_RESIDUES.items()
                       if h_start <= p <= h_end}
        path = render_clustal_figure(
            alignment, col_start, col_end, col_to_human,
            title, f"{prefix}_{suffix}.png", block_width=60,
            highlight_positions=region_subs if region_subs else None,
            highlight_positions_2=region_ring if region_ring else None,
            for_presentation=True)
        panel_a_paths.append(path)
        panel_a_paths_ppt.append(path)

    print(f"\n--- Domain Architecture ---")
    path_b = plot_domain_architecture(
        alignment, col_to_human, human_to_col, human_idx, prefix)

    if panel_a_paths_ppt and path_b:
        print(f"\n--- Slide composites (large-font PNGs) ---")
        try:
            create_combined_figure_slide(panel_a_paths_ppt, path_b, prefix)
        except Exception as exc:
            print(f"  Slide-layout combined figure failed: {exc}")
        try:
            create_slide_figure_splits(panel_a_paths_ppt, path_b, prefix)
        except Exception as exc:
            print(f"  Split slide PNGs failed: {exc}")

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
