"""
Deep analysis of UBR3 conserved regions:
  - Annotated alignment figures with actual AA letters for conserved blocks
  - Zoomed-in per-region alignment panels
  - Conserved-only heatmap (non-conserved columns stripped)
  - Annotated CSV with human UBR3 residue mapping
"""

import os
import csv
import textwrap
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

from Bio import AlignIO

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))
ALIGNED_FASTA = os.path.join(OUTPUT_DIR, "ubr3_aligned.fasta")

AA_COLORS = {
    "D": "#E53935", "E": "#E53935",
    "R": "#1E88E5", "K": "#1E88E5", "H": "#42A5F5",
    "S": "#43A047", "T": "#43A047", "N": "#66BB6A", "Q": "#66BB6A",
    "A": "#FFA726", "V": "#FFA726", "I": "#FFA726", "L": "#FFA726",
    "M": "#FFA726", "F": "#FF7043", "W": "#FF7043", "Y": "#FF7043",
    "P": "#AB47BC", "G": "#78909C",
    "C": "#FFEE58",
    "-": "#FFFFFF",
}
DEFAULT_AA_COLOR = "#BDBDBD"


def load_alignment():
    alignment = AlignIO.read(ALIGNED_FASTA, "fasta")
    return alignment


def compute_conservation(alignment) -> np.ndarray:
    n_seqs = len(alignment)
    n_cols = alignment.get_alignment_length()
    scores = np.zeros(n_cols)
    for col in range(n_cols):
        residues = [alignment[row, col] for row in range(n_seqs)]
        non_gap = [r for r in residues if r != "-"]
        if not non_gap:
            continue
        counts = Counter(non_gap)
        scores[col] = counts.most_common(1)[0][1] / len(non_gap)
    return scores


def get_consensus(alignment, col: int) -> str:
    residues = [alignment[row, col] for row in range(len(alignment))]
    non_gap = [r for r in residues if r != "-"]
    if not non_gap:
        return "-"
    counts = Counter(non_gap)
    return counts.most_common(1)[0][0]


def find_conserved_regions(scores, threshold=0.9, min_len=5):
    regions = []
    in_region = False
    start = 0
    for i, s in enumerate(scores):
        if s >= threshold and not in_region:
            start = i
            in_region = True
        elif s < threshold and in_region:
            if i - start >= min_len:
                regions.append((start, i - 1, float(np.mean(scores[start:i]))))
            in_region = False
    if in_region and len(scores) - start >= min_len:
        regions.append((start, len(scores) - 1,
                        float(np.mean(scores[start:]))))
    return regions


def alignment_col_to_human_pos(alignment):
    """Map alignment column index -> human residue number (1-based), or None for gaps."""
    human_idx = None
    for i, rec in enumerate(alignment):
        if rec.id.lower().startswith("human"):
            human_idx = i
            break
    if human_idx is None:
        human_idx = 0

    col_to_human = {}
    res_num = 0
    for col in range(alignment.get_alignment_length()):
        aa = alignment[human_idx, col]
        if aa != "-":
            res_num += 1
            col_to_human[col] = res_num
        else:
            col_to_human[col] = None
    return col_to_human, human_idx


def nearest_human_pos(col_to_human, col, direction="both"):
    """Find nearest non-None human position for a gap column."""
    pos = col_to_human.get(col)
    if pos is not None:
        return pos
    max_col = max(col_to_human.keys())
    for offset in range(1, max_col + 1):
        if direction in ("both", "left") and col - offset >= 0:
            p = col_to_human.get(col - offset)
            if p is not None:
                return p
        if direction in ("both", "right") and col + offset <= max_col:
            p = col_to_human.get(col + offset)
            if p is not None:
                return p
    return None


def human_range_label(col_to_human, rstart, rend):
    """Get a human-position label for an alignment region, handling gaps."""
    h_start = col_to_human.get(rstart) or nearest_human_pos(col_to_human, rstart, "right")
    h_end = col_to_human.get(rend) or nearest_human_pos(col_to_human, rend, "left")
    if h_start and h_end:
        return h_start, h_end, f"Human {h_start}-{h_end}"
    elif h_start:
        return h_start, h_start, f"Human ~{h_start}"
    elif h_end:
        return h_end, h_end, f"Human ~{h_end}"
    return None, None, f"Aln {rstart}-{rend}"


# ── Figure 1: Zoomed alignment blocks for each conserved region ───────────
def plot_region_alignment_blocks(alignment, scores, regions, col_to_human,
                                  human_idx, prefix):
    """For each conserved region, draw a colored alignment block with AA letters."""
    n_seqs = len(alignment)
    labels = [rec.id for rec in alignment]

    COLS_PER_PAGE = 80
    region_groups = group_regions_for_pages(regions, COLS_PER_PAGE, gap=5)

    paths = []
    for page_idx, group in enumerate(region_groups):
        all_cols = []
        region_boundaries = []
        for rstart, rend, rmean in group:
            boundary_start = len(all_cols)
            for c in range(rstart, rend + 1):
                all_cols.append(c)
            region_boundaries.append((boundary_start, len(all_cols) - 1,
                                      rstart, rend, rmean))
            all_cols.append(None)
        if all_cols and all_cols[-1] is None:
            all_cols.pop()

        n_display = len(all_cols)
        cell_w, cell_h = 0.28, 0.38
        fig_w = max(8, n_display * cell_w + 3.5)
        fig_h = max(4, (n_seqs + 3) * cell_h + 1.5)

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        ax.set_xlim(-0.5, n_display - 0.5)
        ax.set_ylim(-0.5, n_seqs + 2.5)
        ax.invert_yaxis()
        ax.set_aspect("equal")
        ax.axis("off")

        for rb_start, rb_end, rstart, rend, rmean in region_boundaries:
            _, _, h_label = human_range_label(col_to_human, rstart, rend)
            ax.text((rb_start + rb_end) / 2, -0.3,
                    f"pos {rstart}-{rend}  ({h_label})\ncons={rmean:.2f}",
                    ha="center", va="bottom", fontsize=6, fontweight="bold",
                    color="#333333")

        for x_idx, col in enumerate(all_cols):
            if col is None:
                for y in range(n_seqs + 2):
                    ax.add_patch(plt.Rectangle((x_idx - 0.5, y - 0.5), 1, 1,
                                                fc="#F5F5F5", ec="none"))
                continue

            cons = scores[col]
            consensus_aa = get_consensus(alignment, col)

            # Row 0: consensus
            bg = AA_COLORS.get(consensus_aa, DEFAULT_AA_COLOR)
            ax.add_patch(plt.Rectangle((x_idx - 0.5, -0.5), 1, 1,
                                        fc=bg, ec="white", lw=0.3, alpha=0.9))
            ax.text(x_idx, 0, consensus_aa, ha="center", va="center",
                    fontsize=5.5, fontweight="bold", color="black",
                    family="monospace")

            # Row 1: conservation bar
            bar_color = "#B71C1C" if cons >= 1.0 else (
                "#E53935" if cons >= 0.9 else "#FF9800")
            ax.add_patch(plt.Rectangle((x_idx - 0.5, 0.5), 1, 1,
                                        fc=bar_color, ec="white", lw=0.3,
                                        alpha=cons))
            ax.text(x_idx, 1, f"{cons:.0%}" if cons < 1 else "*",
                    ha="center", va="center", fontsize=4, color="white",
                    family="monospace")

            # Rows 2..n_seqs+1: per-species
            for seq_i in range(n_seqs):
                y = seq_i + 2
                aa = alignment[seq_i, col]
                if aa == "-":
                    fc = "#F5F5F5"
                    txt_color = "#BDBDBD"
                elif aa == consensus_aa:
                    fc = AA_COLORS.get(aa, DEFAULT_AA_COLOR)
                    txt_color = "black"
                else:
                    fc = "#FFFFFF"
                    txt_color = "#E53935"
                ax.add_patch(plt.Rectangle((x_idx - 0.5, y - 0.5), 1, 1,
                                            fc=fc, ec="#E0E0E0", lw=0.2))
                ax.text(x_idx, y, aa, ha="center", va="center",
                        fontsize=5, color=txt_color, family="monospace",
                        fontweight="bold" if aa == consensus_aa else "normal")

        # Y-axis labels
        ax.text(-1, 0, "Consensus", ha="right", va="center", fontsize=7,
                fontweight="bold")
        ax.text(-1, 1, "Conserv.", ha="right", va="center", fontsize=6,
                fontstyle="italic", color="#666")
        for seq_i in range(n_seqs):
            weight = "bold" if seq_i == human_idx else "normal"
            ax.text(-1, seq_i + 2, labels[seq_i], ha="right", va="center",
                    fontsize=7, fontweight=weight)

        title_regions = ", ".join(
            f"{rs}-{re}" for _, _, rs, re, _ in region_boundaries)
        ax.set_title(f"UBR3 Conserved Regions (>=90%) - Alignment Block {page_idx+1}\n"
                     f"Positions: {title_regions}",
                     fontsize=10, fontweight="bold", pad=20)

        fig.tight_layout()
        path = os.path.join(OUTPUT_DIR,
                            f"{prefix}_annotated_block_{page_idx+1:02d}.png")
        fig.savefig(path, dpi=250, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)
        print(f"  Saved {os.path.basename(path)}")

    return paths


def group_regions_for_pages(regions, max_cols=80, gap=1):
    """Group conserved regions into pages that fit within max_cols."""
    pages = []
    current_page = []
    current_cols = 0
    for rstart, rend, rmean in regions:
        region_len = rend - rstart + 1
        needed = region_len + (gap if current_page else 0)
        if current_cols + needed > max_cols and current_page:
            pages.append(current_page)
            current_page = []
            current_cols = 0
        current_page.append((rstart, rend, rmean))
        current_cols += region_len + (gap if len(current_page) > 1 else 0)
    if current_page:
        pages.append(current_page)
    return pages


# ── Figure 2: Conserved-only heatmap ─────────────────────────────────────
def plot_conserved_only_heatmap(alignment, scores, regions, col_to_human,
                                 prefix):
    """Heatmap showing ONLY conserved columns, with consensus AA annotated."""
    n_seqs = len(alignment)
    labels = [rec.id for rec in alignment]

    conserved_cols = []
    region_labels = []
    for rstart, rend, rmean in regions:
        for c in range(rstart, rend + 1):
            conserved_cols.append(c)
        region_labels.append((len(conserved_cols) - (rend - rstart + 1),
                              len(conserved_cols) - 1, rstart, rend))

    n_cons = len(conserved_cols)
    if n_cons == 0:
        print("  No conserved columns to plot.")
        return None

    aa_to_num = {
        "A": 1, "R": 2, "N": 3, "D": 4, "C": 5, "E": 6, "Q": 7,
        "G": 8, "H": 9, "I": 10, "L": 11, "K": 12, "M": 13, "F": 14,
        "P": 15, "S": 16, "T": 17, "W": 18, "Y": 19, "V": 20, "-": 0,
    }

    matrix = np.zeros((n_seqs, n_cons))
    for x, col in enumerate(conserved_cols):
        consensus = get_consensus(alignment, col)
        for seq_i in range(n_seqs):
            aa = alignment[seq_i, col]
            if aa == "-":
                matrix[seq_i, x] = -1
            elif aa == consensus:
                matrix[seq_i, x] = 1.0
            else:
                matrix[seq_i, x] = 0.3

    cmap = LinearSegmentedColormap.from_list(
        "cons_match", ["#EEEEEE", "#FFCDD2", "#1565C0"])
    cmap.set_under("#F5F5F5")

    fig_w = max(14, n_cons * 0.06 + 3)
    fig, ax = plt.subplots(figsize=(fig_w, max(4, n_seqs * 0.5 + 1)))
    im = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=-0.5, vmax=1.0,
                   interpolation="none")
    ax.set_yticks(range(n_seqs))
    ax.set_yticklabels(labels, fontsize=9)

    for rl_start, rl_end, rstart, rend in region_labels:
        ax.axvline(rl_start - 0.5, color="#E0E0E0", lw=0.5, ls="--")

    tick_positions = []
    tick_labels_text = []
    for rl_start, rl_end, rstart, rend in region_labels:
        mid = (rl_start + rl_end) / 2
        h_s, h_e, _ = human_range_label(col_to_human, rstart, rend)
        if h_s and h_e:
            tick_labels_text.append(f"{h_s}-{h_e}")
        else:
            tick_labels_text.append(f"aln {rstart}-{rend}")
        tick_positions.append(mid)

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels_text, fontsize=5, rotation=90)
    ax.set_xlabel("Human UBR3 residue ranges (conserved regions only)", fontsize=10)
    ax.set_title("UBR3 Conserved Regions Only - Sequence Match Heatmap",
                 fontsize=13, fontweight="bold")

    legend_items = [
        mpatches.Patch(color="#1565C0", label="Matches consensus"),
        mpatches.Patch(color="#FFCDD2", label="Differs from consensus"),
        mpatches.Patch(color="#F5F5F5", label="Gap"),
    ]
    ax.legend(handles=legend_items, fontsize=8, loc="lower right",
              bbox_to_anchor=(1, -0.18), ncol=3)

    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, f"{prefix}_conserved_only_heatmap.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {os.path.basename(path)}")
    return path


# ── Figure 3: Top conserved regions bar chart with AA sequences ──────────
def plot_top_regions_annotated(alignment, scores, regions, col_to_human,
                                human_idx, prefix, top_n=20):
    """Bar chart of longest/most conserved regions with human AA sequence written out."""
    sorted_regions = sorted(regions, key=lambda r: r[1] - r[0], reverse=True)
    top = sorted_regions[:top_n]
    top = sorted(top, key=lambda r: r[0])

    fig, ax = plt.subplots(figsize=(18, max(6, len(top) * 0.7 + 2)))

    y_positions = []
    for i, (rstart, rend, rmean) in enumerate(top):
        length = rend - rstart + 1
        y = i

        human_seq = ""
        for c in range(rstart, rend + 1):
            human_seq += alignment[human_idx, c]

        _, _, pos_label = human_range_label(col_to_human, rstart, rend)

        color = "#B71C1C" if rmean >= 1.0 else (
            "#E53935" if rmean >= 0.95 else "#FF5722")
        ax.barh(y, length, color=color, alpha=0.85, edgecolor="white", height=0.6)

        ax.text(-1, y, pos_label, ha="right", va="center", fontsize=7,
                fontweight="bold")

        display_seq = human_seq.replace("-", "")
        if len(display_seq) > 60:
            display_seq = display_seq[:28] + "..." + display_seq[-28:]
        ax.text(length + 0.5, y, display_seq, ha="left", va="center",
                fontsize=5.5, family="monospace", color="#333")

        cons_label = f"{rmean:.0%}" if rmean < 1.0 else "100%"
        ax.text(length / 2, y, f"{length}aa ({cons_label})",
                ha="center", va="center", fontsize=6, color="white",
                fontweight="bold")

        y_positions.append(y)

    ax.set_yticks([])
    ax.set_xlabel("Region length (amino acids)", fontsize=11)
    ax.set_title(f"Top {len(top)} Conserved Regions in UBR3 (>=90%)\n"
                 f"with Human AA Sequence",
                 fontsize=13, fontweight="bold")
    ax.invert_yaxis()

    legend_items = [
        mpatches.Patch(color="#B71C1C", label="100% identical"),
        mpatches.Patch(color="#E53935", label="95-99% conserved"),
        mpatches.Patch(color="#FF5722", label="90-95% conserved"),
    ]
    ax.legend(handles=legend_items, fontsize=9, loc="lower right")

    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, f"{prefix}_top_regions_annotated.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {os.path.basename(path)}")
    return path


# ── Figure 4: Conservation map with consensus AA on top ──────────────────
def plot_conservation_with_consensus(alignment, scores, col_to_human,
                                      human_idx, prefix):
    """Full-length conservation plot with consensus AA letters printed above peaks."""
    n_cols = len(scores)

    window = 20
    kernel = np.ones(window) / window
    smoothed = np.convolve(scores, kernel, mode="same")

    fig, axes = plt.subplots(2, 1, figsize=(24, 8),
                              gridspec_kw={"height_ratios": [1, 4]})
    ax_aa, ax_cons = axes

    # Top panel: consensus AA at highly conserved positions
    ax_aa.set_xlim(0, n_cols)
    ax_aa.set_ylim(0, 1)
    ax_aa.axis("off")

    for col in range(n_cols):
        if scores[col] >= 1.0:
            aa = get_consensus(alignment, col)
            color = AA_COLORS.get(aa, DEFAULT_AA_COLOR)
            ax_aa.axvline(col, color=color, alpha=0.6, lw=0.4)

    step = max(1, n_cols // 200)
    for col in range(0, n_cols, step):
        if scores[col] >= 1.0:
            aa = get_consensus(alignment, col)
            ax_aa.text(col, 0.5, aa, ha="center", va="center",
                       fontsize=3, family="monospace", fontweight="bold",
                       color=AA_COLORS.get(aa, "#333"))

    ax_aa.set_title("UBR3 Conservation with Consensus Amino Acids (100% identical shown above)",
                    fontsize=13, fontweight="bold")

    # Bottom panel: conservation curve
    ax_cons.fill_between(range(n_cols), smoothed, alpha=0.3, color="#2196F3")
    ax_cons.plot(range(n_cols), smoothed, lw=0.8, color="#0D47A1")

    for col in range(n_cols):
        if scores[col] >= 1.0:
            ax_cons.axvline(col, color="#B71C1C", alpha=0.08, lw=0.5)

    ax_cons.axhline(0.9, ls="--", color="red", alpha=0.5, lw=1, label="90%")
    ax_cons.axhline(0.7, ls="--", color="orange", alpha=0.5, lw=1, label="70%")
    ax_cons.set_xlim(0, n_cols)
    ax_cons.set_ylim(0, 1.05)
    ax_cons.set_xlabel("Alignment position", fontsize=12)
    ax_cons.set_ylabel("Conservation", fontsize=12)
    ax_cons.legend(fontsize=9)

    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, f"{prefix}_conservation_with_consensus.png")
    fig.savefig(path, dpi=250, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {os.path.basename(path)}")
    return path


# ── Figure 5: Per-region detailed alignment panels ───────────────────────
def plot_individual_region_panels(alignment, scores, regions, col_to_human,
                                   human_idx, prefix, min_region_len=10):
    """One detailed figure per large conserved region."""
    labels = [rec.id for rec in alignment]
    n_seqs = len(alignment)
    large_regions = [(s, e, m) for s, e, m in regions if e - s + 1 >= min_region_len]

    paths = []
    for reg_i, (rstart, rend, rmean) in enumerate(large_regions):
        length = rend - rstart + 1
        cell_w = 0.32
        cell_h = 0.36
        fig_w = max(6, length * cell_w + 3)
        fig_h = max(3, (n_seqs + 2) * cell_h + 1.5)

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        ax.set_xlim(-0.5, length - 0.5)
        ax.set_ylim(-0.5, n_seqs + 1.5)
        ax.invert_yaxis()
        ax.set_aspect("equal")
        ax.axis("off")

        h_start, h_end, h_label = human_range_label(col_to_human, rstart, rend)

        # Column headers: human residue numbers
        for x, col in enumerate(range(rstart, rend + 1)):
            h_pos = col_to_human.get(col)
            if h_pos and (x == 0 or x == length - 1 or x % 10 == 0):
                ax.text(x, -0.3, str(h_pos), ha="center", va="bottom",
                        fontsize=4.5, color="#666", rotation=90)

        for x, col in enumerate(range(rstart, rend + 1)):
            cons = scores[col]
            consensus_aa = get_consensus(alignment, col)

            # Row 0: consensus
            bg = AA_COLORS.get(consensus_aa, DEFAULT_AA_COLOR)
            ax.add_patch(plt.Rectangle((x - 0.5, -0.5), 1, 1,
                                        fc=bg, ec="white", lw=0.3, alpha=0.85))
            ax.text(x, 0, consensus_aa, ha="center", va="center",
                    fontsize=6, fontweight="bold", family="monospace")

            # Rows 1..n_seqs: species
            for seq_i in range(n_seqs):
                y = seq_i + 1
                aa = alignment[seq_i, col]
                if aa == "-":
                    fc = "#F5F5F5"
                    tc = "#CCC"
                elif aa == consensus_aa:
                    fc = AA_COLORS.get(aa, DEFAULT_AA_COLOR)
                    tc = "black"
                else:
                    fc = "#FFF"
                    tc = "#D32F2F"
                ax.add_patch(plt.Rectangle((x - 0.5, y - 0.5), 1, 1,
                                            fc=fc, ec="#E0E0E0", lw=0.15))
                ax.text(x, y, aa, ha="center", va="center",
                        fontsize=5.5, color=tc, family="monospace",
                        fontweight="bold" if aa == consensus_aa else "normal")

        ax.text(-1, 0, "Consensus", ha="right", va="center", fontsize=7,
                fontweight="bold")
        for seq_i in range(n_seqs):
            w = "bold" if seq_i == human_idx else "normal"
            ax.text(-1, seq_i + 1, labels[seq_i], ha="right", va="center",
                    fontsize=7, fontweight=w)

        ax.set_title(
            f"Conserved Region: {h_label}  "
            f"({length} aa, {rmean:.0%} conserved)",
            fontsize=10, fontweight="bold", pad=15)

        fig.tight_layout()
        fn_start = h_start or rstart
        fn_end = h_end or rend
        path = os.path.join(OUTPUT_DIR,
                            f"{prefix}_region_{reg_i+1:02d}_h{fn_start}-{fn_end}.png")
        fig.savefig(path, dpi=250, bbox_inches="tight")
        plt.close(fig)
        paths.append(path)
        print(f"  Saved {os.path.basename(path)}")

    return paths


# ── Annotated CSV ─────────────────────────────────────────────────────────
def save_annotated_csv(alignment, scores, regions, col_to_human, human_idx,
                       prefix):
    """CSV with region details, human AA sequence, and consensus."""
    path = os.path.join(OUTPUT_DIR, f"{prefix}_conserved_regions_annotated.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "region_num", "aln_start", "aln_end", "length",
            "human_start", "human_end", "mean_conservation",
            "pct_100_identical", "human_AA_sequence", "consensus_sequence",
        ])
        for i, (rstart, rend, rmean) in enumerate(regions, 1):
            length = rend - rstart + 1
            h_start, h_end, _ = human_range_label(col_to_human, rstart, rend)

            human_seq = ""
            consensus_seq = ""
            n_identical = 0
            for c in range(rstart, rend + 1):
                human_seq += alignment[human_idx, c]
                consensus_seq += get_consensus(alignment, c)
                if scores[c] >= 1.0:
                    n_identical += 1

            pct_identical = n_identical / length * 100

            writer.writerow([
                i, rstart, rend, length,
                h_start or "", h_end or "", f"{rmean:.4f}",
                f"{pct_identical:.1f}",
                human_seq.replace("-", ""),
                consensus_seq,
            ])

    print(f"  Saved {os.path.basename(path)}")
    return path


# ── Figure 6: Summary overview - only conserved regions on the protein ───
def plot_conserved_regions_map(scores, regions, col_to_human, prefix):
    """Linear protein map showing only where conserved regions fall."""
    n_cols = len(scores)

    max_human = max((v for v in col_to_human.values() if v), default=n_cols)

    fig, ax = plt.subplots(figsize=(20, 3))

    ax.barh(0, max_human, height=0.3, color="#E0E0E0", edgecolor="none")

    for rstart, rend, rmean in regions:
        h_start, h_end, _ = human_range_label(col_to_human, rstart, rend)
        if h_start and h_end:
            width = h_end - h_start + 1
            color = "#B71C1C" if rmean >= 1.0 else (
                "#E53935" if rmean >= 0.95 else "#FF5722")
            ax.barh(0, width, left=h_start, height=0.3,
                    color=color, edgecolor="none", alpha=0.9)

    ax.set_xlim(0, max_human + 10)
    ax.set_ylim(-0.5, 0.5)
    ax.set_yticks([])
    ax.set_xlabel("Human UBR3 residue position", fontsize=12)
    ax.set_title("Conserved Regions Mapped onto Human UBR3 Protein",
                 fontsize=13, fontweight="bold")

    legend_items = [
        mpatches.Patch(color="#B71C1C", label="100% identical"),
        mpatches.Patch(color="#E53935", label="95-99%"),
        mpatches.Patch(color="#FF5722", label="90-95%"),
        mpatches.Patch(color="#E0E0E0", label="Not conserved (>=90%)"),
    ]
    ax.legend(handles=legend_items, fontsize=8, loc="upper right", ncol=4)

    for pos in [1, 500, 1000, 1500, max_human]:
        ax.axvline(pos, color="#999", lw=0.3, ls=":")
        ax.text(pos, -0.25, str(pos), ha="center", fontsize=7, color="#666")

    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, f"{prefix}_conserved_regions_map.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {os.path.basename(path)}")
    return path


# ── Main ──────────────────────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("UBR3 Conserved Regions - Deep Analysis")
    print("=" * 60)

    print("\nLoading alignment...")
    alignment = load_alignment()
    n_seqs = len(alignment)
    n_cols = alignment.get_alignment_length()
    print(f"  {n_seqs} sequences x {n_cols} columns")

    scores = compute_conservation(alignment)
    col_to_human, human_idx = alignment_col_to_human_pos(alignment)
    print(f"  Human sequence index: {human_idx} ({alignment[human_idx].id})")

    regions_90 = find_conserved_regions(scores, threshold=0.9, min_len=5)
    print(f"  Found {len(regions_90)} conserved regions (>=90%, >=5 aa)")

    prefix = "ubr3"

    print("\n--- Generating annotated alignment blocks ---")
    plot_region_alignment_blocks(alignment, scores, regions_90, col_to_human,
                                  human_idx, prefix)

    print("\n--- Generating conserved-only heatmap ---")
    plot_conserved_only_heatmap(alignment, scores, regions_90, col_to_human,
                                 prefix)

    print("\n--- Generating top regions bar chart with AA sequences ---")
    plot_top_regions_annotated(alignment, scores, regions_90, col_to_human,
                                human_idx, prefix)

    print("\n--- Generating conservation plot with consensus AA ---")
    plot_conservation_with_consensus(alignment, scores, col_to_human,
                                      human_idx, prefix)

    print("\n--- Generating individual region panels ---")
    plot_individual_region_panels(alignment, scores, regions_90, col_to_human,
                                   human_idx, prefix, min_region_len=10)

    print("\n--- Generating conserved regions protein map ---")
    plot_conserved_regions_map(scores, regions_90, col_to_human, prefix)

    print("\n--- Saving annotated CSV ---")
    save_annotated_csv(alignment, scores, regions_90, col_to_human,
                       human_idx, prefix)

    print(f"\n{'='*60}")
    print(f"DONE - all outputs in: {OUTPUT_DIR}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
