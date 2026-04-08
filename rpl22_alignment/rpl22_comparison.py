"""
RPL22 vs RPL22L1 paralog comparison:
  - Pairwise alignment of human RPL22 vs RPL22L1 (ClustalW-style figure)
  - Cross-species conservation overlay
  - Divergence highlight figure
"""
import os
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

from Bio import AlignIO, SeqIO
from Bio.Align import PairwiseAligner

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

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

STRONG_GROUPS = [
    set("STA"), set("NEQK"), set("NHQK"), set("NDEQ"), set("QHRK"),
    set("MILV"), set("MILF"), set("HY"), set("FYW"),
]
WEAK_GROUPS = [
    set("CSA"), set("ATV"), set("SAG"), set("STNK"), set("STPA"),
    set("SGND"), set("SNDEQK"), set("NDEQHK"), set("NEQHRK"),
    set("FVLIM"), set("HFY"),
]


def conservation_symbol_pair(aa1, aa2):
    if aa1 == "-" or aa2 == "-":
        return " "
    if aa1 == aa2:
        return "*"
    pair = {aa1, aa2}
    for g in STRONG_GROUPS:
        if pair.issubset(g):
            return ":"
    for g in WEAK_GROUPS:
        if pair.issubset(g):
            return "."
    return " "


def get_human_seq(fasta_path):
    alignment = AlignIO.read(fasta_path, "fasta")
    for rec in alignment:
        if "sapiens" in rec.id.lower() or "human" in rec.id.lower():
            return str(rec.seq).replace("-", "")
    return str(alignment[0].seq).replace("-", "")


def pairwise_align(seq1, seq2):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.substitution_matrix = None
    aligner.match_score = 2
    aligner.mismatch_score = -1
    alignments = aligner.align(seq1, seq2)
    best = alignments[0]
    aligned_seqs = best.format("fasta").split("\n")
    aln_seq1 = ""
    aln_seq2 = ""
    current = 1
    for line in aligned_seqs:
        if line.startswith(">"):
            current += 1
            continue
        if current == 2:
            aln_seq1 += line.strip()
        elif current == 3:
            aln_seq2 += line.strip()
    return aln_seq1, aln_seq2


def render_pairwise_figure(seq1, seq2, name1, name2, filename,
                            block_width=60):
    n_cols = len(seq1)
    n_blocks = (n_cols + block_width - 1) // block_width
    labels = [name1, name2]
    max_label = max(len(l) for l in labels)
    lines_per_block = 5
    total_lines = n_blocks * lines_per_block + 3

    char_w_pts = 6.2
    label_space = max_label * 0.52
    fig_w = max(10, (block_width + label_space + 2) * char_w_pts / 72 + 0.5)
    fig_h = max(3, total_lines * 0.18 + 1.0)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, fig_w)
    ax.set_ylim(0, total_lines + 2)
    ax.invert_yaxis()
    ax.axis("off")

    font = {"family": "monospace", "fontsize": 6}
    bold = {"family": "monospace", "fontsize": 6, "fontweight": "bold"}

    char_w = fig_w / (block_width + label_space + 2)
    x_seq = (label_space + 0.1) * char_w
    y = 0.0

    ax.text(fig_w / 2, y, f"Human {name1} vs {name2} Pairwise Alignment",
            ha="center", va="top", fontsize=10, fontweight="bold",
            fontstyle="italic", family="sans-serif")
    y += 1.8

    identity = sum(1 for a, b in zip(seq1, seq2) if a == b and a != "-")
    compared = sum(1 for a, b in zip(seq1, seq2) if a != "-" and b != "-")
    pct = identity / compared * 100 if compared > 0 else 0
    ax.text(fig_w / 2, y, f"Identity: {identity}/{compared} ({pct:.1f}%)",
            ha="center", va="top", fontsize=7, color="#666",
            family="sans-serif")
    y += 1.2

    for blk in range(n_blocks):
        start = blk * block_width
        end = min(start + block_width, n_cols)
        chunk = list(range(start, end))

        pos1 = sum(1 for c in seq1[:start] if c != "-")
        pos2 = sum(1 for c in seq2[:start] if c != "-")

        for seq_idx, (seq, label) in enumerate([(seq1, name1), (seq2, name2)]):
            color = "#B71C1C" if seq_idx == 0 else "#0D47A1"
            ax.text(0.1, y, label, va="center", color=color, **bold)

            for ci, c in enumerate(chunk):
                aa = seq[c]
                x = x_seq + ci * char_w
                if aa == "-":
                    ax.text(x, y, "-", va="center", ha="center",
                            color="#DDD", **font)
                    continue

                other_aa = seq2[c] if seq_idx == 0 else seq1[c]
                is_identical = aa == other_aa and other_aa != "-"

                if is_identical:
                    bg = AA_BG.get(aa, "#E0E0E0")
                    ax.add_patch(plt.Rectangle(
                        (x - char_w * 0.45, y - 0.42),
                        char_w * 0.9, 0.84,
                        fc=bg, ec="none", zorder=0))
                    ax.text(x, y, aa, va="center", ha="center",
                            color=AA_FG.get(aa, "#333"), **bold)
                else:
                    sym = conservation_symbol_pair(aa, other_aa)
                    if sym == ":":
                        bg = AA_BG.get(aa, "#F5F5F5")
                        ax.add_patch(plt.Rectangle(
                            (x - char_w * 0.45, y - 0.42),
                            char_w * 0.9, 0.84,
                            fc=bg, ec="none", alpha=0.4, zorder=0))
                    ax.text(x, y, aa, va="center", ha="center",
                            color=AA_FG.get(aa, "#666"), **font)
            y += 1

        for ci, c in enumerate(chunk):
            sym = conservation_symbol_pair(seq1[c], seq2[c])
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


def compute_per_species_conservation(fasta_path):
    alignment = AlignIO.read(fasta_path, "fasta")
    n_seqs = len(alignment)
    n_cols = alignment.get_alignment_length()

    human_idx = 0
    for i, rec in enumerate(alignment):
        if "sapiens" in rec.id.lower():
            human_idx = i
            break

    scores = np.zeros(n_cols)
    for col in range(n_cols):
        residues = [alignment[row, col] for row in range(n_seqs)]
        non_gap = [r for r in residues if r != "-"]
        if non_gap:
            counts = Counter(non_gap)
            scores[col] = counts.most_common(1)[0][1] / len(non_gap)

    human_pos = []
    res_num = 0
    for col in range(n_cols):
        if alignment[human_idx, col] != "-":
            res_num += 1
        human_pos.append(res_num)

    return np.array(human_pos, dtype=float), scores


def plot_conservation_comparison(rpl22_fasta, rpl22l1_fasta):
    print("\n--- Cross-species conservation comparison ---")

    hx1, sc1 = compute_per_species_conservation(rpl22_fasta)
    hx2, sc2 = compute_per_species_conservation(rpl22l1_fasta)

    window = 5
    kernel = np.ones(window) / window
    sm1 = np.convolve(sc1, kernel, mode="same")
    sm2 = np.convolve(sc2, kernel, mode="same")

    fig, ax = plt.subplots(figsize=(14, 6))
    ax.fill_between(hx1, sm1, alpha=0.2, color="#1565C0")
    ax.plot(hx1, sm1, linewidth=1.2, color="#1565C0", label="RPL22")
    ax.fill_between(hx2, sm2, alpha=0.2, color="#C62828")
    ax.plot(hx2, sm2, linewidth=1.2, color="#C62828", label="RPL22L1")

    ax.axhline(0.9, ls="--", color="gray", alpha=0.4, lw=0.8)
    ax.axhline(0.7, ls="--", color="gray", alpha=0.3, lw=0.8)

    ax.set_xlim(0, max(max(hx1), max(hx2)) + 1)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Human residue position", fontsize=12)
    ax.set_ylabel("Conservation score (across 13 species)", fontsize=12)
    ax.set_title("RPL22 vs RPL22L1: Cross-Species Conservation Comparison",
                 fontsize=14, fontweight="bold")
    ax.legend(fontsize=11, loc="lower right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, "rpl22_vs_rpl22l1_conservation.png")
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {os.path.basename(path)}")
    return path


def plot_divergence_map(seq1, seq2, name1, name2):
    """Highlight identical vs divergent positions between the two paralogs."""
    n = len(seq1)
    identity = []
    similar = []
    different = []
    pos1 = 0
    for i in range(n):
        a, b = seq1[i], seq2[i]
        if a == "-" and b == "-":
            continue
        pos1 += 1
        if a == b:
            identity.append(pos1)
        elif conservation_symbol_pair(a, b) in (":", "."):
            similar.append(pos1)
        else:
            different.append(pos1)

    fig, ax = plt.subplots(figsize=(14, 3))
    if identity:
        ax.bar(identity, [1]*len(identity), width=1, color="#1565C0",
               alpha=0.7, label="Identical")
    if similar:
        ax.bar(similar, [1]*len(similar), width=1, color="#FFB74D",
               alpha=0.7, label="Similar")
    if different:
        ax.bar(different, [1]*len(different), width=1, color="#E53935",
               alpha=0.7, label="Different")

    total_compared = len(identity) + len(similar) + len(different)
    pct_id = len(identity) / total_compared * 100 if total_compared else 0
    pct_sim = (len(identity) + len(similar)) / total_compared * 100 if total_compared else 0

    ax.set_xlim(0, max(pos1 + 1, 1))
    ax.set_ylim(0, 1.2)
    ax.set_yticks([])
    ax.set_xlabel("Alignment position", fontsize=11)
    ax.set_title(f"{name1} vs {name2}: Residue-by-Residue Comparison  "
                 f"({pct_id:.1f}% identical, {pct_sim:.1f}% similar+identical)",
                 fontsize=12, fontweight="bold")
    ax.legend(fontsize=9, loc="upper right", ncol=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, "rpl22_vs_rpl22l1_divergence_map.png")
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {os.path.basename(path)}")
    return path


def main():
    print("=" * 60)
    print("RPL22 vs RPL22L1 Paralog Comparison")
    print("=" * 60)

    rpl22_fasta = os.path.join(OUTPUT_DIR, "rpl22_aligned.fasta")
    rpl22l1_fasta = os.path.join(OUTPUT_DIR, "rpl22l1_aligned.fasta")

    if not os.path.exists(rpl22_fasta) or not os.path.exists(rpl22l1_fasta):
        print("Need both rpl22_aligned.fasta and rpl22l1_aligned.fasta")
        return

    print("\n--- Extracting human sequences ---")
    seq_rpl22 = get_human_seq(rpl22_fasta)
    seq_rpl22l1 = get_human_seq(rpl22l1_fasta)
    print(f"  RPL22:  {len(seq_rpl22)} aa")
    print(f"  RPL22L1: {len(seq_rpl22l1)} aa")

    print("\n--- Pairwise alignment ---")
    aln1, aln2 = pairwise_align(seq_rpl22, seq_rpl22l1)
    print(f"  Aligned length: {len(aln1)}")

    print("\n--- Pairwise alignment figure ---")
    render_pairwise_figure(aln1, aln2, "RPL22", "RPL22L1",
                            "rpl22_vs_rpl22l1_pairwise.png")

    print("\n--- Divergence map ---")
    plot_divergence_map(aln1, aln2, "RPL22", "RPL22L1")

    plot_conservation_comparison(rpl22_fasta, rpl22l1_fasta)

    print(f"\n{'='*60}")
    print("DONE")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
