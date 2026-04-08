"""
Conservation analysis for RPL22 and RPL22L1 alignments:
  - Per-position conservation scoring
  - Conservation line plots, heatmaps, identity matrices
  - Annotated alignment blocks, zoomed region panels
  - Conserved-only heatmap, conserved regions map
  - CSV exports
"""
import os, csv
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap

from Bio import AlignIO

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

DISPLAY_NAMES = {
    "H_sapiens": "H. sapiens", "M_musculus": "M. musculus",
    "X_laevis": "X. laevis", "P_marinus": "P. marinus",
    "C_intestinalis": "C. intestinalis", "B_floridae": "B. floridae",
    "S_purpuratus": "S. purpuratus", "O_sinensis": "O. sinensis",
    "D_melanogaster": "D. melanogaster", "A_mellifera": "A. mellifera",
    "P_tepidariorum": "P. tepidariorum", "H_vulgaris": "H. vulgaris",
    "N_vectensis": "N. vectensis",
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

AA_COLORS = {
    "D": "#E53935", "E": "#E53935",
    "R": "#1E88E5", "K": "#1E88E5", "H": "#42A5F5",
    "S": "#43A047", "T": "#43A047", "N": "#66BB6A", "Q": "#66BB6A",
    "A": "#FFA726", "V": "#FFA726", "I": "#FFA726", "L": "#FFA726",
    "M": "#FFA726", "F": "#FF7043", "W": "#FF7043", "Y": "#FF7043",
    "P": "#AB47BC", "G": "#78909C", "C": "#FFEE58", "-": "#FFFFFF",
}
DEFAULT_AA_COLOR = "#BDBDBD"


def compute_conservation(alignment):
    n_seqs = len(alignment)
    n_cols = alignment.get_alignment_length()
    scores = np.zeros(n_cols)
    for col in range(n_cols):
        residues = [alignment[row, col] for row in range(n_seqs)]
        non_gap = [r for r in residues if r != "-"]
        if non_gap:
            counts = Counter(non_gap)
            scores[col] = counts.most_common(1)[0][1] / len(non_gap)
    return scores


def get_consensus(alignment, col):
    residues = [alignment[row, col] for row in range(len(alignment))]
    non_gap = [r for r in residues if r != "-"]
    if not non_gap:
        return "-"
    return Counter(non_gap).most_common(1)[0][0]


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


def find_conserved_regions(scores, threshold=0.7, min_len=3):
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
        regions.append((start, len(scores) - 1, float(np.mean(scores[start:]))))
    return regions


def nearest_human_pos(col_to_human, col):
    pos = col_to_human.get(col)
    if pos is not None:
        return pos
    max_col = max(col_to_human.keys())
    for offset in range(1, max_col + 1):
        if col - offset >= 0:
            p = col_to_human.get(col - offset)
            if p is not None:
                return p
        if col + offset <= max_col:
            p = col_to_human.get(col + offset)
            if p is not None:
                return p
    return None


def human_range_label(col_to_human, rstart, rend):
    h_start = col_to_human.get(rstart) or nearest_human_pos(col_to_human, rstart)
    h_end = col_to_human.get(rend) or nearest_human_pos(col_to_human, rend)
    if h_start and h_end:
        return h_start, h_end, f"Human {h_start}-{h_end}"
    return None, None, f"Aln {rstart}-{rend}"


# ── Conservation line plot ────────────────────────────────────────────────
def plot_conservation_line(alignment, scores, col_to_human, protein_name, prefix):
    n_cols = len(scores)
    human_x = []
    last_h = 0
    for col in range(n_cols):
        h = col_to_human.get(col)
        if h is not None:
            last_h = h
        human_x.append(last_h)
    human_x = np.array(human_x, dtype=float)

    window = min(10, n_cols // 5) if n_cols > 20 else 3
    kernel = np.ones(window) / window
    smoothed = np.convolve(scores, kernel, mode="same")

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.fill_between(human_x, smoothed, alpha=0.3, color="#2196F3")
    ax.plot(human_x, smoothed, linewidth=1.0, color="#0D47A1")
    ax.axhline(0.9, ls="--", color="red", alpha=0.5, lw=1, label="90%")
    ax.axhline(0.7, ls="--", color="orange", alpha=0.5, lw=1, label="70%")
    ax.set_xlim(0, max(human_x) + 1)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Human residue position", fontsize=12)
    ax.set_ylabel("Conservation score", fontsize=12)
    ax.set_title(f"{protein_name} Conservation Score Along Alignment",
                 fontsize=14, fontweight="bold")
    ax.legend(fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, f"{prefix}_conservation_line.png")
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {os.path.basename(path)}")
    return path


# ── Conservation heatmap ──────────────────────────────────────────────────
def plot_conservation_heatmap(alignment, scores, protein_name, prefix):
    n_seqs = len(alignment)
    n_cols = alignment.get_alignment_length()
    labels = [DISPLAY_NAMES.get(rec.id, rec.id) for rec in alignment]

    matrix = np.zeros((n_seqs, n_cols))
    for col in range(n_cols):
        consensus = get_consensus(alignment, col)
        for seq_i in range(n_seqs):
            aa = alignment[seq_i, col]
            if aa == "-":
                matrix[seq_i, col] = -1
            elif aa == consensus:
                matrix[seq_i, col] = 1.0
            else:
                matrix[seq_i, col] = 0.3

    cmap = LinearSegmentedColormap.from_list(
        "cons", ["#EEEEEE", "#FFCDD2", "#1565C0"])
    cmap.set_under("#F5F5F5")

    fig, ax = plt.subplots(figsize=(max(10, n_cols * 0.08 + 3),
                                     max(4, n_seqs * 0.5 + 1)))
    ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=-0.5, vmax=1.0,
              interpolation="none")
    ax.set_yticks(range(n_seqs))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Alignment position", fontsize=11)
    ax.set_title(f"{protein_name} Conservation Heatmap", fontsize=13, fontweight="bold")

    legend_items = [
        mpatches.Patch(color="#1565C0", label="Matches consensus"),
        mpatches.Patch(color="#FFCDD2", label="Differs"),
        mpatches.Patch(color="#F5F5F5", label="Gap"),
    ]
    ax.legend(handles=legend_items, fontsize=8, loc="lower right",
              bbox_to_anchor=(1, -0.15), ncol=3)

    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, f"{prefix}_conservation_heatmap.png")
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {os.path.basename(path)}")
    return path


# ── Identity matrix ───────────────────────────────────────────────────────
def plot_identity_matrix(alignment, protein_name, prefix):
    n_seqs = len(alignment)
    n_cols = alignment.get_alignment_length()
    labels = [DISPLAY_NAMES.get(rec.id, rec.id) for rec in alignment]

    matrix = np.zeros((n_seqs, n_seqs))
    for i in range(n_seqs):
        for j in range(n_seqs):
            matches = 0
            compared = 0
            for col in range(n_cols):
                a = alignment[i, col]
                b = alignment[j, col]
                if a != "-" and b != "-":
                    compared += 1
                    if a == b:
                        matches += 1
            matrix[i, j] = matches / compared * 100 if compared > 0 else 0

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(matrix, cmap="YlOrRd", vmin=20, vmax=100)
    ax.set_xticks(range(n_seqs))
    ax.set_xticklabels(labels, fontsize=8, rotation=45, ha="right")
    ax.set_yticks(range(n_seqs))
    ax.set_yticklabels(labels, fontsize=8)

    for i in range(n_seqs):
        for j in range(n_seqs):
            ax.text(j, i, f"{matrix[i,j]:.0f}%", ha="center", va="center",
                    fontsize=6, color="white" if matrix[i,j] > 60 else "black")

    ax.set_title(f"{protein_name} Pairwise Sequence Identity (%)",
                 fontsize=13, fontweight="bold")
    fig.colorbar(im, ax=ax, shrink=0.8, label="Identity %")

    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, f"{prefix}_identity_matrix.png")
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {os.path.basename(path)}")

    csv_path = os.path.join(OUTPUT_DIR, f"{prefix}_identity_matrix.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow([""] + labels)
        for i in range(n_seqs):
            writer.writerow([labels[i]] + [f"{matrix[i,j]:.1f}" for j in range(n_seqs)])
    print(f"  Saved {os.path.basename(csv_path)}")

    return path


# ── Domain conservation bar chart ─────────────────────────────────────────
def plot_domain_conservation(alignment, scores, col_to_human, human_to_col,
                              protein_name, prefix, l22e_start, l22e_end):
    n_cols = len(scores)
    overall_mean = float(np.mean(scores[scores > 0]))

    domain_cols = []
    for col in range(n_cols):
        h = col_to_human.get(col)
        if h is not None and l22e_start <= h <= l22e_end:
            domain_cols.append(col)
    domain_mean = float(np.mean(scores[domain_cols])) if domain_cols else 0

    nterm_cols = [c for c in range(n_cols)
                  if col_to_human.get(c) is not None and col_to_human[c] < l22e_start]
    nterm_mean = float(np.mean(scores[nterm_cols])) if nterm_cols else 0

    cterm_cols = [c for c in range(n_cols)
                  if col_to_human.get(c) is not None and col_to_human[c] > l22e_end]
    cterm_mean = float(np.mean(scores[cterm_cols])) if cterm_cols else 0

    regions = ["Overall", f"L22e domain\n({l22e_start}-{l22e_end})",
               f"N-terminal\n(1-{l22e_start-1})", f"C-terminal\n({l22e_end+1}+)"]
    values = [overall_mean, domain_mean, nterm_mean, cterm_mean]
    colors = ["#78909C", "#7E57C2", "#FFB74D", "#FFB74D"]

    valid = [(r, v, c) for r, v, c in zip(regions, values, colors) if v > 0]
    if not valid:
        return None

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar([r for r, _, _ in valid], [v for _, v, _ in valid],
                  color=[c for _, _, c in valid], edgecolor="white", width=0.6)
    for bar, val in zip(bars, [v for _, v, _ in valid]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                f"{val:.1%}", ha="center", fontsize=10, fontweight="bold")

    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Mean conservation score", fontsize=12)
    ax.set_title(f"{protein_name} Domain-Specific Conservation",
                 fontsize=13, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, f"{prefix}_domain_conservation.png")
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {os.path.basename(path)}")
    return path


# ── Conservation CSV ──────────────────────────────────────────────────────
def save_conservation_csv(alignment, scores, col_to_human, human_idx, prefix):
    n_cols = alignment.get_alignment_length()
    path = os.path.join(OUTPUT_DIR, f"{prefix}_conservation_scores.csv")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["alignment_col", "human_pos", "human_aa",
                          "consensus", "conservation_score"])
        for col in range(n_cols):
            h_pos = col_to_human.get(col, "")
            h_aa = alignment[human_idx, col]
            cons = get_consensus(alignment, col)
            writer.writerow([col, h_pos if h_pos else "", h_aa, cons,
                              f"{scores[col]:.4f}"])
    print(f"  Saved {os.path.basename(path)}")
    return path


# ── Conserved regions map ─────────────────────────────────────────────────
def plot_conserved_regions_map(scores, regions, col_to_human, protein_name,
                                prefix, human_len):
    fig, ax = plt.subplots(figsize=(14, 3))
    ax.barh(0, human_len, height=0.3, color="#E0E0E0", edgecolor="none")

    for rstart, rend, rmean in regions:
        h_start, h_end, _ = human_range_label(col_to_human, rstart, rend)
        if h_start and h_end:
            width = h_end - h_start + 1
            color = "#B71C1C" if rmean >= 1.0 else (
                "#E53935" if rmean >= 0.9 else "#FF5722")
            ax.barh(0, width, left=h_start, height=0.3,
                    color=color, edgecolor="none", alpha=0.9)

    ax.set_xlim(0, human_len + 5)
    ax.set_ylim(-0.5, 0.5)
    ax.set_yticks([])
    ax.set_xlabel(f"Human {protein_name} residue position", fontsize=12)
    ax.set_title(f"Conserved Regions Mapped onto Human {protein_name}",
                 fontsize=13, fontweight="bold")

    legend_items = [
        mpatches.Patch(color="#B71C1C", label="100% identical"),
        mpatches.Patch(color="#E53935", label="90-99%"),
        mpatches.Patch(color="#FF5722", label="70-90%"),
        mpatches.Patch(color="#E0E0E0", label="Not conserved (>=70%)"),
    ]
    ax.legend(handles=legend_items, fontsize=8, loc="upper right", ncol=4)

    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, f"{prefix}_conserved_regions_map.png")
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved {os.path.basename(path)}")
    return path


def analyze_protein(fasta_path, protein_name, prefix, l22e_start, l22e_end):
    print(f"\n{'='*60}")
    print(f"Analyzing {protein_name}")
    print(f"{'='*60}")

    alignment = AlignIO.read(fasta_path, "fasta")
    n_seqs = len(alignment)
    n_cols = alignment.get_alignment_length()
    print(f"  {n_seqs} sequences x {n_cols} columns")

    scores = compute_conservation(alignment)
    col_to_human, human_to_col, human_idx = alignment_col_to_human_pos(alignment)
    human_len = max(v for v in col_to_human.values() if v is not None)
    print(f"  Human sequence: {alignment[human_idx].id} ({human_len} aa)")

    regions_70 = find_conserved_regions(scores, threshold=0.7, min_len=3)
    print(f"  Found {len(regions_70)} conserved regions (>=70%)")

    print("\n--- Conservation line plot ---")
    plot_conservation_line(alignment, scores, col_to_human, protein_name, prefix)

    print("\n--- Conservation heatmap ---")
    plot_conservation_heatmap(alignment, scores, protein_name, prefix)

    print("\n--- Identity matrix ---")
    plot_identity_matrix(alignment, protein_name, prefix)

    print("\n--- Domain conservation ---")
    plot_domain_conservation(alignment, scores, col_to_human, human_to_col,
                              protein_name, prefix, l22e_start, l22e_end)

    print("\n--- Conserved regions map ---")
    plot_conserved_regions_map(scores, regions_70, col_to_human, protein_name,
                                prefix, human_len)

    print("\n--- Conservation CSV ---")
    save_conservation_csv(alignment, scores, col_to_human, human_idx, prefix)

    return alignment, scores, col_to_human, human_to_col, human_idx


def main():
    print("=" * 60)
    print("RPL22 / RPL22L1 Conservation Analysis")
    print("=" * 60)

    rpl22_fasta = os.path.join(OUTPUT_DIR, "rpl22_aligned.fasta")
    rpl22l1_fasta = os.path.join(OUTPUT_DIR, "rpl22l1_aligned.fasta")

    if os.path.exists(rpl22_fasta):
        analyze_protein(rpl22_fasta, "RPL22", "rpl22",
                         l22e_start=15, l22e_end=128)
    else:
        print(f"  SKIP RPL22: {rpl22_fasta} not found")

    if os.path.exists(rpl22l1_fasta):
        analyze_protein(rpl22l1_fasta, "RPL22L1", "rpl22l1",
                         l22e_start=16, l22e_end=122)
    else:
        print(f"  SKIP RPL22L1: {rpl22l1_fasta} not found")

    print(f"\n{'='*60}")
    print("DONE")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
