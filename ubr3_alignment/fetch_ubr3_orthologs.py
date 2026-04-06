"""
Fetch human UBR3 protein sequence and orthologs from NCBI,
perform multiple sequence alignment via EBI Clustal Omega,
and visualize conserved regions.
"""

import os
import sys
import time
import re
from io import StringIO
from collections import Counter
from urllib import request, parse, error

from Bio import SeqIO, AlignIO, Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

Entrez.email = "ubr3.alignment@example.com"
Entrez.tool  = "UBR3_ortholog_alignment"

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

# Full-length UBR3 orthologs verified from NCBI Gene orthologs page
# Human UBR3 = Gene ID 130507; ~1888 aa
UBR3_ORTHOLOGS = {
    "Human":       "NP_742067.3",
    "Mouse":       "NP_808451.2",
    "Rat":         "XP_006237133.1",
    "Cow":         "XP_005214862.1",
    "Dog":         "XP_005618987.1",
    "Chicken":     "XP_004935306.1",
    "Frog":        "XP_002935553.2",
    "Zebrafish":   "XP_021332594.1",
    "Fruitfly":    "NP_572428.2",
}


def fetch_sequences(accession_map: dict) -> list[SeqRecord]:
    """Fetch protein sequences from NCBI by accession."""
    records = []
    print(f"Fetching {len(accession_map)} UBR3 sequences from NCBI...")

    for name, acc in accession_map.items():
        for attempt in range(3):
            try:
                handle = Entrez.efetch(
                    db="protein", id=acc, rettype="fasta", retmode="text"
                )
                text = handle.read()
                handle.close()

                rec = SeqIO.read(StringIO(text), "fasta")
                rec.id = name
                rec.name = name
                rec.description = f"{name} | {acc}"
                records.append(rec)
                print(f"  OK  {name:12s}  {acc:22s}  ({len(rec.seq)} aa)")
                time.sleep(0.5)
                break
            except Exception as exc:
                print(f"  RETRY {name} ({acc}): {exc}")
                time.sleep(3)
        else:
            print(f"  FAIL {name} ({acc}) after 3 attempts - skipping")

    return records


def search_additional_orthologs() -> dict:
    """Search NCBI for additional UBR3 orthologs not in our curated list."""
    extra = {}
    query = '"UBR3"[Gene Name] AND refseq[filter] AND ("1500"[SLEN] : "3000"[SLEN])'
    try:
        handle = Entrez.esearch(db="protein", term=query, retmax=50)
        result = Entrez.read(handle)
        handle.close()
        ids = result.get("IdList", [])
        if ids:
            handle = Entrez.efetch(db="protein", id=",".join(ids[:30]),
                                   rettype="fasta", retmode="text")
            text = handle.read()
            handle.close()
            for rec in SeqIO.parse(StringIO(text), "fasta"):
                org_match = re.search(r'\[(.+?)\]', rec.description)
                if org_match and len(rec.seq) > 1400:
                    org = org_match.group(1)
                    short = org.split()[0][:10]
                    if short not in extra:
                        extra[short] = rec
            print(f"  Found {len(extra)} additional orthologs from NCBI search")
    except Exception as exc:
        print(f"  Additional ortholog search failed: {exc}")
    return extra


def save_fasta(records: list[SeqRecord], filename: str) -> str:
    path = os.path.join(OUTPUT_DIR, filename)
    SeqIO.write(records, path, "fasta")
    print(f"Saved {len(records)} sequences to {path}")
    return path


def run_clustalo_ebi(fasta_path: str) -> str:
    """Clustal Omega via EBI REST API."""
    url_run    = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run"
    url_status = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{}"
    url_result = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{}/aln-fasta"

    with open(fasta_path, encoding="utf-8") as fh:
        seq_data = fh.read()

    data = parse.urlencode({
        "email":    Entrez.email,
        "sequence": seq_data,
        "outfmt":   "fa",
    }).encode()

    print("Submitting Clustal Omega alignment to EBI...")
    req = request.Request(url_run, data=data)
    resp = request.urlopen(req)
    job_id = resp.read().decode().strip()
    print(f"  Job ID: {job_id}")

    for _ in range(180):
        time.sleep(10)
        resp = request.urlopen(url_status.format(job_id))
        status = resp.read().decode().strip()
        print(f"  Status: {status}")
        if status == "FINISHED":
            break
        if status in ("FAILURE", "ERROR", "NOT_FOUND"):
            raise RuntimeError(f"EBI Clustal Omega job failed: {status}")
    else:
        raise RuntimeError("EBI Clustal Omega job timed out after 30 min")

    resp = request.urlopen(url_result.format(job_id))
    aligned = resp.read().decode()
    return aligned


def run_muscle_ebi(fasta_path: str) -> str:
    """MUSCLE via EBI REST API."""
    url_run    = "https://www.ebi.ac.uk/Tools/services/rest/muscle/run"
    url_status = "https://www.ebi.ac.uk/Tools/services/rest/muscle/status/{}"
    url_result = "https://www.ebi.ac.uk/Tools/services/rest/muscle/result/{}/aln-fasta"

    with open(fasta_path, encoding="utf-8") as fh:
        seq_data = fh.read()

    data = parse.urlencode({
        "email":    Entrez.email,
        "sequence": seq_data,
        "format":   "fasta",
    }).encode()

    print("Submitting MUSCLE alignment to EBI...")
    req = request.Request(url_run, data=data)
    resp = request.urlopen(req)
    job_id = resp.read().decode().strip()
    print(f"  Job ID: {job_id}")

    for _ in range(180):
        time.sleep(10)
        resp = request.urlopen(url_status.format(job_id))
        status = resp.read().decode().strip()
        print(f"  Status: {status}")
        if status == "FINISHED":
            break
        if status in ("FAILURE", "ERROR", "NOT_FOUND"):
            raise RuntimeError(f"EBI MUSCLE job failed: {status}")
    else:
        raise RuntimeError("EBI MUSCLE job timed out after 30 min")

    resp = request.urlopen(url_result.format(job_id))
    aligned = resp.read().decode()
    return aligned


def align_sequences(fasta_path: str) -> str:
    """Try Clustal Omega first, fall back to MUSCLE."""
    out_path = os.path.join(OUTPUT_DIR, "ubr3_aligned.fasta")

    for aligner_fn, label in [
        (run_clustalo_ebi, "Clustal Omega"),
        (run_muscle_ebi, "MUSCLE"),
    ]:
        try:
            aligned_text = aligner_fn(fasta_path)
            with open(out_path, "w", encoding="utf-8") as fh:
                fh.write(aligned_text)
            print(f"Alignment ({label}) saved to {out_path}")
            return out_path
        except Exception as exc:
            print(f"  {label} failed: {exc} - trying next aligner")

    raise RuntimeError("All remote alignment services failed.")


# ── Conservation analysis ─────────────────────────────────────────────────
def compute_conservation(alignment) -> np.ndarray:
    """Per-column conservation: fraction of most-common residue (ignoring gaps)."""
    n_seqs = len(alignment)
    n_cols = alignment.get_alignment_length()
    scores = np.zeros(n_cols)

    for col in range(n_cols):
        residues = [alignment[row, col] for row in range(n_seqs)]
        non_gap = [r for r in residues if r != "-"]
        if not non_gap:
            scores[col] = 0.0
            continue
        counts = Counter(non_gap)
        scores[col] = counts.most_common(1)[0][1] / len(non_gap)

    return scores


def smooth(arr: np.ndarray, window: int = 30) -> np.ndarray:
    if window < 2:
        return arr
    kernel = np.ones(window) / window
    return np.convolve(arr, kernel, mode="same")


# ── Visualization ─────────────────────────────────────────────────────────
def plot_conservation_line(scores: np.ndarray, prefix: str):
    n_cols = len(scores)
    smoothed = smooth(scores, window=30)

    fig, ax = plt.subplots(figsize=(20, 5))
    ax.fill_between(range(n_cols), smoothed, alpha=0.35, color="#2196F3")
    ax.plot(range(n_cols), smoothed, linewidth=0.8, color="#0D47A1")
    ax.set_xlim(0, n_cols)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("Alignment position", fontsize=12)
    ax.set_ylabel("Conservation score", fontsize=12)
    ax.set_title("UBR3 Ortholog Conservation (smoothed, window=30)",
                 fontsize=14, fontweight="bold")
    ax.axhline(0.9, ls="--", color="red", alpha=0.6, lw=1.2, label="90% conserved")
    ax.axhline(0.7, ls="--", color="orange", alpha=0.6, lw=1.2, label="70% conserved")
    ax.legend(fontsize=10)
    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, f"{prefix}_conservation_line.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"  Saved {path}")
    return path


def plot_conservation_heatmap(alignment, scores: np.ndarray, prefix: str):
    cmap = LinearSegmentedColormap.from_list(
        "cons", ["#FFFFFF", "#BBDEFB", "#1565C0", "#0D47A1", "#B71C1C"]
    )
    n_seqs = len(alignment)
    n_cols = alignment.get_alignment_length()

    matrix = np.zeros((n_seqs, n_cols))
    for r in range(n_seqs):
        for c in range(n_cols):
            aa = alignment[r, c]
            matrix[r, c] = np.nan if aa == "-" else scores[c]

    fig, ax = plt.subplots(figsize=(22, max(4, n_seqs * 0.55)))
    labels = [rec.id for rec in alignment]
    im = ax.imshow(matrix, aspect="auto", cmap=cmap, vmin=0, vmax=1,
                   interpolation="none")
    ax.set_yticks(range(n_seqs))
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel("Alignment position", fontsize=12)
    ax.set_title("UBR3 Conservation Heatmap across Orthologs",
                 fontsize=14, fontweight="bold")
    fig.colorbar(im, ax=ax, shrink=0.6, label="Conservation score")
    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, f"{prefix}_conservation_heatmap.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"  Saved {path}")
    return path


def plot_identity_matrix(alignment, prefix: str):
    n = len(alignment)
    labels = [rec.id for rec in alignment]
    length = alignment.get_alignment_length()
    mat = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            matches = sum(
                1 for c in range(length)
                if alignment[i, c] == alignment[j, c] and alignment[i, c] != "-"
            )
            compared = sum(
                1 for c in range(length)
                if alignment[i, c] != "-" and alignment[j, c] != "-"
            )
            mat[i, j] = (matches / compared * 100) if compared else 0

    fig, ax = plt.subplots(figsize=(9, 8))
    im = ax.imshow(mat, cmap="YlOrRd", vmin=20, vmax=100)
    ax.set_xticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels, fontsize=9)
    for i in range(n):
        for j in range(n):
            ax.text(j, i, f"{mat[i,j]:.0f}", ha="center", va="center",
                    fontsize=7, color="white" if mat[i, j] > 65 else "black")
    ax.set_title("Pairwise Sequence Identity (%) - UBR3 Orthologs",
                 fontsize=13, fontweight="bold")
    fig.colorbar(im, ax=ax, shrink=0.7, label="% Identity")
    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, f"{prefix}_identity_matrix.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"  Saved {path}")

    csv_path = os.path.join(OUTPUT_DIR, f"{prefix}_identity_matrix.csv")
    with open(csv_path, "w", encoding="utf-8") as fh:
        fh.write("," + ",".join(labels) + "\n")
        for i, lab in enumerate(labels):
            fh.write(lab + "," + ",".join(f"{mat[i,j]:.1f}" for j in range(n)) + "\n")
    print(f"  Saved {csv_path}")
    return path


def plot_domain_conservation(alignment, scores: np.ndarray, prefix: str):
    """Conservation with approximate UBR3 domain annotations."""
    n_cols = len(scores)
    smoothed = smooth(scores, window=20)

    # Approximate domain positions based on human UBR3 (UniProt Q6ZT12)
    domains = [
        (30,  170,  "UBR box",      "#4CAF50"),
        (1540, 1620, "RING domain", "#9C27B0"),
    ]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(20, 8),
                                    gridspec_kw={"height_ratios": [3, 1]},
                                    sharex=True)

    ax1.fill_between(range(n_cols), smoothed, alpha=0.3, color="#2196F3")
    ax1.plot(range(n_cols), smoothed, linewidth=0.8, color="#0D47A1")
    ax1.axhline(0.9, ls="--", color="red", alpha=0.5, lw=1)
    ax1.axhline(0.7, ls="--", color="orange", alpha=0.5, lw=1)
    ax1.set_ylim(0, 1.05)
    ax1.set_ylabel("Conservation", fontsize=12)
    ax1.set_title("UBR3 Conservation with Domain Annotations",
                  fontsize=14, fontweight="bold")

    for start, end, name, color in domains:
        if end < n_cols:
            ax1.axvspan(start, end, alpha=0.15, color=color, label=name)
    ax1.legend(fontsize=10, loc="lower right")

    raw_scores = scores.copy()
    colors = np.where(raw_scores >= 0.9, "#B71C1C",
             np.where(raw_scores >= 0.7, "#FF9800",
             np.where(raw_scores >= 0.5, "#FFC107", "#E0E0E0")))
    ax2.bar(range(n_cols), np.ones(n_cols), color=colors, width=1.0, edgecolor="none")
    ax2.set_yticks([])
    ax2.set_xlabel("Alignment position", fontsize=12)
    ax2.set_ylabel("Region type", fontsize=10)

    import matplotlib.patches as mpatches
    legend_items = [
        mpatches.Patch(color="#B71C1C", label=">90% conserved"),
        mpatches.Patch(color="#FF9800", label="70-90% conserved"),
        mpatches.Patch(color="#FFC107", label="50-70% conserved"),
        mpatches.Patch(color="#E0E0E0", label="<50% conserved"),
    ]
    ax2.legend(handles=legend_items, fontsize=8, ncol=4, loc="upper center")

    ax1.set_xlim(0, n_cols)
    fig.tight_layout()
    path = os.path.join(OUTPUT_DIR, f"{prefix}_domain_conservation.png")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"  Saved {path}")
    return path


def find_conserved_regions(scores: np.ndarray, threshold: float = 0.9,
                           min_len: int = 5) -> list[tuple]:
    """Find stretches of high conservation."""
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
                        float(np.mean(scores[start:len(scores)]))))
    return regions


def save_conserved_regions(scores: np.ndarray, prefix: str):
    csv_path = os.path.join(OUTPUT_DIR, f"{prefix}_conserved_regions.csv")
    regions = find_conserved_regions(scores, threshold=0.9, min_len=5)
    with open(csv_path, "w", encoding="utf-8") as fh:
        fh.write("region_start,region_end,length,mean_conservation\n")
        for s, e, m in regions:
            fh.write(f"{s},{e},{e - s + 1},{m:.4f}\n")
    print(f"  Found {len(regions)} highly conserved regions (>=90%, >=5 aa)")
    print(f"  Saved {csv_path}")

    regions70 = find_conserved_regions(scores, threshold=0.7, min_len=10)
    csv70 = os.path.join(OUTPUT_DIR, f"{prefix}_conserved_regions_70pct.csv")
    with open(csv70, "w", encoding="utf-8") as fh:
        fh.write("region_start,region_end,length,mean_conservation\n")
        for s, e, m in regions70:
            fh.write(f"{s},{e},{e - s + 1},{m:.4f}\n")
    print(f"  Found {len(regions70)} conserved regions (>=70%, >=10 aa)")
    print(f"  Saved {csv70}")


def print_summary(alignment, scores: np.ndarray):
    n_cols = len(scores)
    print(f"\n{'='*60}")
    print(f"ALIGNMENT SUMMARY")
    print(f"{'='*60}")
    print(f"  Sequences:              {len(alignment)}")
    print(f"  Alignment length:       {n_cols} columns")
    print(f"  Mean conservation:      {np.mean(scores):.3f}")
    print(f"  Positions >=90% cons.:  {np.sum(scores >= 0.9):,} / {n_cols:,} "
          f"({np.sum(scores >= 0.9)/n_cols*100:.1f}%)")
    print(f"  Positions 100% cons.:   {np.sum(scores >= 1.0):,} / {n_cols:,} "
          f"({np.sum(scores >= 1.0)/n_cols*100:.1f}%)")
    print(f"  Positions >=70% cons.:  {np.sum(scores >= 0.7):,} / {n_cols:,} "
          f"({np.sum(scores >= 0.7)/n_cols*100:.1f}%)")

    species = [rec.id for rec in alignment]
    print(f"\n  Species included:")
    for sp in species:
        print(f"    - {sp}")
    print(f"{'='*60}")


# ── Main ──────────────────────────────────────────────────────────────────
def main():
    print("=" * 60)
    print("UBR3 Ortholog Alignment & Conservation Analysis")
    print("=" * 60)

    # Step 1: Fetch sequences
    records = fetch_sequences(UBR3_ORTHOLOGS)

    # Try to find additional orthologs
    extra = search_additional_orthologs()
    existing_ids = {r.id for r in records}
    for short_name, rec in extra.items():
        if short_name not in existing_ids and len(records) < 20:
            rec.id = short_name
            rec.name = short_name
            records.append(rec)
            print(f"  Added extra: {short_name} ({len(rec.seq)} aa)")

    # Filter to full-length sequences only (UBR3 is ~1800-2200 aa)
    min_len = 1400
    full_records = [r for r in records if len(r.seq) >= min_len]
    partial = [r for r in records if len(r.seq) < min_len]
    if partial:
        print(f"\n  Excluded {len(partial)} partial/predicted isoforms (<{min_len} aa):")
        for r in partial:
            print(f"    - {r.id} ({len(r.seq)} aa)")
    records = full_records

    if len(records) < 2:
        sys.exit("Need at least 2 sequences for alignment.")

    fasta_path = save_fasta(records, "ubr3_orthologs.fasta")

    # Step 2: Align
    print(f"\nStarting alignment of {len(records)} sequences...")
    aligned_path = align_sequences(fasta_path)

    # Step 3: Parse alignment
    alignment = AlignIO.read(aligned_path, "fasta")
    print(f"Alignment: {len(alignment)} sequences x "
          f"{alignment.get_alignment_length()} columns")

    # Step 4: Conservation
    print("\nComputing conservation scores...")
    scores = compute_conservation(alignment)

    # Step 5: Summary
    print_summary(alignment, scores)

    # Step 6: Plots & outputs
    print("\nGenerating visualizations...")
    prefix = "ubr3"
    plot_conservation_line(scores, prefix)
    plot_conservation_heatmap(alignment, scores, prefix)
    plot_identity_matrix(alignment, prefix)
    plot_domain_conservation(alignment, scores, prefix)
    save_conserved_regions(scores, prefix)

    # Save raw conservation scores
    np.savetxt(
        os.path.join(OUTPUT_DIR, f"{prefix}_conservation_scores.csv"),
        scores, delimiter=",", header="conservation_score", comments=""
    )
    print(f"  Saved {prefix}_conservation_scores.csv")

    print(f"\n{'='*60}")
    print(f"DONE - all outputs in: {OUTPUT_DIR}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
