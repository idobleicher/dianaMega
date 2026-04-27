"""
Shared plotting / counting helpers for the MP-project enrichment analyses.

Used by:
  - mp_enrichment_heatmap.py        (vs MP-only library)
  - mp_vs_full_library.py           (vs full library, all N-termini)
  - mp_extra_plots.py               (companion plots, parameterised by bg)
"""
from __future__ import annotations

from pathlib import Path

import logomaker
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import fisher_exact

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")

AA_GROUPS = {
    "Hydrophobic": list("AILMV"),
    "Aromatic":    list("FWY"),
    "Polar":       list("NQST"),
    "Basic":       list("KRH"),
    "Acidic":      list("DE"),
    "Special":     list("CGP"),
}
GROUP_COLORS = {
    "Hydrophobic": "#4C8C40",
    "Aromatic":    "#9C27B0",
    "Polar":       "#1F77B4",
    "Basic":       "#1565C0",
    "Acidic":      "#D32F2F",
    "Special":     "#F57C00",
}
PSEUDOCOUNT = 0.5


# ---------------------------------------------------------------------------
# Counts / frequencies / enrichment
# ---------------------------------------------------------------------------
def count_matrix(seqs, positions):
    counts = np.zeros((len(AA_ORDER), len(positions)), dtype=float)
    aa_idx = {a: i for i, a in enumerate(AA_ORDER)}
    for s in seqs:
        for j, p in enumerate(positions):
            if len(s) >= p:
                a = s[p - 1]
                if a in aa_idx:
                    counts[aa_idx[a], j] += 1
    return counts


def freq_matrix(counts):
    totals = counts.sum(axis=0, keepdims=True)
    totals[totals == 0] = 1.0
    return counts / totals


def log2_enrichment(hit_seqs, lib_seqs, positions):
    """Returns (enrichment_df, hit_counts_df, lib_counts_df).

    Cells where BOTH the hit and library counts are 0 are forced to 0
    (no signal), to avoid pseudocount artefacts in fully empty cells.
    """
    h = count_matrix(hit_seqs, positions)
    l = count_matrix(lib_seqs, positions)
    h_freq = (h + PSEUDOCOUNT) / (h + PSEUDOCOUNT).sum(axis=0, keepdims=True)
    l_freq = (l + PSEUDOCOUNT) / (l + PSEUDOCOUNT).sum(axis=0, keepdims=True)
    enrich = np.log2(h_freq / l_freq)
    enrich[(h == 0) & (l == 0)] = 0.0
    cols = [str(p) for p in positions]
    return (pd.DataFrame(enrich, index=AA_ORDER, columns=cols),
            pd.DataFrame(h, index=AA_ORDER, columns=cols),
            pd.DataFrame(l, index=AA_ORDER, columns=cols))


# ---------------------------------------------------------------------------
# Plot 1: enrichment heatmap
# ---------------------------------------------------------------------------
def plot_enrichment_heatmap(df: pd.DataFrame, title: str, out_path: Path,
                            vmax: float | None = None):
    if vmax is None:
        vmax = max(float(np.nanpercentile(np.abs(df.values), 99)), 0.5)
    fig, ax = plt.subplots(figsize=(13, 8))
    sns.heatmap(
        df, cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax,
        linewidths=0.4, linecolor="white",
        cbar_kws={"label": "log2(hit / library)", "shrink": 0.8},
        ax=ax, xticklabels=True, yticklabels=True,
    )
    ax.set_xlabel("Residue position", fontsize=13, weight="bold")
    ax.set_ylabel("Amino acid", fontsize=13, weight="bold")
    ax.set_title(title, fontsize=14, weight="bold", pad=14)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Plot 2: side-by-side frequency heatmaps
# ---------------------------------------------------------------------------
def plot_frequency_panels(hit_seqs, lib_seqs, positions, out_path,
                          lib_label: str):
    h = count_matrix(hit_seqs, positions)
    l = count_matrix(lib_seqs, positions)
    h_freq = freq_matrix(h)
    l_freq = freq_matrix(l)
    cols = [str(p) for p in positions]
    df_h = pd.DataFrame(h_freq * 100, index=AA_ORDER, columns=cols)
    df_l = pd.DataFrame(l_freq * 100, index=AA_ORDER, columns=cols)
    vmax = max(df_h.values.max(), df_l.values.max())

    fig, axes = plt.subplots(1, 2, figsize=(24, 8), sharey=True)
    for ax, df, sub in (
        (axes[0], df_h, f"A   MP project hits  (n={len(hit_seqs)})"),
        (axes[1], df_l, f"B   {lib_label}  (n={len(lib_seqs)})"),
    ):
        sns.heatmap(
            df, cmap="YlOrRd", vmin=0, vmax=vmax,
            linewidths=0.4, linecolor="white",
            cbar_kws={"label": "Frequency (%)", "shrink": 0.8},
            ax=ax, xticklabels=True, yticklabels=True,
        )
        ax.set_xlabel("Residue position", fontsize=13, weight="bold")
        ax.set_ylabel("Amino acid", fontsize=13, weight="bold")
        ax.set_title(sub, loc="left", fontsize=13, weight="bold", pad=10)
    fig.suptitle("Amino-acid frequencies per position", fontsize=16,
                 weight="bold", y=1.02)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Plot 3: Fisher's exact significance heatmap
# ---------------------------------------------------------------------------
def plot_significance(hit_seqs, lib_seqs, positions, out_path,
                      out_csv: Path, title_lib: str):
    h = count_matrix(hit_seqs, positions)
    l = count_matrix(lib_seqs, positions)
    h_tot = h.sum(axis=0)
    l_tot = l.sum(axis=0)

    sig = np.zeros_like(h)
    pvals = np.ones_like(h)
    for i in range(h.shape[0]):
        for j in range(h.shape[1]):
            a = int(h[i, j])
            b = int(h_tot[j] - a)
            c = int(l[i, j])
            d = int(l_tot[j] - c)
            if (a + c) == 0 or (a + b) == 0 or (c + d) == 0:
                continue
            try:
                _, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
            except ValueError:
                p = 1.0
            pvals[i, j] = p
            obs = a / max(a + b, 1)
            exp = c / max(c + d, 1)
            sig[i, j] = (1.0 if obs >= exp else -1.0) * (-np.log10(max(p, 1e-300)))

    cols = [str(p) for p in positions]
    df_sig = pd.DataFrame(sig, index=AA_ORDER, columns=cols)
    df_sig.to_csv(out_csv)

    vmax = max(float(np.nanpercentile(np.abs(sig), 99)), 1.0)
    annot = np.full(sig.shape, "", dtype=object)
    annot[(pvals < 0.05) & (pvals >= 0.01)] = "*"
    annot[(pvals < 0.01) & (pvals >= 0.001)] = "**"
    annot[pvals < 0.001] = "***"

    fig, ax = plt.subplots(figsize=(13, 8))
    sns.heatmap(
        df_sig, cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax,
        linewidths=0.4, linecolor="white",
        annot=annot, fmt="", annot_kws={"size": 8, "weight": "bold"},
        cbar_kws={"label": "signed -log10(Fisher p)", "shrink": 0.8},
        ax=ax, xticklabels=True, yticklabels=True,
    )
    ax.set_xlabel("Residue position", fontsize=13, weight="bold")
    ax.set_ylabel("Amino acid", fontsize=13, weight="bold")
    ax.set_title(
        f"Statistical significance — Fisher's exact, MP project (n={len(hit_seqs)}) "
        f"vs {title_lib} (n={len(lib_seqs)})\n"
        "Cell text: *p<0.05  **p<0.01  ***p<0.001",
        fontsize=13, weight="bold", pad=12,
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Plot 4: enrichment sequence logo
# ---------------------------------------------------------------------------
def plot_logo(hit_seqs, lib_seqs, positions, out_path, title_lib: str):
    h_raw = count_matrix(hit_seqs, positions)
    l_raw = count_matrix(lib_seqs, positions)
    h = h_raw + PSEUDOCOUNT
    l = l_raw + PSEUDOCOUNT
    h_freq = h / h.sum(axis=0, keepdims=True)
    l_freq = l / l.sum(axis=0, keepdims=True)
    enrich = np.log2(h_freq / l_freq)
    enrich[(h_raw == 0) & (l_raw == 0)] = 0.0
    df = pd.DataFrame(enrich.T, index=positions, columns=AA_ORDER)

    color_scheme = {}
    for grp, aas in AA_GROUPS.items():
        for aa in aas:
            color_scheme[aa] = GROUP_COLORS[grp]

    fig, ax = plt.subplots(figsize=(15, 5))
    logo = logomaker.Logo(
        df, ax=ax, color_scheme=color_scheme, shade_below=0.5, fade_below=0.5,
        flip_below=True,
    )
    logo.style_xticks(anchor=positions[0], spacing=1)
    logo.ax.set_xlabel("Residue position", fontsize=13, weight="bold")
    logo.ax.set_ylabel("log2(hit / library)", fontsize=13, weight="bold")
    logo.ax.set_title(
        f"Enrichment logo — MP project (n={len(hit_seqs)}) vs {title_lib} "
        f"(n={len(lib_seqs)})\nLetter height = log2 enrichment "
        "(positive above zero, negative flipped below)",
        fontsize=13, weight="bold", pad=10,
    )
    logo.ax.axhline(0, color="black", linewidth=0.7)

    handles = [plt.Rectangle((0, 0), 1, 1, color=GROUP_COLORS[g])
               for g in AA_GROUPS]
    labels = [f"{g} ({''.join(AA_GROUPS[g])})" for g in AA_GROUPS]
    logo.ax.legend(handles, labels, loc="upper right", fontsize=8,
                   ncol=2, frameon=True)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Plot 5: chemistry-group heatmap
# ---------------------------------------------------------------------------
def plot_grouped_heatmap(hit_seqs, lib_seqs, positions, out_path, out_csv: Path,
                         title_lib: str):
    cols = [str(p) for p in positions]
    h = count_matrix(hit_seqs, positions)
    l = count_matrix(lib_seqs, positions)
    h_tot = h.sum(axis=0, keepdims=True)
    l_tot = l.sum(axis=0, keepdims=True)
    rows, labels = [], []
    for grp, aas in AA_GROUPS.items():
        idx = [AA_ORDER.index(a) for a in aas]
        h_grp = h[idx, :].sum(axis=0)
        l_grp = l[idx, :].sum(axis=0)
        h_freq = (h_grp + PSEUDOCOUNT) / (h_tot[0] + PSEUDOCOUNT * len(AA_ORDER))
        l_freq = (l_grp + PSEUDOCOUNT) / (l_tot[0] + PSEUDOCOUNT * len(AA_ORDER))
        rows.append(np.log2(h_freq / l_freq))
        labels.append(f"{grp} ({''.join(aas)})")
    df = pd.DataFrame(np.array(rows), index=labels, columns=cols)
    df.to_csv(out_csv)

    vmax = max(float(np.nanpercentile(np.abs(df.values), 99)), 0.5)
    fig, ax = plt.subplots(figsize=(13, 4))
    sns.heatmap(
        df, cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax,
        linewidths=0.4, linecolor="white",
        annot=True, fmt="+.2f", annot_kws={"size": 8},
        cbar_kws={"label": "log2(hit / library)", "shrink": 0.8},
        ax=ax, xticklabels=True, yticklabels=True,
    )
    ax.set_xlabel("Residue position", fontsize=13, weight="bold")
    ax.set_ylabel("AA chemistry group", fontsize=13, weight="bold")
    ax.set_title(
        f"Property-grouped enrichment — MP project (n={len(hit_seqs)}) "
        f"vs {title_lib} (n={len(lib_seqs)})",
        fontsize=13, weight="bold", pad=12,
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Plot 6: per-position top-mover bars
# ---------------------------------------------------------------------------
def plot_top_movers(hit_seqs, lib_seqs, positions, out_path, title_lib: str):
    h_raw = count_matrix(hit_seqs, positions)
    l_raw = count_matrix(lib_seqs, positions)
    h = h_raw + PSEUDOCOUNT
    l = l_raw + PSEUDOCOUNT
    h_freq = h / h.sum(axis=0, keepdims=True)
    l_freq = l / l.sum(axis=0, keepdims=True)
    enrich = np.log2(h_freq / l_freq)
    enrich[(h_raw == 0) & (l_raw == 0)] = 0.0

    aa_to_color = {}
    for grp, aas in AA_GROUPS.items():
        for aa in aas:
            aa_to_color[aa] = GROUP_COLORS[grp]

    n = len(positions)
    fig, axes = plt.subplots(1, n, figsize=(3.2 * n, 6.4))
    if n == 1:
        axes = [axes]

    all_vals = []
    for j in range(n):
        col = enrich[:, j]
        order = np.argsort(col)
        sel = list(order[-4:]) + list(order[:4])
        all_vals.extend(col[sel])
    xpad = max(0.6, max(abs(min(all_vals)), abs(max(all_vals))) * 0.30)
    xlim = (min(all_vals) - xpad, max(all_vals) + xpad)

    for j, (ax, p) in enumerate(zip(axes, positions)):
        col = enrich[:, j]
        order = np.argsort(col)
        sel = list(order[-4:][::-1]) + list(order[:4])
        labels = [AA_ORDER[i] for i in sel]
        values = [col[i] for i in sel]
        colors = [aa_to_color[a] for a in labels]
        ypos = np.arange(len(sel))
        ax.barh(ypos, values, color=colors, edgecolor="black", linewidth=0.6)
        ax.set_yticks(ypos)
        ax.set_yticklabels(labels, fontsize=14, weight="bold")
        ax.invert_yaxis()
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_xlim(*xlim)
        ax.set_title(f"Position {p}", fontsize=12, weight="bold")
        ax.set_xlabel("log2(hit / library)", fontsize=10)
        ax.grid(axis="x", linestyle=":", alpha=0.55)
        ax.tick_params(axis="x", labelsize=9)
        for k, v in enumerate(values):
            offset = xpad * 0.25
            tx, ha = (v + offset, "left") if v >= 0 else (v - offset, "right")
            ax.text(tx, k, f"{v:+.2f}", va="center", ha=ha,
                    fontsize=8.5, color="#222222")
    axes[0].set_ylabel("Amino acid", fontsize=12, weight="bold")

    fig.suptitle(
        f"Top enriched (above 0) & depleted (below 0) residues per position\n"
        f"MP project (n={len(hit_seqs)}) vs {title_lib} (n={len(lib_seqs)})",
        fontsize=13, weight="bold", y=1.02,
    )
    handles = [plt.Rectangle((0, 0), 1, 1, color=GROUP_COLORS[g])
               for g in AA_GROUPS]
    labels = [f"{g}" for g in AA_GROUPS]
    fig.legend(handles, labels, loc="lower center", ncol=len(AA_GROUPS),
               bbox_to_anchor=(0.5, -0.05), frameon=True, fontsize=9)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# One-stop: run all six plots for a given (hits, library, output_dir)
# ---------------------------------------------------------------------------
def run_full_panel_set(hit_seqs, lib_seqs, out_dir: Path,
                       lib_label: str, lib_short: str,
                       positions_main: list[int],
                       positions_top: list[int]):
    out_dir.mkdir(parents=True, exist_ok=True)
    enrich, h_counts, l_counts = log2_enrichment(
        hit_seqs, lib_seqs, positions_main)
    enrich.to_csv(out_dir / "enrichment.csv")
    h_counts.to_csv(out_dir / "counts_hits.csv")
    l_counts.to_csv(out_dir / "counts_library.csv")

    plot_enrichment_heatmap(
        enrich,
        f"MP project (n={len(hit_seqs)}) vs {lib_label} (n={len(lib_seqs)})\n"
        f"Log2 enrichment per residue (positions {positions_main[0]}–{positions_main[-1]})",
        out_dir / "mp_enrichment_heatmap.png",
    )
    plot_frequency_panels(hit_seqs, lib_seqs, positions_main,
                          out_dir / "mp_frequency_heatmaps.png", lib_label)
    plot_significance(hit_seqs, lib_seqs, positions_main,
                      out_dir / "mp_significance_heatmap.png",
                      out_dir / "significance_signed_log10p.csv",
                      lib_short)
    plot_logo(hit_seqs, lib_seqs, positions_main,
              out_dir / "mp_enrichment_logo.png", lib_short)
    plot_grouped_heatmap(hit_seqs, lib_seqs, positions_main,
                         out_dir / "mp_property_group_heatmap.png",
                         out_dir / "enrichment_property_groups.csv",
                         lib_short)
    plot_top_movers(hit_seqs, lib_seqs, positions_top,
                    out_dir / "mp_top_movers.png", lib_short)
    return enrich
