"""
Separate GO / pathway enrichment for the UBR3 co-IP best hits.

Runs g:Profiler (g:GOSt) INDEPENDENTLY on:
  - WT best hits   : significant interactors of wild-type UBR3   (q < 0.05 vs 293T)
  - RING best hits : significant interactors of the RING-mutant UBR3 (q < 0.05 vs 293T)

"Best hits" = the significant co-IP interactors for each genotype (all pass q<0.05
vs the 293T control and are strongly enriched). Each list is analysed on its own so
the enriched terms reflect that genotype alone, not the pooled union.

Outputs (in this folder):
  GO_WT_besthits.csv          full enrichment table, WT only
  GO_RING_besthits.csv        full enrichment table, RING only
  GO_separate_summary.json    counts + parameters
  GO_WT_vs_RING_separate.(png|pdf|svg)   side-by-side top-terms figure
"""
import json
import pathlib
import requests
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = pathlib.Path(__file__).parent
GPROFILER = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
SOURCES = ["GO:BP", "GO:CC", "GO:MF", "KEGG", "REAC", "CORUM"]

# ---------------------------------------------------------------- best hits
hits = pd.read_csv(HERE / "hits_table.csv")

def best_hits(hit_col, diff_col):
    """Significant interactors for one genotype, ranked by enrichment (best first)."""
    sub = hits[hits[hit_col] == True].copy()
    sub = sub.sort_values(diff_col, ascending=False)
    return sub["primary_gene"].dropna().tolist()

groups = {
    "WT":   best_hits("hit_UBR3wt",   "WT_293T_log2diff"),
    "RING": best_hits("hit_UBR3ring", "UBR3ring_293T_log2diff"),
}
for name, genes in groups.items():
    print(f"{name}: {len(genes)} best hits")

# ---------------------------------------------------------------- enrichment
def run_gost(genes):
    payload = {
        "organism": "hsapiens",
        "query": {"best_hits": genes},
        "sources": SOURCES,
        "user_threshold": 0.05,
        "significance_threshold_method": "g_SCS",  # g:Profiler's recommended correction
        "domain_scope": "annotated",               # whole annotated genome background
        "no_evidences": True,
        "ordered": False,
    }
    r = requests.post(GPROFILER, json=payload, timeout=120)
    r.raise_for_status()
    res = r.json()["result"]
    if not res:
        return pd.DataFrame()
    df = pd.DataFrame(res)
    keep = ["source", "native", "name", "p_value", "term_size",
            "query_size", "intersection_size", "precision", "recall"]
    df = df[keep].rename(columns={"native": "term_id", "name": "term_name"})
    df["neg_log10_p"] = -np.log10(df["p_value"])
    return df.sort_values(["source", "p_value"]).reset_index(drop=True)

results = {}
for name, genes in groups.items():
    df = run_gost(genes)
    results[name] = df
    out = HERE / f"GO_{name}_besthits.csv"
    df.to_csv(out, index=False)
    print(f"  {name}: {len(df)} enriched terms -> {out.name}")

# ---------------------------------------------------------------- summary
summary = {
    "genotypes": {
        name: {"n_best_hits": len(genes), "n_enriched_terms": len(results[name])}
        for name, genes in groups.items()
    },
    "best_hits_definition": "significant co-IP interactors per genotype (q<0.05 vs 293T)",
    "parameters": {
        "tool": "g:Profiler g:GOSt REST API",
        "organism": "hsapiens",
        "sources": SOURCES,
        "correction": "g_SCS",
        "user_threshold": 0.05,
        "background": "whole annotated genome",
        "analysed_separately": True,
    },
}
(HERE / "GO_separate_summary.json").write_text(json.dumps(summary, indent=2))

# ---------------------------------------------------------------- figure
SRC_COLOR = {
    "GO:BP": "#4C72B0", "GO:CC": "#55A868", "GO:MF": "#C44E52",
    "KEGG": "#8172B3", "REAC": "#CCB974", "CORUM": "#64B5CD",
}
TOP_N = 12

fig, axes = plt.subplots(1, 2, figsize=(15, 7), sharex=False)
for ax, name in zip(axes, ["WT", "RING"]):
    df = results[name]
    if df.empty:
        ax.set_title(f"UBR3 {name}: no enriched terms")
        ax.axis("off")
        continue
    top = df.sort_values("p_value").head(TOP_N).iloc[::-1]
    colors = [SRC_COLOR.get(s, "#999999") for s in top["source"]]
    labels = [f"{n[:46]}{'…' if len(n) > 46 else ''}" for n in top["term_name"]]
    ax.barh(range(len(top)), top["neg_log10_p"], color=colors, edgecolor="black", linewidth=0.4)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel(r"$-\log_{10}$ adjusted $p$")
    n = len(groups[name])
    ax.set_title(f"UBR3 {name} best hits (n={n})\n{len(df)} enriched terms",
                 fontsize=12, fontweight="bold")
    ax.axvline(-np.log10(0.05), color="grey", ls="--", lw=0.8)

handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in SRC_COLOR.values()]
fig.legend(handles, SRC_COLOR.keys(), loc="lower center", ncol=6,
           frameon=False, bbox_to_anchor=(0.5, -0.02))
fig.suptitle("GO / pathway enrichment — WT vs RING best hits (analysed separately)",
             fontsize=14, fontweight="bold")
fig.tight_layout(rect=[0, 0.03, 1, 0.96])
for ext in ("png", "pdf", "svg"):
    fig.savefig(HERE / f"GO_WT_vs_RING_separate.{ext}", dpi=200, bbox_inches="tight")
print("figure -> GO_WT_vs_RING_separate.png/pdf/svg")
