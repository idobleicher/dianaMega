"""
Mutual (shared) pathways among the UBR3-WT best hits.

Question: do the best hits of wild-type UBR3 converge on common pathways, or is
each enriched term driven by a different protein?

Approach
  1. Re-run g:Profiler on the WT best hits WITH evidences, so every enriched term
     comes back with the query genes that support it (the earlier run in
     go_analysis_separate.py used no_evidences=True and discarded that).
  2. Term x gene membership -> which hits are shared between which terms.
  3. Term-term graph (edge = gene-set overlap) -> modules of mutually connected
     pathways, each named by its most significant member term.
  4. Gene-gene graph (edge weight = number of specific terms two hits share)
     -> the groups of hits that travel together through those pathways.
  5. WT vs RING term comparison -> which WT pathways are also enriched in the
     RING mutant (mutual across genotypes) and which are WT-exclusive.

Outputs (this folder, prefix WT_mutual_):
  WT_mutual_terms.csv           59 WT terms + supporting genes + module
  WT_mutual_gene_matrix.csv     gene x term binary membership matrix
  WT_mutual_genes.csv           per hit: #terms, modules, partners, log2diff
  WT_mutual_modules.csv         pathway module -> core genes, member terms
  WT_mutual_vs_RING.csv         WT terms shared with / exclusive of RING
  WT_mutual_summary.json        counts + parameters
"""
import json
import pathlib
import itertools
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
import requests
import networkx as nx

HERE = pathlib.Path(__file__).parent
GPROFILER = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
SOURCES = ["GO:BP", "GO:CC", "GO:MF", "KEGG", "REAC", "CORUM"]

# A term is "specific" if it is small enough to describe a real pathway/complex
# rather than a whole-cell compartment ("nucleus", "protein binding").
SPECIFIC_MAX_SIZE = 500
# term-term edge: overlap coefficient = |A n B| / min(|A|,|B|)
TERM_EDGE_OVERLAP = 0.5
# gene-gene edge: minimum number of shared specific terms
GENE_EDGE_MIN_SHARED = 2

# ------------------------------------------------------------------ best hits
hits = pd.read_csv(HERE / "hits_table.csv")


def best_hits(hit_col, diff_col):
    sub = hits[hits[hit_col] == True].copy()
    sub = sub.sort_values(diff_col, ascending=False)
    return sub["primary_gene"].dropna().tolist()


wt_genes = best_hits("hit_UBR3wt", "WT_293T_log2diff")
ring_genes = best_hits("hit_UBR3ring", "UBR3ring_293T_log2diff")
wt_rank = {g: i + 1 for i, g in enumerate(wt_genes)}
wt_diff = dict(zip(hits["primary_gene"], hits["WT_293T_log2diff"]))

roles = pd.read_csv(HERE / "interactor_roles.csv")
roles.columns = [c.strip().lstrip("﻿") for c in roles.columns]
role_of = dict(zip(roles["gene"], roles["cellular_role"]))


# ---------------------------------------------------------------- enrichment
def run_gost(genes, label):
    """Enrichment with per-term member genes (query genes supporting the term)."""
    payload = {
        "organism": "hsapiens",
        "query": {label: genes},
        "sources": SOURCES,
        "user_threshold": 0.05,
        "significance_threshold_method": "g_SCS",
        "domain_scope": "annotated",
        "no_evidences": False,
        "ordered": False,
    }
    r = requests.post(GPROFILER, json=payload, timeout=180)
    r.raise_for_status()
    res = r.json()
    if not res["result"]:
        return pd.DataFrame(), {}

    mapping = res["meta"]["genes_metadata"]["query"][label]["mapping"]
    ensgs = res["meta"]["genes_metadata"]["query"][label]["ensgs"]
    sym_of = {}
    for sym, ids in mapping.items():
        for e in ids:
            sym_of[e] = sym
    order = [sym_of.get(e, e) for e in ensgs]

    rows, members = [], {}
    for t in res["result"]:
        inter = t.get("intersections") or []
        if len(inter) != len(order):
            raise RuntimeError(
                f"intersections length {len(inter)} != query length {len(order)}"
            )
        genes_in = [g for g, ev in zip(order, inter) if ev]
        members[t["native"]] = genes_in
        rows.append({
            "source": t["source"],
            "term_id": t["native"],
            "term_name": t["name"],
            "p_value": t["p_value"],
            "term_size": t["term_size"],
            "query_size": t["query_size"],
            "intersection_size": t["intersection_size"],
            "precision": t["precision"],
            "recall": t["recall"],
            "neg_log10_p": -np.log10(t["p_value"]),
            "genes": ";".join(genes_in),
            "n_genes": len(genes_in),
        })
    df = pd.DataFrame(rows).sort_values("p_value").reset_index(drop=True)
    return df, members


wt_df, wt_members = run_gost(wt_genes, "WT")
ring_df, ring_members = run_gost(ring_genes, "RING")
print(f"WT: {len(wt_genes)} best hits -> {len(wt_df)} enriched terms")
print(f"RING: {len(ring_genes)} best hits -> {len(ring_df)} enriched terms")

wt_df["specific"] = wt_df["term_size"] <= SPECIFIC_MAX_SIZE
specific = wt_df[wt_df["specific"]].copy()
print(f"  {len(specific)} of {len(wt_df)} WT terms are specific (term_size <= {SPECIFIC_MAX_SIZE})")

# ------------------------------------------------- histone leave-out control
# 7 histone/histone-variant genes are in the WT hit list. Histones are annotated
# to a very large number of Reactome "gene expression" pathways, so they can
# create pathway overlap between hits that share nothing else. Re-run without
# them: a term that disappears was carried by the histones, not by a genuine
# multi-protein module.
HISTONES = [g for g in wt_genes
            if g.startswith(("H1", "H2A", "H2B", "H3", "H4")) or g == "MACROH2A1"]
no_hist_df, _ = run_gost([g for g in wt_genes if g not in HISTONES], "WT_noHist")
survives = set(no_hist_df["term_id"]) if len(no_hist_df) else set()
wt_df["survives_without_histones"] = wt_df["term_id"].isin(survives)
wt_df["histone_driven"] = ~wt_df["survives_without_histones"]
print(f"  histones in WT hits ({len(HISTONES)}): {', '.join(HISTONES)}")
print(f"  {int(wt_df['survives_without_histones'].sum())}/{len(wt_df)} terms survive their removal")

# ------------------------------------------------- modules of mutual pathways
# Built on the SPECIFIC terms only: the broad terms ("nucleus", "protein
# binding") contain nearly every hit and would bridge all modules into one blob.
G_terms = nx.Graph()
for _, r in specific.iterrows():
    G_terms.add_node(r["term_id"], name=r["term_name"], source=r["source"],
                     p=r["p_value"], size=r["term_size"])

for a, b in itertools.combinations(specific["term_id"], 2):
    A, B = set(wt_members[a]), set(wt_members[b])
    if not A or not B:
        continue
    inter = A & B
    ov = len(inter) / min(len(A), len(B))
    if ov >= TERM_EDGE_OVERLAP:
        G_terms.add_edge(a, b, overlap=ov, shared=len(inter),
                         jaccard=len(inter) / len(A | B))

comms = nx.community.greedy_modularity_communities(G_terms, weight="overlap")
p_of = dict(zip(wt_df["term_id"], wt_df["p_value"]))
name_of = dict(zip(wt_df["term_id"], wt_df["term_name"]))
size_of = dict(zip(wt_df["term_id"], wt_df["term_size"]))

# order modules by best p-value inside; name by the most significant member term
comms = sorted(comms, key=lambda c: min(p_of[t] for t in c))
surv_of = dict(zip(wt_df["term_id"], wt_df["survives_without_histones"]))
module_of, module_rows = {}, []
for i, c in enumerate(comms, start=1):
    terms = sorted(c, key=lambda t: p_of[t])
    gene_counts = Counter()
    for t in terms:
        gene_counts.update(set(wt_members[t]))
    # core genes = present in at least half of the module's terms
    core = [g for g, n in gene_counts.most_common() if n >= max(2, len(terms) / 2)]
    for t in terms:
        module_of[t] = i
    n_surv = sum(bool(surv_of[t]) for t in terms)
    module_rows.append({
        "module": i,
        "module_label": name_of[terms[0]],
        "n_terms": len(terms),
        "n_terms_surviving_histone_removal": n_surv,
        "robust": n_surv > 0,
        "best_p": min(p_of[t] for t in terms),
        "n_genes": len(gene_counts),
        "core_genes": ";".join(core),
        "all_genes": ";".join(sorted(gene_counts)),
        "terms": ";".join(f"{name_of[t]} [{t}]" for t in terms),
    })

wt_df["module"] = wt_df["term_id"].map(module_of).fillna(0).astype(int)
wt_df["module_label"] = wt_df["module"].map(
    {r["module"]: r["module_label"] for r in module_rows}).fillna("(broad term - not clustered)")
modules = pd.DataFrame(module_rows)

# --------------------------------------------------- gene x pathway matrix
matrix = pd.DataFrame(
    0, index=wt_genes,
    columns=[f"{r.term_name} [{r.term_id}]" for r in wt_df.itertuples()], dtype=int)
for r in wt_df.itertuples():
    col = f"{r.term_name} [{r.term_id}]"
    for g in wt_members[r.term_id]:
        if g in matrix.index:
            matrix.at[g, col] = 1

spec_members = {t: set(wt_members[t]) for t in specific["term_id"]}

# gene-gene co-membership over specific terms only
shared_counts = defaultdict(int)
shared_terms = defaultdict(list)
for t, mem in spec_members.items():
    for a, b in itertools.combinations(sorted(mem), 2):
        shared_counts[(a, b)] += 1
        shared_terms[(a, b)].append(name_of[t])

G_genes = nx.Graph()
G_genes.add_nodes_from(wt_genes)
for (a, b), n in shared_counts.items():
    if n >= GENE_EDGE_MIN_SHARED:
        G_genes.add_edge(a, b, weight=n, terms=";".join(shared_terms[(a, b)]))

gene_comms = nx.community.greedy_modularity_communities(
    G_genes.subgraph([n for n in G_genes if G_genes.degree(n) > 0]), weight="weight") \
    if G_genes.number_of_edges() else []
gene_module_of = {}
for i, c in enumerate(sorted(gene_comms, key=len, reverse=True), start=1):
    for g in c:
        gene_module_of[g] = i

gene_rows = []
for g in wt_genes:
    in_terms = [r.term_id for r in wt_df.itertuples() if g in wt_members[r.term_id]]
    in_spec = [t for t in in_terms if t in spec_members]
    mods = sorted({module_of[t] for t in in_spec})
    partners = sorted(
        {b if a == g else a for (a, b), n in shared_counts.items()
         if g in (a, b) and n >= GENE_EDGE_MIN_SHARED})
    gene_rows.append({
        "gene": g,
        "wt_rank": wt_rank[g],
        "WT_293T_log2diff": wt_diff.get(g, np.nan),
        "cellular_role": role_of.get(g, ""),
        "n_terms": len(in_terms),
        "n_specific_terms": len(in_spec),
        "n_pathway_modules": len(mods),
        "pathway_modules": ";".join(map(str, mods)),
        "gene_module": gene_module_of.get(g, 0),
        "n_mutual_partners": len(partners),
        "mutual_partners": ";".join(partners),
        "specific_terms": ";".join(name_of[t] for t in in_spec),
    })
genes_df = pd.DataFrame(gene_rows).sort_values(
    ["n_specific_terms", "n_terms"], ascending=False).reset_index(drop=True)

# ------------------------------------------------- STRING cross-check
# Annotation databases can group two proteins that were never observed together.
# STRING (already computed for this dataset) is independent evidence: if a module
# is real, its members should also be physically/functionally linked there.
string_edges = pd.read_csv(HERE / "STRING_edges.csv")
string_clusters = pd.read_csv(HERE / "STRING_clusters.csv")
cluster_of = dict(zip(string_clusters["gene"], string_clusters["cluster"]))
sedges = {tuple(sorted((a, b)))
          for a, b in zip(string_edges["preferredName_A"], string_edges["preferredName_B"])}

string_rows = []
for r in modules.itertuples():
    gs = [g for g in r.all_genes.split(";") if g]
    pairs = list(itertools.combinations(sorted(gs), 2))
    n_edge = sum(1 for p in pairs if p in sedges)
    cl = Counter(cluster_of[g] for g in gs if g in cluster_of)
    top_cl, top_n = cl.most_common(1)[0] if cl else (None, 0)
    string_rows.append({
        "module": r.module,
        "module_label": r.module_label,
        "n_genes": len(gs),
        "string_edges_within_module": n_edge,
        "possible_pairs": len(pairs),
        "largest_shared_STRING_cluster": top_cl,
        "genes_in_that_cluster": top_n,
    })
string_check = pd.DataFrame(string_rows)
modules = modules.merge(
    string_check.drop(columns=["module_label", "n_genes"]), on="module", how="left")

# ------------------------------------------------------------ WT vs RING
ring_ids = set(ring_df["term_id"]) if len(ring_df) else set()
cmp_rows = []
for r in wt_df.itertuples():
    in_ring = r.term_id in ring_ids
    rr = ring_df[ring_df["term_id"] == r.term_id] if in_ring else None
    cmp_rows.append({
        "term_id": r.term_id,
        "term_name": r.term_name,
        "source": r.source,
        "term_size": r.term_size,
        "specific": r.specific,
        "module": r.module,
        "module_label": r.module_label,
        "survives_without_histones": r.survives_without_histones,
        "WT_p": r.p_value,
        "WT_genes_n": r.n_genes,
        "WT_genes": r.genes,
        "in_RING": in_ring,
        "RING_p": float(rr["p_value"].iloc[0]) if in_ring else np.nan,
        "RING_genes": rr["genes"].iloc[0] if in_ring else "",
        "status": "mutual WT+RING" if in_ring else "WT-exclusive",
    })
cmp_df = pd.DataFrame(cmp_rows).sort_values(["status", "WT_p"])

# ------------------------------------------------------------------ outputs
wt_df.to_csv(HERE / "WT_mutual_terms.csv", index=False)
matrix.to_csv(HERE / "WT_mutual_gene_matrix.csv")
genes_df.to_csv(HERE / "WT_mutual_genes.csv", index=False)
modules.to_csv(HERE / "WT_mutual_modules.csv", index=False)
cmp_df.to_csv(HERE / "WT_mutual_vs_RING.csv", index=False)

hub = genes_df[genes_df["n_specific_terms"] > 0]
summary = {
    "question": "Do the UBR3-WT best hits share mutual pathways?",
    "wt_best_hits": len(wt_genes),
    "wt_enriched_terms": int(len(wt_df)),
    "wt_specific_terms": int(len(specific)),
    "terms_supported_by_multiple_hits": int((wt_df["n_genes"] >= 2).sum()),
    "hits_in_at_least_one_specific_term": int(len(hub)),
    "hits_in_no_specific_term": int(len(genes_df) - len(hub)),
    "hits_sharing_a_specific_term_with_another_hit": int(
        (genes_df["n_mutual_partners"] > 0).sum()),
    "pathway_modules": int(len(modules)),
    "robust_pathway_modules": int(modules["robust"].sum()),
    "gene_modules": int(len(gene_comms)),
    "mutual_with_RING": int(cmp_df["in_RING"].sum()),
    "wt_exclusive": int((~cmp_df["in_RING"]).sum()),
    "histones_in_wt_hits": HISTONES,
    "terms_lost_without_histones": int(wt_df["histone_driven"].sum()),
    "parameters": {
        "tool": "g:Profiler g:GOSt REST API (no_evidences=False)",
        "organism": "hsapiens",
        "sources": SOURCES,
        "correction": "g_SCS",
        "user_threshold": 0.05,
        "background": "whole annotated genome",
        "specific_term_max_size": SPECIFIC_MAX_SIZE,
        "term_edge_overlap_coefficient": TERM_EDGE_OVERLAP,
        "gene_edge_min_shared_specific_terms": GENE_EDGE_MIN_SHARED,
        "community_detection": "greedy modularity (networkx)",
    },
    "modules": [
        {"module": int(r.module), "label": r.module_label, "n_terms": int(r.n_terms),
         "best_p": float(r.best_p), "robust_to_histone_removal": bool(r.robust),
         "string_edges_within_module": int(r.string_edges_within_module),
         "genes_in_one_STRING_cluster": int(r.genes_in_that_cluster),
         "core_genes": r.core_genes.split(";") if r.core_genes else [],
         "all_genes": r.all_genes.split(";") if r.all_genes else []}
        for r in modules.itertuples()
    ],
}
(HERE / "WT_mutual_summary.json").write_text(json.dumps(summary, indent=2))

# ------------------------------------------------------------------ console
print("\n=== pathway modules (mutual pathway groups) ===")
for r in modules.itertuples():
    tag = "robust" if r.robust else "HISTONE-DRIVEN"
    print(f"M{r.module}: {r.module_label}  ({r.n_terms} terms, best p={r.best_p:.2e}, {tag})")
    print(f"    all hits : {r.all_genes}")
    print(f"    core hits: {r.core_genes}")

print("\n=== hits appearing in the most specific pathways ===")
print(hub.head(20)[["gene", "wt_rank", "n_specific_terms", "n_pathway_modules",
                    "n_mutual_partners", "cellular_role"]].to_string(index=False))

print("\n=== STRING cross-check (independent of the annotation databases) ===")
print(string_check.to_string(index=False))

print("\n=== terms mutual between WT and RING ===")
print(cmp_df[cmp_df["in_RING"]][["term_name", "source", "WT_p", "RING_p"]].to_string(index=False))
print(json.dumps({k: v for k, v in summary.items()
                  if k not in ("parameters", "modules")}, indent=2))
