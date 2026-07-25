"""
Protein stability of the UBR3-WT co-IP best hits, using the Yi et al. (2023)
degradomics dataset (Supplementary Table 1).

IMPORTANT — that table contains no half-life column. It is a TMTpro-MS3 time
course of protein remaining in WT / ATG7-KO / FIP200(RB1CC1)-KO extracts, sampled
at 0 h (two replicates), 5 h, 10 h and 15 h, untreated (UT) or Torin-treated
(Tor), expressed as a ratio to the 0 h WT untreated baseline. Half-lives here are
therefore DERIVED, by log-linear regression of ln(ratio) on time:

    ln(ratio_t) = -k t   ->   t_half = ln2 / k

The assay runs for 15 h and the median protein retains ~82% at 15 h, so most
half-lives fall outside the measured window and are extrapolations. Two guards:
  * `frac_remaining_15h` (mean of the two 15 h ratios) is reported alongside and
    is the model-free primary measure — prefer it for ranking.
  * `beyond_window` flags any fitted half-life > 15 h; `no_decay` flags proteins
    whose fitted slope is flat or rising, for which no half-life is defined.

Hit matching is by UniProt accession, then by any accession in the MaxQuant
majority-protein group (this rescues histones, whose shared peptides are assigned
to different paralogs in the two studies), then by gene symbol.

Outputs (this folder):
  WT_besthits_halflife.csv          per hit: decay rate, derived t1/2, 15 h remaining
  WT_besthits_halflife_summary.json group comparisons + parameters
  WT_besthits_halflife.(png|pdf|svg)
"""
import json
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

HERE = pathlib.Path(__file__).parent
DEGRADOMICS = pathlib.Path(r"C:\Users\User\Downloads\Supplementary Table 1 (3).xlsx")

# WT untreated series -> (column, timepoint in hours)
UT_SERIES = [
    ("WT_UT1_ratio", 0.0), ("WT_UT2_ratio", 0.0),
    ("WT_UT5h1_ratio", 5.0), ("WT_UT5h2_ratio", 5.0),
    ("WT_UT10h1_ratio", 10.0), ("WT_UT10h2_ratio", 10.0),
    ("WT_UT15h1_ratio", 15.0), ("WT_UT15h2_ratio", 15.0),
]
TOR_15H = ["WT_Tor15h1_ratio", "WT_Tor15h2_ratio"]
ATG7_15H = ["ATG7_UT15h1_ratio", "ATG7_UT15h2_ratio"]
RB1CC1_15H = ["RB1CC1_UT15h1_ratio", "RB1CC1_UT15h2_ratio"]
WINDOW_H = 15.0

# ------------------------------------------------------------------ load
deg = pd.read_excel(DEGRADOMICS, sheet_name="Degradomics", header=3)
num_cols = [c for c, _ in UT_SERIES] + TOR_15H + ATG7_15H + RB1CC1_15H
for c in num_cols:
    deg[c] = pd.to_numeric(deg[c], errors="coerce")
deg["UniProt ID"] = deg["UniProt ID"].astype(str)
deg["Gene Symbol"] = deg["Gene Symbol"].astype(str)
deg["quantified"] = deg[[c for c, _ in UT_SERIES]].notna().all(axis=1)

# ------------------------------------------------------- derive decay / t1/2
times = np.array([t for _, t in UT_SERIES])
vals = deg[[c for c, _ in UT_SERIES]].to_numpy(dtype=float)


def fit_decay(y):
    """Log-linear fit of ratio vs time. Returns k (1/h), t_half (h), R^2."""
    if np.isnan(y).any() or (y <= 0).any():
        return np.nan, np.nan, np.nan
    ln = np.log(y)
    slope, intercept, r, _, _ = stats.linregress(times, ln)
    k = -slope
    t_half = np.log(2) / k if k > 0 else np.inf
    return k, t_half, r ** 2


fits = np.array([fit_decay(row) for row in vals])
deg["decay_k_per_h"] = fits[:, 0]
deg["halflife_h"] = fits[:, 1]
deg["fit_r2"] = fits[:, 2]
deg["frac_remaining_15h"] = deg[["WT_UT15h1_ratio", "WT_UT15h2_ratio"]].mean(axis=1)
deg["frac_remaining_15h_Torin"] = deg[TOR_15H].mean(axis=1)
deg["frac_remaining_15h_ATG7KO"] = deg[ATG7_15H].mean(axis=1)
deg["frac_remaining_15h_RB1CC1KO"] = deg[RB1CC1_15H].mean(axis=1)
deg["no_decay"] = deg["quantified"] & ~np.isfinite(deg["halflife_h"])
deg["beyond_window"] = deg["quantified"] & (deg["halflife_h"] > WINDOW_H)

# ------------------------------------------------------------------ hits
hits = pd.read_csv(HERE / "hits_table.csv")
wt = hits[hits["hit_UBR3wt"] == True].sort_values(
    "WT_293T_log2diff", ascending=False).reset_index(drop=True)
ring_genes = set(hits.loc[hits["hit_UBR3ring"] == True, "primary_gene"].dropna())

by_acc = {a: i for i, a in enumerate(deg["UniProt ID"])}
by_sym = {}
for i, s in enumerate(deg["Gene Symbol"]):
    by_sym.setdefault(s, i)


def match(acc, majority, sym):
    if acc in by_acc:
        return by_acc[acc], "accession"
    for a in str(majority).split(";"):
        a = a.strip()
        if a in by_acc:
            return by_acc[a], "protein group"
    if sym in by_sym:
        return by_sym[sym], "gene symbol"
    return None, "not found"


modules = pd.read_csv(HERE / "WT_mutual_modules.csv")
gene_module = {}
for r in modules.itertuples():
    for g in str(r.all_genes).split(";"):
        gene_module.setdefault(g, []).append(r.module)
mod_label = dict(zip(modules["module"], modules["module_label"]))

diff = pd.read_csv(HERE / "WTvsRING_differential.csv")
diff.columns = [c.strip().lstrip("\ufeff") for c in diff.columns]
ring_class = dict(zip(diff["gene"], diff["RING_vs_WT_class"]))

rows = []
for r in wt.itertuples():
    idx, how = match(str(r.primary_acc), r.Majority_protein_IDs, str(r.primary_gene))
    d = deg.iloc[idx] if idx is not None else None
    rows.append({
        "gene": r.primary_gene,
        "wt_rank": r.Index + 1,
        "WT_293T_log2diff": r.WT_293T_log2diff,
        "cellular_role": r.cellular_role,
        "Ndegron_class": r.Ndegron_class,
        "pathway_modules": ";".join(map(str, gene_module.get(r.primary_gene, []))),
        "RING_vs_WT_class": ring_class.get(r.primary_gene, ""),
        "matched_by": how,
        "degradomics_gene": d["Gene Symbol"] if d is not None else "",
        "quantified": bool(d["quantified"]) if d is not None else False,
        "WT_quant_peptides": d["WT Quant Peptides"] if d is not None else np.nan,
        "frac_remaining_15h": d["frac_remaining_15h"] if d is not None else np.nan,
        "decay_k_per_h": d["decay_k_per_h"] if d is not None else np.nan,
        "halflife_h": d["halflife_h"] if d is not None else np.nan,
        "fit_r2": d["fit_r2"] if d is not None else np.nan,
        "beyond_window": bool(d["beyond_window"]) if d is not None else False,
        "no_decay": bool(d["no_decay"]) if d is not None else False,
        "frac_remaining_15h_Torin": d["frac_remaining_15h_Torin"] if d is not None else np.nan,
        "frac_remaining_15h_ATG7KO": d["frac_remaining_15h_ATG7KO"] if d is not None else np.nan,
        "frac_remaining_15h_RB1CC1KO": d["frac_remaining_15h_RB1CC1KO"] if d is not None else np.nan,
        "autophagy_annotation": d["Autophagy"] if d is not None else "",
    })
res = pd.DataFrame(rows)
res_q = res[res["quantified"]].copy()

# ----------------------------------------------------- hits vs background
hit_names = set(res_q["degradomics_gene"])
bg = deg[deg["quantified"] & ~deg["Gene Symbol"].isin(hit_names)]


def compare(a, b, col):
    a, b = a[col].dropna(), b[col].dropna()
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {"n_hits": int(len(a)), "n_background": int(len(b)),
            "median_hits": float(np.median(a)), "median_background": float(np.median(b)),
            "mannwhitney_U": float(u), "p_value": float(p)}


cmp_15h = compare(res_q, bg, "frac_remaining_15h")
cmp_k = compare(res_q, bg, "decay_k_per_h")

# finite half-lives only, for a median that means something
fin = res_q[np.isfinite(res_q["halflife_h"])]
bg_fin = bg[np.isfinite(bg["halflife_h"])]

# ----------------------------------------------------- per-module breakdown
mod_rows = []
for m in sorted(modules["module"]):
    sub = res_q[res_q["pathway_modules"].str.split(";").apply(lambda L: str(m) in L)]
    if len(sub) == 0:
        continue
    c = compare(sub, bg, "frac_remaining_15h")
    mod_rows.append({
        "module": m, "label": mod_label[m], "n": len(sub),
        "median_frac_remaining_15h": float(sub["frac_remaining_15h"].median()),
        "median_halflife_h": float(np.median(sub["halflife_h"][np.isfinite(sub["halflife_h"])]))
        if np.isfinite(sub["halflife_h"]).any() else None,
        "p_vs_background": c["p_value"],
        "p_bonferroni": min(1.0, c["p_value"] * len(modules)),
        "genes": ";".join(sub["gene"]),
    })

# ------------------------------------------- N-degron and RING-class splits
def split_compare(col):
    out = {}
    for lvl, sub in res_q.groupby(res_q[col].fillna("(none)")):
        if len(sub) < 3:
            continue
        c = compare(sub, bg, "frac_remaining_15h")
        out[str(lvl)] = {"n": len(sub), "median_frac_remaining_15h": c["median_hits"],
                         "p_vs_background": c["p_value"],
                         "genes": ";".join(sub["gene"])}
    return out


ndegron_split = split_compare("Ndegron_class")
ringclass_split = split_compare("RING_vs_WT_class")

summary = {
    "source": str(DEGRADOMICS),
    "source_note": ("Yi et al. (2023) Supplementary Table 1, 'Degradomics' sheet. "
                    "No half-life column exists; half-lives are derived here by "
                    "log-linear fit to the 0/5/10/15 h WT untreated ratio series."),
    "wt_best_hits": int(len(res)),
    "matched_in_degradomics": int(res["degradomics_gene"].ne("").sum()),
    "quantified": int(len(res_q)),
    "not_quantified": res.loc[~res["quantified"], "gene"].tolist(),
    "matched_by": res["matched_by"].value_counts().to_dict(),
    "background_proteins": int(len(bg)),
    "fraction_remaining_15h": cmp_15h,
    "decay_rate_k": cmp_k,
    "derived_halflife": {
        "hits_median_h": float(fin["halflife_h"].median()) if len(fin) else None,
        "background_median_h": float(bg_fin["halflife_h"].median()) if len(bg_fin) else None,
        "hits_with_no_measurable_decay": int(res_q["no_decay"].sum()),
        "hits_beyond_15h_window": int(res_q["beyond_window"].sum()),
        "background_beyond_15h_window_pct": float(100 * bg["beyond_window"].mean()),
    },
    "modules": mod_rows,
    "by_Ndegron_class": ndegron_split,
    "by_RING_vs_WT_class": ringclass_split,
    "parameters": {
        "timepoints_h": [0, 0, 5, 5, 10, 10, 15, 15],
        "condition": "WT extract, untreated (UT), ratio to 0 h",
        "fit": "linregress of ln(ratio) on time; t_half = ln2 / k",
        "measurement_window_h": WINDOW_H,
        "test": "Mann-Whitney U, two-sided",
    },
}

res.to_csv(HERE / "WT_besthits_halflife.csv", index=False)
(HERE / "WT_besthits_halflife_summary.json").write_text(json.dumps(summary, indent=2))

# raw trajectories, so the figure can show measured points rather than fits only
ratio_cols = [c for c, _ in UT_SERIES]
curves = []
for r in res_q.itertuples():
    d = deg[deg["Gene Symbol"] == r.degradomics_gene].iloc[0]
    curves.append({"gene": r.gene, "degradomics_gene": r.degradomics_gene,
                   **{c: d[c] for c in ratio_cols}})
pd.DataFrame(curves).to_csv(HERE / "WT_besthits_halflife_curves.csv", index=False)

bg_out = bg[ratio_cols + ["frac_remaining_15h", "decay_k_per_h", "halflife_h"]].copy()
bg_out.insert(0, "Gene Symbol", bg["Gene Symbol"])
bg_out.to_csv(HERE / "WT_besthits_halflife_background.csv", index=False)

# ------------------------------------------------------------------ console
print(f"WT best hits: {len(res)} | matched: {summary['matched_in_degradomics']} | "
      f"quantified: {len(res_q)}")
print(f"not quantified: {', '.join(summary['not_quantified'])}")
print(f"\nfraction remaining at 15 h  hits {cmp_15h['median_hits']:.3f} vs "
      f"background {cmp_15h['median_background']:.3f}  p={cmp_15h['p_value']:.3g}")
print(f"decay rate k (1/h)          hits {cmp_k['median_hits']:.4f} vs "
      f"background {cmp_k['median_background']:.4f}  p={cmp_k['p_value']:.3g}")
print(f"derived median half-life    hits {summary['derived_halflife']['hits_median_h']:.1f} h vs "
      f"background {summary['derived_halflife']['background_median_h']:.1f} h")
print(f"  ({summary['derived_halflife']['hits_beyond_15h_window']}/{len(res_q)} hits "
      f"extrapolate beyond the 15 h window)")

print("\n=== most unstable WT best hits ===")
print(res_q.nsmallest(12, "frac_remaining_15h")[
    ["gene", "wt_rank", "frac_remaining_15h", "halflife_h", "fit_r2",
     "pathway_modules", "cellular_role"]].to_string(index=False))
print("\n=== most stable WT best hits ===")
print(res_q.nlargest(8, "frac_remaining_15h")[
    ["gene", "wt_rank", "frac_remaining_15h", "halflife_h", "cellular_role"]].to_string(index=False))
print("\n=== per module (vs proteome background) ===")
for m in mod_rows:
    hl = f"{m['median_halflife_h']:.1f} h" if m["median_halflife_h"] else "n/a"
    print(f"M{m['module']} {m['label'][:34]:36} n={m['n']:2}  "
          f"15h rem {m['median_frac_remaining_15h']:.3f}  t1/2 {hl:>8}  "
          f"p={m['p_vs_background']:.4f}  p_bonf={m['p_bonferroni']:.4f}")

print("\n=== by N-degron class ===")
for k, v in ndegron_split.items():
    print(f"  {k[:34]:36} n={v['n']:2}  15h rem {v['median_frac_remaining_15h']:.3f}  "
          f"p={v['p_vs_background']:.3f}")
print("\n=== by WT/RING class ===")
for k, v in ringclass_split.items():
    print(f"  {k[:44]:46} n={v['n']:2}  15h rem {v['median_frac_remaining_15h']:.3f}  "
          f"p={v['p_vs_background']:.3f}")
