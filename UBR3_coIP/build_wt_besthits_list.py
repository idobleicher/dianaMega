"""
One consolidated list of the 50 UBR3-WT co-IP best hits, joining the three analyses
in this folder: co-IP enrichment (hits_table.csv), the mutual pathway modules
(WT_mutual_*.csv) and the degradomics stability (WT_besthits_halflife.csv).

Output: WT_besthits_master_list.xlsx  (also .csv for the main sheet)
  Sheet "Best hits"       50 hits, ranked by WT/293T enrichment
  Sheet "Pathway modules" the 4 mutual modules + member genes
  Sheet "Mutual pathways" the 29 specific enriched terms + supporting hits
  Sheet "Notes"           what each column means and what it must not be read as
"""
import pathlib

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).parent

MODULE_NAME = {
    1: "M1 mRNA silencing & RNP granules",
    2: "M2 ATP-dependent chromatin remodeling",
    3: "M3 Nucleosome / histone core",
    4: "M4 BAF-LAP2b nuclear-envelope complex",
}

hits = pd.read_csv(HERE / "hits_table.csv")
wt = hits[hits["hit_UBR3wt"] == True].sort_values(
    "WT_293T_log2diff", ascending=False).reset_index(drop=True)

life = pd.read_csv(HERE / "WT_besthits_halflife.csv")
mods = pd.read_csv(HERE / "WT_mutual_modules.csv")
terms = pd.read_csv(HERE / "WT_mutual_terms.csv")
mgenes = pd.read_csv(HERE / "WT_mutual_genes.csv")
diff = pd.read_csv(HERE / "WTvsRING_differential.csv")
diff.columns = [c.strip().lstrip("﻿") for c in diff.columns]

gene_mods = {}
for r in mods.itertuples():
    for g in str(r.all_genes).split(";"):
        gene_mods.setdefault(g, []).append(int(r.module))

life_i = life.set_index("gene")
mg_i = mgenes.set_index("gene")
ring_class = dict(zip(diff["gene"], diff["RING_vs_WT_class"]))

rows = []
for r in wt.itertuples():
    g = r.primary_gene
    mlist = gene_mods.get(g, [])
    L = life_i.loc[g] if g in life_i.index else None
    M = mg_i.loc[g] if g in mg_i.index else None

    quantified = bool(L["quantified"]) if L is not None else False
    frac = L["frac_remaining_15h"] if quantified else np.nan
    t_half = L["halflife_h"] if quantified else np.nan
    if not quantified:
        stability = "not quantified"
    elif not np.isfinite(t_half):
        stability = "no decay in 15 h"
    elif t_half > 15:
        stability = "t1/2 > 15 h assay (extrapolated)"
    else:
        stability = "t1/2 within assay window"

    rows.append({
        "WT_rank": r.Index + 1,
        "Gene": g,
        "Protein": str(r.Protein_names).split(";")[0],
        "UniProt": r.primary_acc,
        "WT/293T log2 enrichment": round(r.WT_293T_log2diff, 2),
        "WT/293T q": r.WT_293T_qval,
        "RING comparison": ring_class.get(g, ""),
        "Cellular role": r.cellular_role,
        "Pathway module(s)": ";".join(f"M{m}" for m in mlist),
        "Module name(s)": " | ".join(MODULE_NAME[m] for m in mlist),
        "Specific terms (n)": int(M["n_specific_terms"]) if M is not None else 0,
        "Shares pathways with": M["mutual_partners"] if M is not None else "",
        "Remaining at 15 h": round(frac, 3) if quantified else np.nan,
        "Derived t1/2 (h)": (round(t_half, 1) if quantified and np.isfinite(t_half)
                             else np.nan),
        "Stability call": stability,
        "Decay fit R2": round(L["fit_r2"], 2) if quantified and np.isfinite(t_half) else np.nan,
        "N-degron class": r.Ndegron_class,
        "MW (kDa)": r.MolWeight_kDa,
        "Razor+unique peptides": r.razor_unique_peptides,
    })

master = pd.DataFrame(rows)

mod_sheet = mods[["module", "module_label", "n_genes", "n_terms",
                  "n_terms_surviving_histone_removal", "robust",
                  "string_edges_within_module", "possible_pairs",
                  "core_genes", "all_genes"]].copy()
mod_sheet.insert(1, "module_name", mod_sheet["module"].map(MODULE_NAME))
mod_sheet = mod_sheet.rename(columns={
    "module_label": "named_after_term", "robust": "survives_histone_removal"})

term_sheet = terms[terms["specific"]][
    ["module", "module_label", "source", "term_id", "term_name", "term_size",
     "n_genes", "genes", "p_value", "survives_without_histones"]].sort_values(
    ["module", "p_value"])

notes = pd.DataFrame({
    "Column / sheet": [
        "WT/293T log2 enrichment",
        "RING comparison",
        "Pathway module(s)",
        "Specific terms (n)",
        "Remaining at 15 h",
        "Derived t1/2 (h)",
        "Stability call",
        "N-degron class",
        "Sheet: Mutual pathways",
        "Overall stability result",
        "M2 stability trend",
    ],
    "What it means / caution": [
        "log2 fold enrichment of the co-IP over the 293T control; hits are q<0.05.",
        "From WTvsRING_differential.csv. 'WT-preferring' = binding needs an intact RING.",
        "Mutual pathway module from wt_mutual_pathways.py. Blank = this hit shares no "
        "specific enriched term with any other hit.",
        "Number of enriched terms of size <=500 containing this hit. 25 of the 50 hits "
        "score 0 - no pathway statement covers them.",
        "Fraction remaining at 15 h in Yi et al. (2023) extract degradomics. "
        "Model-free; RANK ON THIS, not on the half-life.",
        "DERIVED by log-linear fit, not present in the source table. 39 of 45 quantified "
        "hits never reach 50% inside the 15 h assay, so most values are extrapolations.",
        "Flags whether the half-life is inside the measured window.",
        "From hits_table.csv. All 45 quantified hits are 'Stabilizing', so this column "
        "cannot separate them.",
        "The 29 specific enriched terms. survives_without_histones=FALSE means the term "
        "disappears when the 7 histone hits are removed - histone-carried annotation, "
        "do not report as a UBR3-WT pathway.",
        "The hits as a group are NOT less stable than the proteome "
        "(0.829 vs 0.815 remaining at 15 h, Mann-Whitney p=0.80).",
        "M2 turns over ~2x faster (0.576, median derived t1/2 19 h) but p=0.069 raw / "
        "0.28 Bonferroni with n=9. A lead to test, not a result.",
    ],
})

# --- one-line-per-protein half-life sheet ----------------------------------
# A single number per hit. Half-lives longer than the 15 h assay are fitted
# extrapolations, so each one carries the basis it rests on.
def basis(row):
    if row["Stability call"] == "not quantified":
        return "no data - protein not quantified in the degradomics set"
    if row["Stability call"] == "no decay in 15 h":
        return "no decay detected in 15 h - half-life not determinable"
    if row["Derived t1/2 (h)"] <= 15:
        return "measured - drops below 50% inside the 15 h assay"
    return "extrapolated - still above 50% at 15 h, fitted beyond the assay"


halflife = master[["WT_rank", "Gene", "Derived t1/2 (h)", "Remaining at 15 h",
                   "Stability call", "Decay fit R2", "Pathway module(s)"]].copy()
halflife["Basis"] = master.apply(basis, axis=1)
halflife = halflife.rename(columns={"Derived t1/2 (h)": "Half-life (h)"})
halflife = halflife.drop(columns=["Stability call"]).sort_values(
    "Half-life (h)", na_position="last").reset_index(drop=True)
halflife.insert(0, "Order", range(1, len(halflife) + 1))

out = HERE / "WT_besthits_master_list.xlsx"
with pd.ExcelWriter(out, engine="openpyxl") as xl:
    halflife.to_excel(xl, sheet_name="Half-life", index=False)
    master.to_excel(xl, sheet_name="Best hits", index=False)
    mod_sheet.to_excel(xl, sheet_name="Pathway modules", index=False)
    term_sheet.to_excel(xl, sheet_name="Mutual pathways", index=False)
    notes.to_excel(xl, sheet_name="Notes", index=False)

    widths = {"Half-life": {"C": 12, "D": 14, "E": 18, "F": 13, "G": 18, "H": 60},
              "Best hits": {"B": 12, "C": 46, "D": 11, "H": 30, "I": 16, "J": 38,
                            "L": 34, "O": 30},
              "Pathway modules": {"B": 38, "C": 26, "J": 40, "K": 60},
              "Mutual pathways": {"C": 10, "E": 52, "H": 60},
              "Notes": {"A": 26, "B": 100}}
    for sheet, cols in widths.items():
        ws = xl.sheets[sheet]
        ws.freeze_panes = "A2"
        for col, w in cols.items():
            ws.column_dimensions[col].width = w

master.to_csv(HERE / "WT_besthits_master_list.csv", index=False)
halflife.to_csv(HERE / "WT_besthits_halflife_list.csv", index=False)

print(f"-> {out.name}  ({len(master)} hits, {len(mod_sheet)} modules, "
      f"{len(term_sheet)} specific terms)")
print("-> WT_besthits_master_list.csv, WT_besthits_halflife_list.csv")
print(f"\nin a module: {(master['Pathway module(s)'] != '').sum()} / {len(master)}")
print(f"stability quantified: {(master['Stability call'] != 'not quantified').sum()} / {len(master)}")
print("\n=== half-life, one row per hit ===")
print(halflife.to_string(index=False))
