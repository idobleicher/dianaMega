"""
Cross-check the UBR3-WT best-hit half-lives against two independent, published
protein-turnover datasets.

Datasets
  Yi et al. 2023        in-extract degradation, 0-15 h, WT/ATG7-KO/FIP200-KO.
                        Half-lives DERIVED here (see wt_besthits_halflife.py).
                        This is the dataset the current numbers come from.
  Li et al. 2021        cycloheximide chase + TMTpro in living cells, 0-8 h,
  (Mol Cell, Gygi lab)  four lines including HEK293T - the same background as the
                        co-IP. Table S2 reports the authors' own fitted half-lives.
                        IMPORTANT: S2 lists ONLY proteins that lost >15% at 8 h with
                        R2>0.8. Absence from it means "not degraded OR not quantified"
                        - Table S1 would be needed to tell those apart.
  Mathieson et al. 2018 dynamic SILAC in five non-dividing primary cell types
  (Nat Commun)          (B cells, NK cells, monocytes, hepatocytes, mouse neurons).
                        Downloaded automatically.

Two questions:
  1. Proteome-wide, does the in-extract assay track turnover in living cells at all?
     (Spearman on every protein shared between the datasets.)
  2. Per hit, which of the 50 WT best hits are short-lived in more than one system?

Outputs:
  WT_besthits_halflife_crosscheck.csv    per hit, all three datasets side by side
  WT_besthits_halflife_crosscheck.json   correlations + counts
"""
import json
import pathlib
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats

HERE = pathlib.Path(__file__).parent
LI_S2 = pathlib.Path(r"C:\Users\User\Downloads\NIHMS1747401-supplement-Table_S2.xlsx")
MATHIESON_URL = ("https://static-content.springer.com/esm/art%3A10.1038%2F"
                 "s41467-018-03106-1/MediaObjects/41467_2018_3106_MOESM5_ESM.xlsx")
MATHIESON_LOCAL = HERE / "external_Mathieson2018_halflives.xlsx"
LI_LINES = ["HEK293T", "U2OS", "HCT116", "RPE1"]

# ------------------------------------------------------------------ Yi (ours)
yi = pd.read_csv(HERE / "WT_besthits_halflife.csv")
yi_bg = pd.read_csv(HERE / "WT_besthits_halflife_background.csv")

# ------------------------------------------------------------------ Li 2021
li = {}
for cl in LI_LINES:
    d = pd.read_excel(LI_S2, sheet_name=cl)
    d["gene"] = d["Gene symbol"].astype(str)
    d["acc"] = d["Uniprot ID"].astype(str).str.split("-").str[0]
    d["Half-life"] = pd.to_numeric(d["Half-life"], errors="coerce")
    li[cl] = d
li_hek = li["HEK293T"]

# ------------------------------------------------------------------ Mathieson
if not MATHIESON_LOCAL.exists():
    print("downloading Mathieson et al. 2018 supplementary data 2 ...")
    urllib.request.urlretrieve(MATHIESON_URL, MATHIESON_LOCAL)
mat = pd.read_excel(MATHIESON_LOCAL, sheet_name="protein half lives high qual")
mat["gene"] = mat["gene_name"].astype(str)
# human primary cell types only; the mouse neuron columns are a different organism
mat_cols = [c for c in mat.columns if "half_life" in c and "Mouse" not in c]
for c in mat_cols:
    mat[c] = pd.to_numeric(mat[c], errors="coerce")
mat["median_half_life_h"] = mat[mat_cols].median(axis=1)
mat_hl = mat.dropna(subset=["median_half_life_h"]).set_index("gene")["median_half_life_h"]
mat_hl = mat_hl[~mat_hl.index.duplicated()]

# ------------------------------------------- Q1: proteome-wide concordance
# Yi gives a decay RATE k; higher k = faster. Li/Mathieson give half-lives;
# lower = faster. So the expected correlation between k and half-life is NEGATIVE.
# the background frame excludes the hits by construction - add them back, so the
# hits can be located inside the proteome-wide comparison
yi_all = pd.concat([
    yi_bg[["Gene Symbol", "decay_k_per_h", "frac_remaining_15h"]]
    .rename(columns={"Gene Symbol": "gene"}),
    yi.loc[yi["quantified"], ["degradomics_gene", "decay_k_per_h", "frac_remaining_15h"]]
    .rename(columns={"degradomics_gene": "gene"}),
], ignore_index=True)
yi_all = yi_all.replace([np.inf, -np.inf], np.nan).dropna(subset=["decay_k_per_h"])
yi_all = yi_all[~yi_all["gene"].duplicated()].set_index("gene")

corr = {}
for label, series in [("Li2021_HEK293T", li_hek.dropna(subset=["Half-life"])
                       .drop_duplicates("gene").set_index("gene")["Half-life"]),
                      ("Mathieson2018_primary_cells", mat_hl)]:
    shared = yi_all.index.intersection(series.index)
    if len(shared) < 10:
        continue
    a = yi_all.loc[shared, "decay_k_per_h"]
    b = series.loc[shared]
    rho, p = stats.spearmanr(a, b)
    rho2, p2 = stats.spearmanr(yi_all.loc[shared, "frac_remaining_15h"], b)
    corr[label] = {
        "n_shared_proteins": int(len(shared)),
        "spearman_rho_Yi_k_vs_halflife": float(rho),
        "p_value": float(p),
        "spearman_rho_Yi_15h_remaining_vs_halflife": float(rho2),
        "p_value_15h": float(p2),
        "expected_sign": "negative for k (faster decay = shorter half-life), "
                         "positive for 15h remaining",
    }

# ------------------------------------------- Q2: per-hit cross-check
rows = []
for r in yi.itertuples():
    g = r.gene
    rec = {
        "gene": g,
        "WT_rank": r.wt_rank,
        "pathway_modules": r.pathway_modules,
        "cellular_role": r.cellular_role,
        "Yi_extract_t_half_h": (round(r.halflife_h, 1)
                                if r.quantified and np.isfinite(r.halflife_h) else np.nan),
        "Yi_remaining_15h": round(r.frac_remaining_15h, 3) if r.quantified else np.nan,
        "Yi_status": ("not quantified" if not r.quantified
                      else "no decay" if not np.isfinite(r.halflife_h)
                      else "measured" if r.halflife_h <= 15 else "extrapolated"),
    }
    for cl in LI_LINES:
        d = li[cl]
        m = d[(d["gene"] == g) | (d["acc"] == str(r.gene))]
        rec[f"Li2021_{cl}_t_half_h"] = float(m["Half-life"].iloc[0]) if len(m) else np.nan
        if cl == "HEK293T":
            rec["Li2021_HEK293T_R2"] = float(m["R2"].iloc[0]) if len(m) else np.nan
            rec["Li2021_HEK293T_short_lived"] = (
                str(m["Short-lived"].iloc[0]) if len(m) and pd.notna(m["Short-lived"].iloc[0])
                else "")
    rec["Mathieson_primary_median_t_half_h"] = (
        round(float(mat_hl[g]), 1) if g in mat_hl.index else np.nan)
    n_li = sum(pd.notna(rec[f"Li2021_{cl}_t_half_h"]) for cl in LI_LINES)
    rec["n_Li_cell_lines_degraded"] = n_li
    rows.append(rec)

cross = pd.DataFrame(rows)

# concordance calls
cross["in_Li_HEK293T"] = cross["Li2021_HEK293T_t_half_h"].notna()
cross["short_in_both"] = (cross["Yi_status"].isin(["measured", "extrapolated"])
                          & cross["in_Li_HEK293T"])

summary = {
    "question": "Do the UBR3-WT best-hit half-lives reproduce in independent datasets?",
    "datasets": {
        "Yi2023": "in-extract degradation 0-15 h; half-lives derived by us",
        "Li2021": f"cycloheximide chase in living cells 0-8 h; {LI_LINES}; "
                  "Table S2 lists only proteins losing >15% at 8 h with R2>0.8",
        "Mathieson2018": "dynamic SILAC, 4 human primary cell types (mouse neurons excluded)",
    },
    "proteome_wide_concordance": corr,
    "per_hit": {
        "wt_best_hits": int(len(cross)),
        "in_Li_HEK293T_degraded_list": int(cross["in_Li_HEK293T"].sum()),
        "in_any_Li_cell_line": int((cross["n_Li_cell_lines_degraded"] > 0).sum()),
        "with_Mathieson_halflife": int(cross["Mathieson_primary_median_t_half_h"].notna().sum()),
        "short_in_Yi_and_Li_HEK293T": cross.loc[cross["short_in_both"], "gene"].tolist(),
    },
    "caveat": ("Absence from Li Table S2 means the protein either was not degraded "
               ">15% in 8 h or was not quantified; Table S1 is needed to separate "
               "those. Do not read absence as proof of stability."),
}

cross.to_csv(HERE / "WT_besthits_halflife_crosscheck.csv", index=False)
(HERE / "WT_besthits_halflife_crosscheck.json").write_text(json.dumps(summary, indent=2))

# proteome-wide merged frame, for the concordance scatter
li_hek_u = li_hek.dropna(subset=["Half-life"]).drop_duplicates("gene").set_index("gene")
shared = yi_all.index.intersection(li_hek_u.index)
scatter = pd.DataFrame({
    "gene": shared,
    "Yi_remaining_15h": yi_all.loc[shared, "frac_remaining_15h"].to_numpy(),
    "Yi_decay_k": yi_all.loc[shared, "decay_k_per_h"].to_numpy(),
    "Li_HEK293T_t_half_h": li_hek_u.loc[shared, "Half-life"].to_numpy(),
})
scatter["is_WT_best_hit"] = scatter["gene"].isin(
    set(yi.loc[yi["quantified"], "degradomics_gene"]))
scatter.to_csv(HERE / "WT_besthits_halflife_crosscheck_scatter.csv", index=False)

# add to the master workbook
book = HERE / "WT_besthits_master_list.xlsx"
if book.exists():
    out_cols = ["gene", "WT_rank", "pathway_modules", "Yi_extract_t_half_h",
                "Yi_remaining_15h", "Yi_status", "Li2021_HEK293T_t_half_h",
                "Li2021_HEK293T_R2", "Li2021_HEK293T_short_lived",
                "Li2021_U2OS_t_half_h", "Li2021_HCT116_t_half_h",
                "Li2021_RPE1_t_half_h", "n_Li_cell_lines_degraded",
                "Mathieson_primary_median_t_half_h", "cellular_role"]
    with pd.ExcelWriter(book, engine="openpyxl", mode="a",
                        if_sheet_exists="replace") as xl:
        cross[out_cols].to_excel(xl, sheet_name="Cross-dataset", index=False)
        ws = xl.sheets["Cross-dataset"]
        ws.freeze_panes = "B2"
        for col, w in {"A": 12, "C": 14, "D": 18, "E": 16, "F": 14, "G": 20,
                       "H": 16, "I": 22, "J": 18, "K": 20, "L": 18, "M": 22,
                       "N": 28, "O": 38}.items():
            ws.column_dimensions[col].width = w

# ------------------------------------------------------------------ console
print("=== proteome-wide concordance with the Yi extract assay ===")
for k, v in corr.items():
    print(f"{k}: n={v['n_shared_proteins']}  "
          f"rho(k vs t-half)={v['spearman_rho_Yi_k_vs_halflife']:+.3f} "
          f"(p={v['p_value']:.2g})   "
          f"rho(15h remaining vs t-half)={v['spearman_rho_Yi_15h_remaining_vs_halflife']:+.3f} "
          f"(p={v['p_value_15h']:.2g})")

print("\n=== WT best hits found in the Li 2021 HEK293T degraded list ===")
sub = cross[cross["in_Li_HEK293T"]].sort_values("Li2021_HEK293T_t_half_h")
print(sub[["gene", "WT_rank", "Li2021_HEK293T_t_half_h", "Li2021_HEK293T_R2",
           "Li2021_HEK293T_short_lived", "Yi_extract_t_half_h", "Yi_status",
           "pathway_modules"]].to_string(index=False))

print("\n=== hits in any Li cell line ===")
anyli = cross[cross["n_Li_cell_lines_degraded"] > 0]
print(anyli[["gene", "n_Li_cell_lines_degraded"] +
            [f"Li2021_{c}_t_half_h" for c in LI_LINES]].to_string(index=False))

print("\n=== Mathieson primary-cell half-lives, shortest first ===")
mm = cross.dropna(subset=["Mathieson_primary_median_t_half_h"]).nsmallest(
    12, "Mathieson_primary_median_t_half_h")
print(mm[["gene", "WT_rank", "Mathieson_primary_median_t_half_h",
          "Yi_extract_t_half_h", "Yi_status"]].to_string(index=False))
