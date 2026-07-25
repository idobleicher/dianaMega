"""
The UBR3-WT best hits that are short-lived in BOTH the Yi 2023 extract assay and at
least one independent published dataset.

"Short-lived" needs a definition in each dataset, because the assays differ:

  Yi 2023 (in extract, derived half-lives)
      Passes both tests from wt_besthits_halflife_significant.py:
        - decay is real          : BH q < 0.05 on the log-linear slope
        - fast for the proteome  : decay rate k above the 80th percentile of the
                                   8097-protein background
  Li 2021 (cycloheximide chase in living cells)
      Present in Table S2 for at least one cell line. Membership of that table IS
      the short-lived criterion: the authors included only proteins losing >15% at
      8 h with R2 > 0.8.
  Mathieson 2018 (dynamic SILAC, human primary cells)
      Median half-life below the 20th percentile of that dataset's own distribution.

A hit is "validated" when Yi says short-lived AND at least one independent dataset
agrees. Reported alongside: hits that the independent data calls short-lived but Yi
missed - those are false negatives of the extract assay, not disagreements to ignore.

Output: WT_besthits_validated_shortlived.csv + sheet in the master workbook.
"""
import pathlib

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).parent
LI_LINES = ["HEK293T", "U2OS", "HCT116", "RPE1"]
MATHIESON_PCTL = 20

cross = pd.read_csv(HERE / "WT_besthits_halflife_crosscheck.csv")
slopes = pd.read_csv(HERE / "WT_besthits_halflife_slopetests.csv")
sig = pd.read_csv(HERE / "WT_besthits_halflife_significant.csv")
mat = pd.read_excel(HERE / "external_Mathieson2018_halflives.xlsx",
                    sheet_name="protein half lives high qual")

# Mathieson's own distribution, to define "short" on its scale
mat_cols = [c for c in mat.columns if "half_life" in c and "Mouse" not in c]
for c in mat_cols:
    mat[c] = pd.to_numeric(mat[c], errors="coerce")
mat_med = mat[mat_cols].median(axis=1).dropna()
mat_cut = float(np.percentile(mat_med, MATHIESON_PCTL))

yi_sig = set(sig["Gene"])
df = cross.copy()

li_cols = [f"Li2021_{c}_t_half_h" for c in LI_LINES]
df["Li_n_lines"] = df[li_cols].notna().sum(axis=1)
df["Li_min_t_half_h"] = df[li_cols].min(axis=1)
df["short_Yi"] = df["gene"].isin(yi_sig)
df["short_Li"] = df["Li_n_lines"] > 0
df["short_Mathieson"] = df["Mathieson_primary_median_t_half_h"] < mat_cut
df["short_independent"] = df["short_Li"] | df["short_Mathieson"]

validated = df[df["short_Yi"] & df["short_independent"]].copy()
validated["evidence_datasets"] = validated.apply(
    lambda r: 1 + int(r["Li_n_lines"]) + int(bool(r["short_Mathieson"])), axis=1)
validated = validated.sort_values("Li_min_t_half_h").reset_index(drop=True)

missed = df[~df["short_Yi"] & df["short_independent"]].sort_values("Li_min_t_half_h")
yi_only = df[df["short_Yi"] & ~df["short_independent"]]

out_cols = ["gene", "WT_rank", "pathway_modules", "cellular_role",
            "Yi_extract_t_half_h", "Yi_remaining_15h",
            "Li2021_HEK293T_t_half_h", "Li2021_U2OS_t_half_h",
            "Li2021_HCT116_t_half_h", "Li2021_RPE1_t_half_h",
            "Li_n_lines", "Li_min_t_half_h",
            "Mathieson_primary_median_t_half_h", "evidence_datasets"]
val_out = validated[out_cols].copy()
val_out.insert(0, "#", range(1, len(val_out) + 1))
val_out.to_csv(HERE / "WT_besthits_validated_shortlived.csv", index=False)

book = HERE / "WT_besthits_master_list.xlsx"
if book.exists():
    with pd.ExcelWriter(book, engine="openpyxl", mode="a",
                        if_sheet_exists="replace") as xl:
        val_out.to_excel(xl, sheet_name="Validated short-lived", index=False)
        ws = xl.sheets["Validated short-lived"]
        ws.freeze_panes = "C2"
        for col, w in {"B": 12, "D": 16, "E": 38, "F": 18, "G": 17, "H": 20,
                       "I": 18, "J": 19, "K": 17, "L": 11, "M": 17, "N": 28,
                       "O": 16}.items():
            ws.column_dimensions[col].width = w

print(f"Mathieson 'short' cutoff ({MATHIESON_PCTL}th pctile of its own "
      f"distribution): {mat_cut:.1f} h\n")
print(f"=== VALIDATED SHORT-LIVED  (short in Yi AND in independent data): "
      f"{len(validated)} ===")
print(val_out.drop(columns=["cellular_role"]).to_string(index=False))

print(f"\n=== short in independent data but NOT flagged by Yi ({len(missed)}) ===")
print(missed[["gene", "WT_rank", "Yi_extract_t_half_h", "Yi_remaining_15h",
              "Li_n_lines", "Li_min_t_half_h",
              "Mathieson_primary_median_t_half_h", "pathway_modules"]]
      .to_string(index=False))

print(f"\n=== short in Yi only, no independent measurement available ({len(yi_only)}) ===")
print(yi_only[["gene", "WT_rank", "Yi_extract_t_half_h", "Yi_status",
               "pathway_modules"]].to_string(index=False))
