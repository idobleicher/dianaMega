"""
MP Project — enrichment analysis vs all MP-starting peptides in the library.

Background = the 1,263 unique MP-starting peptides from
"Data S3. N-terminome GPS screen data in different genetic backgrounds (1).xlsx"
(union of sheets Data S3A / S3B / S3C, deduplicated on Ensembl Transcript ID).

Outputs in mpProject/results/vs_MP_library/:
  - hits_resolved.csv               155 user peptides matched to library
  - hits_unmatched.csv              0 entries (kept for completeness)
  - library_mp_peptides.csv         the 1,263-peptide MP background
  - enrichment.csv                  log2(hit/lib), AA × position
  - counts_hits.csv                 raw hit counts per (AA, position)
  - counts_library.csv              raw library counts per (AA, position)
  - significance_signed_log10p.csv  Fisher's exact signed -log10 p
  - enrichment_property_groups.csv  enrichment collapsed to chemistry groups
  - mp_enrichment_heatmap.png       primary heatmap (positions 3–24)
  - mp_frequency_heatmaps.png       hits vs library raw frequencies
  - mp_significance_heatmap.png     Fisher's exact significance
  - mp_enrichment_logo.png          logo with letter heights = log2 enrichment
  - mp_property_group_heatmap.png   chemistry-group enrichment
  - mp_top_movers.png               per-position top movers (positions 3–7)
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from mp_plotting import run_full_panel_set

HERE = Path(__file__).parent
RESULTS = HERE / "results" / "vs_MP_library"
RESULTS.mkdir(parents=True, exist_ok=True)

LIB_PATH = Path(
    r"c:\Users\User\Downloads"
    r"\Data S3. N-terminome GPS screen data in different genetic backgrounds (1).xlsx"
)
SHEETS = ["Data S3A", "Data S3B", "Data S3C"]


def load_library() -> pd.DataFrame:
    rows = []
    for s in SHEETS:
        df = pd.read_excel(LIB_PATH, sheet_name=s, header=[0, 1])
        enst = df[("Unnamed: 0_level_0", "Ensembl Transcript ID")].astype(str)
        gene = df[("Unnamed: 1_level_0", "Gene Symbol")].astype(str)
        aa_col = [c for c in df.columns if "Amino Acid Sequence" in str(c[1])][0]
        aa = df[aa_col].astype(str)
        rows.append(pd.DataFrame({
            "ENST_ID": enst, "Gene_Symbol": gene,
            "AA_seq": aa, "source_sheet": s,
        }))
    lib = pd.concat(rows, ignore_index=True)
    lib["ENST_short"] = lib["ENST_ID"].str.split(".").str[0]
    sources = lib.groupby("ENST_short")["source_sheet"].agg(
        lambda v: ",".join(sorted(set(v))))
    lib = lib.drop_duplicates("ENST_short").set_index("ENST_short")
    lib["source_sheets"] = sources
    return lib.drop(columns="source_sheet").reset_index()


def load_user_list() -> list[str]:
    return [ln.strip().replace(" 2", "").strip()
            for ln in (HERE / "peptide_list.txt").read_text().splitlines()
            if ln.strip()]


def resolve_hits(user, lib):
    by_enst = lib.set_index("ENST_short", drop=False)
    by_gene = lib.groupby("Gene_Symbol", sort=False)
    rows, unmatched = [], []
    for entry in user:
        gene = entry.split("|")[0]
        enst = entry.split("|")[1] if "|" in entry else None
        pick = None
        if enst and enst in by_enst.index:
            cand = by_enst.loc[[enst]]
            mp_cand = cand[cand["AA_seq"].str.startswith("MP")]
            if len(mp_cand):
                pick = mp_cand.iloc[0]
        elif gene in by_gene.groups:
            cand = by_gene.get_group(gene)
            mp_cand = cand[cand["AA_seq"].str.startswith("MP")]
            if len(mp_cand):
                pick = mp_cand.iloc[0]
        if pick is None:
            unmatched.append({"user_entry": entry,
                              "reason": "no MP peptide found"})
        else:
            rows.append({
                "user_entry": entry,
                "Gene_Symbol": pick["Gene_Symbol"],
                "ENST_ID": pick["ENST_ID"],
                "AA_seq": pick["AA_seq"],
                "source_sheets": pick["source_sheets"],
            })
    return pd.DataFrame(rows), pd.DataFrame(unmatched)


def main():
    print("=" * 72)
    print("MP project — enrichment vs MP-only library background")
    print("=" * 72)

    lib = load_library()
    print(f"GPS library (union of S3A/B/C): {len(lib)} unique transcripts")
    lib_mp = lib[lib["AA_seq"].str.startswith("MP")].copy()
    print(f"  of which MP-starting        : {len(lib_mp)}")
    lib_mp.to_csv(RESULTS / "library_mp_peptides.csv", index=False)

    user = load_user_list()
    hits, unmatched = resolve_hits(user, lib)
    hits.to_csv(RESULTS / "hits_resolved.csv", index=False)
    unmatched.to_csv(RESULTS / "hits_unmatched.csv", index=False)
    print(f"Resolved {len(hits)}/{len(user)} hits, {len(unmatched)} unmatched")

    if not len(hits):
        return

    run_full_panel_set(
        hit_seqs=hits["AA_seq"].tolist(),
        lib_seqs=lib_mp["AA_seq"].tolist(),
        out_dir=RESULTS,
        lib_label="all MP peptides in GPS library",
        lib_short="MP library",
        positions_main=list(range(3, 25)),   # skip invariant 1=M, 2=P
        positions_top=list(range(3, 8)),     # focused panel for top-movers
    )
    print("Done.")


if __name__ == "__main__":
    main()
