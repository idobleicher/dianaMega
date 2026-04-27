"""
MP Project — enrichment analysis vs ALL peptides in the library.

Background = all 23,408 unique transcripts from the GPS library
(union of Data S3A / S3B / S3C, dedup on Ensembl Transcript ID).
Every transcript starts with M (N-terminal screen), so this background
captures the global N-terminome distribution.

Because the hits all start with M-P but the library has all possible
second residues, position 2 will show a very strong "P enriched, all
others depleted" signal. Position 1 (always M) is invariant in both
sets and is included for completeness but has zero signal.

Outputs in mpProject/results/vs_full_library/:
  - hits_resolved.csv               same 155 hits as the MP-vs-MP run
  - library_full_peptides.csv       background = full library (23,408 rows)
  - enrichment.csv                  log2(hit/lib) — positions 1–24
  - counts_hits.csv                 raw hit counts per (AA, position)
  - counts_library.csv              raw library counts per (AA, position)
  - significance_signed_log10p.csv  Fisher's exact signed -log10 p
  - enrichment_property_groups.csv  enrichment per chemistry group
  - mp_enrichment_heatmap.png       primary heatmap (positions 1–24)
  - mp_frequency_heatmaps.png       hits vs full-library frequencies
  - mp_significance_heatmap.png     Fisher's exact significance
  - mp_enrichment_logo.png          logo with log2-enrichment letter heights
  - mp_property_group_heatmap.png   chemistry-group enrichment
  - mp_top_movers.png               per-position top movers (positions 2–7)
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from mp_plotting import run_full_panel_set

HERE = Path(__file__).parent
RESULTS = HERE / "results" / "vs_full_library"
RESULTS.mkdir(parents=True, exist_ok=True)

LIB_PATH = Path(
    r"c:\Users\User\Downloads"
    r"\Data S3. N-terminome GPS screen data in different genetic backgrounds (1).xlsx"
)
SHEETS = ["Data S3A", "Data S3B", "Data S3C"]


def load_full_library() -> pd.DataFrame:
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


def main():
    print("=" * 72)
    print("MP project — enrichment vs FULL library background")
    print("=" * 72)

    lib = load_full_library()
    print(f"GPS full library: {len(lib)} unique transcripts")
    lib.to_csv(RESULTS / "library_full_peptides.csv", index=False)

    # Reuse the resolved hits already written by the MP-vs-MP run if present;
    # otherwise resolve directly here.
    hits_csv = HERE / "results" / "vs_MP_library" / "hits_resolved.csv"
    if hits_csv.exists():
        hits = pd.read_csv(hits_csv)
        print(f"Reusing resolved hits from {hits_csv}: {len(hits)} entries")
    else:
        from mp_enrichment_heatmap import load_user_list, resolve_hits
        user = load_user_list()
        hits, _ = resolve_hits(user, lib)
        print(f"Resolved {len(hits)} hits in-place")

    if not len(hits):
        return
    hits.to_csv(RESULTS / "hits_resolved.csv", index=False)

    run_full_panel_set(
        hit_seqs=hits["AA_seq"].tolist(),
        lib_seqs=lib["AA_seq"].tolist(),
        out_dir=RESULTS,
        lib_label="full GPS library (all N-termini)",
        lib_short="full library",
        # Show positions 3–24 in every plot (consistent with the MP-vs-MP run).
        # Positions 1 and 2 are excluded: 1 is always M, and 2 trivially shows
        # P enriched simply because every hit is an MP peptide — neither
        # carries new biological signal.
        positions_main=list(range(3, 25)),
        positions_top=list(range(3, 8)),
    )
    print("Done.")


if __name__ == "__main__":
    main()
