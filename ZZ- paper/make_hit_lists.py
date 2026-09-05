"""
The best-stabilised peptides of the ZYG11B / ZER1 GPS screen.

The equivalent of the UBR3 package's `01_UBR3_substrates` sheet, built for this
dataset. That list was handed to us; this one has to be defined, because Data S3
carries a dPSI for every peptide and no author-defined hit list. So the
definition is written down here, with the arithmetic that justifies it.

WHAT COUNTS AS A HIT
  Three conditions, all applied together:

  1. HEADROOM. Baseline PSI <= 4.0 in the experiment being read. A peptide
     already stable in the control cannot be stabilised further, and above 4.0
     the mean dPSI of the whole library is negative -- so an unfiltered
     threshold reports the ceiling rather than the receptor. This one filter is
     what makes the false-discovery estimates below sane: without it the
     negative tail is heavier than the positive one at every threshold and the
     mirror-null returns an FDR over 100%.

  2. dPSI >= 1.0. A full stability bin. Held constant across the lists so the
     genotypes are compared on one ruler, not on a per-genotype quantile.

  3. AGREEMENT ACROSS THE TWO EXPERIMENTS, for the headline list only. S3A
     (ZYG11B KO vs AAVS1) and S3B (double KO vs wild-type) have different
     controls and different library draws, so a peptide passing in both has
     cleared two independent measurements. Within S3A that check is not
     available: ZYG11B and ZER1 are one clone each, and every S3A dPSI is
     measured against the same AAVS1 sample, which correlates them whether or
     not the biology does.

HOW FALSE DISCOVERY IS ESTIMATED
  A mirror null. For threshold t, the peptides at dPSI <= -t are counted the
  same way as those at dPSI >= +t. Strong destabilisation on losing an E3
  receptor has no mechanism behind it here, so the negative tail is an estimate
  of how much of the positive tail is noise. It is an approximation -- the
  tails are not perfectly symmetric even after the headroom filter -- and it is
  quoted, not hidden, on every list.

THE POSITIVE CONTROL
  ZYG11B and ZER1 are the CRL2 receptors for an N-terminal glycine, so a
  correct hit list must be full of peptides whose position 2 is Gly (position 1
  is Met in 100% of the library and is excised when position 2 is small). Gly
  is 7.7% of the library. If a list is not far above that, it is noise, and the
  Gly fraction is printed for every list produced here as the check.

WHAT THE LISTS ARE
  hits_high_confidence   ZYG11B KO AND double KO both >= 1.0, both with
                         headroom. The headline list, and the one to use.
  hits_double_ko         double KO >= 1.0 with headroom. The cleanest single
                         experiment in the workbook.
  hits_zyg11b_ko         ZYG11B KO >= 1.0 with headroom.
  hits_zer1_ko           ZER1 KO >= 1.0 with headroom. Reported for
                         completeness and NOT recommended: see the README.
                         Its Gly fraction sits at the library rate, i.e. this
                         list carries no detectable Gly/N-degron signal.

Outputs:
  data/hits_*.csv                       one file per list, ranked
  ZZ_ZYG11B_ZER1_hit_lists.xlsx         the same four lists, readable
  data/hit_summary.csv                  the counts and FDRs behind the README
"""
import os

import numpy as np
import pandas as pd

import zz_gps_data as zz

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
XLSX = os.path.join(HERE, "ZZ_ZYG11B_ZER1_hit_lists.xlsx")

CUT = 1.0                      # dPSI, one full stability bin
d = zz.load()
LIB_GLY = 100 * d.is_gly.mean()


def mirror_fdr(frame, col, headroom_col, cut=CUT):
    """Peptides at dPSI <= -cut, counted over the same headroom-filtered set,
    as an estimate of how many of the +cut peptides are noise."""
    pool = d[d[headroom_col] & d[col].notna()]
    neg = int((pool[col] <= -cut).sum())
    return neg, (neg / len(frame) if len(frame) else np.nan)


def summarise(name, frame, col, headroom_col):
    neg, fdr = mirror_fdr(frame, col, headroom_col)
    return {"list": name, "n": len(frame),
            "pct_gly": 100 * frame.is_gly.mean() if len(frame) else np.nan,
            "gly_enrichment": (frame.is_gly.mean() / d.is_gly.mean()) if len(frame) else np.nan,
            "mirror_negative": neg, "est_fdr_pct": 100 * fdr,
            "median_dpsi": frame[col].median() if len(frame) else np.nan,
            "n_genes": frame.gene.nunique() if len(frame) else 0}


# ------------------------------------------------------------------ the lists
both = d[d.in_a & d.in_b]

high = both[both.headroom & both.headroom_b
            & (both.dpsi_zyg11b >= CUT) & (both.dpsi_dko >= CUT)].copy()
high["dpsi_mean"] = high[["dpsi_zyg11b", "dpsi_dko"]].mean(axis=1)
high["dpsi_min"] = high[["dpsi_zyg11b", "dpsi_dko"]].min(axis=1)
high = high.sort_values("dpsi_mean", ascending=False).reset_index(drop=True)

dko = d[d.headroom_b & (d.dpsi_dko >= CUT)].sort_values(
    "dpsi_dko", ascending=False).reset_index(drop=True)
zyg = d[d.headroom & (d.dpsi_zyg11b >= CUT)].sort_values(
    "dpsi_zyg11b", ascending=False).reset_index(drop=True)
zer = d[d.headroom & (d.dpsi_zer1 >= CUT)].sort_values(
    "dpsi_zer1", ascending=False).reset_index(drop=True)

LISTS = [
    ("hits_high_confidence", high, "dpsi_mean", "headroom_b",
     "ZYG11B KO and double KO both ≥ 1.0, both with headroom"),
    ("hits_double_ko", dko, "dpsi_dko", "headroom_b",
     "double KO − wild-type ≥ 1.0, baseline PSI ≤ 4.0"),
    ("hits_zyg11b_ko", zyg, "dpsi_zyg11b", "headroom",
     "ZYG11B KO − AAVS1 ≥ 1.0, baseline PSI ≤ 4.0"),
    ("hits_zer1_ko", zer, "dpsi_zer1", "headroom",
     "ZER1 KO − AAVS1 ≥ 1.0, baseline PSI ≤ 4.0 — no Gly signal, see README"),
]

# ------------------------------------------------------------------ readable columns
COLS = [
    ("rank", "Rank"),
    ("gene", "Gene symbol"),
    ("transcript", "Ensembl transcript"),
    ("nterm", "Residue 2 (N-terminus after Met excision)"),
    ("is_gly", "Gly N-degron"),
    ("met_excised", "Met excision predicted"),
    ("dpsi_mean", "Mean ΔPSI (ZYG11B, double KO)"),
    ("dpsi_min", "Lower of the two ΔPSI"),
    ("dpsi_zyg11b", "ΔPSI ZYG11B KO"),
    ("dpsi_zer1", "ΔPSI ZER1 KO"),
    ("dpsi_dko", "ΔPSI double KO"),
    ("psi_ctrl_a", "PSI AAVS1 control"),
    ("psi_zyg11b", "PSI ZYG11B KO"),
    ("psi_zer1", "PSI ZER1 KO"),
    ("psi_ctrl_b", "PSI wild-type control"),
    ("psi_dko", "PSI double KO"),
    ("dpsi_ubr1", "ΔPSI UBR KO #1 (contrast)"),
    ("seq", "N-terminal 24-mer"),
]


def readable(frame):
    f = frame.copy()
    f.insert(0, "rank", np.arange(1, len(f) + 1))
    for c, _ in COLS:
        if c not in f.columns:
            f[c] = np.nan
    f = f[[c for c, _ in COLS]].rename(columns=dict(COLS))
    f["Gly N-degron"] = f["Gly N-degron"].map({True: "yes", False: ""})
    f["Met excision predicted"] = f["Met excision predicted"].map({True: "yes", False: "no"})
    return f.round({k: 3 for k in f.columns if f[k].dtype.kind == "f"})


# ------------------------------------------------------------------ write
print(f"library Gly at position 2: {LIB_GLY:.1f}%   cut = ΔPSI ≥ {CUT}\n")
rows = []
for name, frame, col, hcol, note in LISTS:
    frame.to_csv(os.path.join(DATA, f"{name}.csv"), index=False, float_format="%.4g")
    s = summarise(name, frame, col if col != "dpsi_mean" else "dpsi_dko", hcol)
    s["definition"] = note
    rows.append(s)
    print(f"{name:<22} n={s['n']:>4}  {s['n_genes']:>4} genes   "
          f"Gly {s['pct_gly']:>5.1f}% ({s['gly_enrichment']:>4.1f}x library)   "
          f"mirror-negative {s['mirror_negative']:>3}  est. FDR {s['est_fdr_pct']:>5.1f}%")

summary = pd.DataFrame(rows)
summary.to_csv(os.path.join(DATA, "hit_summary.csv"), index=False, float_format="%.4g")

with pd.ExcelWriter(XLSX, engine="openpyxl") as xl:
    readme = pd.DataFrame({
        "ZYG11B / ZER1 GPS screen — best-stabilised peptides": [
            "Source: Data S3. N-terminome GPS screen data in different genetic backgrounds.",
            "PSI = read-weighted mean stability bin, 1 (unstable) to 6 (stable).",
            "ΔPSI = PSI(knockout) − PSI(control). Positive = STABILISED when the receptor is lost.",
            "",
            f"A hit is ΔPSI ≥ {CUT} (one full stability bin) with baseline PSI ≤ {zz.HEADROOM_MAX} "
            "in the same experiment.",
            "The headroom filter is not cosmetic: a peptide already stable cannot rise, and without "
            "it the false-discovery estimate exceeds 100%.",
            "",
            "ZYG11B and ZER1 are the CRL2 receptors for an N-terminal glycine. Position 1 is Met in "
            "100% of the library and is excised when position 2 is small, so position 2 is the real "
            "N-terminus. Gly there is 7.7% of the library — the Gly column of each sheet is the "
            "positive control.",
            "",
            "01_high_confidence is the list to use: stabilised in the ZYG11B KO AND in the double "
            "KO, two experiments with different controls.",
            "04_zer1_ko carries no Gly enrichment at all and should not be read as a substrate list.",
            "",
            "Every ΔPSI here rests on ONE clone per genotype. There is no within-genotype replicate "
            "in this workbook, which is why agreement between experiments is what the headline list "
            "is built on.",
        ]})
    readme.to_excel(xl, sheet_name="00_README", index=False)
    for i, (name, frame, *_rest) in enumerate(LISTS, start=1):
        sheet = f"{i:02d}_{name.replace('hits_', '')}"[:31]
        readable(frame).to_excel(xl, sheet_name=sheet, index=False)
    summary.to_excel(xl, sheet_name="05_summary", index=False)

    for ws in xl.book.worksheets:
        ws.freeze_panes = "A2"
        if ws.title != "00_README":
            ws.auto_filter.ref = ws.dimensions
        for col in ws.columns:
            width = max((len(str(c.value)) for c in col[:60] if c.value is not None), default=10)
            ws.column_dimensions[col[0].column_letter].width = min(max(width + 2, 10), 46)

print(f"\nwrote {os.path.basename(XLSX)} and {len(LISTS)} CSVs + hit_summary.csv")

# ------------------------------------------------------------------ readout
print("\n" + "=" * 78)
print("the top 20 of the high-confidence list")
print("=" * 78)
top = readable(high).head(20)[
    ["Rank", "Gene symbol", "Residue 2 (N-terminus after Met excision)",
     "Mean ΔPSI (ZYG11B, double KO)", "ΔPSI ZYG11B KO", "ΔPSI double KO",
     "ΔPSI ZER1 KO", "PSI AAVS1 control"]]
print(top.to_string(index=False))

print("\n" + "=" * 78)
print("does the list behave like a Gly/N-degron list?")
print("=" * 78)
for name, frame, *_ in LISTS:
    if not len(frame):
        continue
    print(f"  {name:<22} Gly {100 * frame.is_gly.mean():>5.1f}%   "
          f"non-Gly residues at position 2: "
          f"{', '.join(f'{k}({v})' for k, v in frame[~frame.is_gly].nterm.value_counts().head(6).items())}")

print("\nZER1 on its own, against the library it was drawn from:")
pool = d[d.headroom]
print(f"   library with headroom: Gly {100 * pool.is_gly.mean():.1f}%   "
      f"ZER1 hits: Gly {100 * zer.is_gly.mean():.1f}%   "
      f"ZYG11B hits: Gly {100 * zyg.is_gly.mean():.1f}%   "
      f"double-KO hits: Gly {100 * dko.is_gly.mean():.1f}%")
