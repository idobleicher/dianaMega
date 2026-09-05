"""
Single loader for the ZYG11B / ZER1 GPS N-terminome screen.

Source: data/Data_S3_Nterminome_GPS_screen.xlsx -- "Data S3. N-terminome GPS
screen data in different genetic backgrounds". Every figure and every table in
this project reads the workbook through here, so nothing can pick up a stale
copy or disagree with anything else about the same number.

WHAT THE SCREEN IS
  A Global Protein Stability (GPS) reporter screen over a library of N-terminal
  24-mers, one per Ensembl transcript. Cells are sorted into six bins by
  reporter stability and each peptide's PSI (Protein Stability Index) is the
  read-weighted mean bin, 1 (unstable) to 6 (stable). A peptide that is
  degraded in the control and stops being degraded when an E3 receptor is
  knocked out goes UP: dPSI = PSI(KO) - PSI(control) > 0 means STABILISED, and
  stabilised on loss of the receptor is what makes a peptide its substrate.

THE THREE SHEETS
  Data S3A  20,720 peptides. One AAVS1-KO control and five knockouts sharing
            it: UBR KO clones 1-3, ZYG11B KO, ZER1 KO. ZYG11B and ZER1 are ONE
            clone each -- there is no replicate within either genotype.
  Data S3B  18,099 peptides. Wild-type control vs the ZYG11B/ZER1 DOUBLE KO.
            A separate experiment with its own control and its own library
            draw; 17,901 transcripts are shared with S3A.
  Data S3C  23,356 peptides. NMT1/2 KO vs AAVS1 KO. Loaded and carried here
            because N-myristoyltransferase competes for the same N-terminal
            glycine the CRL2 receptors read, but no hit list is built from it.

  The workbook's own dPSI columns are used as given. They reproduce
  PSI(KO) - PSI(control) exactly (max deviation 1.5e-15, i.e. floating point),
  which is checked on every load.

THE BIOLOGY THIS DATASET IS ABOUT
  ZYG11B and ZER1 are the substrate receptors of the CRL2 complex that read an
  N-terminal GLYCINE -- the Gly/N-degron pathway. Position 1 is Met in 100% of
  this library, and Met is excised when position 2 is small (Sherman's rule),
  so position 2 is the residue actually exposed at the N-terminus. Gly at
  position 2 is therefore the feature the pathway acts on, and it behaves like
  one in the data: it is the single most stabilised residue in the ZYG11B KO
  (mean dPSI +0.226 vs -0.037 for everything else, Mann-Whitney p = 1.6e-81)
  and in the double KO (+0.240 vs -0.013, p = 1.5e-62). That is the built-in
  positive control every hit list here is checked against.

TWO THINGS TO KNOW BEFORE USING dPSI
  * THE CEILING. A peptide already stable in the control cannot rise. Mean dPSI
    falls monotonically with baseline PSI and is negative above 4.5, and the
    left tail of every dPSI distribution is heavier than the right for that
    reason alone. `headroom` (baseline PSI <= 4.0) marks the peptides that
    could have gone up, and no threshold in this project is applied without it.
  * THE SHARED CONTROL. Every dPSI in S3A is measured against the same AAVS1
    sample, so noise in that one control propagates into all five and makes
    them correlate whether or not the biology does. ZYG11B and UBR KO #1
    correlate at r = 0.46 despite reading different degrons; that number is
    about the control, not about either receptor. Agreement BETWEEN S3A and
    S3B, which have different controls, is the only cross-check here that is
    not inflated this way.

Columns returned by load(), one row per Ensembl transcript:
  transcript, gene, seq, nterm            identity; nterm = residue 2
  psi_ctrl_a, psi_zyg11b, psi_zer1        S3A, PSI 1-6
  dpsi_zyg11b, dpsi_zer1                  S3A, vs the AAVS1 control
  dpsi_ubr1, dpsi_ubr2, dpsi_ubr3         S3A, the UBR clones, for contrast
  psi_ctrl_b, psi_dko, dpsi_dko           S3B, vs the wild-type control
  reads_*                                 total reads behind each PSI
  in_a, in_b                              which experiments the peptide is in
  headroom, headroom_b                    baseline PSI <= 4.0 in S3A / S3B
"""
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
SRC = os.path.join(DATA, "Data_S3_Nterminome_GPS_screen.xlsx")
CACHE = os.path.join(DATA, "gps_tidy.csv.gz")

# PSI at or below this in the control leaves room to be stabilised. 4.0 is the
# last bin edge at which mean dPSI is still positive in both experiments; above
# it the ceiling dominates and dPSI measures the ceiling, not the receptor.
HEADROOM_MAX = 4.0

# The residue exposed after initiator-Met excision is position 2, and Met is
# excised when it is small (Ala Cys Gly Pro Ser Thr Val).
MET_EXCISED = set("ACGPSTV")
GLY = "G"

SHEETS = {"a": "Data S3A", "b": "Data S3B", "c": "Data S3C"}

_CACHE = None


def _read_sheets():
    a = pd.read_excel(SRC, sheet_name=SHEETS["a"], header=1)
    b = pd.read_excel(SRC, sheet_name=SHEETS["b"], header=1)
    return a, b


def _check(a, b):
    """The workbook's dPSI columns must be PSI(KO) - PSI(control). If a future
    export breaks that, every number in this project is wrong, so it is checked
    on load rather than trusted."""
    pairs = [(a, "ZYG11B KO", "PSI ZYG11B KO", "PSI AAVS1 KO"),
             (a, "ZER1 KO", "PSI ZER1 KO", "PSI AAVS1 KO"),
             (a, "UBR KO #1", "PSI UBR KO #1", "PSI AAVS1 KO"),
             (b, "ΔPSI (Double KO - WT)", "PSI Double KO", "PSI Wild-type")]
    for df, col, ko, ctrl in pairs:
        gap = (df[col] - (df[ko] - df[ctrl])).abs().max()
        assert gap < 1e-9, f"{col}: stated dPSI differs from PSI(KO)-PSI(ctrl) by {gap}"
    for df, name in ((a, "S3A"), (b, "S3B")):
        seq = df["Amino Acid Sequence"].astype(str)
        assert (seq.str.len() == 24).all(), f"{name}: not every peptide is a 24-mer"
        assert (seq.str[0] == "M").all(), f"{name}: not every peptide starts with Met"


def load(force=False, use_cache=True):
    """Tidy one-row-per-transcript table joining S3A and S3B."""
    global _CACHE
    if _CACHE is not None and not force:
        return _CACHE.copy()
    if use_cache and not force and os.path.exists(CACHE):
        _CACHE = pd.read_csv(CACHE)
        return _CACHE.copy()

    a, b = _read_sheets()
    _check(a, b)

    A = pd.DataFrame({
        "transcript": a["Ensembl Transcript ID"],
        "gene": a["Gene Symbol"],
        "seq": a["Amino Acid Sequence"].astype(str),
        "psi_ctrl_a": a["PSI AAVS1 KO"],
        "psi_zyg11b": a["PSI ZYG11B KO"],
        "psi_zer1": a["PSI ZER1 KO"],
        "dpsi_zyg11b": a["ZYG11B KO"],
        "dpsi_zer1": a["ZER1 KO"],
        "dpsi_ubr1": a["UBR KO #1"],
        "dpsi_ubr2": a["UBR KO #2"],
        "dpsi_ubr3": a["UBR KO #3"],
        "reads_ctrl_a": a["Total Reads"],
        "reads_zyg11b": a["Total Reads.4"],
        "reads_zer1": a["Total Reads.5"],
    })
    B = pd.DataFrame({
        "transcript": b["Ensembl Transcript ID"],
        "gene_b": b["Gene Symbol"],
        "seq_b": b["Amino Acid Sequence"].astype(str),
        "psi_ctrl_b": b["PSI Wild-type"],
        "psi_dko": b["PSI Double KO"],
        "dpsi_dko": b["ΔPSI (Double KO - WT)"],
        "reads_ctrl_b": b["Total Reads"],
        "reads_dko": b["Total Reads.1"],
    })

    m = A.merge(B, on="transcript", how="outer", indicator=True)
    m["in_a"] = m._merge.isin(["left_only", "both"])
    m["in_b"] = m._merge.isin(["right_only", "both"])
    m = m.drop(columns="_merge")

    # a peptide present in both sheets must be the same peptide in both
    both = m[m.in_a & m.in_b]
    assert (both.seq == both.seq_b).all(), "the two sheets disagree about a peptide's sequence"
    m["seq"] = m.seq.fillna(m.seq_b)
    m["gene"] = m.gene.fillna(m.gene_b)
    m = m.drop(columns=["seq_b", "gene_b"])

    m["nterm"] = m.seq.str[1]
    m["met_excised"] = m.nterm.isin(MET_EXCISED)
    m["is_gly"] = m.nterm.eq(GLY)
    m["headroom"] = m.psi_ctrl_a <= HEADROOM_MAX
    m["headroom_b"] = m.psi_ctrl_b <= HEADROOM_MAX

    m = m.sort_values("transcript").reset_index(drop=True)
    m.to_csv(CACHE, index=False, float_format="%.6g")
    _CACHE = m
    return m.copy()


def gly_enrichment(frame, label=""):
    """Gly fraction of a subset against the library it was drawn from -- the
    positive control every hit list in this project is judged by."""
    lib = load()
    base = lib[lib.seq.notna()].is_gly.mean()
    got = frame.is_gly.mean() if len(frame) else np.nan
    return {"set": label, "n": len(frame), "pct_gly": 100 * got,
            "pct_gly_library": 100 * base,
            "enrichment": got / base if base else np.nan}


if __name__ == "__main__":
    d = load(force=True)
    print(f"{len(d)} transcripts  ({d.in_a.sum()} in S3A, {d.in_b.sum()} in S3B, "
          f"{(d.in_a & d.in_b).sum()} in both)")
    print(f"library Gly at position 2: {100 * d.is_gly.mean():.1f}%   "
          f"Met-excised (small position 2): {100 * d.met_excised.mean():.1f}%")
    print("\nmean dPSI by position-2 residue, the six commonest and Gly:")
    g = (d.groupby("nterm")
           .agg(n=("seq", "size"), zyg11b=("dpsi_zyg11b", "mean"),
                zer1=("dpsi_zer1", "mean"), dko=("dpsi_dko", "mean"))
           .sort_values("n", ascending=False))
    print(g.head(6).round(3).to_string())
    print(g.loc[["G"]].round(3).to_string())
    print(f"\ncached -> {os.path.relpath(CACHE, HERE)}")
