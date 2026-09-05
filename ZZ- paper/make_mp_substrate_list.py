"""
Every Met-Pro peptide in the screen, ranked, with the stabilised ones tiered.

WHY TIERS AND NOT ONE CUT. The hit lists in make_hit_lists.py use a single
threshold, dPSI >= 1.0 with headroom, chosen for precision. Of the two peptides
confirmed in the lab, VWA5B1 clears it in ONE genotype (+1.21 in ZER1, +0.42
and +0.63 elsewhere) and SEPTIN4 in two but only just (+1.17, +1.02, and +0.26
in ZYG11B). Neither reaches the project's high-confidence list, which wants 1.0
in two experiments with different controls. With no within-genotype replicate
anywhere in this workbook a real substrate can easily land at 0.6 in a given
readout, so this list keeps the same measurements and shows the gradient
instead of collapsing it:

  Tier A   dPSI >= 1.0 in two or more of ZYG11B / ZER1 / double KO      n = 7
  Tier B   dPSI >= 1.0 in exactly one                                   n = 29
  Tier C   best dPSI 0.5-1.0, none reaching 1.0 -- the watch list      n = 144
  (below)  best dPSI < 0.5 with headroom, and everything with no
           headroom at all, kept in the all-peptides sheet only

Tier is a statement about how much evidence there is, not about how good a
substrate something is. SEPTIN4 lands in tier A and VWA5B1 in tier B, and
VWA5B1's only qualifying readout is the ZER1 knockout -- the one list this
project's README warns carries no glycine signal. It is a real substrate all
the same, which is exactly why tier C exists.

MET-PRO-GLY IS FLAGGED THROUGHOUT. Both confirmed peptides are Met-Pro-Gly, and
that group behaves like the canonical Met-Gly substrates: mean dPSI +0.313 and
a 12.5% hit rate in the double KO, against +0.064 and 2.4% for other Met-Pro
peptides (p = 0.0001). See the README. Met-Pro-Gly is 11.1% of Met-Pro peptides
overall and 34% of tier B.

THE UBR COLUMN IS A CAVEAT, NOT A BONUS. Data S3A also carries three UBR
knockout clones, and Met-Pro is UBR3's territory in the sibling project
(`../ubr3_PG_reviewer/`, the [P/G]-[E/D] motif). A peptide stabilised both when
UBR is lost and when ZYG11B/ZER1 are lost is not cleanly attributable to
either, so `also_UBR` marks the 6 tier-A/B peptides where all three UBR clones
reach dPSI >= 0.5. Note that Met-Pro-acidic peptides show NO UBR enrichment in
this dataset (mean UBR dPSI +0.085 against +0.099 for other Met-Pro peptides),
so the UBR3 motif does not sort these peptides here.

Met excision is inefficient at Met-Pro, so most of these peptides probably keep
their initiator methionine and present Met, not Pro. Whatever stabilises them
need not be an N-degron mechanism.

Outputs:
  ZZ_MetPro_substrates.xlsx     README + the three tiers + all 1,122 peptides
  data/mp_substrates.csv        tiers A-C, 180 peptides, ranked
  data/mp_all_peptides.csv      every Met-Pro peptide, ranked
"""
import os

import numpy as np
import pandas as pd

import zz_gps_data as zz

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
XLSX = os.path.join(HERE, "ZZ_MetPro_substrates.xlsx")

KO = ["dpsi_zyg11b", "dpsi_zer1", "dpsi_dko"]
UBR = ["dpsi_ubr1", "dpsi_ubr2", "dpsi_ubr3"]
STRONG, WATCH = 1.0, 0.5
LAB_CONFIRMED = {"ENST00000485375.6": "VWA5B1", "ENST00000412945.7": "SEPTIN4"}

d = zz.load()
mp = d[d.nterm == "P"].copy()
mp["p3"] = mp.seq.str[2]
mp["is_mpg"] = mp.p3.eq("G")
mp["headroom_any"] = mp.headroom.fillna(False) | mp.headroom_b.fillna(False)
mp["best_dpsi"] = mp[KO].max(axis=1)
mp["n_strong"] = mp[KO].ge(STRONG).sum(axis=1)
mp["n_measured"] = mp[KO].notna().sum(axis=1)
mp["ubr_mean"] = mp[UBR].mean(axis=1)
mp["also_ubr"] = mp[UBR].ge(WATCH).all(axis=1)
mp["lab_confirmed"] = mp.transcript.map(LAB_CONFIRMED).notna()


def tier(r):
    if not r.headroom_any:
        return "— no headroom"
    if r.n_strong >= 2:
        return "A"
    if r.n_strong == 1:
        return "B"
    if r.best_dpsi >= WATCH:
        return "C"
    return "— below 0.5"


mp["tier"] = mp.apply(tier, axis=1)
mp = mp.sort_values(["best_dpsi"], ascending=False).reset_index(drop=True)

subs = mp[mp.tier.isin(list("ABC"))].copy()
subs = subs.sort_values(
    ["tier", "best_dpsi"], ascending=[True, False]).reset_index(drop=True)

COLS = [
    ("tier", "Tier"), ("gene", "Gene symbol"), ("transcript", "Ensembl transcript"),
    ("nterm3", "N-terminal triplet"), ("p3", "Residue 3"),
    ("is_mpg", "Met-Pro-Gly"), ("lab_confirmed", "Confirmed in lab"),
    ("best_dpsi", "Best ΔPSI"), ("n_strong", "Genotypes at ΔPSI ≥ 1.0"),
    ("dpsi_zyg11b", "ΔPSI ZYG11B KO"), ("dpsi_zer1", "ΔPSI ZER1 KO"),
    ("dpsi_dko", "ΔPSI double KO"),
    ("psi_ctrl_a", "PSI AAVS1 control"), ("psi_ctrl_b", "PSI wild-type control"),
    ("headroom", "Headroom in S3A"), ("headroom_b", "Headroom in S3B"),
    ("ubr_mean", "Mean ΔPSI UBR clones"), ("also_ubr", "Also UBR (all 3 ≥ 0.5)"),
    ("seq", "N-terminal 24-mer"),
]


def readable(frame):
    f = frame.copy()
    f.insert(0, "rank", np.arange(1, len(f) + 1))
    f = f[["rank"] + [c for c, _ in COLS]].rename(
        columns={**dict(COLS), "rank": "Rank"})
    for c in ("Met-Pro-Gly", "Confirmed in lab", "Also UBR (all 3 ≥ 0.5)"):
        f[c] = f[c].map({True: "yes", False: ""})
    for c in ("Headroom in S3A", "Headroom in S3B"):
        f[c] = f[c].map({True: "yes", False: "no"})
    return f.round({k: 3 for k in f.columns if f[k].dtype.kind == "f"})


subs.to_csv(os.path.join(DATA, "mp_substrates.csv"), index=False, float_format="%.4g")
mp.to_csv(os.path.join(DATA, "mp_all_peptides.csv"), index=False, float_format="%.4g")

lib_mpg = mp.is_mpg.mean()
readme = pd.DataFrame({"Met-Pro peptides of the ZYG11B / ZER1 GPS screen": [
    f"{len(mp)} peptides start Met-Pro, across {mp.gene.nunique()} genes. "
    f"{int(mp.headroom_any.sum())} have headroom (baseline PSI ≤ 4.0) in at least one experiment "
    "and could therefore show stabilisation at all.",
    "",
    "ΔPSI = PSI(knockout) − PSI(control); positive means stabilised when the receptor is lost.",
    "",
    f"Tier A  ΔPSI ≥ 1.0 in two or more of ZYG11B / ZER1 / double KO — {int((subs.tier=='A').sum())} peptides",
    f"Tier B  ΔPSI ≥ 1.0 in exactly one — {int((subs.tier=='B').sum())} peptides",
    f"Tier C  best ΔPSI 0.5–1.0, the watch list — {int((subs.tier=='C').sum())} peptides",
    "",
    "Tier is how much evidence there is, not how good a substrate something is. There is no "
    "within-genotype replicate anywhere in this workbook, so one clone's ΔPSI of 0.6 is not "
    "evidence of absence. Of the two peptides confirmed in the lab, SEPTIN4 is tier A and "
    "VWA5B1 is tier B — and VWA5B1 qualifies only through the ZER1 knockout, the readout with "
    "the weakest aggregate signal in this screen.",
    "",
    f"Met-Pro-Gly is flagged throughout. It is {lib_mpg:.1%} of Met-Pro peptides overall and "
    f"{(subs[subs.tier=='B'].is_mpg.mean()):.0%} of tier B, and it behaves like the canonical "
    "Met-Gly substrates the receptors are known to read.",
    "",
    "'Also UBR' is a CAVEAT: all three UBR knockout clones stabilise that peptide too, so it is "
    "not cleanly attributable to ZYG11B or ZER1. Met-Pro is UBR3's territory in the sibling "
    "project.",
    "",
    "Met excision is inefficient at Met-Pro, so most of these peptides probably keep their "
    "initiator methionine and present Met rather than Pro at the N-terminus.",
]})

with pd.ExcelWriter(XLSX, engine="openpyxl") as xl:
    readme.to_excel(xl, sheet_name="00_README", index=False)
    for t, nm in (("A", "01_tierA_two_genotypes"), ("B", "02_tierB_one_genotype"),
                  ("C", "03_tierC_watchlist")):
        readable(subs[subs.tier == t]).to_excel(xl, sheet_name=nm, index=False)
    readable(mp).to_excel(xl, sheet_name="04_all_MetPro", index=False)
    for ws in xl.book.worksheets:
        ws.freeze_panes = "A2"
        if ws.title != "00_README":
            ws.auto_filter.ref = ws.dimensions
        for col in ws.columns:
            w = max((len(str(c.value)) for c in col[:80] if c.value is not None), default=10)
            ws.column_dimensions[col[0].column_letter].width = min(max(w + 2, 10), 44)

print(f"{len(mp)} Met-Pro peptides, {int(mp.headroom_any.sum())} with headroom\n")
for t in "ABC":
    k = subs[subs.tier == t]
    print(f"  tier {t}: {len(k):>4} peptides   Met-Pro-Gly {int(k.is_mpg.sum()):>3} "
          f"({k.is_mpg.mean():>4.0%}, library {lib_mpg:.0%})   also UBR {int(k.also_ubr.sum())}")
print(f"\nwrote {os.path.basename(XLSX)}, data/mp_substrates.csv ({len(subs)}), "
      f"data/mp_all_peptides.csv ({len(mp)})")

print("\n" + "=" * 78)
print("tiers A and B in full")
print("=" * 78)
show = ["tier", "gene", "nterm3", "is_mpg", "lab_confirmed", "dpsi_zyg11b",
        "dpsi_zer1", "dpsi_dko", "also_ubr"]
ab = subs[subs.tier.isin(["A", "B"])]
print(ab[show].to_string(index=False, formatters={
    k: "{:+.2f}".format for k in ["dpsi_zyg11b", "dpsi_zer1", "dpsi_dko"]}))
