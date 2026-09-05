"""
Position x residue comparison WITHIN the Met-Pro peptides.

The question this answers: among peptides that start Met-Pro, is anything
BESIDES the proline different between the ones stabilised on loss of
ZYG11B / ZER1 / the double KO and the ones that are not? Every peptide on both
sides of this comparison has the same position 1 (Met) and the same position 2
(Pro), so the proline itself cannot separate them and positions 3-24 are what
is being asked about.

This is the design the UBR3 package uses in `../ubr3_PG_reviewer/`: substrates
against non-substrates CARRYING THE SAME MOTIF, never against the whole
library, because a comparison against the library re-discovers the motif and
says nothing about what else matters.

THE TWO GROUPS
  substrates  Met-Pro, and dPSI >= 1.0 with headroom in at least one of
              ZYG11B KO, ZER1 KO, the double KO.                      n = 36
  controls    Met-Pro, had headroom in at least one experiment, and dPSI < 0.5
              in every genotype measured for it.                      n = 424

  The 662 Met-Pro peptides between those two definitions are in NEITHER group.
  Dropping the middle is deliberate: with 36 substrates, letting peptides at
  dPSI 0.5-1.0 into the control group would put probable weak substrates into
  the background and flatten exactly the signal being looked for. It also means
  the control group is not "all other Met-Pro peptides", and the fold changes
  here are not comparable to fold changes computed against the full library.

TWO WARNINGS THAT BELONG ON EVERY FIGURE BUILT FROM THIS
  * n = 36. A residue must appear in 4-5 substrates against a 2% background
    before it reaches p < 0.01, and a single peptide moves any cell by 2.8
    percentage points. Across 22 positions x 20 residues there are 440 cells,
    so about 22 reach p < 0.05 by chance alone.
  * 14 of the 36 substrates are ZER1-only calls, and the ZER1 hit list carries
    no glycine signal at all (see README) -- it may be largely noise. Every
    figure is therefore also built for the ZYG11B + double-KO subset alone
    (SUBSET_STRICT), and a result that only survives in the union of all three
    should not be believed.

Met excision is inefficient at Met-Pro, so these peptides probably keep their
initiator Met. Their exposed N-terminus is Met, not Pro. Whatever separates the
two groups here is therefore NOT necessarily an N-degron.
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

import zz_gps_data as zz

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

POS = list(range(3, 25))        # 1 is Met, 2 is Pro, both fixed by definition
AA = list("ACDEFGHIKLMNPQRSTVWY")
CUT_HIT, CUT_CTRL = 1.0, 0.5

# Chemical classes and the warm palette are carried over verbatim from
# ../ubr3_PG_reviewer/pg_motif_data.py so the two projects' logos read as one
# visual system. Copied rather than imported: the projects are separate, and a
# folder name with a space in it is not importable.
CATEGORY_MEMBERS = {
    "Basic": list("KHR"), "Acidic": list("DE"), "Aromatic": list("FYW"),
    "Aliphatic": list("GPAM"), "Hydrophobic": list("VLI"), "Polar": list("STCNQ"),
}
CAT_ORDER = ["Basic", "Acidic", "Aromatic", "Aliphatic", "Hydrophobic", "Polar"]
CAT_COLORS = {"Acidic": "#7B241C", "Aromatic": "#A0522D", "Basic": "#C0392B",
              "Hydrophobic": "#D35400", "Aliphatic": "#E8A33D", "Polar": "#F1C40F"}
CAT_OF_AA = {a: c for c, mem in CATEGORY_MEMBERS.items() for a in mem}

CMAP_STOPS = ["#F7F4EC", "#F4EAD0", "#F1DCA6", "#EDC76E", "#E8A33D",
              "#DC7C31", "#C0392B", "#9B2A20", "#71201A"]
INK, MUTED, RULE, GRID = "#1A1A1A", "#6B6B6B", "#C9C4BA", "#FFFFFF"

_CACHE = {}


def stars(p):
    if p is None or not np.isfinite(p):
        return ""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


def groups(strict=False):
    """(substrates, controls) as frames of Met-Pro peptides.

    strict=True drops ZER1 from the hit definition, leaving ZYG11B and the
    double KO -- the two readouts that carry a glycine signal.
    """
    d = zz.load()
    pro = d[d.nterm == "P"].copy()

    hit = (pro.headroom_b.fillna(False) & (pro.dpsi_dko >= CUT_HIT).fillna(False)) \
        | (pro.headroom.fillna(False) & (pro.dpsi_zyg11b >= CUT_HIT).fillna(False))
    if not strict:
        hit = hit | (pro.headroom.fillna(False) & (pro.dpsi_zer1 >= CUT_HIT).fillna(False))

    dp = pro[["dpsi_zyg11b", "dpsi_zer1", "dpsi_dko"]]
    quiet = (dp.max(axis=1) < CUT_CTRL) & (pro.headroom.fillna(False)
                                           | pro.headroom_b.fillna(False))
    return pro[hit].copy(), pro[~hit & quiet].copy()


def cells(strict=False):
    """One row per position x residue: counts, percentages, fold change,
    Fisher exact p and BH-FDR q across all 440 cells."""
    key = ("strict" if strict else "union")
    if key in _CACHE:
        return _CACHE[key].copy()

    sub, ctrl = groups(strict)
    n_sub, n_ctrl = len(sub), len(ctrl)
    s_seq, c_seq = sub.seq.values, ctrl.seq.values

    rows = []
    for p in POS:
        s_col = np.array([q[p - 1] for q in s_seq])
        c_col = np.array([q[p - 1] for q in c_seq])
        for aa in AA:
            a = int((s_col == aa).sum())
            c = int((c_col == aa).sum())
            f_s, f_c = a / n_sub, c / n_ctrl
            fc = np.nan if (a == 0 and c == 0) else (np.inf if c == 0 else f_s / f_c)
            _, p_f = fisher_exact([[a, n_sub - a], [c, n_ctrl - c]])
            rows.append({"position": p, "label": aa, "n_sub": a, "n_ctrl": c,
                         "pct_sub": 100 * f_s, "pct_ctrl": 100 * f_c,
                         "fold_change": fc,
                         "log2FC": np.log2(fc) if np.isfinite(fc) and fc > 0 else np.nan,
                         "direction": 0 if f_s == f_c else (1 if f_s > f_c else -1),
                         "p_fisher": p_f, "defined": not (a == 0 and c == 0)})

    t = pd.DataFrame(rows)
    m = t.defined
    t.loc[m, "q_fisher"] = multipletests(t.loc[m, "p_fisher"].values, method="fdr_bh")[1]
    t["n_group_sub"], t["n_group_ctrl"] = n_sub, n_ctrl
    _CACHE[key] = t
    return t.copy()


def matrix(col, order=None, table=None):
    table = cells() if table is None else table
    return (table.pivot(index="label", columns="position", values=col)
                 .reindex(index=order or AA, columns=POS))


def logo_frame(col="fold_change", table=None, floor=1.0, only=None):
    """positions x residues, shaped for logomaker. Values at or below `floor`
    become 0, so a fold-change logo draws only what is enriched."""
    m = matrix(col, table=table).astype(float)
    m = m.where(np.isfinite(m), 0.0)
    if floor is not None:
        m = m.where(m > floor, 0.0)
    if only is not None:
        m = m.where(only.reindex_like(m).fillna(False), 0.0)
    return m.T.set_axis(POS)


if __name__ == "__main__":
    for strict in (False, True):
        sub, ctrl = groups(strict)
        t = cells(strict)
        k = t[t.defined]
        tag = "ZYG11B + DKO only" if strict else "ZYG11B / ZER1 / DKO"
        print(f"\n=== {tag}:  {len(sub)} substrates vs {len(ctrl)} controls, all Met-Pro")
        print(f"    {len(k)} testable cells of {len(t)};  {0.05 * len(k):.0f} expected at p<0.05 "
              f"by chance;  {int((k.p_fisher < 0.05).sum())} observed;  "
              f"{int((k.q_fisher < 0.05).sum())} survive BH-FDR")
        top = k[k.direction > 0].nsmallest(8, "p_fisher")
        print(top[["position", "label", "n_sub", "n_ctrl", "pct_sub", "pct_ctrl",
                   "fold_change", "p_fisher", "q_fisher"]].to_string(
            index=False, formatters={"pct_sub": "{:.1f}".format, "pct_ctrl": "{:.1f}".format,
                                     "fold_change": "{:.2f}".format,
                                     "p_fisher": "{:.2e}".format, "q_fisher": "{:.3f}".format}))
