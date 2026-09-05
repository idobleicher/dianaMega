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

# ---------------------------------------------------------------- palette
# TWO hues, blue and pink. This project's figures do NOT share the warm ramp of
# ../ubr3_PG_reviewer/ -- they are a different paper and are painted apart on
# purpose.
#
# Purple was a third category here and is not any more. It is the hue that
# cannot be separated from blue for a red-blind reader -- purple is blue plus
# red, and red is exactly what a protanope does not see -- so a blue/purple
# pair sat at dE 6.2, inside the band that is only legal where the mark carries
# its own label. Two hues remove the problem instead of managing it: blue
# against pink is dE 10.5 under every simulated colour vision and dE 27.2 under
# normal vision, comfortably clear of every floor, and it needs no relief and
# no caveat. Purple now appears only as the short passage the sequential ramp
# makes between the two, where nothing depends on telling it from its
# neighbours.
#
# The residues are grouped by CHARGE, which is the axis this paper's result is
# about (acidic content, net charge). The 15 uncharged residues take a neutral
# grey rather than a third hue: grey is not an identity colour, it is the
# absence of the property the other two encode, which is exactly what those
# residues are. A logo glyph also names itself, so no residue depends on its
# colour to be read.
CATEGORY_MEMBERS = {
    "Acidic": list("DE"),
    "Basic": list("KHR"),
    "Uncharged": list("STCNQGPAMVLIFYW"),
}
CAT_ORDER = ["Acidic", "Basic", "Uncharged"]
CAT_COLORS = {"Acidic": "#2E6FD0", "Basic": "#D6408F", "Uncharged": "#A6A4AE"}
CAT_OF_AA = {a: c for c, mem in CATEGORY_MEMBERS.items() for a in mem}

# Sequential ramp for -log10 p: pale blue -> blue -> a brief violet -> magenta
# -> deep pink. Verified monotonic in relative luminance, every step decreasing
# and the smallest step 0.021 -- a magnitude ramp that stalls in lightness
# stops encoding magnitude, and a blue-to-pink ramp will stall at the magenta
# end unless it is pushed dark, which is why the last stops are as deep as they
# are.
CMAP_STOPS = ["#F5F8FD", "#DDE8F8", "#B6CDF2", "#86A8E8", "#5C82DC",
              "#6A5FC8", "#8E3C9E", "#9E2266", "#6E0E3C"]
INK, MUTED, RULE, GRID = "#1A1A1A", "#6B6B6B", "#C6C6CE", "#FFFFFF"

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
