"""
Single loader for the P/G-D/E/T motif position matrices.

Every figure that plots positions 4-24 -- the three logos and the significance
heatmaps -- imports this module, so they cannot drift apart and cannot pick up
a stale copy of the source data.

Source: data/AminoAcids_PG_motif_with_pvalues.csv
  n = 20 P/G-D/E/T substrates  vs  n = 193 non-substrate controls carrying the
  same motif. Positions 4-24; positions 1-3 are the motif itself.

  The sheet carries, for every residue x position and every chemical class x
  position, the full 2 x 2 contingency table (with / without the feature, in
  each group) and a chi-square p-value. Fisher exact is recomputed here from
  the same tables: the workbook's own expected-count block shows most cells
  expect fewer than 2 substrates, well under the >= 5 that chi-square needs,
  so both tests are carried and the caller chooses.

WHY NOT data/AMINO_ACIDS_PG_motif.csv
  That earlier export stores significance only as * / ** / *** bins, and its
  arginine row is byte-identical to the Basic-class row -- it is not arginine
  at all. The other 19 residue rows match this file exactly. Figures 13-15 were
  originally built from it, so their R glyphs carried class-level values
  (R at position 7 drawn at 2.89x; the true value is 5.63x). The file is kept
  in data/ for provenance only. Do not read it.
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "data", "AminoAcids_PG_motif_with_pvalues.csv")

POS = list(range(4, 25))
N_SUB, N_CTRL = 20, 193

# display order = the Colab aa_order, which is already grouped by class
AA = ["K", "H", "R", "D", "E", "F", "Y", "W", "G", "P",
      "A", "M", "V", "L", "I", "S", "T", "C", "N", "Q"]

CATEGORY_MEMBERS = {
    "Basic":       list("KHR"),
    "Acidic":      list("DE"),
    "Aromatic":    list("FYW"),
    "Aliphatic":   list("GPAM"),
    "Hydrophobic": list("VLI"),
    "Polar":       list("STCNQ"),
}
CAT_ORDER = ["Basic", "Acidic", "Aromatic", "Aliphatic", "Hydrophobic", "Polar"]
CAT_COLORS = {
    "Acidic":      "#7B241C",
    "Aromatic":    "#A0522D",
    "Basic":       "#C0392B",
    "Hydrophobic": "#D35400",
    "Aliphatic":   "#E8A33D",
    "Polar":       "#F1C40F",
}
CAT_OF_AA = {a: c for c, mem in CATEGORY_MEMBERS.items() for a in mem}
# the workbook's spellings -> ours
CAT_IN_FILE = {"BASIC": "Basic", "Acidic": "Acidic", "Aromtic": "Aromatic",
               "Alipatic": "Aliphatic", "Hydrophobic": "Hydrophobic",
               "Polar": "Polar"}

assert sorted(sum(CATEGORY_MEMBERS.values(), [])) == sorted(AA), \
    "the class grouping must be a partition of all 20 residues"

_CACHE = None


def _num(x):
    x = str(x).strip()
    if x in ("", "N/A", "-", "#NUM!", "nan", "#DIV/0!"):
        return np.nan
    try:
        return float(x)
    except ValueError:
        return np.nan


def stars(p):
    """* p<0.05, ** p<0.01, *** p<0.001 -- computed from the p value itself,
    never read from the workbook's asterisk matrix, which is not internally
    consistent with its own p values."""
    if p is None or not np.isfinite(p):
        return ""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


def load(force=False):
    """Tidy one-row-per-cell table for all 420 residue and 126 class cells.

    Columns: kind, label, position, n_sub, pct_sub, n_ctrl, pct_ctrl,
             fold_change, log2FC, direction, defined,
             p_chi2, p_fisher, q_chi2, q_fisher

    fold_change is inf where the feature is absent from the controls and NaN
    where it is absent from both groups (`defined` is False there -- no test is
    possible, and the workbook shows N/A). log2FC is NaN whenever fold_change
    is not finite and positive, matching the workbook.

    BH-FDR is applied within each family (residues, classes) and separately for
    each test.
    """
    global _CACHE
    if _CACHE is not None and not force:
        return _CACHE.copy()

    raw = pd.read_csv(SRC, header=None, dtype=str, keep_default_na=False)

    def find_row(col, text):
        for i in range(len(raw)):
            if str(raw.iloc[i, col]).strip() == text:
                return i
        raise LookupError(f"{text!r} not found in column {col} of {SRC}")

    chi_hdr = find_row(2, "CHI test")        # p-value block anchor
    obs_hdr = find_row(3, "Substrates")      # first row of the observed block

    # --- observed 2x2 tables: label col 2, arm col 3, two columns per position
    counts, i = {}, obs_hdr
    while i < chi_hdr:
        lab = str(raw.iloc[i, 2]).strip()
        if str(raw.iloc[i, 3]).strip() == "Substrates" and lab:
            sub = [_num(v) for v in raw.iloc[i, 4:4 + 2 * len(POS)]]
            non = [_num(v) for v in raw.iloc[i + 1, 4:4 + 2 * len(POS)]]
            # a second block repeats the same labels with EXPECTED counts;
            # keep the first (observed) occurrence only
            counts.setdefault(lab, pd.DataFrame(
                {"a": sub[0::2], "b": sub[1::2], "c": non[0::2], "d": non[1::2]},
                index=POS))
            i += 2
        else:
            i += 1

    # --- the workbook's chi-square p-values (header row is chi_hdr + 1)
    chi_p, i = {}, chi_hdr + 2
    while i < len(raw) and str(raw.iloc[i, 2]).strip():
        chi_p[str(raw.iloc[i, 2]).strip()] = pd.Series(
            [_num(v) for v in raw.iloc[i, 3:3 + len(POS)]], index=POS)
        i += 1

    missing = [l for l in list(AA) + list(CAT_IN_FILE) if l not in counts]
    assert not missing, f"no contingency table for: {missing}"
    for lab, t in counts.items():
        ok = ((t.a + t.b == N_SUB) | t.a.isna()) & ((t.c + t.d == N_CTRL) | t.c.isna())
        assert ok.all(), f"{lab}: contingency rows do not sum to {N_SUB}/{N_CTRL}"

    rows = []
    for lab, t in counts.items():
        is_cat = lab in CAT_IN_FILE
        name = CAT_IN_FILE[lab] if is_cat else lab
        for p in POS:
            a, b, c, d = (float(t.at[p, k]) for k in ("a", "b", "c", "d"))
            f_sub, f_non = a / N_SUB, c / N_CTRL
            fc = np.nan if (a == 0 and c == 0) else (np.inf if c == 0 else f_sub / f_non)
            _, p_fisher = fisher_exact([[int(a), int(b)], [int(c), int(d)]])
            rows.append({
                "kind": "category" if is_cat else "residue",
                "label": name, "position": p,
                "n_sub": int(a), "pct_sub": 100 * f_sub,
                "n_ctrl": int(c), "pct_ctrl": 100 * f_non,
                "fold_change": fc,
                "log2FC": np.log2(fc) if np.isfinite(fc) and fc > 0 else np.nan,
                "direction": 0 if f_sub == f_non else (1 if f_sub > f_non else -1),
                "p_chi2": chi_p.get(lab, pd.Series(dtype=float)).get(p, np.nan),
                "p_fisher": p_fisher,
                "defined": not (a == 0 and c == 0),
            })

    cells = pd.DataFrame(rows)
    for kind in ("residue", "category"):
        for col in ("p_chi2", "p_fisher"):
            m = (cells.kind == kind) & cells[col].notna()
            cells.loc[m, col.replace("p_", "q_")] = multipletests(
                cells.loc[m, col].values, method="fdr_bh")[1]

    _CACHE = cells
    return cells.copy()


def matrix(col, kind="residue", order=None, cells=None):
    """`col` as a DataFrame with rows = residues/classes and columns = POS."""
    cells = load() if cells is None else cells
    order = order or (AA if kind == "residue" else CAT_ORDER)
    return (cells[cells.kind == kind]
            .pivot(index="label", columns="position", values=col)
            .reindex(index=order, columns=POS))


def logo_frame(col="fold_change", kind="residue", order=None, cells=None,
               floor=1.0, only=None):
    """Matrix shaped for logomaker: rows = positions, columns = residues.

    Non-finite values become 0 (the workbook shows them as N/A), as do values
    at or below `floor` -- for fold change that drops everything not enriched,
    exactly as the original logo scripts did. `only` is an optional boolean
    matrix in the same shape as `matrix()`; cells outside it are zeroed.
    """
    m = matrix(col, kind, order, cells).astype(float)
    m = m.where(np.isfinite(m), 0.0)
    if floor is not None:
        m = m.where(m > floor, 0.0)
    if only is not None:
        m = m.where(only.reindex_like(m).fillna(False), 0.0)
    return m.T.set_axis(POS)


if __name__ == "__main__":
    c = load()
    print(f"{len(c)} cells  ({int((c.kind == 'residue').sum())} residue + "
          f"{int((c.kind == 'category').sum())} class)")
    for kind in ("residue", "category"):
        k = c[(c.kind == kind) & c.defined]
        print(f"\n{kind}: {len(k)} testable, "
              f"{0.05 * len(k):.0f} expected at p<0.05 by chance")
        for col in ("p_chi2", "p_fisher"):
            q = col.replace("p_", "q_")
            print(f"   {col:<9} p<0.05: {int((k[col] < 0.05).sum()):>3}"
                  f"   q<0.05: {int((k[q] < 0.05).sum()):>3}")
