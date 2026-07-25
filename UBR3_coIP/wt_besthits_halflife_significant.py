"""
Which UBR3-WT best hits have a half-life worth quoting?

Two things have to be true before a half-life means anything here:

  1. The decay is real. The log-linear fit over the 8 WT untreated points must have
     a slope significantly below zero (two-sided p on the slope, Benjamini-Hochberg
     corrected across the 45 quantified hits). A flat noisy line gives a huge
     "half-life" that is an artefact of dividing by a near-zero slope.
  2. The protein is short-lived relative to the proteome. Its decay rate k is
     compared with the 8097-protein background distribution; a protein at the 50th
     percentile is simply behaving like bulk protein.

Reported separately, because it is a different kind of confidence:
  * measured   - the curve drops below 50% inside the 15 h assay, so t1/2 is read
                 from data rather than extrapolated.
  * extrapolated - decay is significant but t1/2 lies beyond 15 h.

Output: WT_besthits_halflife_significant.csv  and a sheet in the master workbook.
"""
import pathlib

import numpy as np
import pandas as pd
from scipy import stats

HERE = pathlib.Path(__file__).parent

TIMES = np.array([0.0, 0.0, 5.0, 5.0, 10.0, 10.0, 15.0, 15.0])
RCOLS = ["WT_UT1_ratio", "WT_UT2_ratio", "WT_UT5h1_ratio", "WT_UT5h2_ratio",
         "WT_UT10h1_ratio", "WT_UT10h2_ratio", "WT_UT15h1_ratio", "WT_UT15h2_ratio"]
FDR_Q = 0.05
FAST_PCTL = 80          # decay rate must beat this percentile of the proteome
                        # (top 20% of proteome turnover)

curves = pd.read_csv(HERE / "WT_besthits_halflife_curves.csv")
bg = pd.read_csv(HERE / "WT_besthits_halflife_background.csv")
master = pd.read_csv(HERE / "WT_besthits_master_list.csv")
info = master.set_index("Gene")


def fit(y):
    """Slope test on ln(ratio) vs time."""
    if np.isnan(y).any() or (y <= 0).any():
        return np.nan, np.nan, np.nan, np.nan
    r = stats.linregress(TIMES, np.log(y))
    k = -r.slope
    t_half = np.log(2) / k if k > 0 else np.inf
    return k, t_half, r.pvalue, r.rvalue ** 2


rows = []
for c in curves.itertuples():
    y = np.array([getattr(c, col) for col in RCOLS], dtype=float)
    k, t_half, p, r2 = fit(y)
    rows.append({"Gene": c.gene, "k": k, "t_half": t_half, "slope_p": p, "r2": r2})
fits = pd.DataFrame(rows)

# background decay-rate distribution, for the "faster than bulk protein" test
bg_k = bg["decay_k_per_h"].replace([np.inf, -np.inf], np.nan).dropna()
k_cut = np.percentile(bg_k, FAST_PCTL)
fits["percentile_vs_proteome"] = fits["k"].apply(
    lambda v: float(stats.percentileofscore(bg_k, v)) if np.isfinite(v) else np.nan)

# Benjamini-Hochberg across the quantified hits
ok = fits["slope_p"].notna()
p = fits.loc[ok, "slope_p"].to_numpy()
order = np.argsort(p)
ranked = p[order]
q = ranked * len(ranked) / (np.arange(len(ranked)) + 1)
q = np.minimum.accumulate(q[::-1])[::-1]
qv = np.empty_like(q)
qv[order] = np.clip(q, 0, 1)
fits.loc[ok, "slope_q"] = qv

fits["decay_is_real"] = (fits["slope_q"] < FDR_Q) & (fits["k"] > 0)
fits["faster_than_proteome"] = fits["percentile_vs_proteome"] >= FAST_PCTL
fits["significant"] = fits["decay_is_real"] & fits["faster_than_proteome"]
fits["evidence"] = np.where(
    fits["t_half"] <= 15, "measured (crosses 50% inside the assay)",
    "extrapolated (still >50% at 15 h)")

sig = fits[fits["significant"]].copy().sort_values("t_half").reset_index(drop=True)
sig["Half-life (h)"] = sig["t_half"].round(1)
sig["Decay rate k (1/h)"] = sig["k"].round(4)
sig["Slope q (BH)"] = sig["slope_q"].apply(lambda v: f"{v:.2g}")
sig["Fit R2"] = sig["r2"].round(2)
sig["Percentile vs proteome"] = sig["percentile_vs_proteome"].round(1)
sig["Remaining at 15 h"] = [info.loc[g, "Remaining at 15 h"] for g in sig["Gene"]]
sig["WT rank"] = [int(info.loc[g, "WT_rank"]) for g in sig["Gene"]]
sig["Module"] = [info.loc[g, "Pathway module(s)"] if isinstance(
    info.loc[g, "Pathway module(s)"], str) else "" for g in sig["Gene"]]
sig["Role"] = [info.loc[g, "Cellular role"] for g in sig["Gene"]]

out_cols = ["WT rank", "Gene", "Half-life (h)", "Remaining at 15 h", "Decay rate k (1/h)",
            "Slope q (BH)", "Fit R2", "Percentile vs proteome", "evidence", "Module", "Role"]
sig_out = sig[out_cols].rename(columns={"evidence": "Evidence"})
sig_out.insert(0, "#", range(1, len(sig_out) + 1))
sig_out.to_csv(HERE / "WT_besthits_halflife_significant.csv", index=False)

fits.to_csv(HERE / "WT_besthits_halflife_slopetests.csv", index=False)

# add / refresh the sheet in the master workbook
book = HERE / "WT_besthits_master_list.xlsx"
if book.exists():
    with pd.ExcelWriter(book, engine="openpyxl", mode="a",
                        if_sheet_exists="replace") as xl:
        sig_out.to_excel(xl, sheet_name="Half-life significant", index=False)
        ws = xl.sheets["Half-life significant"]
        ws.freeze_panes = "A2"
        for col, w in {"C": 12, "D": 14, "E": 18, "F": 18, "G": 13, "H": 9,
                       "I": 22, "J": 40, "K": 10, "L": 38}.items():
            ws.column_dimensions[col].width = w

print(f"quantified hits tested       : {int(ok.sum())}")
print(f"decay significant (BH q<{FDR_Q}) : {int(fits['decay_is_real'].sum())}")
print(f"  and faster than {FAST_PCTL}th pctile of proteome (k > {k_cut:.4f}/h): "
      f"{int(fits['significant'].sum())}")
print(f"  of those, measured inside the 15 h window: "
      f"{int((sig['t_half'] <= 15).sum())}")
print("\n=== SIGNIFICANT HALF-LIVES ===")
print(sig_out.to_string(index=False))

print("\n--- significant decay but NOT faster than the proteome (excluded) ---")
excl = fits[fits["decay_is_real"] & ~fits["faster_than_proteome"]].sort_values("t_half")
print(", ".join(f"{r.Gene} ({r.t_half:.0f} h, p{r.percentile_vs_proteome:.0f})"
                for r in excl.itertuples()))
