"""
What separates stabilised peptides from unstabilised ones, once the N-terminal
residue is held fixed.

The per-position logos (make_pro_logo.py) find nothing: 25 of 440 cells reach
p < 0.05 against 22 expected by chance, none survives FDR, and the two versions
of the figure disagree about which cells are tall. That is the correct answer
to "is there a second motif position" and it is no.

It is the wrong test for "is there anything else". A per-cell test spends its
power on 440 hypotheses to ask about one residue at one position, when the
thing that actually differs is spread over the whole sequence. Asking about
COMPOSITION instead -- 27 tests rather than 440, each pooling all 22 positions
-- finds it immediately, in both motif classes, at effect sizes that are not
subtle:

    D+E content      Met-Pro   13.3% stabilised vs  7.3% not
                     Met-Gly   13.1% stabilised vs  6.2% not
    hydropathy       Met-Pro   -0.20 stabilised vs +0.24 not
    net charge       Met-Pro   -0.75 stabilised vs +0.50 not

THE CONTROL THAT MATTERS. Acidic peptides could be stabilised-looking simply by
sitting lower at baseline with more room to rise. They do not: among controls,
acidic content correlates POSITIVELY with baseline PSI (Spearman r = +0.22,
p = 8e-6), so the ceiling works against this result rather than producing it.
Compared band by band within matched baseline PSI the difference holds anyway --
combined p = 0.0018 for Met-Pro and 7e-19 for Met-Gly. Panels D and E are that
control, and they are the reason the conclusion is stated at all.

WHAT IT IS NOT. This is not a proline-specific second determinant. The same
composition separates stabilised from unstabilised Met-GLY peptides, more
strongly, and those are the canonical ZYG11B substrates. It is a property of
what makes any peptide a substrate in this screen -- acidic, hydrophilic,
negatively charged N-terminal regions -- not of how proline is read.

Colours: #9B2A20 (stabilised) and #D9A441 (not), both from the project ramp,
checked with the dataviz validator -- CVD dE 29, well clear of the floor. The
gold sits below 3:1 against white, so every group is labelled on the axis and
identity is never carried by colour alone.

Outputs:
  figures/FigureZ3_composition_stabilised_vs_not   the five panels
  data/pro_composition_tests.csv                   all 27 tests, both motifs
"""
import os

import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr, combine_pvalues
from statsmodels.stats.multitest import multipletests

import zz_gps_data as zz
import zz_pro_motif as zp

HERE = os.path.dirname(os.path.abspath(__file__))
FIG, DATA = os.path.join(HERE, "figures"), os.path.join(HERE, "data")
os.makedirs(FIG, exist_ok=True)

HIT, QUIET = "#9B2A20", "#D9A441"
INK, MUTED, RULE = zp.INK, zp.MUTED, zp.RULE

KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2,
          I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9,
          Y=-1.3, V=4.2)
TAIL = slice(2, 24)          # positions 3-24: everything the motif does not fix

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8, "axes.linewidth": 0.8,
    "svg.fonttype": "none", "pdf.fonttype": 42,
    "figure.dpi": 600, "savefig.dpi": 600, "savefig.bbox": "tight",
    "figure.facecolor": "white", "savefig.facecolor": "white",
})


def props(seq):
    t = seq[TAIL]
    n = len(t)
    return pd.Series({
        "acidic": 100 * sum(t.count(a) for a in "DE") / n,
        "basic": 100 * sum(t.count(a) for a in "KR") / n,
        "net_charge": sum(t.count(a) for a in "KR") - sum(t.count(a) for a in "DE"),
        "gravy": float(np.mean([KD[c] for c in t if c in KD])),
        "aromatic": 100 * sum(t.count(a) for a in "FWY") / n,
    })


# ------------------------------------------------------------------ the sets
d = zz.load()
pro_sub, pro_ctrl = zp.groups(strict=False)

gly = d[d.is_gly & d.headroom_b.fillna(False) & d.psi_ctrl_b.notna()]
gly_sub, gly_ctrl = gly[gly.dpsi_dko >= 1.0], gly[gly.dpsi_dko < 0.5]

SETS = [("Met-Pro", pro_sub, pro_ctrl, "psi_ctrl_b", "psi_ctrl_a"),
        ("Met-Gly", gly_sub, gly_ctrl, "psi_ctrl_b", "psi_ctrl_b")]

frames = {}
for name, s, c, _, _ in SETS:
    for grp, f in (("stabilised", s), ("not stabilised", c)):
        t = f.seq.apply(props)
        t["motif"], t["group"] = name, grp
        t["base"] = f.psi_ctrl_b.fillna(f.psi_ctrl_a).values
        frames[(name, grp)] = t

# ------------------------------------------------------------------ the tests
rows = []
for name, *_ in SETS:
    S, C = frames[(name, "stabilised")], frames[(name, "not stabilised")]
    for col in ["acidic", "basic", "net_charge", "gravy", "aromatic"]:
        _, p = mannwhitneyu(S[col], C[col])
        rows.append({"motif": name, "test": col, "stabilised": S[col].mean(),
                     "not_stabilised": C[col].mean(), "n_sub": len(S),
                     "n_ctrl": len(C), "p": p})
    sub_seq = (pro_sub if name == "Met-Pro" else gly_sub).seq
    ctl_seq = (pro_ctrl if name == "Met-Gly" and False else
               (pro_ctrl if name == "Met-Pro" else gly_ctrl)).seq
    for aa in zp.AA:
        s = sub_seq.map(lambda q: q[TAIL].count(aa) / 22)
        c_ = ctl_seq.map(lambda q: q[TAIL].count(aa) / 22)
        _, p = mannwhitneyu(s, c_)
        rows.append({"motif": name, "test": f"{aa} content", "stabilised": 100 * s.mean(),
                     "not_stabilised": 100 * c_.mean(), "n_sub": len(s),
                     "n_ctrl": len(c_), "p": p})
tests = pd.DataFrame(rows)
for name in tests.motif.unique():
    m = tests.motif == name
    tests.loc[m, "q"] = multipletests(tests.loc[m, "p"], method="fdr_bh")[1]
tests.to_csv(os.path.join(DATA, "pro_composition_tests.csv"), index=False,
             float_format="%.5g")

# the baseline-matched control
BANDS = [(0, 2.5), (2.5, 3.0), (3.0, 3.5), (3.5, 4.1)]
matched = []
for name, *_ in SETS:
    S, C = frames[(name, "stabilised")], frames[(name, "not stabilised")]
    for lo, hi in BANDS:
        s = S[(S.base >= lo) & (S.base < hi)]
        c = C[(C.base >= lo) & (C.base < hi)]
        if len(s) >= 4 and len(c) >= 15:
            _, p = mannwhitneyu(s.acidic, c.acidic)
            matched.append({"motif": name, "band": f"{lo:g}–{hi:g}", "lo": lo,
                            "n_sub": len(s), "n_ctrl": len(c),
                            "acidic_sub": s.acidic.mean(),
                            "acidic_ctrl": c.acidic.mean(), "p": p})
matched = pd.DataFrame(matched)


# ------------------------------------------------------------------ drawing
def strip_panel(ax, col, label, motifs=("Met-Pro", "Met-Gly")):
    """Box + jittered points, one pair per motif class. Thin marks, one accent
    against one recessive gold; the axis labels carry identity, not colour."""
    rng = np.random.default_rng(0)
    xs, ticks, labels = [], [], []
    x = 0
    for motif in motifs:
        for grp, colour in (("not stabilised", QUIET), ("stabilised", HIT)):
            v = frames[(motif, grp)][col].values
            bp = ax.boxplot([v], positions=[x], widths=0.52, showfliers=False,
                            patch_artist=True, zorder=3)
            for box in bp["boxes"]:
                box.set(facecolor="white", edgecolor=colour, linewidth=1.0)
            for part in ("whiskers", "caps"):
                for ln in bp[part]:
                    ln.set(color=colour, linewidth=0.9)
            for md in bp["medians"]:
                md.set(color=colour, linewidth=1.6)
            n = len(v)
            size, alpha = (5.0, 0.55) if n < 200 else (1.6, 0.18)
            ax.scatter(x + rng.normal(0, 0.085, n), v, s=size, color=colour,
                       alpha=alpha, linewidths=0, zorder=2)
            ticks.append(x); labels.append(f"{grp}\nn = {n}")
            xs.append((motif, grp, x, v))
            x += 1
        x += 0.7

    # the two comparisons, annotated where they sit
    top = max(np.percentile(v, 97) for *_, v in xs)
    bot = min(np.percentile(v, 3) for *_, v in xs)
    span = top - bot
    for i, motif in enumerate(motifs):
        a, b = xs[2 * i][2], xs[2 * i + 1][2]
        p = tests[(tests.motif == motif) & (tests.test == col)].p.iloc[0]
        y = top + span * (0.10 + 0.0)
        ax.plot([a, a, b, b], [y, y + span * 0.03, y + span * 0.03, y],
                color=MUTED, linewidth=0.7, clip_on=False)
        txt = "p < 1e-4" if p < 1e-4 else f"p = {p:.3g}"
        ax.text((a + b) / 2, y + span * 0.055, txt, ha="center", va="bottom",
                fontsize=6.0, color=INK)
        # the motif name goes UNDER the two-line tick labels, in axes
        # fractions, so it never has to be given room inside the data range
        ax.text((a + b) / 2, -0.145, motif, ha="center", va="top",
                transform=ax.get_xaxis_transform(),
                fontsize=7.2, color=INK, fontweight="bold")

    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=5.8, color=MUTED)
    ax.set_ylabel(label, fontsize=7.5, color=INK)
    # a percentage cannot be negative, so the axis does not go there
    floor_at_zero = min(v.min() for *_, v in xs) >= 0
    lo = 0.0 if floor_at_zero else bot - span * 0.12
    ax.set_ylim(lo, top + span * 0.22)
    ax.tick_params(axis="y", labelsize=6.5, length=2.5, width=0.7, colors=INK)
    ax.tick_params(axis="x", length=0, pad=2)
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color(RULE)
    ax.yaxis.grid(True, color="#EFECE6", linewidth=0.6)
    ax.set_axisbelow(True)


def band_panel(ax, motif):
    """Panel D/E: the same comparison inside matched baseline-PSI bands, which
    is what rules the ceiling out as the explanation."""
    m = matched[matched.motif == motif]
    x = np.arange(len(m))
    ax.plot(np.repeat(x, 3).reshape(-1, 3).T[:2],
            np.vstack([m.acidic_ctrl, m.acidic_sub]), color=RULE, linewidth=0.8,
            zorder=1)
    ax.scatter(x, m.acidic_ctrl, s=26, color=QUIET, zorder=3,
               label="not stabilised", linewidths=0)
    ax.scatter(x, m.acidic_sub, s=26, color=HIT, zorder=3,
               label="stabilised", linewidths=0)
    for xi, r in zip(x, m.itertuples()):
        star = "***" if r.p < 0.001 else "**" if r.p < 0.01 else "*" if r.p < 0.05 else "n.s."
        ax.text(xi, max(r.acidic_sub, r.acidic_ctrl) + 0.9, star, ha="center",
                va="bottom", fontsize=6.2, color=INK,
                fontweight="bold" if star != "n.s." else "normal")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r.band}\n{r.n_sub} vs {r.n_ctrl}" for r in m.itertuples()],
                       fontsize=5.8, color=MUTED)
    ax.set_xlim(-0.5, len(m) - 0.5)
    ax.set_ylim(0, max(m.acidic_sub.max(), m.acidic_ctrl.max()) + 4.2)
    ax.set_xlabel("Baseline PSI band   (peptides: stabilised vs not)", fontsize=6.8,
                  color=MUTED, labelpad=2)
    ax.set_ylabel("D+E content (%)", fontsize=7.5, color=INK)
    p_comb = combine_pvalues(m.p, method="stouffer")[1]
    ax.set_title(f"{motif} — matched for baseline stability   "
                 f"(combined $p$ = {p_comb:.2g})", fontsize=7.2, color=INK, pad=4)
    ax.tick_params(axis="y", labelsize=6.5, length=2.5, width=0.7, colors=INK)
    ax.tick_params(axis="x", length=0, pad=2)
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color(RULE)
    ax.yaxis.grid(True, color="#EFECE6", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.legend(fontsize=6, frameon=False, loc="upper left", handletextpad=0.3,
              borderpad=0.2, labelcolor=INK)


fig = plt.figure(figsize=(11.0, 5.9))
gs = fig.add_gridspec(2, 6, height_ratios=[1.0, 0.82], hspace=0.42, wspace=0.55)
axA = fig.add_subplot(gs[0, 0:2]); axB = fig.add_subplot(gs[0, 2:4])
axC = fig.add_subplot(gs[0, 4:6])
axD = fig.add_subplot(gs[1, 0:3]); axE = fig.add_subplot(gs[1, 3:6])

strip_panel(axA, "acidic", "D+E content of positions 3–24 (%)")
strip_panel(axB, "gravy", "Hydropathy (mean Kyte–Doolittle)")
strip_panel(axC, "net_charge", "Net charge, positions 3–24")
band_panel(axD, "Met-Pro")
band_panel(axE, "Met-Gly")

for ax, letter in ((axA, "A"), (axB, "B"), (axC, "C"), (axD, "D"), (axE, "E")):
    ax.text(-0.16, 1.06, letter, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="bottom", ha="left", color=INK)

fig.suptitle("What separates stabilised peptides from unstabilised ones, with the "
             "N-terminal residue held fixed", fontsize=10, fontweight="bold",
             color=INK, y=0.995)
fig.text(0.5, 0.945, "A–C  composition of positions 3–24, the part the motif does not "
         "fix · D–E  the same comparison inside matched baseline-PSI bands, which "
         "rules out the stability ceiling as the explanation",
         ha="center", fontsize=6.6, color=MUTED, style="italic")

r, p = spearmanr(frames[("Met-Pro", "not stabilised")].acidic,
                 frames[("Met-Pro", "not stabilised")].base)
fig.text(0.5, 0.012,
         "Stabilised = ΔPSI ≥ 1.0 with baseline PSI ≤ 4.0 (Met-Pro: in any of ZYG11B KO, ZER1 KO, "
         "double KO; Met-Gly: double KO).  Not stabilised = ΔPSI < 0.5 everywhere measured.  "
         "Boxes are median and quartiles, whiskers 1.5 IQR, every peptide drawn.  "
         f"Mann-Whitney throughout.\nAmong unstabilised Met-Pro peptides acidic content rises WITH "
         f"baseline stability (Spearman r = {r:+.2f}, p = {p:.0e}), so the ceiling cannot "
         "manufacture panels A and D — it works against them.  The same composition separates "
         "Met-Gly substrates more strongly than Met-Pro ones, so this is not proline-specific.",
         ha="center", va="bottom", fontsize=5.8, color=MUTED, linespacing=1.6)

plt.subplots_adjust(top=0.90, bottom=0.13)
for ext in ("png", "pdf", "svg"):
    path = os.path.join(FIG, f"FigureZ3_composition_stabilised_vs_not.{ext}")
    fig.savefig(path, format=ext)
    print(f"  saved: {os.path.basename(path)}  ({os.path.getsize(path)/1024:.1f} KB)")
plt.close(fig)

# ------------------------------------------------------------------ readout
print("\n" + "=" * 78)
print("composition: 27 tests per motif class, BH-FDR within each")
print("=" * 78)
for name in tests.motif.unique():
    t = tests[tests.motif == name].sort_values("p")
    print(f"\n{name}: {int(t.n_sub.iloc[0])} stabilised vs {int(t.n_ctrl.iloc[0])} not   "
          f"({int((t.p < 0.05).sum())} of {len(t)} at p<0.05, {0.05*len(t):.1f} expected; "
          f"{int((t.q < 0.05).sum())} survive FDR)")
    print(t.head(7)[["test", "stabilised", "not_stabilised", "p", "q"]].to_string(
        index=False, formatters={"stabilised": "{:.3f}".format,
                                 "not_stabilised": "{:.3f}".format,
                                 "p": "{:.2e}".format, "q": "{:.4f}".format}))

print("\n" + "=" * 78)
print("baseline-matched control")
print("=" * 78)
print(matched.to_string(index=False, formatters={"acidic_sub": "{:.1f}".format,
                                                 "acidic_ctrl": "{:.1f}".format,
                                                 "p": "{:.4f}".format}))
for name in matched.motif.unique():
    m = matched[matched.motif == name]
    print(f"   {name}: combined p = {combine_pvalues(m.p, method='stouffer')[1]:.3g}")
