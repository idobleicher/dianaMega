# Stability / half-life of the UBR3-WT best hits

**Question.** What are the half-lives of the best hits from the UBR3-WT co-IP MS?

**Source.** `Supplementary Table 1 (3).xlsx`, Yi et al. (2023), "Degradomics" sheet —
fractionated TMTpro-MS3 analysis of WT, ATG7-KO and FIP200(RB1CC1)-KO extracts.

## Read this first: the source table has no half-life column

It is a **time course**, not a half-life table. Each protein has the fraction
remaining at 0 h (two replicates), 5 h, 10 h and 15 h, in WT / ATG7-KO / RB1CC1-KO
extracts, untreated (UT) or Torin-treated (Tor), as a ratio to the 0 h WT untreated
baseline. Half-lives here are **derived** by log-linear regression:

    ln(ratio_t) = -k·t     t½ = ln2 / k

The assay runs 15 h and the median protein still retains 82% at the end, so **39 of
the 45 quantified hits never reach 50% inside the measured window** — their half-lives
are extrapolations and should be treated as ordering information, not as measurements.
`frac_remaining_15h` (mean of the two 15 h replicates) is model-free and is the
column to rank on. `beyond_window` and `no_decay` flag the extrapolated and the
non-decaying proteins respectively.

Also note this is degradation **in extract**, and the study is an autophagy screen
(ATG7 / FIP200 knockouts). It is not a cellular half-life measured by pulse-chase, and
it is not specific to the ubiquitin–proteasome system that UBR3 works through.

## How to reproduce

```
python wt_besthits_halflife.py       # -> CSVs + summary JSON
python plot_wt_besthits_halflife.py  # -> WT_besthits_halflife.png/pdf/svg
```

Matching is by UniProt accession (45 hits), then by any accession in the MaxQuant
majority-protein group (2 more: H3C1→HIST2H3A, H2BC13→HIST1H2BN — histone peptides
are shared, so the two studies assign them to different paralogs), then by gene
symbol. **5 of the 50 WT best hits have no usable data**: H3-7, H2AC6, CLECL1P and
H3C1 are absent, and TNRC6C is present with 0 quantified peptides.

## Result

**As a group, the WT best hits are not unusually short-lived.** Median 15 h remaining
0.829 for the hits vs 0.815 for the 8097-protein background (Mann–Whitney p = 0.80);
decay rate p = 0.78; derived median half-life 41 h in both. Whatever UBR3-WT is
pulling down, it is not enriched for fast-turnover proteins in this assay.

**The signal is inside one module.** Split by the pathway modules from
`WT_mutual_pathways_README.md`:

| Module | n | median 15 h remaining | median derived t½ | p vs background |
|---|---|---|---|---|
| M1 mRNA silencing & RNP granules | 12 | 0.844 | 65 h | 0.48 |
| **M2 ATP-dependent chromatin remodeling** | **9** | **0.576** | **19 h** | **0.069** |
| M3 Nucleosome / histone core | 4 | 0.846 | 64 h | 0.99 |
| M4 BAF–LAP2β nuclear envelope | 2 | 1.019 | 68 h | 0.19 |

M2 — SMARCA5, SMARCA1, CHD7, BLM, POLE, TET1, UFD1, MACROH2A1, H2BC13 — turns over
roughly twice as fast as the proteome. **This is a trend, not a significant result**:
p = 0.069 raw and 0.28 after Bonferroni across the four modules, with n = 9. It is
worth a targeted follow-up (cycloheximide chase on SMARCA5 / CHD7 / TET1 ± UBR3), not
a claim in a figure legend.

**UBR3 itself is among the least stable proteins in its own hit list** — 33% remaining
at 15 h, derived t½ 9.7 h, fit R² = 0.87. Consistent with E3 ligase auto-degradation.

**Least stable hits** (15 h remaining): CCHCR1 0.24, ISL1 0.25, UBR3 0.33, TET1 0.40,
H2AC18 0.45, POLE 0.47, UFD1 0.52, CHD7 0.54, BLM 0.58.
**Most stable** (no measurable decay): H4C1, BANF1, TDRD3, EIF3G, FTSJ3, OTUD4, GCFC2.

Two negative results worth recording: all 45 quantified hits carry a *stabilizing*
N-degron in `hits_table.csv`, so there is no N-degron/stability split to test; and
WT-preferring (RING-dependent) hits are no less stable than the background
(0.830, p = 0.69), so RING-dependent binding does not predict fast turnover here.

## Which half-lives are worth quoting

`wt_besthits_halflife_significant.py` applies two tests, because either alone misleads:

1. **The decay is real** — the slope of ln(ratio) on time is significantly below zero,
   Benjamini-Hochberg corrected across the 45 quantified hits. This alone is weak:
   **31 of 45 hits pass**, because with 8 timepoints the fit is precise enough to call
   even a gentle slope real. Statistically significant decay ≠ a short-lived protein.
2. **It is fast relative to the proteome** — decay rate k above the 80th percentile of
   the 8097-protein background (k > 0.0363/h), i.e. the top 20% of proteome turnover.
   This is the test that discriminates. Threshold is `FAST_PCTL` in the script.

**Nine hits pass both** (was 4 at the stricter 90th-percentile setting):

| Gene | t½ (h) | 15 h remaining | q (BH) | R² | percentile | evidence | module |
|---|---|---|---|---|---|---|---|
| CCHCR1 | 7.6 | 0.239 | 0.003 | 0.85 | 97th | measured | — |
| ISL1 | 7.7 | 0.250 | 5.6e-05 | 0.97 | 97th | measured | — |
| UBR3 | 9.7 | 0.333 | 0.002 | 0.87 | 94th | measured | — |
| TET1 | 11.6 | 0.399 | 4.4e-04 | 0.93 | 91st | measured | M2 |
| H2AC18 | 13.2 | 0.452 | 1.2e-04 | 0.97 | 88th | measured | — |
| POLE | 14.5 | 0.474 | 9.2e-04 | 0.91 | 86th | measured | M2 |
| UFD1 | 16.2 | 0.521 | 2.3e-05 | 0.99 | 84th | extrapolated | M2 |
| CHD7 | 17.0 | 0.541 | 2.9e-05 | 0.98 | 83rd | extrapolated | M2 |
| BLM | 19.0 | 0.576 | 2.1e-05 | 0.99 | 80th | extrapolated | M2 |

"measured" = the curve crosses 50% inside the 15 h assay, so t½ is read from data.
"extrapolated" = decay is significant but t½ falls beyond the assay; the three here
have R² ≥ 0.98 and land just past 15 h, so they are the defensible extrapolations.

Five of the nine are M2 chromatin-remodeling members (TET1, POLE, UFD1, CHD7, BLM),
which is the per-protein counterpart of the group-level M2 trend above. SMARCA5 and
SMARCA1, the two proteins that anchor M2 as a complex, do **not** pass (48th and 66th
percentile) — so the module's instability is not uniform across its members.

## Cross-check against published turnover data

`wt_besthits_halflife_crosscheck.py` compares the Yi-derived numbers with two
independent datasets:

* **Li et al. 2021** (Mol Cell, Gygi lab) — cycloheximide chase + TMTpro in *living
  cells*, 0–8 h, in U2OS, **HEK293T** (the co-IP background), HCT116 and RPE1.
  Table S2 gives the authors' own fitted half-lives. It lists **only** proteins that
  lost >15% at 8 h with R²>0.8, so absence is "not degraded **or** not quantified" —
  Table S1 would be needed to separate those. Never read absence as stability.
* **Mathieson et al. 2018** (Nat Commun) — dynamic SILAC in 5 non-dividing primary
  cell types; we use the median of the 4 human ones. Downloaded automatically to
  `external_Mathieson2018_halflives.xlsx`.

**The in-extract assay is validated.** Across every protein shared with the Li
HEK293T chase, fraction remaining at 15 h in extract correlates with cellular
half-life at **Spearman ρ = +0.64 (n = 1112, p = 2e-130)**; against Mathieson primary
cells, ρ = +0.44 (n = 5249). Both have the expected sign. The Li range is truncated to
degraded proteins, which attenuates rather than inflates ρ. So the Yi-derived ordering
is trustworthy even though the absolute half-lives are extrapolated.

**13 of the 50 hits have an independent half-life in at least one Li cell line**, 6 in
HEK293T, and 31 have a Mathieson value. Concordance per hit:

| Gene | Yi (extract) | Li HEK293T | Li other lines | Mathieson | verdict |
|---|---|---|---|---|---|
| CCHCR1 | 7.6 | **3.65** (short-lived) | 3.63 U2OS, 1.67 HCT116 | — | confirmed, very short |
| N4BP2L2 | 21.8 (R²=0.36) | **5.25** (short-lived) | 8.78 / 8.80 / 7.47 — all 4 lines | — | **confirmed short; Yi missed it** |
| TET1 | 11.6 | 10.25 | — | — | confirmed |
| UBR3 | 9.7 | not listed | 11.44 HCT116 | 21.9 | confirmed short |
| CHD7 | 17.0 | 23.89 | 20.19 HCT116, 13.07 RPE1 | — | confirmed |
| BLM | 19.0 | not listed | 19.81 HCT116 | — | confirmed |
| SMARCA1 | 31.5 | not listed | 24.16 HCT116, 26.29 RPE1 | 17.9 | short in cells |
| NCOA1 | 33.3 | 22.15 | 24.33 / 26.82 / 19.02 — all 4 | 25.6 | consistent ~20–27 h |
| RAI1 | 81.2 | not listed | 15.91 HCT116, 17.92 RPE1 | — | **Yi overestimated badly** |
| CSTF2T | 95.3 | 32.36 | — | 50.7 | Yi overestimated |

Three things this changes:

1. **N4BP2L2 is the clearest short-lived hit in the list** — 5–9 h in all four Li cell
   lines — and the Yi-only analysis ranked it 21.8 h on a bad fit (R² = 0.36) and
   dropped it from the significant set. It belongs at the top.
2. **The M2 instability signal is corroborated.** TET1, CHD7, BLM, UFD1 and SMARCA1 all
   appear in Li degraded lists. SMARCA1 in particular looked unremarkable in Yi (31.5 h,
   66th percentile) but is 18–26 h in three independent systems.
3. **SMARCA5 remains long-lived** (Yi 55.5 h, Mathieson 63.4 h, absent from all Li
   lists), so M2's instability really is confined to the DNA-associated members rather
   than the remodeler ATPases as a class.

## Outputs

| File | Contents |
|---|---|
| `WT_besthits_halflife.csv` | per hit: 15 h remaining, decay k, derived t½, fit R², flags, Torin / ATG7-KO / RB1CC1-KO values, module, N-degron class |
| `WT_besthits_halflife_summary.json` | group comparisons, per-module tests, parameters |
| `WT_besthits_halflife_curves.csv` | raw 0/5/10/15 h ratios for the 45 quantified hits |
| `WT_besthits_halflife_background.csv` | same columns for the 8097 background proteins |
| `WT_besthits_halflife.png/pdf/svg` | 3-panel figure (ECDF vs background; per-hit ranking; measured decay curves) |
| `WT_besthits_halflife_significant.csv` | the 9 hits passing both significance tests |
| `WT_besthits_halflife_slopetests.csv` | per-protein slope test, q-value, percentile — including everything that failed |
| `WT_besthits_master_list.xlsx` | all of the above joined per hit; sheets "Half-life" and "Half-life significant" |
