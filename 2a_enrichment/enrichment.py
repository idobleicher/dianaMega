"""
Is D/E at position 2 of the 2A downstream product (the UBR3 motif) enriched or
depleted relative to what you would expect by chance?

The 2A skip happens between the terminal G and P of the conserved ...PG|P motif, so
every downstream product starts with that proline (P1 = P, always). The residue that
follows it, P2, is what UBR3 reads. The question is whether D/E turns up at P2 more or
less often than background.

Three backgrounds, from 17,517 reviewed (Swiss-Prot) viral proteins:

  BG1  all amino acids, all positions
       The crudest expectation: how common are D and E in viral protein sequence at all?
  BG2  position 2 of viral proteins (residue after the initiator Met)
       Controls for the fact that N-terminal positions have their own composition,
       which is not the same as the bulk average.
  BG3  the residue immediately after ANY proline (P+1)
       The best-matched control: P2 in a 2A product sits directly after a proline, and
       proline constrains what follows it. If D/E were merely rare after prolines in
       general, BG1 and BG2 would overstate the depletion.

Observed foreground, from the combined 2A set (38,180 instances, 2a_combined_UBR3).

THE COUNTING UNIT DECIDES THE ANSWER, so all of them are reported. 84% of instances
come from UniParc, which deposits every sequence variant separately - the same protein
can appear hundreds of times. Counting raw instances therefore counts the same
downstream context over and over, which drags the composition toward whatever the most
duplicated peptides happen to carry. The defensible unit is the DISTINCT DOWNSTREAM
CONTEXT: the residue at P2 is a property of the protein that follows the skip site, not
of the 2A peptide, so each distinct downstream product should count once.

  - all instances                      raw, redundancy included
  - unique source accession            one site per protein
  - unique downstream product (20 aa)  PRIMARY - one distinct downstream context each
  - unique 2A core                     one distinct 2A peptide each
  - viral-lineage subset               raw and deduplicated

Outputs: enrichment_results.csv, enrichment_summary.json,
         enrichment_figure.(png|pdf|svg), 2A_UBR3_enrichment.pptx
"""
import gzip
import json
import pathlib
from collections import Counter

import numpy as np
import pandas as pd
from scipy import stats

HERE = pathlib.Path(__file__).parent
INSTANCES = HERE.parent / "2a_combined_UBR3" / "combined_all_instances.csv"
VIRAL_FASTA = HERE / "viral_swissprot.fasta.gz"
AA = list("ACDEFGHIKLMNPQRSTVWY")
TARGET = ("D", "E")

# ------------------------------------------------------- foreground (2A set)
inst = pd.read_csv(INSTANCES)
inst["lineage"] = inst["lineage"].fillna("")
inst["is_viral"] = inst["lineage"].str.contains("Viruses", case=False)

res = inst[inst["pos2"].notna() & inst["pos2"].isin(AA)].copy()
res_d = res.dropna(subset=["downstream_product"])
vir = res[res["is_viral"]]
vir_d = res_d[res_d["is_viral"]]

# ordered so the figure reads from most-redundant to properly deduplicated
fg = {
    "All 2A instances (raw)": res["pos2"],
    "Unique source accession": res.drop_duplicates("accession")["pos2"],
    "Unique 2A core peptide": res.drop_duplicates("core")["pos2"],
    "Unique downstream context": res_d.drop_duplicates("downstream_product")["pos2"],
    "Viral 2A instances (raw)": vir["pos2"],
    "Viral, unique downstream context": vir_d.drop_duplicates("downstream_product")["pos2"],
}
PRIMARY = "Unique downstream context"

# ------------------------------------------------------- backgrounds (viral)
def read_fasta(path):
    name, buf = None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield "".join(buf)
                name, buf = line, []
            else:
                buf.append(line.strip())
    if name:
        yield "".join(buf)


all_pos = Counter()
pos2 = Counter()
after_p = Counter()
n_prot = n_met = 0
for seq in read_fasta(VIRAL_FASTA):
    if not seq:
        continue
    n_prot += 1
    all_pos.update(c for c in seq if c in AA)
    if seq[0] == "M" and len(seq) > 1 and seq[1] in AA:
        n_met += 1
        pos2[seq[1]] += 1
    for i, c in enumerate(seq[:-1]):
        if c == "P" and seq[i + 1] in AA:
            after_p[seq[i + 1]] += 1

backgrounds = {
    "All amino acids, all positions": all_pos,
    "Position 2 of viral proteins": pos2,
    "Residue after any proline (P+1)": after_p,
}

# ------------------------------------------------------------------ testing
def counts_of(series_or_counter):
    if isinstance(series_or_counter, Counter):
        c = series_or_counter
    else:
        c = Counter(series_or_counter)
    total = sum(c[a] for a in AA)
    hit = sum(c[a] for a in TARGET)
    return hit, total, c


rows = []
for fname, fseries in fg.items():
    f_hit, f_tot, f_c = counts_of(fseries)
    f_pct = 100 * f_hit / f_tot
    lo, hi = stats.binomtest(f_hit, f_tot).proportion_ci(confidence_level=0.95)
    for bname, bc in backgrounds.items():
        b_hit, b_tot, _ = counts_of(bc)
        b_pct = 100 * b_hit / b_tot
        table = [[f_hit, f_tot - f_hit], [b_hit, b_tot - b_hit]]
        odds, p = stats.fisher_exact(table, alternative="two-sided")
        chi2, p_chi, _, _ = stats.chi2_contingency(table)
        rows.append({
            "foreground": fname,
            "foreground_n": f_tot,
            "observed_DE_n": f_hit,
            "observed_DE_pct": round(f_pct, 3),
            "observed_CI95_low_pct": round(100 * lo, 3),
            "observed_CI95_high_pct": round(100 * hi, 3),
            "background": bname,
            "background_n": b_tot,
            "background_DE_pct": round(b_pct, 3),
            "fold_change": round(f_pct / b_pct, 3),
            "depletion_fold": round(b_pct / f_pct, 2) if f_pct else np.nan,
            "odds_ratio": round(odds, 3),
            "fisher_p": p,
            "chi2_p": p_chi,
            "direction": "enriched" if f_pct > b_pct else "depleted",
        })
results = pd.DataFrame(rows)

# per-amino-acid composition table, for the figure
comp = pd.DataFrame({"amino_acid": AA})
for fname, fseries in fg.items():
    _, tot, c = counts_of(fseries)
    comp[fname] = [100 * c[a] / tot for a in AA]
for bname, bc in backgrounds.items():
    _, tot, c = counts_of(bc)
    comp[bname] = [100 * c[a] / tot for a in AA]
comp = comp.round(3)

results.to_csv(HERE / "enrichment_results.csv", index=False)
comp.to_csv(HERE / "enrichment_composition.csv", index=False)

summary = {
    "question": ("Is D/E at position 2 of the 2A downstream product (UBR3 motif) "
                 "enriched or depleted vs background?"),
    "foreground_source": str(INSTANCES.relative_to(HERE.parent)),
    "background_source": "UniProt Swiss-Prot, reviewed viral proteins (taxonomy 10239)",
    "background_proteins": n_prot,
    "background_residues": int(sum(all_pos.values())),
    "background_proteins_starting_with_Met": n_met,
    "foregrounds": {k: int(counts_of(v)[1]) for k, v in fg.items()},
    "headline": {},
    "parameters": {"target_residues": list(TARGET),
                   "test": "Fisher exact (two-sided) and chi-square on 2x2 counts",
                   "ci": "exact binomial 95% CI on the observed proportion"},
}
summary["primary_foreground"] = PRIMARY
summary["redundancy"] = {
    "instances": int(len(res)),
    "unique_accessions": int(res["accession"].nunique()),
    "unique_downstream_contexts": int(res_d["downstream_product"].nunique()),
    "instances_per_unique_downstream_context": round(
        len(res_d) / res_d["downstream_product"].nunique(), 1),
    "pct_instances_from_UniParc": round(100 * (res["db"] == "UniParc").mean(), 1),
}
for unit in (PRIMARY, "All 2A instances (raw)"):
    summary["headline"][unit] = {
        r.background: {"observed_pct": r.observed_DE_pct,
                       "background_pct": r.background_DE_pct,
                       "fold": r.fold_change, "p": r.fisher_p,
                       "direction": r.direction}
        for r in results[results["foreground"] == unit].itertuples()}
(HERE / "enrichment_summary.json").write_text(json.dumps(summary, indent=2))

# ------------------------------------------------------------------ console
print(f"background: {n_prot:,} reviewed viral proteins, "
      f"{sum(all_pos.values()):,} residues ({n_met:,} start with Met)\n")
for fname in fg:
    sub = results[results["foreground"] == fname]
    o = sub.iloc[0]
    print(f"--- {fname}  (n={o.foreground_n:,})  D/E at P2 = {o.observed_DE_pct}% "
          f"[95% CI {o.observed_CI95_low_pct}-{o.observed_CI95_high_pct}]")
    for r in sub.itertuples():
        print(f"      vs {r.background:34} {r.background_DE_pct:6.2f}%   "
              f"fold {r.fold_change:5.2f}  ({r.direction} {r.depletion_fold}x)  "
              f"p={r.fisher_p:.3g}")
    print()
print("D and E separately:")
for unit in ("All 2A instances (raw)", PRIMARY):
    _, tot, c = counts_of(fg[unit])
    print(f"  {unit:34} D {100*c['D']/tot:5.2f}%   E {100*c['E']/tot:5.2f}%")
print(f"\nredundancy: {summary['redundancy']['instances']:,} instances collapse to "
      f"{summary['redundancy']['unique_downstream_contexts']:,} distinct downstream "
      f"contexts ({summary['redundancy']['instances_per_unique_downstream_context']}x); "
      f"{summary['redundancy']['pct_instances_from_UniParc']}% of instances are UniParc")
