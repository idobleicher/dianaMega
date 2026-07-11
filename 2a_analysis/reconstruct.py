#!/usr/bin/env python3
"""Reconstruct downstream products (post-skip) for every 2A row and filter the
UBR3 motif: downstream position 1 = P (skip proline), position 2 = D or E.
Outputs full table + FASTA and the sorted motif-hit subset."""
import json, os, csv
from collections import Counter

ROOT = os.path.dirname(os.path.abspath(__file__))
recs = json.load(open(os.path.join(ROOT, "parsed_rows.json")))
cache = json.load(open(os.path.join(ROOT, "fasta_cache.json"))) if os.path.exists(os.path.join(ROOT,"fasta_cache.json")) else {}

DOWN_LEN = 20  # residues of downstream product to report (including skip-P)

def locate_skipP(src, seq):
    """Locate the skip proline (downstream position 1) in the SOURCE protein.
    Skip site is the G|P of the conserved ...PG|P motif. Returns
    (skipP_index, method) where src[skipP_index] is the terminal/skip P."""
    if not src:
        return None, "no_source"
    # --- anchor the 2A region in the source ---
    idx = src.find(seq)
    if idx != -1:
        a_start, a_end, method = idx, idx + len(seq), "exact"
    else:
        a_start = a_end = -1; method = None
        for k in (14, 13, 12, 11, 10, 9, 8, 7):
            if k >= len(seq):
                continue
            j = src.rfind(seq[-k:])
            if j != -1:
                a_start, a_end, method = j, j + k, f"suffix{k}"
                break
        if method is None:
            return None, "not_found"
    # --- find the skip proline via the ...PGP motif near the anchor end ---
    # window covers the anchor tail plus up to 3 residues past it (handles
    # shown sequences that stop at ...NPG, one residue short of the skip P)
    win_start = max(a_start, a_end - 6)
    window = src[win_start:a_end + 3]
    p = window.rfind("PGP")
    if p != -1:
        return win_start + p + 2, method            # index of the final P in PGP
    # fallbacks when the PGP motif is absent from the shown region
    if a_end - 1 < len(src) and src[a_end - 1] == "P":
        return a_end - 1, method + "_endP"
    if a_end < len(src) and src[a_end] == "P":
        return a_end, method + "_nextP"
    return None, method + "_noPGP"

rows_out = []
stat = Counter()
pos2_counter = Counter()
for r in recs:
    src = cache.get(r["acc"], "") if r["acc"] else ""
    skipP, method = locate_skipP(src, r["seq"])
    pos2 = ""; downstream = ""; found = False
    if skipP is None:
        stat[method] += 1
        downstream = "P"               # skip-P is guaranteed even if source is missing
    else:
        stat[method] += 1
        if src[skipP] != "P":
            stat["skipP_notP"] += 1    # should not happen; diagnostic
        found = True
        downstream = src[skipP:skipP + DOWN_LEN]   # begins AT the skip proline
        if skipP + 1 < len(src):
            pos2 = src[skipP + 1]
            pos2_counter[pos2] += 1
        else:
            stat["at_protein_Cterm"] += 1          # 2A at very C-terminus
    rows_out.append({
        "name": r["name"], "aid": r["aid"], "db": r["db"], "class": r["cls"],
        "origin": r["origin"], "twoA_seq": r["seq"], "found": found,
        "pos2": pos2, "downstream_product": downstream, "lineage": r["lineage"],
    })

# ---- write full table ----
cols = ["name","aid","db","class","origin","twoA_seq","found","pos2","downstream_product","lineage"]
with open(os.path.join(ROOT,"downstream_all.csv"),"w",newline="",encoding="utf-8") as f:
    w=csv.DictWriter(f,fieldnames=cols); w.writeheader()
    for x in rows_out: w.writerow(x)

# ---- motif filter: pos2 in {D,E} ----
hits = [x for x in rows_out if x["pos2"] in ("D","E")]
hits.sort(key=lambda x: (x["downstream_product"], x["aid"]))  # sort by amino-acid sequence

with open(os.path.join(ROOT,"motif_PDE_hits.csv"),"w",newline="",encoding="utf-8") as f:
    w=csv.DictWriter(f,fieldnames=cols); w.writeheader()
    for x in hits: w.writerow(x)

with open(os.path.join(ROOT,"motif_PDE_hits.fasta"),"w",encoding="utf-8") as f:
    for x in hits:
        f.write(f">{x['aid']} | {x['name']} | class{x['class']} | {x['origin']} | pos2={x['pos2']}\n{x['downstream_product']}\n")

# ---- full downstream FASTA (only entries with a resolved downstream residue) ----
with open(os.path.join(ROOT,"downstream_all.fasta"),"w",encoding="utf-8") as f:
    for x in rows_out:
        if x["found"] and x["pos2"]:
            f.write(f">{x['aid']} | {x['name']} | class{x['class']} | {x['origin']} | pos2={x['pos2']}\n{x['downstream_product']}\n")

# ---- summary ----
resolved = [x for x in rows_out if x["found"] and x["pos2"]]
print("TOTAL rows:", len(rows_out))
print("status breakdown:", dict(stat))
print("resolved downstream (pos2 known):", len(resolved))
print("pos2 distribution (top 25):", pos2_counter.most_common(25))
nD = sum(1 for x in hits if x["pos2"]=="D"); nE = sum(1 for x in hits if x["pos2"]=="E")
print(f"MOTIF HITS (P then D/E): {len(hits)}  (PD={nD}, PE={nE})")
byclass = Counter((x["class"],x["origin"]) for x in hits)
print("motif hits by (class,origin):", dict(byclass))
print("\nFirst 15 motif hits (sorted):")
for x in hits[:15]:
    print(f"  {x['downstream_product']:22s} {x['db']:9s} {x['aid'][:38]:38s} {x['origin']}")

# integrity check: every motif hit MUST be P at pos1 and D/E at pos2
bad = [x for x in hits if not (x["downstream_product"][:1] == "P" and x["downstream_product"][1:2] in ("D","E"))]
print("\nINTEGRITY: malformed motif hits (should be 0):", len(bad))
if bad:
    for x in bad[:5]: print("   BAD:", x["downstream_product"], x["aid"])
