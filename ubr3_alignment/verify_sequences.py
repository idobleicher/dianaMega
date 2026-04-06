"""Verify all UBR3 sequences are correct full-length proteins with UBR box."""
from Bio import AlignIO, Entrez, SeqIO
from io import StringIO
import time

Entrez.email = "ubr3.alignment@example.com"

alignment = AlignIO.read("ubr3_aligned.fasta", "fasta")

# Find human index and map columns to human positions
human_idx = 0
for i, rec in enumerate(alignment):
    if "sapiens" in rec.id:
        human_idx = i
        break

col_to_human = {}
res = 0
for c in range(alignment.get_alignment_length()):
    if alignment[human_idx, c] != "-":
        res += 1
        col_to_human[c] = res

human_to_col = {v: k for k, v in col_to_human.items()}

# UBR box = human 118-189
ubr_col_s = human_to_col[118]
ubr_col_e = human_to_col[189]
# RING = human 1306-1364
ring_col_s = human_to_col[1306]
ring_col_e = human_to_col[1364]

print(f"Alignment: {len(alignment)} seqs x {alignment.get_alignment_length()} cols")
print(f"UBR box cols: {ubr_col_s}-{ubr_col_e}")
print(f"RING cols:    {ring_col_s}-{ring_col_e}")
print()

print(f"{'Species':<25s} {'Length':>6s}  {'UBR gaps':>9s}  {'RING gaps':>10s}  {'UBR box?':>8s}  {'RING?':>6s}  Status")
print("-" * 100)

problems = []
for seq_i in range(len(alignment)):
    rec = alignment[seq_i]
    seq_len = sum(1 for c in range(alignment.get_alignment_length())
                  if rec.seq[c] != "-")

    ubr_gaps = sum(1 for c in range(ubr_col_s, ubr_col_e + 1)
                   if rec.seq[c] == "-")
    ubr_total = ubr_col_e - ubr_col_s + 1
    ubr_pct = (ubr_total - ubr_gaps) / ubr_total * 100

    ring_gaps = sum(1 for c in range(ring_col_s, ring_col_e + 1)
                    if rec.seq[c] == "-")
    ring_total = ring_col_e - ring_col_s + 1
    ring_pct = (ring_total - ring_gaps) / ring_total * 100

    has_ubr = ubr_pct > 50
    has_ring = ring_pct > 50

    status = "OK" if (has_ubr and has_ring) else "PROBLEM"
    if not has_ubr:
        status += " - MISSING UBR BOX"
        problems.append((rec.id, "UBR box", ubr_pct))
    if not has_ring:
        status += " - MISSING RING"
        problems.append((rec.id, "RING", ring_pct))

    print(f"{rec.id:<25s} {seq_len:>6d}  {ubr_gaps:>4d}/{ubr_total:<4d}  "
          f"{ring_gaps:>4d}/{ring_total:<4d}   "
          f"{'YES' if has_ubr else 'NO':>6s}   "
          f"{'YES' if has_ring else 'NO':>5s}  {status}")

if problems:
    print(f"\n{'='*60}")
    print(f"PROBLEMS FOUND: {len(problems)}")
    for sp, domain, pct in problems:
        print(f"  {sp}: {domain} only {pct:.0f}% present")
    print(f"{'='*60}")
else:
    print(f"\nAll sequences have both UBR box and RING domain.")

# Also check key conserved cysteines in UBR box
print(f"\n=== Key UBR box residues (human C120, C133, C142, C148, C151) ===")
key_positions = [120, 133, 142, 148, 151]
for hpos in key_positions:
    col = human_to_col[hpos]
    residues = []
    for seq_i in range(len(alignment)):
        aa = alignment[seq_i, col]
        residues.append(aa)
    all_same = len(set(r for r in residues if r != "-")) <= 1
    print(f"  Human pos {hpos} (col {col}): {' '.join(residues)}  "
          f"{'CONSERVED' if all_same else 'VARIABLE'}")
