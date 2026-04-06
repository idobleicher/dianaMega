"""Debug: check why D. melanogaster UBR box is missing from the alignment."""
from Bio import AlignIO

alignment = AlignIO.read("ubr3_aligned.fasta", "fasta")

fly_idx = human_idx = None
for i, rec in enumerate(alignment):
    if "melanogaster" in rec.id:
        fly_idx = i
    if "sapiens" in rec.id:
        human_idx = i

print(f"Human index: {human_idx} ({alignment[human_idx].id})")
print(f"Fly index:   {fly_idx} ({alignment[fly_idx].id})")
print(f"Alignment length: {alignment.get_alignment_length()}")

# Count actual residues (non-gap)
fly_len = sum(1 for c in range(alignment.get_alignment_length())
              if alignment[fly_idx, c] != "-")
human_len = sum(1 for c in range(alignment.get_alignment_length())
                if alignment[human_idx, c] != "-")
print(f"Human length: {human_len} aa")
print(f"Fly length:   {fly_len} aa")

# Map human pos -> alignment col
human_to_col = {}
res = 0
for c in range(alignment.get_alignment_length()):
    if alignment[human_idx, c] != "-":
        res += 1
        human_to_col[res] = c

# UBR box = human 118-189
col_s = human_to_col[118]
col_e = human_to_col[189]
print(f"\n=== UBR Box (Human 118-189) -> Alignment cols {col_s}-{col_e} ===")

human_seq = ""
fly_seq = ""
for c in range(col_s, col_e + 1):
    human_seq += alignment[human_idx, c]
    fly_seq += alignment[fly_idx, c]

print(f"Human: {human_seq}")
print(f"Fly:   {fly_seq}")
gap_count = fly_seq.count("-")
print(f"Fly gaps in this region: {gap_count} / {len(fly_seq)}")

# Where does the fly sequence START in the alignment?
first_fly_col = None
for c in range(alignment.get_alignment_length()):
    if alignment[fly_idx, c] != "-":
        first_fly_col = c
        break
print(f"\nFly first residue at alignment col: {first_fly_col}")

# What human position does that correspond to?
res = 0
for c in range(alignment.get_alignment_length()):
    if alignment[human_idx, c] != "-":
        res += 1
    if c == first_fly_col:
        print(f"At fly start col {first_fly_col}, human position = {res}")
        break

# Show fly residues 1-200 and their alignment positions
print(f"\n=== Fly residues 1-150 mapped to human ===")
fly_res = 0
for c in range(alignment.get_alignment_length()):
    if alignment[fly_idx, c] != "-":
        fly_res += 1
        if fly_res <= 150:
            h_res = None
            # find closest human position
            hr = 0
            for cc in range(c + 1):
                if alignment[human_idx, cc] != "-":
                    hr += 1
            h_res = hr if alignment[human_idx, c] != "-" else f"~{hr}"
            if fly_res % 10 == 1 or fly_res <= 5:
                print(f"  Fly {fly_res:4d}  ->  col {c:5d}  "
                      f"fly={alignment[fly_idx, c]}  "
                      f"human={alignment[human_idx, c]}  hpos={h_res}")

# Show the raw fly FASTA header to check accession
print(f"\n=== Fly record details ===")
print(f"ID: {alignment[fly_idx].id}")
print(f"Name: {alignment[fly_idx].name}")
print(f"Description: {alignment[fly_idx].description}")
