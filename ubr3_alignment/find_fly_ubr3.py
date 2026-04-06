"""Find the full-length D. melanogaster UBR3 with UBR box."""
from Bio import Entrez, SeqIO
from io import StringIO

Entrez.email = "ubr3.alignment@example.com"

# Search for all Drosophila UBR3 RefSeq isoforms
h = Entrez.esearch(db="protein",
                   term='"Ubr3" "Drosophila melanogaster" refseq',
                   retmax=20)
res = Entrez.read(h); h.close()
ids = res.get("IdList", [])
print(f"Found {len(ids)} RefSeq results for Dmel Ubr3")

if ids:
    h2 = Entrez.efetch(db="protein", id=",".join(ids),
                       rettype="fasta", retmode="text")
    text = h2.read(); h2.close()
    for rec in SeqIO.parse(StringIO(text), "fasta"):
        print(f"  {rec.id:25s} {len(rec.seq):5d} aa  {rec.description[:90]}")

# Also search by gene name
print("\n--- Searching by gene ID ---")
h = Entrez.esearch(db="protein",
                   term='"Drosophila melanogaster"[Organism] AND Ubr3[Gene Name]',
                   retmax=20)
res = Entrez.read(h); h.close()
ids2 = res.get("IdList", [])
print(f"Found {len(ids2)} results by gene name")

if ids2:
    h2 = Entrez.efetch(db="protein", id=",".join(ids2),
                       rettype="fasta", retmode="text")
    text = h2.read(); h2.close()
    for rec in SeqIO.parse(StringIO(text), "fasta"):
        print(f"  {rec.id:25s} {len(rec.seq):5d} aa  {rec.description[:90]}")
