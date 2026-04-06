"""Search for UBR3 orthologs across diverse animal phyla."""
from Bio import Entrez, SeqIO
from io import StringIO
import re, time

Entrez.email = "ubr3.alignment@example.com"

queries = [
    ("C. elegans",  '"UBR" "Caenorhabditis elegans" refseq', 1200),
    ("Sea urchin",  '"ubr3" "Strongylocentrotus" refseq', 1200),
    ("Amphioxus",   '"ubr3" "Branchiostoma" refseq', 1200),
    ("Lamprey",     '"ubr3" "Petromyzon" refseq', 1200),
    ("Hydra",       '"ubr3" "Hydra vulgaris" refseq', 800),
    ("Sea anemone", '"ubr3" "Nematostella" refseq', 800),
    ("Octopus",     '"ubr3" "Octopus" refseq', 1200),
    ("Oyster",      '"ubr3" "Crassostrea" refseq', 1200),
    ("Spider",      '"ubr3" "Parasteatoda" refseq', 1200),
    ("Honeybee",    '"ubr3" "Apis mellifera" refseq', 1200),
    ("Mosquito",    '"ubr3" "Anopheles" refseq', 1200),
    ("Turtle",      '"ubr3" "Chrysemys" refseq', 1400),
    ("Lizard",      '"ubr3" "Anolis" refseq', 1400),
    ("Xenopus",     '"ubr3" "Xenopus laevis" refseq', 1400),
    ("Ciona",       '"ubr3" "Ciona" refseq', 800),
    ("Sponge",      '"ubr" "Amphimedon" refseq', 800),
    ("Planarian",   '"ubr" "Schmidtea" refseq', 800),
]

for name, query, min_len in queries:
    try:
        h = Entrez.esearch(db="protein", term=query, retmax=10)
        res = Entrez.read(h); h.close()
        ids = res.get("IdList", [])
        if not ids:
            print(f"{name:15s}  -- no results --")
            time.sleep(0.4)
            continue
        h2 = Entrez.efetch(db="protein", id=",".join(ids[:5]),
                           rettype="fasta", retmode="text")
        text = h2.read(); h2.close()
        best = None
        for rec in SeqIO.parse(StringIO(text), "fasta"):
            if len(rec.seq) >= min_len:
                if best is None or len(rec.seq) > len(best.seq):
                    best = rec
        if best:
            org = re.search(r'\[(.+?)\]', best.description)
            org_name = org.group(1) if org else "?"
            print(f"{name:15s} {len(best.seq):5d} aa  {best.id:25s}  [{org_name}]")
        else:
            # show what we found even if short
            for rec in SeqIO.parse(StringIO(text), "fasta"):
                org = re.search(r'\[(.+?)\]', rec.description)
                org_name = org.group(1) if org else "?"
                print(f"{name:15s} {len(rec.seq):5d} aa  {rec.id:25s}  [{org_name}]  (short)")
                break
        time.sleep(0.5)
    except Exception as e:
        print(f"{name:15s}  ERROR: {e}")
        time.sleep(1)
