"""
Broader search for missing RPL22 and RPL22L1 orthologs.
"""
import time
from io import StringIO
from Bio import SeqIO, Entrez

Entrez.email = "rpl22.alignment@example.com"
Entrez.tool  = "RPL22_ortholog_search"


def search_and_print(query, label, min_len=80, max_len=200):
    try:
        handle = Entrez.esearch(db="protein", term=query, retmax=10)
        result = Entrez.read(handle)
        handle.close()
        if not result["IdList"]:
            print(f"  {label}: no results")
            return
        ids = ",".join(result["IdList"][:10])
        handle = Entrez.efetch(db="protein", id=ids, rettype="fasta", retmode="text")
        text = handle.read()
        handle.close()
        for rec in SeqIO.parse(StringIO(text), "fasta"):
            sl = len(rec.seq)
            marker = " <-- GOOD" if min_len <= sl <= max_len else ""
            print(f"  {label}: {rec.id} ({sl} aa) {rec.description[:60]}{marker}")
        time.sleep(0.5)
    except Exception as e:
        print(f"  {label}: Error - {e}")
        time.sleep(2)


def main():
    print("=== Drosophila RPL22 (need RefSeq, not PDB) ===")
    search_and_print(
        '"RPL22"[Gene Name] AND "Drosophila melanogaster"[Organism]',
        "Drosophila RPL22")
    search_and_print(
        '"RpL22"[Gene Name] AND "Drosophila melanogaster"[Organism]',
        "Drosophila RpL22")
    search_and_print(
        '"60S ribosomal protein L22" AND "Drosophila melanogaster"[Organism] AND refseq[filter]',
        "Drosophila L22 refseq")

    print("\n=== Missing RPL22 species ===")
    for org, name in [
        ("Ciona intestinalis", "C_intestinalis"),
        ("Branchiostoma floridae", "B_floridae"),
        ("Parasteatoda tepidariorum", "P_tepidariorum"),
    ]:
        search_and_print(
            f'"60S ribosomal protein L22" AND "{org}"[Organism]',
            f"{name} L22")
        search_and_print(
            f'"ribosomal protein L22" AND "{org}"[Organism]',
            f"{name} RPL22")
        search_and_print(
            f'"eL22" AND "{org}"[Organism]',
            f"{name} eL22")

    print("\n=== RPL22L1 broader search ===")
    for org, name in [
        ("Petromyzon marinus", "P_marinus"),
        ("Danio rerio", "D_rerio"),
        ("Gallus gallus", "G_gallus"),
        ("Branchiostoma floridae", "B_floridae"),
        ("Strongylocentrotus purpuratus", "S_purpuratus"),
        ("Drosophila melanogaster", "D_melanogaster"),
        ("Octopus sinensis", "O_sinensis"),
    ]:
        search_and_print(
            f'"RPL22L1"[Gene Name] AND "{org}"[Organism]',
            f"{name} RPL22L1")
        search_and_print(
            f'"L22-like" AND "{org}"[Organism] AND refseq[filter]',
            f"{name} L22-like")
        search_and_print(
            f'"eL22-like" AND "{org}"[Organism]',
            f"{name} eL22-like")


if __name__ == "__main__":
    main()
