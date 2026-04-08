"""
Search NCBI for correct RPL22 and RPL22L1 orthologs.
Filter by expected size range to avoid multi-domain fusion proteins.
"""
import time
from io import StringIO
from Bio import SeqIO, Entrez

Entrez.email = "rpl22.alignment@example.com"
Entrez.tool  = "RPL22_ortholog_search"

ORGANISMS = {
    "H_sapiens": "Homo sapiens",
    "M_musculus": "Mus musculus",
    "X_laevis": "Xenopus laevis",
    "P_marinus": "Petromyzon marinus",
    "C_intestinalis": "Ciona intestinalis",
    "B_floridae": "Branchiostoma floridae",
    "S_purpuratus": "Strongylocentrotus purpuratus",
    "O_sinensis": "Octopus sinensis",
    "D_melanogaster": "Drosophila melanogaster",
    "A_mellifera": "Apis mellifera",
    "P_tepidariorum": "Parasteatoda tepidariorum",
    "H_vulgaris": "Hydra vulgaris",
    "N_vectensis": "Nematostella vectensis",
}


def search_protein(gene_name, organism, min_len=80, max_len=250):
    """Search for a protein by gene name and organism, filtering by size."""
    queries = [
        f'"{gene_name}"[Gene Name] AND "{organism}"[Organism] AND refseq[filter]',
        f'"{gene_name}"[Gene Name] AND "{organism}"[Organism]',
        f'"ribosomal protein L22" AND "{organism}"[Organism] AND refseq[filter]',
        f'"60S ribosomal protein L22" AND "{organism}"[Organism]',
    ]

    if gene_name == "RPL22L1":
        queries = [
            f'"RPL22L1"[Gene Name] AND "{organism}"[Organism] AND refseq[filter]',
            f'"RPL22L1"[Gene Name] AND "{organism}"[Organism]',
            f'"RPL22-like" AND "{organism}"[Organism] AND refseq[filter]',
            f'"ribosomal protein L22-like" AND "{organism}"[Organism]',
            f'"eL22-like" AND "{organism}"[Organism]',
        ]

    for query in queries:
        try:
            handle = Entrez.esearch(db="protein", term=query, retmax=20)
            result = Entrez.read(handle)
            handle.close()

            if not result["IdList"]:
                continue

            ids = ",".join(result["IdList"][:20])
            handle = Entrez.efetch(db="protein", id=ids,
                                   rettype="fasta", retmode="text")
            text = handle.read()
            handle.close()

            candidates = []
            for rec in SeqIO.parse(StringIO(text), "fasta"):
                seq_len = len(rec.seq)
                if min_len <= seq_len <= max_len:
                    candidates.append((rec, seq_len))

            if candidates:
                candidates.sort(key=lambda x: abs(x[1] - 128))
                return candidates[0][0]

            time.sleep(0.5)
        except Exception as e:
            print(f"    Error: {e}")
            time.sleep(2)

    return None


def main():
    for gene in ["RPL22", "RPL22L1"]:
        expected_size = 128 if gene == "RPL22" else 122
        min_len = 80
        max_len = 200
        print(f"\n{'='*60}")
        print(f"Searching for {gene} orthologs (expected ~{expected_size} aa)")
        print(f"{'='*60}")

        results = {}
        for sp_id, organism in ORGANISMS.items():
            print(f"\n  {sp_id} ({organism}):")
            rec = search_protein(gene, organism, min_len, max_len)
            if rec:
                acc = rec.id.split("|")[-1] if "|" in rec.id else rec.id
                print(f"    FOUND: {acc} ({len(rec.seq)} aa)")
                print(f"    Desc: {rec.description[:80]}")
                results[sp_id] = acc
            else:
                print(f"    NOT FOUND")
            time.sleep(0.5)

        print(f"\n\n{gene}_ORTHOLOGS = {{")
        for sp_id, acc in results.items():
            print(f'    "{sp_id}": "{acc}",')
        print("}")


if __name__ == "__main__":
    main()
