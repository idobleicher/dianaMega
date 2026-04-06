"""
Rebuild UBR3 alignment with a broad phylogenetic panel,
then regenerate improved publication figures.
"""
import os, sys, time, re
from io import StringIO
from collections import Counter

from Bio import SeqIO, AlignIO, Entrez
from Bio.SeqRecord import SeqRecord
from urllib import request, parse

import numpy as np

Entrez.email = "ubr3.alignment@example.com"
Entrez.tool  = "UBR3_ortholog_alignment"

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

CURATED_ORTHOLOGS = {
    # Mammals
    "H_sapiens":         "NP_742067.3",        # Human
    "M_musculus":        "NP_808451.2",         # Mouse
    # Amphibian
    "X_laevis":          "XP_041433091.1",      # Frog
    # Fish (jawless)
    "P_marinus":         "XP_075927777.1",      # Lamprey
    # Tunicate
    "C_intestinalis":    "XP_078489238.1",      # Sea squirt
    # Cephalochordate
    "B_floridae":        "XP_078702957.1",      # Amphioxus
    # Echinoderm
    "S_purpuratus":      "XP_030846700.1",      # Sea urchin
    # Mollusc
    "O_sinensis":        "XP_036362376.1",      # Octopus
    # Insects
    "D_melanogaster":    "NP_572428.3",         # Fruit fly (full-length isoform A)
    "A_mellifera":       "XP_016766865.1",      # Honeybee
    # Arachnid
    "P_tepidariorum":    "XP_015917471.2",      # Spider
    # Cnidaria
    "H_vulgaris":        "XP_065665179.1",      # Hydra
    "N_vectensis":       "XP_048579132.1",      # Sea anemone
}

DISPLAY_NAMES = {
    "H_sapiens":       "H. sapiens (Human)",
    "M_musculus":      "M. musculus (Mouse)",
    "X_laevis":        "X. laevis (Frog)",
    "P_marinus":       "P. marinus (Lamprey)",
    "C_intestinalis":  "C. intestinalis (Sea squirt)",
    "B_floridae":      "B. floridae (Amphioxus)",
    "S_purpuratus":    "S. purpuratus (Sea urchin)",
    "O_sinensis":      "O. sinensis (Octopus)",
    "D_melanogaster":  "D. melanogaster (Fruit fly)",
    "A_mellifera":     "A. mellifera (Honeybee)",
    "P_tepidariorum":  "P. tepidariorum (Spider)",
    "H_vulgaris":      "H. vulgaris (Hydra)",
    "N_vectensis":     "N. vectensis (Sea anemone)",
}


def fetch_sequences(accession_map):
    records = []
    print(f"Fetching {len(accession_map)} UBR3 sequences from NCBI...")
    for name, acc in accession_map.items():
        for attempt in range(3):
            try:
                handle = Entrez.efetch(db="protein", id=acc,
                                       rettype="fasta", retmode="text")
                text = handle.read(); handle.close()
                rec = SeqIO.read(StringIO(text), "fasta")
                rec.id = name
                rec.name = name
                rec.description = f"{name} | {acc}"
                records.append(rec)
                print(f"  OK  {name:22s}  {acc:22s}  ({len(rec.seq)} aa)")
                time.sleep(0.5)
                break
            except Exception as exc:
                print(f"  RETRY {name}: {exc}")
                time.sleep(3)
        else:
            print(f"  FAIL {name} after 3 attempts")
    return records


def run_clustalo_ebi(fasta_path):
    url_run    = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/run"
    url_status = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/status/{}"
    url_result = "https://www.ebi.ac.uk/Tools/services/rest/clustalo/result/{}/aln-fasta"

    with open(fasta_path, encoding="utf-8") as fh:
        seq_data = fh.read()

    data = parse.urlencode({
        "email": Entrez.email, "sequence": seq_data, "outfmt": "fa",
    }).encode()

    print("Submitting Clustal Omega alignment to EBI...")
    req = request.Request(url_run, data=data)
    try:
        resp = request.urlopen(req)
    except Exception as exc:
        if hasattr(exc, 'read'):
            body = exc.read().decode()
            print(f"  Server response: {body[:500]}")
        raise
    job_id = resp.read().decode().strip()
    print(f"  Job ID: {job_id}")

    for _ in range(200):
        time.sleep(10)
        resp = request.urlopen(url_status.format(job_id))
        status = resp.read().decode().strip()
        print(f"  Status: {status}")
        if status == "FINISHED":
            break
        if status in ("FAILURE", "ERROR", "NOT_FOUND"):
            raise RuntimeError(f"Job failed: {status}")
    else:
        raise RuntimeError("Timed out")

    resp = request.urlopen(url_result.format(job_id))
    return resp.read().decode()


def main():
    print("=" * 60)
    print("Rebuilding UBR3 alignment with broad phylogenetic panel")
    print("=" * 60)

    records = fetch_sequences(CURATED_ORTHOLOGS)
    if len(records) < 2:
        sys.exit("Need >= 2 sequences")

    fasta_path = os.path.join(OUTPUT_DIR, "ubr3_orthologs_broad.fasta")
    SeqIO.write(records, fasta_path, "fasta")
    print(f"\nSaved {len(records)} sequences to {fasta_path}")

    print(f"\nAligning {len(records)} sequences...")
    aligned_text = run_clustalo_ebi(fasta_path)

    aligned_path = os.path.join(OUTPUT_DIR, "ubr3_aligned.fasta")
    with open(aligned_path, "w", encoding="utf-8") as fh:
        fh.write(aligned_text)
    print(f"Alignment saved to {aligned_path}")

    alignment = AlignIO.read(aligned_path, "fasta")
    print(f"\nAlignment: {len(alignment)} sequences x "
          f"{alignment.get_alignment_length()} columns")

    print("\n" + "=" * 60)
    print("DONE - ready for publication_figures.py")
    print("=" * 60)


if __name__ == "__main__":
    main()
