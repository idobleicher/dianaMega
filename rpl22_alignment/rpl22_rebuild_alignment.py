"""
Fetch RPL22 and RPL22L1 orthologs from 13 species and align via Clustal Omega (EBI).
Produces rpl22_aligned.fasta and rpl22l1_aligned.fasta.
"""
import os, sys, time, re
from io import StringIO
from collections import Counter

from Bio import SeqIO, AlignIO, Entrez
from Bio.SeqRecord import SeqRecord
from urllib import request, parse

import numpy as np

Entrez.email = "rpl22.alignment@example.com"
Entrez.tool  = "RPL22_ortholog_alignment"

OUTPUT_DIR = os.path.dirname(os.path.abspath(__file__))

# ── Curated accessions ────────────────────────────────────────────────────
# RPL22 orthologs (60S ribosomal protein L22 / eL22)
RPL22_ORTHOLOGS = {
    "H_sapiens":         "NP_000974.1",          # Human (128 aa)
    "M_musculus":        "NP_033105.1",           # Mouse (128 aa)
    "X_laevis":          "NP_001081541.1",        # Frog (128 aa)
    "P_marinus":         "XP_032817946.1",        # Lamprey (135 aa)
    "C_intestinalis":    "XP_002125710.2",        # Sea squirt (131 aa, eL22-like best hit)
    "B_floridae":        "XP_035663727.1",        # Amphioxus (132 aa, L22-like)
    "S_purpuratus":      "XP_782165.2",           # Sea urchin (130 aa)
    "O_sinensis":        "XP_029648559.1",        # Octopus (136 aa)
    "D_melanogaster":    "NP_477134.1",           # Fruit fly (299 aa, has long N-term extension)
    "A_mellifera":       "XP_625009.2",           # Honeybee (139 aa)
    "P_tepidariorum":    "XP_015921618.1",        # Spider (135 aa)
    "H_vulgaris":        "CDG71062.1",            # Hydra (124 aa)
    "N_vectensis":       "XP_001626833.2",        # Sea anemone (123 aa)
}

# RPL22L1 orthologs (60S ribosomal protein L22-like 1)
# RPL22L1 is primarily a vertebrate paralog; fewer species have it
RPL22L1_ORTHOLOGS = {
    "H_sapiens":         "NP_001093115.1",        # Human (122 aa)
    "M_musculus":        "NP_080793.1",            # Mouse (122 aa)
    "X_laevis":          "NP_001091257.1",         # Frog (123 aa)
    "G_gallus":          "NP_001264646.1",         # Chicken (122 aa)
    "D_rerio":           "NP_001038800.1",         # Zebrafish (126 aa)
    "C_intestinalis":    "XP_002125710.2",         # Sea squirt (131 aa)
    "B_floridae":        "XP_035663727.1",         # Amphioxus (132 aa)
}

DISPLAY_NAMES = {
    "H_sapiens":       "H. sapiens (Human)",
    "M_musculus":      "M. musculus (Mouse)",
    "X_laevis":        "X. laevis (Frog)",
    "G_gallus":        "G. gallus (Chicken)",
    "D_rerio":         "D. rerio (Zebrafish)",
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


def fetch_sequences(accession_map, protein_name):
    records = []
    print(f"Fetching {len(accession_map)} {protein_name} sequences from NCBI...")
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


def search_ortholog(gene_name, organism, db="protein"):
    """Search NCBI for an ortholog when curated accession fails."""
    query = f"{gene_name}[Gene Name] AND {organism}[Organism] AND refseq[filter]"
    try:
        handle = Entrez.esearch(db=db, term=query, retmax=5)
        result = Entrez.read(handle)
        handle.close()
        if result["IdList"]:
            acc_id = result["IdList"][0]
            handle = Entrez.efetch(db=db, id=acc_id, rettype="fasta", retmode="text")
            text = handle.read(); handle.close()
            rec = SeqIO.read(StringIO(text), "fasta")
            return rec
    except Exception:
        pass
    return None


def fetch_with_fallback(accession_map, protein_name, gene_name):
    """Fetch sequences, falling back to NCBI search if accession fails."""
    records = []
    organism_map = {
        "H_sapiens": "Homo sapiens", "M_musculus": "Mus musculus",
        "X_laevis": "Xenopus laevis", "P_marinus": "Petromyzon marinus",
        "C_intestinalis": "Ciona intestinalis", "B_floridae": "Branchiostoma floridae",
        "S_purpuratus": "Strongylocentrotus purpuratus",
        "O_sinensis": "Octopus sinensis",
        "D_melanogaster": "Drosophila melanogaster",
        "A_mellifera": "Apis mellifera",
        "P_tepidariorum": "Parasteatoda tepidariorum",
        "H_vulgaris": "Hydra vulgaris", "N_vectensis": "Nematostella vectensis",
    }
    print(f"Fetching {len(accession_map)} {protein_name} sequences from NCBI...")
    for name, acc in accession_map.items():
        success = False
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
                success = True
                break
            except Exception as exc:
                print(f"  RETRY {name}: {exc}")
                time.sleep(3)
        if not success:
            print(f"  Accession failed for {name}, searching NCBI...")
            organism = organism_map.get(name, name.replace("_", " "))
            rec = search_ortholog(gene_name, organism)
            if rec:
                rec.id = name
                rec.name = name
                rec.description = f"{name} | search"
                records.append(rec)
                print(f"  SEARCH OK  {name:22s}  ({len(rec.seq)} aa)")
            else:
                print(f"  FAIL {name} - no sequence found")
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


def align_protein(accession_map, protein_name, gene_name, output_prefix):
    print(f"\n{'='*60}")
    print(f"Processing {protein_name}")
    print(f"{'='*60}")

    records = fetch_with_fallback(accession_map, protein_name, gene_name)
    if len(records) < 2:
        print(f"Need >= 2 sequences for {protein_name}, got {len(records)}")
        return None

    fasta_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_orthologs.fasta")
    SeqIO.write(records, fasta_path, "fasta")
    print(f"\nSaved {len(records)} sequences to {fasta_path}")

    print(f"\nAligning {len(records)} sequences...")
    aligned_text = run_clustalo_ebi(fasta_path)

    aligned_path = os.path.join(OUTPUT_DIR, f"{output_prefix}_aligned.fasta")
    with open(aligned_path, "w", encoding="utf-8") as fh:
        fh.write(aligned_text)
    print(f"Alignment saved to {aligned_path}")

    alignment = AlignIO.read(aligned_path, "fasta")
    print(f"Alignment: {len(alignment)} sequences x "
          f"{alignment.get_alignment_length()} columns")

    return aligned_path


def main():
    print("=" * 60)
    print("RPL22 / RPL22L1 Ortholog Alignment Pipeline")
    print("=" * 60)

    rpl22_path = align_protein(
        RPL22_ORTHOLOGS, "RPL22", "RPL22", "rpl22")

    rpl22l1_path = align_protein(
        RPL22L1_ORTHOLOGS, "RPL22L1", "RPL22L1", "rpl22l1")

    print(f"\n{'='*60}")
    print("DONE")
    if rpl22_path:
        print(f"  RPL22 alignment:   {rpl22_path}")
    if rpl22l1_path:
        print(f"  RPL22L1 alignment: {rpl22l1_path}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
