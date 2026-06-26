import os
import json
import time
from Bio import Entrez, SeqIO

# NCBI credentials
Entrez.email = "shahram.saghaei@uni-jena.de"
Entrez.api_key = "e9d71c34d90febd81e75a45cd31146ce1208"

INPUT_DIR = "ncbi"
OUTPUT_FASTA = "fasta"
OUTPUT_METADATA = "metadata"

os.makedirs(OUTPUT_FASTA, exist_ok=True)
os.makedirs(OUTPUT_METADATA, exist_ok=True)


def fetch_fasta(accession):
    try:
        handle = Entrez.efetch(
            db="nucleotide",
            id=accession,
            rettype="fasta",
            retmode="text"
        )

        fasta_text = handle.read()
        handle.close()

        return fasta_text

    except Exception as e:
        print(f"FASTA error for {accession}: {e}")
        return None


def fetch_metadata(accession):
    try:
        handle = Entrez.efetch(
            db="nucleotide",
            id=accession,
            rettype="gb",
            retmode="text"
        )

        record = SeqIO.read(handle, "genbank")
        handle.close()

        metadata = {
            "accession": accession,
            "name": record.name,
            "description": record.description,
            "sequence_length": len(record.seq),
            "organism": record.annotations.get("organism"),
            "taxonomy": record.annotations.get("taxonomy"),
            "date": record.annotations.get("date"),
            "topology": record.annotations.get("topology"),
            "molecule_type": record.annotations.get("molecule_type"),
            "keywords": record.annotations.get("keywords"),
            "references": [
                {
                    "title": ref.title,
                    "authors": ref.authors,
                    "journal": ref.journal
                }
                for ref in record.annotations.get("references", [])
            ],
            "features": []
        }

        for feature in record.features:
            metadata["features"].append({
                "type": feature.type,
                "qualifiers": feature.qualifiers
            })

        return metadata

    except Exception as e:
        print(f"Metadata error for {accession}: {e}")
        return None


import re

ACCESSION_PATTERN = re.compile(
    r'\b([A-Z]{1,4}_?\d{5,9}(?:\.\d+)?)\b'
)


def extract_accessions(row):

    accessions = set()

    for key, value in row.items():

        if "accession" not in key.lower():
            continue

        if not value:
            continue

        text = str(value)

        matches = ACCESSION_PATTERN.findall(text)

        for acc in matches:
            accessions.add(acc.strip())

    return accessions



def process_json_files():

    all_accessions = set()

    # Read all output JSON files
    for filename in os.listdir(INPUT_DIR):

        if not filename.endswith(".json"):
            continue

        path = os.path.join(INPUT_DIR, filename)

        try:
            with open(path, "r", encoding="utf-8") as f:
                data = json.load(f)

            for row in data:

                row_accessions = extract_accessions(row)

                all_accessions.update(row_accessions)

        except Exception as e:
            print(f"Error reading {filename}: {e}")

    print(f"\nFound {len(all_accessions)} unique accession numbers\n")

    # Download sequences + metadata
    for i, accession in enumerate(sorted(all_accessions), 1):

        print(f"[{i}/{len(all_accessions)}] {accession}")

        # FASTA
        fasta = fetch_fasta(accession)

        if fasta:

            fasta_file = os.path.join(
                OUTPUT_FASTA,
                f"{accession}.fasta"
            )

            with open(fasta_file, "w") as f:
                f.write(fasta)

        # METADATA
        metadata = fetch_metadata(accession)

        if metadata:

            metadata_file = os.path.join(
                OUTPUT_METADATA,
                f"{accession}.json"
            )

            with open(metadata_file, "w", encoding="utf-8") as f:
                json.dump(metadata, f, indent=2)

        # NCBI rate limiting
        time.sleep(0.12)


if __name__ == "__main__":
    process_json_files()

