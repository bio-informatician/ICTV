import os
import json
import time
import re
from Bio import Entrez, SeqIO
# ============================================================
# NCBI credentials
# ============================================================
Entrez.email = "shahram.saghaei@uni-jena.de"
Entrez.api_key = "e9d71c34d90febd81e75a45cd31146ce1208"
# ============================================================
# Directories
# ============================================================
INPUT_DIR = "ncbi"
OUTPUT_FASTA = "fasta"
OUTPUT_METADATA = "metadata"
os.makedirs(INPUT_DIR, exist_ok=True)
os.makedirs(OUTPUT_FASTA, exist_ok=True)
os.makedirs(OUTPUT_METADATA, exist_ok=True)
# ============================================================
# Accession pattern
#
# Examples matched:
#   MT241517
#   AY225133
#   ABC12345
#   NC_001234
#
# Segment names such as:
#   L
#   M
#   S1
#   RNA10
#   SEG10
#   SEGA1
#
# are NOT matched because they do not contain enough digits.
# ============================================================
ACCESSION_PATTERN = re.compile(
    r"\b[A-Z]+_?\d{3,}\b"
)
# ============================================================
# Extract accessions from one JSON row
# ============================================================
def extract_accessions(row):
    accessions = set()
    if not isinstance(row, dict):
        return accessions
    for key, value in row.items():
        # Only examine fields whose name contains "accession"
        if "accession" not in key.lower():
            continue
        if value is None:
            continue
        text = str(value).upper()
        # Find all accession-like strings in the field
        for accession in ACCESSION_PATTERN.findall(text):
            # Keep only accessions with at least 5 characters
            if len(accession) >= 5:
                accessions.add(accession)
    return accessions
# ============================================================
# Fetch FASTA from NCBI
# ============================================================
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
        print(f"  FASTA error for {accession}: {e}")
        return None
# ============================================================
# Fetch GenBank metadata from NCBI
# ============================================================
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
        print(f"  Metadata error for {accession}: {e}")
        return None
# ============================================================
# Create accession list for ONE JSON file
# ============================================================
def create_accession_list(json_file):
    json_path = os.path.join(INPUT_DIR, json_file)
    base_name = os.path.splitext(json_file)[0]
    accession_file = os.path.join(
        INPUT_DIR,
        f"{base_name}.txt"
    )
    print(f"\nReading: {json_file}")
    try:
        with open(
            json_path,
            "r",
            encoding="utf-8"
        ) as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error reading {json_file}: {e}")
        return None
    all_accessions = set()
    # JSON is expected to contain a list of rows
    if isinstance(data, list):
        for row in data:
            all_accessions.update(
                extract_accessions(row)
            )
    # In case the JSON contains a single object
    elif isinstance(data, dict):
        all_accessions.update(
            extract_accessions(data)
        )
    else:
        print(
            f"Unsupported JSON structure in {json_file}"
        )
        return None
    accessions = sorted(all_accessions)
    # --------------------------------------------------------
    # IMPORTANT:
    # Write the accession list BEFORE downloading anything.
    # --------------------------------------------------------
    with open(
        accession_file,
        "w",
        encoding="utf-8"
    ) as f:
        for accession in accessions:
            f.write(accession + "\n")
    print(
        f"Found {len(accessions)} unique accessions"
    )
    print(
        f"Accession list saved to: {accession_file}"
    )
    return accession_file
# ============================================================
# Read accession list from .txt
# ============================================================
def read_accession_list(accession_file):
    accessions = []
    with open(
        accession_file,
        "r",
        encoding="utf-8"
    ) as f:
        for line in f:
            accession = line.strip()
            if not accession:
                continue
            accessions.append(accession)
    return accessions
# ============================================================
# Download/process one accession list
# ============================================================
def process_accession_list(accession_file):
    accessions = read_accession_list(
        accession_file
    )
    print(
        f"\nProcessing {len(accessions)} accessions "
        f"from {os.path.basename(accession_file)}"
    )
    for i, accession in enumerate(
        accessions,
        1
    ):
        print(
            f"[{i}/{len(accessions)}] {accession}"
        )
        fasta_file = os.path.join(
            OUTPUT_FASTA,
            f"{accession}.fasta"
        )
        metadata_file = os.path.join(
            OUTPUT_METADATA,
            f"{accession}.json"
        )
        # ----------------------------------------------------
        # Both files already exist
        # ----------------------------------------------------
        if (
            os.path.exists(fasta_file)
            and
            os.path.exists(metadata_file)
        ):
            print(
                "  FASTA + metadata already exist. Skipping."
            )
            continue
        # ----------------------------------------------------
        # FASTA
        # ----------------------------------------------------
        if not os.path.exists(fasta_file):
            print("  Downloading FASTA...")
            fasta = fetch_fasta(accession)
            if fasta:
                with open(
                    fasta_file,
                    "w",
                    encoding="utf-8"
                ) as f:
                    f.write(fasta)
                print("  FASTA saved.")
            else:
                print("  FASTA download failed.")
        else:
            print(
                "  FASTA already exists. Skipping."
            )
        # ----------------------------------------------------
        # Metadata
        # ----------------------------------------------------
        if not os.path.exists(metadata_file):
            print("  Downloading metadata...")
            metadata = fetch_metadata(accession)
            if metadata:
                with open(
                    metadata_file,
                    "w",
                    encoding="utf-8"
                ) as f:
                    json.dump(
                        metadata,
                        f,
                        indent=2,
                        ensure_ascii=False
                    )
                print("  Metadata saved.")
            else:
                print("  Metadata download failed.")
        else:
            print(
                "  Metadata already exists. Skipping."
            )
        # ----------------------------------------------------
        # NCBI rate limiting
        # 0.12 sec ≈ 8 requests/sec
        # ----------------------------------------------------
        time.sleep(0.12)
# ============================================================
# Main processing
# ============================================================
def process_json_files():
    json_files = sorted(
        filename
        for filename in os.listdir(INPUT_DIR)
        if filename.endswith(".json")
    )
    if not json_files:
        print(
            f"No JSON files found in {INPUT_DIR}/"
        )
        return
    print(
        f"Found {len(json_files)} JSON file(s)"
    )
    # ========================================================
    # STEP 1
    #
    # Create ALL accession .txt files first.
    #
    # NO DOWNLOADS happen in this step.
    # ========================================================
    accession_files = []
    print("\n" + "=" * 60)
    print("STEP 1: Creating accession lists")
    print("=" * 60)
    for json_file in json_files:
        accession_file = create_accession_list(
            json_file
        )
        if accession_file:
            accession_files.append(
                accession_file
            )
    # ========================================================
    # STEP 2
    #
    # Now download/process the accession lists.
    # ========================================================
    print("\n" + "=" * 60)
    print("STEP 2: Downloading missing data")
    print("=" * 60)
    for accession_file in accession_files:
        process_accession_list(
            accession_file
        )
# ============================================================
# Run
# ============================================================
if __name__ == "__main__":
    process_json_files()
