import os
import json

# ==========================================
# INPUTS
# ==========================================

ACCESSION_MAP_FILE = "accession_to_files.json"

FASTA_DIR = "fasta"
METADATA_DIR = "metadata"

OUTPUT_DIR = "vmr_multifasta"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ==========================================
# LOAD ACCESSION MAP
# ==========================================

with open(
    ACCESSION_MAP_FILE,
    "r",
    encoding="utf-8"
) as f:

    accession_entries = json.load(f)

print(
    f"Loaded {len(accession_entries)} accession mappings"
)

# ==========================================
# CREATE EMPTY MULTIFASTA FILES
# ==========================================

all_vmr_files = set()

for entry in accession_entries:

    for vmr_file in entry["files"]:
        all_vmr_files.add(vmr_file)

# Create empty multifasta files
for vmr_file in sorted(all_vmr_files):

    out_name = vmr_file.replace(".json", ".fasta")

    out_path = os.path.join(
        OUTPUT_DIR,
        out_name
    )

    # Empty/create file
    open(out_path, "w").close()

print(
    f"Created {len(all_vmr_files)} empty multifasta files"
)

# ==========================================
# PROCESS ACCESSIONS
# ==========================================

written = 0
missing_fasta = 0
missing_metadata = 0

for i, entry in enumerate(accession_entries, 1):

    accession = entry["accession"]

    fasta_file = os.path.join(
        FASTA_DIR,
        f"{accession}.fasta"
    )

    metadata_file = os.path.join(
        METADATA_DIR,
        f"{accession}.json"
    )

    # --------------------------------------
    # CHECK FILES EXIST
    # --------------------------------------

    if not os.path.exists(fasta_file):

        missing_fasta += 1
        continue

    if not os.path.exists(metadata_file):

        missing_metadata += 1
        continue

    # --------------------------------------
    # LOAD METADATA
    # --------------------------------------

    try:

        with open(
            metadata_file,
            "r",
            encoding="utf-8"
        ) as f:

            metadata = json.load(f)

    except Exception as e:

        print(
            f"Metadata error {accession}: {e}"
        )

        continue

    # --------------------------------------
    # BUILD FASTA HEADER
    # --------------------------------------

    organism = metadata.get(
        "organism",
        "Unknown_organism"
    )

    description = metadata.get(
        "description",
        ""
    )

    sequence_length = metadata.get(
        "sequence_length",
        ""
    )

    header = (
        f">{accession}"
        f" | organism={organism}"
        f" | length={sequence_length}"
        f" | description={description}"
    )

    # --------------------------------------
    # LOAD GENOME SEQUENCE
    # --------------------------------------

    try:

        with open(
            fasta_file,
            "r",
            encoding="utf-8"
        ) as f:

            lines = f.readlines()

        # Remove original FASTA header
        sequence = "".join(
            line.strip()
            for line in lines
            if not line.startswith(">")
        )

    except Exception as e:

        print(
            f"FASTA error {accession}: {e}"
        )

        continue

    # --------------------------------------
    # CREATE FINAL FASTA ENTRY
    # --------------------------------------

    fasta_entry = (
        f"{header}\n"
        f"{sequence}\n"
    )

    # --------------------------------------
    # WRITE TO ALL CORRESPONDING VMR FILES
    # --------------------------------------

    for vmr_file in entry["files"]:

        out_name = vmr_file.replace(
            ".json",
            ".fasta"
        )

        out_path = os.path.join(
            OUTPUT_DIR,
            out_name
        )

        with open(
            out_path,
            "a",
            encoding="utf-8"
        ) as out:

            out.write(fasta_entry)

    written += 1

    if written % 1000 == 0:

        print(
            f"Processed {written} accessions"
        )

# ==========================================
# SUMMARY
# ==========================================

print("\nDONE\n")

print(f"Written accessions: {written}")

print(f"Missing FASTA files: {missing_fasta}")

print(f"Missing metadata files: {missing_metadata}")

print(
    f"Output directory: {OUTPUT_DIR}"
)
