import os
import re
import json
from collections import defaultdict

INPUT_DIR = "ncbi"
OUTPUT_FILE = "accession_to_files.json"

# Valid accession patterns:
# KP710246
# NC_011268
# AB123456
# etc.
ACCESSION_PATTERN = re.compile(
    r'([A-Z]{1,4}_?\d{5,9})'
)


def clean_accession(raw_value):
    """
    Clean accession strings.

    Handles:
        L:KP710246
        M: KP710264
        KP710246.1
        LK928904 (nts 2253-10260)
        (nts 2253-10260):LK928904
        , NC_011268
    """

    value = str(raw_value).strip()

    # Remove segment labels
    # L:KP710246 -> KP710246
    if ":" in value:
        value = value.split(":")[-1].strip()

    # Remove parenthesis annotations
    # LK928904 (nts 2253-10260)
    if "(" in value:
        value = value.split("(")[0].strip()

    # Remove version number
    # NC_011268.1 -> NC_011268
    if "." in value:
        left, right = value.rsplit(".", 1)

        if right.isdigit():
            value = left

    # Remove leading/trailing junk
    value = value.strip(" ,;")

    # Extract only valid accession
    match = ACCESSION_PATTERN.search(value)

    if match:
        return match.group(1)

    return None


# accession -> set(files)
accession_map = defaultdict(set)

# =========================
# PROCESS FILES
# =========================

for filename in os.listdir(INPUT_DIR):

    if not filename.endswith(".json"):
        continue

    filepath = os.path.join(INPUT_DIR, filename)

    print(f"Processing {filename}")

    try:

        with open(filepath, "r", encoding="utf-8") as f:
            data = json.load(f)

        for row in data:

            for key, value in row.items():

                # Only accession columns
                if "accession" not in key.lower():
                    continue

                # Skip empty fields
                if not value:
                    continue

                text = str(value).strip()

                if not text:
                    continue

                # Split multiple accessions
                # Example:
                # L:KP710246; M:KP710264; S:KP710267
                parts = text.split(";")

                for part in parts:

                    accession = clean_accession(part)

                    if not accession:
                        continue

                    accession_map[accession].add(
                        filename
                    )

    except Exception as e:

        print(f"ERROR {filename}: {e}")

# =========================
# BUILD OUTPUT
# =========================

output = []

for accession in sorted(accession_map):

    output.append({
        "accession": accession,
        "files": sorted(
            accession_map[accession]
        )
    })

# =========================
# WRITE OUTPUT
# =========================

with open(
    OUTPUT_FILE,
    "w",
    encoding="utf-8"
) as f:

    json.dump(
        output,
        f,
        indent=4
    )

print(
    f"\nDONE\n"
    f"Unique accessions: {len(output)}\n"
    f"Output written to: {OUTPUT_FILE}"
)
