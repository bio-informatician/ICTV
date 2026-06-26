import os
import json
from collections import defaultdict

INPUT_DIR = "ncbi"
OUTPUT_FILE = "accession_to_files.json"


def clean_accession(raw_value):
    """
    Clean accession strings.

    Examples:
        "L:KP710246" -> KP710246
        " KP710246 " -> KP710246
        "KP710246.1" -> KP710246
        "LK928904 (nts 2253-10260)" -> LK928904
        "(nts 2253-10260):LK928904" -> LK928904
    """

    value = raw_value.strip()

    # Remove segment name before :
    if ":" in value:
        value = value.split(":")[-1].strip()

    # Remove parenthesis annotation
    if "(" in value:
        value = value.split("(")[0].strip()

    # Remove version number
    if "." in value:
        left, right = value.rsplit(".", 1)

        # Only remove if right side is numeric
        if right.isdigit():
            value = left

    return value.strip()


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

                # Only accession fields
                if "accession" not in key.lower():
                    continue

                # Skip empty values
                if not value:
                    continue

                text = str(value).strip()

                if not text:
                    continue

                # Split multiple entries
                # Example:
                # L:KP710246; M:KP710264; S:KP710267
                parts = text.split(";")

                for part in parts:

                    clean_acc = clean_accession(part)

                    if not clean_acc:
                        continue

                    accession_map[clean_acc].add(
                        filename
                    )

    except Exception as e:

        print(f"ERROR {filename}: {e}")

# =========================
# CONVERT TO OUTPUT
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
# WRITE JSON
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
