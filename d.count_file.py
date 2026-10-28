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
        KP710246
        NC_011268
        KP710246.1
        LK928904 (nts 2253-10260)
        (nts 2253-10260):LK928904
        , NC_011268
    """

    value = str(raw_value).strip()

    if not value:
        return None

    # Remove segment labels
    # L:KP710246 -> KP710246
    if ":" in value:
        value = value.split(":")[-1].strip()

    # Remove parenthesis annotations
    # LK928904 (nts 2253-10260) -> LK928904
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

    # Extract valid accession
    match = ACCESSION_PATTERN.search(value)

    if match:
        return match.group(1)

    return None


# accession -> set(VMR files)
accession_map = defaultdict(set)

# =========================
# PROCESS TXT FILES
# =========================

txt_files = [
    filename
    for filename in os.listdir(INPUT_DIR)
    if filename.endswith(".txt")
]

print(f"Found {len(txt_files)} TXT files")

for filename in sorted(txt_files):

    filepath = os.path.join(INPUT_DIR, filename)

    # Corresponding VMR file name
    # Example:
    # VMR_MSL41.v1.20260320.txt
    # ->
    # VMR_MSL41.v1.20260320.json
    vmr_filename = filename[:-4] + ".json"

    print(f"Processing {filename}")

    try:

        with open(filepath, "r", encoding="utf-8") as f:

            for line in f:

                line = line.strip()

                if not line:
                    continue

                # In case a line contains multiple accessions
                # separated by semicolon
                parts = line.split(";")

                for part in parts:

                    accession = clean_accession(part)

                    if not accession:
                        continue

                    accession_map[accession].add(
                        vmr_filename
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
        "files": sorted(accession_map[accession])
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


# =========================
# SUMMARY
# =========================

print()
print("DONE")
print(f"TXT files processed: {len(txt_files)}")
print(f"Unique accessions:   {len(output)}")
print(f"Output written to:   {OUTPUT_FILE}")
