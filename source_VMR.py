
import os
import re
import json

INPUT_DIR = "ncbi"

ACCESSION_PATTERN = re.compile(
    r'\b([A-Z]{1,4}_?\d{5,9}(?:\.\d+)?)\b'
)


def determine_isolate_type(exemplar_value):

    text = str(exemplar_value).strip()

    # Common VMR patterns
    if text.endswith("(A)"):
        return "A"

    if text.endswith("(E)"):
        return "E"

    # fallback heuristics
    lowered = text.lower()

    if "additional" in lowered:
        return "A"

    return "E"


def extract_accessions(text):

    return ACCESSION_PATTERN.findall(str(text))


all_records = []

for filename in os.listdir(INPUT_DIR):

    if not filename.endswith(".json"):
        continue

    filepath = os.path.join(INPUT_DIR, filename)

    with open(filepath, "r", encoding="utf-8") as f:
        data = json.load(f)

    for row in data:

        exemplar_value = row.get("Exemplar", "")

        isolate_type = determine_isolate_type(
            exemplar_value
        )

        for key, value in row.items():

            if "accession" not in key.lower():
                continue

            accession_type = (
                "RefSeq"
                if "refseq" in key.lower()
                else "GenBank"
            )

            accessions = extract_accessions(value)

            for accession in accessions:

                all_records.append({
                    "accession": accession,
                    "accession_type": accession_type,
                    "isolate_type": isolate_type,
                    "source_file": filename,
                    "exemplar_value": exemplar_value
                })


print(f"\nTotal accession records: {len(all_records)}")

# Deduplicated accession count
unique_accessions = {
    r["accession"]
    for r in all_records
}

print(
    f"Unique accessions: "
    f"{len(unique_accessions)}"
)

# Example
print("\nExample:\n")
print(all_records[0])

