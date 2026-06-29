#!/usr/bin/env python3

import os
import json
import requests
from pathlib import Path

INPUT_DIR = "vmr_multifasta"
EXISTS_DIR = "existsindb"
NOT_EXISTS_DIR = "notindb"

API_URL = "https://api2.virjendb.org/v2/search"

# create output directories
os.makedirs(EXISTS_DIR, exist_ok=True)
os.makedirs(NOT_EXISTS_DIR, exist_ok=True)


def parse_fasta(filepath):
    """
    Simple FASTA parser preserving sequence formatting.
    Returns list of tuples: (header_without_>, sequence_string)
    """
    records = []

    with open(filepath, "r") as f:
        header = None
        seq_lines = []

        for line in f:
            line = line.rstrip("\n")

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "\n".join(seq_lines)))

                header = line[1:]  # remove >
                seq_lines = []
            else:
                seq_lines.append(line)

        if header is not None:
            records.append((header, "\n".join(seq_lines)))

    return records


def get_insdc_accession(header):
    """
    Extract accession from beginning of header.
    Example:
    AB000923 | organism=...
    """
    return header.split("|")[0].strip()


def query_virjendb(accession):
    payload = {
        "query": f"insdc_accession:{accession}"
    }

    try:
        response = requests.post(
            API_URL,
            headers={
                "accept": "application/json",
                "Content-Type": "application/json"
            },
            json=payload,
            timeout=30
        )

        response.raise_for_status()

        data = response.json()

        # adjust according to API response structure
        if isinstance(data, list):
            results = data
        elif "results" in data:
            results = data["results"]
        elif "data" in data:
            results = data["data"]
        else:
            results = []

        return results

    except Exception as e:
        print(f"ERROR querying {accession}: {e}")
        return None


def extract_vj_accession(record):
    """
    Extract VirJenDB accession from API record.
    """
    return record.get("VirJenDB Accession")


for fasta_file in Path(INPUT_DIR).glob("*"):

    if not fasta_file.is_file():
        continue

    print(f"Processing {fasta_file.name}")

    records = parse_fasta(fasta_file)

    exists_output = []
    not_exists_output = []

    for header, sequence in records:

        accession = get_insdc_accession(header)

        results = query_virjendb(accession)

        if results is None:
            # API failure -> treat as not found
            not_exists_output.append((header, sequence))
            continue

        if len(results) == 1:

            vj_accession = extract_vj_accession(results[0])

            if vj_accession:

                new_header = f"{vj_accession} | {header}"

                exists_output.append((new_header, sequence))

            else:
                not_exists_output.append((header, sequence))

        else:
            # zero or multiple results
            not_exists_output.append((header, sequence))

    # write existsindb output
    exists_path = Path(EXISTS_DIR) / fasta_file.name

    with open(exists_path, "w") as out:
        for header, sequence in exists_output:
            out.write(f">{header}\n")
            out.write(f"{sequence}\n")

    # write notindb output
    not_exists_path = Path(NOT_EXISTS_DIR) / fasta_file.name

    with open(not_exists_path, "w") as out:
        for header, sequence in not_exists_output:
            out.write(f">{header}\n")
            out.write(f"{sequence}\n")

print("Done.")