
import os
import time
import requests

from pathlib import Path
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


INPUT_DIR = "vmr_multifasta"
EXISTS_DIR = "existsindb"
NOT_EXISTS_DIR = "notindb"

API_URL = "https://api2.virjendb.org/v2/search"


# ------------------------------------------------------------------
# create output directories
# ------------------------------------------------------------------

os.makedirs(EXISTS_DIR, exist_ok=True)
os.makedirs(NOT_EXISTS_DIR, exist_ok=True)


# ------------------------------------------------------------------
# requests session with retry
# ------------------------------------------------------------------

session = requests.Session()

retries = Retry(
    total=5,
    backoff_factor=2,
    status_forcelist=[429, 500, 502, 503, 504],
    allowed_methods=["POST"]
)

adapter = HTTPAdapter(max_retries=retries)

session.mount("https://", adapter)
session.mount("http://", adapter)


# ------------------------------------------------------------------
# cache
# ------------------------------------------------------------------

cache = {}


# ------------------------------------------------------------------
# fasta parser
# ------------------------------------------------------------------

def parse_fasta(filepath):

    with open(filepath, "r") as f:

        header = None
        seq_lines = []

        for line in f:

            line = line.rstrip("\n")

            if line.startswith(">"):

                if header is not None:
                    yield header, "\n".join(seq_lines)

                header = line[1:]   # remove >
                seq_lines = []

            else:
                seq_lines.append(line)

        if header is not None:
            yield header, "\n".join(seq_lines)


# ------------------------------------------------------------------
# accession extraction
# ------------------------------------------------------------------

def get_insdc_accession(header):

    return header.split("|")[0].strip()


# ------------------------------------------------------------------
# query virjendb
# ------------------------------------------------------------------

def query_virjendb(accession):

    if accession in cache:
        return cache[accession]

    payload = {
        "query": f"insdc_accession:{accession}"
    }

    try:

        response = session.post(
            API_URL,
            headers={
                "accept": "application/json",
                "Content-Type": "application/json"
            },
            json=payload,
            timeout=(10, 60)
        )

        response.raise_for_status()

        data = response.json()

        # adjust for API structure
        if isinstance(data, list):
            results = data

        elif "results" in data:
            results = data["results"]

        elif "data" in data:
            results = data["data"]

        else:
            results = []

        cache[accession] = results

        # small delay to avoid hammering API
        time.sleep(0.2)

        return results

    except Exception as e:

        print(f"ERROR querying {accession}: {e}")

        cache[accession] = None

        return None


# ------------------------------------------------------------------
# extract VirJenDB accession
# ------------------------------------------------------------------

def extract_vj_accession(record):

    return record.get("VirJenDB Accession")


# ------------------------------------------------------------------
# processing
# ------------------------------------------------------------------

total_found = 0
total_not_found = 0
total_errors = 0


fasta_files = sorted(Path(INPUT_DIR).glob("*"))


for fasta_file in fasta_files:

    if not fasta_file.is_file():
        continue

    print("\n==================================================")
    print(f"Processing file: {fasta_file.name}")
    print("==================================================")

    exists_path = Path(EXISTS_DIR) / fasta_file.name
    not_exists_path = Path(NOT_EXISTS_DIR) / fasta_file.name

    exists_handle = open(exists_path, "w")
    not_exists_handle = open(not_exists_path, "w")

    record_count = 0

    for header, sequence in parse_fasta(fasta_file):

        record_count += 1

        accession = get_insdc_accession(header)

        print(f"[{record_count}] Checking {accession}")

        results = query_virjendb(accession)

        # ----------------------------------------------------------
        # API failure
        # ----------------------------------------------------------

        if results is None:

            total_errors += 1

            not_exists_handle.write(f">{header}\n")
            not_exists_handle.write(f"{sequence}\n")

            not_exists_handle.flush()

            print("  -> ERROR")

            continue

        # ----------------------------------------------------------
        # exactly one result
        # ----------------------------------------------------------

        if len(results) == 1:

            vj_accession = extract_vj_accession(results[0])

            if vj_accession:

                new_header = f"{vj_accession} | {header}"

                exists_handle.write(f">{new_header}\n")
                exists_handle.write(f"{sequence}\n")

                exists_handle.flush()

                total_found += 1

                print(f"  -> FOUND ({vj_accession})")

            else:

                not_exists_handle.write(f">{header}\n")
                not_exists_handle.write(f"{sequence}\n")

                not_exists_handle.flush()

                total_not_found += 1

                print("  -> NO VirJenDB accession")

        # ----------------------------------------------------------
        # zero or multiple results
        # ----------------------------------------------------------

        else:

            not_exists_handle.write(f">{header}\n")
            not_exists_handle.write(f"{sequence}\n")

            not_exists_handle.flush()

            total_not_found += 1

            print(f"  -> NOT FOUND ({len(results)} results)")

    exists_handle.close()
    not_exists_handle.close()

    print(f"Finished {fasta_file.name}")
    print(f"Records processed: {record_count}")


# ------------------------------------------------------------------
# summary
# ------------------------------------------------------------------

print("\n==================================================")
print("DONE")
print("==================================================")

print(f"Found in DB     : {total_found}")
print(f"Not found       : {total_not_found}")
print(f"Errors          : {total_errors}")
print(f"Cache size      : {len(cache)}")
