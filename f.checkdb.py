#!/usr/bin/env python3

import os
import time
import requests

from pathlib import Path
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


INPUT_DIR = "vmr_multifasta"
EXISTS_DIR = "existsindb"
NOT_EXISTS_DIR = "notindb"

API_URL = "https://api.dev.virjendb.org/v2/search"


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

    accession = header.split("|")[0].strip()

    # remove version suffix if present
    # AB000403.1 -> AB000403
    accession = accession.split(".")[0]

    return accession


# ------------------------------------------------------------------
# query virjendb
# ------------------------------------------------------------------

def query_virjendb(accession):

    # cache hit
    if accession in cache:

        print("  -> CACHE HIT")

        return cache[accession]

    payload = {
        "query": f"insdc_accession:{accession}"
    }

    print("  -> Querying API...", end="", flush=True)

    start_time = time.time()

    try:

        response = session.post(
            API_URL,
            headers={
                "accept": "application/json",
                "Content-Type": "application/json"
            },
            json=payload,

            # (connect timeout, read timeout)
            timeout=(5, 30)
        )

        elapsed = round(time.time() - start_time, 2)

        print(f" done ({elapsed}s)")

        response.raise_for_status()

        data = response.json()

        results = data.get("results", [])

        cache[accession] = results

        # polite delay
        time.sleep(0.2)

        return results

    except Exception as e:

        elapsed = round(time.time() - start_time, 2)

        print(f" FAILED ({elapsed}s)")
        print(f"  -> ERROR: {e}")

        cache[accession] = None

        return None


# ------------------------------------------------------------------
# extract VirJenDB accession
# ------------------------------------------------------------------

def extract_vj_accession(record):

    return record.get("source", {}).get("VirJenDB Accession")


# ------------------------------------------------------------------
# validate accession match
# ------------------------------------------------------------------

def validate_match(results, accession):

    for r in results:

        source = r.get("source", {})

        insdc_list = source.get("INSDC Accession", [])

        normalized_db_accessions = [
            x.split(".")[0] for x in insdc_list
        ]

        if accession in normalized_db_accessions:

            return r

    return None


# ------------------------------------------------------------------
# processing
# ------------------------------------------------------------------

total_found = 0
total_not_found = 0
total_errors = 0

fasta_files = sorted( Path(INPUT_DIR).glob("*"), reverse=True )

for fasta_file in fasta_files:

    if not fasta_file.is_file():
        continue

    print("\n==================================================")
    print(f"Processing file: {fasta_file.name}")
    print("==================================================")

    exists_path = Path(EXISTS_DIR) / fasta_file.name
    not_exists_path = Path(NOT_EXISTS_DIR) / fasta_file.name
    exists_txt_path = Path(EXISTS_DIR) / f"{fasta_file.stem}.txt"
    not_exists_txt_path = Path(NOT_EXISTS_DIR) / f"{fasta_file.stem}.txt"

    fasta_exists = exists_path.exists() and not_exists_path.exists()
    txt_exists = exists_txt_path.exists() and not_exists_txt_path.exists()

    if fasta_exists and txt_exists:
        print(f"Skipping {fasta_file.name}: ALL output files already exist.")
        continue
        
    elif fasta_exists and not txt_exists:
        print(f"FASTA exists but TXT missing for {fasta_file.name}. Generating TXT only...")
        
        # 1. Generate TXT from existsindb FASTA
        with open(exists_txt_path, "w") as txt_out:
            for out_header, _ in parse_fasta(exists_path):
                # The exists FASTA header looks like: VJ_ACC | INSDC_ACC | ...
                # So we split by "|" and take the SECOND item (index 1)
                parts = out_header.split("|")
                if len(parts) >= 2:
                    insdc = parts[1].strip().split(".")[0]
                    txt_out.write(f"{insdc}\n")
                    
        # 2. Generate TXT from notindb FASTA
        with open(not_exists_txt_path, "w") as txt_out:
            for out_header, _ in parse_fasta(not_exists_path):
                # The not_exists FASTA header looks like: INSDC_ACC | ...
                # So we split by "|" and take the FIRST item (index 0)
                insdc = out_header.split("|")[0].strip().split(".")[0]
                txt_out.write(f"{insdc}\n")
                
        print(" -> TXT files created successfully. Skipping API queries.")
        continue

    exists_handle = open(exists_path, "w")
    not_exists_handle = open(not_exists_path, "w")

    exists_txt_handle = open(exists_txt_path, "w")
    not_exists_txt_handle = open(not_exists_txt_path, "w")

    record_count = 0

    for header, sequence in parse_fasta(fasta_file):

        record_count += 1

        accession = get_insdc_accession(header)

        print(f"\n[{record_count}] Checking {accession}")

        results = query_virjendb(accession)

        # ----------------------------------------------------------
        # API failure
        # ----------------------------------------------------------

        if results is None:

            total_errors += 1

            not_exists_handle.write(f">{header}\n")
            not_exists_handle.write(f"{sequence}\n")

            not_exists_handle.flush()

            not_exists_txt_handle.write(f"{accession}\n")
            not_exists_txt_handle.flush()

            print("  -> STORED IN notindb (API ERROR)")

            continue

        # ----------------------------------------------------------
        # validate accession match
        # ----------------------------------------------------------

        matched_record = validate_match(results, accession)

        # ----------------------------------------------------------
        # valid match found
        # ----------------------------------------------------------

        if matched_record is not None:

            vj_accession = extract_vj_accession(matched_record)

            if vj_accession:

                new_header = f"{vj_accession} | {header}"

                exists_handle.write(f">{new_header}\n")
                exists_handle.write(f"{sequence}\n")

                exists_handle.flush()

                exists_txt_handle.write(f"{accession}\n")
                exists_txt_handle.flush()

                total_found += 1

                print(f"  -> FOUND ({vj_accession})")
                print("  -> STORED IN existsindb")

            else:

                not_exists_handle.write(f">{header}\n")
                not_exists_handle.write(f"{sequence}\n")

                not_exists_handle.flush()

                not_exists_txt_handle.write(f"{accession}\n")
                not_exists_txt_handle.flush()

                total_not_found += 1

                print("  -> MATCHED but NO VirJenDB accession")
                print("  -> STORED IN notindb")

        # ----------------------------------------------------------
        # no validated match
        # ----------------------------------------------------------

        else:

            not_exists_handle.write(f">{header}\n")
            not_exists_handle.write(f"{sequence}\n")

            not_exists_handle.flush()

            not_exists_txt_handle.write(f"{accession}\n")
            not_exists_txt_handle.flush()
            
            total_not_found += 1

            print("  -> NOT FOUND")
            print("  -> STORED IN notindb")

    exists_handle.close()
    not_exists_handle.close()
    exists_txt_handle.close()       
    not_exists_txt_handle.close()

    print(f"\nFinished {fasta_file.name}")
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
