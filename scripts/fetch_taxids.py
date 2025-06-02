# scripts/fetch_taxids.py

import json
import os
from time import sleep
from Bio import Entrez

# Use your GitHub secret email for NCBI API
Entrez.email = os.getenv("ENTREZ_EMAIL", "fallback@example.com")

json_path = "converted_files/merged_ictv.json"

# Load JSON data
with open(json_path, "r") as f:
    data = json.load(f)

def get_taxon_id_from_accession(accession):
    try:
        handle = Entrez.esearch(db="nuccore", term=accession)
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            print(f"Warning: No NCBI records found for accession: {accession}")
            return None

        uid = record["IdList"][0]
        handle = Entrez.esummary(db="nuccore", id=uid)
        summary = Entrez.read(handle)
        handle.close()

        if not summary or len(summary) == 0:
            print(f"Warning: No summary data found for accession: {accession}")
            return None

        taxid = summary[0].get("TaxId")
        if not taxid:
            print(f"Warning: TaxId not found in summary for accession: {accession}")
            return None

        return int(taxid)

    except Exception as e:
        print(f"Error fetching Taxon ID for {accession}: {e}")
        return None

# Enrich entries with Taxon_ID where possible
for entry in data:
    accession = entry.get("Virus_GENBANK_accession", "").strip()
    if accession:
        taxid = get_taxon_id_from_accession(accession)
        if taxid is not None:
            entry["Taxon_ID"] = taxid
            print(f"{accession} => Taxon ID: {taxid}")
        else:
            print(f"{accession} => No Taxon ID found, skipping.")
        sleep(0.34)  # ~3 requests per second (NCBI limit)

# Save updated JSON back to file
with open(json_path, "w") as f:
    json.dump(data, f, indent=2)
