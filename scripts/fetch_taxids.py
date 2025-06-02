import json
import requests
from time import sleep
from Bio import Entrez

import os
Entrez.email = os.getenv("ENTREZ_EMAIL", "fallback@example.com")

# Download JSON file from URL
url = "https://raw.githubusercontent.com/bio-informatician/ICTV/refs/heads/main/converted_files/merged_ictv.json"
response = requests.get(url)
data = response.json()

def get_taxon_id_from_accession(accession):
    try:
        # Search the nucleotide database using the accession number
        handle = Entrez.esearch(db="nuccore", term=accession)
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            return None

        # Fetch the full record to extract Taxonomy info
        uid = record["IdList"][0]
        handle = Entrez.esummary(db="nuccore", id=uid)
        summary = Entrez.read(handle)
        handle.close()

        taxid = summary[0].get("TaxId", None)
        return taxid
    except Exception as e:
        print(f"Error processing accession {accession}: {e}")
        return None

# Iterate over each entry and collect Taxon IDs
results = []
for entry in data:
    accession = entry.get("Virus_GENBANK_accession", "").strip()
    if accession:
        taxid = get_taxon_id_from_accession(accession)
        print(f"{accession} => Taxon ID: {taxid}")
        results.append({"accession": accession, "taxon_id": taxid})
        sleep(0.34)  # Respect NCBI rate limit (~3 requests/sec)

# Optional: save results to JSON
with open("accession_to_taxid.json", "w") as f:
    json.dump(results, f, indent=2)
