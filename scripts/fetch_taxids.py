import json
import os
from time import sleep
from Bio import Entrez

# Use your GitHub secret email for NCBI API
Entrez.email = os.getenv("ENTREZ_EMAIL", "fallback@example.com")

json_path = "converted_files/merged_ictv.json"
backup_path = "converted_files/merged_ictv_backup.json"

# Load JSON data
with open(json_path, "r") as f:
    data = json.load(f)

# Create a set of accessions we've already processed
processed_accessions = {
    entry.get("Virus_GENBANK_accession")
    for entry in data
    if entry and "Taxon_ID" in entry
}

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

# Enrich entries with Taxon_ID where possible and write progress
for i, entry in enumerate(data):
    if not entry:
        print(f"Skipping null entry at index {i}")
        continue

    accession = entry.get("Virus_GENBANK_accession")
    if not accession or not isinstance(accession, str):
        print(f"Skipping entry with invalid accession at index {i}")
        continue

    accession = accession.strip()

    if accession in processed_accessions:
        print(f"Already processed {accession}, skipping.")
        continue

    taxid = get_taxon_id_from_accession(accession)
    if taxid is not None:
        entry["Taxon_ID"] = taxid
        print(f"{accession} => Taxon ID: {taxid}")
    else:
        print(f"{accession} => No Taxon ID found, skipping.")

    # Write after every update to ensure resuming is possible
    with open(json_path, "w") as f:
        json.dump(data, f, indent=2)
    with open(backup_path, "w") as f:
        json.dump(data, f, indent=2)

    sleep(0.34)  # NCBI limit
