import json
import os
from time import sleep
from Bio import Entrez

Entrez.email = os.getenv("ENTREZ_EMAIL", "fallback@example.com")

json_path = "converted_files/merged_ictv.json"
with open(json_path, "r") as f:
    data = json.load(f)

def get_taxon_id_from_accession(accession):
    try:
        handle = Entrez.esearch(db="nuccore", term=accession)
        record = Entrez.read(handle)
        handle.close()

        if not record["IdList"]:
            raise ValueError(f"No NCBI records found for accession: {accession}")

        uid = record["IdList"][0]
        handle = Entrez.esummary(db="nuccore", id=uid)
        summary = Entrez.read(handle)
        handle.close()

        taxid = summary[0].get("TaxId")
        if not taxid:
            raise ValueError(f"TaxId not found in summary for accession: {accession}")
        return int(taxid)
    except Exception as e:
        print(f"Error fetching Taxon ID for {accession}: {e}")
        raise  # Re-raise to stop execution

for entry in data:
    accession = entry.get("Virus_GENBANK_accession", "").strip()
    if accession:
        taxid = get_taxon_id_from_accession(accession)
        entry["Taxon_ID"] = taxid
        print(f"{accession} => Taxon ID: {taxid}")
        sleep(0.34)  # respect NCBI rate limits

with open(json_path, "w") as f:
    json.dump(data, f, indent=2)
