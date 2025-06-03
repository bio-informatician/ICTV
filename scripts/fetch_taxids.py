import json
import requests
import time

INPUT_FILE = 'merged_ictv.json'
API_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy"
HEADERS = {
    'accept': 'application/json',
    'content-type': 'application/json'
}

def query_taxonomy(taxon):
    try:
        response = requests.post(API_URL, headers=HEADERS, json={"taxons": [str(taxon)]})
        if response.status_code == 200:
            return response.json()
    except Exception:
        pass
    return None

def extract_taxonomy_object(api_response):
    if not api_response:
        return None
    try:
        return api_response.get("taxonomy_nodes", [{}])[0].get("taxonomy", None)
    except (IndexError, AttributeError):
        return None

def extract_species_tax_id(api_response):
    taxonomy = extract_taxonomy_object(api_response)
    if taxonomy and taxonomy.get("rank") == "SPECIES":
        return taxonomy.get("tax_id")
    return None

def extract_tax_id_if_name_matches(api_response, virus_name):
    taxonomy = extract_taxonomy_object(api_response)
    if taxonomy and taxonomy.get("organism_name", "").strip() == virus_name.strip():
        return taxonomy.get("tax_id")
    return None

def extract_lineage(api_response):
    taxonomy = extract_taxonomy_object(api_response)
    if taxonomy:
        return taxonomy.get("lineage", [])
    return []

def search_species_in_lineage(lineage):
    for tax_id in reversed(lineage):
        resp = query_taxonomy(tax_id)
        if resp:
            species_tax_id = extract_species_tax_id(resp)
            if species_tax_id:
                return species_tax_id
        time.sleep(0.5)
    return None

def main():
    with open(INPUT_FILE, 'r', encoding='utf-8') as f:
        entries = json.load(f)

    for entry in entries:
        species_tax_id = None
        response_from_taxid = None
        response_from_name = None

        taxid = entry.get("TaxID")
        species_name = entry.get("Species")
        virus_name = entry.get("Virus_name_s_")

        # Step 1: Try TaxID
        if taxid:
            response_from_taxid = query_taxonomy(taxid)
            species_tax_id = extract_species_tax_id(response_from_taxid)

        # Step 2: Try Species name
        if not species_tax_id and species_name:
            response_from_name = query_taxonomy(species_name)
            species_tax_id = extract_species_tax_id(response_from_name)

        # Step 3: Match Virus_name_s_
        if not species_tax_id and response_from_name and virus_name:
            species_tax_id = extract_tax_id_if_name_matches(response_from_name, virus_name)

        # Step 4: Use lineage from TaxID response
        if not species_tax_id and response_from_taxid:
            lineage = extract_lineage(response_from_taxid)
            if lineage:
                species_tax_id = search_species_in_lineage(lineage)

        if species_tax_id:
            entry["Species_Tax_ID"] = species_tax_id
        else:
            print(f"❌ Could not determine Species_Tax_ID for entry: {species_name or taxid}")

        # Live update: write updated JSON back to input file after each entry processed
        with open(INPUT_FILE, 'w', encoding='utf-8') as f:
            json.dump(entries, f, indent=2)

        time.sleep(0.5)

if __name__ == '__main__':
    main()
