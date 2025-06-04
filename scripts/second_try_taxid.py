import json
import requests
import time

INPUT_FILE = "converted_files/merged_ictv.json"
API_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy"
HEADERS = {
    'accept': 'application/json',
    'content-type': 'application/json'
}

taxonomy_cache = {}

def query_taxonomy(taxon):
    taxon_str = str(taxon).strip()
    if taxon_str in taxonomy_cache:
        return taxonomy_cache[taxon_str]
    try:
        response = requests.post(API_URL, headers=HEADERS, json={"taxons": [taxon_str]})
        if response.status_code == 200:
            data = response.json()
            taxonomy_cache[taxon_str] = data
            return data
    except Exception:
        pass
    taxonomy_cache[taxon_str] = None
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
        taxonomy = extract_taxonomy_object(resp)
        if taxonomy and taxonomy.get("rank") == "SPECIES":
            return taxonomy.get("tax_id"), taxonomy
        time.sleep(0.1)
    return None, None

def split_values(value):
    if not value:
        return []
    return [v.strip() for v in value.split(";") if v.strip()]

def main():
    with open(INPUT_FILE, 'r', encoding='utf-8') as f:
        entries = json.load(f)

    for entry in entries:
        species_tax_id = None
        species_taxonomy = None
        response_from_taxid = None
        response_from_name = None

        taxid = entry.get("TaxID")
        species_name = entry.get("Species")
        virus_names = split_values(entry.get("Virus_name_s_"))
        genbank_accessions = split_values(entry.get("Virus_GENBANK_accession"))

        # Step 1: Try TaxID
        if taxid:
            response_from_taxid = query_taxonomy(taxid)
            species_tax_id = extract_species_tax_id(response_from_taxid)
            if species_tax_id:
                species_taxonomy = extract_taxonomy_object(query_taxonomy(species_tax_id))

        # Step 2: Try Species name
        if not species_tax_id and species_name:
            response_from_name = query_taxonomy(species_name)
            species_tax_id = extract_species_tax_id(response_from_name)
            if species_tax_id:
                species_taxonomy = extract_taxonomy_object(query_taxonomy(species_tax_id))

        # Step 3: Try querying each virus name directly
        if not species_tax_id:
            for virus_name in virus_names:
                response_from_virus = query_taxonomy(virus_name)
                matched_tax_id = extract_tax_id_if_name_matches(response_from_virus, virus_name)
                if matched_tax_id:
                    matched_tax = extract_taxonomy_object(response_from_virus)
                    if matched_tax and matched_tax.get("rank") == "SPECIES":
                        species_tax_id = matched_tax.get("tax_id")
                        species_taxonomy = matched_tax
                        break  # Stop on first valid match

        # Step 4: Use lineage from TaxID response
        if not species_tax_id and response_from_taxid:
            lineage = extract_lineage(response_from_taxid)
            if lineage:
                species_tax_id, species_taxonomy = search_species_in_lineage(lineage)

        if species_tax_id and species_taxonomy and species_taxonomy.get("rank") == "SPECIES":
            entry["Species_Tax_ID"] = species_tax_id
        else:
            print(f"❌ {species_name or taxid or virus_names}")

    with open(INPUT_FILE, 'w', encoding='utf-8') as f:
        json.dump(entries, f, indent=2)

if __name__ == '__main__':
    main()

