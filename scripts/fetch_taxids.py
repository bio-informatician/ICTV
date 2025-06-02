import json
from Bio import Entrez

Entrez.email = "your.email@example.com"  # replace with your email

def get_all_accessions_from_json(file_path, key="Virus_GENBANK_accession"):
    def extract_accessions(val):
        result = []
        if isinstance(val, list):
            values = val
        else:
            values = [val]
        for v in values:
            if not v:
                continue
            if ";" in v or ":" in v:
                parts = v.split(";")
                for part in parts:
                    if ":" in part:
                        acc = part.strip().split(":")[-1].strip()
                        result.append(acc)
            else:
                result.append(v.strip())
        return result

    with open(file_path, "r") as f:
        data = json.load(f)

    accessions = []
    if isinstance(data, list):
        for entry in data:
            val = entry.get(key)
            if val is not None:
                accessions.extend(extract_accessions(val))
    elif isinstance(data, dict):
        for k, entry in data.items():
            if isinstance(entry, dict):
                val = entry.get(key)
                if val is not None:
                    accessions.extend(extract_accessions(val))

    unique_accessions = list(dict.fromkeys(accessions))
    return unique_accessions, data


def fetch_taxid_for_accessions(accessions):
    import json

    taxid_map = {}
    batch_size = 50
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i+batch_size]
        try:
            handle = Entrez.esummary(db="nuccore", id=",".join(batch), retmode="json")
            response_text = handle.read()
            handle.close()

            records = json.loads(response_text)
            summaries = records.get("result", {})
            uids = summaries.get("uids", [])
            for uid in uids:
                summary = summaries.get(uid, {})
                acc = summary.get("caption") or summary.get("Caption")
                taxid = summary.get("taxid") or summary.get("TaxId")
                if acc and taxid:
                    taxid_map[acc] = taxid

        except Exception as e:
            print(f"Error fetching taxid for batch {batch}: {e}")

    return taxid_map


def update_json_with_taxid(data, taxid_map, accession_key="Virus_GENBANK_accession", taxid_key="TaxID"):
    def get_taxids_for_accessions(accessions):
        if isinstance(accessions, list):
            taxids = [taxid_map.get(acc) for acc in accessions if acc in taxid_map]
            taxids = list(dict.fromkeys(filter(None, taxids)))
            if len(taxids) == 1:
                return taxids[0]
            elif len(taxids) > 1:
                return taxids
            else:
                return None
        else:
            return taxid_map.get(accessions)

    if isinstance(data, list):
        for entry in data:
            if accession_key in entry:
                taxid_value = get_taxids_for_accessions(entry[accession_key])
                if taxid_value is not None:
                    entry[taxid_key] = taxid_value

    elif isinstance(data, dict):
        for k, entry in data.items():
            if isinstance(entry, dict) and accession_key in entry:
                taxid_value = get_taxids_for_accessions(entry[accession_key])
                if taxid_value is not None:
                    entry[taxid_key] = taxid_value
    return data


if __name__ == "__main__":
    input_json_path = "converted_files/merged_ictv.json"

    # Step 1: Extract all accession numbers and load JSON data
    accessions, json_data = get_all_accessions_from_json(input_json_path)
    print(f"Got {len(accessions)} total accessions from JSON.")

    # Step 2: Fetch TaxIDs for all accessions (in batches)
    accession_taxid_map = fetch_taxid_for_accessions(accessions)
    print(f"Fetched TaxIDs for {len(accession_taxid_map)} accessions.")

    # Step 3: Update JSON with TaxID info
    updated_json = update_json_with_taxid(json_data, accession_taxid_map)

    # Step 4: Save updated JSON (overwrite original)
    with open(input_json_path, "w") as f_out:
        json.dump(updated_json, f_out, indent=2)

    print(f"Finished updating JSON file: {input_json_path}")
