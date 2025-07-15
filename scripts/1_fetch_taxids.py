import json
import time
import requests
from tqdm import tqdm
import re
import openpyxl  # pip install openpyxl

EMAIL = "shahram4dev@gmail.com"
API_KEY = "a1a61a96906ca0d589efef9e91541019b808"
BATCH_SIZE = 100
BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

def fetch_taxid_via_efetch(accession):
    try:
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=gb&retmode=text"
        response = requests.get(url)
        response.raise_for_status()
        text = response.text
        match = re.search(r'/db_xref="taxon:(\d+)"', text)
        if match:
            return match.group(1)
    except Exception as e:
        print(f"❌ Error fetching TaxID via efetch for {accession}: {e}")
    return None

def extract_accessions(entry):
    val = entry.get("Virus_GENBANK_accession")
    if not val:
        return []
    parts = [v.strip() for v in val.split(';')]
    accs = []
    for part in parts:
        if ':' in part:
            accs.append(part.split(':', 1)[1].strip())
        else:
            accs.append(part)
    return accs

def get_all_accessions_from_json(file_path):
    with open(file_path, "r") as f:
        data = json.load(f)

    acc_map = {}
    if isinstance(data, list):
        for entry in data:
            if "TaxID" in entry and entry["TaxID"]:
                continue
            for acc in extract_accessions(entry):
                acc_map[acc] = entry
    elif isinstance(data, dict):
        for entry in data.values():
            if isinstance(entry, dict):
                if "TaxID" in entry and entry["TaxID"]:
                    continue
                for acc in extract_accessions(entry):
                    acc_map[acc] = entry
    return list(acc_map.keys()), acc_map, data

def fetch_summary_by_accession(acc):
    try:
        response = requests.get(
            BASE_URL + "esummary.fcgi",
            params={
                "db": "nuccore",
                "id": acc,
                "retmode": "json",
                "email": EMAIL,
                "api_key": API_KEY
            }
        )
        response.raise_for_status()
        data = response.json()
        result = data.get("result", {})
        uid = next((k for k in result if k != "uids"), None)
        if uid:
            entry = result.get(uid)
            if entry and "taxid" in entry:
                return entry["taxid"]
    except Exception as e:
        print(f"❌ Error fetching summary for {acc}: {e}")
    return None

def fetch_and_update_taxids(accessions, acc_map, json_data, json_path):
    failed = []

    for i in tqdm(range(0, len(accessions), BATCH_SIZE), desc="Fetching TaxIDs"):
        batch = accessions[i:i + BATCH_SIZE]
        taxid_map = {}

        for acc in batch:
            taxid = fetch_summary_by_accession(acc)
            if taxid:
                taxid_map[acc] = taxid
            else:
                failed.append(acc)

        # Update original data with new TaxIDs
        for acc, taxid in taxid_map.items():
            entry = acc_map.get(acc)
            if entry and "TaxID" not in entry:
                entry["TaxID"] = taxid

        with open(json_path, "w") as f:
            json.dump(json_data, f, indent=2)

        time.sleep(0.1 if API_KEY else 0.34)

    if failed:
        print("\n⚠️ Failed to retrieve TaxID for:")
        for acc in failed:
            print(f"  - {acc}")
        print("\n🔄 Retrying failed accessions with efetch method...")
        retry_failed = []
        for acc in failed:
            taxid = fetch_taxid_via_efetch(acc)
            if taxid:
                entry = acc_map.get(acc)
                if entry and "TaxID" not in entry:
                    entry["TaxID"] = taxid
            else:
                retry_failed.append(acc)

        # Save updates after retry
        with open(json_path, "w") as f:
            json.dump(json_data, f, indent=2)

        if retry_failed:
            print("\n⚠️ Still failed to retrieve TaxID for:")
            for acc in retry_failed:
                print(f"  - {acc}")
        else:
            print("✅ All previously failed accessions processed successfully on retry.")
    else:
        print("✅ All accessions processed successfully.")

def extract_accessions_from_excel(excel_path, sheet_name=None):
    """
    Extract accession IDs and their full URLs from hyperlink URLs in Excel file.
    Returns dict: { accession: full_url }
    """
    wb = openpyxl.load_workbook(excel_path, data_only=True)
    ws = wb[sheet_name] if sheet_name else wb.active

    accession_url_map = {}

    for row in ws.iter_rows():
        for cell in row:
            if cell.hyperlink:
                url = cell.hyperlink.target
                if url:
                    parts = url.strip().split('/')
                    if parts:
                        last_part = parts[-1]
                        for acc in last_part.split(','):
                            acc = acc.strip()
                            if acc:
                                accession_url_map[acc] = url
    return accession_url_map


if __name__ == "__main__":
    input_json_path = "converted_files/merged_ictv.json"
    input_excel_path = "your_excel_file.xlsx"  # <-- change this to your Excel filename

    # Extract accessions from JSON
    accessions_json, acc_map, json_data = get_all_accessions_from_json(input_json_path)
    print(f"🔍 Found {len(accessions_json)} accessions in JSON.")

    # Extract accessions and URLs from Excel
    accession_url_map = extract_accessions_from_excel(input_excel_path)
    print(f"🔍 Found {len(accession_url_map)} accessions in Excel hyperlinks.")

    # Combine unique accessions
    all_accessions = list(set(accessions_json) | set(accession_url_map.keys()))
    print(f"🔍 Total unique accessions to process: {len(all_accessions)}")

    # Add or update "from_url" field for accessions from Excel
    for acc, url in accession_url_map.items():
        if acc in acc_map:
            if acc_map[acc].get("from_url") != url:
                acc_map[acc]["from_url"] = url
        else:
            acc_map[acc] = {"from_url": url}

    fetch_and_update_taxids(all_accessions, acc_map, json_data, input_json_path)
    print(f"✅ JSON file updated: {input_json_path}")
