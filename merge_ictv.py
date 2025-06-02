import requests
import json
from collections import OrderedDict

# URLs of the JSON files
msl_url = "https://raw.githubusercontent.com/bio-informatician/ICTV/refs/heads/main/converted_files/MSL.json"
vmr_url = "https://raw.githubusercontent.com/bio-informatician/ICTV/refs/heads/main/converted_files/VMR.json"

# Download JSON data
def download_json(url):
    response = requests.get(url)
    response.raise_for_status()
    return response.json()

# Merge two entries by ICTV_ID, removing duplicate keys (MSL fields take precedence)
def merge_entries(msl_entry, vmr_entry):
    merged = {}
    merged.update(vmr_entry)
    merged.update(msl_entry)  # MSL overwrites duplicate keys
    return merged

# Main script
def main():
    print("Downloading MSL...")
    msl_data = download_json(msl_url)
    print("Downloading VMR...")
    vmr_data = download_json(vmr_url)

    # Index entries by ICTV_ID
    msl_index = {entry["ICTV_ID"]: entry for entry in msl_data if "ICTV_ID" in entry}
    vmr_index = {entry["ICTV_ID"]: entry for entry in vmr_data if "ICTV_ID" in entry}

    # Merge where ICTV_ID matches
    merged_data = []
    for ictv_id, msl_entry in msl_index.items():
        if ictv_id in vmr_index:
            merged = merge_entries(msl_entry, vmr_index[ictv_id])
            merged_data.append(merged)

    # Output result
    with open("merged_ictv.json", "w") as f:
        json.dump(merged_data, f, indent=2)

    print(f"Merged {len(merged_data)} entries. Saved to 'converted_files/merged_ictv.json'.")

if __name__ == "__main__":
    main()
