import json
import os
from collections import defaultdict

INPUT_FOLDER = "converted_files"
MSL_FILE = os.path.join(INPUT_FOLDER, "MSL.json")
VMR_FILE = os.path.join(INPUT_FOLDER, "VMR.json")
OUTPUT_FILE = os.path.join(INPUT_FOLDER, "merged_ictv.json")

def merge_vmr_entries(vmr_entries):
    merged = {}
    if not vmr_entries:
        return merged

    keys = vmr_entries[0].keys()
    for key in keys:
        values = [entry[key] for entry in vmr_entries if key in entry and entry[key]]

        # Deduplicate and concatenate string values with ";"
        if all(isinstance(v, str) for v in values):
            unique_vals = sorted(set(values))
            merged[key] = "; ".join(unique_vals)
        else:
            # Assume other types can be overridden or handled case-by-case
            merged[key] = values[0]

    return merged

def merge_entries(msl_entry, vmr_entry):
    merged = {}
    merged.update(vmr_entry)
    merged.update(msl_entry)  # MSL overwrites duplicates
    return merged

def main():
    if not os.path.exists(INPUT_FOLDER):
        raise FileNotFoundError(f"Input folder '{INPUT_FOLDER}' does not exist")

    print(f"Loading {MSL_FILE} ...")
    with open(MSL_FILE, "r", encoding="utf-8") as f:
        msl_data = json.load(f)

    print(f"Loading {VMR_FILE} ...")
    with open(VMR_FILE, "r", encoding="utf-8") as f:
        vmr_data = json.load(f)

    # Group VMR entries by "Species_Sort"
    vmr_grouped = defaultdict(list)
    for entry in vmr_data:
        key = entry.get("Species_Sort")
        if key:
            vmr_grouped[key].append(entry)

    # Merge VMR entries by "Species_Sort"
    merged_vmr = {key: merge_vmr_entries(entries) for key, entries in vmr_grouped.items()}

    # Index MSL by "Sort"
    msl_index = {entry["Sort"]: entry for entry in msl_data if "Sort" in entry}

    # Merge VMR and MSL based on matching "Species_Sort" and "Sort"
    merged_data = []
    for sort_key, msl_entry in msl_index.items():
        if sort_key in merged_vmr:
            merged = merge_entries(msl_entry, merged_vmr[sort_key])
            merged_data.append(merged)

    print(f"Writing merged data to {OUTPUT_FILE} ...")
    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        json.dump(merged_data, f, indent=2)

    print(f"Merged {len(merged_data)} entries. Saved to '{OUTPUT_FILE}'.")

if __name__ == "__main__":
    main()
