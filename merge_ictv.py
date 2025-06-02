import json
import os

INPUT_FOLDER = "converted_files"
MSL_FILE = os.path.join(INPUT_FOLDER, "MSL.json")
VMR_FILE = os.path.join(INPUT_FOLDER, "VMR.json")
OUTPUT_FILE = os.path.join(INPUT_FOLDER, "merged_ictv.json")

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

    msl_index = {entry["ICTV_ID"]: entry for entry in msl_data if "ICTV_ID" in entry}
    vmr_index = {entry["ICTV_ID"]: entry for entry in vmr_data if "ICTV_ID" in entry}

    merged_data = []
    for ictv_id, msl_entry in msl_index.items():
        if ictv_id in vmr_index:
            merged = merge_entries(msl_entry, vmr_index[ictv_id])
            merged_data.append(merged)

    print(f"Writing merged data to {OUTPUT_FILE} ...")
    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        json.dump(merged_data, f, indent=2)

    print(f"Merged {len(merged_data)} entries. Saved to '{OUTPUT_FILE}'.")

if __name__ == "__main__":
    main()
