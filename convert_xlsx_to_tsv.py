import pandas as pd
import os

INPUT_XLSX = "VMR.xlsx"  # Replace with your actual XLSX file name
OUTPUT_DIR = "raw_files"  # This should match the INPUT_DIR in your existing script

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Define a list where each item specifies a sheet and the filename for output
WANTED_FILES = [
    {"sheet_name": "VMR MSL40", "filename": "VMR"},
]

# Load the Excel file
xlsx = pd.ExcelFile(INPUT_XLSX)

# Convert each sheet to a TSV file with all columns
for wanted_file in WANTED_FILES:
    print(f"Processing: {wanted_file}")
    if wanted_file["sheet_name"] in xlsx.sheet_names:
        df = pd.read_excel(xlsx, sheet_name=wanted_file["sheet_name"], skiprows=3)
        tsv_path = os.path.join(OUTPUT_DIR, f"{wanted_file['filename']}.tsv")
        df.to_csv(tsv_path, sep="\t", index=False)
    else:
        print("No such sheet: ", wanted_file["sheet_name"])

# Export full sheet as "VJDB_catalogue.tsv"
if "VMR MSL40" in xlsx.sheet_names:
    full_df = pd.read_excel(xlsx, sheet_name="VMR MSL40", skiprows=3)
    full_tsv_path = os.path.join(OUTPUT_DIR, "VMR_catalogue.tsv")
    full_df.to_csv(full_tsv_path, sep="\t", index=False)
    print("Full sheet saved as VMR_catalogue.tsv")
else:
    print("No such sheet: VMR MSL40")

print("TSV files have been generated!")
