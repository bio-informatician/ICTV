import pandas as pd
import os

INPUT_XLSX = "ICTV_Master_Species_List_2024_MSL40.v1.xlsx"
OUTPUT_DIR = "raw_files"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load the Excel file
xlsx = pd.ExcelFile(INPUT_XLSX)

# Export the full "MSL" sheet with all content (no rows or columns skipped)
if "MSL" in xlsx.sheet_names:
    df = pd.read_excel(xlsx, sheet_name="MSL")  # No skiprows
    tsv_path = os.path.join(OUTPUT_DIR, "MSL_full_export.tsv")
    df.to_csv(tsv_path, sep="\t", index=False)
    print("MSL sheet saved as MSL_full_export.tsv")
else:
    print("No such sheet: MSL")
