import pandas as pd
import os

INPUT_XLSX = "VMR_MSL40.v1.20250307.xlsx"
OUTPUT_DIR = "raw_files"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load the Excel file
xlsx = pd.ExcelFile(INPUT_XLSX)

# Export the full "VMR" sheet with all content (no rows or columns skipped)
if "MSL" in xlsx.sheet_names:
    df = pd.read_excel(xlsx, sheet_name="VMR MSL40")  # No skiprows
    tsv_path = os.path.join(OUTPUT_DIR, "VMR_full_export.tsv")
    df.to_csv(tsv_path, sep="\t", index=False)
    print("VMR sheet saved as VMR_full_export.tsv")
else:
    print("No such sheet: VMR")
