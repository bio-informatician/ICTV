import pandas as pd
import os

# Input Excel file name
input_file = "VMR.xlsx"

# Output directories
output_dirs = {
    "json": "json",
    "csv": "csv",
    "tsv": "tsv",
    "xml": "xml"
}

# Create output directories if they don't exist
for dir_path in output_dirs.values():
    os.makedirs(dir_path, exist_ok=True)

# Read the Excel file
try:
    df = pd.read_excel(input_file, engine='openpyxl')
except Exception as e:
    print(f"Error reading Excel file: {e}")
    exit(1)

# Export to JSON
json_path = os.path.join(output_dirs["json"], "output.json")
df.to_json(json_path, orient="records", indent=2)

# Export to CSV
csv_path = os.path.join(output_dirs["csv"], "output.csv")
df.to_csv(csv_path, index=False)

# Export to TSV
tsv_path = os.path.join(output_dirs["tsv"], "output.tsv")
df.to_csv(tsv_path, sep="\t", index=False)

# Export to XML
xml_path = os.path.join(output_dirs["xml"], "output.xml")
try:
    df.to_xml(xml_path, index=False)
except AttributeError:
    print("Pandas version may be too old for .to_xml(). Update pandas if needed.")
    exit(1)

print("Conversion completed successfully.")
