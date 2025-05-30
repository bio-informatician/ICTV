import pandas as pd
import os
import re

input_file = 'MSL.xlsx'
sheet_name = 'MSL'
output_dir = 'converted_files'

os.makedirs(output_dir, exist_ok=True)

# Load the specific sheet
df = pd.read_excel(input_file, sheet_name=sheet_name)

# Drop unwanted index column if present
if 'Unnamed: 0' in df.columns:
    df = df.drop(columns=['Unnamed: 0'])

# Clean column names for XML tags
df.columns = [re.sub(r'\W+', '_', str(col)) for col in df.columns]

# Export to TSV
df.to_csv(os.path.join(output_dir, f'{sheet_name}.tsv'), sep='\t', index=False)

# Export to JSON
df.to_json(os.path.join(output_dir, f'{sheet_name}.json'), orient='records', indent=2)

# Export to XML
df.to_xml(os.path.join(output_dir, f'{sheet_name}.xml'), index=False)

print("Conversion completed successfully.")
