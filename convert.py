import pandas as pd
import os
import re

input_file = 'VMR.xlsx'
sheet_name = 'VMR MSL40'
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
df.to_csv(os.path.join(output_dir, f'{input_file}.tsv'), sep='\t', index=False)

# Export to JSON
df.to_json(os.path.join(output_dir, f'{input_file}.json'), orient='records', indent=2)

# Export to XML
df.to_xml(os.path.join(output_dir, f'{input_file}.xml'), index=False)

print("Conversion completed successfully.")
