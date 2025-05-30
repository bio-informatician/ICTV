import pandas as pd
import os

input_file = 'VMR.xlsx'
sheet_name = 'VMR MSL40'
output_dir = 'converted_files'

os.makedirs(output_dir, exist_ok=True)

# Read only the specified sheet
df = pd.read_excel(input_file, sheet_name=sheet_name)

# Drop common unwanted index column if present
if 'Unnamed: 0' in df.columns:
    df = df.drop(columns=['Unnamed: 0'])

# Replace invalid XML tag characters in column names
df.columns = [str(col).replace(' ', '_').replace(':', '_') for col in df.columns]

# Export to TSV
df.to_csv(os.path.join(output_dir, f'{sheet_name}.tsv'), sep='\t', index=False)

# Export to JSON
df.to_json(os.path.join(output_dir, f'{sheet_name}.json'), orient='records', indent=2)

# Export to XML
df.to_xml(os.path.join(output_dir, f'{sheet_name}.xml'), index=False)

print("Conversion completed successfully.")
