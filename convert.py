import pandas as pd
import os

input_file = 'VMR.xlsx'
output_dir = 'converted_files'

os.makedirs(output_dir, exist_ok=True)

# Load Excel file (all sheets)
xls = pd.ExcelFile(input_file)

for sheet_name in xls.sheet_names:
    df = pd.read_excel(xls, sheet_name=sheet_name)

    # Save TSV
    tsv_path = os.path.join(output_dir, f'{sheet_name}.tsv')
    df.to_csv(tsv_path, sep='\t', index=False)

    # Save JSON
    json_path = os.path.join(output_dir, f'{sheet_name}.json')
    df.to_json(json_path, orient='records', indent=2)

    # Save XML (Pandas 1.3+ has to_xml)
    xml_path = os.path.join(output_dir, f'{sheet_name}.xml')
    df.to_xml(xml_path, index=False)

print("Conversion completed!")
