import pandas as pd
import json
import os
import warnings

warnings.filterwarnings(
    "ignore",
    category=UserWarning,
    module="openpyxl"
)

# Configuration: Map the file and tab to the EXACT column header names
# The script will look for these strings in the first row.
FILE_CONFIG = {
    "VMR_MSL41.v1.20260320.xlsx": {"sheet": "VMR MSL41", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession"]},
    "VMR_MSL40.v2.20260223.xlsx": {"sheet": "VMR MSL40", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession"]},
    "VMR_MSL40.v2.20251013.xlsx": {"sheet": "VMR MSL40", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession"]},
    "VMR_MSL40.v1.20250307.xlsx": {"sheet": "VMR MSL40", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession"]},
    "VMR_MSL39.v4_20241106.xlsx": {"sheet": "VMR MSL39", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession"]},
    "VMR_MSL39.v2_20240920.xlsx": {"sheet": "VMR MSL39", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_MSL39.v1_20240912.xlsx": {"sheet": "VMR MSL39 v1", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_MSL38_v3.xlsx": {"sheet": "VMR MSL38 v3", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_MSL38_v2.xlsx": {"sheet": "VMR MSL38 v2", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_MSL38_v1.xlsx": {"sheet": "VMR MSL38 v1", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_21-221122_MSL37.xlsx": {"sheet": "VMRb37", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_20-190822_MSL37.3.xlsx": {"sheet": "VMRb37", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession"]},
    "VMR_20-190822_MSL37.2.xlsx": {"sheet": "VMRb37", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession"]},
    "VMR_20-190822_MSL37.1.xlsx": {"sheet": "VMRb37", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession"]},
    "VMR_19-250422_MSL37.xlsx": {"sheet": "VMR 19 MSL37", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession"]},
    "VMR_18-191021_MSL36.xlsx": {"sheet": "VMR 18 MSL36", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_17-200721_MSL36.xlsx": {"sheet": "VMRb36", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_16-180521_MSL36.xlsx": {"sheet": "VMRb36", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_15-010820_MSL35.xlsx": {"sheet": "VMRb35 010820 MSL35", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_14-010520_MSL35.xlsx": {"sheet": "VMRb35", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_13-271119_MSL34.xlsx": {"sheet": "VMR 271119 MSL34", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_12-160919.2_MSL34.xlsx": {"sheet": "VMR 160919 MSL34", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_11-010619_MSL34.xlsx": {"sheet": "VMRb34", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_10-190519_MSL34.xlsx": {"sheet": "VMRb34", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_09-220319_MSL34.xlsx": {"sheet": "VMRb34v1", "cols": ["Exemplar or additional isolate", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_08-201218_MSL33.xlsx": {"sheet": "VMRb33", "cols": ["Exemplar (E) or additional isolate (A)", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_07-200718_MSL32.xlsx": {"sheet": "VMR", "cols": ["Exemplar or additional isolate ", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_06-110518.MSL32.xlsx": {"sheet": "VMR 210318", "cols": ["Exemplar or Additional ", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_05-210318.MSL32.xlsx": {"sheet": "VMR 210318", "cols": ["Exemplar ", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_04-290118_MSL31.xlsx": {"sheet": "VMR 081117", "cols": ["Exemplar ", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_03-030118_MSL31.xlsx": {"sheet": "VMR 081117", "cols": ["Exemplar ", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_02-081117_MSL31.xlsx": {"sheet": "VMR 081117", "cols": ["Exemplar ", "Virus GENBANK accession", "Virus REFSEQ accession"]},
    "VMR_01-010817_MSL31.xlsx": {"sheet": "VMR 010817", "cols": ["Exemplar ", "Virus GENBANK accession", "Virus REFSEQ accession"]},
}

def process_files(input_dir, output_dir):
    if not os.path.exists(output_dir): os.makedirs(output_dir)

    for filename, config in FILE_CONFIG.items():
        file_path = os.path.join(input_dir, filename)
        if not os.path.exists(file_path): continue

        try:
            # Load sheet using header=0 (first row)
            df = pd.read_excel(file_path, sheet_name=config['sheet'], header=0)
            
            # Clean column names (strip spaces to ensure matching)
            df.columns = df.columns.astype(str).str.strip()
            
            output_rows = []
            for _, row in df.iterrows():
                entry = {}
                for col_name in config['cols']:
                    clean_col = col_name.strip()

                    if clean_col in df.columns:
                        val = row[clean_col]

                        if pd.notna(val):
                            val_str = str(val).strip()

                            if val_str:
                                # Standardize ALL Exemplar column names
                                if "Exemplar" in clean_col:
                                    entry["Exemplar"] = val_str
                                else:
                                    entry[clean_col] = val_str
                
                # Check for the primary identifier (look for keys starting with 'Exemplar')
                if any("Exemplar" in k for k in entry.keys()):
                    output_rows.append(entry)

            out_name = os.path.splitext(filename)[0] + ".json"
            with open(os.path.join(output_dir, out_name), 'w', encoding='utf-8') as f:
                json.dump(output_rows, f, indent=4)
            print(f"Processed: {filename}")
        except Exception as e:
            print(f"Error processing {filename}: {e}")

if __name__ == "__main__":
    process_files("VMR_Files", "ncbi")