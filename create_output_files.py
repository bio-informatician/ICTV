import os
import pandas as pd
import json
import xml.etree.ElementTree as ET

INPUT_DIR = "raw_files"  # Folder containing TSV files

# Ensure top-level output directories exist
for folder in ["json", "csv", "tsv", "xml"]:
    os.makedirs(folder, exist_ok=True)


def dict_to_xml(tag, d):
    """Convert a dictionary to an XML element"""
    elem = ET.Element(tag)
    for key, val in d.items():
        child = ET.SubElement(elem, key.replace(" ", "_"))  # Replace spaces with underscores
        child.text = str(val)
    return elem


def create_xml_structure(data):
    """Convert flat data into XML format"""
    root = ET.Element("Root")

    for item in data:
        entry_elem = dict_to_xml("Entry", item)
        root.append(entry_elem)

    return ET.ElementTree(root)


for file in os.listdir(INPUT_DIR):
    if file.endswith(".tsv"):
        input_path = os.path.join(INPUT_DIR, file)
        base_name = os.path.splitext(file)[0]

        # Read TSV file
        df = pd.read_csv(input_path, sep='\t')

        # Handle NaN values before processing
        df = df.fillna('N/A')  # Replace NaN values with "N/A" or another placeholder

        # Ensure we have enough columns
        if df.shape[1] < 3:
            print(f"Skipping {file} (needs at least 3 columns)")
            continue

        # Check for the column containing the array (adjust column name accordingly)
        if 'vjdbv0.3_tags' in df.columns:
            # Convert the column containing the string representation of a list into an actual list
            df['vjdbv0.3_tags'] = df['vjdbv0.3_tags'].apply(lambda x: json.loads(x) if isinstance(x, str) else x)

        # Flatten the structure, no hierarchical nesting
        flat_data = []

        for _, row in df.iterrows():
            entry = row.to_dict()  # Flatten the row into a single dictionary
            flat_data.append(entry)

        # Save JSON file
        json_path = os.path.join("json", f"{base_name}.json")
        with open(json_path, "w", encoding="utf-8") as json_file:
            json.dump(flat_data, json_file, indent=4)

        # Save XML file
        xml_path = os.path.join("xml", f"{base_name}.xml")
        xml_tree = create_xml_structure(flat_data)
        xml_tree.write(xml_path, encoding="utf-8", xml_declaration=True)

        # Save CSV and TSV (no change needed as they are already flat)
        df.to_csv(os.path.join("csv", f"{base_name}.csv"), index=False)
        df.to_csv(os.path.join("tsv", f"{base_name}.tsv"), sep='\t', index=False)

print("Conversion complete!")








# import os
# import pandas as pd
# import json
# import xml.etree.ElementTree as ET

# INPUT_DIR = "raw_files"  # Folder containing TSV files

# # Ensure top-level output directories exist
# for folder in ["json", "csv", "tsv", "xml"]:
#     os.makedirs(folder, exist_ok=True)


# def dict_to_xml(tag, d):
#     #"""Convert a dictionary to an XML element"""
#     elem = ET.Element(tag)
#     for key, val in d.items():
#         child = ET.SubElement(elem, key.replace(" ", "_"))  # Replace spaces with underscores
#         child.text = str(val)
#     return elem


# def create_xml_structure(data):
#     #"""Convert hierarchical dictionary to XML format"""
#     root = ET.Element("Root")

#     for category, sub_categories in data.items():
#         category_elem = ET.SubElement(root, category)
#         for sub_category, items in sub_categories.items():
#             sub_category_elem = ET.SubElement(category_elem, sub_category)
#             for item in items:
#                 entry_elem = dict_to_xml("Entry", item)
#                 sub_category_elem.append(entry_elem)

#     return ET.ElementTree(root)


# for file in os.listdir(INPUT_DIR):
#     if file.endswith(".tsv"):
#         input_path = os.path.join(INPUT_DIR, file)
#         base_name = os.path.splitext(file)[0]

#         # Read TSV file
#         df = pd.read_csv(input_path, sep='\t')

#         # Handle NaN values before processing
#         df = df.fillna('N/A')  # Replace NaN values with "N/A" or another placeholder

#         # Ensure we have enough columns
#         if df.shape[1] < 3:
#             print(f"Skipping {file} (needs at least 3 columns)")
#             continue

#         # Extract column names
#         col1, col2 = df.columns[:2]  # First two columns for grouping
#         data_columns = df.columns[2:]  # Remaining columns for data

#         # Convert DataFrame into hierarchical JSON-like structure
#         json_structure = {}

#         for _, row in df.iterrows():
#             category = row[col1]
#             sub_category = row[col2]
#             entry = row[data_columns].to_dict()  # Convert remaining columns to dict

#             if category not in json_structure:
#                 json_structure[category] = {}
#             if sub_category not in json_structure[category]:
#                 json_structure[category][sub_category] = []

#             json_structure[category][sub_category].append(entry)

#         # Save JSON file
#         json_path = os.path.join("json", f"{base_name}.json")
#         with open(json_path, "w", encoding="utf-8") as json_file:
#             json.dump(json_structure, json_file, indent=4)

#         # Save XML file
#         xml_path = os.path.join("xml", f"{base_name}.xml")
#         xml_tree = create_xml_structure(json_structure)
#         xml_tree.write(xml_path, encoding="utf-8", xml_declaration=True)

#         # Save CSV and TSV
#         df.to_csv(os.path.join("csv", f"{base_name}.csv"), index=False)
#         df.to_csv(os.path.join("tsv", f"{base_name}.tsv"), sep='\t', index=False)

# print("Conversion complete!")
