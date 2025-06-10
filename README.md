[![Convert MSL Excel](https://github.com/bio-informatician/ICTV/actions/workflows/convert-msl.yml/badge.svg)](https://github.com/bio-informatician/ICTV/actions/workflows/convert-msl.yml) - [![Convert VMR Excel](https://github.com/bio-informatician/ICTV/actions/workflows/convert-vmr.yml/badge.svg)](https://github.com/bio-informatician/ICTV/actions/workflows/convert-vmr.yml)  -  [![Merge ICTV JSON files](https://github.com/bio-informatician/ICTV/actions/workflows/merge_ictv.yml/badge.svg)](https://github.com/bio-informatician/ICTV/actions/workflows/merge_ictv.yml)  -  [![Fetch Taxon IDs](https://github.com/bio-informatician/ICTV/actions/workflows/fetch_taxids.yml/badge.svg)](https://github.com/bio-informatician/ICTV/actions/workflows/fetch_taxids.yml)
# ICTV: Virus Metadata Conversion and Automation

This repository automates the conversion of virus metadata stored in `VMR.xlsx` into multiple machine-readable formats using GitHub Actions. When a new version of `VMR.xlsx` is committed, the repository automatically generates the following:

- TSV (Tab-Separated Values)
- JSON (JavaScript Object Notation)
- XML (Extensible Markup Language)

These files are output to a dedicated folder named `converted_files/`.

---

## 📦 File Overview

### 📂 Source

- `VMR.xlsx`: Source Excel file containing virus metadata
  - The conversion targets the sheet named **"VMR MSL40"**.

### 🧮 Script

- `convert.py`: Python script that performs the actual conversion of Excel data to `.tsv`, `.json`, and `.xml`.

### 📁 Output

- `converted_files/`: Directory containing converted files:
  - `VMR MSL40.tsv`
  - `VMR MSL40.json`
  - `VMR MSL40.xml`

---

## ⚙️ GitHub Actions Automation

This repository includes a GitHub Actions workflow defined in:

```

.github/workflows/convert-vmr.yml

````

### 🔁 Trigger

The workflow runs automatically when `VMR.xlsx` is updated (pushed to the repository).

### 🔨 Workflow Process

1. Checks out the repository.
2. Sets up Python 3.9 and installs dependencies.
3. Converts `VMR.xlsx` sheet `"VMR MSL40"` into TSV, JSON, and XML.
4. Commits and pushes the converted files back to the repository.

### ✅ Permissions

Ensure GitHub Actions has **write access** to the repository under:
> Settings → Actions → General → Workflow permissions → Enable "Read and write permissions"

---

## Steps 
1. merge ICTV files, based on the ICTV ID by combineing the informationo f Exampler and Addisionals (uf exists) into one record. [0_merge_ictv_files.py]
2. convert GenBank Accession ID to Taxonomy ID using NCBI resources. [1_fetch_taxids.py]
3. find the species taxon ID using the TaxID, species name or virus name. [2_species_taxid.py]
4. insert the taxonomy information using the species taxid into the database. [3_ictv2db.py]
