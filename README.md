[![Convert MSL Excel](https://github.com/bio-informatician/ICTV/actions/workflows/convert-msl.yml/badge.svg)](https://github.com/bio-informatician/ICTV/actions/workflows/convert-msl.yml)   [![Convert VMR Excel](https://github.com/bio-informatician/ICTV/actions/workflows/convert-vmr.yml/badge.svg)](https://github.com/bio-informatician/ICTV/actions/workflows/convert-vmr.yml)   

[![Fetch Taxon IDs](https://github.com/bio-informatician/ICTV/actions/workflows/fetch_taxids.yml/badge.svg)](https://github.com/bio-informatician/ICTV/actions/workflows/fetch_taxids.yml)
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

## 🚀 How to Use

1. **Clone the repo** or edit `VMR.xlsx` directly in GitHub.
2. Commit the updated Excel file.
3. GitHub Actions will run automatically.
4. Check the `converted_files/` directory for the updated TSV, JSON, and XML files.

---

## 📚 Requirements (for local development)

If you want to run the script locally:

### Install Python packages:

```bash
pip install -r requirements.txt
````

### Requirements include:

* `pandas`
* `openpyxl`
* `lxml`

---

## ✍️ Author & License

This project is maintained by the bioinformatics team. Contributions and improvements are welcome.

> Licensed under MIT. See `LICENSE` file for more information.

---

## 🧠 Notes

* If your conversion fails due to XML tag naming errors, the script will sanitize column names to valid XML-safe formats automatically.
* If pushing fails during automation, ensure Actions bot has appropriate permissions or consider using GitHub artifacts instead of committing outputs.
