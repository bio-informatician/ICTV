"""Load VMR JSON and map its columns to the configured VJDB schema."""

import json
import re


# These are vjdbv1.1 schema field IDs, not the similarly named legacy CSV
# field IDs or their numeric field indices. The importer resolves each ID
# through the configured schema before writing.
VMR_COLUMN_FIELD_IDS = {
    "Realm": "realm_ictv",
    "Subrealm": "subrealm_ictv",
    "Kingdom": "kingdom_ictv",
    "Subkingdom": "subkingdom_ictv",
    "Phylum": "phylum_ictv",
    "Subphylum": "subphylum_ictv",
    "Class": "class_ictv",
    "Subclass": "subclass_ictv",
    "Order": "order_ictv",
    "Suborder": "suborder_ictv",
    "Family": "family_ictv",
    "Subfamily": "subfamily_ictv",
    "Genus": "genus_ictv",
    "Subgenus": "subgenus_ictv",
    "Species": "species_ictv",
    "ICTV_ID": "ictv_id",
    "Exemplar_or_additional_isolate": "is_ictv_exemplar",
    "Virus_name_s_": "virus_name_ictv",
    "Virus_name_abbreviation_s_": "abbreviation",
    "Virus_isolate_designation": "isolate",
    "Genome_coverage": "genome_completeness",
    "Host_source": "host_group_ictv",
}

ACCESSION_COLUMN = "Virus_GENBANK_accession"
_ACCESSION_RE = re.compile(
    r"(?<![A-Z0-9])([A-Z]{1,6}_?\d{5,9}(?:\.\d+)?)(?![A-Z0-9])",
    re.IGNORECASE,
)


def _is_empty(value):
    return value is None or (isinstance(value, str) and not value.strip()) or value == []


def _accessions(value):
    if not isinstance(value, str):
        return []
    return [match.group(1).upper() for match in _ACCESSION_RE.finditer(value)]


def _base_accession(accession):
    return re.sub(r"\.\d+$", "", accession.strip().upper())


def load_vmr_json(path):
    with open(path, "r", encoding="utf-8") as vmr_file:
        records = json.load(vmr_file)
    if not isinstance(records, list):
        raise ValueError(f"{path}: expected a JSON array of VMR record objects")
    for row_number, record in enumerate(records, start=1):
        if not isinstance(record, dict):
            raise ValueError(f"{path}: VMR record {row_number} is not a JSON object")
    return records


def build_vmr_accession_index(vmr_records):
    index = {}
    with_accession = 0
    for record in vmr_records:
        accessions = {
            _base_accession(accession)
            for accession in _accessions(record.get(ACCESSION_COLUMN))
        }
        if accessions:
            with_accession += 1
        for accession in accessions:
            index.setdefault(accession, []).append(record)
    stats = {
        "total_records": len(vmr_records),
        "records_with_accession_field": with_accession,
        "records_without_accession_field": len(vmr_records) - with_accession,
        "distinct_base_accessions": len(index),
    }
    return index, stats


def report_unrecognized_keys(vmr_records):
    recognized = set(VMR_COLUMN_FIELD_IDS) | {ACCESSION_COLUMN}
    columns = sorted(
        {key for record in vmr_records for key in record} - recognized
    )
    return columns, len(columns)


def _field_value(record, column, warnings, line_label):
    value = record.get(column)
    if _is_empty(value):
        return None
    if isinstance(value, str):
        value = value.strip()
    if column == "Exemplar_or_additional_isolate":
        if isinstance(value, bool):
            return value
        if isinstance(value, str):
            normalized = value.strip().casefold()
            if normalized in ("e", "true", "yes"):
                return True
            if normalized in ("a", "false", "no"):
                return False
        warnings.append(
            f"{line_label}: VMR column {column!r} has unsupported value "
            f"{value!r}; it was not mapped to 'is_ictv_exemplar'."
        )
        return None
    return value


def _row_values(record, schema, warnings, line_label):
    values = {}
    for column, field_id in VMR_COLUMN_FIELD_IDS.items():
        value = _field_value(record, column, warnings, line_label)
        if value is None:
            continue
        if field_id not in schema:
            warnings.append(
                f"{line_label}: VMR column {column!r} maps to schema field ID "
                f"{field_id!r}, but that ID is absent from the configured schema; "
                "the value was not stored."
            )
            continue
        values[field_id] = value
    return values


def merge_vmr_into_normalized(normalized, vmr_records, old_doc, schema,
                              line_label, warnings):
    """Merge non-empty VMR fields; the main importer handles coercion/updates."""
    records = list(vmr_records)
    if len(records) > 1 and old_doc is not None:
        ictv_id_entry = schema.get("ictv_id")
        existing_ictv_id = (
            old_doc.get(ictv_id_entry["index"]) if ictv_id_entry else None
        )
        if not _is_empty(existing_ictv_id):
            matching_records = [
                record for record in records
                if str(record.get("ICTV_ID", "")).strip()
                == str(existing_ictv_id).strip()
            ]
            if matching_records:
                records = matching_records
                warnings.append(
                    f"{line_label}: accession occurs in multiple VMR rows; "
                    f"selected the row matching existing ICTV_ID "
                    f"{existing_ictv_id!r}."
                )

    field_values = {}
    for record in records:
        row_values = _row_values(record, schema, warnings, line_label)
        for field_id, value in row_values.items():
            field_values.setdefault(field_id, [])
            if value not in field_values[field_id]:
                field_values[field_id].append(value)

    for field_id, values in field_values.items():
        if len(values) > 1:
            warnings.append(
                f"{line_label}: matched VMR rows have conflicting non-empty "
                f"values for schema field ID {field_id!r}; using the first "
                "value in the VMR JSON."
            )
        value = values[0]
        existing_value = normalized.get(field_id)
        if not _is_empty(existing_value) and existing_value != value:
            warnings.append(
                f"{line_label}: VMR value for schema field ID {field_id!r} "
                "differs from NCBI/FASTA metadata and takes precedence."
            )
        normalized[field_id] = value