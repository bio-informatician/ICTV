#!/usr/bin/env python3
"""
insert_jsonl_to_vjdb.py

Maps records from a VirJenDB-style jsonl file to the vjdbv1.1 DB field-index
scheme (0aa, 0ab, ...) and inserts them into the Elasticsearch index
`metadata_index_2026_08_05` on https://virjendb.org:9201/.

Workflow:
    1. Fetches the current vjdbv1.1 DB_Scheme.json fresh from GitHub as soon
       as the script starts (pass --schema to use a local file instead).
    2. Always runs a DRY RUN first: builds every mapped document in memory,
       validates every field against the schema's declared type, and
       reports any problems WITHOUT touching the database.
    3. Only if the dry run is 100% clean (no hard errors) does it proceed
       to the real bulk insert -- and only if --commit was passed.

Usage:
    # Dry run only (default, safe, no network writes; schema pulled from GitHub)
    python3 insert_jsonl_to_vjdb.py --input 41.jsonl

    # Dry run + real insertion if dry run is clean
    python3 insert_jsonl_to_vjdb.py --input 41.jsonl --commit

    # Use a local schema file instead of fetching from GitHub
    python3 insert_jsonl_to_vjdb.py --input 41.jsonl --schema DB_Scheme.json
"""

import argparse
import csv
import json
import os
import sys
from datetime import date, datetime

import requests
from requests.auth import AuthBase

# --------------------------------------------------------------------------
# Fixed / constant configuration
# --------------------------------------------------------------------------

ES_URL_DEFAULT = "https://virjendb.org:9201"
INDEX_NAME = "metadata_index_2026_08_05"
AUTH_TOKEN_DEFAULT = (
    "OERyN3JKOEJPaThtV3FrMkUtWVE6Q3FqSVFidU02NF9tT29VWXloQ1RjQQ=="
)

# Schema is fetched fresh from GitHub every run by default (see load_schema()).
SCHEMA_URL_DEFAULT = (
    "https://raw.githubusercontent.com/VirJenDB/virjendb_metadata/"
    "refs/heads/v1.1/json/DB_Scheme.json"
)

VJDB_ACCESSION_FIELD = "0dn"          # virjendb_accession, auto-incremented
# Fallback only -- normally the real latest accession is looked up live from
# the index at run time (see fetch_latest_accession()) so this never goes
# stale. Only used if --offline-latest-accession is passed explicitly.
LATEST_KNOWN_ACCESSION = "vj000016771170"
ACCESSION_PREFIX = "vj"
ACCESSION_DIGITS = 12                 # vj + 12 digits = matches example above

INSDC_ACCESSION_FIELD = "0aa"         # insdc_accession -- the original
                                       # GenBank/RefSeq accession, kept
                                       # alongside the generated vjid so
                                       # genome sequences can be attached
                                       # to the right record later.

INSERTION_DATE_FIELD = "0fg"          # insertion_date
VMR_VERSION_FIELD = "0bg"             # ictv_vmr_version

METADATA_VERSION_FIELD = "0do"        # virjendb_metadata_version
METADATA_VERSION_VALUE = 1            # always 1

REFSEQ_SOURCE_FIELD = "is_annotated_refseq"   # jsonl key (special-cased, not
                                               # a plain 1:1 override -- see
                                               # map_record())
IS_REFSEQ_FIELD_ID = "is_refseq_sequence"     # -> 0ai
IS_GENBANK_FIELD_ID = "is_genbank_sequence"   # -> 0ah
# A RefSeq record is always also a GenBank record: if is_annotated_refseq is
# true, both is_refseq_sequence and is_genbank_sequence are set true; if not
# true, only is_genbank_sequence is set true (is_refseq_sequence -> false).

# --------------------------------------------------------------------------
# jsonl_field_name -> schema field_id overrides.
# Only needed where the jsonl key does NOT literally match a
# `vjdbv1.1_field_id` in the schema file. Everything else is resolved
# automatically by exact field_id match.
# --------------------------------------------------------------------------
FIELD_ID_OVERRIDES = {
    "virus_organism_name": "organism_name",
    "virus_tax_id_ncbi": "tax_id_ncbi",
    "sequence_hash_refseq": "genome_sequence_hash",
    "sra_accessions": "sra_accession",
    "geographic_region": "collection_region",
    "update_date_refseq": "bvbrc_update_date",   # 0cw, date
    "source_database": "sequence_source",         # 0dk, keyword
}

# jsonl fields with NO corresponding DB column at all (confirmed by
# scanning the full 233-entry schema -- there is genuinely no field for
# these). They are dropped with a warning rather than guessed at.
UNMAPPED_FIELDS = set()

# jsonl fields that are intentionally dropped with NO warning at all (by
# default). Add field names here to silence their "no corresponding DB
# column" warning entirely.
SILENT_DROP_FIELDS = {
    "protein_count",
    "is_lab_host",
    "gene_count",
}

# jsonl fields whose value is a ';'-delimited string that must be split
# into a list before type validation/insertion (schema types them as
# scalar integer/keyword but the intent is clearly a list).
SEMICOLON_LIST_FIELDS = {
    "lineage_ids_ncbi": int,     # split -> list[int]
    "lineage_names_ncbi": str,   # split -> list[str]
}


# --------------------------------------------------------------------------
# Schema loading
# --------------------------------------------------------------------------

def load_schema(schema_path=None, schema_url=SCHEMA_URL_DEFAULT):
    """Return dict: field_id -> {"index": "0xx", "type": "keyword"/"text"/...}

    By default the schema is fetched fresh from GitHub (schema_url) every
    run, so the script always maps against the current v1.1 scheme. Pass
    schema_path to force loading a local file instead (e.g. offline testing).
    """
    if schema_path:
        print(f"Loading schema from local file: {schema_path}")
        with open(schema_path, "r", encoding="utf-8") as f:
            raw = json.load(f)
    else:
        print(f"Fetching schema from GitHub: {schema_url}")
        resp = requests.get(schema_url, timeout=30)
        resp.raise_for_status()
        raw = resp.json()

    schema = {}
    for entry in raw:
        field_id = entry["vjdbv1.1_field_id"]
        schema[field_id] = {
            "index": entry["vjdbv1.1_field_index"],
            "type": entry["vjdbv1.1_fields_type"],
            "name": entry.get("vjdbv1.1_name", field_id),
        }
    print(f"Schema loaded: {len(schema)} field(s).")
    return schema


def resolve_field(jsonl_field, schema):
    """Return (schema_field_id, index, type) or (None, None, None) if unmapped."""
    if jsonl_field in UNMAPPED_FIELDS:
        return None, None, None

    schema_field_id = FIELD_ID_OVERRIDES.get(jsonl_field, jsonl_field)
    entry = schema.get(schema_field_id)
    if entry is None:
        return None, None, None
    return schema_field_id, entry["index"], entry["type"]


# --------------------------------------------------------------------------
# Type coercion / validation
# --------------------------------------------------------------------------

DATE_FMT = "%Y-%m-%d"


import re

_RANGE_DATE_RE = re.compile(
    r"^\[\s*(\d{4}(?:-\d{2}(?:-\d{2})?)?)\s+TO\s+(\d{4}(?:-\d{2}(?:-\d{2})?)?)\s*\]$",
    re.IGNORECASE,
)


_MISORDERED_DATE_RE = re.compile(r"^(\d{2})[-/](\d{2})[-/](\d{4})$")


def _normalize_date(v, field_label):
    """Return (normalized_date_str, warning_or_None). Raises ValueError if
    the value cannot be interpreted as a date at all.

    Partial dates (YYYY-MM, YYYY) are passed through UNCHANGED -- they are
    not padded into full YYYY-MM-DD dates. Only genuinely mis-formatted
    full dates (e.g. DD-MM-YYYY, DD/MM/YYYY) are corrected, since those are
    a formatting error rather than a legitimate partial date."""
    v = v.strip()

    # GenBank uncertain-date range syntax, e.g. "[2018 TO 2019]"
    m = _RANGE_DATE_RE.match(v)
    if m:
        start_raw = m.group(1)
        # Recurse on the range start (itself may be YYYY / YYYY-MM / YYYY-MM-DD)
        normalized_start, _ = _normalize_date(start_raw, field_label)
        return normalized_start, (
            f"{field_label}: uncertain date range {v!r} collapsed to range "
            f"start {normalized_start!r} (earliest possible date used)"
        )

    # Full YYYY-MM-DD -- already correct, use as-is
    try:
        datetime.strptime(v, "%Y-%m-%d")
        return v, None
    except ValueError:
        pass

    # Partial: YYYY-MM (common for GenBank collection_date / release dates)
    # -- kept as-is, NOT padded to a full date.
    try:
        datetime.strptime(v, "%Y-%m")
        return v, None
    except ValueError:
        pass

    # Partial: YYYY only -- kept as-is, NOT padded to a full date.
    try:
        datetime.strptime(v, "%Y")
        return v, None
    except ValueError:
        pass

    # Mis-formatted full date, e.g. DD-MM-YYYY or DD/MM/YYYY -- this is a
    # genuine formatting error (wrong field order), so it gets corrected
    # to YYYY-MM-DD.
    m2 = _MISORDERED_DATE_RE.match(v)
    if m2:
        first, second, year = m2.groups()
        for day, month, order_label in (
            (first, second, "DD-MM-YYYY"),
            (second, first, "MM-DD-YYYY"),
        ):
            candidate = f"{year}-{month}-{day}"
            try:
                datetime.strptime(candidate, "%Y-%m-%d")
                return candidate, (
                    f"{field_label}: mis-formatted date {v!r} interpreted as "
                    f"{order_label} and corrected to {candidate!r}"
                )
            except ValueError:
                continue

    raise ValueError(f"{field_label}: date {v!r} is not a recognizable YYYY[-MM[-DD]] value")


def _coerce_scalar(value, vjdb_type, field_label):
    """Coerce/validate a single scalar value against its declared schema type.
    Returns (coerced_value, warning_or_None).
    Raises ValueError with a descriptive message on failure."""
    if vjdb_type == "integer":
        if isinstance(value, bool):  # bool is a subclass of int -- reject
            raise ValueError(f"{field_label}: expected integer, got boolean")
        if isinstance(value, int):
            return value, None
        if isinstance(value, str) and value.strip().lstrip("-").isdigit():
            return int(value.strip()), None
        raise ValueError(f"{field_label}: expected integer, got {value!r}")

    if vjdb_type == "float":
        if isinstance(value, bool):
            raise ValueError(f"{field_label}: expected float, got boolean")
        if isinstance(value, (int, float)):
            return float(value), None
        try:
            return float(value), None
        except (TypeError, ValueError):
            raise ValueError(f"{field_label}: expected float, got {value!r}")

    if vjdb_type == "boolean":
        if isinstance(value, bool):
            return value, None
        if isinstance(value, str) and value.strip().lower() in ("true", "false"):
            return value.strip().lower() == "true", None
        raise ValueError(f"{field_label}: expected boolean, got {value!r}")

    if vjdb_type == "date":
        if isinstance(value, str):
            return _normalize_date(value, field_label)
        raise ValueError(f"{field_label}: expected date string, got {value!r}")

    # text / keyword / geo_point / anything else -> pass through as string
    if isinstance(value, (int, float, bool)):
        return value, None  # numbers stored under a text/keyword field: keep as-is
    return value, None


def coerce_value(value, vjdb_type, field_label):
    """Handle lists transparently: coerce each element, preserve list-ness.
    Returns (coerced_value, list_of_warnings)."""
    if isinstance(value, list):
        out = []
        warnings = []
        for v in value:
            coerced, warn = _coerce_scalar(v, vjdb_type, field_label)
            out.append(coerced)
            if warn:
                warnings.append(warn)
        return out, warnings

    coerced, warn = _coerce_scalar(value, vjdb_type, field_label)
    return coerced, ([warn] if warn else [])


# --------------------------------------------------------------------------
# Accession ID generator
# --------------------------------------------------------------------------

def fetch_latest_accession(es_url, index_name, token, verify_ssl=True):
    """Query the index for the current highest virjendb_accession (0dn) so
    the counter always starts from the real DB state instead of a constant
    that goes stale as soon as anything else gets inserted. Equivalent to:

        GET <index_name>/_search
        {
          "size": 1,
          "sort": [{"0dn": {"order": "desc"}}],
          "_source": ["0dn"]
        }
    """
    session = requests.Session()
    session.auth = ApiKeyAuth(token)
    session.headers.update({"Content-Type": "application/json"})

    url = f"{es_url.rstrip('/')}/{index_name}/_search"
    query = {
        "size": 1,
        "sort": [{VJDB_ACCESSION_FIELD: {"order": "desc"}}],
        "_source": [VJDB_ACCESSION_FIELD],
    }

    resp = session.post(url, json=query, verify=verify_ssl, timeout=30)
    resp.raise_for_status()
    hits = resp.json().get("hits", {}).get("hits", [])
    if not hits:
        raise RuntimeError(
            f"No documents found in index '{index_name}' when looking up "
            "the latest accession -- cannot determine where to start the "
            "counter. Pass --offline-latest-accession to override."
        )

    latest = hits[0]["_source"][VJDB_ACCESSION_FIELD]
    print(f"Latest accession in '{index_name}': {latest}")
    return latest


def make_accession_counter(latest_known=LATEST_KNOWN_ACCESSION):
    numeric_part = int(latest_known[len(ACCESSION_PREFIX):])
    counter = {"n": numeric_part}

    def next_id():
        counter["n"] += 1
        return f"{ACCESSION_PREFIX}{counter['n']:0{ACCESSION_DIGITS}d}"

    return next_id

def vmr_version_from_filename(input_path):
    """'ncbi/VMR_MSL41.v1.20260320.txt' -> 'VMR_MSL41.v1.20260320'."""
    name = os.path.basename(input_path)
    for ext in (".txt", ".json", ".jsonl"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name

# --------------------------------------------------------------------------
# Record mapping
# --------------------------------------------------------------------------

def map_record(
    raw_record,
    schema,
    line_no,
    insertion_date,
    next_accession_id,
    vmr_version,
):
    """
    Returns (es_document: dict, errors: list[str], warnings: list[str])
    es_document is the fully mapped 0xx-keyed document ready for indexing.
    """
    errors = []
    warnings = []
    doc = {}

    for jsonl_field, value in raw_record.items():
        if value is None:
            continue

        # Fields intentionally dropped with no warning at all (e.g. protein_count).
        if jsonl_field in SILENT_DROP_FIELDS:
            continue

        # RefSeq/GenBank relationship: a RefSeq record is always also a
        # GenBank record. is_annotated_refseq has no 1:1 schema column --
        # it fans out into two boolean fields.
        if jsonl_field == REFSEQ_SOURCE_FIELD:
            refseq_label = f"line {line_no} field '{jsonl_field}'"
            try:
                is_refseq, refseq_warn = _coerce_scalar(value, "boolean", refseq_label)
            except ValueError as e:
                errors.append(str(e))
                continue
            if refseq_warn:
                warnings.append(refseq_warn)

            for field_id, bool_value in (
                (IS_REFSEQ_FIELD_ID, is_refseq),
                (IS_GENBANK_FIELD_ID, True),
            ):
                entry = schema.get(field_id)
                if entry is None:
                    errors.append(
                        f"line {line_no}: schema field '{field_id}' "
                        f"(derived from '{jsonl_field}') not found in schema"
                    )
                    continue
                index = entry["index"]
                if index in doc:
                    errors.append(
                        f"line {line_no}: DB field '{index}' already set "
                        f"(collision from derived field '{field_id}')"
                    )
                    continue
                doc[index] = bool_value
            continue

        # Split ';'-delimited fields into real lists first
        if jsonl_field in SEMICOLON_LIST_FIELDS and isinstance(value, str):
            caster = SEMICOLON_LIST_FIELDS[jsonl_field]
            parts = [p.strip() for p in value.split(";") if p.strip() != ""]
            try:
                value = [caster(p) for p in parts]
            except ValueError:
                errors.append(
                    f"line {line_no}: field '{jsonl_field}' could not split/cast "
                    f"';'-list value {value!r} to {caster.__name__}"
                )
                continue

        schema_field_id, index, vjdb_type = resolve_field(jsonl_field, schema)

        if index is None:
            if jsonl_field in UNMAPPED_FIELDS:
                warnings.append(
                    f"line {line_no}: field '{jsonl_field}' has no corresponding "
                    f"DB column -- value {value!r} was dropped"
                )
            else:
                errors.append(
                    f"line {line_no}: field '{jsonl_field}' not found in schema "
                    f"and has no override -- unmapped"
                )
            continue

        try:
            coerced, coerce_warnings = coerce_value(
                value, vjdb_type, f"line {line_no} field '{jsonl_field}' ({index})"
            )
            warnings.extend(coerce_warnings)
        except ValueError as e:
            errors.append(str(e))
            continue

        if index in doc:
            errors.append(
                f"line {line_no}: DB field '{index}' already set (collision "
                f"from jsonl field '{jsonl_field}')"
            )
            continue

        doc[index] = coerced

    # Inject the fixed / auto fields
    doc[VJDB_ACCESSION_FIELD] = next_accession_id()
    doc[INSERTION_DATE_FIELD] = insertion_date
    doc[VMR_VERSION_FIELD] = vmr_version
    doc[METADATA_VERSION_FIELD] = METADATA_VERSION_VALUE

    return doc, errors, warnings


# --------------------------------------------------------------------------
# Dry run / main pipeline
# --------------------------------------------------------------------------

def run_pipeline(input_path, schema_path=None, schema_url=SCHEMA_URL_DEFAULT,
                  latest_known_accession=LATEST_KNOWN_ACCESSION):
    schema = load_schema(schema_path=schema_path, schema_url=schema_url)
    next_accession_id = make_accession_counter(latest_known_accession)
    insertion_date = date.today().strftime(DATE_FMT)

    documents = []
    all_errors = []
    all_warnings = []
    vmr_version = vmr_version_from_filename(input_path)
    print(f"VMR version (from filename): {vmr_version}")

    with open(input_path, "r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            try:
                record = json.loads(line)
            except json.JSONDecodeError as e:
                all_errors.append(f"line {line_no}: invalid JSON ({e})")
                continue

            doc, errors, warnings = map_record(
		   record,
		   schema,
		   line_no,
		   insertion_date,
 		   next_accession_id,
 		   vmr_version,
		)
            all_errors.extend(errors)
            all_warnings.extend(warnings)
            documents.append((line_no, doc))

    return documents, all_errors, all_warnings


def print_report(documents, errors, warnings, sample_n=2):
    print(f"Processed {len(documents)} record(s).")
    print(f"Errors:   {len(errors)}")
    print(f"Warnings: {len(warnings)} (dropped/unmapped fields -- non-fatal)")

    if warnings:
        print("\n--- Sample warnings (first 10) ---")
        for w in warnings[:10]:
            print("  WARN:", w)
        if len(warnings) > 10:
            print(f"  ... and {len(warnings) - 10} more")

    if errors:
        print("\n--- Errors (first 20) ---")
        for e in errors[:20]:
            print("  ERROR:", e)
        if len(errors) > 20:
            print(f"  ... and {len(errors) - 20} more")

    if documents:
        print(f"\n--- Sample mapped document(s) ---")
        for line_no, doc in documents[:sample_n]:
            print(f"  line {line_no}:")
            print(json.dumps(doc, indent=2, ensure_ascii=False))


def write_accession_csv(documents, csv_path):
    """Write a csv with columns: vjid, accession.

    vjid      = the generated virjendb_accession (vj + 12 digits)
    accession = the original INSDC/GenBank accession from the source record

    This gives a simple lookup table to attach genome sequences to the
    right record later, keyed by either column.
    """
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["vjid", "accession"])
        for _, doc in documents:
            vjid = doc.get(VJDB_ACCESSION_FIELD, "")
            accession = doc.get(INSDC_ACCESSION_FIELD, "")
            writer.writerow([vjid, accession])
    print(f"Wrote {len(documents)} vjid/accession pair(s) to: {csv_path}")


class Tee:
    """Duplicates writes to both the terminal and a log file, so the full
    dry-run/insert report is always saved somewhere, not just printed."""

    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            s.write(data)
            s.flush()

    def flush(self):
        for s in self.streams:
            s.flush()


# --------------------------------------------------------------------------
# Elasticsearch bulk insert
# --------------------------------------------------------------------------

class ApiKeyAuth(AuthBase):
    def __init__(self, token):
        self.token = token

    def __call__(self, r):
        r.headers["Authorization"] = f"ApiKey {self.token}"
        return r


def bulk_insert(documents, es_url, index_name, token, verify_ssl=True, batch_size=500):
    """Insert documents via the Elasticsearch _bulk API."""
    session = requests.Session()
    session.auth = ApiKeyAuth(token)
    session.headers.update({"Content-Type": "application/x-ndjson"})

    url = f"{es_url.rstrip('/')}/_bulk"

    total = len(documents)
    inserted = 0
    failures = []

    for start in range(0, total, batch_size):
        batch = documents[start:start + batch_size]
        lines = []
        for line_no, doc in batch:
            accession = doc.get(VJDB_ACCESSION_FIELD)
            action = {"index": {"_index": index_name, "_id": accession}}
            lines.append(json.dumps(action))
            lines.append(json.dumps(doc, ensure_ascii=False))
        payload = "\n".join(lines) + "\n"

        resp = session.post(url, data=payload.encode("utf-8"), verify=verify_ssl, timeout=120)
        if resp.status_code >= 300:
            failures.append((start, resp.status_code, resp.text[:2000]))
            continue

        result = resp.json()
        if result.get("errors"):
            for item in result["items"]:
                op = item.get("index", {})
                if op.get("status", 200) >= 300:
                    failures.append((op.get("_id"), op.get("status"), op.get("error")))
                else:
                    inserted += 1
        else:
            inserted += len(batch)

        print(f"  Batch {start}-{start + len(batch)}: submitted "
              f"({resp.status_code}), running total inserted={inserted}")

    return inserted, failures


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Map & insert VirJenDB jsonl into ES.")
    parser.add_argument("--input", required=True, help="Path to the source .jsonl file")
    parser.add_argument("--schema", default=None,
                         help="Path to a local DB_Scheme.json to use INSTEAD of "
                              "fetching from GitHub. By default the schema is "
                              "always fetched fresh from --schema-url.")
    parser.add_argument("--schema-url", default=SCHEMA_URL_DEFAULT,
                         help="URL to fetch DB_Scheme.json from (ignored if "
                              "--schema is given).")
    parser.add_argument("--es-url", default=ES_URL_DEFAULT, help="Elasticsearch base URL")
    parser.add_argument("--index", default=INDEX_NAME, help="Target index name")
    parser.add_argument("--token", default=AUTH_TOKEN_DEFAULT,
                         help="Elasticsearch API key (sent as 'Authorization: ApiKey <token>')")
    parser.add_argument("--no-verify-ssl", action="store_true",
                         help="Disable TLS certificate verification (test servers)")
    parser.add_argument("--offline-latest-accession", default=None,
                         help="Skip the live ES lookup and start the vjid counter "
                              "from this accession instead (e.g. for offline "
                              "testing). By default the latest accession is "
                              "always fetched fresh from --index.")
    parser.add_argument("--commit", action="store_true",
                         help="Perform the real insertion if the dry run is clean. "
                              "Without this flag, only the dry run report is produced.")
    parser.add_argument("--log-file", default=None,
                         help="Path to write the full run log to (everything that's "
                              "printed to the console). Default: "
                              "<input>_<timestamp>.log next to the input file.")
    parser.add_argument("--csv-output", default=None,
                         help="Path to write the vjid/accession csv to (columns: "
                              "vjid, accession). Default: <input>_vjid_accession.csv "
                              "next to the input file.")
    args = parser.parse_args()

    input_base, _ = os.path.splitext(args.input)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    log_path = args.log_file or f"{input_base}_{timestamp}.log"
    csv_path = args.csv_output or f"{input_base}_vjid_accession.csv"

    log_fh = open(log_path, "w", encoding="utf-8")
    sys.stdout = Tee(sys.__stdout__, log_fh)

    try:
        _run(args, log_path, csv_path)
    finally:
        sys.stdout = sys.__stdout__
        log_fh.close()


def _run(args, log_path, csv_path):
    print(f"Logging full output to: {log_path}")

    if args.no_verify_ssl:
        import urllib3
        urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

    if args.offline_latest_accession:
        latest_known_accession = args.offline_latest_accession
        print(f"Using offline override for latest accession: {latest_known_accession}")
    else:
        latest_known_accession = fetch_latest_accession(
            args.es_url, args.index, args.token,
            verify_ssl=not args.no_verify_ssl,
        )

    print("=== DRY RUN ===")
    documents, errors, warnings = run_pipeline(
        args.input, schema_path=args.schema, schema_url=args.schema_url,
        latest_known_accession=latest_known_accession,
    )
    print_report(documents, errors, warnings)

    if errors:
        print(f"\nDry run FAILED with {len(errors)} error(s). Aborting -- nothing was inserted.")
        sys.exit(1)

    print("\nDry run PASSED (no hard errors).")

    write_accession_csv(documents, csv_path)

    if not args.commit:
        print("Run again with --commit to actually insert these "
              f"{len(documents)} record(s) into '{args.index}'.")
        print("(Note: the vjid/accession csv above reflects this dry run's "
              "assignment -- if you don't --commit now, re-running later "
              "may assign different vjids to the same rows.)")
        return

    print(f"\n=== INSERTING {len(documents)} record(s) into "
          f"{args.es_url}/{args.index} ===")

    inserted, failures = bulk_insert(
        documents, args.es_url, args.index, args.token,
        verify_ssl=not args.no_verify_ssl,
    )

    print(f"\nInserted: {inserted}/{len(documents)}")
    if failures:
        print(f"Failures: {len(failures)}")
        for f in failures[:20]:
            print("  ", f)
        print(f"\nWARNING: {len(failures)} record(s) failed to insert but are "
              f"still listed in {csv_path} -- cross-check against the "
              "failures above before using that csv to attach sequences.")


if __name__ == "__main__":
    main()
