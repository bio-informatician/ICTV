#!/usr/bin/env python3
"""
k.update_existing_records.py

Updates records that ALREADY EXIST in the vjdbv1.1 metadata index
(default `metadata_index_2026_08_05` on https://virjendb.org:9201/), and --
new in this version -- pushes fresh sequences to `genome_sequences` when an
accession's version changes. It never creates a record for a brand-new vjid
that isn't already in the index (those are skipped and reported).

-------------------------------------------------------------------------
EVENT MODEL (per accession, per run)
-------------------------------------------------------------------------
Every accession resolves to exactly ONE of:

1. INSERTION  (no log, no versioning)
   - Filling fields that are currently empty on the active record, and/or
   - The recomputed "always-set" fields (0bg ictv_vmr_version, 0ai/0ah
     refseq/genbank flags) are (re)written in place.
   - Applied as a partial update to the existing active document. 0do/0cx/
     0hd/0fh are untouched.

2. CONTENT UPDATE  (0aa / accession+version unchanged, but >=1 field's VALUE
   differs from what's already stored, non-empty -> non-empty-and-different)
   - The active document is flagged 0hd=false.
   - A new document is indexed: a full duplicate of the old one, with the
     changed field(s) applied, empty fields filled too, 0do incremented,
     0hd=true, 0fh renewed to today. 0dn and 0cx are unchanged.
   - 0dk (sequence_source) append ("GenBank" added if missing) counts as a
     content change if the list was non-empty before.
   - One document is written to the change_log index (see below).

3. VERSION-DIFF  (0aa's accession is the same but its NCBI version number
   changed, e.g. AB000403.1 -> AB000403.2)
   - The active document is flagged 0hd=false.
   - A brand new document is indexed: NOTHING is copied from the old one
     except 0dn. 0do and 0cx are both incremented. 0hd=true, 0fh=today.
     Every field the current NCBI report supplies is written fresh
     (unlogged -- see change_log note below).
   - One SUMMARY document is written to the change_log index: a single
     {field: "insdc_accession", old_value, new_value} entry -- the rest of
     the fresh record is NOT diffed/logged field-by-field.
   - The new sequence is pushed to genome_sequences (0dn + 0cx as the
     composite key -- see GENOME SEQUENCES below).

Only CONTENT UPDATE and VERSION-DIFF ever touch 0do, 0cx, 0hd, 0fh, or the
change_log index. INSERTION never does.

-------------------------------------------------------------------------
CHANGE LOG
-------------------------------------------------------------------------
Written to a SEPARATE index (default `change_log`, see --change-log-index).
One document per CONTENT UPDATE or VERSION-DIFF event:

    {
      "0dn": "<vjid>",
      "0do": <new metadata version>,
      "time": "YYYY-MM-DD",
      "method": "<from config.json>",
      "reason": "<from config.json>",
      "changes": [
        {"field": "<schema field id, e.g. host_species>",
         "old_value": "<JSON-encoded>",
         "new_value": "<JSON-encoded>"},
        ...
      ]
    }

"method"/"reason" come from --config (JSON), with independent values for
"content_update" and "version_change" event types. Default config used if
--config is omitted or the file doesn't exist:

    {
      "version_change": {"method": "NCBI_Version_Bump",
                          "reason": "Accession version incremented at source"},
      "content_update": {"method": "Metadata_Refresh",
                          "reason": "NCBI metadata content updated"}
    }

Bookkeeping fields (0hd, 0do, 0cx, 0fh) and the recomputed 0bg/0ai/0ah flags
are NEVER logged, even though they change on every content-update/version-
diff event -- only genuine sourced-content changes are (0dk append included).

-------------------------------------------------------------------------
GENOME SEQUENCES  (only on VERSION-DIFF)
-------------------------------------------------------------------------
Same field aliases as j.insert_sequences_to_vjdb.py, PLUS the new 0cx
(sequence_version) composite-key field:

    0aa -> insdc_accession_id  (full accession WITH version)
    0dn -> virjen_id
    0cx -> sequence_version    (NEW -- same value just computed for the
                                 corresponding metadata record)
    0hh -> genome_sequence
    0hi -> genome_hash

This requires the sequence body itself, so the FASTA parser in this script
now captures the sequence lines under each header (not just the header),
unlike the previous version which only read headers.

ASSUMPTIONS FLAGGED FOR CONFIRMATION (see chat reply):
  - change_log index name/field names (0dn/0do reused, but "time"/"method"/
    "reason"/"changes"/"field"/"old_value"/"new_value" are plain, unmapped
    field names -- change_log is not part of DB_Scheme.json).
  - 0fh (renewed date) field id, reused from the original script.
  - 0bg/0ai/0ah are always (re)written but never trigger an event by
    themselves and are never logged, even if their value changes.
  - If 0do or 0cx is missing on an old doc, it is treated as 1 before
    incrementing (per your answer).

-------------------------------------------------------------------------
INPUT
-------------------------------------------------------------------------
Multi-FASTA file, headers like:

    >vj000013520911 | AB000403 | organism=Tenuivirus oryzabrevis | length=2704 | description=...
    ACGTACGT...
    ACGTACGT...

Usage:
    # Dry run only (default, safe, no network writes)
    python3 k.update_existing_records.py --fasta VMR_MSL41.v1.20260320.fasta

    # Dry run + real writes if dry run is clean
    python3 k.update_existing_records.py --fasta VMR_MSL41.v1.20260320.fasta --commit
"""

import argparse
import concurrent.futures
import hashlib
import json
import logging
import os
import re
import subprocess
import sys
import threading
import time
from copy import deepcopy
from datetime import datetime
from pathlib import Path

import requests
from tqdm import tqdm

# --------------------------------------------------------------------------
# Fixed / constant configuration
# --------------------------------------------------------------------------

ES_URL_DEFAULT = "https://virjendb.org:9201"
INDEX_NAME_DEFAULT = "metadata_index_2026_08_05"
GENOME_INDEX_DEFAULT = "genome_sequences"
CHANGE_LOG_INDEX_DEFAULT = "change_log"
AUTH_TOKEN_DEFAULT = "OERyN3JKOEJPaThtV3FrMkUtWVE6Q3FqSVFidU02NF9tT29VWXloQ1RjQQ=="

SCHEMA_URL_DEFAULT = (
    "https://raw.githubusercontent.com/VirJenDB/virjendb_metadata/"
    "refs/heads/v1.1/json/DB_Scheme.json"
)

# Hardcoded field ids -- NOT resolved via schema (kept consistent with the
# rest of the pipeline's convention of treating these as fixed control
# fields rather than ordinary content fields).
VJID_FIELD = "0dn"                 # virjen id (composite key part 1)
VERSION_FIELD = "0do"              # metadata version (composite key part 2)
SEQ_VERSION_FIELD = "0cx"          # sequence_version -- NEW field
ACTIVE_FLAG_FIELD = "0hd"          # true = current active record
RENEWED_DATE_FIELD = "0fh"         # date the record was last (re)written

VMR_VERSION_FIELD = "0bg"          # ictv_vmr_version -- always overwritten
GENBANK_SOURCE_VALUE = "GenBank"

IS_REFSEQ_FIELD_ID = "is_refseq_sequence"      # -> 0ai (schema-resolved)
IS_GENBANK_FIELD_ID = "is_genbank_sequence"    # -> 0ah (schema-resolved)
REFSEQ_SOURCE_KEY = "is_annotated_refseq"      # key produced by normalize_report()

INSDC_ACCESSION_KEY = "insdc_accession"        # normalize_report() target key
INSDC_ACCESSION_FIELD_ID = "insdc_accession"   # schema field_id (-> 0aa)

FIELD_ID_OVERRIDES = {
    "virus_organism_name": "organism_name",
    "virus_tax_id_ncbi": "tax_id_ncbi",
    "sequence_hash_refseq": "genome_sequence_hash",
    "sra_accessions": "sra_accession",
    "geographic_region": "collection_region",
    "update_date_refseq": "bvbrc_update_date",
    "source_database": "sequence_source",   # handled specially (0dk append)
}

# Fields intentionally dropped with no warning, or handled specially and
# never run through the generic fill/update comparison loop.
SILENT_DROP_FIELDS = {
    "protein_count",
    "gene_count",
    "is_lab_host",
    "source_database",       # -> 0dk, handled specially (append rule)
    REFSEQ_SOURCE_KEY,       # -> 0ai/0ah, handled specially (always-set)
}

SEMICOLON_LIST_FIELDS = {
    "lineage_ids_ncbi": int,
    "lineage_names_ncbi": str,
}

REPORT_SCHEMA_MAP = {
    "accession": "insdc_accession",
    "bioprojects[]": "bioproject_accession",
    "completeness": "genome_completeness",
    "length": "sequence_length",
    "release_date": "release_date",
    "update_date": "update_date_refseq",
    "gene_count": "gene_count",
    "is_annotated": "is_annotated_refseq",
    "protein_count": "protein_count",
    "source_database": "source_database",
    "is_lab_host": "is_lab_host",
    "lab_host": "lab_host_species",
    "segment": "segment_name",
    "submitter.names[]": "submitter",
    "submitter.affiliation": "submitter_organization",
    "submitter.country": "submitter_country",
    "isolate.name": "isolate_ncbi",
    "nucleotide.sequence_hash": "sequence_hash_refseq",
    "virus.organism_name": "virus_organism_name",
    "virus.tax_id": "virus_tax_id_ncbi",
    "host.common_name": "host_common_name",
    "host.organism_name": "host_species",
    "host.tax_id": "host_ncbi_tax_id",
    "host.species": "host_species_ncbi",
    "host.genus": "host_genus_ncbi",
    "host.family": "host_family_ncbi",
    "isolate.collection_date": "collection_date",
    "isolate.geo_location": "collection_geo_location",
    "isolate.source": "collection_source_host_tissue",
    "location.geographic_location": "collection_geo_location",
    "location.geographic_region": "geographic_region",
    "location.usa_state": "ncbi_usa",
    "biosample": "biosample_accession",
    "sra_accessions[]": "sra_accessions",
}

TAXONOMIC_SUFFIXES = (
    ("_viria", "realm_ncbi"),
    ("_virae", "kingdom_ncbi"),
    ("_viricota", "phylum_ncbi"),
    ("_viricetes", "class_ncbi"),
    ("_virales", "order_ncbi"),
    ("_viridae", "family_ncbi"),
    ("_virus", "genus_ncbi"),
)


# --------------------------------------------------------------------------
# Config (method/reason per event type)
# --------------------------------------------------------------------------

DEFAULT_CONFIG = {
    "version_change": {
        "method": "NCBI_Version_Bump",
        "reason": "Accession version incremented at source",
    },
    "content_update": {
        "method": "Metadata_Refresh",
        "reason": "NCBI metadata content updated",
    },
}


def load_config(config_path):
    if not config_path or not os.path.exists(config_path):
        return deepcopy(DEFAULT_CONFIG)
    with open(config_path, "r", encoding="utf-8") as f:
        cfg = json.load(f)
    for key in ("version_change", "content_update"):
        if key not in cfg or "method" not in cfg[key] or "reason" not in cfg[key]:
            raise ValueError(
                f"config file missing required '{key}.method' / '{key}.reason'"
            )
    return cfg


# --------------------------------------------------------------------------
# schema.py logic, inlined
# --------------------------------------------------------------------------

def get_nested(obj, path, default=None):
    current = obj
    for part in path.split("."):
        if part.endswith("[]"):
            part = part[:-2]
            if not isinstance(current, dict):
                return default
            current = current.get(part)
            if current is None:
                return default
            if not isinstance(current, list):
                return default
            return deepcopy(current)
        else:
            if not isinstance(current, dict):
                return default
            if part not in current:
                return default
            current = current[part]
    return deepcopy(current)


def set_if_not_none(result, key, value):
    if value is None or value == "" or value == []:
        return
    result[key] = value


def iso_date(value):
    if value is None:
        return None
    return value.split("T")[0]


def lineage_names(lineage):
    if lineage is None:
        return []
    return [n["name"] for n in lineage if n.get("name") is not None]


def lineage_ids(lineage):
    if lineage is None:
        return []
    return [str(n["tax_id"]) for n in lineage if n.get("tax_id") is not None]


def virus_taxonomy(lineage):
    taxonomy = {}
    if lineage is None:
        return taxonomy
    valid = [n for n in lineage if "name" in n]
    if not valid:
        return taxonomy
    for node in valid:
        name = node["name"]
        if name == "Viruses":
            continue
        for suffix, output in TAXONOMIC_SUFFIXES:
            if name.endswith(suffix):
                taxonomy[output] = name
                break
    used = set(taxonomy.values())
    for node in reversed(valid):
        name = node["name"]
        if name == "Viruses":
            continue
        if name not in used:
            taxonomy["species_ncbi"] = name
            break
    return taxonomy


def normalize_report(report):
    result = {}
    for source, target in REPORT_SCHEMA_MAP.items():
        value = get_nested(report, source)
        if source in ("release_date", "update_date"):
            value = iso_date(value)
        set_if_not_none(result, target, value)

    isolate = report.get("isolate", {}) or {}
    location = report.get("location", {}) or {}
    geo = isolate.get("geo_location") or location.get("geographic_location")
    if geo is not None:
        result["collection_geo_location"] = geo
    region = location.get("geographic_region")
    if region is not None:
        result["geographic_region"] = region

    host = report.get("host")
    if host is not None:
        h_lineage = host.get("lineage")
        names = lineage_names(h_lineage)
        ids = lineage_ids(h_lineage)
        if names:
            result["lineage_names_ncbi"] = ";".join(names)
        if ids:
            result["lineage_ids_ncbi"] = ";".join(ids)

    virus = report.get("virus")
    if virus is not None:
        result.update(virus_taxonomy(virus.get("lineage")))

    for key in [k for k, v in result.items() if v is None or v == "" or v == []]:
        del result[key]

    return result


# --------------------------------------------------------------------------
# type coercion
# --------------------------------------------------------------------------

_RANGE_DATE_RE = re.compile(
    r"^\[\s*(\d{4}(?:-\d{2}(?:-\d{2})?)?)\s+TO\s+(\d{4}(?:-\d{2}(?:-\d{2})?)?)\s*\]$",
    re.IGNORECASE,
)
_MISORDERED_DATE_RE = re.compile(r"^(\d{2})-(\d{2})-(\d{4})$")


def _normalize_date(v, field_label):
    v = v.strip()
    m = _RANGE_DATE_RE.match(v)
    if m:
        normalized_start, _ = _normalize_date(m.group(1), field_label)
        return normalized_start, (
            f"{field_label}: uncertain date range {v!r} collapsed to range "
            f"start {normalized_start!r} (earliest possible date used)"
        )
    for fmt in ("%Y-%m-%d", "%Y-%m", "%Y"):
        try:
            datetime.strptime(v, fmt)
            return v, None
        except ValueError:
            continue
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
    if vjdb_type == "integer":
        if isinstance(value, bool):
            raise ValueError(f"{field_label}: expected integer, got boolean")
        if isinstance(value, int):
            return value, None
        if isinstance(value, str) and value.strip().lstrip("-").isdigit():
            return int(value.strip()), None
        raise ValueError(f"{field_label}: expected integer, got {value!r}")
    if vjdb_type == "float":
        if isinstance(value, bool):
            raise ValueError(f"{field_label}: expected float, got boolean")
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
    if isinstance(value, (int, float, bool)):
        return value, None
    return value, None


def coerce_value(value, vjdb_type, field_label):
    if isinstance(value, list):
        out, warnings = [], []
        for v in value:
            c, w = _coerce_scalar(v, vjdb_type, field_label)
            out.append(c)
            if w:
                warnings.append(w)
        return out, warnings
    c, w = _coerce_scalar(value, vjdb_type, field_label)
    return c, ([w] if w else [])


def load_schema(schema_path=None, schema_url=SCHEMA_URL_DEFAULT):
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
            "type": entry.get("vjdbv1.1_fields_type", "keyword"),
            "name": entry.get("vjdbv1.1_name", field_id),
        }
    print(f"Schema loaded: {len(schema)} field(s).")
    if SEQ_VERSION_FIELD not in {e["index"] for e in schema.values()} and "sequence_version" not in schema:
        print(
            f"WARNING: 'sequence_version' ({SEQ_VERSION_FIELD}) not found in schema -- "
            f"per your confirmation it should already exist in DB_Scheme.json. "
            f"Double-check the schema source if this warning appears.",
            file=sys.stderr,
        )
    return schema


def resolve_field(jsonl_field, schema):
    schema_field_id = FIELD_ID_OVERRIDES.get(jsonl_field, jsonl_field)
    entry = schema.get(schema_field_id)
    if entry is None:
        return None, None, None
    return schema_field_id, entry["index"], entry["type"]


def is_empty(value):
    return value is None or value == "" or value == []


def strip_version(accession: str) -> str:
    accession = accession.strip()
    return accession.rsplit(".", 1)[0] if "." in accession else accession


def sha256_hex(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def to_bool(value):
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() == "true"
    return bool(value)


# --------------------------------------------------------------------------
# FASTA parsing (now captures the sequence body, needed for genome_sequences)
# --------------------------------------------------------------------------

def vmr_version_from_filename(input_path):
    name = os.path.basename(input_path)
    for ext in (".fasta", ".txt", ".json", ".jsonl"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def parse_fasta(fasta_path):
    """Yield {'vjid', 'accession', 'organism', 'length', 'description',
    'sequence'} for each '>' header. Headers look like:
        >vj000013520911 | AB000403 | organism=... | length=... | description=...
    followed by 0+ raw sequence lines, concatenated (no newlines, upper-cased).
    """
    logger = logging.getLogger("update_vjdb")
    records = []
    current_rec = None
    current_seq = []

    def flush():
        if current_rec is not None:
            current_rec["sequence"] = "".join(current_seq).replace(" ", "").upper()
            records.append(current_rec)

    with open(fasta_path, "r") as fh:
        for line_no, line in enumerate(fh, start=1):
            line = line.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue
            if line.startswith(">"):
                flush()
                current_rec = None
                current_seq = []

                header = line[1:].strip()
                parts = [p.strip() for p in header.split("|")]
                if len(parts) < 2:
                    logger.warning(
                        "line %d: header does not have at least vjid|accession -- skipped: %r",
                        line_no, header,
                    )
                    continue  # current_rec stays None -- no stale carry-over
                vjid, accession = parts[0], parts[1]
                extra = {}
                for kv in parts[2:]:
                    if "=" in kv:
                        k, v = kv.split("=", 1)
                        extra[k.strip()] = v.strip()
                current_rec = {
                    "vjid": vjid,
                    "accession": accession,
                    "organism": extra.get("organism"),
                    "length": extra.get("length"),
                    "description": extra.get("description"),
                }
            else:
                if current_rec is not None:
                    current_seq.append(line.strip())
        flush()

    return records


# --------------------------------------------------------------------------
# NCBI Datasets resolution
# --------------------------------------------------------------------------

class FetchNotFound(Exception):
    pass


class RateLimiter:
    def __init__(self, rate_per_sec):
        self.min_interval = 1.0 / rate_per_sec
        self.lock = threading.Lock()
        self._next_time = 0.0

    def wait(self):
        with self.lock:
            now = time.monotonic()
            target = max(now, self._next_time)
            self._next_time = target + self.min_interval
            sleep_for = target - now
        if sleep_for > 0:
            time.sleep(sleep_for)


_eutils_limiter = None
_datasets_limiter = None
_NCBI_API_KEY = None


def latest_accession(accession: str) -> str:
    _eutils_limiter.wait()
    cmd = r'''
set -o pipefail
esearch -db nucleotide -query "{accession}" 2>/dev/null |
efetch -format docsum 2>/dev/null |
xtract -pattern DocumentSummary -element AccessionVersion |
sort -V |
tail -n1
'''.format(accession=accession)

    env = os.environ.copy()
    if _NCBI_API_KEY:
        env["NCBI_API_KEY"] = _NCBI_API_KEY

    p = subprocess.run(
        cmd, shell=True, executable="/bin/bash", env=env,
        capture_output=True, text=True,
    )
    if p.returncode != 0:
        raise RuntimeError(
            f"esearch/efetch failed for {accession}\nSTDOUT:\n{p.stdout}\nSTDERR:\n{p.stderr}"
        )
    version = p.stdout.strip()
    if not version:
        raise FetchNotFound(f"{accession}: not resolvable via esearch/efetch\nSTDERR:\n{p.stderr}")
    return version


def fetch_datasets_report(accession: str):
    version = latest_accession(accession)
    cmd = ["datasets", "summary", "virus", "genome", "accession", version]
    if _NCBI_API_KEY:
        cmd += ["--api-key", _NCBI_API_KEY]
    _datasets_limiter.wait()
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(p.stderr.strip())
    obj = json.loads(p.stdout)
    if obj.get("total_count", 0) == 0:
        return version, None
    return version, obj["reports"][0]


def fetch_datasets_report_with_retries(accession, retries, base_backoff):
    last_err = None
    for attempt in range(retries):
        try:
            return fetch_datasets_report(accession)
        except FetchNotFound as e:
            last_err = e
            if attempt < retries - 1:
                time.sleep(base_backoff * (2 ** attempt))
    raise last_err


_cache_write_lock = threading.Lock()
_accession_locks = {}
_accession_locks_guard = threading.Lock()


def _lock_for_accession(accession):
    with _accession_locks_guard:
        lock = _accession_locks.get(accession)
        if lock is None:
            lock = threading.Lock()
            _accession_locks[accession] = lock
        return lock


def _cache_path(cache_dir: Path, accession: str) -> Path:
    safe = accession.replace("/", "_").replace(" ", "_")
    return cache_dir / f"{safe}.json"


def load_cached(cache_dir: Path, accession: str):
    path = _cache_path(cache_dir, accession)
    if not path.exists():
        return None
    with path.open() as f:
        return json.load(f)


def save_cached(cache_dir: Path, accession: str, version: str, report):
    cache_dir.mkdir(parents=True, exist_ok=True)
    with _cache_path(cache_dir, accession).open("w") as f:
        json.dump({"version": version, "report": report}, f)


def resolve_one(vjid, accession, cache_dir, use_cache, retries, retry_backoff):
    with _lock_for_accession(accession):
        try:
            cached = load_cached(cache_dir, accession) if use_cache else None
            if cached is not None:
                return vjid, accession, cached["version"], cached["report"], True, None
            version, report = fetch_datasets_report_with_retries(accession, retries, retry_backoff)
            if use_cache:
                with _cache_write_lock:
                    save_cached(cache_dir, accession, version, report)
            return vjid, accession, version, report, False, None
        except FetchNotFound as e:
            return vjid, accession, None, None, False, str(e)
        except Exception as e:
            return vjid, accession, None, None, False, str(e)


# --------------------------------------------------------------------------
# Elasticsearch helpers
# --------------------------------------------------------------------------

def build_session(token: str) -> requests.Session:
    session = requests.Session()
    session.headers.update({"Authorization": f"ApiKey {token}"})
    return session


def check_auth(session, es_url, index_name, verify_ssl):
    logger = logging.getLogger("update_vjdb")
    url = f"{es_url.rstrip('/')}/{index_name}/_count"
    resp = session.get(url, verify=verify_ssl, timeout=30)
    if resp.status_code == 401:
        raise SystemExit(f"Auth failed (401) against {url} -- check the API key.")
    if resp.status_code == 403:
        raise SystemExit(f"Auth OK but forbidden (403) against {url} -- key lacks permission on this index.")
    resp.raise_for_status()
    count = resp.json().get("count", 0)
    logger.info("Auth OK -- %s/%s has %d document(s).", es_url, index_name, count)


def find_active_by_vjid(session, es_url, index_name, vjids, verify_ssl, batch_size=200):
    """Look up the ACTIVE (0hd == true) document for each vjid. Returns
    {vjid: (es_id, source) or (None, None)}.

    Filtering on 0hd is required now that content-update/version-diff events
    can leave multiple documents sharing the same 0dn (old inactive + new
    active) -- a plain terms-on-0dn lookup would risk returning a stale
    duplicate."""
    out = {v: (None, None) for v in vjids}
    url = f"{es_url.rstrip('/')}/{index_name}/_search"
    logger = logging.getLogger("update_vjdb")
    n_batches = (len(vjids) + batch_size - 1) // batch_size
    progress = tqdm(
        range(0, len(vjids), batch_size), total=n_batches,
        unit="batch", desc="Looking up active docs", file=sys.__stderr__,
    )
    for start in progress:
        batch = vjids[start:start + batch_size]
        payload = {
            "query": {"bool": {"must": [
                {"terms": {VJID_FIELD: batch}},
                {"term": {ACTIVE_FLAG_FIELD: True}},
            ]}},
            "size": len(batch),
            "_source": True,
        }
        try:
            resp = session.post(url, json=payload, verify=verify_ssl, timeout=30)
            resp.raise_for_status()
        except requests.exceptions.RequestException as e:
            logger.error("Lookup batch %d-%d failed: %s", start, start + len(batch), e)
            raise
        hits = resp.json().get("hits", {}).get("hits", [])
        for hit in hits:
            source = hit.get("_source", {})
            vjid = source.get(VJID_FIELD)
            if vjid in out:
                out[vjid] = (hit["_id"], source)
        progress.set_postfix(found=sum(1 for eid, _ in out.values() if eid is not None))
    return out


def execute_bulk(session, es_url, index_name, actions, verify_ssl, timeout=180, chunk_size=2000):
    """actions: flat list of NDJSON-ready dicts (action-line, source-line, ...).
    chunk_size is set to an even number to keep action/document pairs safely together."""
    if not actions:
        return {"errors": False, "items": []}
        
    url = f"{es_url.rstrip('/')}/{index_name}/_bulk"
    
    all_items = []
    has_errors = False
    
    # Process in chunks to avoid HTTP 413 Request Entity Too Large
    for i in range(0, len(actions), chunk_size):
        chunk = actions[i:i + chunk_size]
        payload = "\n".join(json.dumps(a) for a in chunk) + "\n"
        
        resp = session.post(
            url, data=payload.encode("utf-8"),
            headers={"Content-Type": "application/x-ndjson"},
            verify=verify_ssl, timeout=timeout,
        )
        resp.raise_for_status()
        resp_json = resp.json()
        
        # Accumulate results across all batches so check_bulk_errors still works
        if resp_json.get("errors"):
            has_errors = True
        if "items" in resp_json:
            all_items.extend(resp_json["items"])
            
    return {"errors": has_errors, "items": all_items}


def check_bulk_errors(resp_json, logger, label):
    if not resp_json.get("errors"):
        return True
    failed = []
    for item in resp_json.get("items", []):
        for op, detail in item.items():
            if detail.get("error"):
                failed.append(detail["error"])
    logger.error("%s: %d/%d item(s) failed", label, len(failed), len(resp_json.get("items", [])))
    for err in failed[:20]:
        logger.error("  %s", err)
    return False


def set_index_write_block(session, es_url, index_name, blocked, verify_ssl):
    url = f"{es_url.rstrip('/')}/{index_name}/_settings"
    payload = {"index.blocks.write": True if blocked else None}
    resp = session.put(url, json=payload, verify=verify_ssl, timeout=60)
    resp.raise_for_status()
    return resp.json()


# --------------------------------------------------------------------------
# Per-record event evaluation
# --------------------------------------------------------------------------

def build_record_actions(vjid, es_id, old_doc, normalized, resolved_version,
                          sequence, schema, config, vmr_version, line_label):
    """Returns a dict describing exactly one event:
        {event, errors, warnings,
         old_doc_update: {"_id","doc"} | None,
         new_meta_doc: {...} | None,
         fill_update: {"_id","doc"} | None,
         change_log_doc: {...} | None,
         genome_doc: {...} | None}
    """
    errors, warnings = [], []
    today_str = datetime.utcnow().strftime("%Y-%m-%d")

    old_version = old_doc.get(VERSION_FIELD) or 1
    old_seq_version = old_doc.get(SEQ_VERSION_FIELD) or 1

    insdc_entry = schema.get(INSDC_ACCESSION_FIELD_ID)
    if insdc_entry is None:
        errors.append(f"{line_label}: schema field '{INSDC_ACCESSION_FIELD_ID}' (0aa) not found in schema")
        return {"event": "error", "errors": errors, "warnings": warnings,
                "old_doc_update": None, "new_meta_doc": None,
                "fill_update": None, "change_log_doc": None, "genome_doc": None}
    acc_index = insdc_entry["index"]
    old_accessions = old_doc.get(acc_index) or []
    if isinstance(old_accessions, str):
        old_accessions = [old_accessions]

    # --- always-set / recomputed fields (never logged, never trigger an event alone)
    is_annotated_refseq = normalized.get(REFSEQ_SOURCE_KEY, False)
    is_refseq = to_bool(is_annotated_refseq)
    always_set = {}
    refseq_entry = schema.get(IS_REFSEQ_FIELD_ID)
    genbank_entry = schema.get(IS_GENBANK_FIELD_ID)
    if refseq_entry:
        always_set[refseq_entry["index"]] = is_refseq
    else:
        errors.append(f"{line_label}: schema field '{IS_REFSEQ_FIELD_ID}' not found")
    if genbank_entry:
        always_set[genbank_entry["index"]] = True
    else:
        errors.append(f"{line_label}: schema field '{IS_GENBANK_FIELD_ID}' not found")
    if vmr_version:
        always_set[VMR_VERSION_FIELD] = vmr_version

    result = {"event": None, "errors": errors, "warnings": warnings,
              "old_doc_update": None, "new_meta_doc": None,
              "fill_update": None, "change_log_doc": None, "genome_doc": None}

    # ---------------------------------------------------------------
    # 1. VERSION-DIFF: same base accession, different version suffix
    # ---------------------------------------------------------------
    if old_accessions and old_accessions[0]:
        old_acc = old_accessions[0]
        if strip_version(old_acc) == strip_version(resolved_version) and old_acc != resolved_version:
            new_version = old_version + 1
            new_seq_version = old_seq_version + 1

            new_doc = {
                VJID_FIELD: vjid,
                VERSION_FIELD: new_version,
                SEQ_VERSION_FIELD: new_seq_version,
                ACTIVE_FLAG_FIELD: True,
                RENEWED_DATE_FIELD: today_str,
                acc_index: [resolved_version],
            }
            for jsonl_field, value in normalized.items():
                if jsonl_field in SILENT_DROP_FIELDS:
                    continue
                if jsonl_field in SEMICOLON_LIST_FIELDS and isinstance(value, str):
                    caster = SEMICOLON_LIST_FIELDS[jsonl_field]
                    parts = [p.strip() for p in value.split(";") if p.strip()]
                    try:
                        value = [caster(p) for p in parts]
                    except ValueError:
                        errors.append(f"{line_label}: field '{jsonl_field}' could not cast ';'-list {value!r}")
                        continue
                schema_field_id, idx, vtype = resolve_field(jsonl_field, schema)
                if idx is None:
                    errors.append(
                        f"{line_label}: field '{jsonl_field}' not found in schema, no override, "
                        f"not in SILENT_DROP_FIELDS -- unmapped."
                    )
                    continue
                try:
                    coerced, w = coerce_value(value, vtype, f"{line_label} field '{jsonl_field}' ({idx})")
                except ValueError as e:
                    errors.append(str(e))
                    continue
                warnings.extend(w)
                new_doc[idx] = coerced
            # sequence_source: fresh record -> just "GenBank" plus whatever normalize gave (rare)
            dk_entry = schema.get("sequence_source")
            if dk_entry and dk_entry["index"] not in new_doc:
                new_doc[dk_entry["index"]] = [GENBANK_SOURCE_VALUE]
            new_doc.update(always_set)

            result["event"] = "version-diff"
            result["old_doc_update"] = {"_id": es_id, "doc": {ACTIVE_FLAG_FIELD: False}}
            result["new_meta_doc"] = new_doc
            result["change_log_doc"] = {
                VJID_FIELD: vjid,
                VERSION_FIELD: new_version,
                "time": today_str,
                "method": config["version_change"]["method"],
                "reason": config["version_change"]["reason"],
                "changes": [{
                    "field": INSDC_ACCESSION_FIELD_ID,
                    "old_value": json.dumps(old_acc),
                    "new_value": json.dumps(resolved_version),
                }],
            }
            if sequence:
                result["genome_doc"] = {
                    "0aa": resolved_version,
                    "0dn": vjid,
                    SEQ_VERSION_FIELD: new_seq_version,
                    "0hh": sequence,
                    "0hi": sha256_hex(sequence),
                }
            else:
                warnings.append(
                    f"{line_label}: version-diff event but no sequence body found in the "
                    f"FASTA record -- genome_sequences was NOT updated for this accession."
                )
            return result

    # ---------------------------------------------------------------
    # 2. CONTENT UPDATE / FILL (0aa unchanged)
    # ---------------------------------------------------------------
    fill_dict = {}
    updates_dict = {}
    changes_log = []

    for jsonl_field, value in normalized.items():
        if jsonl_field in SILENT_DROP_FIELDS:
            continue
        if jsonl_field in SEMICOLON_LIST_FIELDS and isinstance(value, str):
            caster = SEMICOLON_LIST_FIELDS[jsonl_field]
            parts = [p.strip() for p in value.split(";") if p.strip()]
            try:
                value = [caster(p) for p in parts]
            except ValueError:
                errors.append(f"{line_label}: field '{jsonl_field}' could not cast ';'-list {value!r}")
                continue
        schema_field_id, idx, vtype = resolve_field(jsonl_field, schema)
        if idx is None:
            errors.append(
                f"{line_label}: field '{jsonl_field}' not found in schema, no override, "
                f"not in SILENT_DROP_FIELDS -- unmapped."
            )
            continue
        try:
            coerced, w = coerce_value(value, vtype, f"{line_label} field '{jsonl_field}' ({idx})")
        except ValueError as e:
            errors.append(str(e))
            continue
        warnings.extend(w)

        old_val = old_doc.get(idx)
        if is_empty(old_val):
            fill_dict[idx] = coerced
        elif old_val != coerced:
            updates_dict[idx] = coerced
            changes_log.append({
                "field": schema_field_id,
                "old_value": json.dumps(old_val),
                "new_value": json.dumps(coerced),
            })

    # sequence_source (0dk): always append "GenBank" if missing from the list.
    # Empty-before -> fill (unlogged). Non-empty-before and value actually
    # changes -> logged content update.
    dk_entry = schema.get("sequence_source")
    if dk_entry is None:
        errors.append(f"{line_label}: schema field 'sequence_source' (0dk) not found")
    else:
        dk_idx = dk_entry["index"]
        existing_dk = old_doc.get(dk_idx)
        if existing_dk is None:
            dk_list = []
        elif isinstance(existing_dk, list):
            dk_list = list(existing_dk)
        else:
            dk_list = [existing_dk]
        if GENBANK_SOURCE_VALUE not in dk_list:
            new_dk_list = dk_list + [GENBANK_SOURCE_VALUE]
            if not dk_list:
                fill_dict[dk_idx] = new_dk_list
            else:
                updates_dict[dk_idx] = new_dk_list
                changes_log.append({
                    "field": "sequence_source",
                    "old_value": json.dumps(dk_list),
                    "new_value": json.dumps(new_dk_list),
                })

    if updates_dict:
        new_version = old_version + 1
        new_doc = deepcopy(old_doc)
        new_doc[VERSION_FIELD] = new_version
        new_doc[ACTIVE_FLAG_FIELD] = True
        new_doc[RENEWED_DATE_FIELD] = today_str
        new_doc.update(fill_dict)
        new_doc.update(updates_dict)
        new_doc.update(always_set)

        result["event"] = "content-update"
        result["old_doc_update"] = {"_id": es_id, "doc": {ACTIVE_FLAG_FIELD: False}}
        result["new_meta_doc"] = new_doc
        result["change_log_doc"] = {
            VJID_FIELD: vjid,
            VERSION_FIELD: new_version,
            "time": today_str,
            "method": config["content_update"]["method"],
            "reason": config["content_update"]["reason"],
            "changes": changes_log,
        }
        return result

    # ---------------------------------------------------------------
    # 3. INSERTION (fill-only and/or always-set, no genuine content change)
    # ---------------------------------------------------------------
    fill_dict.update(always_set)
    if fill_dict:
        result["event"] = "insertion"
        result["fill_update"] = {"_id": es_id, "doc": fill_dict}
    else:
        result["event"] = "no-change"
    return result


# --------------------------------------------------------------------------
# Main pipeline
# --------------------------------------------------------------------------

class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            s.write(data)
            s.flush()

    def flush(self):
        for s in self.streams:
            s.flush()


def run(args):
    global _NCBI_API_KEY, _eutils_limiter, _datasets_limiter
    _NCBI_API_KEY = args.ncbi_api_key or os.environ.get("NCBI_API_KEY") or os.environ.get("API_KEY")

    logger = logging.getLogger("update_vjdb")
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s", "%Y-%m-%d %H:%M:%S"))
    logger.addHandler(handler)

    default_rate = 8.0 if _NCBI_API_KEY else 2.5
    eutils_rate = args.eutils_rate if args.eutils_rate is not None else default_rate
    datasets_rate = args.datasets_rate if args.datasets_rate is not None else default_rate
    _eutils_limiter = RateLimiter(eutils_rate)
    _datasets_limiter = RateLimiter(datasets_rate)
    logger.info("Rate limits: eutils=%.1f req/sec, datasets=%.1f req/sec", eutils_rate, datasets_rate)

    config = load_config(args.config)
    logger.info("Config: version_change=%s, content_update=%s",
                config["version_change"]["method"], config["content_update"]["method"])

    fasta_records = parse_fasta(args.fasta)
    logger.info("Parsed %d record(s) from %s", len(fasta_records), args.fasta)
    n_with_seq = sum(1 for r in fasta_records if r.get("sequence"))
    logger.info("%d/%d record(s) include a sequence body (needed for version-diff -> genome_sequences)",
                n_with_seq, len(fasta_records))

    if args.vmr_version is None:
        args.vmr_version = vmr_version_from_filename(args.fasta)
    logger.info("VMR version: %s", args.vmr_version)

    schema = load_schema(schema_path=args.schema, schema_url=args.schema_url)
    cache_dir = Path(args.cache_dir)
    use_cache = not args.no_cache

    verify_ssl = not args.no_verify_ssl
    session = build_session(args.token)
    check_auth(session, args.es_url, args.index, verify_ssl)

    # --- Step 1: resolve accessions ----------------------------------------
    resolved = {}
    fasta_by_vjid = {r["vjid"]: r for r in fasta_records}
    resolve_errors = []
    no_report_count = 0
    cache_hits = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [
            executor.submit(resolve_one, rec["vjid"], rec["accession"], cache_dir, use_cache,
                             args.retries, args.retry_backoff)
            for rec in fasta_records
        ]
        progress = tqdm(concurrent.futures.as_completed(futures), total=len(futures),
                         unit="acc", desc="Resolving accessions", file=sys.__stderr__)
        for fut in progress:
            vjid, accession, version, report, was_cached, err = fut.result()
            if err:
                resolve_errors.append(f"{vjid} ({accession}): {err}")
            else:
                resolved[vjid] = (accession, version, report)
                if report is None:
                    no_report_count += 1
                if was_cached:
                    cache_hits += 1
            progress.set_postfix(cached=cache_hits, failed=len(resolve_errors), no_report=no_report_count)

    if resolve_errors:
        logger.warning("%d accession(s) could not be resolved:", len(resolve_errors))
        for e in resolve_errors[:20]:
            logger.warning("  %s", e)

    # --- Step 2: look up ACTIVE existing docs (0dn + 0hd:true) -------------
    all_vjids = list(resolved.keys())
    lookup = find_active_by_vjid(session, args.es_url, args.index, all_vjids, verify_ssl, batch_size=args.batch_size)
    missing_vjids = [v for v, (eid, _) in lookup.items() if eid is None]
    if missing_vjids:
        logger.warning("%d vjid(s) have no ACTIVE document -- SKIPPED (never created):", len(missing_vjids))
        for v in missing_vjids[:20]:
            logger.warning("  %s", v)

    # --- Step 3: build per-record actions -----------------------------------
    all_errors, all_warnings = [], []
    meta_bulk = []          # flat NDJSON action list for the metadata index
    change_log_bulk = []    # flat NDJSON action list for the change_log index
    genome_bulk = []        # flat NDJSON action list for genome_sequences
    event_counts = {"insertion": 0, "content-update": 0, "version-diff": 0, "no-change": 0, "error": 0}

    for vjid, (accession, version, report) in resolved.items():
        es_id, old_doc = lookup.get(vjid, (None, None))
        if es_id is None:
            continue
        normalized = normalize_report(report) if report is not None else {}
        sequence = fasta_by_vjid.get(vjid, {}).get("sequence") or ""
        line_label = f"vjid {vjid} ({accession})"

        r = build_record_actions(vjid, es_id, old_doc, normalized, version, sequence,
                                  schema, config, args.vmr_version, line_label)
        all_errors.extend(r["errors"])
        all_warnings.extend(r["warnings"])
        event_counts[r["event"]] = event_counts.get(r["event"], 0) + 1

        if r["old_doc_update"]:
            meta_bulk.append({"update": {"_id": r["old_doc_update"]["_id"]}})
            meta_bulk.append({"doc": r["old_doc_update"]["doc"]})
        if r["new_meta_doc"] is not None:
            meta_bulk.append({"index": {}})
            meta_bulk.append(r["new_meta_doc"])
        if r["fill_update"]:
            meta_bulk.append({"update": {"_id": r["fill_update"]["_id"]}})
            meta_bulk.append({"doc": r["fill_update"]["doc"]})
        if r["change_log_doc"] is not None:
            change_log_bulk.append({"index": {}})
            change_log_bulk.append(r["change_log_doc"])
        if r["genome_doc"] is not None:
            genome_bulk.append({"index": {}})
            genome_bulk.append(r["genome_doc"])

    # --- Report --------------------------------------------------------------
    print(f"\nParsed:              {len(fasta_records)}")
    print(f"Resolved via NCBI:   {len(resolved)}  (failed: {len(resolve_errors)}, no-report: {no_report_count})")
    print(f"Active doc found:    {len(resolved) - len(missing_vjids)}  (missing/skipped: {len(missing_vjids)})")
    print(f"Events: insertion={event_counts.get('insertion', 0)}  "
          f"content-update={event_counts.get('content-update', 0)}  "
          f"version-diff={event_counts.get('version-diff', 0)}  "
          f"no-change={event_counts.get('no-change', 0)}")
    print(f"Errors:              {len(all_errors)}")
    print(f"Warnings:            {len(all_warnings)}")

    if all_warnings:
        print("\n--- Sample warnings (first 10) ---")
        for w in all_warnings[:10]:
            print("  WARN:", w)

    if all_errors:
        print("\n--- Errors (first 20) ---")
        for e in all_errors[:20]:
            print("  ERROR:", e)
        print(f"\nDry run FAILED with {len(all_errors)} error(s). Aborting -- nothing was written.")
        sys.exit(1)

    print("\nDry run PASSED (no hard errors).")
    if not args.commit:
        print(f"Run again with --commit to write: "
              f"{len(meta_bulk)//2} metadata action(s), "
              f"{len(change_log_bulk)//2} change_log doc(s), "
              f"{len(genome_bulk)//2} genome_sequences doc(s).")
        return

    if meta_bulk:
        print(f"\n=== Writing {len(meta_bulk)//2} metadata action(s) to {args.index} ===")
        resp = execute_bulk(session, args.es_url, args.index, meta_bulk, verify_ssl)
        check_bulk_errors(resp, logger, "metadata bulk")

    if change_log_bulk:
        print(f"=== Writing {len(change_log_bulk)//2} change_log doc(s) to {args.change_log_index} ===")
        resp = execute_bulk(session, args.es_url, args.change_log_index, change_log_bulk, verify_ssl)
        check_bulk_errors(resp, logger, "change_log bulk")

    if genome_bulk:
        print(f"=== Writing {len(genome_bulk)//2} doc(s) to {args.genome_index} ===")
        try:
            set_index_write_block(session, args.es_url, args.genome_index, False, verify_ssl)
            resp = execute_bulk(session, args.es_url, args.genome_index, genome_bulk, verify_ssl)
            check_bulk_errors(resp, logger, "genome_sequences bulk")
        finally:
            set_index_write_block(session, args.es_url, args.genome_index, True, verify_ssl)

    print("\nDone.")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True, help="Multi-FASTA file with '>vjid | accession | ...' headers + sequence bodies")
    parser.add_argument("--config", default="config.json", help="JSON file with version_change/content_update method+reason")
    parser.add_argument("--schema", default=None, help="Local DB_Scheme.json (default: fetch fresh from GitHub)")
    parser.add_argument("--schema-url", default=SCHEMA_URL_DEFAULT)
    parser.add_argument("--es-url", default=ES_URL_DEFAULT)
    parser.add_argument("--index", default=INDEX_NAME_DEFAULT)
    parser.add_argument("--genome-index", default=GENOME_INDEX_DEFAULT)
    parser.add_argument("--change-log-index", default=CHANGE_LOG_INDEX_DEFAULT)
    parser.add_argument("--token", default=AUTH_TOKEN_DEFAULT)
    parser.add_argument("--no-verify-ssl", action="store_true")
    parser.add_argument("--commit", action="store_true")
    parser.add_argument("--batch-size", type=int, default=200)
    parser.add_argument("--cache-dir", default=".cache/reports")
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--ncbi-api-key", default=None)
    parser.add_argument("--eutils-rate", type=float, default=None)
    parser.add_argument("--datasets-rate", type=float, default=None)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--retry-backoff", type=float, default=1.0)
    parser.add_argument("--vmr-version", default=None)
    parser.add_argument("--log-file", default=None)
    args = parser.parse_args()

    input_base = Path(args.fasta).stem
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = args.log_file or f"{input_base}_update_{timestamp}.log"
    log_fh = open(log_path, "w", encoding="utf-8")
    sys.stdout = Tee(sys.__stdout__, log_fh)
    print(f"Logging full output to: {log_path}")
    try:
        run(args)
    finally:
        sys.stdout = sys.__stdout__
        log_fh.close()


if __name__ == "__main__":
    main()
