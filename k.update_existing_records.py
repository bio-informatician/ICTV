#!/usr/bin/env python3
"""
update_existing_vjdb.py

Updates records that ALREADY EXIST in the vjdbv1.1 metadata index
(`metadata_index_2026_08_05` on https://virjendb.org:9201/). It does NOT create
new records and does NOT touch the `genome_sequences` index -- for genome
sequence ingestion see j.insert_sequences_to_vjdb.py.

-------------------------------------------------------------------------
INPUT
-------------------------------------------------------------------------
The ONLY input is a multi-FASTA file whose headers look like:

    >vj000013520911 | AB000403 | organism=Tenuivirus oryzabrevis | length=2704 | description=...

    field 1 (before first "|")  -> existing vjid           (0dn, used as ES _id)
    field 2 (after first "|")   -> INSDC accession, NO version

No json/csv/jsonl input is used or produced.

-------------------------------------------------------------------------
IMPORTANT: vjid is NOT the Elasticsearch _id
-------------------------------------------------------------------------
Documents in this index are keyed by an internal numeric _id (e.g. "9780307").
The vjid (e.g. "vj000013900745") is stored as a regular field, "0dn". This
script therefore looks records up via a `_search` with a `terms` query on
"0dn" (NOT `_mget` by _id), and uses the REAL `_id` returned by that search
for the actual `_bulk` update -- never the vjid itself.

-------------------------------------------------------------------------
WORKFLOW (per record)
-------------------------------------------------------------------------
1. Parse the FASTA header -> (vjid, accession_no_version).
2. Resolve the LATEST versioned accession for accession_no_version via
   esearch | efetch | xtract (same approach as NCBIDatasetsResolver in
   fetchers.py), then fetch the full report via
   `datasets summary virus genome accession <version>`.
   If NCBI Datasets has no virus-genome report for this accession (common
   for partial/isolate GenBank submissions), this is NOT treated as a
   resolution failure: the resolved version is kept and no NCBI metadata
   is used for that record (see step 6).
3. Fetch the current vjdbv1.1 DB_Scheme.json fresh from GitHub (pass
   --schema to use a local file instead) to resolve field_id -> (index, type).
4. Normalize the NCBI Datasets report into the same flat, jsonl-shaped
   fields used by insert_jsonl_to_vjdb.py (schema.py logic, inlined below).
   If no report was available, this is an empty dict -- no content to
   compare or fill in from this step.
5. Fetch the EXISTING document from Elasticsearch by searching "0dn" = vjid
   (batched terms query), which also gives us the real ES _id to update.
6. Build a PARTIAL update document:
     - 0dk (sequence_source): ALWAYS append "GenBank" to the existing list
       (existing value is normalized to a list first; "NCBI Virus" etc. is
       left untouched).
     - 0bg (ictv_vmr_version): ALWAYS set/overwritten to the current VMR
       version constant.
     - 0ai (is_refseq_sequence) / 0ah (is_genbank_sequence): ALWAYS
       recomputed from the report's `is_annotated_refseq` (a RefSeq record
       is always also a GenBank record -- see REFSEQ logic below) and
       overwritten. If there was no report, is_refseq defaults to False and
       is_genbank is still True (it's a resolved GenBank accession either way).
     - insdc_accession (0aa): filled with the resolved version (as a list,
       matching the field's stored shape) if currently empty.
     - Every OTHER mapped field (country/collection info, host info,
       taxonomy, etc.): only added to the update if it is MISSING or EMPTY
       in the existing document, AND normalized report data is available.
       If no report was available, no other field is touched at all --
       content comparison is skipped entirely for these records.
7. Records whose vjid is not found in the index are SKIPPED (never
   created) and reported at the end.
8. Always runs a DRY RUN first (prints the per-record update payloads,
   warnings, and errors). Only sends the real `_bulk` update requests to
   Elasticsearch if --commit is passed AND the dry run had 0 hard errors.

Usage:
    # Dry run only (default, safe, no network writes to ES)
    python3 update_existing_vjdb.py --fasta VMR_MSL41.v1.20260320.fasta

    # Dry run + real partial updates if dry run is clean
    python3 update_existing_vjdb.py --fasta VMR_MSL41.v1.20260320.fasta --commit
"""

import argparse
import concurrent.futures
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
# Fixed / constant configuration (kept consistent with i.insert_jsonl_to_vjdb.py
# and j.insert_sequences_to_vjdb.py)
# --------------------------------------------------------------------------

ES_URL_DEFAULT = "https://virjendb.org:9201"
INDEX_NAME_DEFAULT = "metadata_index_2026_08_05"
AUTH_TOKEN_DEFAULT = (
    "OERyN3JKOEJPaThtV3FrMkUtWVE6Q3FqSVFidU02NF9tT29VWXloQ1RjQQ=="
)

SCHEMA_URL_DEFAULT = (
    "https://raw.githubusercontent.com/VirJenDB/virjendb_metadata/"
    "refs/heads/v1.1/json/DB_Scheme.json"
)

VMR_VERSION_FIELD = "0bg"          # ictv_vmr_version

SOURCE_FIELD = "0dk"                # sequence_source -- a list; existing values
                                     # (e.g. "NCBI Virus") are kept, "GenBank" is
                                     # appended if not already present.
GENBANK_SOURCE_VALUE = "GenBank"

IS_REFSEQ_FIELD_ID = "is_refseq_sequence"     # -> 0ai
IS_GENBANK_FIELD_ID = "is_genbank_sequence"   # -> 0ah
REFSEQ_SOURCE_KEY = "is_annotated_refseq"      # key produced by normalize_report()
# A RefSeq record is always also a GenBank record: if is_annotated_refseq is
# true, both is_refseq_sequence and is_genbank_sequence are set true; if not
# true (including "no report at all"), only is_genbank_sequence is set true
# (is_refseq_sequence -> false).
# These two fields are ALWAYS (re)computed and overwritten, never fill-if-empty.

INSDC_ACCESSION_KEY = "insdc_accession"       # normalize_report() target key
INSDC_ACCESSION_FIELD_ID = "insdc_accession"  # schema field_id (-> 0aa)
# NOTE: 0aa is stored as a LIST in the DB (e.g. ["AJ867584.2"]), and its ES
# mapping uses a lowercase normalizer (ES itself lowercases for indexing), so
# the fill-if-empty value below is wrapped in a list to match the existing shape.

# vjid lives in this field -- it is NOT the Elasticsearch _id.
VJID_FIELD = "0dn"

# jsonl_field_name -> schema field_id overrides. Same table as
# i.insert_jsonl_to_vjdb.py (kept in sync manually).
FIELD_ID_OVERRIDES = {
    "virus_organism_name": "organism_name",
    "virus_tax_id_ncbi": "tax_id_ncbi",
    "sequence_hash_refseq": "genome_sequence_hash",
    "sra_accessions": "sra_accession",
    "geographic_region": "collection_region",
    "update_date_refseq": "bvbrc_update_date",
    "source_database": "sequence_source",   # handled specially, popped before
                                             # generic mapping -- see below.
}

# Fields intentionally dropped with no warning (confirmed no DB column, or
# handled specially and not meant to be fill-if-empty candidates).
SILENT_DROP_FIELDS = {
    "protein_count",
    "gene_count",            # same kind of report-only stat as protein_count;
                              # no DB column for it (see 2026-07-31 dry run)
    "is_lab_host",
    "source_database",       # handled specially (0dk append logic)
    REFSEQ_SOURCE_KEY,       # handled specially (0ai/0ah logic)
}

SEMICOLON_LIST_FIELDS = {
    "lineage_ids_ncbi": int,
    "lineage_names_ncbi": str,
}

# --------------------------------------------------------------------------
# schema.py logic, inlined (per instructions: no external file dependency --
# this script is fully self-contained apart from live HTTP/CLI calls).
# --------------------------------------------------------------------------

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
    """NCBI Datasets `virus genome accession` report -> flat jsonl-shaped dict
    (same target key names as main.py/schema.py's normalize())."""
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

    # remove empties
    for key in [k for k, v in result.items() if v is None or v == "" or v == []]:
        del result[key]

    return result


# --------------------------------------------------------------------------
# insert_jsonl_to_vjdb.py logic, inlined: field resolution + type coercion
# --------------------------------------------------------------------------

DATE_FMT = "%Y-%m-%d"

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
            "type": entry["vjdbv1.1_fields_type"],
            "name": entry.get("vjdbv1.1_name", field_id),
        }
    print(f"Schema loaded: {len(schema)} field(s).")
    return schema


def resolve_field(jsonl_field, schema):
    schema_field_id = FIELD_ID_OVERRIDES.get(jsonl_field, jsonl_field)
    entry = schema.get(schema_field_id)
    if entry is None:
        return None, None, None
    return schema_field_id, entry["index"], entry["type"]


# --------------------------------------------------------------------------
# FASTA parsing
# --------------------------------------------------------------------------

def strip_version(accession: str) -> str:
    accession = accession.strip()
    return accession.rsplit(".", 1)[0] if "." in accession else accession


def vmr_version_from_filename(input_path):
    """'VMR_MSL41.v1.20260320.fasta' -> 'VMR_MSL41.v1.20260320'.

    Same approach as vmr_version_from_filename() in i.insert_jsonl_to_vjdb.py,
    extended to also strip a .fasta extension (this script's input is a
    multi-FASTA file rather than a jsonl file)."""
    name = os.path.basename(input_path)
    for ext in (".fasta", ".txt", ".json", ".jsonl"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def parse_fasta_headers(fasta_path):
    """Yield {'vjid':..., 'accession':..., 'organism':..., 'length':...,
    'description':...} for each '>' header. Headers look like:
        >vj000013520911 | AB000403 | organism=... | length=... | description=...
    """
    records = []
    with open(fasta_path, "r") as fh:
        for line_no, line in enumerate(fh, start=1):
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            parts = [p.strip() for p in header.split("|")]
            if len(parts) < 2:
                logging.getLogger("update_vjdb").warning(
                    "line %d: header does not have at least vjid|accession -- skipped: %r",
                    line_no, header,
                )
                continue
            vjid, accession = parts[0], parts[1]
            extra = {}
            for kv in parts[2:]:
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    extra[k.strip()] = v.strip()
            records.append(
                {
                    "vjid": vjid,
                    "accession": accession,
                    "organism": extra.get("organism"),
                    "length": extra.get("length"),
                    "description": extra.get("description"),
                }
            )
    return records


# --------------------------------------------------------------------------
# NCBI Datasets resolution (same approach as fetchers.NCBIDatasetsResolver)
# --------------------------------------------------------------------------

class FetchNotFound(Exception):
    pass


class RateLimiter:
    """Paces calls to at most `rate_per_sec` per second, GLOBALLY across all
    threads (not per-thread) -- this is what actually keeps NCBI happy, not
    the --workers count. Threads block in .wait() until their turn."""

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


# Two separate budgets: eutils (esearch/efetch) and the NCBI Datasets API are
# different services with their own quotas. Each is set (in run(), once the
# API key is known) to a bit under NCBI's documented ~10 req/sec-with-key /
# ~3 req/sec-without-key ceiling, so this script self-throttles instead of
# relying on --workers to (accidentally) do it via failure.
_eutils_limiter = None
_datasets_limiter = None

# Set once at the start of run() from --ncbi-api-key or the NCBI_API_KEY /
# API_KEY environment variable, then passed explicitly on every esearch /
# efetch / datasets command line below -- rather than just relying on the
# subprocess inheriting the parent's environment and hoping each CLI tool
# picks it up on its own.
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
        cmd,
        shell=True,
        executable="/bin/bash",
        env=env,
        capture_output=True,
        text=True,
    )

    if p.returncode != 0:
        raise RuntimeError(
            f"esearch/efetch failed for {accession}\n"
            f"STDOUT:\n{p.stdout}\n"
            f"STDERR:\n{p.stderr}"
        )

    version = p.stdout.strip()
    if not version:
        raise FetchNotFound(
            f"{accession}: not resolvable via esearch/efetch\n"
            f"STDERR:\n{p.stderr}"
        )

    return version


def fetch_datasets_report(accession: str):
    """Returns (version, report_or_None).

    If NCBI Datasets has no virus-genome report for this accession (common
    for partial/isolate GenBank submissions -- confirmed to happen even when
    esearch/efetch resolves the accession just fine), this is NOT an error:
    we return (version, None) and the caller skips content comparison for
    this record, writing only the always-set fields.
    """
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
    """A FetchNotFound here now only comes from latest_accession() (esearch/
    efetch genuinely can't resolve the accession) -- a missing Datasets
    report is handled separately as (version, None), not as an exception.
    Retrying with backoff distinguishes a transient throttle (succeeds on
    retry) from a real miss (keeps failing, reported as such at the end)."""
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
    """One lock per distinct accession, so if the same accession appears on
    multiple FASTA records it's only ever fetched from NCBI once, and
    concurrent workers don't race to write the same cache file."""
    with _accession_locks_guard:
        lock = _accession_locks.get(accession)
        if lock is None:
            lock = threading.Lock()
            _accession_locks[accession] = lock
        return lock


def resolve_one(vjid, accession, cache_dir, use_cache, retries, retry_backoff):
    """Thread-worker: resolve+fetch a single accession, using the on-disk
    cache. Returns (vjid, accession, version, report_or_None, was_cached, error_str).
    report is None when the accession resolved fine but NCBI Datasets has no
    virus-genome report for it -- that is NOT an error."""
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


# --------------------------------------------------------------------------
# Elasticsearch helpers
# --------------------------------------------------------------------------

def build_session(token: str) -> requests.Session:
    """ES API-key auth is ALWAYS 'Authorization: ApiKey <base64(id:key)>' --
    never Bearer, never Basic. Set it directly as a session header (rather
    than via requests.auth.AuthBase) so there is no ambiguity about what is
    actually sent on the wire."""
    session = requests.Session()
    session.headers.update({"Authorization": f"ApiKey {token}"})
    return session


def check_auth(session, es_url, index_name, verify_ssl):
    """Fail loudly on bad auth/wrong index instead of silently producing a
    0-matches result later on that looks identical to 'no records need
    updating'."""
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
    if count == 0:
        logger.warning(
            "Index '%s' reports 0 documents -- every vjid lookup will miss. "
            "Check --index / --es-url point at the right place.", index_name
        )


def is_empty(value):
    return value is None or value == "" or value == []


def find_existing_by_vjid(session, es_url, index_name, vjids, verify_ssl, batch_size=200):
    """Documents are looked up by field '0dn' (the vjid) via a `_search` with
    a `terms` query -- NOT by the Elasticsearch-internal _id, which is some
    other unrelated identifier (e.g. a legacy numeric id) that this script
    intentionally never relies on. Returns {vjid: source_or_None}.

    NOTE: '0dn' is confirmed mapped as 'keyword', so this terms query is
    fast and exact-match; no analyzed-text concerns here.
    """
    out = {v: None for v in vjids}
    url = f"{es_url.rstrip('/')}/{index_name}/_search"
    logger = logging.getLogger("update_vjdb")
    n_batches = (len(vjids) + batch_size - 1) // batch_size
    progress = tqdm(
        range(0, len(vjids), batch_size), total=n_batches,
        unit="batch", desc="Looking up existing docs", file=sys.__stderr__,
    )
    for start in progress:
        batch = vjids[start:start + batch_size]
        payload = {
            "query": {"terms": {VJID_FIELD: batch}},
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
                out[vjid] = source
        progress.set_postfix(found=sum(1 for s in out.values() if s is not None))
    return out


def _painless_set_script(doc: dict):
    """Build a Painless script + params that sets each field in `doc` on
    ctx._source. Values are passed as params (not interpolated into the
    script source) so this is safe regardless of value content/type."""
    lines = []
    params = {}
    for i, (field, value) in enumerate(doc.items()):
        pname = f"p{i}"
        lines.append(f"ctx._source['{field}'] = params.{pname};")
        params[pname] = value
    return " ".join(lines), params


def update_by_vjid(session, es_url, index_name, vjid, doc, verify_ssl, timeout=60):
    """Apply a partial update to the document(s) matching 0dn == vjid, via
    Elasticsearch's _update_by_query + a Painless script -- this NEVER uses
    the Elasticsearch-internal _id, only the vjid field, per project
    convention. Returns (updated_count, response_json)."""
    url = f"{es_url.rstrip('/')}/{index_name}/_update_by_query"
    script_source, params = _painless_set_script(doc)
    payload = {
        "query": {"term": {VJID_FIELD: vjid}},
        "script": {"source": script_source, "params": params, "lang": "painless"},
        "conflicts": "proceed",
    }
    resp = session.post(url, json=payload, verify=verify_ssl, timeout=timeout, params={"refresh": "false"})
    resp.raise_for_status()
    result = resp.json()
    return result.get("updated", 0), result


def bulk_update(session, es_url, index_name, actions, verify_ssl, workers=8):
    """actions: list of (vjid, doc_partial_dict). Each is applied with its
    own _update_by_query call matched on '0dn' (vjid), run across a thread
    pool since these are independent network round trips (like the NCBI
    resolution step) -- not a single ES _bulk request, since _bulk requires
    a concrete _id and this project intentionally matches by vjid instead."""
    total_ok = 0
    failures = []
    logger = logging.getLogger("update_vjdb")
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        future_to_vjid = {
            executor.submit(update_by_vjid, session, es_url, index_name, vjid, doc, verify_ssl): vjid
            for vjid, doc in actions
        }
        progress = tqdm(
            concurrent.futures.as_completed(future_to_vjid), total=len(future_to_vjid),
            unit="doc", desc="Writing updates", file=sys.__stderr__,
        )
        for fut in progress:
            vjid = future_to_vjid[fut]
            try:
                updated, result = fut.result()
                if updated == 0:
                    failures.append((vjid, "matched 0 documents", result))
                elif updated > 1:
                    # Should not happen if 0dn is unique, but don't silently
                    # assume -- flag it so it can be investigated rather than
                    # trusting a possibly-wrong update fanned out to duplicates.
                    logger.warning("vjid %s: _update_by_query updated %d documents (expected 1)", vjid, updated)
                    total_ok += 1
                else:
                    total_ok += 1
            except requests.exceptions.RequestException as e:
                failures.append((vjid, str(e), None))
            progress.set_postfix(ok=total_ok, failed=len(failures))
    return total_ok, failures


# --------------------------------------------------------------------------
# Per-record payload building
# --------------------------------------------------------------------------

def build_update_doc(normalized, resolved_version, existing_source, schema, vmr_version, line_label):
    """Returns (doc_partial, errors, warnings).

    `normalized` is {} when no NCBI Datasets report was available for this
    accession -- in that case only the ALWAYS-SET fields below get written,
    plus insdc_accession is filled if empty; no other content comparison
    happens (there is nothing to compare against)."""
    errors, warnings = [], []
    existing_source = existing_source or {}
    doc = {}

    # --- ALWAYS-SET fields -------------------------------------------------
    doc[VMR_VERSION_FIELD] = vmr_version

    is_annotated_refseq = normalized.get(REFSEQ_SOURCE_KEY, False)
    is_refseq = bool(is_annotated_refseq) if isinstance(is_annotated_refseq, (bool, str)) else False
    if isinstance(is_annotated_refseq, str):
        is_refseq = is_annotated_refseq.strip().lower() == "true"
    for field_id, value in (
        (IS_REFSEQ_FIELD_ID, is_refseq),
        (IS_GENBANK_FIELD_ID, True),
    ):
        entry = schema.get(field_id)
        if entry is None:
            errors.append(f"{line_label}: schema field '{field_id}' not found in schema")
            continue
        doc[entry["index"]] = value

    dk_entry = schema.get("sequence_source")
    if dk_entry is None:
        errors.append(f"{line_label}: schema field 'sequence_source' (0dk) not found in schema")
    else:
        existing_dk = existing_source.get(dk_entry["index"])
        if existing_dk is None:
            dk_list = []
        elif isinstance(existing_dk, list):
            dk_list = list(existing_dk)
        else:
            dk_list = [existing_dk]
        if GENBANK_SOURCE_VALUE not in dk_list:
            dk_list.append(GENBANK_SOURCE_VALUE)
        doc[dk_entry["index"]] = dk_list

    # --- FILL-IF-EMPTY fields -----------------------------------------------
    # Use the freshly-resolved, versioned accession as the candidate value for
    # insdc_accession (0aa) instead of whatever the report itself reports.
    # 0aa is stored as a LIST in the DB (e.g. ["AJ867584.2"]) -- wrap to match.
    fill_candidates = dict(normalized)
    fill_candidates.pop(REFSEQ_SOURCE_KEY, None)
    fill_candidates.pop("source_database", None)
    fill_candidates[INSDC_ACCESSION_KEY] = [resolved_version]

    for jsonl_field, value in fill_candidates.items():
        if jsonl_field in SILENT_DROP_FIELDS:
            continue

        if jsonl_field in SEMICOLON_LIST_FIELDS and isinstance(value, str):
            caster = SEMICOLON_LIST_FIELDS[jsonl_field]
            parts = [p.strip() for p in value.split(";") if p.strip() != ""]
            try:
                value = [caster(p) for p in parts]
            except ValueError:
                errors.append(
                    f"{line_label}: field '{jsonl_field}' could not split/cast "
                    f"';'-list value {value!r} to {caster.__name__}"
                )
                continue

        schema_field_id, index, vjdb_type = resolve_field(jsonl_field, schema)
        if index is None:
            # Closed-world policy: a field is only ever allowed through if it is
            # either mapped to a schema field_id (directly or via
            # FIELD_ID_OVERRIDES) or explicitly listed in SILENT_DROP_FIELDS
            # above. Anything else is an unrecognized/new field (e.g. a field
            # NCBI Datasets started returning that nobody has decided a policy
            # for yet) -- this is a HARD ERROR that aborts the dry run rather
            # than silently skipping it, so nothing not in the mapping can ever
            # be written to Elasticsearch "on the fly".
            errors.append(
                f"{line_label}: field '{jsonl_field}' not found in schema and has "
                f"no override, and is not in SILENT_DROP_FIELDS -- unmapped. "
                f"Add it to FIELD_ID_OVERRIDES or SILENT_DROP_FIELDS explicitly."
            )
            continue

        # fill-if-empty: only touch fields missing/empty in the existing doc
        if not is_empty(existing_source.get(index)):
            continue
        # a field already staged above (always-set fields) is not overwritten here
        if index in doc:
            continue

        try:
            coerced, coerce_warnings = coerce_value(
                value, vjdb_type, f"{line_label} field '{jsonl_field}' ({index})"
            )
            warnings.extend(coerce_warnings)
        except ValueError as e:
            errors.append(str(e))
            continue

        doc[index] = coerced

    return doc, errors, warnings


# --------------------------------------------------------------------------
# Logging tee
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


# --------------------------------------------------------------------------
# Main pipeline
# --------------------------------------------------------------------------

def run(args):
    global _NCBI_API_KEY
    _NCBI_API_KEY = args.ncbi_api_key or os.environ.get("NCBI_API_KEY") or os.environ.get("API_KEY")

    logger = logging.getLogger("update_vjdb")
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s", "%Y-%m-%d %H:%M:%S"))
    logger.addHandler(handler)

    if _NCBI_API_KEY:
        logger.info("NCBI API key detected -- using it (~10 req/sec limit instead of ~3 req/sec).")
    else:
        logger.info(
            "No NCBI API key found (checked --ncbi-api-key, $NCBI_API_KEY, $API_KEY) -- "
            "rate-limited to ~3 req/sec."
        )

    global _eutils_limiter, _datasets_limiter
    default_rate = 8.0 if _NCBI_API_KEY else 2.5
    eutils_rate = args.eutils_rate if args.eutils_rate is not None else default_rate
    datasets_rate = args.datasets_rate if args.datasets_rate is not None else default_rate
    _eutils_limiter = RateLimiter(eutils_rate)
    _datasets_limiter = RateLimiter(datasets_rate)
    logger.info(
        "Rate limits: eutils=%.1f req/sec, datasets=%.1f req/sec (self-throttled "
        "globally, independent of --workers)", eutils_rate, datasets_rate,
    )
    logger.info(
        "Each new (non-cached) accession costs 2-3 requests (esearch+efetch[+datasets]), "
        "so expect roughly min(eutils_rate/2, datasets_rate) accessions/sec of real "
        "throughput once the cache is warm."
    )

    fasta_records = parse_fasta_headers(args.fasta)
    logger.info("Parsed %d header(s) from %s", len(fasta_records), args.fasta)

    if args.vmr_version is None:
        args.vmr_version = vmr_version_from_filename(args.fasta)
        logger.info("VMR version (from filename): %s", args.vmr_version)
    else:
        logger.info("VMR version (from --vmr-version): %s", args.vmr_version)

    schema = load_schema(schema_path=args.schema, schema_url=args.schema_url)

    cache_dir = Path(args.cache_dir)
    use_cache = not args.no_cache

    verify_ssl = not args.no_verify_ssl
    session = build_session(args.token)
    check_auth(session, args.es_url, args.index, verify_ssl)

    # --- Step 1: resolve accessions + fetch NCBI Datasets reports ----------
    # Each iteration shells out (esearch/efetch/xtract, then `datasets summary`),
    # which is dominated by network latency to NCBI, not CPU -- so this is run
    # across a thread pool (--workers, default 8) instead of one-at-a-time.
    resolved = {}       # vjid -> (accession_no_version, version, report_or_None)
    resolve_errors = []
    no_report_count = 0
    cache_hits = 0
    logger.info(
        "Resolving %d accession(s) via NCBI Datasets using %d worker thread(s) "
        "(this is the slow step)", len(fasta_records), args.workers,
    )
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [
            executor.submit(
                resolve_one, rec["vjid"], rec["accession"], cache_dir, use_cache,
                args.retries, args.retry_backoff,
            )
            for rec in fasta_records
        ]
        progress = tqdm(
            concurrent.futures.as_completed(futures),
            total=len(futures), unit="acc", desc="Resolving accessions",
            file=sys.__stderr__,
        )
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
    if cache_hits:
        logger.info("Served %d/%d accession(s) from cache", cache_hits, len(fasta_records))
    if no_report_count:
        logger.info(
            "%d accession(s) resolved fine but have no NCBI Datasets virus-genome "
            "report -- these will still be updated with version/source/booleans "
            "only, no content comparison.", no_report_count,
        )

    non_cached_total = len(fasta_records) - cache_hits
    if non_cached_total > 0 and len(resolve_errors) / non_cached_total > 0.10:
        logger.warning(
            "%d/%d (%.0f%%) of non-cached lookups failed even after retries -- "
            "this is higher than expected and may mean the rate limit is still "
            "too aggressive for this network, or NCBI is having issues. Consider "
            "lowering --eutils-rate/--datasets-rate and/or raising --retries.",
            len(resolve_errors), non_cached_total,
            100.0 * len(resolve_errors) / non_cached_total,
        )

    if resolve_errors:
        logger.warning("%d accession(s) could not be resolved at all (esearch/efetch failure):", len(resolve_errors))
        for e in resolve_errors[:20]:
            logger.warning("  %s", e)

    # --- Step 2: fetch existing documents from Elasticsearch (by 0dn, not _id) ----
    all_vjids = list(resolved.keys())
    logger.info("Looking up %d existing document(s) in %s/%s (by field '%s')",
                len(all_vjids), args.es_url, args.index, VJID_FIELD)
    lookup = find_existing_by_vjid(session, args.es_url, args.index, all_vjids, verify_ssl, batch_size=args.batch_size)
    missing_vjids = [v for v, src in lookup.items() if src is None]
    if missing_vjids:
        logger.warning(
            "%d vjid(s) not found in the index -- these will be SKIPPED (never created):",
            len(missing_vjids),
        )
        for v in missing_vjids[:20]:
            logger.warning("  %s", v)

    # --- Step 3: build per-record update payloads ---------------------------
    all_errors, all_warnings = [], []
    actions = []
    for vjid, (accession, version, report) in resolved.items():
        existing_source = lookup.get(vjid)
        if existing_source is None:
            continue  # not in db -- never create
        # No NCBI Datasets report available: skip content comparison entirely --
        # only the always-set fields (version, source append, refseq/genbank
        # booleans) plus the resolved accession itself get written.
        normalized = normalize_report(report) if report is not None else {}
        line_label = f"vjid {vjid} ({accession})"
        doc, errors, warnings = build_update_doc(
            normalized, version, existing_source, schema, args.vmr_version, line_label
        )
        all_errors.extend(errors)
        all_warnings.extend(warnings)
        actions.append((vjid, doc))

    # --- Report ---------------------------------------------------------------
    print(f"\nParsed:              {len(fasta_records)}")
    print(f"Resolved via NCBI:   {len(resolved)}  (failed: {len(resolve_errors)}, no-report: {no_report_count})")
    print(f"Found in index:      {len(actions)}  (missing/skipped: {len(missing_vjids)})")
    print(f"Errors:              {len(all_errors)}")
    print(f"Warnings:            {len(all_warnings)}")

    if all_warnings:
        print("\n--- Sample warnings (first 10) ---")
        for w in all_warnings[:10]:
            print("  WARN:", w)
        if len(all_warnings) > 10:
            print(f"  ... and {len(all_warnings) - 10} more")

    if all_errors:
        print("\n--- Errors (first 20) ---")
        for e in all_errors[:20]:
            print("  ERROR:", e)
        if len(all_errors) > 20:
            print(f"  ... and {len(all_errors) - 20} more")

    if actions:
        print("\n--- Sample update payload(s) ---")
        for vjid, doc in actions[:3]:
            print(f"  vjid={vjid}:")
            print(json.dumps(doc, indent=2, ensure_ascii=False))

    if all_errors:
        print(f"\nDry run FAILED with {len(all_errors)} error(s). Aborting -- nothing was updated.")
        sys.exit(1)

    print("\nDry run PASSED (no hard errors).")

    if not args.commit:
        print(f"Run again with --commit to actually update these {len(actions)} record(s) in '{args.index}'.")
        return

    print(f"\n=== UPDATING {len(actions)} record(s) in {args.es_url}/{args.index} (matched by '{VJID_FIELD}') ===")
    updated, failures = bulk_update(
        session, args.es_url, args.index, actions, verify_ssl, workers=args.workers
    )
    print(f"\nUpdated: {updated}/{len(actions)}")
    if failures:
        print(f"Failures: {len(failures)}")
        for f in failures[:20]:
            print("  ", f)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True, help="Multi-FASTA file with '>vjid | accession | ...' headers")
    parser.add_argument("--schema", default=None, help="Local DB_Scheme.json (default: fetch fresh from GitHub)")
    parser.add_argument("--schema-url", default=SCHEMA_URL_DEFAULT)
    parser.add_argument("--es-url", default=ES_URL_DEFAULT)
    parser.add_argument("--index", default=INDEX_NAME_DEFAULT)
    parser.add_argument("--token", default=AUTH_TOKEN_DEFAULT)
    parser.add_argument("--no-verify-ssl", action="store_true")
    parser.add_argument("--commit", action="store_true", help="Perform the real partial updates if the dry run is clean")
    parser.add_argument("--batch-size", type=int, default=200)
    parser.add_argument("--cache-dir", default=".cache/reports")
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument(
        "--workers", type=int, default=8,
        help=(
            "number of parallel threads used to resolve/fetch accessions from "
            "NCBI (default: 8). Actual request pacing is controlled separately "
            "by --eutils-rate/--datasets-rate, NOT by this number -- --workers just "
            "needs to be at least high enough to keep the rate limiter fed (e.g. "
            "an in-flight request taking ~0.3-1s round trip means a handful of "
            "workers is already enough to saturate an 8 req/sec limiter). Raising "
            "this well past that won't go faster and just adds idle/blocked threads."
        ),
    )
    parser.add_argument(
        "--ncbi-api-key",
        default=None,
        help=(
            "NCBI API key. Passed to the 'datasets' CLI via "
            "--api-key and exported as the NCBI_API_KEY environment "
            "variable for EDirect (esearch/efetch). Defaults to "
            "$NCBI_API_KEY or $API_KEY if set."
        ),
    )
    parser.add_argument(
        "--eutils-rate", type=float, default=None,
        help="max esearch/efetch requests per second, globally across all workers "
             "(default: 8.0 with an API key, 2.5 without)",
    )
    parser.add_argument(
        "--datasets-rate", type=float, default=None,
        help="max `datasets summary` requests per second, globally across all "
             "workers (default: 8.0 with an API key, 2.5 without)",
    )
    parser.add_argument(
        "--retries", type=int, default=3,
        help="retry attempts for an accession before giving up (distinguishes a "
             "transient throttle from a genuine miss; default: 3)",
    )
    parser.add_argument(
        "--retry-backoff", type=float, default=1.0,
        help="base seconds to wait before a retry, doubled each attempt (default: 1.0)",
    )
    parser.add_argument(
        "--vmr-version", default=None,
        help="ictv_vmr_version (0bg) to write. Defaults to the value derived "
             "from --fasta's filename, e.g. 'VMR_MSL41.v1.20260320.fasta' -> "
             "'VMR_MSL41.v1.20260320' (same convention as "
             "i.insert_jsonl_to_vjdb.py).",
    )
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
