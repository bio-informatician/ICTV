#!/usr/bin/env python3
"""
m.upsert_vjdb.py

Merged replacement for the i / j / l pipeline. This is l.update_and_insert_w_genome.py
(the well-structured "current" script) with the missing "brand new accession"
path folded back in from i.insert_jsonl_to_vjdb.py, and the existing-record
lookup switched from vjid-keyed to accession-keyed, per your instructions:

    "The input will have records that may exist in database or may not.
     For a new record an insertion event happens (vjid generated) followed
     by adding the genome sequence, with NO change log.
     If it already exists in the database (matched by ACCESSION NUMBER, not
     virjendb id) and the resolved version is higher than what's stored (or
     the stored accession has no version at all), do a VERSION-DIFF using
     the SAME virjendb id, per the rules already defined in l.py.
     If the version is the same, also follow l.py's existing rules
     (content-update / fill-only insertion)."

-------------------------------------------------------------------------
WHAT CHANGED VS l.update_and_insert_w_genome.py
-------------------------------------------------------------------------
1. INPUT FORMAT
   FASTA headers no longer need to carry a vjid. The record is now looked
   up purely by accession:

       >AB000403 | organism=Tenuivirus oryzabrevis | length=2704 | description=...
       ACGTACGT...

   Backward compatibility: if a header still uses the old
   ">vjid | accession | ..." shape (first field matches ``vj\\d+``), the
   vjid field is parsed but IGNORED for lookup purposes (only logged) --
   matching is always done by accession, never by a vjid claimed in the
   input, per your instruction.

2. LOOKUP: ACCESSION-KEYED, NOT VJID-KEYED, CASE/WHITESPACE NORMALIZED
   find_active_by_accession() replaces find_active_by_vjid(). It searches
   the metadata index for an ACTIVE (0hd=true) document whose 0aa list
   contains either the bare (unversioned) accession or accession+".<any
   version>". This correctly matches the "stored accession has no version
   at all" case you called out, in addition to ordinary version bumps.
   Every accession is passed through normalize_accession() (trim + upper)
   BEFORE it is ever compared, cached, or queried, and the ES query itself
   uses case_insensitive wildcard matching as a second line of defense --
   this closes the "record silently treated as brand-new because casing
   differed" bug that caused a duplicate active record for AY038069.

3. NEW EVENT: BRAND-NEW INSERTION (accession not found in DB at all)
   Ported from i.insert_jsonl_to_vjdb.py:
     - fetch_latest_accession(): looks up the current max 0dn across BOTH
       metadata_index AND genome_sequences (not metadata alone), so a
       vjid already used in either index -- even if the other write never
       happened because of a mid-run crash -- is never reissued to an
       unrelated accession.
     - make_accession_counter(): hands out vj000000000001-style ids.
   A brand-new document is built (0dn = freshly minted vjid, 0do=1, 0cx=1,
   0hd=true, 0fh=today, 0fg=today (insertion_date, NEW/immutable field,
   never touched again after this), 0bg=vmr_version, 0ai/0ah recomputed,
   0dk=["GenBank"], 0aa=[resolved_version], plus every other field the
   NCBI report / FASTA header supplies).
   This event NEVER writes to change_log (per your instruction), but DOES
   push the sequence to genome_sequences (0cx=1), same field shape as the
   version-diff case.

4. GENOME-GATED WRITES (brand-new-insertion + version-diff only)
   These are the only two event types that ever produce a genome_doc.
   Their metadata/change_log actions are now held back and only released
   AFTER the genome bulk write's PER-ITEM result (not just the aggregate
   "errors" flag) confirms that specific record's genome write actually
   succeeded. This prevents the metadata/genome desync where an active
   metadata record could exist while its sequence silently failed to
   write, or vice versa.

5. NCBI DATASETS "NO VIRUS REPORT" FALLBACK
   If `datasets summary virus genome accession <version>` returns zero
   results (e.g. the accession isn't classified as a virus genome by NCBI
   Datasets, or hasn't been indexed there yet -- see Z27119.1), the code
   no longer silently inserts a bare-header stub NOR silently skips it.
   It falls back to fetching the plain GenBank flat file for that SAME
   resolved, versioned accession via efetch/nuccore, extracts organism /
   length / country via simple regex, and proceeds with THAT data --
   always flagged with a warning for manual review. If the fallback also
   yields nothing, the accession is a hard error (skipped), never a
   silent stub insert.

6. EXISTING-RECORD EVENTS: UNCHANGED FROM l.py (aside from genome-gating
   above)
   Once an active document is found by accession, the vjid is read off
   that document and everything downstream (version-diff / content-update
   / fill-only insertion / no-change, change_log writing, bookkeeping
   fields) is exactly build_record_actions() from l.py, untouched. The
   "stored accession has no version" case is already handled by that
   function's existing comparison:
       strip_version(old_acc) == strip_version(resolved_version) and old_acc != resolved_version
   which is True whenever old_acc has no version suffix at all (its
   stripped form still equals its own base), so no logic change was
   needed there -- only the *lookup* needed to change.

-------------------------------------------------------------------------
ASSUMPTIONS FLAGGED FOR CONFIRMATION
-------------------------------------------------------------------------
  - Header "organism="/"length=" extras are mapped onto the SAME
    normalized-report field names used elsewhere (virus_organism_name,
    sequence_length) and used ONLY to fill gaps when no NCBI Datasets
    report is available (or the report is missing that value). Report
    data always wins when both exist. "description=" has no DB column
    (same as in the original scripts) and is dropped.
  - The GenBank flat-file fallback (see #5 above) is a LESSER source than
    a real NCBI Datasets virus report: no taxonomy, no host lineage, no
    isolate data. It is only ever used when the real report is absent,
    and its use is always flagged with a warning.
  - vjid minting reuses ACCESSION_PREFIX/ACCESSION_DIGITS = "vj"/12 as in
    i.py. These are read from config.yml if present
    (values.accession_prefix / values.accession_digits), else fall back
    to "vj"/12 with a printed warning -- add them to config.yml to silence
    the warning.
  - insertion_date field id (0fg) is read from config.yml
    (field_ids.insertion_date_field) if present, else falls back to "0fg"
    with a printed warning.
  - A batch of accession lookups is done as ONE bool/should query per
    batch (term for bare accession OR wildcard "<accession>.*" for any
    versioned form), not one query per accession, to keep this reasonably
    fast on large inputs. Wildcard "ACC.*" cannot cross-match a
    different-prefixed accession like "ACC10" because the pattern
    requires a literal "." right after "ACC". Wildcard queries use
    "case_insensitive": true (requires Elasticsearch >= 7.10 -- flag if
    your cluster is older, a script_score/match fallback would be needed
    instead).
  - If two FASTA records resolve to the SAME base accession (shouldn't
    happen in a well-formed VMR export) the second one to be processed
    will see the first one's freshly-written state only if you re-run
    with an updated lookup -- this script, like l.py, resolves and looks
    up existing docs in one batch pass, so such intra-run collisions are
    NOT specially handled. Flagging in case this matters for your data.
  - MAX_CHUNK_BYTES (see --max-chunk-bytes) now defaults to 900000, well
    under nginx's common 1MB client_max_body_size default, to avoid the
    wasted 413-retry churn seen when the old 8MB default was miles above
    the server's real limit. A single oversized genome document (bigger
    than this cap on its own) is now excluded from the run at the
    pre-flight stage instead of only surfacing as a crash at --commit
    time -- see the "genome" pre-flight check in run().

-------------------------------------------------------------------------
USAGE
-------------------------------------------------------------------------
    # Dry run only (default, safe, no network writes)
    python3 m.upsert_vjdb.py --fasta VMR_MSL41.v1.20260320.fasta

    # Dry run + real writes if dry run is clean
    python3 m.upsert_vjdb.py --fasta VMR_MSL41.v1.20260320.fasta --commit
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
import yaml
from tqdm import tqdm

import vmr_ictv_merge

# --------------------------------------------------------------------------
# Settings (config.yml) loading -- unchanged from l.py
# --------------------------------------------------------------------------

SETTINGS_PATH_DEFAULT = "config.yml"


def load_settings(path=None):
    resolved_path = path or os.environ.get("VJDB_SETTINGS_PATH") or SETTINGS_PATH_DEFAULT
    with open(resolved_path, "r", encoding="utf-8") as f:
        settings = yaml.safe_load(f)
    return settings


_SETTINGS = load_settings()


def _settings_get(path_str, default, warn_label):
    """Nested dict lookup with a dotted path, falling back to `default`
    (and printing a one-time warning) if the key is absent from
    config.yml. Used for the handful of NEW keys this merged script needs
    that config.yml (written for l.py) may not yet contain."""
    node = _SETTINGS
    for part in path_str.split("."):
        if not isinstance(node, dict) or part not in node:
            print(
                f"WARNING: config.yml missing '{path_str}' -- falling back to "
                f"default {default!r} for {warn_label}. Add it to config.yml "
                f"to silence this warning.",
                file=sys.stderr,
            )
            return default
        node = node[part]
    return node


# --------------------------------------------------------------------------
# Fixed / constant configuration (values sourced from config.yml, same as l.py)
# --------------------------------------------------------------------------

ES_URL_DEFAULT = _SETTINGS["elasticsearch"]["es_url"]
INDEX_NAME_DEFAULT = _SETTINGS["elasticsearch"]["index_name"]
GENOME_INDEX_DEFAULT = _SETTINGS["elasticsearch"]["genome_index"]
CHANGE_LOG_INDEX_DEFAULT = _SETTINGS["elasticsearch"]["change_log_index"]
AUTH_TOKEN_DEFAULT = _SETTINGS["elasticsearch"]["auth_token"]

SCHEMA_URL_DEFAULT = _SETTINGS["schema"]["schema_url"]

VJID_FIELD = _SETTINGS["field_ids"]["vjid_field"]
VERSION_FIELD = _SETTINGS["field_ids"]["version_field"]
SEQ_VERSION_FIELD = _SETTINGS["field_ids"]["seq_version_field"]
ACTIVE_FLAG_FIELD = _SETTINGS["field_ids"]["active_flag_field"]
RENEWED_DATE_FIELD = _SETTINGS["field_ids"]["renewed_date_field"]

VMR_VERSION_FIELD = _SETTINGS["field_ids"]["vmr_version_field"]
GENBANK_SOURCE_VALUE = _SETTINGS["values"]["genbank_source_value"]

IS_REFSEQ_FIELD_ID = _SETTINGS["field_ids"]["is_refseq_field_id"]
IS_GENBANK_FIELD_ID = _SETTINGS["field_ids"]["is_genbank_field_id"]
REFSEQ_SOURCE_KEY = _SETTINGS["field_ids"]["refseq_source_key"]

INSDC_ACCESSION_KEY = _SETTINGS["field_ids"]["insdc_accession_key"]
INSDC_ACCESSION_FIELD_ID = _SETTINGS["field_ids"]["insdc_accession_field_id"]

FIELD_ID_OVERRIDES = _SETTINGS["field_id_overrides"]
SILENT_DROP_FIELDS = set(_SETTINGS["silent_drop_fields"])

_CASTER_BY_NAME = {"int": int, "str": str}
SEMICOLON_LIST_FIELDS = {
    field: _CASTER_BY_NAME[caster_name]
    for field, caster_name in _SETTINGS["semicolon_list_fields"].items()
}

REPORT_SCHEMA_MAP = _SETTINGS["report_schema_map"]
TAXONOMIC_SUFFIXES = tuple(tuple(pair) for pair in _SETTINGS["taxonomic_suffixes"])

# --- NEW constants needed for the brand-new-insertion path (ported from i.py) ---
INSERTION_DATE_FIELD = _settings_get(
    "field_ids.insertion_date_field", "0fg", "INSERTION_DATE_FIELD"
)
ACCESSION_PREFIX = _settings_get("values.accession_prefix", "vj", "ACCESSION_PREFIX")
ACCESSION_DIGITS = _settings_get("values.accession_digits", 12, "ACCESSION_DIGITS")
# Only used if --offline-latest-accession is passed; otherwise the real
# latest vjid is always looked up live from the index (fetch_latest_accession).
LATEST_KNOWN_ACCESSION_FALLBACK = _settings_get(
    "values.latest_known_accession_fallback", "vj000000000000",
    "LATEST_KNOWN_ACCESSION_FALLBACK",
)

# Schema field_ids whose value should be compared CASE-INSENSITIVELY when
# deciding whether a content-update event has actually occurred. A change
# that's only a casing difference (e.g. genome_completeness "COMPLETE" vs
# "complete") is not a real content change and should not be logged to
# change_log or bump the metadata version. Configurable via config.yml
# (values.case_insensitive_fields: [...]); falls back to just
# genome_completeness if not present.
CASE_INSENSITIVE_FIELDS = set(_settings_get(
    "values.case_insensitive_fields", ["genome_completeness"],
    "CASE_INSENSITIVE_FIELDS",
))


# --------------------------------------------------------------------------
# Config (method/reason per event type) -- unchanged from l.py
# --------------------------------------------------------------------------

DEFAULT_CONFIG = _SETTINGS["change_log_default_config"]


def load_config(config_path):
    if not config_path or not os.path.exists(config_path):
        return deepcopy(DEFAULT_CONFIG)
    with open(config_path, "r", encoding="utf-8") as f:
        cfg = json.load(f)
    for key in ("version_change", "content_update"):
        if key not in cfg or "method" not in cfg[key] or "reason" not in cfg[key]:
            raise ValueError(f"config file missing required '{key}.method' / '{key}.reason'")
    return cfg


# --------------------------------------------------------------------------
# schema.py logic, inlined -- unchanged from l.py
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


def header_extras_as_normalized_fields(rec):
    """Map the FASTA header's 'organism=' / 'length=' extras onto the same
    target field names normalize_report() would use, so they can serve as
    a fallback when no NCBI Datasets report is available (or the report
    doesn't supply that particular value). 'description=' has no DB column
    and is intentionally dropped, same as in the original scripts."""
    out = {}
    if rec.get("organism"):
        out["virus_organism_name"] = rec["organism"]
    if rec.get("length"):
        try:
            out["sequence_length"] = int(rec["length"])
        except (TypeError, ValueError):
            out["sequence_length"] = rec["length"]
    return out


def merge_normalized(report_normalized, header_normalized, fallback_normalized=None):
    """Priority, weakest to strongest: header extras < GenBank flat-file
    fallback < real NCBI Datasets virus report. Report data always wins
    when present; the fallback only fills gaps a real report would have
    covered, and only when there IS no real report."""
    merged = dict(header_normalized)
    if fallback_normalized:
        merged.update(fallback_normalized)
    merged.update(report_normalized)
    return merged


# --------------------------------------------------------------------------
# type coercion -- unchanged from l.py
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
            f"double-check the schema source if this warning appears.",
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


def normalize_accession(acc: str) -> str:
    """Canonical form used for EVERY accession comparison/lookup/cache-key
    in this pipeline -- trims stray whitespace and forces uppercase.
    NCBI accessions are canonically uppercase; any drift here (a lowercase
    letter from a hand-edited FASTA header, a stray tab, etc.) silently
    breaks term/wildcard matching against 0aa and causes an existing
    record to be misclassified as brand-new -- which is exactly what
    caused a duplicate active record for AY038069.2."""
    return acc.strip().upper()


def strip_version(accession: str) -> str:
    accession = normalize_accession(accession)
    return accession.rsplit(".", 1)[0] if "." in accession else accession


def sha256_hex(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def to_bool(value):
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() == "true"
    return bool(value)


def values_equal(schema_field_id, old_val, new_val):
    """Equality check used to decide whether a stored field genuinely
    changed. For fields listed in CASE_INSENSITIVE_FIELDS, a difference
    that is ONLY casing (e.g. "COMPLETE" vs "complete") does not count as
    a change. For everything else this is a plain equality check, same as
    the original l.py behaviour."""
    if schema_field_id not in CASE_INSENSITIVE_FIELDS:
        return old_val == new_val

    def norm(v):
        if isinstance(v, str):
            return v.strip().lower()
        if isinstance(v, list):
            return [x.strip().lower() if isinstance(x, str) else x for x in v]
        return v

    return norm(old_val) == norm(new_val)


# --------------------------------------------------------------------------
# vjid minting -- ported from i.insert_jsonl_to_vjdb.py
# --------------------------------------------------------------------------

def fetch_latest_accession(session, es_url, index_name, genome_index_name, verify_ssl):
    """Take the max 0dn seen in EITHER metadata_index OR genome_sequences,
    so a vjid that's already used in one (even if the other write never
    happened, e.g. due to a crash between the two) is never reissued."""
    def _max_dn(idx):
        url = f"{es_url.rstrip('/')}/{idx}/_search"
        query = {"size": 1, "sort": [{VJID_FIELD: {"order": "desc"}}], "_source": [VJID_FIELD]}
        resp = session.post(url, json=query, verify=verify_ssl, timeout=30)
        resp.raise_for_status()
        hits = resp.json().get("hits", {}).get("hits", [])
        return hits[0]["_source"][VJID_FIELD] if hits else None

    candidates = [d for d in (_max_dn(index_name), _max_dn(genome_index_name)) if d]
    if not candidates:
        raise RuntimeError(
            "No documents found in metadata or genome index -- cannot determine "
            "vjid counter start. Pass --offline-latest-accession to override."
        )
    latest = max(candidates)  # string comparison works since vjid is zero-padded
    print(f"Latest vjid across metadata+genome indices: {latest}")
    return latest


def make_accession_counter(latest_known):
    numeric_part = int(latest_known[len(ACCESSION_PREFIX):])
    counter = {"n": numeric_part}

    def next_id():
        counter["n"] += 1
        return f"{ACCESSION_PREFIX}{counter['n']:0{ACCESSION_DIGITS}d}"

    return next_id


# --------------------------------------------------------------------------
# FASTA parsing -- accession-first, vjid optional/ignored (NEW)
# --------------------------------------------------------------------------

_LEGACY_VJID_RE = re.compile(r"^vj\d+$", re.IGNORECASE)


def vmr_version_from_filename(input_path):
    name = os.path.basename(input_path)
    for ext in (".fasta", ".txt", ".json", ".jsonl"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def parse_fasta(fasta_path):
    """Yield {'accession', 'organism', 'length', 'description', 'sequence',
    'legacy_vjid'} for each '>' header.

    Primary (expected) header shape:
        >ACCESSION | organism=... | length=... | description=...

    Legacy shape (still parsed, vjid captured but NEVER used for lookup --
    only logged for visibility):
        >vjNNNNNNNNNNNN | ACCESSION | organism=... | length=... | description=...

    The accession is passed through normalize_accession() (trim + upper)
    the moment it's parsed, so every downstream comparison, cache key, and
    ES lookup is already working from a canonical form.
    """
    logger = logging.getLogger("upsert_vjdb")
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
                if not parts or not parts[0]:
                    logger.warning("line %d: empty header -- skipped: %r", line_no, header)
                    continue

                legacy_vjid = None
                if _LEGACY_VJID_RE.match(parts[0]) and len(parts) >= 2:
                    legacy_vjid = parts[0]
                    accession = normalize_accession(parts[1])
                    extra_parts = parts[2:]
                    logger.info(
                        "line %d: legacy 'vjid | accession' header detected "
                        "(vjid=%s) -- vjid IGNORED for lookup, matching by "
                        "accession=%s only", line_no, legacy_vjid, accession,
                    )
                else:
                    accession = normalize_accession(parts[0])
                    extra_parts = parts[1:]

                extra = {}
                for kv in extra_parts:
                    if "=" in kv:
                        k, v = kv.split("=", 1)
                        extra[k.strip()] = v.strip()

                current_rec = {
                    "accession": accession,
                    "organism": extra.get("organism"),
                    "length": extra.get("length"),
                    "description": extra.get("description"),
                    "legacy_vjid": legacy_vjid,
                }
            else:
                if current_rec is not None:
                    current_seq.append(line.strip())
        flush()

    return records


# --------------------------------------------------------------------------
# NCBI Datasets resolution -- unchanged from l.py, plus the GenBank
# flat-file fallback for accessions with no virus-genome report.
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


def fetch_genbank_flatfile_fallback(accession_version):
    """Fallback used ONLY when 'datasets summary virus genome' returned zero
    results for an accession (e.g. NCBI Datasets doesn't classify it as a
    virus genome, or hasn't indexed it there yet, but the accession still
    exists as an ordinary GenBank record -- e.g. Z27119.1).

    IMPORTANT: `accession_version` here is the RESOLVED, FULLY VERSIONED
    accession returned by latest_accession() (i.e. the SAME string that
    ends up in 0aa) -- never the raw, possibly version-less accession as
    typed in the input FASTA. This guarantees the flat file fetched here
    is for the exact sequence version being inserted, not just whatever
    version happened to be in the input file.

    Fetches the flat-file record directly via efetch/nuccore and extracts
    a SMALL, conservative set of fields using plain regex (no biopython
    dependency, matching the esearch/efetch/xtract style already used
    elsewhere in this script). Returns a normalized-field dict using ONLY
    target field names already proven valid elsewhere in this file
    (virus_organism_name, sequence_length, collection_geo_location) -- or
    {} if nothing could be parsed.

    This is intentionally a LESSER fallback: it can never contain taxonomy,
    host lineage, or isolate data the way a real Datasets report can, it
    never overrides an actual report, and its use is always flagged as a
    warning at the call site for manual review.
    """
    _eutils_limiter.wait()
    cmd = ["efetch", "-db", "nuccore", "-id", accession_version, "-format", "gb"]
    env = os.environ.copy()
    if _NCBI_API_KEY:
        env["NCBI_API_KEY"] = _NCBI_API_KEY

    p = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if p.returncode != 0 or not p.stdout.strip():
        return {}

    text = p.stdout
    out = {}

    m = re.search(r"^LOCUS\s+\S+\s+(\d+)\s+bp", text, re.MULTILINE)
    if m:
        out["sequence_length"] = int(m.group(1))

    m = re.search(r"^SOURCE\s+(.+)$", text, re.MULTILINE)
    if m:
        # Drop a trailing common name in parens, e.g. "(fruit fly)".
        out["virus_organism_name"] = re.sub(r"\s*\([^)]*\)\s*$", "", m.group(1).strip())

    m = re.search(r'/country="([^"]+)"', text)
    if m:
        out["collection_geo_location"] = m.group(1)

    return out


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
        data = json.load(f)
    # Backward compatible with caches written before the fallback field existed.
    data.setdefault("fallback_normalized", {})
    return data


def save_cached(cache_dir: Path, accession: str, version: str, report, fallback_normalized):
    cache_dir.mkdir(parents=True, exist_ok=True)
    with _cache_path(cache_dir, accession).open("w") as f:
        json.dump({"version": version, "report": report,
                    "fallback_normalized": fallback_normalized}, f)


def resolve_one(accession, cache_dir, use_cache, retries, retry_backoff):
    with _lock_for_accession(accession):
        try:
            cached = load_cached(cache_dir, accession) if use_cache else None
            if cached is not None:
                return (accession, cached["version"], cached["report"],
                        cached["fallback_normalized"], True, None)

            version, report = fetch_datasets_report_with_retries(accession, retries, retry_backoff)

            fallback_normalized = {}
            if report is None:
                # Only spend the extra efetch call when the virus-specific
                # report genuinely came back empty. Uses the RESOLVED,
                # versioned accession (see fetch_genbank_flatfile_fallback
                # docstring), not the raw input accession.
                fallback_normalized = fetch_genbank_flatfile_fallback(version)

            if use_cache:
                with _cache_write_lock:
                    save_cached(cache_dir, accession, version, report, fallback_normalized)
            return accession, version, report, fallback_normalized, False, None
        except FetchNotFound as e:
            return accession, None, None, {}, False, str(e)
        except Exception as e:
            return accession, None, None, {}, False, str(e)


# --------------------------------------------------------------------------
# Elasticsearch helpers
# --------------------------------------------------------------------------

def build_session(token: str) -> requests.Session:
    session = requests.Session()
    session.headers.update({"Authorization": f"ApiKey {token}"})
    return session


def check_auth(session, es_url, index_name, verify_ssl):
    logger = logging.getLogger("upsert_vjdb")
    url = f"{es_url.rstrip('/')}/{index_name}/_count"
    resp = session.get(url, verify=verify_ssl, timeout=30)
    if resp.status_code == 401:
        raise SystemExit(f"Auth failed (401) against {url} -- check the API key.")
    if resp.status_code == 403:
        raise SystemExit(f"Auth OK but forbidden (403) against {url} -- key lacks permission on this index.")
    resp.raise_for_status()
    count = resp.json().get("count", 0)
    logger.info("Auth OK -- %s/%s has %d document(s).", es_url, index_name, count)


def find_active_by_accession(session, es_url, index_name, schema, accessions_no_version,
                              verify_ssl, batch_size=100):
    """NEW (replaces find_active_by_vjid). Looks up the ACTIVE (0hd==true)
    document for each accession, matched against 0aa (stored as a list) --
    either a bare/unversioned entry equal to the accession, or a versioned
    entry "<accession>.<anything>". Returns
    {accession_no_version: (es_id, source) or (None, None)}.

    One ES query per batch: a bool/should of (term + wildcard) pairs, one
    pair per accession in the batch. Wildcard queries use
    "case_insensitive": true as a second line of defense (requires
    Elasticsearch >= 7.10) on top of normalize_accession() already having
    upper-cased everything on the input side -- so a stored value with any
    unexpected casing still matches. Hits are then matched back to the
    right query accession client-side by stripping the versions found in
    each hit's own 0aa list (also normalized before comparison).
    """
    insdc_entry = schema.get(INSDC_ACCESSION_FIELD_ID)
    if insdc_entry is None:
        raise RuntimeError(
            f"schema field '{INSDC_ACCESSION_FIELD_ID}' (0aa) not found in schema -- "
            "cannot look up existing records by accession."
        )
    acc_index = insdc_entry["index"]

    out = {a: (None, None) for a in accessions_no_version}
    url = f"{es_url.rstrip('/')}/{index_name}/_search"
    logger = logging.getLogger("upsert_vjdb")
    n_batches = (len(accessions_no_version) + batch_size - 1) // batch_size
    progress = tqdm(
        range(0, len(accessions_no_version), batch_size), total=n_batches,
        unit="batch", desc="Looking up existing accessions", file=sys.__stderr__,
    )
    for start in progress:
        batch = accessions_no_version[start:start + batch_size]
        should = []
        for acc in batch:
            should.append({"term": {acc_index: acc}})
            should.append({"wildcard": {acc_index: {"value": f"{acc}.*", "case_insensitive": True}}})
            should.append({"wildcard": {acc_index: {"value": acc, "case_insensitive": True}}})
        payload = {
            "query": {"bool": {"filter": [{"term": {ACTIVE_FLAG_FIELD: True}}],
                                "must": [{"bool": {"should": should, "minimum_should_match": 1}}]}},
            "size": len(batch) * 3,
            "_source": True,
        }
        try:
            resp = session.post(url, json=payload, verify=verify_ssl, timeout=30)
            resp.raise_for_status()
        except requests.exceptions.RequestException as e:
            logger.error("Lookup batch %d-%d failed: %s", start, start + len(batch), e)
            raise
        hits = resp.json().get("hits", {}).get("hits", [])
        batch_set = set(batch)
        for hit in hits:
            source = hit.get("_source", {})
            if not to_bool(source.get(ACTIVE_FLAG_FIELD)):
                continue
            stored = source.get(acc_index) or []
            if isinstance(stored, str):
                stored = [stored]
            for stored_acc in stored:
                base = strip_version(stored_acc)
                if base in batch_set:
                    if out[base][0] is not None and out[base][0] != hit["_id"]:
                        logger.warning(
                            "accession %s matched more than one active document "
                            "(%s and %s) -- keeping the first one seen",
                            base, out[base][0], hit["_id"],
                        )
                        continue
                    out[base] = (hit["_id"], source)
        progress.set_postfix(found=sum(1 for eid, _ in out.values() if eid is not None))
    return out


def _iter_byte_aware_chunks(actions, chunk_size, max_chunk_bytes):
    """Yield sub-lists of `actions` (kept in whole action-pairs -- every
    op in this script is exactly 2 list elements: meta line + doc/source
    line) that respect BOTH a max element count (chunk_size) AND a max
    serialized NDJSON byte size (max_chunk_bytes)."""
    chunk = []
    chunk_bytes = 0
    # step=2 because every action is a (meta, doc) pair -- never split a pair.
    for i in range(0, len(actions), 2):
        pair = actions[i:i + 2]
        pair_bytes = sum(len(json.dumps(a).encode("utf-8")) + 1 for a in pair)  # +1 per newline

        # A single pair larger than the whole budget still has to go out
        # alone -- there's nothing smaller to split it into here.
        if chunk and (
            len(chunk) + len(pair) > chunk_size
            or chunk_bytes + pair_bytes > max_chunk_bytes
        ):
            yield chunk
            chunk, chunk_bytes = [], 0

        chunk.extend(pair)
        chunk_bytes += pair_bytes

    if chunk:
        yield chunk


def execute_bulk(session, es_url, index_name, actions, verify_ssl, timeout=180,
                  chunk_size=2000, max_chunk_bytes=900_000, _depth=0):
    """
    1. Chunks are capped by BOTH element count AND serialized byte size
       (max_chunk_bytes, default 900000 bytes -- safely under nginx's
       common 1MB client_max_body_size default; pass a larger value via
       --max-chunk-bytes once your proxy/ES limits are confirmed higher).
    2. On a 413 specifically, instead of letting raise_for_status() blow
       up the whole run, the offending chunk is halved and retried
       recursively (bounded depth) until it fits, rather than aborting.
    """
    if not actions:
        return {"errors": False, "items": []}

    logger = logging.getLogger("upsert_vjdb")
    url = f"{es_url.rstrip('/')}/{index_name}/_bulk"
    all_items = []
    has_errors = False

    for chunk in _iter_byte_aware_chunks(actions, chunk_size, max_chunk_bytes):
        payload = "\n".join(json.dumps(a) for a in chunk) + "\n"
        resp = session.post(
            url, data=payload.encode("utf-8"),
            headers={"Content-Type": "application/x-ndjson"},
            verify=verify_ssl, timeout=timeout,
        )

        if resp.status_code == 413:
            n_pairs = len(chunk) // 2
            if n_pairs <= 1 or _depth >= 8:
                raise RuntimeError(
                    f"413 Request Entity Too Large against {url} even after "
                    f"shrinking to {n_pairs} document(s) "
                    f"({len(payload.encode('utf-8'))} bytes). The proxy/server "
                    f"body-size limit is smaller than a single document can "
                    f"fit -- check nginx client_max_body_size / ES "
                    f"http.max_content_length."
                )
            half = (n_pairs // 2) * 2  # keep it even -- whole pairs only
            logger.warning(
                "413 from %s on a %d-document chunk (%d bytes) -- retrying "
                "as two smaller chunks of ~%d document(s) each.",
                url, n_pairs, len(payload.encode("utf-8")), half // 2,
            )
            sub_result = execute_bulk(
                session, es_url, index_name, chunk[:half], verify_ssl,
                timeout=timeout, chunk_size=chunk_size,
                max_chunk_bytes=max_chunk_bytes // 2, _depth=_depth + 1,
            )
            sub_result2 = execute_bulk(
                session, es_url, index_name, chunk[half:], verify_ssl,
                timeout=timeout, chunk_size=chunk_size,
                max_chunk_bytes=max_chunk_bytes // 2, _depth=_depth + 1,
            )
            has_errors = has_errors or sub_result["errors"] or sub_result2["errors"]
            all_items.extend(sub_result["items"])
            all_items.extend(sub_result2["items"])
            continue

        resp.raise_for_status()
        resp_json = resp.json()
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
# Shared field-mapping helper (used by both brand-new insertion and
# version-diff, factored out to avoid duplicating the loop)
# --------------------------------------------------------------------------

def map_normalized_fields(normalized, schema, line_label, errors, warnings):
    """Map a normalized {jsonl_field: value} dict onto {schema_index: coerced_value}.
    Appends to the given errors/warnings lists in place. Unmapped fields
    (not in schema, no override, not silently dropped) are HARD ERRORS,
    same closed-world policy as l.py."""
    doc = {}
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
        doc[idx] = coerced
    return doc


def compute_always_set_fields(normalized, schema, vmr_version, errors, line_label):
    """0ai/0ah (recomputed every time) + 0bg (vmr version). Shared by
    brand-new insertion, version-diff, and content-update/fill."""
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
    return always_set


# --------------------------------------------------------------------------
# Event builders
# --------------------------------------------------------------------------

def build_brand_new_insertion(vjid, normalized, resolved_version, sequence,
                                schema, vmr_version, line_label):
    """NEW EVENT (ported from i.insert_jsonl_to_vjdb.py + j.insert_sequences_to_vjdb.py):
    the accession does not exist in the metadata index at all. Mint a vjid,
    build a full fresh document, and push the sequence -- no change_log
    entry, ever, for this event type."""
    errors, warnings = [], []
    today_str = datetime.utcnow().strftime("%Y-%m-%d")

    insdc_entry = schema.get(INSDC_ACCESSION_FIELD_ID)
    if insdc_entry is None:
        errors.append(f"{line_label}: schema field '{INSDC_ACCESSION_FIELD_ID}' (0aa) not found in schema")
        return {"event": "error", "errors": errors, "warnings": warnings,
                "new_meta_doc": None, "genome_doc": None}
    acc_index = insdc_entry["index"]

    doc = map_normalized_fields(normalized, schema, line_label, errors, warnings)
    always_set = compute_always_set_fields(normalized, schema, vmr_version, errors, line_label)
    doc.update(always_set)

    doc[VJID_FIELD] = vjid
    doc[VERSION_FIELD] = 1
    doc[SEQ_VERSION_FIELD] = 1
    doc[ACTIVE_FLAG_FIELD] = True
    doc[RENEWED_DATE_FIELD] = today_str
    doc[INSERTION_DATE_FIELD] = today_str  # immutable, set once, never touched again
    doc[acc_index] = [resolved_version]

    dk_entry = schema.get("sequence_source")
    if dk_entry is None:
        errors.append(f"{line_label}: schema field 'sequence_source' (0dk) not found")
    elif dk_entry["index"] not in doc:
        doc[dk_entry["index"]] = [GENBANK_SOURCE_VALUE]

    genome_doc = None
    if sequence:
        genome_doc = {
            "0aa": resolved_version,
            "0dn": vjid,
            SEQ_VERSION_FIELD: 1,
            "0hh": sequence,
            "0hi": sha256_hex(sequence),
        }
    else:
        warnings.append(
            f"{line_label}: brand-new insertion but no sequence body found in the "
            f"FASTA record -- genome_sequences was NOT written for this accession."
        )

    return {
        "event": "brand-new-insertion" if not errors else "error",
        "errors": errors, "warnings": warnings,
        "new_meta_doc": doc if not errors else None,
        "genome_doc": genome_doc if not errors else None,
    }


def build_record_actions(vjid, es_id, old_doc, normalized, resolved_version,
                          sequence, schema, config, vmr_version, line_label):
    """UNCHANGED from l.py: version-diff / content-update / fill-only
    insertion / no-change for a record that ALREADY has an active document
    (found via accession lookup upstream; vjid is read off that document)."""
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

    always_set = compute_always_set_fields(normalized, schema, vmr_version, errors, line_label)

    result = {"event": None, "errors": errors, "warnings": warnings,
              "old_doc_update": None, "new_meta_doc": None,
              "fill_update": None, "change_log_doc": None, "genome_doc": None}

    # ---------------------------------------------------------------
    # 1. VERSION-DIFF: same base accession, different version suffix
    #    (also covers "stored accession has no version at all", since
    #    strip_version(old_acc) == old_acc in that case, and old_acc will
    #    still differ from the fully-versioned resolved_version).
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
            mapped = map_normalized_fields(normalized, schema, line_label, errors, warnings)
            new_doc.update(mapped)
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
        elif not values_equal(schema_field_id, old_val, coerced):
            updates_dict[idx] = coerced
            changes_log.append({
                "field": schema_field_id,
                "old_value": json.dumps(old_val),
                "new_value": json.dumps(coerced),
            })

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
        if INSERTION_DATE_FIELD not in new_doc:
            new_doc[INSERTION_DATE_FIELD] = today_str
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
        if INSERTION_DATE_FIELD not in old_doc:
            fill_dict[INSERTION_DATE_FIELD] = today_str
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

    logger = logging.getLogger("upsert_vjdb")
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
    logger.info("%d/%d record(s) include a sequence body", n_with_seq, len(fasta_records))

    if args.vmr_version is None:
        args.vmr_version = vmr_version_from_filename(args.fasta)
    logger.info("VMR version: %s", args.vmr_version)

    schema = load_schema(schema_path=args.schema, schema_url=args.schema_url)
    cache_dir = Path(args.cache_dir)
    use_cache = not args.no_cache

    # --- NEW: VMR/ICTV metadata JSON enrichment (optional) -----------------
    vmr_index = {}
    if args.vmr_json:
        vmr_records = vmr_ictv_merge.load_vmr_json(args.vmr_json)
        vmr_index, vmr_stats = vmr_ictv_merge.build_vmr_accession_index(vmr_records)
        logger.info(
            "VMR/ICTV JSON: %d record(s) loaded from %s (%d with a parseable "
            "accession field, %d without; %d distinct base accession(s) indexed)",
            vmr_stats["total_records"], args.vmr_json,
            vmr_stats["records_with_accession_field"],
            vmr_stats["records_without_accession_field"],
            vmr_stats["distinct_base_accessions"],
        )
        unrecognized, n_unrecognized = vmr_ictv_merge.report_unrecognized_keys(vmr_records)
        if n_unrecognized:
            logger.warning(
                "VMR/ICTV JSON: %d raw JSON key(s) across the file never matched "
                "any mapped field/alias/accession/lineage label -- sample: %s "
                "(if one of these is actually a field you expect to be mapped, "
                "its spelling differs from what vmr_ictv_merge.py recognizes; "
                "add it as an alias there).",
                n_unrecognized, unrecognized,
            )

    verify_ssl = not args.no_verify_ssl
    session = build_session(args.token)
    check_auth(session, args.es_url, args.index, verify_ssl)

    # --- vjid counter setup (for brand-new-insertion events) ---------------
    if args.offline_latest_accession:
        latest_known = args.offline_latest_accession
        logger.info("Using offline override for latest vjid: %s", latest_known)
    else:
        latest_known = fetch_latest_accession(session, args.es_url, args.index, args.genome_index, verify_ssl)
    next_vjid = make_accession_counter(latest_known)

    # --- Step 1: resolve accessions via NCBI --------------------------------
    resolved = {}  # accession_no_version -> (version, report, fallback_normalized)
    fasta_by_accession = {}
    for rec in fasta_records:
        fasta_by_accession[strip_version(rec["accession"])] = rec

    resolve_errors = []
    no_report_count = 0
    no_fallback_count = 0
    cache_hits = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = [
            executor.submit(resolve_one, rec["accession"], cache_dir, use_cache,
                             args.retries, args.retry_backoff)
            for rec in fasta_records
        ]
        progress = tqdm(concurrent.futures.as_completed(futures), total=len(futures),
                         unit="acc", desc="Resolving accessions", file=sys.__stderr__)
        for fut in progress:
            accession, version, report, fallback_normalized, was_cached, err = fut.result()
            acc_no_version = strip_version(accession)
            if err:
                resolve_errors.append(f"{accession}: {err}")
            else:
                resolved[acc_no_version] = (version, report, fallback_normalized)
                if report is None:
                    no_report_count += 1
                    if not fallback_normalized:
                        no_fallback_count += 1
                if was_cached:
                    cache_hits += 1
            progress.set_postfix(cached=cache_hits, failed=len(resolve_errors), no_report=no_report_count)

    if resolve_errors:
        logger.warning("%d accession(s) could not be resolved:", len(resolve_errors))
        for e in resolve_errors[:20]:
            logger.warning("  %s", e)
    if no_report_count:
        logger.info(
            "%d accession(s) had no NCBI Datasets virus report; %d of those also had "
            "no GenBank flat-file fallback data (will be skipped as errors).",
            no_report_count, no_fallback_count,
        )

    # --- Step 2: look up ACTIVE existing docs -- BY ACCESSION, not vjid ----
    all_accessions_no_version = list(resolved.keys())
    lookup = find_active_by_accession(
        session, args.es_url, args.index, schema, all_accessions_no_version,
        verify_ssl, batch_size=args.batch_size,
    )
    brand_new_accessions = [a for a, (eid, _) in lookup.items() if eid is None]
    existing_accessions = [a for a, (eid, _) in lookup.items() if eid is not None]
    logger.info(
        "%d accession(s) already exist in the index (active doc found), "
        "%d are brand new (no active doc found -- will mint a new vjid)",
        len(existing_accessions), len(brand_new_accessions),
    )

    # --- Step 3: build per-record actions -----------------------------------
    # brand-new-insertion and version-diff are the only two event types that
    # ever produce a genome_doc. For those, their metadata/change_log
    # actions are NOT flattened into the global bulk lists yet -- they're
    # held in `pending_genome_gated` and only released after we've verified
    # (via the genome bulk response's PER-ITEM status) that THIS record's
    # genome write actually succeeded.
    all_errors, all_warnings = [], []
    non_gated_meta_actions = []        # flat [action, doc, action, doc, ...] -- always written
    non_gated_change_log_actions = []  # same shape
    pending_genome_gated = []          # [{"label", "genome_pair", "meta_pairs", "change_log_pairs"}, ...]
    event_counts = {
        "brand-new-insertion": 0, "insertion": 0, "content-update": 0,
        "version-diff": 0, "no-change": 0, "error": 0,
    }
    vmr_matched_count = 0

    for acc_no_version, (version, report, fallback_normalized) in resolved.items():
        rec = fasta_by_accession.get(acc_no_version, {})
        report_normalized = normalize_report(report) if report is not None else {}
        header_normalized = header_extras_as_normalized_fields(rec)
        normalized = merge_normalized(report_normalized, header_normalized, fallback_normalized)
        sequence = rec.get("sequence") or ""
        line_label = f"accession {acc_no_version}"

        es_id, old_doc = lookup.get(acc_no_version, (None, None))

        # --- NEW: merge matched VMR/ICTV JSON fields into `normalized`,
        # BEFORE building either event, so build_brand_new_insertion() /
        # build_record_actions() see them exactly like any other normalized
        # field (no changes needed to either function).
        vmr_matches = vmr_index.get(acc_no_version, [])
        if vmr_matches:
            vmr_ictv_merge.merge_vmr_into_normalized(
                normalized, vmr_matches, old_doc, schema, line_label, all_warnings
            )
            vmr_matched_count += 1

        if es_id is None:
            # --- BRAND-NEW INSERTION ---
            if report is None and not fallback_normalized:
                all_errors.append(
                    f"accession {acc_no_version}: no NCBI Datasets VIRUS report AND no "
                    f"GenBank flat-file fallback data available -- refusing to insert "
                    f"a metadata-only stub. Skipping."
                )
                event_counts["error"] += 1
                continue
            if report is None and fallback_normalized:
                all_warnings.append(
                    f"accession {acc_no_version}: no NCBI Datasets VIRUS report -- used "
                    f"the GenBank flat-file fallback (organism/length/country only, no "
                    f"taxonomy/host/isolate data). VERIFY this accession belongs in the "
                    f"database before trusting this record."
                )

            vjid = next_vjid()
            line_label = f"accession {acc_no_version} (new vjid {vjid})"
            r = build_brand_new_insertion(
                vjid, normalized, version, sequence, schema, args.vmr_version, line_label
            )
            all_errors.extend(r["errors"])
            all_warnings.extend(r["warnings"])
            event_counts[r["event"]] = event_counts.get(r["event"], 0) + 1

            if r["genome_doc"] is not None:
                pending_genome_gated.append({
                    "label": line_label,
                    "genome_pair": ({"index": {}}, r["genome_doc"]),
                    "meta_pairs": [({"index": {}}, r["new_meta_doc"])] if r["new_meta_doc"] is not None else [],
                    "change_log_pairs": [],  # brand-new NEVER writes change_log, per spec
                })
            elif r["new_meta_doc"] is not None:
                # No sequence body was available at all (already warned inside
                # build_brand_new_insertion) -- nothing to gate against, so the
                # metadata doc is written unconditionally, same as before.
                non_gated_meta_actions.append({"index": {}})
                non_gated_meta_actions.append(r["new_meta_doc"])
            continue

        # --- EXISTING RECORD: same rules as l.py, vjid taken from old_doc ---
        vjid = old_doc.get(VJID_FIELD)
        if not vjid:
            all_errors.append(f"{line_label}: matched active document {es_id} has no vjid ({VJID_FIELD}) set")
            event_counts["error"] += 1
            continue
        line_label = f"vjid {vjid} ({acc_no_version})"

        r = build_record_actions(
            vjid, es_id, old_doc, normalized, version, sequence,
            schema, config, args.vmr_version, line_label,
        )
        all_errors.extend(r["errors"])
        all_warnings.extend(r["warnings"])
        event_counts[r["event"]] = event_counts.get(r["event"], 0) + 1

        if r["event"] == "version-diff" and r["genome_doc"] is not None:
            # Gate ALL THREE actions together: if the genome write for the
            # new version fails, do NOT deactivate the old (still-working)
            # record, and do NOT create a new "active" record whose sequence
            # never actually landed.
            meta_pairs = []
            if r["old_doc_update"]:
                meta_pairs.append(({"update": {"_id": r["old_doc_update"]["_id"]}}, {"doc": r["old_doc_update"]["doc"]}))
            if r["new_meta_doc"] is not None:
                meta_pairs.append(({"index": {}}, r["new_meta_doc"]))
            change_log_pairs = [({"index": {}}, r["change_log_doc"])] if r["change_log_doc"] is not None else []
            pending_genome_gated.append({
                "label": line_label,
                "genome_pair": ({"index": {}}, r["genome_doc"]),
                "meta_pairs": meta_pairs,
                "change_log_pairs": change_log_pairs,
            })
            continue

        # content-update / fill-only insertion / no-change / version-diff
        # WITHOUT a sequence body (already warned) -- nothing to gate on.
        if r["old_doc_update"]:
            non_gated_meta_actions.append({"update": {"_id": r["old_doc_update"]["_id"]}})
            non_gated_meta_actions.append({"doc": r["old_doc_update"]["doc"]})
        if r["new_meta_doc"] is not None:
            non_gated_meta_actions.append({"index": {}})
            non_gated_meta_actions.append(r["new_meta_doc"])
        if r["fill_update"]:
            non_gated_meta_actions.append({"update": {"_id": r["fill_update"]["_id"]}})
            non_gated_meta_actions.append({"doc": r["fill_update"]["doc"]})
        if r["change_log_doc"] is not None:
            non_gated_change_log_actions.append({"index": {}})
            non_gated_change_log_actions.append(r["change_log_doc"])

    # --- Preview totals for reporting (genome-gated records not yet flattened) ---
    genome_gated_meta_count = sum(len(g["meta_pairs"]) for g in pending_genome_gated)
    genome_gated_cl_count = sum(len(g["change_log_pairs"]) for g in pending_genome_gated)

    # --- Pre-flight payload-size check (this is what makes 413s visible in
    # dry run, instead of only surfacing at commit time). We're not sending
    # anything -- just running the same chunk-size math execute_bulk() would,
    # to catch any SINGLE document big enough to 413 even after chunking.
    #
    # For "genome" specifically, a document that alone exceeds the cap is
    # now EXCLUDED from the run (not just warned about) -- a >max_chunk_bytes
    # single genome doc has nowhere smaller to be split into, so leaving it
    # in would just reproduce the 413 crash at --commit time. Its vjid stays
    # reserved/safe to retry once the underlying data issue is fixed. ---
    def _preview_oversized(label, pairs_flat, max_chunk_bytes, exclude_ids=None):
        found = []
        for i in range(0, len(pairs_flat), 2):
            pair = pairs_flat[i:i + 2]
            size = sum(len(json.dumps(a).encode("utf-8")) + 1 for a in pair)
            if size > max_chunk_bytes:
                doc = pair[1] if len(pair) > 1 else pair[0]
                if label == "genome":
                    # Genome docs always carry BOTH "0aa" (accession) and "0dn"
                    # (vjid) -- show accession first since that's what you'd
                    # actually look up in NCBI, plus the raw (uncompressed,
                    # un-JSON-escaped) sequence length so you can immediately
                    # tell whether this is a parsing artifact (multiple records
                    # collapsed into one) vs a genuinely huge sequence.
                    accession = doc.get("0aa", "?")
                    vjid = doc.get("0dn", "?")
                    raw_seq_len = len(doc.get("0hh") or "")
                    ident = f"accession={accession} vjid={vjid} raw_sequence_length={raw_seq_len}bp"
                else:
                    ident = doc.get(VJID_FIELD) or doc.get("0aa") or doc.get("0dn") or "?"
                found.append(
                    f"{label}: document ({ident}) serializes to {size} byte(s), which "
                    f"ALONE exceeds the {max_chunk_bytes}-byte chunk cap"
                    + (" -- EXCLUDED from this run for manual review."
                       if exclude_ids is not None else
                       " -- this will still be sent as an oversized request and "
                       "may 413 regardless of chunking. Consider trimming this "
                       "record's fields.")
                )
                if exclude_ids is not None:
                    exclude_ids.add(id(doc))
        return found

    oversized_genome_doc_ids = set()
    genome_preview_flat = []
    for g in pending_genome_gated:
        genome_preview_flat.extend(g["genome_pair"])
    genome_size_warnings = _preview_oversized(
        "genome", genome_preview_flat, args.max_chunk_bytes,
        exclude_ids=oversized_genome_doc_ids,
    )
    size_warnings = (
        _preview_oversized("metadata", non_gated_meta_actions, args.max_chunk_bytes)
        + _preview_oversized("change_log", non_gated_change_log_actions, args.max_chunk_bytes)
        + genome_size_warnings
    )
    all_warnings.extend(size_warnings)

    n_oversized_skipped = 0
    if oversized_genome_doc_ids:
        n_before = len(pending_genome_gated)
        pending_genome_gated = [
            g for g in pending_genome_gated
            if id(g["genome_pair"][1]) not in oversized_genome_doc_ids
        ]
        n_oversized_skipped = n_before - len(pending_genome_gated)
        event_counts["skipped-oversized"] = event_counts.get("skipped-oversized", 0) + n_oversized_skipped
        # SOFT skip, not a hard error: this does NOT go into all_errors and does
        # NOT abort the dry run / block --commit. The rest of the batch proceeds
        # normally; only these specific records (and their vjids) are held back
        # for a future run once the underlying data/limit issue is resolved.
        all_warnings.append(
            f"{n_oversized_skipped} genome document(s) exceed --max-chunk-bytes "
            f"({args.max_chunk_bytes} bytes) on their own and were SKIPPED for this "
            f"run -- see warnings above for identifiers. Their metadata/change_log "
            f"writes were also skipped (nothing to gate them on); their vjids stay "
            f"reserved and safe to retry later. Fix/exclude the underlying record(s), "
            f"or raise --max-chunk-bytes together with your nginx client_max_body_size "
            f"/ ES http.max_content_length if this data is legitimate, then re-run."
        )
        genome_gated_meta_count = sum(len(g["meta_pairs"]) for g in pending_genome_gated)
        genome_gated_cl_count = sum(len(g["change_log_pairs"]) for g in pending_genome_gated)

    total_meta_preview = len(non_gated_meta_actions) // 2 + genome_gated_meta_count
    total_cl_preview = len(non_gated_change_log_actions) // 2 + genome_gated_cl_count
    total_genome_preview = len(pending_genome_gated)

    # --- Report --------------------------------------------------------------
    print(f"\nParsed:               {len(fasta_records)}")
    print(f"Resolved via NCBI:    {len(resolved)}  (failed: {len(resolve_errors)}, no-report: {no_report_count})")
    print(f"Brand new accessions: {len(brand_new_accessions)}")
    print(f"Existing accessions:  {len(existing_accessions)}")
    if args.vmr_json:
        print(f"VMR/ICTV JSON matches: {vmr_matched_count}/{len(resolved)} accession(s) "
              f"matched a VMR record (see vmr_ictv_merge.py assumptions if this looks off)")
    print(f"Events: brand-new-insertion={event_counts.get('brand-new-insertion', 0)}  "
          f"insertion(fill)={event_counts.get('insertion', 0)}  "
          f"content-update={event_counts.get('content-update', 0)}  "
          f"version-diff={event_counts.get('version-diff', 0)}  "
          f"no-change={event_counts.get('no-change', 0)}")
    print(f"Errors:               {len(all_errors)}")
    print(f"Warnings:             {len(all_warnings)}")
    if n_oversized_skipped:
        print(f"Skipped (oversized):  {n_oversized_skipped}  (genome doc(s) over "
              f"--max-chunk-bytes={args.max_chunk_bytes}; vjids reserved, safe to retry)")
        print("--- Skipped record(s), full list (accession / vjid / raw sequence length) ---")
        for w in genome_size_warnings:
            print("  ", w)

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
              f"{total_meta_preview} metadata action(s) "
              f"({genome_gated_meta_count} gated on their genome write succeeding), "
              f"{total_cl_preview} change_log doc(s), "
              f"{total_genome_preview} genome_sequences doc(s).")
        print("(Note: vjids for brand-new accessions above are provisional -- "
              "if you don't --commit now, re-running later may mint different "
              "vjids for the same accessions, since the counter is re-derived "
              "from live DB state each run.)")
        return

    meta_bulk = list(non_gated_meta_actions)
    change_log_bulk = list(non_gated_change_log_actions)

    # --- Genome writes FIRST, gated: only release a record's metadata/
    # change_log actions into the bulk lists once we've confirmed (via the
    # bulk response's PER-ITEM status) that ITS SPECIFIC genome write
    # succeeded. ---
    if pending_genome_gated:
        genome_bulk = []
        for g in pending_genome_gated:
            genome_bulk.extend(g["genome_pair"])

        print(f"\n=== Writing {len(genome_bulk)//2} gated genome doc(s) "
              f"(brand-new + version-diff) to {args.genome_index} ===")
        try:
            set_index_write_block(session, args.es_url, args.genome_index, False, verify_ssl)
            resp = execute_bulk(session, args.es_url, args.genome_index, genome_bulk, verify_ssl,
                                 max_chunk_bytes=args.max_chunk_bytes)
        finally:
            set_index_write_block(session, args.es_url, args.genome_index, True, verify_ssl)

        items = resp.get("items", [])
        if len(items) != len(pending_genome_gated):
            # Response shape didn't line up 1:1 with what we sent -- we can no
            # longer trust positional matching, so refuse to guess: skip ALL
            # gated metadata/change_log writes this run rather than risk
            # attaching the wrong result to the wrong record.
            logger.error(
                "genome bulk response item count (%d) != expected record count (%d) "
                "-- cannot safely match results to records. Skipping ALL gated "
                "metadata/change_log writes this run; their vjids remain reserved "
                "and safe to retry on a future run.",
                len(items), len(pending_genome_gated),
            )
            gated_ok = [False] * len(pending_genome_gated)
        else:
            gated_ok = []
            for g, item in zip(pending_genome_gated, items):
                op_result = next(iter(item.values())) if item else {"error": "empty item"}
                ok = not op_result.get("error")
                if not ok:
                    logger.error(
                        "%s: genome write FAILED (%s) -- its metadata/change_log "
                        "will NOT be written this run; vjid stays reserved and "
                        "safe to retry later.", g["label"], op_result.get("error"),
                    )
                gated_ok.append(ok)

        n_released, n_blocked = 0, 0
        for g, ok in zip(pending_genome_gated, gated_ok):
            if ok:
                for action, doc in g["meta_pairs"]:
                    meta_bulk.append(action)
                    meta_bulk.append(doc)
                for action, doc in g["change_log_pairs"]:
                    change_log_bulk.append(action)
                    change_log_bulk.append(doc)
                n_released += 1
            else:
                n_blocked += 1
        print(f"Genome-gated records: {n_released} released to metadata/change_log, "
              f"{n_blocked} blocked (genome write failed).")

    if meta_bulk:
        print(f"\n=== Writing {len(meta_bulk)//2} metadata action(s) to {args.index} ===")
        resp = execute_bulk(session, args.es_url, args.index, meta_bulk, verify_ssl,
                             max_chunk_bytes=args.max_chunk_bytes)
        check_bulk_errors(resp, logger, "metadata bulk")

    if change_log_bulk:
        print(f"=== Writing {len(change_log_bulk)//2} change_log doc(s) to {args.change_log_index} ===")
        resp = execute_bulk(session, args.es_url, args.change_log_index, change_log_bulk, verify_ssl,
                             max_chunk_bytes=args.max_chunk_bytes)
        check_bulk_errors(resp, logger, "change_log bulk")

    print("\nDone.")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True,
                         help="Multi-FASTA file, headers '>ACCESSION | organism=... | length=... | "
                              "description=...' (legacy 'vjid | accession | ...' also parsed but "
                              "vjid is ignored for lookup) + sequence bodies")
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
    parser.add_argument("--batch-size", type=int, default=100,
                         help="accessions per lookup query batch (default: 100; each accession "
                              "contributes 3 'should' clauses, so keep this reasonable)")
    parser.add_argument("--max-chunk-bytes", type=int, default=900_000,
                         help="Max serialized _bulk request size in bytes per chunk "
                              "(default: 900000, safely under nginx's common 1MB "
                              "client_max_body_size default). A single document larger "
                              "than this on its own cannot be sent regardless of chunking; "
                              "for genome_sequences such a document is excluded from the "
                              "run at the dry-run stage instead of crashing at --commit. "
                              "Raise this only after confirming your proxy/ES limits "
                              "(client_max_body_size / http.max_content_length) allow it.")
    parser.add_argument("--cache-dir", default=".cache/reports")
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--ncbi-api-key", default=None)
    parser.add_argument("--eutils-rate", type=float, default=None)
    parser.add_argument("--datasets-rate", type=float, default=None)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--retry-backoff", type=float, default=1.0)
    parser.add_argument("--vmr-version", default=None)
    parser.add_argument("--vmr-json", default=None,
                         help="Path to a VMR/ICTV metadata JSON file (array of record "
                              "objects with a 'Virus GENBANK accession'-style field). "
                              "When given, records are matched to insertion accessions "
                              "(with or without version) and their mapped ICTV taxonomy/ "
                              "isolate fields are merged in. See vmr_ictv_merge.py.")
    parser.add_argument("--offline-latest-accession", default=None,
                         help="skip the live ES lookup for the vjid counter seed and start from "
                              "this vjid instead (for offline testing)")
    parser.add_argument("--log-file", default=None)
    args = parser.parse_args()

    input_base = Path(args.fasta).stem
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = args.log_file or f"{input_base}_upsert_{timestamp}.log"
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
