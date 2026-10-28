#!/usr/bin/env python3
"""Replace changelog field aliases with their DB_Scheme.json field indices.

The script is read-only by default. Pass --commit to update matching records
in place; updates replace only the document's ``changes`` array.

Usage:
    python3 scripts/repair_change_log_field_ids.py
    python3 scripts/repair_change_log_field_ids.py --commit
    python3 scripts/repair_change_log_field_ids.py --schema DB_Scheme.json
"""

import argparse
import json
import re
from collections import Counter

import requests
import yaml


FIELD_INDEX_RE = re.compile(r"^0[a-z]{2}$")
DEFAULT_CONFIG = "config.yml"
DEFAULT_SCROLL_KEEPALIVE = "2m"
DEFAULT_BATCH_SIZE = 500
REQUEST_TIMEOUT = 60
BULK_TIMEOUT = 180


def load_config(path):
    with open(path, "r", encoding="utf-8") as config_file:
        config = yaml.safe_load(config_file)
    if not isinstance(config, dict):
        raise ValueError(f"Configuration file {path!r} must contain a YAML mapping.")
    return config


def load_schema(schema_path, schema_url):
    if schema_path:
        with open(schema_path, "r", encoding="utf-8") as schema_file:
            entries = json.load(schema_file)
    else:
        response = requests.get(schema_url, timeout=REQUEST_TIMEOUT)
        response.raise_for_status()
        entries = response.json()

    if not isinstance(entries, list):
        raise ValueError("DB_Scheme.json must contain a JSON array of field definitions.")

    candidates = {}
    for position, entry in enumerate(entries):
        if not isinstance(entry, dict):
            raise ValueError(f"Schema entry {position} must be a JSON object.")
        field_id = entry.get("vjdbv1.1_field_id")
        field_index = entry.get("vjdbv1.1_field_index")
        if not isinstance(field_id, str) or not isinstance(field_index, str):
            raise ValueError(
                f"Schema entry {position} must define string "
                "'vjdbv1.1_field_id' and 'vjdbv1.1_field_index' values."
            )

        for alias in (field_id, entry.get("vjdbv1.1_name")):
            if isinstance(alias, str) and alias:
                candidates.setdefault(alias, set()).add(field_index)

    aliases = {}
    ambiguous = set()
    for alias, indices in candidates.items():
        if len(indices) == 1:
            aliases[alias] = next(iter(indices))
        else:
            ambiguous.add(alias)
    return aliases, ambiguous


def field_index_for(field, aliases, ambiguous):
    if not isinstance(field, str) or FIELD_INDEX_RE.fullmatch(field):
        return None, None
    if field in ambiguous:
        return None, "ambiguous"
    if field not in aliases:
        return None, "unmapped"
    return aliases[field], None


def rewrite_changes(changes, aliases, ambiguous):
    """Return a corrected copy and counts of resolved/unresolved field names."""
    if not isinstance(changes, list):
        return changes, Counter(), Counter()

    corrected = []
    replacements = Counter()
    unresolved = Counter()
    for change in changes:
        if not isinstance(change, dict):
            corrected.append(change)
            continue
        copied_change = dict(change)
        if "field" in copied_change:
            new_field, issue = field_index_for(copied_change["field"], aliases, ambiguous)
            if new_field is not None:
                replacements[copied_change["field"]] += 1
                copied_change["field"] = new_field
            elif issue:
                unresolved[(copied_change["field"], issue)] += 1
        corrected.append(copied_change)
    return corrected, replacements, unresolved


def make_session(config):
    try:
        token = config["elasticsearch"]["auth_token"]
    except (KeyError, TypeError) as error:
        raise ValueError("Configuration must define elasticsearch.auth_token.") from error
    session = requests.Session()
    session.headers.update({"Authorization": f"ApiKey {token}"})
    return session


def send_bulk_updates(session, es_url, index_name, updates, verify_ssl):
    if not updates:
        return 0

    lines = []
    for document_id, changes in updates:
        lines.append(json.dumps({"update": {"_id": document_id}}))
        lines.append(json.dumps({"doc": {"changes": changes}}))
    payload = "\n".join(lines) + "\n"
    response = session.post(
        f"{es_url.rstrip('/')}/{index_name}/_bulk",
        data=payload.encode("utf-8"),
        headers={"Content-Type": "application/x-ndjson"},
        verify=verify_ssl,
        timeout=BULK_TIMEOUT,
    )
    response.raise_for_status()
    result = response.json()
    items = result.get("items", [])
    failures = []
    for item in items:
        update_result = item.get("update")
        if (
            update_result is None
            or update_result.get("error")
            or not 200 <= update_result.get("status", 0) < 300
        ):
            failures.append(update_result or item)
    if result.get("errors") or failures:
        details = [
            f"{item.get('_id')}: {item.get('error')}"
            for item in failures[:10]
        ]
        raise RuntimeError(
            f"Elasticsearch bulk update reported {len(failures)} failed item(s)"
            + (f": {'; '.join(details)}" if details else ".")
        )
    if len(items) != len(updates):
        raise RuntimeError(
            f"Elasticsearch bulk response returned {len(items)} item(s) for "
            f"{len(updates)} requested update(s)."
        )
    return len(items)


def repair_change_log(session, es_url, index_name, aliases, ambiguous,
                      verify_ssl, batch_size, scroll_keepalive, commit):
    search_url = f"{es_url.rstrip('/')}/{index_name}/_search"
    response = session.post(
        search_url,
        params={"scroll": scroll_keepalive},
        json={
            "size": batch_size,
            "sort": ["_doc"],
            "_source": ["changes"],
            "query": {"match_all": {}},
        },
        verify=verify_ssl,
        timeout=REQUEST_TIMEOUT,
    )
    response.raise_for_status()
    result = response.json()
    scroll_id = result.get("_scroll_id")
    if result.get("hits", {}).get("hits") and not scroll_id:
        raise RuntimeError(
            "Elasticsearch returned records without a scroll ID; refusing to "
            "continue with an incomplete changelog scan."
        )

    scanned = 0
    changed_documents = 0
    updated_documents = 0
    replacement_counts = Counter()
    unresolved_counts = Counter()
    samples = []

    try:
        hits = result.get("hits", {}).get("hits", [])
        while hits:
            scanned += len(hits)
            batch_updates = []
            for hit in hits:
                source = hit.get("_source", {})
                corrected, replacements, unresolved = rewrite_changes(
                    source.get("changes"), aliases, ambiguous
                )
                replacement_counts.update(replacements)
                unresolved_counts.update(unresolved)
                if replacements:
                    changed_documents += 1
                    batch_updates.append((hit["_id"], corrected))
                    if len(samples) < 20:
                        samples.append(
                            {
                                "_id": hit["_id"],
                                "fields": dict(replacements),
                            }
                        )

            if commit:
                updated_documents += send_bulk_updates(
                    session, es_url, index_name, batch_updates, verify_ssl
                )

            if not scroll_id:
                break
            response = session.post(
                f"{es_url.rstrip('/')}/_search/scroll",
                json={"scroll": scroll_keepalive, "scroll_id": scroll_id},
                verify=verify_ssl,
                timeout=REQUEST_TIMEOUT,
            )
            response.raise_for_status()
            result = response.json()
            scroll_id = result.get("_scroll_id", scroll_id)
            hits = result.get("hits", {}).get("hits", [])
    finally:
        if scroll_id:
            response = session.delete(
                f"{es_url.rstrip('/')}/_search/scroll",
                json={"scroll_id": [scroll_id]},
                verify=verify_ssl,
                timeout=REQUEST_TIMEOUT,
            )
            response.raise_for_status()

    return {
        "scanned": scanned,
        "changed_documents": changed_documents,
        "updated_documents": updated_documents,
        "replacement_counts": replacement_counts,
        "unresolved_counts": unresolved_counts,
        "samples": samples,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=DEFAULT_CONFIG, help="YAML settings file (default: config.yml)")
    parser.add_argument("--schema", help="Local DB_Scheme.json; otherwise fetch the configured schema URL")
    parser.add_argument("--schema-url", help="Override the schema URL from the config file")
    parser.add_argument("--es-url", help="Override elasticsearch.es_url from the config file")
    parser.add_argument("--change-log-index", help="Override elasticsearch.change_log_index from the config file")
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    parser.add_argument("--scroll-keepalive", default=DEFAULT_SCROLL_KEEPALIVE)
    parser.add_argument("--no-verify-ssl", action="store_true")
    parser.add_argument(
        "--commit",
        action="store_true",
        help="Apply updates; without this option the script only reports what would change.",
    )
    args = parser.parse_args()

    if args.batch_size < 1:
        parser.error("--batch-size must be at least 1")

    try:
        config = load_config(args.config)
        elasticsearch = config["elasticsearch"]
        schema_config = config["schema"]
        es_url = args.es_url or elasticsearch["es_url"]
        index_name = args.change_log_index or elasticsearch["change_log_index"]
        schema_url = args.schema_url or schema_config["schema_url"]
        aliases, ambiguous = load_schema(args.schema, schema_url)
        session = make_session(config)
        result = repair_change_log(
            session=session,
            es_url=es_url,
            index_name=index_name,
            aliases=aliases,
            ambiguous=ambiguous,
            verify_ssl=not args.no_verify_ssl,
            batch_size=args.batch_size,
            scroll_keepalive=args.scroll_keepalive,
            commit=args.commit,
        )
    except (KeyError, TypeError) as error:
        raise SystemExit(f"Invalid configuration: {error}") from error
    except (OSError, ValueError, requests.exceptions.RequestException, RuntimeError) as error:
        raise SystemExit(f"Changelog field repair failed: {error}") from error

    print(f"Scanned {result['scanned']} change_log record(s).")
    print(f"Records with mapped field aliases: {result['changed_documents']}.")
    print("Field replacements:")
    if result["replacement_counts"]:
        for alias, count in sorted(result["replacement_counts"].items()):
            print(f"  {alias} -> {aliases[alias]}: {count}")
    else:
        print("  none")

    if result["unresolved_counts"]:
        print("Unresolved field values (left unchanged):")
        for (field, issue), count in sorted(result["unresolved_counts"].items()):
            print(f"  {field!r} ({issue}): {count}")

    if result["samples"]:
        print("Sample affected records:")
        for sample in result["samples"]:
            print(f"  {sample['_id']}: {sample['fields']}")

    if args.commit:
        print(f"Updated {result['updated_documents']} change_log record(s) in place.")
    else:
        print("Dry run only; no Elasticsearch records were changed. Use --commit to apply.")


if __name__ == "__main__":
    main()
