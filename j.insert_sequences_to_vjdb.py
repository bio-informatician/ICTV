#!/usr/bin/env python3
"""
insert_genome_sequences.py

Reads a CSV file of (vjid, accession) pairs and a FASTA file whose headers look like:

    >AB000906 | organism=Infectious flacherie virus | length=9650 | description=...

For each CSV row, finds the FASTA record whose accession (version-stripped) matches
the CSV accession (also version-stripped), computes a SHA256 hash of the sequence,
and bulk-inserts the result into an Elasticsearch index using the field aliases:

    0aa -> insdc_accession_id  (full accession WITH version, as given in the CSV)
    0dn -> virjen_id           (vjid, as given in the CSV)
    0hh -> genome_sequence     (the raw sequence, concatenated, no newlines)
    0hi -> genome_hash         (sha256 hex digest of the sequence)

Usage:
    python3 insert_genome_sequences.py \
        --csv mapping.csv \
        --fasta VMR_MSL41.v1.20260320.fasta \
        [--es-url https://virjendb.org:9201] \
        [--index genome_sequences] \
        [--auth-token OERyN3JKOEJPaThtV3FrMkUtWVE6Q3FqSVFidU02NF9tT29VWXloQ1RjQQ==] \
        [--batch-size 500] \
        [--dry-run] \
        [--insecure] \
        [--report-missing missing.csv]

Notes:
- --auth-token is expected to already be a base64-encoded "user:pass" string
  (i.e. exactly what goes after "Basic " in the Authorization header).
- Accession matching strips any ".<version>" suffix from BOTH the CSV accession
  and the FASTA header accession before comparing, but the FULL, versioned
  accession from the CSV is what gets stored in insdc_accession_id (0aa).
- Only standard library is used (csv, hashlib, json, argparse, urllib), so no
  extra pip installs are required.
"""
import gzip
import argparse
import csv
import hashlib
import json
import ssl
import sys
import urllib.error
import urllib.request

ES_URL_DEFAULT = "https://virjendb.org:9201"
INDEX_NAME_DEFAULT = "genome_sequences"
AUTH_TOKEN_DEFAULT = "OERyN3JKOEJPaThtV3FrMkUtWVE6Q3FqSVFidU02NF9tT29VWXloQ1RjQQ=="


def strip_version(accession: str) -> str:
    """Remove a trailing '.<version>' suffix, e.g. 'AB000906.1' -> 'AB000906'."""
    accession = accession.strip()
    if "." in accession:
        return accession.rsplit(".", 1)[0]
    return accession


def parse_fasta(fasta_path: str) -> dict:
    """
    Parse a FASTA file into {accession_without_version: sequence}.
    Header format assumed: '>ACCESSION | organism=... | length=... | description=...'
    Only the first whitespace-delimited token after '>' is used as the accession.
    """
    seqs = {}
    current_acc = None
    current_chunks = []

    def flush():
        if current_acc is not None:
            seq = "".join(current_chunks).replace(" ", "").upper()
            seqs[current_acc] = seq

    with open(fasta_path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                # first token before whitespace or '|' is the accession
                token = header.split("|", 1)[0].strip()
                token = token.split()[0] if token.split() else token
                current_acc = strip_version(token)
                current_chunks = []
            else:
                current_chunks.append(line.strip())
        flush()

    return seqs


def read_csv_rows(csv_path: str):
    with open(csv_path, "r", newline="") as fh:
        reader = csv.DictReader(fh)
        # normalize header names just in case of case/whitespace differences
        fieldmap = {k.strip().lower(): k for k in reader.fieldnames}
        if "vjid" not in fieldmap or "accession" not in fieldmap:
            raise ValueError(
                f"CSV must have 'vjid' and 'accession' columns, found: {reader.fieldnames}"
            )
        vjid_key = fieldmap["vjid"]
        acc_key = fieldmap["accession"]
        for row in reader:
            yield row[vjid_key].strip(), row[acc_key].strip()


def sha256_hex(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def update_index_settings(es_url, index_name, auth_token, settings, insecure=False):
    """
    Updates the settings of an Elasticsearch index via PUT /<index>/_settings.
    """
    url = f"{es_url.rstrip('/')}/{index_name}/_settings"
    payload = json.dumps(settings).encode("utf-8")

    req = urllib.request.Request(url, data=payload, method="PUT")
    req.add_header("Content-Type", "application/json")
    req.add_header("Authorization", f"ApiKey {auth_token}")

    ctx = None
    if insecure:
        ctx = ssl.create_default_context()
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE

    try:
        with urllib.request.urlopen(req, context=ctx, timeout=60) as resp:
            body = resp.read().decode("utf-8")
            return json.loads(body)
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"HTTP {e.code} updating settings for {index_name}: {body}") from e
    except urllib.error.URLError as e:
        raise RuntimeError(f"Could not reach Elasticsearch at {url}: {e}") from e


def bulk_insert(es_url, index_name, auth_token, actions, insecure=False):
    """
    actions: list of (vjid, accession_full, sequence, seq_hash) tuples
    Sends one _bulk request.
    """
    if not actions:
        return {"took": 0, "errors": False, "items": []}

    lines = []
    for vjid, accession_full, sequence, seq_hash in actions:
        lines.append(json.dumps({"index": {"_index": index_name}}))
        lines.append(
            json.dumps(
                {
                    "0aa": accession_full,
                    "0dn": vjid,
                    "0hh": sequence,
                    "0hi": seq_hash,
                }
            )
        )
    payload = ("\n".join(lines) + "\n").encode("utf-8")

    url = f"{es_url.rstrip('/')}/_bulk"

    def _post(data: bytes):
        compressed = gzip.compress(data)
        req = urllib.request.Request(url, data=compressed, method="POST")
        req.add_header("Content-Type", "application/x-ndjson")
        req.add_header("Content-Encoding", "gzip")
        req.add_header("Authorization", f"ApiKey {auth_token}")
        ctx = None
        if insecure:
            ctx = ssl.create_default_context()
            ctx.check_hostname = False
            ctx.verify_mode = ssl.CERT_NONE

        with urllib.request.urlopen(req, context=ctx, timeout=180) as resp:
            body = resp.read().decode("utf-8")
            return json.loads(body)

    try:
        return _post(payload)
    except urllib.error.HTTPError as e:
        if e.code == 413 and len(actions) > 1:
            e.read()  # drain response body
            mid = len(actions) // 2
            print(
                f"    413 Request Entity Too Large for {len(actions)} docs "
                f"({len(payload)} bytes) -> splitting into {mid} + {len(actions) - mid} and retrying",
                file=sys.stderr,
            )
            first = bulk_insert(es_url, index_name, auth_token, actions[:mid], insecure=insecure)
            second = bulk_insert(es_url, index_name, auth_token, actions[mid:], insecure=insecure)
            return {
                "errors": first.get("errors", False) or second.get("errors", False),
                "items": first.get("items", []) + second.get("items", []),
            }
        body = e.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"HTTP {e.code} from Elasticsearch: {body}") from e
    except urllib.error.URLError as e:
        raise RuntimeError(f"Could not reach Elasticsearch at {url}: {e}") from e


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", required=True, help="Path to CSV file with columns vjid,accession")
    ap.add_argument("--fasta", required=True, help="Path to FASTA file")
    ap.add_argument("--es-url", default=ES_URL_DEFAULT, help="Elasticsearch base URL")
    ap.add_argument("--index", default=INDEX_NAME_DEFAULT, help="Elasticsearch index name")
    ap.add_argument("--auth-token", default=AUTH_TOKEN_DEFAULT, help="Base64 'user:pass' token for Basic auth")
    ap.add_argument("--batch-size", type=int, default=50, help="Number of docs per _bulk request (auto-splits further on HTTP 413)")
    ap.add_argument("--dry-run", action="store_true", help="Parse and match only, do not call Elasticsearch")
    ap.add_argument("--insecure", action="store_true", help="Disable TLS certificate verification")
    ap.add_argument("--report-missing", default=None, help="Optional path to write unmatched rows to a CSV")
    args = ap.parse_args()

    print(f"Parsing FASTA file: {args.fasta}", file=sys.stderr)
    fasta_index = parse_fasta(args.fasta)
    print(f"  -> {len(fasta_index)} FASTA records indexed by accession (no version)", file=sys.stderr)

    matched = []
    missing = []

    for vjid, accession_full in read_csv_rows(args.csv):
        acc_no_version = strip_version(accession_full)
        seq = fasta_index.get(acc_no_version)
        if seq is None:
            missing.append((vjid, accession_full))
            continue
        seq_hash = sha256_hex(seq)
        matched.append((vjid, accession_full, seq, seq_hash))

    print(f"Matched: {len(matched)}   Missing: {len(missing)}", file=sys.stderr)

    if missing:
        print("Unmatched accessions (no FASTA record found):", file=sys.stderr)
        for vjid, acc in missing[:20]:
            print(f"  - {vjid}: {acc}", file=sys.stderr)
        if len(missing) > 20:
            print(f"  ... and {len(missing) - 20} more", file=sys.stderr)
        if args.report_missing:
            with open(args.report_missing, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["vjid", "accession"])
                w.writerows(missing)
            print(f"Wrote missing rows to {args.report_missing}", file=sys.stderr)

    if args.dry_run:
        print("\n--dry-run set: no data was sent to Elasticsearch.", file=sys.stderr)
        for vjid, accession_full, seq, seq_hash in matched[:5]:
            print(
                f"  SAMPLE -> vjid={vjid} accession={accession_full} "
                f"hash={seq_hash} seq_len={len(seq)} seq_preview={seq[:40]}...",
                file=sys.stderr,
            )
        return

    # Unlock the index before writing
    if matched:
        print(f"\nUnlocking index '{args.index}' (index.blocks.write=null)...", file=sys.stderr)
        update_index_settings(
            args.es_url, args.index, args.auth_token,
            {"index.blocks.write": None},
            insecure=args.insecure
        )

    total_indexed = 0
    total_errors = 0

    try:
        for i in range(0, len(matched), args.batch_size):
            batch = matched[i : i + args.batch_size]
            result = bulk_insert(args.es_url, args.index, args.auth_token, batch, insecure=args.insecure)
            items = result.get("items", [])
            errors_in_batch = [it for it in items if "index" in it and it["index"].get("status", 200) >= 300]
            total_indexed += len(items) - len(errors_in_batch)
            total_errors += len(errors_in_batch)
            print(
                f"Batch {i // args.batch_size + 1}: sent {len(batch)}, "
                f"errors {len(errors_in_batch)}",
                file=sys.stderr,
            )
            for err in errors_in_batch[:5]:
                print(f"    error: {json.dumps(err)}", file=sys.stderr)
    finally:
        # Guarantee the index is set back to read-only even if the script fails or exits early
        if matched:
            print(f"\nRestoring read-only state for index '{args.index}'...", file=sys.stderr)
            update_index_settings(
                args.es_url, args.index, args.auth_token,
                {"index.blocks.write": True},
                insecure=args.insecure
            )

    print(f"\nDone. Indexed OK: {total_indexed}, Errors: {total_errors}, Unmatched (skipped): {len(missing)}", file=sys.stderr)


if __name__ == "__main__":
    main()
