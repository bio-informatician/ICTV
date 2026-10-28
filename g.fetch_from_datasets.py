#!/usr/bin/env python3
"""Resolve viral GenBank accessions to the latest version, retrieve NCBI Datasets metadata, and stream normalized JSONL."""

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path

from tqdm import tqdm


def setup_logger(log_file: Path) -> logging.Logger:
    logger = logging.getLogger("virus_fetch")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    fmt = logging.Formatter(
        "%(asctime)s | %(levelname)-8s | %(message)s",
        "%Y-%m-%d %H:%M:%S",
    )

    sh = logging.StreamHandler(sys.stderr)
    sh.setLevel(logging.INFO)
    sh.setFormatter(fmt)

    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.WARNING)
    fh.setFormatter(fmt)

    logger.addHandler(sh)
    logger.addHandler(fh)

    return logger


def run(cmd):
    p = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
    )

    if p.returncode != 0:
        raise RuntimeError(p.stderr.strip())

    return p.stdout


def latest_accession(acc: str) -> str:
    cmd = f"""
esearch -db nucleotide -query "{acc}" 2>/dev/null |
efetch -format docsum 2>/dev/null |
xtract -pattern DocumentSummary -element AccessionVersion |
sort -V |
tail -n1
"""

    p = subprocess.run(
        cmd,
        shell=True,
        executable="/bin/bash",
        capture_output=True,
        text=True,
    )

    version = p.stdout.strip()

    if not version:
        raise RuntimeError("unable to resolve accession")

    return version


def lineage_lookup(lineage):
    return {x["name"]: x for x in lineage}


def normalize(report):
    L = lineage_lookup(report["virus"]["lineage"])

    return {
        "insdc_accession": report.get("accession"),
        "bioproject_accession": report.get("bioprojects", []),
        "genome_completeness": report.get("completeness"),
        "sequence_length": report.get("length"),
        "release_date": report.get("release_date", "").split("T")[0],
        "update_date_refseq": report.get("update_date", "").split("T")[0],
        "gene_count": report.get("gene_count"),
        "is_annotated_refseq": report.get("is_annotated"),
        "protein_count": report.get("protein_count"),
        "source_database": report.get("source_database"),
        "collection_date": report.get("isolate", {}).get("collection_date"),
        "collection_country": report.get("isolate", {}).get("geo_location"),
        "host_common_name": report.get("host", {}).get("common_name"),
        "host_ncbi_tax_id": report.get("host", {}).get("tax_id"),
        "virus_organism_name": report["virus"].get("organism_name"),
        "submitter": report.get("submitter", {}).get("names"),
        "submitter_organization": report.get("submitter", {}).get("affiliation"),
        "geographic_region": report.get("isolate", {}).get("geographic_region"),
        "sequence_hash_refseq": report.get("nucleotide", {}).get("sequence_hash"),
        "submitter_country": report.get("submitter", {}).get("country"),
        "virus_tax_id_ncbi": report["virus"].get("tax_id"),
        "realm_ncbi": L.get("Riboviria", {}).get("name"),
        "kingdom_ncbi": L.get("Orthornavirae", {}).get("name"),
        "phylum_ncbi": L.get("Pisuviricota", {}).get("name"),
        "class_ncbi": L.get("Pisoniviricetes", {}).get("name"),
        "order_ncbi": L.get("Picornavirales", {}).get("name"),
        "family_ncbi": L.get("Iflaviridae", {}).get("name"),
        "genus_ncbi": L.get("Iflavirus", {}).get("name"),
        "species_ncbi": L.get("Iflavirus flacherie", {}).get("name"),
        "host_species_ncbi": report.get("host", {}).get("species"),
        "host_genus_ncbi": report.get("host", {}).get("genus"),
        "host_family_ncbi": report.get("host", {}).get("family"),
    }


def fetch(acc: str):
    version = latest_accession(acc)

    raw = run(
        [
            "datasets",
            "summary",
            "virus",
            "genome",
            "accession",
            version,
        ]
    )

    obj = json.loads(raw)

    if obj.get("total_count", 0) == 0:
        raise RuntimeError("datasets returned no report")

    return normalize(obj["reports"][0])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("-o", "--output", default="virus_metadata.jsonl")
    parser.add_argument("--log", default="virus_metadata.log")

    args = parser.parse_args()

    logger = setup_logger(Path(args.log))

    if Path(args.input).is_file():
        accessions = [
            x.strip()
            for x in Path(args.input).read_text().splitlines()
            if x.strip()
        ]
    else:
        accessions = [args.input]

    logger.info("Processing %d accession(s)", len(accessions))

    success = 0
    failed = 0

    with open(args.output, "w") as fout:
        for acc in tqdm(accessions, unit="acc"):
            try:
                result = fetch(acc)

                fout.write(json.dumps(result, separators=(",", ":")) + "\n")
                fout.flush()

                success += 1

            except Exception as e:
                failed += 1
                logger.warning("%s | %s", acc, e)

    logger.info(
        "Finished | success=%d failed=%d output=%s",
        success,
        failed,
        args.output,
    )


if __name__ == "__main__":
    main()
