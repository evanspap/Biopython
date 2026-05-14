#!/usr/bin/env python3
"""
Add UniProt entry names and primary gene names to a CSV/TSV file.

Recommended runtime: Ubuntu/WSL with python3. This avoids Windows PowerShell
redirect encoding issues and writes clean UTF-8 text to standard output.

The script reads UniProt accessions from a selected input column, queries the
UniProt REST API, and writes the enriched table to standard output so it can be
redirected to a new file. Progress messages and warnings are written to standard
error. The added columns are UniProt_entry_name and Gene_name.

If no command-line arguments are provided, the script prints help and exits.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
DEFAULT_INPUT = (
    r"C:\Users\geras\Partners HealthCare Dropbox\Evangelos Papadopoulos"
    r"\LAB-Evan\Manuscript MMSEG pockets\AF2BIND\analysis_files\analysis_files"
    r"\af2bind_p2rank_human_proteome_predictions_revision_with_sitemap.csv"
)
DEFAULT_INPUT_WSL = (
    "/mnt/c/Users/geras/Partners HealthCare Dropbox/Evangelos Papadopoulos"
    "/LAB-Evan/Manuscript MMSEG pockets/AF2BIND/analysis_files/analysis_files"
    "/af2bind_p2rank_human_proteome_predictions_revision_with_sitemap.csv"
)


def parse_args() -> argparse.Namespace:
    default_input = DEFAULT_INPUT if os.name == "nt" else DEFAULT_INPUT_WSL
    parser = argparse.ArgumentParser(
        description=(
            "Read UniProt accessions from a CSV/TSV column and add UniProt "
            "entry name plus primary gene name columns to standard output."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        default=default_input,
        help="Input CSV/TSV file. Defaults to the AF2BIND revision CSV.",
    )
    parser.add_argument(
        "-c",
        "--column",
        default="Uniprot",
        help="Column containing UniProt accessions. Matching is case-insensitive.",
    )
    parser.add_argument(
        "--delimiter",
        default=None,
        help="Input delimiter. By default it is inferred from the file extension.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=80,
        help="Number of unique UniProt accessions to query per API request.",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.15,
        help="Seconds to wait between UniProt API requests.",
    )
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        raise SystemExit(0)
    return parser.parse_args()


def infer_delimiter(path: Path, explicit_delimiter: str | None) -> str:
    if explicit_delimiter:
        return explicit_delimiter
    if path.suffix.lower() == ".tsv":
        return "\t"
    return ","


def find_column(fieldnames: list[str], requested_column: str) -> str:
    normalized = {name.strip().lower(): name for name in fieldnames}
    key = requested_column.strip().lower()
    if key not in normalized:
        available = ", ".join(fieldnames)
        raise ValueError(
            f"Could not find column '{requested_column}'. Available columns: {available}"
        )
    return normalized[key]


def clean_accession(value: str | None) -> str:
    if not value:
        return ""
    return value.strip()


def chunks(values: list[str], size: int) -> list[list[str]]:
    return [values[index : index + size] for index in range(0, len(values), size)]


def fetch_uniprot_batch(accessions: list[str], timeout: int = 60) -> dict[str, tuple[str, str]]:
    query = " OR ".join(f"accession:{accession}" for accession in accessions)
    params = {
        "query": f"({query})",
        "fields": "accession,id,gene_primary",
        "format": "tsv",
        "size": str(len(accessions)),
    }
    url = f"{UNIPROT_SEARCH_URL}?{urllib.parse.urlencode(params)}"
    request = urllib.request.Request(url, headers={"User-Agent": "add-uniprot-entry-gene-names/1.0"})

    with urllib.request.urlopen(request, timeout=timeout) as response:
        text = response.read().decode("utf-8")

    mapping: dict[str, tuple[str, str]] = {}
    reader = csv.DictReader(text.splitlines(), delimiter="\t")
    for row in reader:
        accession = row.get("Entry", "").strip()
        entry_name = row.get("Entry Name", "").strip()
        gene_name = row.get("Gene Names (primary)", "").strip()
        if accession:
            mapping[accession] = (entry_name, gene_name)
    return mapping


def fetch_uniprot_mapping(
    accessions: list[str], batch_size: int, sleep_seconds: float
) -> dict[str, tuple[str, str]]:
    mapping: dict[str, tuple[str, str]] = {}
    total_batches = len(chunks(accessions, batch_size))

    for batch_index, batch in enumerate(chunks(accessions, batch_size), start=1):
        try:
            mapping.update(fetch_uniprot_batch(batch))
        except urllib.error.URLError as exc:
            print(
                f"WARNING: UniProt request failed for batch {batch_index}/{total_batches}: {exc}",
                file=sys.stderr,
            )
        if sleep_seconds and batch_index < total_batches:
            time.sleep(sleep_seconds)

    return mapping


def read_rows(input_path: Path, delimiter: str) -> tuple[list[dict[str, str]], list[str]]:
    with input_path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Input file has no header: {input_path}")
        return list(reader), list(reader.fieldnames)


def write_rows_to_stdout(
    delimiter: str, rows: list[dict[str, str]], fieldnames: list[str]
) -> None:
    writer = csv.DictWriter(
        sys.stdout, fieldnames=fieldnames, delimiter=delimiter, lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    delimiter = infer_delimiter(input_path, args.delimiter)

    rows, fieldnames = read_rows(input_path, delimiter)
    uniprot_column = find_column(fieldnames, args.column)

    accessions = sorted(
        {
            clean_accession(row.get(uniprot_column))
            for row in rows
            if clean_accession(row.get(uniprot_column))
        }
    )
    print(f"Found {len(accessions)} unique UniProt accessions.", file=sys.stderr)

    mapping = fetch_uniprot_mapping(accessions, args.batch_size, args.sleep)
    print(f"Retrieved UniProt metadata for {len(mapping)} accessions.", file=sys.stderr)

    for row in rows:
        accession = clean_accession(row.get(uniprot_column))
        entry_name, gene_name = mapping.get(accession, ("", ""))
        row["UniProt_entry_name"] = entry_name
        row["Gene_name"] = gene_name

    output_fieldnames = [
        name
        for name in fieldnames
        if name not in {"UniProt_entry_name", "Gene_name"}
    ] + ["UniProt_entry_name", "Gene_name"]
    write_rows_to_stdout(delimiter, rows, output_fieldnames)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
