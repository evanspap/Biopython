#!/usr/bin/env python3
"""
===============================================================
Pocket-to-GenBank Converter
===============================================================
Author : Evangelos Papadopoulos
Date   : 2025-11-06
Version: 1.4

Description:
    Reads a CSV/TSV file of protein pocket data and generates
    one GenBank (.gb) file per protein.

New in v1.4:
    • Added --outdir <folder> argument to choose output directory.

Features:
    • Auto-detects CSV/TSV delimiters
    • Cleans header names (removes BOM/spaces)
    • Groups pockets by protein
    • Optional numbering of pocket features (--pocket-names=numbered)
    • Creates join() features from residue positions
    • Adds a misc_feature covering the full protein
    • Writes each protein's GenBank file to chosen folder
    • Supports --dry-run mode

Usage:
    python pocket_to_genbank.py <input_file.csv> [--dry-run]
                                [--pocket-names=numbered]
                                [--outdir <output_folder>]

Example:
    python pocket_to_genbank.py top_pockets.csv
    python pocket_to_genbank.py top_pockets.csv --pocket-names=numbered --outdir ./GenBanks
"""

import sys, os, csv
from collections import defaultdict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

# ------------------------------------------------------------
# --- Parse arguments ---
# ------------------------------------------------------------
if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(1)

input_file = sys.argv[1]
dry_run = "--dry-run" in sys.argv
numbered_flag = any("--pocket-names=numbered" in arg for arg in sys.argv)

# Detect optional output directory
if "--outdir" in sys.argv:
    try:
        outdir_index = sys.argv.index("--outdir")
        outdir = sys.argv[outdir_index + 1]
    except IndexError:
        print("Error: --outdir flag requires a folder path.")
        sys.exit(1)
else:
    outdir = "."

# Ensure directory exists
if not os.path.exists(outdir):
    os.makedirs(outdir, exist_ok=True)

if not os.path.exists(input_file):
    print(f"Error: Input file '{input_file}' not found.")
    sys.exit(1)

print(f"Running: python pocket_to_genbank.py {input_file}"
      + (" --dry-run" if dry_run else "")
      + (" --pocket-names=numbered" if numbered_flag else "")
      + (f" --outdir {outdir}" if outdir else ""))

# ------------------------------------------------------------
# --- Read CSV data (auto-detect delimiter) ---
# ------------------------------------------------------------
proteins = defaultdict(list)
with open(input_file, newline='', encoding='utf-8-sig') as f:
    sample = f.read(4096)
    f.seek(0)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
    except csv.Error:
        dialect = csv.get_dialect("excel")
    reader = csv.DictReader(f, dialect=dialect)
    if reader.fieldnames:
        reader.fieldnames = [h.strip() for h in reader.fieldnames]

    for row in reader:
        if not row:
            continue
        if "Protein" not in row:
            print("⚠️  Warning: 'Protein' column not found in this row, skipping.")
            continue
        prot = row["Protein"].strip()
        proteins[prot].append(row)

print(f"Found {len(proteins)} unique proteins in {input_file}")

# ------------------------------------------------------------
# --- Generate GenBank files ---
# ------------------------------------------------------------
for prot, entries in proteins.items():
    record = SeqRecord(
        Seq(""),
        id=prot,
        name=prot,
        description=f"Predicted pockets for {prot}"
    )
    record.annotations["molecule_type"] = "protein"

    # --- Add top-level misc_feature covering full sequence (placeholder) ---
    seq_len = len(record.seq) if len(record.seq) > 0 else 1
    full_feature = SeqFeature(
        FeatureLocation(0, seq_len),
        type="misc_feature",
        qualifiers={"note": f"Full protein: {prot}"}
    )
    record.features.append(full_feature)

    # --- Add pocket features ---
    for i, row in enumerate(entries, start=1):
        pocket_name = row.get("Pocket", f"{prot}_P{i}")
        residues_str = row.get("Residue list", "").strip("[] ")
        residues = [int(x) for x in residues_str.replace(",", " ").split() if x.isdigit()]

        qualifiers = {
            "product": pocket_name,
            "note": f"Residues: {residues_str}",
        }
        for key in [
            "max%",
            "Product AS*FP*DG",
            "Product FP*DG",
            "FP DScore",
            "DG DScore",
            "ASe N_CDF",
            "common residues",
        ]:
            if key in row and row[key].strip():
                clean_key = key.replace(" ", "_").replace("*", "_")
                qualifiers[clean_key] = row[key].strip()

        # Build join() feature
        if residues:
            locations = [FeatureLocation(r - 1, r) for r in residues]
            location = sum(locations[1:], locations[0]) if len(locations) > 1 else locations[0]
        else:
            location = FeatureLocation(0, 0)

        feature_type = f"pocket{str(i).zfill(2)}" if numbered_flag else "pocket"
        feature = SeqFeature(location=location, type=feature_type, qualifiers=qualifiers)
        record.features.append(feature)

    # --- Write file ---
    if dry_run:
        print(f"Would create: {os.path.join(outdir, f'{prot}.gb')} with {len(entries)} pockets (+ full misc_feature)")
    else:
        outfile = os.path.join(outdir, f"{prot}.gb")
        with open(outfile, "w") as out_handle:
            SeqIO.write(record, out_handle, "genbank")
        print(f"Wrote {outfile} with {len(entries)} pockets (+ full misc_feature)")

print("✅ Done.")

