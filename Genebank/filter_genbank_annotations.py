#!/usr/bin/env python3
"""
filter_genbank_annotations.py
Date: 2025-05-27
Version: 1.3
Author: EP

Script to retain only GenBank features whose feature type, qualifier names, or qualifier values
start with specified prefixes, and write filtered files to an output directory.

Usage:
    python filter_genbank_annotations.py <gb_directory> <output_directory> <prefix1> [<prefix2> ...] [--dry-run]

Example:
    python filter_genbank_annotations.py ./genbank_files ./filtered ASpock Fpock --dry-run
"""
import os
import sys
from Bio import SeqIO

def print_usage_and_exit():
    print(__doc__)
    sys.exit(1)


def filter_features(record, prefixes):
    """
    Keep only features whose feature.type, qualifier names, or qualifier values start with any of the prefixes.
    """
    kept = []
    for feature in record.features:
        match = False
        for prefix in prefixes:
            # Check feature type
            if feature.type.startswith(prefix):
                match = True
            # Check qualifier names
            if not match:
                for qual_name in feature.qualifiers:
                    if qual_name.startswith(prefix):
                        match = True
                        break
            # Check qualifier values
            if not match:
                for values in feature.qualifiers.values():
                    for v in values:
                        if v.startswith(prefix):
                            match = True
                            break
                    if match:
                        break
            if match:
                kept.append(feature)
                break
        # end for prefixes
    return kept


def main():
    args = sys.argv[1:]
    dry_run = False

    if '--dry-run' in args:
        dry_run = True
        args.remove('--dry-run')

    if len(args) < 3:
        print_usage_and_exit()

    gb_dir = args[0]
    out_dir = args[1]
    prefixes = args[2:]

    if not os.path.isdir(gb_dir):
        print(f"Error: '{gb_dir}' is not a directory.")
        sys.exit(1)

    if not dry_run:
        os.makedirs(out_dir, exist_ok=True)

    for fname in os.listdir(gb_dir):
        if not fname.lower().endswith('.gb'):
            continue
        in_path = os.path.join(gb_dir, fname)
        print(f"Running: Processing {in_path}")
        try:
            record = SeqIO.read(in_path, 'genbank')
        except Exception as e:
            print(f"Error reading {in_path}: {e}")
            continue

        kept = filter_features(record, prefixes)
        print(f"Running: Kept {len(kept)} of {len(record.features)} features in {fname}")

        if dry_run:
            print(f"Running: Dry run enabled, not writing filtered file for {fname}\n")
            continue

        record.features = kept
        out_path = os.path.join(out_dir, fname)
        print(f"Running: Writing filtered record to {out_path}")
        try:
            SeqIO.write(record, out_path, 'genbank')
            print(f"Running: Wrote filtered record to {out_path}\n")
        except Exception as e:
            print(f"Error writing {out_path}: {e}\n")

if __name__ == '__main__':
    main()
