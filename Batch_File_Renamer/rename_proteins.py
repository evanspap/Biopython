#!/usr/bin/env python3

"""
Batch File Renamer for Protein Files

This script renames all files with a specified extension in a given folder to a sequential format:
    protein_001.ext, protein_002.ext, etc.

Usage:
    python rename_proteins.py <folder_path> <extension> [--dry-run]

Arguments:
    <folder_path>  - Path to the folder containing files to rename.
    <extension>    - File extension (without the dot) to filter files for renaming.
    --dry-run      - Optional flag to preview renaming without applying changes.

Example:
    python rename_proteins.py /path/to/folder txt --dry-run
    python rename_proteins.py /path/to/folder pdb
"""

import os
import sys

def rename_files(folder_path, extension, dry_run):
    if not os.path.isdir(folder_path):
        print(f"Error: Folder '{folder_path}' does not exist.")
        sys.exit(1)

    files = sorted([f for f in os.listdir(folder_path) if f.endswith(f".{extension}")])
    
    if not files:
        print(f"No files with extension '.{extension}' found in '{folder_path}'.")
        sys.exit(1)

    print("Renaming files...")
    
    for i, old_name in enumerate(files, start=1):
        new_name = f"protein_{i:03d}.{extension}"
        old_path = os.path.join(folder_path, old_name)
        new_path = os.path.join(folder_path, new_name)

        if dry_run:
            print(f"Dry Run: {old_name} -> {new_name}")
        else:
            print(f"Renaming: {old_name} -> {new_name}")
            os.rename(old_path, new_path)

    print("Done.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    folder_path = sys.argv[1]
    extension = sys.argv[2]
    dry_run = "--dry-run" in sys.argv

    rename_files(folder_path, extension, dry_run)
