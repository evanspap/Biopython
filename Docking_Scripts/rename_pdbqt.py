#!/usr/bin/env python3

import sys
import os
import glob

def header():
    print('''
Script to rename all PDBQT files in a folder based on the first "REMARK  Name =" line.

Usage:
    python rename_pdbqt.py <input_folder> [--dry-run | --no-dry-run]

Example:
    python rename_pdbqt.py ./pdbqt_files
    python rename_pdbqt.py ./pdbqt_files --dry-run
    python rename_pdbqt.py ./pdbqt_files --no-dry-run
''')

if len(sys.argv) < 2 or len(sys.argv) > 3:
    header()
    sys.exit(1)

input_folder = sys.argv[1]
DRY_RUN = True  # Default

if len(sys.argv) == 3:
    if sys.argv[2] == '--no-dry-run':
        DRY_RUN = False
    elif sys.argv[2] != '--dry-run':
        header()
        sys.exit(1)

if not os.path.isdir(input_folder):
    print(f"Error: Folder '{input_folder}' does not exist.")
    sys.exit(1)

pdbqt_files = glob.glob(os.path.join(input_folder, '*.pdbqt'))

if not pdbqt_files:
    print("No '.pdbqt' files found in the specified folder.")
    sys.exit(1)

for input_file in pdbqt_files:
    new_name = None
    with open(input_file, 'r') as file:
        for line in file:
            if line.startswith('REMARK  Name ='):
                new_name = line.split('=', 1)[1].strip().replace(' ', '_') + '.pdbqt'
                break

    if not new_name:
        print(f"Warning: No 'REMARK  Name =' line found in file '{input_file}'. Skipping.")
        continue

    target_path = os.path.join(input_folder, new_name)

    if os.path.isfile(target_path):
        print(f"Warning: Target filename '{new_name}' already exists. Skipping '{input_file}'.")
        continue

    print(f"Running: Renaming '{input_file}' to '{target_path}'")

    if DRY_RUN:
        print("Dry run: no actual file operation performed.")
    else:
        os.rename(input_file, target_path)
        print("File renamed successfully.")
