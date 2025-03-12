#!/usr/bin/env python3

"""
Script to extract numerical values from *.pdbqt files.

The script scans a specified directory for pdbqt files, extracts the numeric value from the line containing
"REMARK INTER + INTRA:" (typically representing docking scores), and compiles the results into a CSV file.

Usage:
    python extract_pdbqt_docking_values.py <directory>

Arguments:
    <directory>  Path to the directory containing *.pdbqt files.

Example:
    extract_pdbqt_docking_values.py ./pdbqt_files
"""

import os
import sys
import glob
import csv

def print_header():
    header = """
Extract numerical values from *.pdbqt files.

Usage:
    python extract_pdbqt_values.py <directory>

Arguments:
    <directory>  Path to the directory containing *.pdbqt files.

Example:
    python extract_pdbqt_values.py ./pdbqt_files
"""
    print(header)

def extract_value_from_pdbqt(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            if len(lines) > 0:
                for line in lines:
                    if "REMARK INTER + INTRA:" in line:
                        value = line.strip().split()[-1]
                        return float(value)
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
    return None

def main():
    if len(sys.argv) != 2:
        print_header()
        sys.exit(1)

    directory = sys.argv[1]

    if not os.path.isdir(directory):
        print(f"Error: Directory '{directory}' does not exist.")
        sys.exit(1)

    pdbqt_files = glob.glob(os.path.join(directory, "*.pdbqt"))

    if not pdbqt_files:
        print(f"No *.pdbqt files found in '{directory}'.")
        sys.exit(1)

    csv_output = os.path.join(directory, "pdbqt_summary.csv")

    with open(csv_output, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Filename", "Value"])

        for file_path in pdbqt_files:
            value = extract_value_from_pdbqt(file_path)
            filename = os.path.basename(file_path)

            if value is not None:
                writer.writerow([filename, value])
                print(f"Processed {filename}: {value}")
            else:
                print(f"Value not found or invalid format in {filename}")

    print(f"\nCSV summary saved to {csv_output}")

if __name__ == "__main__":
    main()
