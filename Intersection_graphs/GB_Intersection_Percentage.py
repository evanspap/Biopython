import pandas as pd
import argparse
import os

# ==============================================================
# Script: Intersection Percentage Table Generator
# Description: This script processes a GenBank (.gb) file, extracts
# Version: 1.0.1
# Date: 2025-05-13
#              residue positions from all annotations, and generates
#              an intersection percentage table based on shared residues.
# Usage: python GB_Intersection_Percentage.py <input_file> <output_file>
# Example: python GB_Intersection_Percentage.py ADSS2.gb results/Intersection_Percentage_Table.csv
# ==============================================================

# Function to extract residue positions from a join() string
def extract_positions(join_text):
    # Remove leading 'join(' and trailing ')' then split by commas
    inner = join_text[join_text.find('join(') + 5:]
    # Remove any trailing chars after the last ')'
    if ')' in inner:
        inner = inner[:inner.rfind(')')]
    parts = inner.replace('(', '').replace(')', '').split(',')
    positions = set()
    for part in parts:
        # Extract digits
        nums = ''.join(ch for ch in part if ch.isdigit())
        if nums:
            positions.add(int(nums))
    return positions

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Generate an intersection percentage table for all annotation residues.")
parser.add_argument("input_file", help="Path to the GenBank input file")
parser.add_argument("output_file", help="Full path to the output CSV file (including filename)")
args = parser.parse_args()

genbank_file = args.input_file
output_file = args.output_file

# Ensure the output directory exists
output_dir = os.path.dirname(output_file)
if output_dir:
    os.makedirs(output_dir, exist_ok=True)
# Check if output_file is a directory
if os.path.isdir(output_file):
    raise ValueError(f"Error: The specified output path '{output_file}' is a directory, not a file.")

pockets = {}
print("Processing file:", genbank_file)

# Read file lines for multi-line join handling
with open(genbank_file, 'r') as fh:
    lines = fh.readlines()
i = 0
while i < len(lines):
    raw = lines[i]
    if 'join(' in raw:
        line = raw.strip()
        # The pocket name is the first token
        pocket_name = line.split()[0]
        # Accumulate join text until parentheses close
        join_text = line
        open_par = join_text.count('(')
        close_par = join_text.count(')')
        j = i
        while open_par > close_par and j + 1 < len(lines):
            j += 1
            next_line = lines[j].strip()
            join_text += next_line
            open_par += next_line.count('(')
            close_par += next_line.count(')')
        # Extract positions
        positions = extract_positions(join_text)
        pockets[pocket_name] = positions
        i = j + 1
    else:
        i += 1

# Create an intersection percentage table
pocket_names = list(pockets.keys())
total_residues = {p: len(pockets[p]) for p in pocket_names}

intersection_percentage_matrix = {p1: {p2: 0 for p2 in pocket_names} for p1 in pocket_names}

# Compute intersection percentages (dividing by pockets[p2])
for p1 in pocket_names:
    for p2 in pocket_names:
        if p1 == p2:
            intersection_percentage_matrix[p1][p2] = 100
        else:
            common = pockets[p1] & pockets[p2]
            if pockets[p2]:
                intersection_percentage_matrix[p1][p2] = int(len(common) / len(pockets[p2]) * 100)

# Convert matrix into a DataFrame
import pandas as pd
intersection_percentage_df = pd.DataFrame(intersection_percentage_matrix).T

# Add total residues column
intersection_percentage_df.insert(0, "Total Residues", [total_residues[p] for p in pocket_names])

# Add the "Total Residues" row at the top
total_row = ["-"] + [total_residues[p] for p in pocket_names]
intersection_percentage_df.loc["Total Residues"] = total_row
# Reorder to move it to first row
intersection_percentage_df = intersection_percentage_df.iloc[[-1] + list(range(len(intersection_percentage_df)-1))]

# Save to CSV
intersection_percentage_df.to_csv(output_file, index=True)
print(f"Intersection percentage table saved to {output_file}")
