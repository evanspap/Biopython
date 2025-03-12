import pandas as pd
import argparse
import os

# ==============================================================
# Script: Intersection Percentage Table Generator
# Description: This script processes a GenBank (.gb) file, extracts
#              residue positions from all annotations, and generates
#              an intersection percentage table based on shared residues.
# Usage: python GB_Intersection_Percentage.py <input_file> <output_file>
# Example: python GB_Intersection_Percentage.py ADSS2.gb results/Intersection_Percentage_Table.csv
# ==============================================================

# Function to extract residue positions from a join() field
def extract_positions(line):
    positions = line.split("join(")[-1].replace(")", "").split(",")
    return set(map(int, filter(str.isdigit, positions)))

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

# Check if output_file is actually a file (not a directory)
if os.path.isdir(output_file):
    raise ValueError(f"Error: The specified output path '{output_file}' is a directory, not a file.")

pockets = {}
print("Processing file:", genbank_file)

with open(genbank_file, "r") as file:
    for line in file:
        line = line.strip()
        if "join(" in line:
            pocket_name = line.split()[0]  # First word is the pocket name
            positions = extract_positions(line)
            pockets[pocket_name] = positions

# Create an intersection percentage table
pocket_names = list(pockets.keys())
total_residues = {p: len(pockets[p]) for p in pocket_names}

intersection_percentage_matrix = {p1: {p2: 0 for p2 in pocket_names} for p1 in pocket_names}

# Compute intersection percentages (dividing by pockets[p2])
for p1 in pocket_names:
    for p2 in pocket_names:
        if p1 == p2:
            intersection_percentage_matrix[p1][p2] = 100  # Self-comparison should be 100%
        else:
            common_residues = pockets[p1] & pockets[p2]
            if len(pockets[p2]) > 0:
                intersection_percentage_matrix[p1][p2] = int((len(common_residues) / len(pockets[p2])) * 100)

# Convert matrix into a DataFrame
intersection_percentage_df = pd.DataFrame(intersection_percentage_matrix).T

# Add total residues column
intersection_percentage_df.insert(0, "Total Residues", [total_residues[p] for p in pocket_names])

# Add the "Total Residues" row at the top
total_residue_row = ["-"] + [total_residues[p] for p in pocket_names]
intersection_percentage_df.loc["Total Residues"] = total_residue_row

# Reorder to move "Total Residues" to the first row
intersection_percentage_df = intersection_percentage_df.iloc[[-1] + list(range(len(intersection_percentage_df) - 1))]

# Save the output to the specified file
intersection_percentage_df.to_csv(output_file, index=True)

print(f"Intersection percentage table saved to {output_file}")
