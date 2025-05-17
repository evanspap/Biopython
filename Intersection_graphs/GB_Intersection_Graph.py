import argparse
import os
import csv

# =============================================================
# Script: Intersection Graph Table Generator
# Description:
#   Processes a GenBank (.gb) file, extracts residues from multi-line
#   join() annotations, computes intersection counts, and dumps
#   metadata qualifiers as a single comma-separated field at end of each row.
# Version: 1.0.10
# Date: 2025-05-13
# Usage: python GB_Intersection_Graph.py <input_file> <output_file>
# Example: python GB_Intersection_Graph.py TLE3.gb results/Intersection_Graph_Table.csv
# =============================================================

def extract_positions(join_text):
    # Extract integers from join(...) possibly spanning multiple lines
    inner = join_text[join_text.find('join(') + 5:]
    if ')' in inner:
        inner = inner[:inner.rfind(')')]
    parts = inner.replace('(', '').replace(')', '').split(',')
    positions = set()
    for part in parts:
        nums = ''.join(ch for ch in part if ch.isdigit())
        if nums:
            positions.add(int(nums))
    return positions


def extract_metadata(lines, idx):
    # After a join(...) block, capture lines starting with '/'
    qualifiers = []
    j = idx
    while j + 1 < len(lines):
        nxt = lines[j+1].strip()
        if nxt.startswith('/'):
            # strip leading slash, parse key=value without quotes
            key_val = nxt.lstrip('/').split('=', 1)
            key = key_val[0]
            val = key_val[1].strip().strip('"') if len(key_val) > 1 else ''
            qualifiers.append(f"{key}={val}")
            j += 1
        else:
            break
    return qualifiers, j


def main():
    parser = argparse.ArgumentParser(description="Generate intersection graph table with metadata dumped per pocket.")
    parser.add_argument("input_file", help="Path to GenBank input file (.gb)")
    parser.add_argument("output_file", help="Path to CSV output file")
    args = parser.parse_args()

    # Ensure output directory exists
    out_dir = os.path.dirname(args.output_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Read file
    with open(args.input_file) as fh:
        lines = fh.readlines()

    # Parse join() and metadata
    pockets = {}       # pocket_name -> set of residues
    pocket_meta = {}   # pocket_name -> list of key=value strings

    i = 0
    while i < len(lines):
        line = lines[i]
        if 'join(' in line:
            pocket = line.strip().split()[0]
            join_txt = line.strip()
            # accumulate multi-line join
            open_p = join_txt.count('(')
            close_p = join_txt.count(')')
            j = i
            while open_p > close_p and j + 1 < len(lines):
                j += 1
                nxt = lines[j].strip()
                join_txt += nxt
                open_p += nxt.count('(')
                close_p += nxt.count(')')
            # extract residues
            pockets[pocket] = extract_positions(join_txt)
            # extract metadata qualifiers
            qualifiers, end_idx = extract_metadata(lines, j)
            pocket_meta[pocket] = qualifiers
            i = end_idx + 1
        else:
            i += 1

    # Build intersection counts
    pocket_names = list(pockets)
    total = {p: len(pockets[p]) for p in pocket_names}
    # prepare matrix p1->p2 count
    matrix = {p1: {p2: (total[p1] if p1==p2 else len(pockets[p1] & pockets[p2]))
                   for p2 in pocket_names}
              for p1 in pocket_names}

        # Write CSV: header + rows using raw comma-joined lines
    header = ['PocketName', 'Total Residues'] + pocket_names + ['Metadata']
    with open(args.output_file, 'w') as out_file:
        # Write header line
        out_file.write(','.join(header) + '\n')
        for p1 in pocket_names:
            row = [p1, str(total[p1])]
            row.extend(str(matrix[p1][p2]) for p2 in pocket_names)
            # metadata qualifiers as comma-separated
            meta_str = ', '.join(pocket_meta.get(p1, []))
            row.append(meta_str)
            # write each data row
            out_file.write(','.join(row) + '\n')

        # Finish
    print(f"Intersection graph table saved to {args.output_file}")

if __name__ == '__main__':
    main()
