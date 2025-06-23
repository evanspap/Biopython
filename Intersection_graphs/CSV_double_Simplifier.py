# ==============================================================
# Script: CSV_double_Simplifier.py
# Version: 1.5
# Date: 2025-05-27
# Description:
#   Reads a CSV of pairwise pocket intersections, selects rows and columns matching
#   two pocket-base substrings (case-insensitive), produces a simplified symmetric
#   matrix of their intersections, and appends the 'Metadata' column if present.
#
# Usage:
#   python CSV_double_Simplifier.py <input_csv> <pocket_base1> <pocket_base2> <output_csv>
#
# Example:
#   python CSV_double_Simplifier.py ADSS2.csv ASpock F_Pock ADSS2_simple.csv
# ==============================================================

import csv
import os
import sys
import argparse
from datetime import datetime
import pandas as pd

__version__ = "1.5"
__date__ = datetime.now().strftime('%Y-%m-%d')

def main():
    parser = argparse.ArgumentParser(
        description="Simplify a pairwise pocket CSV to two bases and metadata."
    )
    parser.add_argument("input_csv", help="Input CSV file path")
    parser.add_argument("pocket_base1", help="First pocket-base substring (e.g. 'ASpock')")
    parser.add_argument("pocket_base2", help="Second pocket-base substring (e.g. 'F_Pock')")
    parser.add_argument("output_csv", help="Output CSV file path")
    args = parser.parse_args()

    # Validate input
    if not os.path.isfile(args.input_csv):
        print(f"Error: '{args.input_csv}' not found.")
        sys.exit(1)

    # Read CSV robustly
    with open(args.input_csv, newline='') as f:
        reader = csv.reader(f)
        rows = [r for r in reader if r]
    if not rows:
        print("Error: input CSV is empty or malformed.")
        sys.exit(1)

    header = rows[0]
    ncols = len(header)
    data = []
    for row in rows[1:]:
        if len(row) > ncols:
            row = row[:ncols-1] + [','.join(row[ncols-1:])]
        elif len(row) < ncols:
            row = row + ['']*(ncols - len(row))
        data.append(row)

    df = pd.DataFrame(data, columns=header)
    idx_col = header[0]
    df.set_index(idx_col, inplace=True)

    # Prepare case-insensitive bases
    bases = [args.pocket_base1.lower(), args.pocket_base2.lower()]

    # Select pocket rows
    pocket_rows = [r for r in df.index if any(b in r.lower() for b in bases)]
    if not pocket_rows:
        print("No rows match the given pocket bases.")
        sys.exit(1)

    # Select pocket columns
    pocket_cols = [c for c in df.columns if any(b in c.lower() for b in bases)]
    if not pocket_cols:
        print("No columns match the given pocket bases.")
        sys.exit(1)

    # Check for 'Metadata' column
    metadata_col = 'Metadata' if 'Metadata' in df.columns else None

    # Build simplified records
    records = []
    for r in pocket_rows:
        rec = {idx_col: r}
        # Diagonalized intersections
        for c in pocket_cols:
            v1 = df.at[r, c] if c in df.columns else ''
            v2 = df.at[c, r] if r in df.columns else ''
            # parse numbers
            try:
                num = max(int(v1), int(v2))
            except:
                try:
                    num = max(float(v1) if v1 else 0, float(v2) if v2 else 0)
                except:
                    num = 0
            rec[c] = num
        # Append metadata
        if metadata_col:
            rec['Metadata'] = df.at[r, metadata_col]
        records.append(rec)

    # Construct output DataFrame
    out_cols = [idx_col] + pocket_cols + (['Metadata'] if metadata_col else [])
    out_df = pd.DataFrame(records, columns=out_cols)

    # Save
    out_df.to_csv(args.output_csv, index=False)
    print(f"Output written to {args.output_csv} (v{__version__} - {__date__})")

if __name__ == '__main__':
    main()
