#!/usr/bin/env python3

"""
CSV Heatmap Generator
=====================
This script reads a CSV file containing data for a heatmap, with the following:
 1. The first data row (index=0) is renamed to "Total Residues" and displayed
    with a black background, but still shows numeric annotations.
 2. The second row (index=1) remains unchanged (e.g., "Pocket14, 11, 0, 100...").
 3. The rest of the rows and columns are color-mapped normally.
 4. The second column is not blacked out.
 5. The x-axis (column labels) is placed at the top, oriented vertically.

Usage:
    python csv_heatmap.py <input_csv> <config_json> <output_pdf>

Example:
    python csv_heatmap.py ADSS2_simple.csv config.json ADSS2_simple_heatmap.pdf

Explanation:
    This script uses two DataFrames:
      - df for the heatmap color data,
      - annot_df for the numeric text annotations.
    By converting row=0's data in df to NaN, that row is painted black.
    Meanwhile, row=0 in annot_df still contains the numbers.
    The second row (index=1) is left as-is, so if it says
    "Pocket14, 11, 0, 100,...", it will appear normally.
"""

import sys
import json
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np


def usage():
    print(__doc__)
    sys.exit(1)


def morpheus_json_to_colormap(json_scheme):
    fractions = json_scheme['fractions']
    colors = json_scheme['colors']

    # Normalize if needed
    fractions = np.array(fractions)
    if fractions.max() > 1:
        fractions /= fractions.max()

    cmap = mcolors.LinearSegmentedColormap.from_list('morpheus', list(zip(fractions, colors)))
    cmap.set_bad(color='black')
    return cmap


def main():
    if len(sys.argv) != 4:
        usage()

    csv_file = sys.argv[1]
    json_file = sys.argv[2]
    output_file = sys.argv[3]

    print(f"Running: {os.path.basename(__file__)} {csv_file} {json_file} {output_file}")

    # Load CSV
    df = pd.read_csv(csv_file, index_col=0, na_values=["-", "NA", ""])
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0)

    # Copy for annotations
    annot_df = df.copy()

    # If there's at least one row, rename row=0 to "Total Residues" & set df row=0 -> NaN
    if df.shape[0] > 0:
        row0_label = df.index[0]
        # rename that label in both
        df.rename(index={row0_label: "Total Residues"}, inplace=True)
        annot_df.rename(index={row0_label: "Total Residues"}, inplace=True)

        # black out row=0 in df so it appears black, but keep numeric in annot
        row0_idx = df.index.get_loc("Total Residues")
        df.iloc[row0_idx, :] = np.nan

    # Do not touch row=1 => it remains as is, e.g. "Pocket14"

    # Load config
    with open(json_file, 'r') as jf:
        config = json.load(jf)

    title = config.get('title', 'Heatmap')
    xlabel = config.get('xlabel', '')
    ylabel = config.get('ylabel', '')
    figsize = tuple(config.get('figsize', [10, 8]))

    # Prepare colormap
    if 'valueToColorScheme' in config:
        # Morpheus style
        scheme = next(iter(config['valueToColorScheme'].values()))
        cmap = morpheus_json_to_colormap(scheme)
    else:
        cmap_name = config.get('colormap', 'viridis')
        base = plt.get_cmap(cmap_name)
        cmap = mcolors.ListedColormap(base(np.linspace(0, 1, 256)))
        cmap.set_bad(color='black')

    plt.figure(figsize=figsize)

    ax = sns.heatmap(
        df,
        cmap=cmap,
        annot=annot_df,
        fmt="g",
        cbar=True,
        linewidths=0.5,
        linecolor='black'
    )

    # X-axis on top, vertical labels
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    plt.setp(ax.get_xticklabels(), rotation=90)

    plt.title(title, pad=20)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.tight_layout()
    plt.savefig(output_file)
    print(f"Heatmap saved to '{output_file}'.")

if __name__ == "__main__":
    main()
