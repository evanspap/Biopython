#!/usr/bin/env python3

"""
CSV Heatmap Generator (Simplified)
==================================
This script reads a CSV file containing data for a heatmap, then displays
all rows and columns with standard color coding and numeric annotations.
No rows or columns are renamed or blacked out.

Usage:
    python csv_heatmap_simple.py <input_csv> <config_json> <output_pdf>

Example:
    python csv_heatmap_simple.py ADSS2_simple.csv config.json ADSS2_simple_heatmap.pdf

Explanation:
    1. Load the CSV, converting all entries to numeric.
    2. Create a heatmap with seaborn, using numeric annotations.
    3. Place the X-axis (column labels) on top, oriented vertically.
    4. Save the resulting heatmap as a PDF (or other supported format).
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

    cmap = mcolors.LinearSegmentedColormap.from_list(
        'morpheus',
        list(zip(fractions, colors))
    )
    cmap.set_bad(color='black')
    return cmap

def main():
    if len(sys.argv) != 4:
        usage()

    csv_file = sys.argv[1]
    json_file = sys.argv[2]
    output_file = sys.argv[3]

    print(f"Running: {os.path.basename(__file__)} {csv_file} {json_file} {output_file}")

    # 1. Load CSV, convert to numeric
    df = pd.read_csv(csv_file, index_col=0, na_values=["-", "NA", ""])
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0)

    # We'll annotate with the same data
    annot_df = df.copy()

    # 2. Load config
    with open(json_file, 'r') as jf:
        config = json.load(jf)

    title = config.get('title', 'Heatmap')
    xlabel = config.get('xlabel', '')
    ylabel = config.get('ylabel', '')
    figsize = tuple(config.get('figsize', [10, 8]))

    # 3. Prepare colormap
    if 'valueToColorScheme' in config:
        # Morpheus style
        scheme = next(iter(config['valueToColorScheme'].values()))
        cmap = morpheus_json_to_colormap(scheme)
    else:
        cmap_name = config.get('colormap', 'viridis')
        base = plt.get_cmap(cmap_name)
        cmap = mcolors.ListedColormap(base(np.linspace(0, 1, 256)))
        cmap.set_bad(color='black')  # NaN => black

    plt.figure(figsize=figsize)

    # 4. Create standard heatmap with numeric annotations
    ax = sns.heatmap(
        df,
        cmap=cmap,
        annot=annot_df,
        fmt="g",
        cbar=True,
        linewidths=0.5,
        linecolor='black'
    )

    # Place X-axis on top, orient vertically
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
