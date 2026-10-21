#!/usr/bin/env bash
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
python_bin="${PYTHON:-$HOME/.venvs/mmseg-overlap-v2/bin/python}"
evan='/mnt/c/Partners HealthCare Dropbox/Evangelos Papadopoulos/MMSEGs Discussion/GeneDeepDives/Evangelos'
raw="$evan/Pocket_Overlap_Complete_v2_WSL_all_20261013_141015"
tag="${RUN_TAG:-$(date +%Y%m%d_%H%M%S)}"
[[ "$tag" =~ ^[A-Za-z0-9_-]+$ ]] || exit 2
ordered="$evan/Pocket_Overlap_Ordered_v3_ZNF592_$tag"
figures="$evan/Figures/Pocket_Overlap_Morpheus_v3_ZNF592_$tag"
"$python_bin" "$script_dir/three_method_overlap_order_v3.py" \
  --batch-root "$raw" --counts-root "$raw" --percentage-policy exact \
  --drop-globally-below 50 --protein ZNF592 --out-dir "$ordered"
"$python_bin" "$script_dir/../Morpheus/csv_heatmap_ordered.py" \
  --batch-dir "$ordered" --out-dir "$figures" --protein ZNF592 \
  --mode both --formats pdf svg png --composite \
  --structure-reference-dir '/mnt/c/Partners HealthCare Dropbox/Evangelos Papadopoulos/LAB-Evan/Manuscript MMSEG pockets/Figures/Comparison_Plots_all'
printf 'Completed ZNF592 only. Figures: %s\n' "$figures"