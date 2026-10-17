#!/usr/bin/env bash
# All proteins: common ordering, global <50% filter, Morpheus colors, A-D panels.
# Uses the existing mmseg-overlap-v2 Python environment. No installation/activation.
# Usage: bash run_all_pocket_overlap_morpheus_v2.sh
# Optional overrides: PYTHON, EVAN_ROOT, REFERENCE_DIR, RUN_TAG, SELECTED_TSV, PROTEIN.
set -euo pipefail

if [[ "$#" -ne 0 ]]; then
  echo "Usage: bash $0 (all proteins; no positional arguments)" >&2
  exit 2
fi

script_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
order_script="$script_dir/three_method_overlap_order.py"
render_script="$script_dir/../Morpheus/csv_heatmap_ordered.py"
generate_script="$script_dir/genbank_pairwise_overlap_complete_v2.py"
python_bin="${PYTHON:-$HOME/.venvs/mmseg-overlap-v2/bin/python}"
evan="${EVAN_ROOT:-/mnt/c/Partners HealthCare Dropbox/Evangelos Papadopoulos/MMSEGs Discussion/GeneDeepDives/Evangelos}"
references="${REFERENCE_DIR:-/mnt/c/Partners HealthCare Dropbox/Evangelos Papadopoulos/LAB-Evan/Manuscript MMSEG pockets/Figures/Comparison_Plots_all}"
pair_root="$evan/Genebank_Annotations/merge_AFB_Fpock_DGSite_ASpock_45nm/Fpockt"
run_tag="${RUN_TAG:-$(date +%Y%m%d_%H%M%S)}"
selected="${SELECTED_TSV:-/mnt/c/Partners HealthCare Dropbox/Evangelos Papadopoulos/LAB-Evan/Manuscript MMSEG pockets/Excel/combined_triple_pocket_analysis_20260720_v2_selected.tsv}"

if [[ ! "$run_tag" =~ ^[A-Za-z0-9_-]+$ ]]; then
  echo "RUN_TAG must contain only letters, digits, underscores or hyphens." >&2
  exit 2
fi
scope="${PROTEIN:-all}"
protein_args=()
if [[ -n "${PROTEIN:-}" ]]; then
  if [[ ! "$PROTEIN" =~ ^[A-Za-z0-9_-]+$ ]]; then
    echo "Invalid PROTEIN identifier." >&2
    exit 2
  fi
  protein_args=(--protein "$PROTEIN")
fi
raw="$evan/Pocket_Overlap_Complete_v2_WSL_${scope}_$run_tag"
ordered="$evan/Pocket_Overlap_Ordered_v2_WSL_${scope}_fullresidues_min50_diagonal_$run_tag"
heatmaps="$evan/Figures/Pocket_Overlap_Morpheus_v2_WSL_${scope}_fullresidues_min50_diagonal_$run_tag"

if [[ ! -x "$python_bin" ]]; then
  echo "Existing Python environment not found: $python_bin" >&2
  echo "Set PYTHON to your existing Linux Python executable." >&2
  exit 2
fi
for path in "$generate_script" "$order_script" "$render_script" "$selected"; do
  if [[ ! -r "$path" ]]; then
    echo "Cannot read script: $path. Check the Google Drive mount." >&2
    exit 2
  fi
done
for path in "$pair_root" "$references"; do
  if [[ ! -d "$path" ]]; then
    echo "Input directory not found: $path. Check the C: mount." >&2
    exit 2
  fi
done
if [[ -e "$raw" || -e "$ordered" || -e "$heatmaps" ]]; then
  echo "An output directory already exists; choose a different RUN_TAG." >&2
  exit 2
fi

printf 'Python: %s\nOrdered matrices/QC: %s\nHeatmaps/composites: %s\n' \
  "$python_bin" "$ordered" "$heatmaps"

printf 'Complete-residue matrices: %s\n' "$raw"
if ! "$python_bin" "$generate_script" \
  --genbank-root "$pair_root" \
  --out-dir "$raw" \
  --selected-tsv "$selected" "${protein_args[@]}"; then
  echo "Complete-residue generation/selected validation failed; ordering has not started." >&2
  echo "Inspect: $raw/batch_complete_pairwise_manifest.tsv" >&2
  exit 1
fi

if ! "$python_bin" "$order_script" \
  --batch-root "$raw" \
  --counts-root "$raw" \
  --percentage-policy exact \
  --out-dir "$ordered" \
  --drop-globally-below 50 \
  "${protein_args[@]}"; then
  echo "Ordering failed; rendering has not started." >&2
  echo "Inspect: $ordered/batch_overlap_order_manifest.tsv" >&2
  exit 1
fi

if ! "$python_bin" "$render_script" \
  --batch-dir "$ordered" \
  --out-dir "$heatmaps" \
  --mode both \
  --formats pdf svg png \
  --composite \
  --structure-reference-dir "$references" "${protein_args[@]}"; then
  echo "Rendering failed; inspect the terminal output and any generated manifest:" >&2
  echo "$heatmaps/batch_heatmap_manifest.tsv" >&2
  exit 1
fi

printf 'Completed requested proteins.\nFigures: %s\n' "$heatmaps"
