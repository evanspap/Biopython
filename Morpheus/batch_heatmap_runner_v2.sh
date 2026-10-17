#!/usr/bin/env bash
# Usage: bash batch_heatmap_runner_v2.sh ORDERED_ROOT NEW_HEATMAP_ROOT [extra CLI options]
# Defaults: percentage plots, PDF + SVG. Add --mode both --formats pdf svg png.
set -euo pipefail
if [[ "$#" -lt 2 ]]; then
  echo "Usage: $0 <ordered-root> <new-heatmap-dir> [--mode both --formats pdf svg png ...]" >&2
  exit 2
fi
script_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
python_bin="${PYTHON:-python3}"
input_root="$1"
output_root="$2"
shift 2
exec "$python_bin" "$script_dir/csv_heatmap_ordered.py" \
  --batch-dir "$input_root" --out-dir "$output_root" "$@"
