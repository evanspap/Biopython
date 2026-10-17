#!/usr/bin/env bash
# Usage: bash batch_three_method_overlap_order.sh ROOT NEW_OUTPUT [extra Python CLI options]
# ROOT contains Fpockt_ASpock, ASpock_DGsite and Fpockt_DGsite.
# Override the interpreter with PYTHON=/path/to/python (one executable, no arguments).
set -euo pipefail
if [[ "$#" -lt 2 ]]; then
  echo "Usage: $0 <Fpockt-root> <new-output-dir> [--drop-globally-zero ...]" >&2
  exit 2
fi
script_dir="$(cd -- "$(dirname -- "$0")" && pwd)"
python_bin="${PYTHON:-python3}"
input_root="$1"
output_root="$2"
shift 2
exec "$python_bin" "$script_dir/three_method_overlap_order.py" \
  --batch-root "$input_root" --out-dir "$output_root" "$@"
