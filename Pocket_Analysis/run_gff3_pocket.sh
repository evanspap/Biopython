# Check if the script path is provided
if [[ -z "$1" ]]; then
  echo "Usage: bash run_gff3_pocket.sh /path/to/dynamic_gff3_pocket.py"
  exit 1
fi

# Get the script path from the argument
SCRIPT_PATH="$1"

# Check if the Python script exists
if [[ ! -f "$SCRIPT_PATH" ]]; then
  echo "Error: Python script '$SCRIPT_PATH' not found!"
  exit 1
fi

# Iterate through all *_summary.csv files in subdirectories
for summary_csv in */*_summary.csv; do
  # Ensure the file exists
  if [[ -f "$summary_csv" ]]; then
    echo "Processing $summary_csv..."
    
    # Run the Python script
    python "$SCRIPT_PATH" "$summary_csv" -t 3
  else
    echo "No *_summary.csv files found in subdirectories."
    exit 1
  fi
done

echo "All files processed."
