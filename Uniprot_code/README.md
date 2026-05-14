# Add UniProt Entry And Gene Names

This folder contains `add_uniprot_entry_gene_names.py`, a Python script that reads a CSV or TSV file, takes UniProt accessions from a chosen column, queries the UniProt REST API, and writes the enriched table to standard output with two added columns:

- `UniProt_entry_name`
- `Gene_name`

## Recommended Runtime

Ubuntu/WSL is recommended for running this script on Windows. It gives cleaner command-line behavior for standard output redirection and avoids Windows PowerShell 5.1 encoding issues.

Use `python3` from Ubuntu/WSL when possible.

The default input is:

```text
C:\Users\geras\Partners HealthCare Dropbox\Evangelos Papadopoulos\LAB-Evan\Manuscript MMSEG pockets\AF2BIND\analysis_files\analysis_files\af2bind_p2rank_human_proteome_predictions_revision_with_sitemap.csv
```

In that file, the UniProt accession column is the first column and is named `uniprot`. The script matches the requested column name case-insensitively, so the default `Uniprot` argument works with `uniprot`.

## Usage

Run with defaults in WSL/Ubuntu:

```bash
cd "/mnt/h/My Drive/VSCode_Github/BioPython/Uniprot_code"
python3 add_uniprot_entry_gene_names.py \
  --input "/mnt/c/Users/geras/Partners HealthCare Dropbox/Evangelos Papadopoulos/LAB-Evan/Manuscript MMSEG pockets/AF2BIND/analysis_files/analysis_files/af2bind_p2rank_human_proteome_predictions_revision_with_sitemap.csv" \
  > output_with_uniprot_entry_gene_names.csv
```

This assumes the Windows `H:` drive is mounted in WSL as `/mnt/h`. If `/mnt/h` is not available in your WSL setup, copy or clone this script into a Linux-visible folder, for example under your WSL home directory, and run it from there. The default input file is still read from `/mnt/c/...`, which is usually available in WSL.

Run with explicit input in WSL/Ubuntu:

```bash
python3 add_uniprot_entry_gene_names.py \
  --input "/mnt/c/Users/geras/Partners HealthCare Dropbox/Evangelos Papadopoulos/LAB-Evan/Manuscript MMSEG pockets/AF2BIND/analysis_files/analysis_files/af2bind_p2rank_human_proteome_predictions_revision_with_sitemap.csv" \
  > output_with_uniprot_entry_gene_names.csv
```

Run with defaults in Windows PowerShell:

```powershell
python .\add_uniprot_entry_gene_names.py `
  --input "C:\Users\geras\Partners HealthCare Dropbox\Evangelos Papadopoulos\LAB-Evan\Manuscript MMSEG pockets\AF2BIND\analysis_files\analysis_files\af2bind_p2rank_human_proteome_predictions_revision_with_sitemap.csv" `
  > output_with_uniprot_entry_gene_names.csv
```

Run with explicit input in Windows PowerShell and redirect the output to a file:

```powershell
python .\add_uniprot_entry_gene_names.py `
  --input "C:\path\to\input.csv" `
  > "C:\path\to\output.csv"
```

Use a different UniProt column:

```powershell
python .\add_uniprot_entry_gene_names.py --input "C:\path\to\input.csv" --column Uniprot
```

Running the script without arguments prints the help message and exits:

```bash
python3 add_uniprot_entry_gene_names.py
```

For TSV files, the delimiter is inferred from the `.tsv` extension. For other delimiters, pass `--delimiter`.

## Output

The enriched CSV/TSV is written to standard output. Use shell redirection (`>`) to save it to a file.

Progress messages and warnings are written to standard error, so they do not become part of the redirected CSV/TSV output.

The original rows and columns are preserved. The two new columns are appended at the end of the table.

WSL/Ubuntu is recommended for the simplest CSV redirect behavior because `>` writes normal UTF-8 text output.

On Windows PowerShell 5.1, `>` may save redirected text as UTF-16. If a downstream tool expects UTF-8, use WSL/Ubuntu, PowerShell 7, or write with explicit encoding:

```powershell
python .\add_uniprot_entry_gene_names.py | Set-Content -Encoding utf8 "output_with_uniprot_entry_gene_names.csv"
```

## Dependencies And Requirements

- Python 3.9 or newer is recommended.
- WSL/Ubuntu with `python3` is recommended on Windows.
- No third-party Python packages are required.
- The script uses only the Python standard library: `argparse`, `csv`, `pathlib`, `sys`, `time`, and `urllib`.
- Internet access is required.
- The input file must have a header row.
- The UniProt accession column defaults to `Uniprot` and is matched case-insensitively.

The script queries the UniProt REST API:

```text
https://rest.uniprot.org/uniprotkb/search
```
