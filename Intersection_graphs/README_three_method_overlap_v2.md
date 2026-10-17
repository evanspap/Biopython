# Coordinated three-method pocket overlaps, v2

Current final workflow: use the complete-GenBank Bash runner described in the FINAL workflow section below. Direct examples using legacy *_relative.csv inputs are retained for reference and should not be used for final figures.

Scripts prepared on 2026-10-11. No dependency installation, pipeline, syntax check,
test, or rendering command was run when these scripts were created.
Legacy scripts, CSVs, figures, and the manuscript are unchanged.

Python 3.10 or newer is required. Ordering uses only the standard library.
The renderer needs numpy and matplotlib. Windows paths, including long Dropbox
filenames, are supported. The four scripts are new files beside their canonical
legacy predecessors in the Google Drive repository.

## Windows PowerShell: complete batch

Run these commands yourself. Google Drive stays G:; Dropbox uses the new location.
If python resolves to the wrong interpreter, use your environment's python.exe
or replace python with py -3. Python batch mode requires no Bash installation.
The dependency-installation command can be omitted if these packages are installed.

~~~powershell
$repo = 'G:\My Drive\VSCode_Github'
$evan = 'C:\Partners HealthCare Dropbox\Evangelos Papadopoulos\MMSEGs Discussion\GeneDeepDives\Evangelos'
$pairRoot = Join-Path $evan 'Genebank_Annotations\merge_AFB_Fpock_DGSite_ASpock_45nm\Fpockt'
$ordered = Join-Path $evan 'Pocket_Overlap_Ordered_v2_20261011'
$heatmaps = Join-Path $evan 'Figures\Pocket_Overlap_Ordered_v2_20261011'

python -m pip install -r "$repo\BioPython\Intersection_graphs\requirements_overlap_v2.txt"

python "$repo\BioPython\Intersection_graphs\three_method_overlap_order.py" `
  --batch-root "$pairRoot" `
  --out-dir "$ordered" `
  --drop-globally-zero
if ($LASTEXITCODE -ne 0) { throw 'Ordering failed; inspect the error and batch manifest before rendering.' }

python "$repo\BioPython\Morpheus\csv_heatmap_ordered.py" `
  --batch-dir "$ordered" `
  --out-dir "$heatmaps" `
  --mode both `
  --formats pdf svg png
if ($LASTEXITCODE -ne 0) { throw 'Rendering failed; inspect batch_heatmap_manifest.tsv.' }
~~~

Use --formats pdf svg to omit PNG; --mode percentage to omit count plots.
Omit --drop-globally-zero to retain all pockets. Output roots must be NEW.
The scripts refuse to overwrite any earlier run. Choose new names for retries.

## Optional ADSS2 check before the full batch

Use the same $repo, $evan and $pairRoot definitions above. Separate output roots
avoid any conflict with the complete-batch example.

~~~powershell
$orderedAdss2 = Join-Path $evan 'Pocket_Overlap_Ordered_v2_ADSS2_20261011'
$heatmapsAdss2 = Join-Path $evan 'Figures\Pocket_Overlap_Ordered_v2_ADSS2_20261011'

python "$repo\BioPython\Intersection_graphs\three_method_overlap_order.py" `
  --batch-root "$pairRoot" --protein ADSS2 `
  --out-dir "$orderedAdss2" --drop-globally-zero
if ($LASTEXITCODE -ne 0) { throw 'ADSS2 ordering failed.' }

python "$repo\BioPython\Morpheus\csv_heatmap_ordered.py" `
  --batch-dir "$orderedAdss2" --out-dir "$heatmapsAdss2" `
  --mode both --formats pdf svg png
~~~

The earlier context proposed dropping ASpock11, ASpock15, Fpockt02, Fpockt04,
Fpockt10 and no DoGSite pockets for ADSS2. This list was not reverified by running
v2. Inspect the removed records in ADSS2_overlap_order_manifest.tsv after your run.

## Direct single-protein ordering

~~~powershell
python "$repo\BioPython\Intersection_graphs\three_method_overlap_order.py" `
  --protein ADSS2 `
  --as-fp "$pairRoot\Fpockt_ASpock\Intersection_Percentage\ADSS2_relative.csv" `
  --as-dg "$pairRoot\ASpock_DGsite\Intersection_Percentage\ADSS2_relative.csv" `
  --fp-dg "$pairRoot\Fpockt_DGsite\Intersection_Percentage\ADSS2_relative.csv" `
  --out-dir "$evan\Pocket_Overlap_Ordered_v2_ADSS2_direct_20261011" `
  --drop-globally-zero
~~~

Batch defaults to *_relative.csv in each of Fpockt_ASpock, ASpock_DGsite,
Fpockt_DGsite, under Intersection_Percentage. Add --input-kind simplified to select
*_relative_simple.csv in their Simplified directories instead. Single-protein
mode recognizes either layout, including transposed matrices.

## Bash / WSL: mounts, Python environment, and execution

Commands below are supplied for you to run; none was executed during preparation.
Keep Google Drive for desktop and Dropbox running in Windows. The expected paths
are G:\My Drive\VSCode_Github and the new C:\Partners HealthCare Dropbox location.

First open a normal Windows PowerShell terminal and enter WSL:

~~~powershell
wsl --list --verbose
wsl
~~~

This opens your default distribution. For a specific distribution, use
wsl -d NAME, with NAME exactly as shown in the list. The package-installation
commands below assume Ubuntu or Debian. Python 3.10 or newer is required.

All remaining commands run INSIDE WSL, not in PowerShell.
Mount C: and G: only if their expected mountpoints are not already mounted:

~~~bash
sudo mkdir -p /mnt/c /mnt/g
mountpoint -q /mnt/c || sudo mount -t drvfs C: /mnt/c
mountpoint -q /mnt/g || sudo mount -t drvfs G: /mnt/g

repo='/mnt/g/My Drive/VSCode_Github'
evan='/mnt/c/Partners HealthCare Dropbox/Evangelos Papadopoulos/MMSEGs Discussion/GeneDeepDives/Evangelos'
pair_root="$evan/Genebank_Annotations/merge_AFB_Fpock_DGSite_ASpock_45nm/Fpockt"

ls -l "$repo/BioPython/Intersection_graphs/three_method_overlap_order.py"
ls -l "$repo/BioPython/Morpheus/csv_heatmap_ordered.py"

for comparison in Fpockt_ASpock ASpock_DGsite Fpockt_DGsite; do
  ls -l "$pair_root/$comparison/Intersection_Percentage/ADSS2_relative.csv"
done
~~~

Proceed after these checks show the expected files. DrvFs is the Windows-drive
mount mechanism documented by [Microsoft](https://learn.microsoft.com/en-us/windows/wsl/wsl2-mount-disk#differences-between-mounting-an-external-drive-with-windows-formatting-versus-linux-formatting).
If mounting G: fails, first confirm that G:\My Drive is accessible in Windows
and Google Drive for desktop is running. Existing mounts are left in place.
These manual mount commands can be repeated in a later WSL session when needed.

Create a Linux Python environment in your WSL home directory:

~~~bash
sudo apt update
sudo apt install -y python3 python3-venv python3-pip

python3 --version
mkdir -p "$HOME/.venvs"
python3 -m venv "$HOME/.venvs/mmseg-overlap-v2"
source "$HOME/.venvs/mmseg-overlap-v2/bin/activate"

python -m pip install -r "$repo/BioPython/Intersection_graphs/requirements_overlap_v2.txt"
export PYTHON="$VIRTUAL_ENV/bin/python"
~~~

On later sessions, only activate the environment and set PYTHON again; installation
is needed once. The [Python venv documentation](https://docs.python.org/3/library/venv.html)
describes creating and activating a Linux virtual environment.

Optional ADSS2 run first, using new output directories:

~~~bash
ordered_adss2="$evan/Pocket_Overlap_Ordered_v2_WSL_ADSS2_20261011"
heatmaps_adss2="$evan/Figures/Pocket_Overlap_Ordered_v2_WSL_ADSS2_20261011"

bash "$repo/BioPython/Intersection_graphs/batch_three_method_overlap_order.sh" \
  "$pair_root" "$ordered_adss2" --protein ADSS2 --drop-globally-zero &&
bash "$repo/BioPython/Morpheus/batch_heatmap_runner_v2.sh" \
  "$ordered_adss2" "$heatmaps_adss2" --mode both --formats pdf svg png
~~~

Complete batch, using different new output directories:

~~~bash
ordered="$evan/Pocket_Overlap_Ordered_v2_WSL_20261011"
heatmaps="$evan/Figures/Pocket_Overlap_Ordered_v2_WSL_20261011"

bash "$repo/BioPython/Intersection_graphs/batch_three_method_overlap_order.sh" \
  "$pair_root" "$ordered" --drop-globally-zero &&
bash "$repo/BioPython/Morpheus/batch_heatmap_runner_v2.sh" \
  "$ordered" "$heatmaps" --mode both --formats pdf svg png
~~~

The && prevents rendering if ordering returns an error. Inspect
batch_overlap_order_manifest.tsv in that case. Existing output directories are
never overwritten; choose new output names for another run. Use --formats pdf svg
to omit PNG, or --mode percentage to omit count plots. Both .sh files have LF line
endings; invoking them through bash requires no chmod.

To view the full-batch heatmaps in Windows Explorer from WSL:

~~~bash
explorer.exe "$(wslpath -w "$heatmaps")"
~~~

These paths are for WSL. Git Bash uses /g and /c instead of /mnt/g and /mnt/c.


## Exact counts and legacy percentages

The legacy generator stores int((shared / denominator_pocket_total) * 100).
The simplifier takes the larger of the two directional percentages. v2 preserves
that definition and the input percentages.

Legacy matrices contain pocket totals. v2 enumerates compatible integer counts in
a bounded interval and verifies the original floating-point formula in both
directions. It accepts a count only if EXACTLY ONE integer is compatible, recording
uniquely_inferred_from_percentage_and_totals as its source. For example, a
max-directional 100% with minimum pocket total 25 implies 25 shared residues; with
minimum total 3 it implies 3. This is exact recovery under the legacy definition,
not a rounded percentage-times-size estimate.

With larger totals, truncation may allow several counts. Missing totals or
inconsistent/ambiguous values produce an error for that protein and nonzero batch
exit. Supply count matrices using --as-fp-counts, --as-dg-counts, --fp-dg-counts
in single-protein mode. Count matrices must contain the same pocket sets as their
percentage matrices; orientation/order may differ. Blank/NA is unknown, never zero.
Supplied counts are validated against percentages and totals when available.

Batch --counts-root ROOT selects:
ROOT/Fpockt_ASpock/PROTEIN_counts.csv,
ROOT/ASpock_DGsite/PROTEIN_counts.csv,
ROOT/Fpockt_DGsite/PROTEIN_counts.csv.
All these files are required for each selected protein when this option is used.

For full-precision max-directional percentages, use --percentage-policy exact.
The legacy-floor default rejects noninteger percentages.

The existing GB_Intersection_Percentage.py reads only the first line containing
join(...); general/multiline GenBank feature parsing is not robust. Recovered
counts describe the legacy matrices' PARSED residue sets; they do not validate
complete biological feature locations. Repairing that upstream parser and
recomputing overlaps remains a separate workflow, as requested.

## Ordering and zero policy

Canonical outputs: AutoSite rows / FPocket columns; AutoSite rows / DoGSite columns;
FPocket rows / DoGSite columns. Each method's sequence is reused in both matrices.

1. Seed edge: descending percentage, shared count, summed percentage support of
   both endpoints to the third method, then ascending identifiers. Third pocket:
   maximize summed triad percentages, counts, total incident support, then IDs.
2. Remaining triads: descending sum of three percentages, sum of three counts,
   total incident percentage support, then ascending IDs. Incident support uses
   input pockets. When a method is exhausted, apply the same rule to pairs/singles.
3. Up to --swap-passes improving swaps (default 30), exchanging two pockets within
   one method. The seed triad stays locked at position 1. Each accepted swap
   improves the lexicographic global objective: (summed diagonal percentages,
   summed diagonal shared counts). Unequal method lengths keep contiguous axes.

This deterministic greedy/local heuristic does not guarantee a global optimum or
improvement over the original baseline. QC records the baseline comparison and
prints a warning if worse. The greedy score is recorded separately so that swap
improvement can be verified. This diagonal score is not the biological three-way
residue-intersection score used in the separate consensus workflow.
Retained identifiers are checked as exact permutations before output.

--drop-globally-zero requires EVERY cell to BOTH other methods to have a known zero
percentage AND known zero count. Unknown cells keep their pocket. Rounded 0% with
positive supplied counts also keeps its pocket. Filtering uses the full graph
before any removal.

Missing contributes no percentage/count to the objective and remains blank in CSV.
Absent pockets at the tail have empty method fields in QC; no fake pockets or
interior gaps are inserted into axes.

## Outputs and QC

Each protein directory contains six ordered CSVs (percentage/count for each pair),
three pocket_order.tsv files (position, identifier, total residues),
PROTEIN_overlap_order_manifest.tsv and PROTEIN_overlap_order_qc.json.

The manifest has summary, diagonal and removed records; all three percentages,
counts and count sources; per-position scores; before/after diagonal scores.
QC includes version, input paths/SHA-256 hashes, baseline/greedy/final scores,
method orders, seed, swap count and removal lists.

The baseline is the common input order after filtering, using each method's first
appearance in canonical pair order. It is not the unrelated independent ordering
of each old heatmap.

Batch preflight matches exact protein filename stems and checks all three
collections before processing. Missing comparisons abort. Per-protein data errors
are logged in batch_overlap_order_manifest.tsv; processing continues for other
proteins and the batch exits nonzero if any failed. Only successfully completed
proteins get a complete QC JSON. Partial outputs after I/O failure are preserved;
choose a new output directory when retrying.

Renderer batch mode reads completed v2 QC and listed output CSVs only. It writes
batch_heatmap_manifest.tsv and one render_qc.json per heatmap. Percentage range
is fixed to 0..100; counts have a separate data-dependent scale or --count-vmax.
A too-small count maximum raises an error instead of clipping. PDF/SVG cells and
text stay vector; optional PNG defaults to 300 dpi. Figure size scales with pockets.
Grey hatching/NA means missing; numeric 0 and the colormap endpoint mean true zero.
--no-annotate omits cell text; hatching remains.

--color-config PATH accepts a Morpheus JSON valueToColorScheme (fractions/colors).
Default colors use the exact legacy black/blue/white/red Morpheus scheme. No rows/columns are sorted by the renderer.
Totals are metadata and are excluded from plot matrices. Empty retained axes
render a clear no-retained-pockets panel. Four-panel composition is now supported with --composite; see the ADSS2 commands below. The manuscript remains unchanged.

## Optional verification commands for you

These commands were supplied only, not executed. The compile check writes Python
cache files but does not run the pipeline.

~~~powershell
python -m py_compile `
  "$repo\BioPython\Intersection_graphs\three_method_overlap_order.py" `
  "$repo\BioPython\Morpheus\csv_heatmap_ordered.py"
~~~

After the ADSS2 run, inspect QC, the manifest, and percentage/count plots. Verify
common orders in both matrices containing each method, removal records versus
the earlier expected list, and a nondecreasing greedy-to-final score.

## Morpheus colors and four-panel ADSS2 figure

The renderer now defaults to the exact legacy Morpheus configuration stored beside
the script in morpheus_pocket_overlap_legacy.json. It preserves the intentional
discontinuity at 50%: black below that boundary, blue above it, white at 75%, and
red at 100%. Duplicate color fractions are supported. Missing cells use grey
hatching and NA, distinct from black true-zero cells. --color-config still allows
an explicit alternate configuration.

--composite adds an A–D figure: structure at top left, AutoSite/FPocket at top
right, AutoSite/DoGSite at bottom left, and FPocket/DoGSite at bottom right.
Heatmaps are drawn directly from ordered CSVs, with matching method orders checked
before composition. PDF/SVG heatmaps and colorbars remain vector. Panel A reuses
the existing structure image and is raster; it is not a new molecular rendering.

Run the following INSIDE WSL to rerender only ADSS2, reusing its existing ordering.
No ordering command is needed. These commands were supplied, not executed here.

~~~bash
repo='/mnt/g/My Drive/VSCode_Github'
evan='/mnt/c/Partners HealthCare Dropbox/Evangelos Papadopoulos/MMSEGs Discussion/GeneDeepDives/Evangelos'
reference='/mnt/c/Partners HealthCare Dropbox/Evangelos Papadopoulos/LAB-Evan/Manuscript MMSEG pockets/Figures/Comparison_Plots_all/ADSS2_structure_and_pairwise_pocket_overlap_comparisons.png'
source "$HOME/.venvs/mmseg-overlap-v2/bin/activate"

python "$repo/BioPython/Morpheus/csv_heatmap_ordered.py" \
  --batch-dir "$evan/Pocket_Overlap_Ordered_v2_WSL_ADSS2_20261011" \
  --protein ADSS2 \
  --out-dir "$evan/Figures/Pocket_Overlap_Morpheus_v2_WSL_ADSS2_20261011" \
  --mode both --formats pdf svg png \
  --composite --structure-reference "$reference"
~~~

Outputs are in the new output root's ADSS2 subdirectory: six separate heatmaps
plus two four-panel composites (percentage and counts), each in PDF/SVG/PNG, and QC.
Composite filenames end in comparisons_ordered_percentage or
comparisons_ordered_counts. Existing figures and ordered matrices are retained.
Use a new output root on retries. --mode percentage produces only percentage
heatmaps and one percentage composite.

The default --structure-crop .055 .145 .475 .495 is a normalized rectangle matched
to the supplied ADSS2 legacy figure. It excludes the old panel title/border/letter.
For differently laid-out figures, supply another rectangle. A standalone structure
PNG can be supplied with --structure-image (no crop), and a folder of legacy
figures can be supplied with --structure-reference-dir for a full batch; filenames
must be PROTEIN_structure_and_pairwise_pocket_overlap_comparisons.png.

## Global 50% filter (ADSS2)

Use --drop-globally-below 50 in the ordering step to remove a pocket only when
every percentage to both other methods is known and strictly below 50%.
A pocket with any percentage >=50%, or any missing percentage, is retained.
This removes axes consistently across the three comparisons, not separately for
each heatmap. It is a separate option from --drop-globally-zero and preserves
that earlier zero-only option. All filtering is evaluated against the original
graph before removals. QC records removed_pockets, removal_reasons and the threshold.
The renderer no longer displays a Missing (NA) legend; real missing cells remain
hatched/annotated NA and are counted in QC.

Run inside your activated WSL environment. These commands were not run here.
New directories are required because this filter changes the retained pockets.

~~~bash
repo='/mnt/g/My Drive/VSCode_Github'
evan='/mnt/c/Partners HealthCare Dropbox/Evangelos Papadopoulos/MMSEGs Discussion/GeneDeepDives/Evangelos'
pair_root="$evan/Genebank_Annotations/merge_AFB_Fpock_DGSite_ASpock_45nm/Fpockt"
reference='/mnt/c/Partners HealthCare Dropbox/Evangelos Papadopoulos/LAB-Evan/Manuscript MMSEG pockets/Figures/Comparison_Plots_all/ADSS2_structure_and_pairwise_pocket_overlap_comparisons.png'
ordered="$evan/Pocket_Overlap_Ordered_v2_WSL_ADSS2_min50_20261011"
heatmaps="$evan/Figures/Pocket_Overlap_Morpheus_v2_WSL_ADSS2_min50_20261011"

source "$HOME/.venvs/mmseg-overlap-v2/bin/activate"

python "$repo/BioPython/Intersection_graphs/three_method_overlap_order.py" \
  --batch-root "$pair_root" --protein ADSS2 \
  --out-dir "$ordered" --drop-globally-below 50 &&
python "$repo/BioPython/Morpheus/csv_heatmap_ordered.py" \
  --batch-dir "$ordered" --protein ADSS2 --out-dir "$heatmaps" \
  --mode both --formats pdf svg png \
  --composite --structure-reference "$reference"
~~~


## One Bash script for the complete final workflow

Run inside WSL, using the already installed environment:

~~~bash
bash "/mnt/g/My Drive/VSCode_Github/BioPython/Intersection_graphs/run_all_pocket_overlap_morpheus_v2.sh"
~~~

The script processes all proteins with --drop-globally-below 50 and renders
percentage/count Morpheus heatmaps plus A-D composites in PDF/SVG/PNG.
It directly uses the existing $HOME/.venvs/mmseg-overlap-v2/bin/python, or PYTHON
if supplied. It does not install or activate an environment. The script itself
was not executed during preparation.

Output folders include a runtime timestamp, so earlier results are retained.
Their full paths are printed before processing and on completion.
A failed ordering stops execution before rendering. The script assumes the
C: and G: mounts are available. Optional overrides are EVAN_ROOT, REFERENCE_DIR,
PYTHON and RUN_TAG. RUN_TAG permits only letters, digits, underscores and hyphens.


## Empty-method input correction

Legacy square matrices containing only one method are valid when the other method
has no detected pockets. The parser now preserves empty method axes instead of
raising an identification error. QC records methods_without_input_pockets. No
zero-valued overlaps or artificial pockets are inserted. This addresses MAF,
POU2AF1, UBALD2 and XBP1 in the first full batch: 112 completed and 4 errored,
which prevented the runner from entering its rendering stage. Rerun the same Bash
runner after this correction; it creates new timestamped output directories.
This correction was inspected as source only; the pipeline was not rerun here.


## Selected TSV ordering priority (2026-10-12)

The all-protein runner now defaults to --selected-priority using
combined_triple_pocket_analysis_20260720_v2_selected.tsv from the manuscript Excel
folder. Override SELECTED_TSV to use another exported selected table. No Excel
reader is required; the TSV is the authoritative supplied row order.

For each protein, unique constituent pocket IDs are placed first in the order of
their first appearance in the TSV. These per-method prefixes are fixed during
swap optimization and reused across both relevant pairwise matrices. Unselected
pockets follow using the existing greedy/swap method. The 50% filter remains in
force: a selected row with a filtered constituent is skipped and recorded in QC.
A selected ID absent from the actual input is an error rather than a guessed ID.
Proteins absent from the selected TSV use the existing default ordering.

Repeated constituents make simultaneous diagonal placement of every selected
triad impossible without duplicating pockets. For ALG3, P1 and P2 share FPocket4
and DoGSite1, so the unique axes cannot place both entire triads at separate
diagonal positions. QC selected_records records original TSV lines, constituent
IDs, filtered constituents, final positions and aligned_on_diagonal. This affects
display order only: overlaps and counts are not increased or recalculated.

The TSV ranking is not a decreasing-overlap ranking: ALG3_P1 has max%=0.586,
ALG3_P2 has 0.75, and ALG3_P3 has 1. The first selected triad therefore need not
have 100% pairwise overlaps. A three-way directional max% of 1 does not imply
all three max-directional pairwise percentages are 100%.

Test ALG3 with the existing environment, without running all proteins:

~~~bash
PROTEIN=ALG3 bash "/mnt/g/My Drive/VSCode_Github/BioPython/Intersection_graphs/run_all_pocket_overlap_morpheus_v2.sh"
~~~

Then run all proteins:

~~~bash
bash "/mnt/g/My Drive/VSCode_Github/BioPython/Intersection_graphs/run_all_pocket_overlap_morpheus_v2.sh"
~~~

Output names include min50_selected and a fresh timestamp. The script prints exact
paths. Commands were supplied only; no pipeline or rendering was executed here.


## FINAL workflow: complete GenBank residue sets (2026-10-12)

The runner now executes these stages in sequence:
1. genbank_pairwise_overlap_complete_v2.py reads the existing merged GenBank files
   from Genebank_Annotations/merge_AFB_Fpock_DGSite_ASpock_45nm/Fpockt.
2. It writes new exact-count and full-precision percentage matrices, plus per-pocket
   complete residue lists and QC, into Pocket_Overlap_Complete_v2_WSL_SCOPE_TIMESTAMP.
3. Ordering reads ONLY these new matrices, with --counts-root and
   --percentage-policy exact. The selected TSV prefix order and global <50% filter
   are applied to the complete-residue comparisons.
4. Morpheus individual heatmaps and structure/comparison composites are rendered
   in PNG, SVG and PDF, for percentages and counts.

Original pocket-detection outputs are not needed or reread. Existing merged
GenBank annotations are the chosen source. No pocket detection is performed.
Legacy scripts, raw matrices, selected TSV, workbook and existing figures are
not overwritten. Older figures based on first-line parsing remain provisional.

Biopython parses wrapped join() expressions, all location parts, and entire
inclusive ranges. Every residue is converted to 1-based numbering and deduplicated.
Remote/fuzzy locations, duplicate pocket IDs, absent parsed locations, empty
pocket residue sets and out-of-sequence positions are rejected. Absent whole
methods are recorded explicitly; no zero-filled invented comparisons are made.
See [Biopython's location iteration documentation](https://biopython.org/docs/latest/api/Bio.SeqFeature.html).

Percentages are max(100*shared/left_size, 100*shared/right_size), without integer
truncation. Shared counts are the exact pairwise residue-set intersection.
The existing old GB_Intersection_Percentage.py truncates multiline annotations,
which explains ALG3 ASpock01's 13-residue legacy total versus 29 in the full feature.
The new per-pocket residue TSV should report 29 for that feature when you run it.

Before writing a protein's matrices, selected TSV rows are checked against the
complete GenBank: constituent IDs and sizes, exact triple shared count, full triple
shared residue list, three directional fractions and max%. Counts/lists must agree
exactly; fractions allow 0.00051 absolute tolerance for the TSV's displayed
three-decimal rounding. Score columns, including spreadsheet #VALUE! errors, are
not interpreted as overlap data. A mismatch fails that protein and the generation
batch returns a nonzero exit, stopping ordering and rendering.

The legacy triple builder itself takes only each location part's start; for
singleton joins this is complete, but genuine ranges may require upstream triple
recomputation. This generator expands ranges fully and will not silently reproduce
that limitation: selected validation catches incompatible outcomes.

No pipeline, dependency installation, tests, compilation or rendering was run
during preparation. The following commands are for you to execute inside WSL.
Use the existing Python environment; only Biopython is an added dependency.
Installation is needed once if it is not already installed.

~~~bash
"$HOME/.venvs/mmseg-overlap-v2/bin/python" -m pip install "biopython>=1.83"

PROTEIN=ALG3 bash "/mnt/g/My Drive/VSCode_Github/BioPython/Intersection_graphs/run_all_pocket_overlap_morpheus_v2.sh"
~~~

Inspect ALG3_P1 validation and ASpock01 total (29) in complete-matrix QC. It should
now pass the 50% threshold based on its full residue set, provided the supplied
selected TSV matches the source annotations. No inferred legacy counts are used.

Then execute all proteins:

~~~bash
bash "/mnt/g/My Drive/VSCode_Github/BioPython/Intersection_graphs/run_all_pocket_overlap_morpheus_v2.sh"
~~~

Exact output paths are printed. New ordered and figure directory names contain
fullresidues_min50_selected, distinguishing them from earlier incomplete matrices.
On generation failure, inspect batch_complete_pairwise_manifest.tsv in the printed
complete-matrix directory. Per-protein QC is under QC/PROTEIN and includes source
and selected TSV hashes and all matched selected records. No silent fallback to
legacy matrices is performed.


## Final figure layout and diagonal ordering — 2026-10-13

The runner uses complete GenBank residue matrices and the common greedy/swap
ordering again. It no longer supplies --selected-priority: selected TSV rows
are validated against complete GenBank data in the generation stage, but their
ranking does not pin the heatmap axes. The global <50% filter remains based on
exact unrounded percentages. Legacy matrices are not reused.

New composite layout:
A (top left): structure.
B (top right): AutoSite–DoGSite.
C (bottom left): FPocket–DoGSite.
D (bottom right): AutoSite–FPocket.

Percentage cell text rounds UP to an integer (ceil): 58.62 -> 59, 93.10 -> 94.
Near-integer floating-point noise is stabilized to 10 decimal places first.
Count cell text remains exact integer counts. CSVs, colors, filtering, ordering
and diagonal objectives retain full precision: a displayed 100 can therefore
represent an actual percentage between 99 and 100. Rounding is for labels only.
The missing legend stays removed. Rendering QC records rounding and panel layout.

Inside WSL, using the existing environment:
~~~bash
PROTEIN=ALG3 bash "/mnt/g/My Drive/VSCode_Github/BioPython/Intersection_graphs/run_all_pocket_overlap_morpheus_v2.sh"
~~~

For all proteins:
~~~bash
bash "/mnt/g/My Drive/VSCode_Github/BioPython/Intersection_graphs/run_all_pocket_overlap_morpheus_v2.sh"
~~~

Outputs are new timestamped directories containing fullresidues_min50_diagonal.
All percentage/count heatmaps and both composites are generated as PNG/SVG/PDF.
Commands were supplied only; no pipeline, tests or rendering were executed here.


## Resume interrupted rendering (2026-10-13)

--resume reuses an existing rendering output directory, skipping completed figures
only when their QC matches input hashes, matrix axes, mode, scale, rounding, color
config path and composite layout/source/crop. All requested output files must exist
and be nonempty. Missing/partial generated figures are rebuilt; completion markers
are removed before rebuilding. Changed inputs require a new output directory.
New results append to the batch manifest; skipped results are not duplicated.
Do not run two renderer processes against the same output directory.

The renderer now caps effective raster DPI at a default 40-megapixel budget with
a margin for tight bounding boxes. This also limits embedded structure-image
rasterization inside PDF/SVG. Heatmap cells/text stay vector in PDF/SVG.
Large figures may use less than 300 dpi; used DPI is recorded in new QC files.
--max-raster-megapixels can adjust the budget. Completed earlier files are retained.

Observed interrupted run 20261013_141015:
all 116 complete-residue generations and all 116 orderings completed.
696 standalone plots completed; 54 of 232 composite plots completed through DLST.
The next protein DMXL1 would produce an approximately 882-million-pixel composite
at the previous unconditional 300 dpi. This is a likely memory bottleneck, not a
confirmed cause of process termination. No renderer process was found running
during inspection. The new resume/memory changes were not executed here.

Run inside WSL with the existing environment:

~~~bash
evan='/mnt/c/Partners HealthCare Dropbox/Evangelos Papadopoulos/MMSEGs Discussion/GeneDeepDives/Evangelos'

"$HOME/.venvs/mmseg-overlap-v2/bin/python" \
  "/mnt/g/My Drive/VSCode_Github/BioPython/Morpheus/csv_heatmap_ordered.py" \
  --batch-dir "$evan/Pocket_Overlap_Ordered_v2_WSL_all_fullresidues_min50_diagonal_20261013_141015" \
  --out-dir "$evan/Figures/Pocket_Overlap_Morpheus_v2_WSL_all_fullresidues_min50_diagonal_20261013_141015" \
  --mode both --formats pdf svg png \
  --composite \
  --structure-reference-dir "/mnt/c/Partners HealthCare Dropbox/Evangelos Papadopoulos/LAB-Evan/Manuscript MMSEG pockets/Figures/Comparison_Plots_all" \
  --resume
~~~

Generation and ordering are not repeated. Resume writes into the exact existing
figure directory, not a new timestamped directory. SKIP messages are expected.


## Resume fix: equivalent color-config paths (2026-10-13)

The first resume version incorrectly compared color-config path strings:
Intersection_graphs/../Morpheus/... versus Morpheus/... . These identify the same
file. Resume now resolves both paths before comparing and new QC includes a
color-file content hash. Existing QC without color hashes remains supported.
Input hashes, layout, axes and other resume checks are still enforced.

The latest pasted output ended with Killed. The WSL diagnostic did not confirm
the termination cause, so memory exhaustion remains a likely explanation.
Raster budget is now 20 MP, applied to initial figure creation as well as saving;
closed figures/axes are released and garbage collection runs between figures.
Vector cells/text retain full detail. Existing completed outputs are preserved.
Rerun the same --resume command above. No rendering was executed during this fix.


## Large-composite fix and WSL memory (2026-10-13)

A further adjustment caps composites at 32 x 24 inches. For composite panels above
4000 cells, only nonzero cell values are annotated; zero cells remain in the matrix
and use the zero color. Axis labels and nonzero text scale with matrix size. Detailed
individual plots and all CSV values are retained. This avoids allocating tens of
thousands of redundant zero text artists. Small figures retain all annotations.
Existing completed files are skipped by --resume; partial outputs are rebuilt.
These changes were not executed or benchmarked during preparation.

The inspected Windows machine has about 128 GB physical RAM, while its existing
C:\Users\INFORMATICS\.wslconfig sets memory=12GB. Increasing to memory=48GB and adding
swap=16GB is an optional setting change for this machine. Preserve other settings.
Use Windows PowerShell to open the file:
notepad "$env:USERPROFILE\.wslconfig"

Under the existing [wsl2] section, change memory and add/update swap:
memory=48GB
swap=16GB

After other WSL jobs are stopped, apply with wsl --shutdown and restart WSL.
Shutdown terminates all WSL distributions/jobs. Recheck the C:/G: mounts before
running the existing exact --resume command above. No settings were changed here.
Microsoft documentation:
https://learn.microsoft.com/en-us/windows/wsl/wsl-config#main-wsl-settings

Using Windows Python from PowerShell is possible but is not required for this
resume; the current command uses the existing Linux Python environment. Switching
hosts also changes recorded path forms, so keep WSL for the current output run.
