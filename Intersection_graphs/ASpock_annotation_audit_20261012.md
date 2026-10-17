# Audit of missing AutoSite / triple-consensus annotations — 2026-10-12

Scope: MAF, POU2AF1, UBALD2, and XBP1. Read-only inspection of stored data;
no pocket detection, overlap analysis, ordering, or rendering was executed.

## Findings

The stored B-filtered workflow contains FPocket and DoGSite annotations for all
four proteins, but no AutoSite annotations. This agrees with the AutoSite run logs
(0 clusters), the empty AutoSite GenBank FEATURES sections, and the merged GenBank
files. Their individual triple-analysis .out files exist but contain only a CSV
header, with no result rows. The absence therefore predates Excel aggregation.

| Protein | FPocket features in merged GB | DoGSite features in merged GB | AutoSite features in merged GB | Triple .out data rows | AutoSite pockets in separate Autosite_2 run |
|---|---:|---:|---:|---:|---:|
| MAF | 7 | 7 | 0 | 0 | 1 |
| POU2AF1 | 2 | 9 | 0 | 0 | 2 |
| UBALD2 | 5 | 4 | 0 | 0 | 1 |
| XBP1 | 12 | 16 | 0 | 0 | 0 |

This supports "no triple combinations formed in the stored B-filtered run",
rather than annotations being lost during the final merging/Excel step.
It does not establish that no other historical run exists anywhere.

## Paths

All paths below are relative to:

C:\Partners HealthCare Dropbox\Evangelos Papadopoulos\MMSEGs Discussion\GeneDeepDives\Evangelos

Replace PROTEIN with MAF, POU2AF1, UBALD2 or XBP1.

B-filtered AutoSite run log:
Pockets_from_B-filtered_AF\Autosite\PROTEIN_AutoSiteSummary.log

All four logs explicitly say:
AutoSite identified and characterized 0 clusters.

B-filtered AutoSite annotation files:
Pockets_from_B-filtered_AF\Autosite\pockets_B-filtered_45nm\PROTEIN_corr.gb
Genebank_Annotations\ASpock_45nm\PROTEIN.gb
Genebank_Annotations\ASpock_45nm\ASpock\PROTEIN.gb
Genebank_Annotations\Autosite_Pockets\PROTEIN.gb

The inspected standardized AutoSite GenBank files have no pocket FEATURES.

Merged input containing method-specific features:
Genebank_Annotations\merge_AFB_Fpock_DGSite_ASpock_45nm\PROTEIN.gb
Genebank_Annotations\merge_AFB_Fpock_DGSite_ASpock_45nm\Fpockt\PROTEIN.gb

Existing triple-analysis output:
Genebank_Annotations\merge_AFB_Fpock_DGSite_ASpock_45nm\triple_pocket_analysis\PROTEIN.out

These four .out files contain only:
Protein,max%,Pocket_A,Nr_Fp,Drug Score,Pocket Score,Pocket_B,Nr_DG,drugScore,depth,Pocket_C,Nr_Asp,e,v*buriedness^2/rg,common residues,%Fpockt,%DGsite,%ASpock,list

No matches for these proteins were found in:
Genebank_Annotations\merge_AFB_Fpock_DGSite_ASpock_45nm\triple_pocket_analysis\combined_triple_pocket_analysis.out
Genebank_Annotations\merge_AFB_Fpock_DGSite_ASpock_45nm\triple_pocket_analysis\filtered_combined_triple_pocket_analysis_20251107_v01.csv
Genebank_Annotations\merge_AFB_Fpock_DGSite_ASpock_45nm\filtered_20251107\filtered_combined_triple_pocket_analysis_20251107_v02.csv

The previous read-only inspection of combined_triple_pocket_analysis_20260720_v2.xlsx
also found no rows for these four proteins across its five sheets.

## Separate Autosite_2 run

Pocket definitions and scores are available in:
Autosite_2\PROTEIN_reduce\PROTEIN_reduce_45nm.gb
Autosite_2\PROTEIN_reduce\PROTEIN_reduce_pockets.gbk
Autosite_2\PROTEIN_reduce\PROTEIN_reduce_summary.csv

MAF: Pocket01, residues 103–110.
POU2AF1: Pocket01 (125,127,133,139,140,141,142,143,144,145,146,147);
Pocket02 (123,124,125,137,140,141).
UBALD2: Pocket01, residues 48–55.
XBP1: empty summary/FEATURES; the Autosite_2 log explicitly reports 0 clusters.

These annotations belong to a separate structure/run version. They must not be
silently substituted for the B-filtered workflow's missing pockets. A separate
recomputed comparison would need verified coordinate/residue-numbering and source
compatibility. No such substitution or recomputation was performed.

## Implication for existing figures

The B-filtered data support FPocket–DoGSite comparisons for these proteins.
AutoSite comparisons should show that no pockets were detected in that run;
they are not zero-overlap comparisons between existing AutoSite pockets.
A structure + three comparison figure can still be generated with an explanatory
empty AutoSite panel. The ordering parser's empty-method support is appropriate
for this stored workflow, and QC records methods_without_input_pockets.
