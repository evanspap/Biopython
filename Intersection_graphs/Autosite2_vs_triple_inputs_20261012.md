# Autosite_2 versus stored triple-analysis inputs — 2026-10-12

## Conclusion

The AutoSite annotations in the stored merged GenBank inputs match the corrected
B-filtered AutoSite annotations for all 116 proteins. They do not represent the
Autosite_2 run. This establishes the source of the existing triple-analysis inputs;
it does not independently decide which structure-filtering protocol is scientifically
preferable.

No detection, triple analysis, overlap ordering, or figure pipeline was rerun.
The comparison reads stored GenBank FEATURES and ORIGIN, normalizes PocketNN and
ASpockNN to numeric pocket IDs, and compares complete residue sets including wrapped
locations. The attached CSV contains per-protein results.

## Full comparison

- Merged input GenBank files inspected: 116.
- Corrected B-filtered source files available: 116.
- Exact agreement of pocket IDs and residue sets, merged versus B-filtered: 116/116.
- Proteins with nonempty merged AutoSite annotations: 112; all 112 match B-filtered.
- Autosite_2 source files at the expected path: 115/116 (DMXL1 absent at that path).
- Exact merged versus Autosite_2 agreement: only XBP1, where both have zero pockets.
- Nonempty merged versus Autosite_2 agreement: 0.
- Full ORIGIN sequence equality: merged/B-filtered 115/116; merged/Autosite_2 114/115 available files.
- MTHFD1 has the same sequence length (935) but different sequence text in these
  versions, while its merged/B-filtered pocket residue sets match exactly. This
  additional sequence discrepancy is recorded for follow-up.

"Exact agreement" here means the complete pocket-ID/residue-set collection agrees,
not that no individual pockets overlap between separate prediction runs.

## Selected proteins

| Protein | Merged triple input AS pockets | Corrected B-filtered AS pockets | Autosite_2 AS pockets |
|---|---:|---:|---:|
| ADSS2 | 15 | 15 | 16 |
| AK2 | 11 | 11 | 10 |
| EP300 | 48 | 48 | 80 |
| MAF | 0 | 0 | 1 |
| POU2AF1 | 0 | 0 | 2 |
| UBALD2 | 0 | 0 | 1 |
| XBP1 | 0 | 0 | 0 |

For these seven proteins, ORIGIN sequences agree across the compared versions.
That does not make their independently predicted pocket sets interchangeable.

Example ADSS2 Pocket01 location begins:
Merged/B-filtered:
join(38,39,40,41,42,43,44,45,46,48,49,65,66,67,68,69,70,71,...)
Autosite_2:
join(38,39,40,42,43,65,66,67,68,80,144,147,151,159,160,161,...)

This is a real residue-set difference, not just different folder names.

## Run-log chronology (directly checked examples)

For ADSS2, AK2 and EP300:
- Autosite_2 logs: 2024-10-14; respectively 16, 10 and 80 clusters.
- B-filtered logs: 2025-03-05; respectively 15, 11 and 48 clusters.

The same four problem proteins' B-filtered logs explicitly report zero clusters.
Their triple_pocket_analysis/PROTEIN.out files exist but contain only the header.
They are absent from the inspected combined triple outputs and July-2026 workbook.
This is consistent with no possible three-method combinations in that stored run,
not missing rows introduced by the final Excel aggregation.

## Exact source patterns

Root:
C:\Partners HealthCare Dropbox\Evangelos Papadopoulos\MMSEGs Discussion\GeneDeepDives\Evangelos

Stored merged input:
Genebank_Annotations\merge_AFB_Fpock_DGSite_ASpock_45nm\PROTEIN.gb

Matching corrected B-filtered AutoSite source:
Pockets_from_B-filtered_AF\Autosite\pockets_B-filtered_45nm\PROTEIN_corr.gb

Different Autosite_2 source:
Autosite_2\PROTEIN_reduce\PROTEIN_reduce_45nm.gb

Stored triple output:
Genebank_Annotations\merge_AFB_Fpock_DGSite_ASpock_45nm\triple_pocket_analysis\PROTEIN.out

The canonical triple_pocket_analysis.py reads a supplied merged GenBank and
enumerates Fpockt x DGsite x ASpock pockets. It requires an AutoSite member to
form a triple and writes combinations only when at least one directional
three-way overlap is strictly greater than 0.5.

## Interpretation

For reproducing the existing triple results, keep the B-filtered source.
Introducing Autosite_2 pockets for MAF/POU2AF1/UBALD2 would constitute a different
analysis, not recovery of the existing final annotations. Any decision to use
Autosite_2 should apply a documented consistent workflow and recompute the affected
comparisons/consensus. Do not silently replace only selected missing pockets.

Data:
Autosite2_vs_triple_inputs_20261012.csv
