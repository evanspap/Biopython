"""Read-only scientific invariant checks plus trial report and flowchart export."""
import csv
import itertools
import json
import random
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from three_method_overlap_order_v3 import assignment

b=Path('/mnt/c/Partners HealthCare Dropbox/Evangelos Papadopoulos/MMSEGs Discussion/GeneDeepDives/Evangelos')
old=b/'Pocket_Overlap_Ordered_v2_WSL_all_fullresidues_min50_diagonal_20261013_141015/ZNF592'
new=b/'Pocket_Overlap_Ordered_v3_ZNF592_20261021_trial1/ZNF592'
out=b/'Figures/Pocket_Overlap_Morpheus_v3_ZNF592_20261021_trial1'
rng=random.Random(11)
for n in range(1,7):
    for _ in range(8):
        cost=[[rng.randrange(-100,100) for j in range(n)] for i in range(n)]
        a=assignment(cost)
        actual=sum(cost[a[j]][j] for j in range(n))
        expected=min(sum(cost[p[j]][j] for j in range(n)) for p in itertools.permutations(range(n)))
        assert actual==expected

def load(path):
    with path.open() as f:
        rows=list(csv.reader(f))
    return {(r[0],c):float(v) for r in rows[1:] for c,v in zip(rows[0][1:],r[1:])}

def stats(path):
    with path.open() as f:
        rows=list(csv.reader(f))
    values=[float(rows[i][i]) for i in range(1,min(len(rows),len(rows[0])))]
    return sum(values),sum(abs(v-100)<1e-8 for v in values),len(values)

report=['# ZNF592 V3 trial audit','',
'Full GenBank residue data and the global <50% filter are unchanged. Only ordering changes.',
'Pair percentage = 100 * |A intersection B| / min(|A|, |B|). A 100% score means containment, not necessarily identical pockets.',
'Primary objective: sum of the fixed-denominator three-pair mean across diagonal positions. Secondary objective: total pairwise shared-residue counts; these are not unique triple residues.',
'Percentages are compared at 1e-10 percentage-point precision. Counts never outweigh one primary score unit.',
'32 deterministic starts (input, greedy, 30 shuffled), seed 592; exact integer Hungarian assignment for each axis with the others fixed. No seed is locked. Up to 100 coordinate sweeps per start.',
'Incomplete positions use zero contributions for absent pairs in the optimization only. Missing remains missing in exported matrices.',
'Final positions are ranked by mean percentage, then shared counts, only within groups with the same active methods to preserve contiguous axes.',
'This is a multistart coordinate heuristic, not a proven global optimum and not a maximization of the number of 100% cells.',
'', '| Comparison | Old diagonal percentage sum | New diagonal percentage sum | Old exact-100 cells | New exact-100 cells |',
'|---|---:|---:|---:|---:|']
for pair in ('AutoSite_FPocket','AutoSite_DoGSite','FPocket_DoGSite'):
    for metric in ('percentage','counts'):
        name=f'ZNF592_{pair}_ordered_{metric}.csv'
        assert load(old/name)==load(new/name),name
    name=f'ZNF592_{pair}_ordered_percentage.csv'
    a,c,_=stats(old/name);d,e,_=stats(new/name)
    report.append(f'| {pair} | {a:.6f} | {d:.6f} | {c} | {e} |')
a=json.loads((old/'ZNF592_overlap_order_qc.json').read_text())
z=json.loads((new/'ZNF592_overlap_order_qc.json').read_text())
assert a['removed_pockets']==z['removed_pockets']
for m in a['orders']:
    assert set(a['orders'][m])==set(z['orders'][m])
assert z['diagonal_after'][0]>=a['diagonal_after'][0]
report+=['','Validation: all six matrices contain exactly the same pocket-keyed values; retained/removed IDs unchanged. Hungarian solver matched exhaustive enumeration on 48 random matrices of sizes 1-6.',
         f"Old objective: {a['diagonal_after']}; new objective: {z['diagonal_after']}.",
         'Fpockt20 and DGsite04 remain off-diagonal: old positions 11/16; new positions 16/11. This example is not forced; multiple FPocket pockets compete for each DoGSite position.']
(out/'ZNF592_v3_method_and_audit.md').write_text('\n'.join(report)+'\n')

boxes=[
('Validated inputs','Full GenBank pairwise percentages + exact shared-residue counts'),
('Keep the same pockets','Remove a pocket only when all comparisons to both other methods are <50%'),
('Define priorities','Primary: total diagonal mean of AS-FP, AS-DG, FP-DG percentages\nTie: total shared-residue counts; deterministic identifiers'),
('Initialize 32 arrangements','Input order + greedy order + 30 deterministic shuffled orders'),
('Optimize each method in turn','Hold two axes fixed; Hungarian assignment optimizes the third\nAccept strict improvement; repeat until stable or 100 sweeps'),
('Choose and rank','Keep the best start; rank positions by mean then counts\nReorder together only within equal active-method groups'),
('Validate and export','Same pocket IDs and cell values; new common orders, QC and figures\nOne-to-many matches can remain off-diagonal; no global optimum guarantee')]
fig,ax=plt.subplots(figsize=(11,12));ax.set_xlim(0,1);ax.set_ylim(0,1);ax.axis('off')
for i,(title,body) in enumerate(boxes):
    y=.94-i*.135
    ax.add_patch(FancyBboxPatch((.025,y-.09),.95,.105,boxstyle='round,pad=0.01',facecolor='#edf3fa',edgecolor='#274967'))
    ax.text(.5,y-.006,title,ha='center',va='center',fontsize=12,fontweight='bold')
    ax.text(.5,y-.052,body,ha='center',va='center',fontsize=9)
    if i<6: ax.annotate('',xy=(.5,y-.116),xytext=(.5,y-.092),arrowprops={'arrowstyle':'->','color':'#274967'})
fig.suptitle('ZNF592: common three-method ordering (V3)',fontsize=15,fontweight='bold')
for ext in ('pdf','svg','png'):
    fig.savefig(out/f'ZNF592_v3_ordering_flowchart.{ext}',dpi=180,bbox_inches='tight')
plt.close(fig)
print('\n'.join(report))