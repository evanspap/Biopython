"""V3: deterministic multistart coordinate assignment; percentage first, counts second.
Reuses validated V2 parsing/export. No original data or legacy outputs are modified.
"""
import itertools
import json
import random
import sys
from pathlib import Path
import three_method_overlap_order as base

STARTS = 32
SCALE = 10**10
RUN = {}


def assignment(cost):
    """Exact integer Hungarian minimization, square matrix, deterministic ties."""
    n = len(cost)
    u, v, p, way = [0]*(n+1), [0]*(n+1), [0]*(n+1), [0]*(n+1)
    for i in range(1, n+1):
        p[0] = i
        j0 = 0
        mins, used = [None]*(n+1), [False]*(n+1)
        while True:
            used[j0] = True
            i0, delta, j1 = p[j0], None, 0
            for j in range(1, n+1):
                if used[j]:
                    continue
                cur = cost[i0-1][j-1] - u[i0] - v[j]
                if mins[j] is None or cur < mins[j]:
                    mins[j], way[j] = cur, j0
                if delta is None or mins[j] < delta:
                    delta, j1 = mins[j], j
            for j in range(n+1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                elif mins[j] is not None:
                    mins[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while True:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            if j0 == 0:
                break
    result = [0]*n
    for j in range(1, n+1):
        result[j-1] = p[j]-1
    return result


def key(graph, order):
    pct = count = 0
    for i in range(max(map(len, order.values()), default=0)):
        triad = base.at(order, i)
        for a, b in itertools.combinations(triad, 2):
            e = graph.edge(a, b)
            pct += round((e.pct or 0)*SCALE)
            count += e.count or 0
    return pct, count


def rank_slots(graph, order):
    # Reorder only positions with the same active methods: no internal gaps.
    groups = {}
    for i in range(max(map(len, order.values()), default=0)):
        t = base.at(order, i)
        groups.setdefault(tuple(p is not None for p in t), []).append(t)
    result = {m: [] for m in base.METHODS}
    for mask, triads in groups.items():
        triads.sort(key=lambda t: (-round(graph.score(t)[0]*SCALE),
                                  -graph.score(t)[1], tuple(p or '' for p in t)))
        for t in triads:
            for m, pocket in zip(base.METHODS, t):
                if pocket is not None:
                    result[m].append(pocket)
    return result


old_greedy = base.greedy_order

def arrange(graph, ids):
    baseline = {m: list(ids[m]) for m in base.METHODS}
    greedy, _ = old_greedy(graph, ids)
    # Bounds total residue tie contribution across every possible diagonal.
    bound = 1 + sum(e.count or 0 for (a,b), e in graph.edges.items() if a < b)
    rng = random.Random(592)
    starts = [baseline, greedy]
    for _ in range(STARTS-2):
        start = {m: sorted(ids[m]) for m in base.METHODS}
        for values in start.values():
            rng.shuffle(values)
        starts.append(start)
    best, best_key, trials = baseline, key(graph, baseline), []
    for start_index, start in enumerate(starts):
        order = {m: list(v) for m,v in start.items()}
        sweeps = 0
        for sweep in range(100):
            before = key(graph, order)
            for m in base.METHODS:
                pockets = order[m]
                cost = []
                for pocket in pockets:
                    row = []
                    for i in range(len(pockets)):
                        pct = count = 0
                        for other in base.METHODS:
                            if other == m:
                                continue
                            q = order[other][i] if i < len(order[other]) else None
                            e = graph.edge(pocket, q)
                            pct += round((e.pct or 0)*SCALE)
                            count += e.count or 0
                        row.append(-(pct*bound + count))
                    cost.append(row)
                candidate = [pockets[j] for j in assignment(cost)]
                prior = order[m]
                old_key = key(graph, order)
                order[m] = candidate
                if key(graph, order) <= old_key:
                    order[m] = prior
            sweeps = sweep+1
            if key(graph, order) == before:
                break
        order = rank_slots(graph, order)
        score = key(graph, order)
        trials.append({'start': start_index, 'sweeps': sweeps, 'objective': score,
                       'sweep_limit_reached': sweeps == 100})
        if score > best_key:
            best, best_key = order, score
    best = rank_slots(graph, best)
    RUN.update(trials=trials, baseline_objective=key(graph, baseline),
               best_objective=key(graph, best), algorithm='v3 multistart coordinate Hungarian assignment',
               starts=STARTS, random_seed=592, percentage_comparison_precision=1/SCALE,
               position_score='(AS-FP percentage + AS-DG percentage + FP-DG percentage)/3',
               absent_policy='absent pairs contribute zero to fixed-denominator mean; never written as measured zeros',
               tie_break='sum of three shared-residue counts, then deterministic pocket IDs',
               slot_sort='descending mean, then counts, within equal active-method masks',
               guarantee='nondecreasing from supplied baseline; coordinate local optimum unless sweep limit reached; no global guarantee')
    return best, None


old_process = base.process

def process(protein, paths, count_paths, output, args):
    if args.selected_priority:
        raise ValueError('V3 does not support selected-priority constraints')
    old_process(protein, paths, count_paths, output, args)
    path = Path(output) / f'{protein}_overlap_order_qc.json'
    qc = json.loads(path.read_text(encoding='utf-8'))
    qc.update(ordering_algorithm=RUN, optimization=RUN['algorithm'],
              score_definition='lexicographic total diagonal percentage (equivalent to fixed-denominator mean), then total shared residues')
    path.write_text(json.dumps(qc, indent=2)+'\n', encoding='utf-8')


base.greedy_order = arrange
base.optimize = lambda *args, **kwargs: 0
base.process = process
if __name__ == '__main__':
    sys.exit(base.main())