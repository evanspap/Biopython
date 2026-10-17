#!/usr/bin/env python3
"""Coordinated pocket ordering, v2. Read inputs only; write a new output directory.

Accepts legacy *_relative.csv (two directions + residue totals), legacy
*_relative_simple.csv, or plain rectangular matrices. Percentages are never
recomputed. Counts must be supplied or uniquely recoverable from percentages
and pocket totals. Missing cells stay missing. Python >= 3.10, standard library.
"""
from __future__ import annotations

import argparse
import csv
import functools
import hashlib
import itertools
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path



def local_path(value):
    """Use extended Windows paths for long Dropbox filenames; ordinary paths elsewhere."""
    value = str(Path(value).expanduser())
    if sys.platform == "win32":
        import ntpath
        value = ntpath.abspath(value)
        if not value.startswith("\\\\?\\"):
            value = "\\\\?\\UNC\\" + value[2:] if value.startswith("\\\\") else "\\\\?\\" + value
    return Path(value)

VERSION = "2.0.0"
METHODS = ("AutoSite", "FPocket", "DoGSite")
PREFIXES = {"AutoSite": "ASpock", "FPocket": "Fpockt", "DoGSite": "DGsite"}
PAIRS = (("AutoSite", "FPocket"), ("AutoSite", "DoGSite"), ("FPocket", "DoGSite"))
FOLDERS = ("Fpockt_ASpock", "ASpock_DGsite", "Fpockt_DGsite")
MISSING = {"", "-", "na", "nan", "null", "none"}
META = "totalresidues"


def meta(label):
    return re.sub(r"[\s_]", "", label).lower() == META


def method(label):
    for name, prefix in PREFIXES.items():
        if label.lower().startswith(prefix.lower()):
            return name
    return None


def number(text, *, count=False, context=""):
    text = text.strip()
    if text.lower() in MISSING:
        return None
    try:
        value = float(text.removesuffix("%"))
    except ValueError as exc:
        raise ValueError(f"{context}: invalid number {text!r}") from exc
    if not math.isfinite(value) or value < 0 or (not count and value > 100):
        raise ValueError(f"{context}: out-of-range value {text!r}")
    if count:
        if not value.is_integer():
            raise ValueError(f"{context}: residue count must be an integer")
        return int(value)
    return value


@dataclass
class Table:
    rows: list[str]
    cols: list[str]
    values: dict
    totals: dict
    path: Path


def read_table(path, *, count=False):
    path = Path(path)
    with path.open(encoding="utf-8-sig", newline="") as stream:
        records = [r for r in csv.reader(stream) if any(c.strip() for c in r)]
    if not records:
        raise ValueError(f"{path}: empty CSV")
    headers = [c.strip() for c in records[0]]
    if len(headers) < 1 or len(headers[1:]) != len(set(headers[1:])):
        raise ValueError(f"{path}: duplicate column identifiers")
    cols = [c for c in headers[1:] if not meta(c)]
    if any(not c for c in cols):
        raise ValueError(f"{path}: empty pocket column identifier")
    rows, values, totals = [], {}, {}
    seen = set()
    for line, record in enumerate(records[1:], 2):
        if len(record) != len(headers):
            raise ValueError(f"{path}:{line}: expected {len(headers)} CSV fields")
        row = record[0].strip()
        if not row or row in seen:
            raise ValueError(f"{path}:{line}: empty/duplicate row identifier {row!r}")
        seen.add(row)
        is_total = meta(row)
        if not is_total:
            rows.append(row)
        for col, text in zip(headers[1:], record[1:]):
            if is_total and meta(col):
                continue
            if is_total or meta(col):
                value = number(text, count=True, context=f"{path}:{line}:{col}")
                pocket = col if is_total else row
                if value is not None:
                    if pocket in totals and totals[pocket] != value:
                        raise ValueError(f"{path}: inconsistent total for {pocket}")
                    totals[pocket] = value
            else:
                values[row, col] = number(
                    text, count=count, context=f"{path}:{line}:{col}"
                )
    return Table(rows, cols, values, totals, path)


def orient(table, left, right, *, count=False):
    rl = [p for p in table.rows if method(p) == left]
    rr = [p for p in table.rows if method(p) == right]
    cl = [p for p in table.cols if method(p) == left]
    cr = [p for p in table.cols if method(p) == right]
    values = {}
    if not table.rows and not table.cols:
        return Table([], [], {}, table.totals, table.path)
    # A legacy square matrix can contain only one method (e.g. no AutoSite pockets).
    # Classify by the full square axes, not by requiring nonempty pocket groups.
    if set(table.rows) == set(table.cols) and (rl or rr or cl or cr):
        if set(rl) != set(cl) or set(rr) != set(cr):
            raise ValueError(f"{table.path}: inconsistent square-matrix axes")
        rows, cols = rl, cr
        for a, b in itertools.product(rows, cols):
            x, y = table.values[a, b], table.values[b, a]
            if count:
                if x != y:
                    raise ValueError(f"{table.path}: asymmetric count for {a}, {b}")
                values[a, b] = x
            else:
                # An absent direction cannot be silently replaced by zero.
                values[a, b] = None if x is None or y is None else max(x, y)
    elif (not rr and not cl and (rl or cr)
          and all(method(p) == left for p in table.rows)
          and all(method(p) == right for p in table.cols)):
        rows, cols = rl, cr
        values = {(a, b): table.values[a, b] for a, b in itertools.product(rows, cols)}
    elif (not rl and not cr and (rr or cl)
          and all(method(p) == right for p in table.rows)
          and all(method(p) == left for p in table.cols)):
        rows, cols = cl, rr
        values = {(a, b): table.values[b, a] for a, b in itertools.product(rows, cols)}
    else:
        raise ValueError(
            f"{table.path}: cannot identify {left}/{right} axes; "
            "expected ASpock*, Fpockt*, DGsite* identifiers"
        )
    return Table(rows, cols, values, table.totals, table.path)


def expected_percentage(shared, na, nb, policy):
    if min(na, nb) == 0:
        return 0.0
    x, y = 100 * (shared / na), 100 * (shared / nb)
    return max(int(x), int(y)) if policy == "legacy-floor" else max(x, y)


def matching_counts(pct, na, nb, policy):
    """Enumerate a bounded integer interval, then verify with the actual formula.

    int((shared / total) * 100) matches the legacy generator's floating-point
    operation order, including its behavior at integer boundaries.
    """
    n = min(na, nb)
    if n == 0:
        return [0] if pct == 0 else []
    if policy == "legacy-floor" and not float(pct).is_integer():
        return []
    lo = max(0, math.floor(pct * n / 100) - 2)
    hi = min(n, math.ceil((pct + (1 if policy == "legacy-floor" else 0)) * n / 100) + 2)
    return [
        c for c in range(lo, hi + 1)
        if math.isclose(expected_percentage(c, na, nb, policy), pct,
                        rel_tol=0, abs_tol=1e-8)
    ]


@dataclass
class Edge:
    pct: float | None
    count: int | None
    source: str


class Graph:
    def __init__(self, tables, count_paths, policy):
        self.tables = tables
        self.ids, self.totals, self.edges = {}, {}, {}
        self.sources = [t.path for t in tables]
        self.policy = policy
        for pair, table, count_path in zip(PAIRS, tables, count_paths):
            for name, ids in zip(pair, (table.rows, table.cols)):
                if name in self.ids and set(ids) != set(self.ids[name]):
                    raise ValueError(f"Pocket identifiers for {name} differ between matrices")
                self.ids.setdefault(name, list(ids))
            for pocket, total in table.totals.items():
                if method(pocket) not in pair:
                    continue
                if pocket in self.totals and self.totals[pocket] != total:
                    raise ValueError(f"Conflicting residue totals for {pocket}")
                self.totals[pocket] = total
        ambiguous = []
        for pair, table, count_path in zip(PAIRS, tables, count_paths):
            counts = None
            if count_path:
                counts = orient(read_table(count_path, count=True), *pair, count=True)
                if set(counts.rows) != set(table.rows) or set(counts.cols) != set(table.cols):
                    raise ValueError(f"{counts.path}: count axes differ from percentage axes")
                self.sources.append(counts.path)
            for a, b in itertools.product(table.rows, table.cols):
                pct = table.values[a, b]
                na, nb = self.totals.get(a), self.totals.get(b)
                count, source = None, "missing"
                if counts:
                    count = counts.values[a, b]
                    if count is not None:
                        source = "input_count_matrix"
                if pct is not None and count is None:
                    if na is None or nb is None:
                        ambiguous.append(f"{a}/{b}: missing totals or count")
                    else:
                        candidates = matching_counts(pct, na, nb, policy)
                        if len(candidates) != 1:
                            ambiguous.append(f"{a}/{b}: {len(candidates)} compatible counts")
                        else:
                            count, source = candidates[0], "uniquely_inferred_from_percentage_and_totals"
                if count is not None:
                    if na is not None and count > na or nb is not None and count > nb:
                        raise ValueError(f"{a}/{b}: count exceeds pocket size")
                    if pct is not None and na is not None and nb is not None:
                        expected = expected_percentage(count, na, nb, policy)
                        if not math.isclose(expected, pct, rel_tol=0, abs_tol=1e-8):
                            raise ValueError(
                                f"{a}/{b}: count {count} implies {expected}%, input is {pct}%"
                            )
                self.edges[a, b] = Edge(pct, count, source)
                self.edges[b, a] = self.edges[a, b]
        if ambiguous:
            raise ValueError(
                "Exact shared-residue counts cannot be determined. Supply the corresponding "
                "--as-fp-counts/--as-dg-counts/--fp-dg-counts matrices.\n"
                + "\n".join(ambiguous[:12])
                + (f"\n... {len(ambiguous)} affected cells" if len(ambiguous) > 12 else "")
            )

    def edge(self, a, b):
        if a is None or b is None:
            return Edge(None, None, "absent_pocket")
        return self.edges[a, b]

    @functools.lru_cache(maxsize=None)
    def support(self, pocket, target=None):
        return sum(
            self.edge(pocket, other).pct or 0
            for name in METHODS if name != method(pocket) and (target is None or name == target)
            for other in self.ids[name]
        )

    def score(self, triad):
        edges = [self.edge(triad[i], triad[j]) for i, j in ((0, 1), (0, 2), (1, 2))]
        return (sum(e.pct or 0 for e in edges),
                sum(e.count or 0 for e in edges))

    def global_zero(self, pocket):
        edges = [self.edge(pocket, other)
                 for name in METHODS if name != method(pocket)
                 for other in self.ids[name]]
        return bool(edges) and all(e.pct == 0 and e.count == 0 for e in edges)

    def globally_below(self, pocket, threshold):
        edges = [self.edge(pocket, other)
                 for name in METHODS if name != method(pocket)
                 for other in self.ids[name]]
        return bool(edges) and all(e.pct is not None and e.pct < threshold for e in edges)


def triad_key(graph, triad):
    pct, count = graph.score(triad)
    return (-pct, -count, -sum(graph.support(p) for p in triad if p is not None),
            tuple(p or "" for p in triad))


def greedy_order(graph, ids):
    remaining = {m: sorted(ids[m]) for m in METHODS}
    result = {m: [] for m in METHODS}
    candidates = []
    for left, right in PAIRS:
        third = next(m for m in METHODS if m not in (left, right))
        for a, b in itertools.product(remaining[left], remaining[right]):
            edge = graph.edge(a, b)
            if edge.pct is not None:
                key = (-edge.pct, -(edge.count or 0),
                       -(graph.support(a, third) + graph.support(b, third)), a, b)
                candidates.append((key, left, right, a, b))
    seed = None
    if candidates:
        _, left, right, a, b = min(candidates)
        third = next(m for m in METHODS if m not in (left, right))
        fixed = {left: a, right: b}
        options = remaining[third] or [None]
        triads = [tuple({**fixed, third: p}.get(m) for m in METHODS) for p in options]
        seed = min(triads, key=lambda triad: triad_key(graph, triad))

    def append(triad):
        for m, pocket in zip(METHODS, triad):
            if pocket is not None:
                result[m].append(pocket)
                remaining[m].remove(pocket)

    if seed is not None:
        append(seed)
    while any(remaining.values()):
        options = [remaining[m] or [None] for m in METHODS]
        best = min(itertools.product(*options), key=lambda t: triad_key(graph, t))
        append(best)
    return result, seed


def at(order, position):
    return tuple(order[m][position] if position < len(order[m]) else None for m in METHODS)


def objective(graph, order):
    positions = max((len(v) for v in order.values()), default=0)
    scores = [graph.score(at(order, i)) for i in range(positions)]
    return (sum(s[0] for s in scores), sum(s[1] for s in scores))


def optimize(graph, order, passes, lock_seed, locked_prefix=None):
    """Deterministic best improving swap per pass; keep each axis contiguous."""
    swaps = 0
    for _ in range(passes):
        best, move = (0.0, 0), None
        for m in METHODS:
            for i in range(max(1 if lock_seed else 0, (locked_prefix or {}).get(m, 0)), len(order[m])):
                for j in range(i + 1, len(order[m])):
                    old = [graph.score(at(order, k)) for k in (i, j)]
                    order[m][i], order[m][j] = order[m][j], order[m][i]
                    new = [graph.score(at(order, k)) for k in (i, j)]
                    order[m][i], order[m][j] = order[m][j], order[m][i]
                    delta = (sum(s[0] for s in new) - sum(s[0] for s in old),
                             sum(s[1] for s in new) - sum(s[1] for s in old))
                    if delta > best:
                        best, move = delta, (m, i, j)
        if move is None:
            break
        m, i, j = move
        order[m][i], order[m][j] = order[m][j], order[m][i]
        swaps += 1
    return swaps


def write_csv(path, graph, order, pair, metric):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["PocketID", *order[pair[1]]])
        for a in order[pair[0]]:
            writer.writerow([a, *[
                getattr(graph.edge(a, b), metric) if getattr(graph.edge(a, b), metric) is not None else ""
                for b in order[pair[1]]
            ]])



def selected_orders(path, protein, graph, retained):
    """TSV row order -> unique first-appearance order for each method."""
    prefix = {m: [] for m in METHODS}
    records = []
    if path is None:
        return prefix, records
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"Protein", "Pocket_A", "Pocket_B", "Pocket_C"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError(f"{path}: requires TSV columns {sorted(required)}")
        for line, row in enumerate(reader, 2):
            if row["Protein"].strip() != protein:
                continue
            triad = {}
            for name, column in (("AutoSite", "Pocket_C"), ("FPocket", "Pocket_A"), ("DoGSite", "Pocket_B")):
                token = row[column].strip()
                match = re.fullmatch(re.escape(PREFIXES[name]) + r"_?0*(\d+)", token, re.IGNORECASE)
                if not match:
                    raise ValueError(f"{path}:{line}: invalid {column} identifier {token!r}")
                candidates = [p for p in graph.ids[name]
                              if re.fullmatch(re.escape(PREFIXES[name]) + r"0*" + str(int(match[1])), p, re.IGNORECASE)]
                if len(candidates) != 1:
                    raise ValueError(f"{path}:{line}: {token} does not identify exactly one input pocket")
                triad[name] = candidates[0]
            filtered = [p for name, p in triad.items() if p not in retained[name]]
            record = {"tsv_line": line, "selected_pocket": row.get("Pocket", ""),
                      "pockets": triad, "filtered_constituents": filtered}
            if not filtered:
                for name in METHODS:
                    if triad[name] not in prefix[name]:
                        prefix[name].append(triad[name])
            records.append(record)
    return prefix, records


def process(protein, paths, count_paths, output, args):
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", protein) or protein in (".", ".."):
        raise ValueError(f"Unsafe protein identifier {protein!r}")
    output = Path(output)
    # Refuse output reuse: input/legacy files and earlier runs stay reviewable.
    if output.exists():
        raise ValueError(f"Output already exists: {output}; choose a new --out-dir")
    tables = [orient(read_table(p), *pair) for p, pair in zip(paths, PAIRS)]
    graph = Graph(tables, count_paths, args.percentage_policy)
    removal_reasons = {}
    for m in METHODS:
        for pocket in graph.ids[m]:
            if args.drop_globally_zero and graph.global_zero(pocket):
                removal_reasons[pocket] = "known_zero_in_both_comparisons"
            elif args.drop_globally_below is not None and graph.globally_below(pocket, args.drop_globally_below):
                removal_reasons[pocket] = f"all_known_percentages_below_{args.drop_globally_below:g}_in_both_comparisons"
    removed = {m: [p for p in graph.ids[m] if p in removal_reasons] for m in METHODS}
    baseline = {m: [p for p in graph.ids[m] if p not in removed[m]] for m in METHODS}
    before = objective(graph, baseline)
    priority, selected_records = selected_orders(args.selected_priority, protein, graph, baseline)
    remainder = {m: [p for p in baseline[m] if p not in priority[m]] for m in METHODS}
    tail, seed = greedy_order(graph, remainder)
    order = {m: priority[m] + tail[m] for m in METHODS}
    has_priority = any(priority.values())
    after_greedy = objective(graph, order)
    swaps = optimize(graph, order, args.swap_passes, seed is not None and not has_priority,
                     {m: len(priority[m]) for m in METHODS} if has_priority else None)
    for record in selected_records:
        record['final_positions'] = {m: order[m].index(p) + 1 if p in order[m] else None
                                     for m, p in record['pockets'].items()}
        positions = list(record['final_positions'].values())
        record['aligned_on_diagonal'] = None not in positions and len(set(positions)) == 1
    after = objective(graph, order)
    # Internal invariants are checked when the user runs this script.
    for m in METHODS:
        if len(order[m]) != len(set(order[m])) or set(order[m]) != set(baseline[m]):
            raise RuntimeError(f"Invalid {m} permutation")
    if after < after_greedy:
        raise RuntimeError("Swap optimization reduced the objective")
    output.mkdir(parents=True)
    try:
        outputs = []
        for pair in PAIRS:
            tag = "_".join(pair)
            for label, metric in (("percentage", "pct"), ("counts", "count")):
                filename = f"{protein}_{tag}_ordered_{label}.csv"
                write_csv(output / filename, graph, order, pair, metric)
                outputs.append(filename)
        for m in METHODS:
            with (output / f"{protein}_{m}_pocket_order.tsv").open("w", encoding="utf-8", newline="") as stream:
                writer = csv.writer(stream, delimiter="\t")
                writer.writerow(["position", "pocket_id", "total_residues"])
                writer.writerows((i, p, graph.totals.get(p, "")) for i, p in enumerate(order[m], 1))
        manifest = output / f"{protein}_overlap_order_manifest.tsv"
        fields = ["record_type", "diagonal_position", *METHODS]
        fields += [f"{'_'.join(pair)}_{metric}" for pair in PAIRS
                   for metric in ("percentage", "shared_residues", "count_source")]
        fields += ["position_percentage_score", "position_count_score",
                   "removed_method", "removed_pocket", "removal_reason",
                   "diagonal_percentage_before", "diagonal_percentage_after",
                   "diagonal_count_before", "diagonal_count_after"]
        with manifest.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fields, delimiter="\t")
            writer.writeheader()
            summary = dict(diagonal_percentage_before=before[0], diagonal_percentage_after=after[0],
                           diagonal_count_before=before[1], diagonal_count_after=after[1])
            writer.writerow({"record_type": "summary", **summary})
            for i in range(max((len(v) for v in order.values()), default=0)):
                triad = at(order, i)
                record = {"record_type": "diagonal", "diagonal_position": i + 1,
                          **dict(zip(METHODS, (p or "" for p in triad))), **summary}
                for pair in PAIRS:
                    edge = graph.edge(triad[METHODS.index(pair[0])], triad[METHODS.index(pair[1])])
                    tag = "_".join(pair)
                    record[f"{tag}_percentage"] = "" if edge.pct is None else edge.pct
                    record[f"{tag}_shared_residues"] = "" if edge.count is None else edge.count
                    record[f"{tag}_count_source"] = edge.source
                score = graph.score(triad)
                record.update(position_percentage_score=score[0], position_count_score=score[1])
                writer.writerow(record)
            for m in METHODS:
                for p in removed[m]:
                    writer.writerow({"record_type": "removed", "removed_method": m,
                                     "removed_pocket": p, "removal_reason": removal_reasons[p],
                                     **summary})
        qc = {
            "version": VERSION, "protein": protein, "percentage_policy": graph.policy,
            "score_definition": "lexicographic (sum of diagonal percentages, sum of diagonal counts)",
            "baseline_definition": "common per-method input order after globally-zero filtering",
            "baseline_orders": baseline, "orders": order,
            "methods_without_input_pockets": [m for m in METHODS if not graph.ids[m]],
            "removed_globally_zero": {m: [p for p in removed[m] if graph.global_zero(p)] for m in METHODS},
            "removed_pockets": removed, "removal_reasons": removal_reasons,
            "drop_globally_below": args.drop_globally_below,
            "seed": seed, "seed_locked_during_swaps": seed is not None and not has_priority,
            "selected_priority_source": str(args.selected_priority) if args.selected_priority else None,
            "selected_priority_sha256": hashlib.sha256(args.selected_priority.read_bytes()).hexdigest() if args.selected_priority else None,
            "selected_prefix_orders": priority, "selected_records": selected_records,
            "diagonal_before": before, "diagonal_after_greedy": after_greedy,
            "diagonal_after": after, "swap_count": swaps, "swap_pass_limit": args.swap_passes,
            "baseline_improved_or_equal": after >= before,
            "optimization": "greedy with bounded local swaps; no global optimality guarantee",
            "missing_score_policy": "missing contributes no percentage/count to the objective; stays blank in CSV",
            "count_sources": sorted({e.source for e in graph.edges.values()}),
            "inputs": [{"path": str(p.resolve()), "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
                       for p in graph.sources],
            "outputs": outputs, "status": "complete",
        }
        (output / f"{protein}_overlap_order_qc.json").write_text(
            json.dumps(qc, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    except Exception:
        # Do not delete partial artifacts; no complete QC JSON is written on failure.
        raise
    print(f"{protein}: diagonal {before} -> {after}; swaps={swaps}; output={output}")
    if after < before:
        print(f"WARNING {protein}: seeded local heuristic scored below the input baseline; inspect QC.",
              file=sys.stderr)


def batch(args):
    root = args.batch_root
    suffix = "_relative.csv" if args.input_kind == "raw" else "_relative_simple.csv"
    dirs = [root / folder / "Intersection_Percentage" /
            ("Simplified" if args.input_kind == "simplified" else "") for folder in FOLDERS]
    for directory in dirs:
        if not directory.is_dir():
            raise ValueError(f"Input directory not found: {directory}")
    maps = []
    for directory in dirs:
        maps.append({p.name[:-len(suffix)]: p for p in directory.glob(f"*{suffix}")})
    proteins = sorted(set().union(*(set(mapping) for mapping in maps)))
    if args.protein:
        if args.protein not in proteins:
            raise ValueError(f"No matrices found for {args.protein}")
        proteins = [args.protein]
    if not proteins:
        raise ValueError("No matching protein CSVs found")
    incomplete = {p: [FOLDERS[i] for i, mapping in enumerate(maps) if p not in mapping]
                  for p in proteins if not all(p in mapping for mapping in maps)}
    if incomplete:
        raise ValueError("Batch preflight failed; missing comparison matrices: " +
                         json.dumps(incomplete, ensure_ascii=False))
    if args.out_dir.exists():
        raise ValueError(f"Batch output already exists: {args.out_dir}; choose a new directory")
    args.out_dir.mkdir(parents=True)
    errors = []
    with (args.out_dir / "batch_overlap_order_manifest.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["protein", "status", "output_directory", "message"])
        for protein in proteins:
            destination = args.out_dir / protein
            paths = [mapping[protein] for mapping in maps]
            counts = [args.counts_root / folder / f"{protein}_counts.csv"
                      if args.counts_root else None for folder in FOLDERS]
            try:
                process(protein, paths, counts, destination, args)
                writer.writerow([protein, "complete", str(destination), ""])
            except (OSError, ValueError, RuntimeError) as exc:
                errors.append(protein)
                writer.writerow([protein, "error", str(destination), str(exc)])
                print(f"ERROR {protein}: {exc}", file=sys.stderr)
            stream.flush()
    return 1 if errors else 0


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--as-fp", type=local_path, help="AutoSite/FPocket percentage matrix")
    p.add_argument("--as-dg", type=local_path, help="AutoSite/DoGSite percentage matrix")
    p.add_argument("--fp-dg", type=local_path, help="FPocket/DoGSite percentage matrix")
    p.add_argument("--as-fp-counts", type=local_path)
    p.add_argument("--as-dg-counts", type=local_path)
    p.add_argument("--fp-dg-counts", type=local_path)
    p.add_argument("--selected-priority", type=local_path,
                   help="Selected TSV: prioritize unique constituent pockets in file row order, with fixed per-method prefixes")
    p.add_argument("--protein", help="Protein identifier; optional batch filter")
    p.add_argument("--out-dir", type=local_path, required=True, help="New output directory")
    p.add_argument("--drop-globally-zero", action="store_true")
    p.add_argument("--drop-globally-below", type=float, metavar="PERCENT",
                   help="Remove pockets whose known percentages to BOTH other methods are all strictly below this threshold")
    p.add_argument("--percentage-policy", choices=("legacy-floor", "exact"), default="legacy-floor")
    p.add_argument("--swap-passes", type=int, default=30, help="Maximum improving swaps (default 30)")
    p.add_argument("--batch-root", type=local_path, help="Parent Fpockt directory containing the three comparisons")
    p.add_argument("--input-kind", choices=("raw", "simplified"), default="raw")
    p.add_argument("--counts-root", type=local_path, help="Batch counts: ROOT/comparison_folder/PROTEIN_counts.csv")
    return p


def main():
    p = parser()
    args = p.parse_args()
    if args.drop_globally_below is not None and (not math.isfinite(args.drop_globally_below) or not 0 <= args.drop_globally_below <= 100):
        p.error("--drop-globally-below must be finite and between 0 and 100")
    if args.swap_passes < 0:
        p.error("--swap-passes must be nonnegative")
    try:
        if args.batch_root:
            if any((args.as_fp, args.as_dg, args.fp_dg, args.as_fp_counts, args.as_dg_counts, args.fp_dg_counts)):
                p.error("Batch mode cannot be combined with individual matrix paths")
            return batch(args)
        if args.counts_root:
            p.error("--counts-root requires --batch-root")
        if not all((args.as_fp, args.as_dg, args.fp_dg, args.protein)):
            p.error("Single-protein mode requires --protein, --as-fp, --as-dg, --fp-dg")
        process(args.protein, (args.as_fp, args.as_dg, args.fp_dg),
                (args.as_fp_counts, args.as_dg_counts, args.fp_dg_counts), args.out_dir, args)
        return 0
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
