#!/usr/bin/env python3
"""Exact pairwise matrices from complete merged GenBank locations, plus selected-TSV QC.

Reads existing annotations only. No pocket detection. Python >=3.10; Biopython.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import re
import sys
from Bio import SeqIO
from Bio.SeqFeature import ExactPosition
from three_method_overlap_order import METHODS, PAIRS, FOLDERS, PREFIXES, local_path

VERSION = "2.1.0"


def read_pockets(path):
    record = SeqIO.read(str(path), "genbank")
    pockets = {m: {} for m in METHODS}
    ranges = []
    for feature in record.features:
        match = re.fullmatch(r"(ASpock|Fpockt|DGsite)_?(\d+)", feature.type, re.IGNORECASE)
        if not match:
            continue
        name = next(m for m in METHODS if PREFIXES[m].lower() == match[1].lower())
        pocket = PREFIXES[name] + f"{int(match[2]):02d}"
        if pocket in pockets[name]:
            raise ValueError(f"{path}: duplicate identifier {pocket}")
        if feature.location is None:
            raise ValueError(f"{path}: no parsed location for {pocket}")
        for part in feature.location.parts:
            if part.ref is not None or type(part.start) is not ExactPosition or type(part.end) is not ExactPosition:
                raise ValueError(f"{path}: remote/fuzzy location for {pocket}")
            if len(part) > 1:
                ranges.append(pocket)
        residues = {int(position) + 1 for position in feature.location}
        if not residues or min(residues) < 1 or max(residues) > len(record.seq):
            raise ValueError(f"{path}: empty/out-of-sequence location for {pocket}")
        pockets[name][pocket] = residues
    if not any(pockets.values()):
        raise ValueError(f"{path}: no recognized pocket features")
    return pockets, len(record.seq), sorted(set(ranges))


def selected_rows(path):
    result = {}
    if path is None:
        return result
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not {"Protein", "Pocket_A", "Pocket_B", "Pocket_C"}.issubset(reader.fieldnames or []):
            raise ValueError("Selected TSV requires Protein, Pocket_A, Pocket_B, Pocket_C")
        for line, row in enumerate(reader, 2):
            row["_line"] = line
            result.setdefault(row["Protein"].strip(), []).append(row)
    return result


def validate_selected(protein, pockets, rows):
    checks = []
    for row in rows:
        resolved = {}
        for name, column in (("FPocket", "Pocket_A"), ("DoGSite", "Pocket_B"), ("AutoSite", "Pocket_C")):
            match = re.fullmatch(re.escape(PREFIXES[name]) + r"_?(\d+)", row[column].strip(), re.IGNORECASE)
            if not match:
                raise ValueError(f"Selected TSV line {row['_line']}: invalid {column}")
            pocket = PREFIXES[name] + f"{int(match[1]):02d}"
            if pocket not in pockets[name]:
                raise ValueError(f"{protein}: selected {pocket} absent from GenBank")
            resolved[name] = pocket
        sets = {m: pockets[m][resolved[m]] for m in METHODS}
        shared = set.intersection(*(sets[m] for m in METHODS))
        observed = {
            "Nr_Fp": len(sets["FPocket"]), "Nr_DG": len(sets["DoGSite"]),
            "Nr_Asp": len(sets["AutoSite"]), "common residues": len(shared),
            "%Fpockt": len(shared) / len(sets["FPocket"]),
            "%DGsite": len(shared) / len(sets["DoGSite"]),
            "%ASpock": len(shared) / len(sets["AutoSite"]),
        }
        observed["max%"] = max(observed["%Fpockt"], observed["%DGsite"], observed["%ASpock"])
        for key, actual in observed.items():
            text = (row.get(key) or "").strip()
            if not text:
                continue
            expected = float(text)
            tolerance = 0.00051 if key.startswith("%") or key == "max%" else 0
            if not math.isfinite(expected) or not math.isclose(actual, expected, rel_tol=0, abs_tol=tolerance):
                raise ValueError(f"Selected line {row['_line']} {row.get('Pocket', '')}: "
                                 f"{key}={expected}, full GenBank implies {actual}")
        text = (row.get("Residue list") or row.get("list") or "").strip()
        if text and set(map(int, re.findall(r"\d+", text))) != shared:
            raise ValueError(f"Selected line {row['_line']}: triple residue list differs from GenBank")
        checks.append({"tsv_line": row["_line"], "selected_pocket": row.get("Pocket", ""),
                       "pockets": resolved, "common_residues": sorted(shared),
                       "sizes": {m: len(sets[m]) for m in METHODS}, "status": "matched"})
    return checks


def write_matrix(path, pockets, pair, metric):
    left, right = pair
    rows = sorted(pockets[left], key=lambda p: int(re.search(r"\d+$", p)[0]))
    cols = sorted(pockets[right], key=lambda p: int(re.search(r"\d+$", p)[0]))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["PocketID", "TotalResidues", *cols])
        writer.writerow(["Total Residues", "", *(len(pockets[right][p]) for p in cols)])
        for a in rows:
            values = []
            for b in cols:
                common = len(pockets[left][a] & pockets[right][b])
                value = common if metric == "counts" else max(
                    100 * (common / len(pockets[left][a])), 100 * (common / len(pockets[right][b])))
                values.append(value)
            writer.writerow([a, len(pockets[left][a]), *values])


def process(path, output, selected, selected_path):
    protein = path.stem
    pockets, length, ranges = read_pockets(path)
    checks = validate_selected(protein, pockets, selected.get(protein, []))
    files = []
    for pair, folder in zip(PAIRS, FOLDERS):
        for metric in ("percentage", "counts"):
            relative = (f"Intersection_Percentage/{protein}_relative.csv" if metric == "percentage"
                        else f"{protein}_counts.csv")
            target = output / folder / relative
            write_matrix(target, pockets, pair, metric)
            files.append(str(target.relative_to(output)))
    qc_dir = output / "QC" / protein
    qc_dir.mkdir(parents=True, exist_ok=True)
    with (qc_dir / f"{protein}_complete_pocket_residues.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["method", "pocket_id", "residue_count", "residues"])
        for name in METHODS:
            for pocket, residues in pockets[name].items():
                writer.writerow([name, pocket, len(residues), ",".join(map(str, sorted(residues)))])
    qc = {
        "version": VERSION, "protein": protein, "status": "complete",
        "source": str(path), "source_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "sequence_length": length, "pocket_counts": {m: len(pockets[m]) for m in METHODS},
        "methods_without_pockets": [m for m in METHODS if not pockets[m]],
        "range_features_expanded": ranges,
        "percentage_definition": "max(100*common/left_size,100*common/right_size), no truncation",
        "counts_definition": "exact intersection of all positions in complete feature locations",
        "selected_tsv": str(selected_path) if selected_path else None,
        "selected_tsv_sha256": hashlib.sha256(selected_path.read_bytes()).hexdigest() if selected_path else None,
        "selected_checks": checks, "outputs": files,
        "legacy_triple_range_caveat": "Legacy triple builder uses each part.start only. "
                                     "Here ranges expand fully; selected validation detects differences.",
    }
    (qc_dir / f"{protein}_complete_pairwise_qc.json").write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    print(f"{protein}: complete counts {qc['pocket_counts']}; matched selected rows={len(checks)}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genbank-root", type=local_path, required=True)
    parser.add_argument("--out-dir", type=local_path, required=True)
    parser.add_argument("--protein")
    parser.add_argument("--selected-tsv", type=local_path)
    args = parser.parse_args()
    try:
        paths = sorted(args.genbank_root.glob("*.gb"))
        if args.protein:
            paths = [p for p in paths if p.stem == args.protein]
        if not paths:
            raise ValueError("No matching merged GenBank files")
        selected = selected_rows(args.selected_tsv)
        if not args.protein:
            absent = set(selected) - {p.stem for p in paths}
            if absent:
                raise ValueError(f"Selected proteins missing from GenBank input: {sorted(absent)}")
        if args.out_dir.exists():
            raise ValueError(f"Output exists: {args.out_dir}; choose a new directory")
        args.out_dir.mkdir(parents=True)
        failures = 0
        with (args.out_dir / "batch_complete_pairwise_manifest.tsv").open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t")
            writer.writerow(["protein", "status", "message"])
            for path in paths:
                try:
                    process(path, args.out_dir, selected, args.selected_tsv)
                    writer.writerow([path.stem, "complete", ""])
                except (OSError, ValueError, RuntimeError) as exc:
                    failures += 1
                    writer.writerow([path.stem, "error", str(exc)])
                    print(f"ERROR {path.stem}: {exc}", file=sys.stderr)
                stream.flush()
        return 1 if failures else 0
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
