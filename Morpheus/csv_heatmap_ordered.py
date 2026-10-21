#!/usr/bin/env python3
"""Render v2 ordered CSVs exactly as stored, with vector PDF/SVG and optional PNG.

Single file: --input CSV --output-prefix PATH --mode percentage|counts.
Batch: --batch-dir ORDERED_ROOT --out-dir NEW_ROOT --mode percentage|counts|both.
Batch --composite adds structure + three CSV-based heatmap panels.
Python >= 3.10; numpy, matplotlib and Pillow. No clustering, sorting, or zero filling.
"""
from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import json
import math
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as colors
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
from PIL import Image



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
MISSING = {"", "-", "na", "nan", "none", "null"}


def load_matrix(path, mode):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        records = [r for r in csv.reader(stream) if any(c.strip() for c in r)]
    if not records:
        raise ValueError(f"{path}: empty CSV")
    cols = [c.strip() for c in records[0][1:]]
    if any(not c for c in cols) or len(cols) != len(set(cols)):
        raise ValueError(f"{path}: empty/duplicate column labels")
    rows, data = [], []
    if any(c.replace(" ", "").replace("_", "").lower() == "totalresidues" for c in cols):
        raise ValueError("Renderer expects v2 pocket-only matrices; residue totals are metadata")
    for line, record in enumerate(records[1:], 2):
        if len(record) != len(cols) + 1:
            raise ValueError(f"{path}:{line}: inconsistent CSV width")
        pocket = record[0].strip()
        if not pocket or pocket in rows:
            raise ValueError(f"{path}:{line}: empty/duplicate pocket identifier")
        if pocket.replace(" ", "").replace("_", "").lower() == "totalresidues":
            raise ValueError("Renderer expects v2 output without a Total Residues row")
        rows.append(pocket)
        values = []
        for text in record[1:]:
            if text.strip().lower() in MISSING:
                values.append(float("nan"))
                continue
            try:
                value = float(text)
            except ValueError as exc:
                raise ValueError(f"{path}:{line}: invalid numeric value {text!r}") from exc
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"{path}:{line}: negative or nonfinite value")
            if mode == "percentage" and value > 100:
                raise ValueError(f"{path}:{line}: percentage exceeds 100")
            if mode == "counts" and not value.is_integer():
                raise ValueError(f"{path}:{line}: noninteger residue count")
            values.append(value)
        data.append(values)
    return rows, cols, np.asarray(data, dtype=float).reshape((len(rows), len(cols)))


def colormap(config_path):
    # Default is the original Morpheus scheme, including the intentional 50% jump.
    config_path = config_path or Path(__file__).with_name("morpheus_pocket_overlap_legacy.json")
    config = json.loads(config_path.read_text(encoding="utf-8-sig"))
    missing = "#c0c0c0"
    if "valueToColorScheme" in config:
        scheme = next(iter(config["valueToColorScheme"].values()))
        fractions = np.asarray(scheme["fractions"], dtype=float)
        color_list = scheme["colors"]
        if len(fractions) != len(color_list) or len(fractions) < 2:
            raise ValueError("Invalid Morpheus color fractions/colors")
        if not np.isfinite(fractions).all() or fractions.max() <= 0:
            raise ValueError("Invalid Morpheus color fractions")
        if fractions.max() > 1:
            fractions = fractions / fractions.max()
        if fractions[0] != 0 or fractions[-1] != 1 or np.any(np.diff(fractions) < 0):
            raise ValueError("Morpheus fractions must be nondecreasing from 0 to 1")
        # A repeated position encodes left/right colors of a discontinuity.
        groups = []
        for position, color in zip(fractions, color_list):
            rgb = colors.to_rgb(color)
            if groups and position == groups[-1][0]:
                groups[-1][2] = rgb
            else:
                groups.append([float(position), rgb, rgb])
        channels = {
            name: [(p, left[i], right[i]) for p, left, right in groups]
            for i, name in enumerate(("red", "green", "blue"))
        }
        cmap = colors.LinearSegmentedColormap("morpheus_ordered", channels, N=256)
        missing = scheme.get("missingColor", missing)
    else:
        cmap = plt.get_cmap(config.get("colormap", "Blues")).copy()
    cmap.set_bad(missing)
    return cmap


def method(label):
    for prefix, name in (("ASpock", "AutoSite"), ("Fpockt", "FPocket"), ("DGsite", "DoGSite")):
        if label.lower().startswith(prefix.lower()):
            return name
    return "Pocket"


def draw_heatmap(fig, ax, rows, cols, data, mode, args, title, composite=False):
    large_panel = composite and data.size > 4000
    cell_fontsize = max(1.5, min(7, 350 / max(len(rows), len(cols), 1))) if large_panel else 7
    finite = data[np.isfinite(data)]
    vmax = 100.0 if mode == "percentage" else (
        args.count_vmax if args.count_vmax is not None
        else max(1.0, float(finite.max()) if finite.size else 1.0)
    )
    if finite.size and finite.max() > vmax:
        raise ValueError(f"Data exceeds --count-vmax {vmax}")
    cmap = colormap(args.color_config)
    missing_color = colors.to_hex(cmap(np.nan))
    if rows and cols:
        mesh = ax.pcolormesh(
            np.ma.masked_invalid(data), cmap=cmap, vmin=0, vmax=vmax,
            edgecolors="#202020", linewidth=0.35, shading="flat", rasterized=False)
        ax.set_xlim(0, len(cols))
        ax.set_ylim(len(rows), 0)
        ax.set_aspect("equal")
        ax.set_xticks(np.arange(len(cols)) + 0.5, labels=cols, rotation=90)
        ax.set_yticks(np.arange(len(rows)) + 0.5, labels=rows)
        ax.tick_params(length=0)
        if large_panel:
            ax.tick_params(labelsize=max(3, min(9, 500 / max(len(rows), len(cols), 1))))
        for i, j in zip(*np.where(~np.isfinite(data))):
            ax.add_patch(Rectangle((j, i), 1, 1, facecolor=missing_color,
                                   hatch="///", edgecolor="#909090", linewidth=0.35))
            if not args.no_annotate:
                ax.text(j + .5, i + .5, "NA", ha="center", va="center", fontsize=6)
        if not args.no_annotate:
            annotation_mask = np.isfinite(data) & (data != 0) if large_panel else np.isfinite(data)
            for i, j in zip(*np.where(annotation_mask)):
                value = data[i, j]
                rgba = cmap(value / vmax)
                luminance = .2126 * rgba[0] + .7152 * rgba[1] + .0722 * rgba[2]
                text = (str(math.ceil(round(float(value), 10)))
                        if mode == "percentage" else f"{value:g}")
                ax.text(j + .5, i + .5, text, ha="center", va="center",
                        fontsize=cell_fontsize, color="white" if luminance < .45 else "black")
        bar = fig.colorbar(mesh, ax=ax, fraction=.035, pad=.025)
        if bar.solids is not None:
            bar.solids.set_rasterized(False)
        bar.set_label("Overlap (%)" if mode == "percentage" else "Shared residues (count)")
    else:
        ax.text(.5, .5, "No retained pockets for this comparison",
                transform=ax.transAxes, ha="center", va="center")
        ax.set_xticks([])
        ax.set_yticks([])
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position("top")
    ax.set_xlabel(args.xlabel or (method(cols[0]) if cols else "Pocket"))
    ax.set_ylabel(args.ylabel or (method(rows[0]) if rows else "Pocket"))
    ax.set_title(title, pad=18)
    return vmax




def raster_dpi(width, height, args):
    return min(args.dpi, max(1, math.floor(math.sqrt(
        args.max_raster_megapixels * 1_000_000 / (1.2 * width * height)))))


def save_figures(fig, files, args):
    """Bound raster memory, including embedded structure images in PDF/SVG."""
    width, height = fig.get_size_inches()
    safe_dpi = raster_dpi(width, height, args)
    if safe_dpi < args.dpi:
        print(f"Large figure: raster dpi capped at {safe_dpi} (requested {args.dpi}); vector cells unchanged")
    for path in files:
        fig.savefig(path, format=path.suffix[1:], dpi=safe_dpi, bbox_inches="tight")
    return safe_dpi


def resume_ready(qc_path, files, expected, resume):
    """Skip only matching completed outputs; incomplete generated files may be rebuilt."""
    if not resume or not qc_path.exists():
        return False
    try:
        saved = json.loads(qc_path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    if saved.get("status") != "complete":
        return False
    for key, value in expected.items():
        if key == "color_config":
            # Older runs used Intersection_graphs/../Morpheus; it is the same file.
            if saved.get(key) is not None and local_path(saved[key]).resolve() == local_path(value).resolve():
                continue
        if key == "color_config_sha256" and saved.get(key) is None:
            # Backward compatibility for QC created before color hashes were recorded.
            continue
        if saved.get(key) != value:
            raise ValueError(f"Cannot resume {qc_path}: QC field {key} differs; choose a new output directory")
    if not all(path.is_file() and path.stat().st_size > 0 for path in files):
        return False
    print(f"SKIP completed: {qc_path.stem}")
    return True


def render(input_path, output_prefix, mode, args, title=None):
    rows, cols, data = load_matrix(input_path, mode)
    finite = data[np.isfinite(data)]
    vmax = 100.0 if mode == "percentage" else (
        args.count_vmax if args.count_vmax is not None else max(1.0, float(finite.max()) if finite.size else 1.0)
    )
    if finite.size and finite.max() > vmax:
        raise ValueError(f"{input_path}: data exceeds --count-vmax {vmax}; choose a larger range")
    files = [Path(str(output_prefix) + "." + fmt) for fmt in args.formats]
    qc_path = Path(str(output_prefix) + "_render_qc.json")
    expected = {
        "input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "mode": mode, "vmin": 0, "vmax": vmax, "rows": rows, "columns": cols,
        "percentage_label_rounding": "ceil, after 10-decimal floating-point stabilization",
        "color_config": str((args.color_config or Path(__file__).with_name("morpheus_pocket_overlap_legacy.json")).resolve()),
        "color_config_sha256": hashlib.sha256((args.color_config or Path(__file__).with_name("morpheus_pocket_overlap_legacy.json")).read_bytes()).hexdigest(),
    }
    if resume_ready(qc_path, files, expected, args.resume):
        return False
    if not args.resume and any(path.exists() for path in [*files, qc_path]):
        raise ValueError(f"Heatmap output already exists: {output_prefix}")
    # Remove the completion marker before rebuilding incomplete generated outputs.
    if args.resume and qc_path.exists():
        qc_path.unlink()
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    width = max(5.0, 3.0 + args.cell_inches * len(cols))
    height = max(4.0, 2.8 + args.cell_inches * len(rows))
    plt.rcParams.update({"pdf.fonttype": 42, "svg.fonttype": "none", "font.size": 9})
    fig, ax = plt.subplots(figsize=(width, height), dpi=min(100, raster_dpi(width, height, args)))
    try:
        vmax = draw_heatmap(fig, ax, rows, cols, data, mode, args,
                            title or args.title or input_path.stem)
        fig.tight_layout()
        used_raster_dpi = save_figures(fig, files, args)
    finally:
        plt.close(fig)
        fig = ax = None
        gc.collect()
    qc = {
        "version": VERSION, "color_config": str((args.color_config or Path(__file__).with_name("morpheus_pocket_overlap_legacy.json")).resolve()),
        "color_config_sha256": hashlib.sha256((args.color_config or Path(__file__).with_name("morpheus_pocket_overlap_legacy.json")).read_bytes()).hexdigest(),
        "input": str(input_path.resolve()),
        "input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "mode": mode, "vmin": 0, "vmax": vmax, "rows": rows, "columns": cols,
        "missing_cells": int(np.isnan(data).sum()), "zero_cells": int((data == 0).sum()),
        "data_order": "exact CSV order", "pdf_svg_cells": "vector",
        "percentage_label_rounding": "ceil, after 10-decimal floating-point stabilization",
        "color_values": "unrounded CSV values",
        "raster_dpi_used": used_raster_dpi, "dpi_requested": args.dpi,
        "outputs": [str(path) for path in files], "status": "complete",
    }
    qc_path.write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    print(f"Rendered {input_path.name}: " + ", ".join(str(p) for p in files))
    return True



def render_composite(protein, folder, mode, output_dir, args):
    pairs = (("AutoSite", "FPocket"), ("AutoSite", "DoGSite"), ("FPocket", "DoGSite"))
    inputs = [folder / f"{protein}_{a}_{b}_ordered_{mode}.csv" for a, b in pairs]
    matrices = [load_matrix(path, mode) for path in inputs]
    # Verify common orders before composition, without sorting any axes.
    if (matrices[0][0] != matrices[1][0] or matrices[0][1] != matrices[2][0]
            or matrices[1][1] != matrices[2][1]):
        raise ValueError(f"{protein}: inconsistent method ordering across panels")
    source = args.structure_image or args.structure_reference
    if source is None and args.structure_reference_dir:
        source = args.structure_reference_dir / f"{protein}_structure_and_pairwise_pocket_overlap_comparisons.png"
    if source is None:
        raise ValueError("Composite needs --structure-image, --structure-reference or --structure-reference-dir")
    if not source.is_file():
        raise ValueError(f"Structure reference not found: {source}")
    with Image.open(source) as image:
        image = image.convert("RGB")
        if args.structure_image:
            crop_box = None
        else:
            x0, y0, x1, y1 = args.structure_crop
            w, h = image.size
            crop_box = (round(x0 * w), round(y0 * h), round(x1 * w), round(y1 * h))
            image = image.crop(crop_box)
        structure = np.asarray(image).copy()
    prefix = output_dir / protein / f"{protein}_structure_and_pairwise_pocket_overlap_comparisons_ordered_{mode}"
    files = [Path(str(prefix) + "." + fmt) for fmt in args.formats]
    qc_path = Path(str(prefix) + "_composite_qc.json")
    expected = {
        "protein": protein, "mode": mode,
        "structure_source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "structure_crop_pixels": list(crop_box) if crop_box else None,
        "percentage_label_rounding": "ceil, after 10-decimal floating-point stabilization",
        "color_config": str((args.color_config or Path(__file__).with_name("morpheus_pocket_overlap_legacy.json")).resolve()),
        "color_config_sha256": hashlib.sha256((args.color_config or Path(__file__).with_name("morpheus_pocket_overlap_legacy.json")).read_bytes()).hexdigest(),
        "panel_layout": {"A": "structure top left", "B": "AutoSite_DoGSite top right",
                         "C": "FPocket_DoGSite bottom left", "D": "AutoSite_FPocket bottom right"},
    }
    if args.resume and qc_path.exists():
        try:
            saved = json.loads(qc_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            saved = {}
        previous = {p.get("input"): p.get("sha256") for p in saved.get("panels", [])}
        if saved.get("status") == "complete" and previous != {
                str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}:
            raise ValueError(f"Cannot resume {qc_path}: composite CSV inputs changed")
    if resume_ready(qc_path, files, expected, args.resume):
        return False
    if not args.resume and any(path.exists() for path in [*files, qc_path]):
        raise ValueError(f"Composite output already exists: {prefix}")
    if args.resume and qc_path.exists():
        qc_path.unlink()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    width = min(32, max(18, 6 + .65 * max(len(cols) for _, cols, _ in matrices)))
    height = min(24, max(14, 4 + .50 * max(len(rows) for rows, _, _ in matrices)))
    plt.rcParams.update({"pdf.fonttype": 42, "svg.fonttype": "none", "font.size": 9})
    fig, axes = plt.subplots(2, 2, figsize=(width, height),
                             dpi=min(100, raster_dpi(width, height, args)), layout="constrained")
    try:
        axes[0, 0].imshow(structure)
        axes[0, 0].set_axis_off()
        axes[0, 0].set_title("A   Filtered AlphaFold model used for pocket prediction",
                           fontweight="bold", pad=18)
        titles = ("B   AutoSite–DoGSite pocket overlap",
                  "C   FPocket–DoGSite pocket overlap",
                  "D   AutoSite–FPocket pocket overlap")
        panels = []
        for ax, title, index in zip(
                (axes[0, 1], axes[1, 0], axes[1, 1]), titles, (1, 2, 0)):
            path = inputs[index]
            rows, cols, data = matrices[index]
            vmax = draw_heatmap(fig, ax, rows, cols, data, mode, args, title, composite=True)
            panels.append({"panel": title[0], "comparison": "_".join(pairs[index]),
                           "cell_annotations": "nonzero values only" if data.size > 4000 else "all values",
                           "input": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                           "rows": rows, "columns": cols, "vmin": 0, "vmax": vmax})
        label = "percentage overlaps" if mode == "percentage" else "shared-residue counts"
        fig.suptitle(f"{protein} structural model and pairwise pocket comparisons — {label}",
                     fontsize=18, fontweight="bold")
        used_raster_dpi = save_figures(fig, files, args)
    finally:
        plt.close(fig)
        fig = ax = axes = None
        gc.collect()
    qc = {"version": VERSION, "protein": protein, "mode": mode, "status": "complete",
          "structure_source": str(source), "structure_source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
          "structure_crop_pixels": crop_box, "figure_inches": [width, height],
          "large_composite_policy": "bounded to 32x24 inches; above 4000 cells omit only zero text, retain nonzero values and all data",
          "structure_panel": "raster from supplied existing image",
          "heatmap_panels": "vector, drawn directly from ordered CSVs", "panels": panels,
          "panel_layout": {"A": "structure top left", "B": "AutoSite_DoGSite top right",
                           "C": "FPocket_DoGSite bottom left", "D": "AutoSite_FPocket bottom right"},
          "percentage_label_rounding": "ceil, after 10-decimal floating-point stabilization",
          "color_values": "unrounded CSV values",
          "color_config": str((args.color_config or Path(__file__).with_name("morpheus_pocket_overlap_legacy.json")).resolve()),
        "color_config_sha256": hashlib.sha256((args.color_config or Path(__file__).with_name("morpheus_pocket_overlap_legacy.json")).read_bytes()).hexdigest(),
          "raster_dpi_used": used_raster_dpi, "dpi_requested": args.dpi,
          "outputs": [str(path) for path in files]}
    qc_path.write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    print(f"Composite {protein} ({mode}): " + ", ".join(str(p) for p in files))
    return True


def batch(args):
    qcs = sorted(args.batch_dir.glob("**/*_overlap_order_qc.json"))
    if not qcs:
        raise ValueError(f"No completed ordering QC files found in {args.batch_dir}")
    jobs = []
    seen = set()
    for path in qcs:
        qc = json.loads(path.read_text(encoding="utf-8"))
        if qc.get("status") != "complete" or qc.get("version") != VERSION:
            raise ValueError(f"Unsupported or incomplete ordering QC: {path}")
        protein = qc["protein"]
        if not isinstance(protein, str) or not re.fullmatch(r"[A-Za-z0-9_.-]+", protein) or protein in (".", ".."):
            raise ValueError(f"Unsafe protein identifier in {path}")
        if args.protein and protein != args.protein:
            continue
        if protein in seen:
            raise ValueError(f"Duplicate protein in ordering QC: {protein}")
        seen.add(protein)
        # Output names come from QC, never sort pocket axes.
        for name in qc["outputs"]:
            if Path(name).name != name:
                raise ValueError(f"Unsafe output filename in {path}")
            mode = "percentage" if name.endswith("_ordered_percentage.csv") else (
                "counts" if name.endswith("_ordered_counts.csv") else None)
            if mode is None or args.mode not in ("both", mode):
                continue
            source = path.parent / name
            if not source.is_file():
                raise ValueError(f"Ordered input not found: {source}")
            prefix = args.out_dir / protein / (source.stem + "_heatmap")
            jobs.append((protein, source, prefix, mode))
    if not jobs:
        raise ValueError("No matrices match the requested --mode")
    if args.composite:
        if (args.structure_image or args.structure_reference) and len(seen) != 1:
            raise ValueError("A single structure image/reference requires one protein; use --protein")
        for protein in sorted(seen):
            source = args.structure_image or args.structure_reference
            if source is None and args.structure_reference_dir:
                source = args.structure_reference_dir / f"{protein}_structure_and_pairwise_pocket_overlap_comparisons.png"
            if source is None or not source.is_file():
                raise ValueError(f"Missing structure source for {protein}")
    if args.out_dir.exists() and not args.resume:
        raise ValueError(f"Heatmap output already exists: {args.out_dir}; use --resume or a new directory")
    if args.resume and not args.out_dir.is_dir():
        raise ValueError("--resume requires an existing output directory")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    failures = 0
    manifest = args.out_dir / "batch_heatmap_manifest.tsv"
    existing_manifest = manifest.is_file() and manifest.stat().st_size > 0
    with manifest.open("a" if args.resume else "w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        if not args.resume or not existing_manifest:
            writer.writerow(["protein", "input_csv", "mode", "status", "message"])
        for protein, source, prefix, mode in jobs:
            try:
                created = render(source, prefix, mode, args, title=f"{protein} - {source.stem[len(protein)+1:].replace('_ordered_', ' - ').replace('_', ' ')}")
                if created:
                    writer.writerow([protein, str(source), mode, "complete", ""])
            except (OSError, ValueError, RuntimeError) as exc:
                failures += 1
                writer.writerow([protein, str(source), mode, "error", str(exc)])
                print(f"ERROR: {exc}", file=sys.stderr)
            stream.flush()
        if args.composite:
            modes = ("percentage", "counts") if args.mode == "both" else (args.mode,)
            for protein in sorted(seen):
                folder = next(source.parent for p, source, _, _ in jobs if p == protein)
                for mode in modes:
                    try:
                        created = render_composite(protein, folder, mode, args.out_dir, args)
                        if created:
                            writer.writerow([protein, "four-panel composite", mode, "complete", ""])
                    except (OSError, ValueError, RuntimeError) as exc:
                        failures += 1
                        writer.writerow([protein, "four-panel composite", mode, "error", str(exc)])
                        print(f"ERROR: {exc}", file=sys.stderr)
                    stream.flush()
    return 1 if failures else 0


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--max-raster-megapixels", type=float, default=20,
                   help="Raster memory budget per figure, including embedded images (default 20 MP); vector cells retain full detail")
    p.add_argument("--resume", action="store_true", help="Skip verified completed outputs and rebuild missing/partial figures in the existing directory")
    p.add_argument("--input", type=local_path)
    p.add_argument("--output-prefix", type=local_path, help="Output basename without extension")
    p.add_argument("--batch-dir", type=local_path)
    p.add_argument("--out-dir", type=local_path, help="New heatmap directory for batch mode")
    p.add_argument("--protein", help="Select one protein from a batch directory")
    p.add_argument("--composite", action="store_true", help="Also draw structure + three overlap panels")
    structures = p.add_mutually_exclusive_group()
    structures.add_argument("--structure-image", type=local_path, help="Standalone structure image")
    structures.add_argument("--structure-reference", type=local_path, help="Existing four-panel figure; reuse panel A")
    structures.add_argument("--structure-reference-dir", type=local_path, help="Directory of per-protein legacy four-panel PNGs")
    p.add_argument("--structure-crop", nargs=4, type=float, default=[.055, .145, .475, .495],
                   metavar=("X0", "Y0", "X1", "Y1"), help="Normalized reference-image crop for panel A")
    p.add_argument("--mode", choices=("percentage", "counts", "both"), default="percentage")
    p.add_argument("--formats", nargs="+", choices=("pdf", "svg", "png"), default=["pdf", "svg"])
    p.add_argument("--color-config", type=local_path, help="Optional existing Morpheus JSON color scheme")
    p.add_argument("--count-vmax", type=float, help="Fixed maximum for count plots; never clips silently")
    p.add_argument("--cell-inches", type=float, default=.38)
    p.add_argument("--dpi", type=int, default=300, help="PNG resolution; PDF/SVG cells remain vector")
    p.add_argument("--no-annotate", action="store_true")
    p.add_argument("--title")
    p.add_argument("--xlabel")
    p.add_argument("--ylabel")
    args = p.parse_args()
    if not math.isfinite(args.max_raster_megapixels) or args.max_raster_megapixels <= 0:
        p.error("--max-raster-megapixels must be finite and positive")
    if args.cell_inches <= 0 or not math.isfinite(args.cell_inches) or args.dpi <= 0:
        p.error("--cell-inches and --dpi must be positive and finite")
    if args.count_vmax is not None and (args.count_vmax <= 0 or not math.isfinite(args.count_vmax)):
        p.error("--count-vmax must be positive and finite")
    if len(args.formats) != len(set(args.formats)):
        p.error("Duplicate --formats")
    x0, y0, x1, y1 = args.structure_crop
    if not (0 <= x0 < x1 <= 1 and 0 <= y0 < y1 <= 1):
        p.error("--structure-crop must satisfy 0 <= X0 < X1 <= 1 and 0 <= Y0 < Y1 <= 1")
    if args.composite and not args.batch_dir:
        p.error("--composite requires --batch-dir")
    if (args.structure_image or args.structure_reference or args.structure_reference_dir) and not args.composite:
        p.error("Structure options require --composite")
    if args.protein and not args.batch_dir:
        p.error("--protein requires --batch-dir")
    try:
        if args.batch_dir:
            if not args.out_dir or args.input or args.output_prefix:
                p.error("Batch mode requires --out-dir and cannot use single-file arguments")
            return batch(args)
        if not args.input or not args.output_prefix or args.out_dir or args.mode == "both":
            p.error("Single-file mode requires --input, --output-prefix, --mode percentage|counts")
        render(args.input, args.output_prefix, args.mode, args)
        return 0
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
