#!/usr/bin/env python3
"""Validate the complete artifact contract for one synthetic container run."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


TILE_COLUMNS = [
    "file", "sensor", "country", "tile_id", "tile_xmin", "tile_ymin",
    "edge_tile", "n_tiles_plot", "vox_filled", "vox_total", "metric_height_upper_m", "veg_density",
    "zsd", "zskew", "zkurt", "zq90", "box_dim_fixed", "vci", "rumple",
    "gap_fraction", "chm_sd", "chm_cv", "height_max", "height_mean",
    "n_seg_total", "n_trees", "tree_height_max", "tree_height_mean",
    "tree_height_gini", "tree_crownarea_mean", "tree_crownarea_max",
    "tree_crownarea_gini", "tree_volume_mean", "tree_volume_max",
    "tree_volume_gini",
]
DIAGNOSTIC_COLUMNS = [
    "file", "sensor", "country", "tile_id", "PredInstance", "n_vox",
    "apex_z", "pc_ext1", "pc_ext2", "pc_ext3", "n_zlayer", "pass_vox",
    "pass_apex", "pass_thick", "pass_zlayer", "is_tree", "apex_in_tile",
]
TREE_COLUMNS = {
    "n_seg_total", "n_trees", "tree_height_max", "tree_height_mean",
    "tree_height_gini", "tree_crownarea_mean", "tree_crownarea_max",
    "tree_crownarea_gini", "tree_volume_mean", "tree_volume_max",
    "tree_volume_gini",
}
EXPECTED_TREE_COUNTS = {
    "PredInstance": 2,
    "PredInstance_SAT": 1,
    "TreeAlias": 1,
}


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames or [], list(reader)


def verify_geojson(
    path: Path,
    rows: list[dict[str, str]],
    expected_tiles: int,
) -> None:
    with path.open(encoding="utf-8") as handle:
        features = json.load(handle).get("features", [])
    if len(features) != expected_tiles:
        raise SystemExit(f"{path}: feature count differs")

    for row, feature in zip(rows, features):
        properties = feature["properties"]
        for column, csv_value in row.items():
            geojson_value = properties.get(column)
            if csv_value == "NA":
                matches = geojson_value is None
            elif csv_value in {"TRUE", "FALSE"}:
                matches = geojson_value is (csv_value == "TRUE")
            elif isinstance(geojson_value, (int, float)):
                matches = abs(float(csv_value) - float(geojson_value)) <= 1e-9
            else:
                matches = str(geojson_value) == csv_value
            if not matches:
                raise SystemExit(f"{path}: {column} differs")


def verify_dimension(
    output_dir: Path,
    dataset_id: str,
    dimension: str | None,
    expected_tiles: int,
) -> Path:
    suffix = "" if dimension is None else f"_{dimension}"
    result_path = output_dir / f"{dataset_id}{suffix}_results.csv"
    diagnostic_path = output_dir / (
        f"{dataset_id}{suffix}_segment_diagnostics.csv"
    )
    geojson_path = output_dir / f"{dataset_id}{suffix}_tiles.geojson"

    fields, rows = read_csv(result_path)
    if fields != TILE_COLUMNS:
        raise SystemExit(f"tile CSV schema differs from current contract: {result_path}")
    if len(rows) != expected_tiles:
        raise SystemExit(
            f"{result_path}: expected {expected_tiles} rows, got {len(rows)}"
        )
    if rows:
        height_domains = {float(row["metric_height_upper_m"]) for row in rows}
        if len(height_domains) != 1:
            raise SystemExit(f"{result_path}: tiles do not share one height domain")
        height_domain = next(iter(height_domains))
        if not math.isfinite(height_domain) or height_domain < 1:
            raise SystemExit(f"{result_path}: invalid metric height domain")
        if any(float(row["height_max"]) > height_domain for row in rows if row["height_max"] != "NA"):
            raise SystemExit(f"{result_path}: metric height domain is below observed canopy")
    expected_ids = [str(index) for index in range(1, expected_tiles + 1)]
    if [row["tile_id"] for row in rows] != expected_ids:
        raise SystemExit(f"{result_path}: tile IDs are not deterministic")

    expected_basename = (
        "point_cloud_segmented.laz" if dimension else "point_cloud.laz"
    )
    if any(row["file"] != expected_basename for row in rows):
        raise SystemExit(f"{result_path}: wrong point-cloud basename")
    if any(row["sensor"] != "NA" or row["country"] != "NA" for row in rows):
        raise SystemExit(f"{result_path}: automatic metadata must be NA")
    if any(row["n_tiles_plot"] != str(expected_tiles) for row in rows):
        raise SystemExit(f"{result_path}: n_tiles_plot differs")

    if dimension is None:
        if any(row[column] != "NA" for row in rows for column in TREE_COLUMNS):
            raise SystemExit("tree fields must be NA without an instance dimension")
    elif sum(int(row["n_trees"]) for row in rows) != EXPECTED_TREE_COUNTS[dimension]:
        raise SystemExit(f"{dimension}: trees were not counted exactly once")

    diagnostic_fields, diagnostics = read_csv(diagnostic_path)
    if diagnostic_fields != DIAGNOSTIC_COLUMNS:
        raise SystemExit(f"segment schema differs: {diagnostic_path}")
    if dimension is None and diagnostics:
        raise SystemExit("missing dimensions must produce header-only diagnostics")
    if dimension is not None:
        if not diagnostics or any(not row["PredInstance"] for row in diagnostics):
            raise SystemExit(f"{dimension}: missing segment IDs")
        if any(row["file"] != "point_cloud_segmented.laz" for row in diagnostics):
            raise SystemExit(f"{dimension}: wrong segment source basename")

    verify_geojson(geojson_path, rows, expected_tiles)
    if output_dir.name == "inclusion_only":
        actual_rumple = [float(row["rumple"]) for row in rows]
        expected_rumple = [12.1396, 12.0178, 12.2369, 12.72]
        if actual_rumple != expected_rumple:
            raise SystemExit(
                "tile rumple values differ from direct-TIN regression: "
                f"expected={expected_rumple}, actual={actual_rumple}"
            )
    return result_path


def verify_performance(
    output_dir: Path,
    dataset_id: str,
    dimensions: list[str],
    footprint_source: str,
) -> None:
    path = output_dir / f"{dataset_id}_performance.csv"
    expected = output_dir.name in {"aliased", "whole_cloud"}
    if path.exists() != expected:
        raise SystemExit("performance-report opt-in contract was not honored")
    if not expected:
        return

    _, rows = read_csv(path)
    if len(rows) != 1:
        raise SystemExit("performance report must contain exactly one row")
    performance = rows[0]
    required = {
        "footprint_source", "point_count", "tile_count", "instance_dimension",
        "peak_rss_mib", "grid_seconds", "dtm_seconds", "chm_seconds",
        "segment_seconds", "tile_seconds", "output_seconds", "total_seconds",
        "threads_requested", "threads_effective", "catalog_workers",
        "dtm_chunk_size", "chunk_size", "dtm_buffer", "dtm_storage_mode",
        "chm_storage_mode", "metric_height_upper_m", "metric_height_policy", "maximum_height", "observed_maximum_voxel_height_m", "observed_maximum_voxel_height_m_exact",
    }
    if required.difference(performance):
        raise SystemExit("performance report is missing required measurements")
    if performance["metric_height_policy"] != "observed_voxel_height_v1" or float(performance["maximum_height"]) != math.inf:
        raise SystemExit("performance report has the wrong normalized-height policy")
    exact_maximum = float(performance["observed_maximum_voxel_height_m_exact"])
    if float(performance["metric_height_upper_m"]) != max(1, math.ceil(exact_maximum)):
        raise SystemExit("performance height domain differs from the exact observed maximum")
    for result_path in output_dir.glob(f"{dataset_id}*_results.csv"):
        _, tile_rows = read_csv(result_path)
        if any(float(row["metric_height_upper_m"]) != float(performance["metric_height_upper_m"]) for row in tile_rows):
            raise SystemExit("performance and tile height domains differ")
    expected_dimensions = "NA" if not dimensions else "|".join(dimensions)
    if performance["instance_dimension"] != expected_dimensions:
        raise SystemExit("performance report has the wrong dimensions")
    if performance["footprint_source"] != footprint_source:
        raise SystemExit("performance report has the wrong footprint source")
    if performance["dtm_chunk_size"] != "60" or performance["chunk_size"] != "60":
        raise SystemExit("performance report has the wrong chunk defaults")
    if performance["dtm_buffer"] != "20":
        raise SystemExit("performance report has the wrong DTM buffer")
    expected_workers = min(2, int(performance["threads_effective"]))
    if int(performance["catalog_workers"]) != expected_workers:
        raise SystemExit("performance report has the wrong worker count")
    for artifact in ("dtm_storage_mode", "chm_storage_mode"):
        if performance[artifact] != "disk_backed_vrt":
            raise SystemExit(f"{artifact} was not produced through a disk-backed VRT")


def verify_binary_artifacts(output_dir: Path, dataset_id: str) -> None:
    png_path = output_dir / f"{dataset_id}_tiles.png"
    dtm_path = output_dir / f"{dataset_id}_dtm.tif"
    chm_path = output_dir / f"{dataset_id}_chm.tif"
    if png_path.read_bytes()[:8] != b"\x89PNG\r\n\x1a\n":
        raise SystemExit(f"missing or invalid PNG: {png_path}")
    for path in (dtm_path, chm_path):
        if path.read_bytes()[:4] not in (b"II*\x00", b"MM\x00*"):
            raise SystemExit(f"missing or invalid GeoTIFF: {path}")
    if (output_dir / "chm").exists():
        raise SystemExit("per-tile CHM directory must not be created")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("expected_tiles", type=int)
    parser.add_argument("expected_dimensions")
    parser.add_argument("expected_footprint_source")
    parser.add_argument("dataset_id")
    args = parser.parse_args()

    dimensions = (
        [] if args.expected_dimensions == "NA"
        else args.expected_dimensions.split("|")
    )
    expected_results = {
        verify_dimension(
            args.output_dir,
            args.dataset_id,
            dimension,
            args.expected_tiles,
        )
        for dimension in (dimensions or [None])
    }
    observed_results = set(args.output_dir.glob(f"{args.dataset_id}*_results.csv"))
    if observed_results != expected_results:
        raise SystemExit(
            "dimension-specific tile CSV set differs: "
            f"{observed_results} != {expected_results}"
        )
    verify_performance(
        args.output_dir,
        args.dataset_id,
        dimensions,
        args.expected_footprint_source,
    )
    verify_binary_artifacts(args.output_dir, args.dataset_id)


if __name__ == "__main__":
    main()
