#!/usr/bin/env python3
"""Compare candidate tables/vectors with a Julia-faithful oracle exactly.

Numeric tolerances are intentionally absent: the ordered-COPC acceptance goal is
reproduction, not statistical similarity.  Only the input filename metadata is
ignored because the candidate reads ``input.copc.laz`` while the oracle names
the ordered LAS/LAZ companion.
"""

from __future__ import annotations

import csv
import json
import pathlib
import sys
from typing import Any


# Storage filename is expected to differ; every scientific field must remain.
IGNORED_METADATA_FIELDS = {"file"}


def normalized_csv(path: pathlib.Path) -> tuple[list[str], list[list[str]]]:
    """Read a CSV and remove only explicitly non-scientific metadata columns."""

    # Preserve values as strings so equality is byte-level after CSV parsing and
    # no floating-point coercion can hide representation differences.
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    # Empty output is always invalid, including datasets with zero result rows,
    # because even those files must contain the expected header.
    if not rows:
        raise AssertionError(f"{path} is empty")
    # Resolve retained column ordinals from the file's own header.
    keep = [
        index
        for index, name in enumerate(rows[0])
        if name not in IGNORED_METADATA_FIELDS
    ]
    # Return the filtered header and every data row in original order.
    return (
        [rows[0][index] for index in keep],
        [[row[index] for index in keep] for row in rows[1:]],
    )


def compare_csv(expected: pathlib.Path, actual: pathlib.Path) -> dict[str, Any]:
    """Compare headers, row counts, and all retained cell strings exactly."""

    # Normalize expected and actual independently so reordered headers are caught.
    expected_header, expected_rows = normalized_csv(expected)
    actual_header, actual_rows = normalized_csv(actual)
    # Retain only the first 20 details to keep failure reports bounded.
    first_differences: list[dict[str, Any]] = []
    differing_values = 0
    # Compare corresponding rows and columns; separate length/header checks below
    # catch extra rows or columns that ``zip`` intentionally does not traverse.
    for row_index, (expected_row, actual_row) in enumerate(
        zip(expected_rows, actual_rows), start=1
    ):
        for column, expected_value, actual_value in zip(
            expected_header, expected_row, actual_row
        ):
            # Equal string values need no further reporting.
            if expected_value == actual_value:
                continue
            # Count every mismatch even after the detail sample is full.
            differing_values += 1
            if len(first_differences) < 20:
                first_differences.append(
                    {
                        "row": row_index,
                        "column": column,
                        "expected": expected_value,
                        "actual": actual_value,
                    }
                )
    # Require schema, cardinality, and every compared cell to agree.
    valid = (
        expected_header == actual_header
        and len(expected_rows) == len(actual_rows)
        and differing_values == 0
    )
    # Return structured evidence for the combined JSON acceptance report.
    return {
        "status": "valid" if valid else "invalid",
        "expected_rows": len(expected_rows),
        "actual_rows": len(actual_rows),
        "headers_equal": expected_header == actual_header,
        "differing_values": differing_values,
        "first_differences": first_differences,
    }


def remove_ignored_metadata(value: Any) -> Any:
    """Recursively remove allowed metadata without disturbing array order."""

    # GeoJSON metadata may occur at collection, feature, or property depth.
    if isinstance(value, dict):
        return {
            key: remove_ignored_metadata(child)
            for key, child in value.items()
            if key not in IGNORED_METADATA_FIELDS
        }
    # Preserve list ordering because feature/coordinate order is part of parity.
    if isinstance(value, list):
        return [remove_ignored_metadata(child) for child in value]
    # Scalars, including all numeric coordinates/properties, remain unchanged.
    return value


def compare_geojson(
    expected: pathlib.Path, actual: pathlib.Path
) -> dict[str, Any]:
    """Compare parsed GeoJSON structures after removing allowed metadata."""

    # Parse JSON to ignore insignificant whitespace while retaining exact values.
    expected_value = remove_ignored_metadata(json.loads(expected.read_text()))
    actual_value = remove_ignored_metadata(json.loads(actual.read_text()))
    # Python structural equality checks keys, values, feature order, and arrays.
    valid = expected_value == actual_value
    return {
        "status": "valid" if valid else "invalid",
        "expected_features": len(expected_value.get("features", [])),
        "actual_features": len(actual_value.get("features", [])),
        "equal": valid,
    }


COLLECTION_FOLDERS = {
    "results": "tile_metrics",
    "segment_diagnostics": "segment_diagnostics",
    "tiles": "tiles",
}


def artifact_path(directory: pathlib.Path, stem: str, dimension: str, kind: str,
                  extension: str) -> pathlib.Path:
    """Return the collection-folder file, or the flat legacy name."""

    folder_path = directory / COLLECTION_FOLDERS[kind] / f"{dimension}.{extension}"
    return folder_path if folder_path.exists() else directory / f"{stem}_{kind}.{extension}"


def main() -> None:
    """Compare all scientific artifacts and write one machine-readable verdict."""

    # The shell harness passes exactly one dataset, dimension, oracle, candidate,
    # and report path per process.
    if len(sys.argv) != 6:
        raise SystemExit(
            "usage: compare_valid_updated_science.py DATASET_ID DIMENSION "
            "VALID_UPDATED_DIR CANDIDATE_DIR REPORT_JSON"
        )
    # Keep the dataset text for filenames and convert it only in report metadata.
    dataset_id, dimension = sys.argv[1:3]
    # Resolve path objects without requiring the output report to pre-exist.
    expected_dir = pathlib.Path(sys.argv[3])
    actual_dir = pathlib.Path(sys.argv[4])
    report_path = pathlib.Path(sys.argv[5])
    # ForestStructure uses one shared stem for each instance-dimension artifact set.
    stem = f"{dataset_id}_{dimension}"
    # Compare the two CSV products and the vector product independently so the
    # report identifies the exact failing artifact class.
    comparisons = {
        "results_csv": compare_csv(
            artifact_path(expected_dir, stem, dimension, "results", "csv"),
            artifact_path(actual_dir, stem, dimension, "results", "csv"),
        ),
        "segment_diagnostics_csv": compare_csv(
            artifact_path(expected_dir, stem, dimension, "segment_diagnostics", "csv"),
            artifact_path(actual_dir, stem, dimension, "segment_diagnostics", "csv"),
        ),
        "tiles_geojson": compare_geojson(
            artifact_path(expected_dir, stem, dimension, "tiles", "geojson"),
            artifact_path(actual_dir, stem, dimension, "tiles", "geojson"),
        ),
    }
    # The dataset passes only when every scientific artifact passes.
    valid = all(item["status"] == "valid" for item in comparisons.values())
    # Persist summary and detailed component evidence in one durable JSON object.
    report = {
        "status": "valid" if valid else "invalid",
        "dataset_id": int(dataset_id),
        "dimension": dimension,
        "comparison": "zero-tolerance scientific equality",
        "ignored_metadata_fields": sorted(IGNORED_METADATA_FIELDS),
        **comparisons,
    }
    # Pretty JSON plus a newline is readable in both PR artifacts and shell tools.
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    # Nonzero exit prevents the shell harness from accepting later raster success
    # after any tabular/vector mismatch.
    if not valid:
        raise SystemExit(
            f"valid_updated alignment failed for dataset {dataset_id}; "
            f"see {report_path}"
        )
    # Emit one concise success message; exact counts remain in the report.
    print(
        f"valid_updated alignment passed for dataset {dataset_id}: "
        "zero CSV/GeoJSON scientific differences"
    )


if __name__ == "__main__":
    # Keep helper functions importable without triggering CLI processing.
    main()
