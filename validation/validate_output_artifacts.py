#!/usr/bin/env python3
"""Validate the complete publishable artifact set for one dataset.

Only the four known GDAL PAM sidecars for the DTM and CHM are ignored. They
contain cached raster metadata, are not scientific products, and are excluded
from ``publish_artifacts`` in the resulting report. Every other unexpected file
remains a validation error.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SUPPORTED_DIMENSIONS = (
    "PredInstance",
    "PredInstance_SAT",
    "PredInstance_FM",
)
VALIDATION_REPORT = "forest_structure.validation.json"


def as_list(value: Any) -> list[Any]:
    """Normalize scalar or array provenance fields to a list."""

    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def expected_artifact_names(dataset_id: int | str, processed: list[str], voxel_counts: bool = False,
                            denoising: bool = False) -> set[str]:
    """Return the exact scientific and provenance artifact contract."""

    names = {
        f"{dataset_id}_dtm.tif",
        f"{dataset_id}_chm.tif",
        f"{dataset_id}_tiles.png",
        f"{dataset_id}_performance.csv",
        f"{dataset_id}_aoi_conversion.json",
        f"{dataset_id}_julia_memory_safe_run.json",
    }
    for dimension in processed or [""]:
        prefix = f"{dataset_id}_{dimension}" if dimension else str(dataset_id)
        names.update(
            {
                f"{prefix}_results.csv",
                f"{prefix}_segment_diagnostics.csv",
                f"{prefix}_tiles.geojson",
            }
        )
    if voxel_counts:
        names.update(f"{dataset_id}_voxel_counts{suffix}" for suffix in (".csv", ".tif", ".json", "_validity.tif"))
    if denoising:
        names.update((f"{dataset_id}_noise_points.csv", f"{dataset_id}_denoising.json"))
    return names


def allowed_gdal_sidecar_names(dataset_id: int | str) -> set[str]:
    """Return only the known non-scientific DTM/CHM PAM filenames."""

    return {
        f"{dataset_id}_{raster}.tif{suffix}.aux.xml"
        for raster in ("dtm", "chm")
        for suffix in ("", ".partial.tif")
    }


def validate_output_directory(dataset_id: int | str, output: Path) -> dict[str, Any]:
    """Return a machine-readable artifact validation report."""

    errors: list[str] = []
    provenance_path = output / f"{dataset_id}_julia_memory_safe_run.json"
    provenance: dict[str, Any] = {}
    try:
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        errors.append(f"cannot read provenance: {error}")

    processed = [
        str(value)
        for value in as_list(provenance.get("processed_instance_dimensions"))
    ]
    missing = [
        str(value)
        for value in as_list(provenance.get("missing_instance_dimensions"))
    ]
    requested = [
        str(value)
        for value in as_list(provenance.get("requested_instance_dimensions"))
    ] or list(SUPPORTED_DIMENSIONS)
    if (
        len(processed + missing) != len(set(processed + missing))
        or sorted(processed + missing) != sorted(requested)
    ):
        errors.append("processed and missing dimensions do not partition the request")
    unsupported = sorted(set(requested) - set(SUPPORTED_DIMENSIONS))
    if unsupported:
        errors.append(f"unsupported requested dimensions: {unsupported}")
    if provenance and provenance.get("status") != "valid":
        errors.append("provenance status is not valid")

    expected = expected_artifact_names(dataset_id, processed,
        provenance.get("voxel_counts", {}).get("enabled", False),
        provenance.get("denoising", {}).get("enabled", False))
    present = {path.name for path in output.iterdir() if path.is_file()}
    allowed_sidecars = allowed_gdal_sidecar_names(dataset_id)
    ignored_sidecars = sorted(present & allowed_sidecars)
    optional_metadata = {VALIDATION_REPORT}
    missing_artifacts = sorted(expected - present)
    unexpected_artifacts = sorted(
        present - expected - allowed_sidecars - optional_metadata
    )
    empty_artifacts = sorted(
        name
        for name in expected & present
        if (output / name).stat().st_size == 0
    )
    if missing_artifacts:
        errors.append(f"missing artifacts: {missing_artifacts}")
    if unexpected_artifacts:
        errors.append(f"unexpected artifacts: {unexpected_artifacts}")
    if empty_artifacts:
        errors.append(f"empty artifacts: {empty_artifacts}")

    return {
        "validated_at_utc": datetime.now(timezone.utc).isoformat(),
        "dataset_id": dataset_id,
        "status": "valid" if not errors else "invalid",
        "errors": errors,
        "processed_instance_dimensions": processed,
        "missing_instance_dimensions": missing,
        "publish_artifacts": sorted(expected),
        "ignored_gdal_sidecars": ignored_sidecars,
    }


def write_report(report: dict[str, Any], path: Path) -> None:
    """Write a report atomically so interrupted validation cannot look valid."""

    temporary = path.with_name(f".{path.name}.partial")
    temporary.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dataset_id", help="Output prefix: foreststructure, or a legacy dataset ID")
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--report", type=Path)
    arguments = parser.parse_args()

    report = validate_output_directory(
        arguments.dataset_id,
        arguments.output_directory.resolve(),
    )
    if arguments.report is not None:
        write_report(report, arguments.report.resolve())
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["status"] == "valid" else 1


if __name__ == "__main__":
    raise SystemExit(main())
