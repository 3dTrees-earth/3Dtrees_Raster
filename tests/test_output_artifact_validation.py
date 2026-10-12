#!/usr/bin/env python3
"""Focused tests for the reusable output artifact validator."""

from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
VALIDATOR_PATH = REPOSITORY_ROOT / "validation" / "validate_output_artifacts.py"


def load_validator():
    spec = importlib.util.spec_from_file_location("output_validator", VALIDATOR_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {VALIDATOR_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class OutputArtifactValidationTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.output = Path(self.temporary_directory.name)
        self.dataset_id = 2001
        provenance = {
            "status": "valid",
            "requested_instance_dimensions": [
                "PredInstance",
                "PredInstance_SAT",
                "PredInstance_FM",
            ],
            "processed_instance_dimensions": [
                "PredInstance_SAT",
                "PredInstance_FM",
            ],
            "missing_instance_dimensions": "PredInstance",
        }
        (self.output / "2001_julia_memory_safe_run.json").write_text(
            json.dumps(provenance),
            encoding="utf-8",
        )
        shared = [
            "2001_dtm.tif",
            "2001_chm.tif",
            "2001_tiles.png",
            "2001_performance.csv",
            "2001_aoi_conversion.json",
        ]
        dimension_files = [
            f"2001_{dimension}_{suffix}"
            for dimension in ("PredInstance_SAT", "PredInstance_FM")
            for suffix in (
                "results.csv",
                "segment_diagnostics.csv",
                "tiles.geojson",
            )
        ]
        for name in shared + dimension_files:
            (self.output / name).write_bytes(b"non-empty")

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def test_known_gdal_sidecars_are_reported_and_ignored(self) -> None:
        validator = load_validator()
        sidecars = {
            "2001_dtm.tif.aux.xml",
            "2001_dtm.tif.partial.tif.aux.xml",
            "2001_chm.tif.aux.xml",
            "2001_chm.tif.partial.tif.aux.xml",
        }
        for name in sidecars:
            (self.output / name).write_text("<PAMDataset/>", encoding="utf-8")
        (self.output / "forest_structure.validation.json").write_text(
            "{}",
            encoding="utf-8",
        )

        report = validator.validate_output_directory(self.dataset_id, self.output)

        self.assertEqual(report["status"], "valid")
        self.assertEqual(set(report["ignored_gdal_sidecars"]), sidecars)

    def test_unknown_extra_file_is_rejected(self) -> None:
        validator = load_validator()
        (self.output / "unexpected.txt").write_text("unexpected", encoding="utf-8")

        report = validator.validate_output_directory(self.dataset_id, self.output)

        self.assertEqual(report["status"], "invalid")
        self.assertTrue(
            any("unexpected artifacts" in error for error in report["errors"])
        )

    def test_enabled_voxel_product_requires_all_four_files(self) -> None:
        validator = load_validator()
        path = self.output / "2001_julia_memory_safe_run.json"
        provenance = json.loads(path.read_text())
        provenance["voxel_counts"] = {"enabled": True}
        path.write_text(json.dumps(provenance))
        self.assertEqual(validator.validate_output_directory(2001, self.output)["status"], "invalid")
        for suffix in (".csv", ".tif", ".json", "_validity.tif"):
            (self.output / f"2001_voxel_counts{suffix}").write_bytes(b"non-empty")
        self.assertEqual(validator.validate_output_directory(2001, self.output)["status"], "valid")

    def test_unsegmented_source_uses_existing_unsuffixed_outputs(self) -> None:
        names = load_validator().expected_artifact_names(2001, [], True)
        self.assertIn("2001_results.csv", names)
        self.assertIn("2001_tiles.geojson", names)
        self.assertIn("2001_voxel_counts.tif", names)
        self.assertNotIn("2001__results.csv", names)


if __name__ == "__main__":
    unittest.main()
