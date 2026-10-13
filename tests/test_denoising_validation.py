"""Quality-product corruption must fail before artifact publication."""

import importlib.util
import json
import sqlite3
import tempfile
import unittest
from pathlib import Path

SPEC = importlib.util.spec_from_file_location("denoising_validation",
    Path(__file__).resolve().parents[1] / "src/denoising_validation.py")
VALIDATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VALIDATOR)


class DenoisingValidationTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        root = Path(self.directory.name)
        self.metadata = root / "denoising.json"
        self.csv = root / "noise_points.csv"
        self.count = 2_000_000_000
        self.summary = dict(enabled=True, method="fixed_window_local_sor",
            point_order_dimension="OriginalPointIndex", coordinate_changes=False,
            raw_voxel_counts_filtered=False, point_count=self.count,
            retained_points=self.count - 2, rejected_points=2,
            insufficient_context_points=0, partition_count=2,
            partitions=[dict(ix=0, iy=0, core_points=self.count - 1,
                rejected_points=1, accepted_points=self.count - 2,
                insufficient_context_points=0),
                dict(ix=1, iy=0, core_points=1, rejected_points=1,
                     accepted_points=0, insufficient_context_points=0)])
        self.csv.write_text("OriginalPointIndex\n1999999999\n0\n")

    def validate(self):
        self.metadata.write_text(json.dumps(self.summary))
        return VALIDATOR.validate(self.metadata, self.csv, self.count)

    def test_large_source_ids_without_a_whole_source_memory_table(self):
        self.assertEqual(self.validate()["rejected_points"], 2)

    def test_repeated_id_across_partitions_fails(self):
        self.csv.write_text("OriginalPointIndex\n0\n0\n")
        with self.assertRaises(sqlite3.IntegrityError):
            self.validate()

    def test_bad_or_missing_csv_identity_fails(self):
        for content in ("OriginalPointIndex\n2000000000\n0\n",
                        "OriginalPointIndex\n1\n",
                        "wrong_column\n1\n2\n"):
            with self.subTest(content=content):
                self.csv.write_text(content)
                with self.assertRaises(ValueError):
                    self.validate()

    def test_inconsistent_partition_accounting_fails(self):
        self.summary["partitions"][0]["accepted_points"] -= 1
        with self.assertRaises(ValueError):
            self.validate()

    def test_filtered_raw_counts_cannot_be_labelled_as_this_product(self):
        self.summary["raw_voxel_counts_filtered"] = True
        with self.assertRaises(ValueError):
            self.validate()


if __name__ == "__main__":
    unittest.main()
