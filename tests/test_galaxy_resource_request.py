"""Resource boundary, cap and Galaxy expression checks without point-cloud I/O."""

import ast
import json
from pathlib import Path
import random
import subprocess
import sys
import unittest
import xml.etree.ElementTree as ET
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from galaxy_resource_request import MIB, request_resources

try:
    import yaml
except ImportError:
    yaml = None


class ResourceRequestTests(unittest.TestCase):
    def test_decimal_boundaries_and_cpu_transition(self):
        cases = [(1, 8, 28610), (1_000_000_000, 8, 28610),
                 (1_000_000_001, 10, 28610), (5_000_000_000, 10, 28610),
                 (5_000_000_001, 10, 57220), (10_000_000_000, 10, 57220),
                 (10_000_000_001, 10, 95367), (10**30, 10, 95367)]
        for size, cpus, mib in cases:
            with self.subTest(size=size):
                request = request_resources(size)
                self.assertEqual((request.cpus, request.memory_mib), (cpus, mib))
                self.assertEqual(request.memory_gib * 1024, mib)

    def test_recorded_merged_copc_sizes(self):
        # These are the merged/prepared COPCs, not original upload sizes.
        cases = {647: (879_019_410, 8, 28610), 2080: (795_154_736, 8, 28610),
                 2551: (1_093_956_560, 10, 28610), 1998: (3_910_169_162, 10, 28610),
                 1805: (13_332_084_830, 10, 95367)}
        for dataset, (size, cpus, mib) in cases.items():
            with self.subTest(dataset=dataset):
                request = request_resources(size)
                self.assertEqual((request.cpus, request.memory_mib), (cpus, mib))

    def test_monotone_and_never_above_decimal_ceiling(self):
        rng = random.Random(91)
        sizes = sorted({1, 10**20, *(rng.randrange(1, 20_000_000_001) for _ in range(1000))})
        requests = [request_resources(size) for size in sizes]
        self.assertEqual([r.cpus for r in requests], sorted(r.cpus for r in requests))
        self.assertEqual([r.memory_mib for r in requests], sorted(r.memory_mib for r in requests))
        for request in requests:
            self.assertLessEqual(request.cpus, 10)
            self.assertLessEqual(request.memory_mib * MIB, 100_000_000_000)

    def test_missing_invalid_metadata_fails_without_file_reads(self):
        with patch("builtins.open", side_effect=AssertionError("unexpected file read")), \
             patch("os.stat", side_effect=AssertionError("unexpected file stat")):
            self.assertEqual(request_resources(13_332_084_830).memory_mib, 95367)
            for value in [None, 0, -1, True, False, 1.0, float("inf"), float("nan"), "12"]:
                with self.subTest(value=value), self.assertRaisesRegex(ValueError, "positive integer"):
                    request_resources(value)

    def test_cli_reports_request_and_rejects_invalid_metadata(self):
        command = [sys.executable, str(ROOT / "src/galaxy_resource_request.py"), "--input-bytes"]
        valid = subprocess.run(command + ["13332084830"], text=True, capture_output=True, check=True)
        self.assertEqual(json.loads(valid.stdout)["memory_mib"], 95367)
        for value in ["0", "-1", "nan", "1.5"]:
            invalid = subprocess.run(command + [value], text=True, capture_output=True)
            self.assertNotEqual(invalid.returncode, 0)
            self.assertFalse(invalid.stdout)

    @unittest.skipIf(yaml is None, "PyYAML needed to parse the TPV integration example")
    def test_tpv_example_uses_only_merged_input_and_preserves_units(self):
        config = yaml.safe_load((ROOT / "galaxy/tpv.example.yml").read_text())
        pattern = ".*/3dtrees_foreststructure/3dtrees_foreststructure/.*"
        rule = config["tools"][pattern]
        # The rule must select the wrapper's Tool Shed GUID and sized input.
        guid = "toolshed.g2.bx.psu.edu/repos/bgruening/3dtrees_foreststructure/3dtrees_foreststructure/1.0.0+galaxy0"
        self.assertRegex(guid, f"^{pattern}$")
        wrapper = ET.parse(ROOT / "galaxy/3dtrees_foreststructure/foreststructure.xml").getroot()
        self.assertEqual(wrapper.get("id"), "3dtrees_foreststructure")
        self.assertIsNotNone(wrapper.find("inputs/param[@name='point_cloud']"))

        class Helpers:
            def __init__(self, size):
                self.size = size

            def get_input_size(self, job, parameter, *, estimate_uncompressed_size):
                self.assertions = (job == "test-job", parameter == "point_cloud",
                                   estimate_uncompressed_size is False)
                return self.size / 1024**3

        def evaluate(block, helpers):
            # TPV expression blocks evaluate their final expression as the value.
            tree = ast.parse(block)
            scope = {"helpers": helpers, "job": "test-job"}
            exec(compile(ast.Module(body=tree.body[:-1], type_ignores=[]), "tpv-example", "exec"), scope)
            return eval(compile(ast.Expression(tree.body[-1].value), "tpv-example", "eval"), scope)

        for size in [1, 1_000_000_000, 1_000_000_001, 5_000_000_001, 13_332_084_830]:
            helper = Helpers(size)
            expected = request_resources(size)
            self.assertEqual(evaluate(rule["cores"], helper), expected.cpus)
            self.assertEqual(evaluate(rule["mem"], helper), expected.memory_gib)
            self.assertTrue(all(helper.assertions))
        self.assertEqual(rule["gpus"], 0)
        self.assertEqual(rule["max_cores"], 10)
        self.assertLessEqual(evaluate(rule["max_mem"], Helpers(1)) * 1024**3, 100_000_000_000)
        with self.assertRaisesRegex(ValueError, "positive integer"):
            evaluate(rule["mem"], Helpers(0))


if __name__ == "__main__":
    unittest.main()
