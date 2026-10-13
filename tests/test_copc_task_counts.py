"""Metadata-only fixtures exercise COPC scheduling estimates without lazrs."""

import importlib.util
import json
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


SCRIPT = Path(__file__).resolve().parents[1] / "src/copc_task_counts.py"
spec = importlib.util.spec_from_file_location("copc_task_counts", SCRIPT)
counts = importlib.util.module_from_spec(spec)
spec.loader.exec_module(counts)


def entry(key, offset, size, count):
    return counts.ENTRY.pack(*key, offset, size, count)


def fixture(path, *, pages=None, point_count=18, legacy_point_count=0, copc=True):
    """Write header/hierarchy and deliberately undecodable fake point chunks."""
    raw = bytearray(1500)
    raw[:4] = b"LASF"
    raw[24:26] = bytes((1, 4))
    struct.pack_into("<HII", raw, 94, 375, 589, int(copc))
    raw[104] = 0x86 if copc else 0x83
    struct.pack_into("<H", raw, 105, 30 if copc else 34)
    struct.pack_into("<I", raw, 107, legacy_point_count)
    struct.pack_into("<6d", raw, 179, 8, 0, 8, 0, 8, 0)
    struct.pack_into("<Q", raw, 247, point_count)
    if copc:
        raw[377:381] = b"copc"
        struct.pack_into("<HH", raw, 393, 1, 160)
        struct.pack_into("<5dQQ2d", raw, 429, 4, 4, 4, 4, 1, 1000, 96, 0, 1)
        if pages is None:
            pages = {1000: entry((0, 0, 0, 0), 600, 10, 3)
                     + entry((1, 0, 0, 0), 1200, 64, -1)
                     + entry((1, 1, 1, 1), 620, 10, 7),
                     1200: entry((1, 0, 0, 0), 610, 10, 5)
                     + entry((2, 0, 0, 0), 630, 10, 3)}
        for offset, data in pages.items():
            raw[offset:offset + len(data)] = data
        struct.pack_into("<Q", raw, 477, len(pages[1000]))
    path.write_bytes(raw)
    return path


class CopcTaskCounts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = fixture(Path(self.tmp.name) / "input.copc.laz")

    def test_all_depths_paginated_and_window_order(self):
        result = counts.estimate(self.path, [[6, 6, 7, 7], [0, 0, 8, 8], [0, 0, 1, 1], [2.5, 2.5, 3, 3]])
        self.assertEqual(result["point_count_upper_bounds"], [10, 18, 11, 8])
        self.assertEqual(result["hierarchy_page_count"], 2)
        self.assertEqual(result["hierarchy_point_count"], 18)
        self.assertEqual(result["hierarchy_node_count"], 4)
        self.assertEqual(result["metadata_bytes_read"], 589 + 160)

    def test_inclusive_edges_and_outside_header(self):
        result = counts.estimate(self.path, [[4, 4, 4, 4], [8, 8, 8, 8], [8.1, 0, 9, 1], [-2, -2, -1, -1]])
        self.assertEqual(result["point_count_upper_bounds"], [15, 10, 0, 0])

    def test_window_tree_matches_simple_intersections(self):
        windows = [[x / 10, y / 10, x / 10 + .2, y / 10 + .2] for x in range(-10, 100, 3) for y in range(-10, 100, 5)]
        result = counts.estimate(self.path, windows)
        boxes = [((0, 0, 8, 8), 3), ((0, 0, 4, 4), 5), ((4, 4, 8, 8), 7), ((0, 0, 2, 2), 3)]
        expected = [sum(n for bounds, n in boxes if counts._intersects(window, bounds)) for window in windows]
        self.assertEqual(result["point_count_upper_bounds"], expected)

    def test_fallback_explicit_and_conservative(self):
        fixture(self.path, copc=False)
        with self.assertRaisesRegex(counts.MetadataError, "explicit opt-in"):
            counts.estimate(self.path, [[0, 0, 1, 1]])
        result = counts.estimate(self.path, [[0, 0, 1, 1], [10, 10, 11, 11]], allow_header_fallback=True)
        self.assertEqual(result["mode"], "header_fallback")
        self.assertEqual(result["point_count_upper_bounds"], [18, 0])

    def test_las14_legacy_count_fallback(self):
        fixture(self.path, copc=False, point_count=0, legacy_point_count=18)
        result = counts.estimate(self.path, [[0, 0, 1, 1], [10, 10, 11, 11]], allow_header_fallback=True)
        self.assertEqual(result["point_count"], 18)
        self.assertEqual(result["point_count_upper_bounds"], [18, 0])

    def test_las14_conflicting_counts_reserve_larger_count(self):
        for legacy, extended in ((18, 12), (12, 18)):
            with self.subTest(legacy=legacy, extended=extended):
                fixture(self.path, copc=False, point_count=extended, legacy_point_count=legacy)
                result = counts.estimate(self.path, [[0, 0, 1, 1]], allow_header_fallback=True)
                self.assertEqual(result["point_count_upper_bounds"], [18])

    def test_copc_never_uses_legacy_count(self):
        for legacy in (12, 18, 24):
            with self.subTest(legacy=legacy):
                fixture(self.path, point_count=18, legacy_point_count=legacy)
                result = counts.estimate(self.path, [[0, 0, 8, 8]], allow_header_fallback=True)
                self.assertEqual(result["point_count"], 18)
                self.assertEqual(result["point_count_upper_bounds"], [18])
        for extended in (0, 12, 24):
            with self.subTest(extended=extended):
                fixture(self.path, point_count=extended, legacy_point_count=18)
                with self.assertRaisesRegex(counts.MetadataError, "exceeds LAS header|differs from LAS header"):
                    counts.estimate(self.path, [[0, 0, 8, 8]], allow_header_fallback=True)

    def test_invalid_windows(self):
        for windows in (None, [[1, 0, 0, 1]], [[0, 0, float("nan"), 1]], [[0, 0, 1]], [[False, 0, 1, 1]], [["0", 0, 1, 1]], [[0, 0, 10 ** 400, 1]]):
            with self.subTest(windows=windows), self.assertRaises(counts.MetadataError):
                counts.estimate(self.path, windows)

    def test_count_mismatch_even_for_no_windows(self):
        fixture(self.path, point_count=19)
        with self.assertRaisesRegex(counts.MetadataError, "differs"):
            counts.estimate(self.path, [])

    def test_malformed_hierarchy(self):
        cases = {
            "cycle": ({1000: entry((0, 0, 0, 0), 1000, 32, -1)}, "cycle"),
            "out_of_file": ({1000: entry((0, 0, 0, 0), 1490, 32, -1)}, "bounds"),
            "misaligned": ({1000: b"X" * 33}, "alignment"),
            "negative_count": ({1000: entry((0, 0, 0, 0), 600, 10, -2)}, "negative"),
            "duplicate_key": ({1000: entry((0, 0, 0, 0), 600, 10, 9) + entry((0, 0, 0, 0), 620, 10, 9)}, "duplicate.*key"),
            "overlap": ({1000: entry((0, 0, 0, 0), 600, 20, 9) + entry((1, 0, 0, 0), 610, 20, 9)}, "overlapping"),
            "invalid_key": ({1000: entry((1, 2, 0, 0), 600, 10, 18)}, "voxel key"),
            "empty_range": ({1000: entry((0, 0, 0, 0), 600, 0, 18)}, "empty"),
            "empty_node_data": ({1000: entry((0, 0, 0, 0), 600, 10, 0)}, "empty hierarchy"),
        }
        for name, (pages, message) in cases.items():
            with self.subTest(name=name):
                fixture(self.path, pages=pages)
                with self.assertRaisesRegex(counts.MetadataError, message):
                    counts.estimate(self.path, [[0, 0, 8, 8]], allow_header_fallback=True)

    def test_metadata_budget_and_truncation(self):
        with self.assertRaisesRegex(counts.MetadataError, "limit"):
            counts.estimate(self.path, [], max_metadata_bytes=600)
        self.path.write_bytes(self.path.read_bytes()[:1240])
        with self.assertRaisesRegex(counts.MetadataError, "bounds"):
            counts.estimate(self.path, [])

    def test_invalid_copc_info_never_falls_back(self):
        for offset, data, message in (
            (453, struct.pack("<d", 0), "size"),
            (453, struct.pack("<d", 1e308), "range"),
            (429, struct.pack("<d", 50), "enclose"),
            (501, b"\1", "reserved"),
            (104, b"\x83", "point format"),
        ):
            with self.subTest(offset=offset):
                fixture(self.path)
                raw = bytearray(self.path.read_bytes())
                raw[offset:offset + len(data)] = data
                self.path.write_bytes(raw)
                with self.assertRaisesRegex(counts.MetadataError, message):
                    counts.estimate(self.path, [], allow_header_fallback=True)

    def test_header_bounds_override_broader_octree(self):
        raw = bytearray(self.path.read_bytes())
        struct.pack_into("<6d", raw, 179, 7, 1, 7, 1, 7, 1)
        self.path.write_bytes(raw)
        result = counts.estimate(self.path, [[0, 0, .5, .5], [1, 1, 1, 1]])
        self.assertEqual(result["point_count_upper_bounds"], [0, 11])

    def test_large_page_reads_stay_bounded(self):
        page = entry((0, 0, 0, 0), 600, 10, 18)
        page += b"".join(entry((12, i, 0, 0), 0, 0, 0) for i in range(3000))
        fixture(self.path, pages={1000: page})
        requested = []
        read = counts._MetadataReader.read

        def observed(reader, offset, size):
            requested.append(size)
            return read(reader, offset, size)

        with mock.patch.object(counts._MetadataReader, "read", observed):
            result = counts.estimate(self.path, [[0, 0, 8, 8]])
        self.assertEqual(result["point_count_upper_bounds"], [18])
        self.assertEqual(result["metadata_bytes_read"], 589 + len(page))
        self.assertLessEqual(max(requested), 64 * 1024)

    def test_cli_atomic_result(self):
        windows = Path(self.tmp.name) / "windows.json"
        output = Path(self.tmp.name) / "counts.json"
        windows.write_text(json.dumps([[0, 0, 8, 8]]))
        command = [sys.executable, str(SCRIPT), "--input", str(self.path), "--windows", str(windows), "--output", str(output)]
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(output.read_text())["point_count_upper_bounds"], [18])
        before = output.read_bytes()
        fixture(self.path, point_count=19)
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(output.read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
