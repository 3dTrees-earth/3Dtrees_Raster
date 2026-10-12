"""Regression coverage for pre-quantization source header bounds."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import importlib.util
import laspy
import numpy as np

VALIDATOR = Path(__file__).resolve().parents[1] / 'validation/validate_indexed_copc_streaming.py'

class BoundsValidation(unittest.TestCase):
    def test_crs_serializations_preserve_meaning_and_axis_order(self):
        from pyproj import CRS
        spec = importlib.util.spec_from_file_location('indexed_validator', VALIDATOR)
        validator = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(validator)
        crs = CRS.from_epsg(32632)
        modern = {'crs_wkt': crs.to_wkt()}
        legacy = {'crs_wkt': crs.to_wkt(version='WKT1_GDAL')}
        self.assertNotEqual(modern, legacy)
        self.assertTrue(validator.crs_compatible(modern, legacy))
        self.assertFalse(validator.crs_compatible(modern, {'crs_wkt': CRS.from_epsg(32633).to_wkt()}))
        self.assertFalse(validator.crs_compatible(modern, {'crs_wkt': None}))
        self.assertFalse(validator.crs_compatible(
            {'crs_wkt': CRS.from_epsg(4326).to_wkt()},
            {'crs_wkt': CRS.from_user_input('OGC:CRS84').to_wkt()}))

    def test_bounds_correction_preserves_strict_point_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            h = laspy.LasHeader(point_format=3, version='1.2')
            h.scales = np.full(3, .0001)
            from pyproj import CRS
            h.add_crs(CRS.from_epsg(32633))
            source = laspy.LasData(h)
            source.X = [100, 200, 300]
            source.Y = [300, 100, 200]
            source.Z = [100, 300, 200]
            source.classification = [1, 2, 3]
            original = root/'source.laz'
            source.write(original)
            candidate = laspy.read(original)
            candidate.add_extra_dim(laspy.ExtraBytesParams(name='OriginalPointIndex', type=np.uint64))
            candidate.OriginalPointIndex = [0, 1, 2]
            path = root/'candidate.laz'
            candidate.write(path)
            def run(expected, candidate_path=path):
                result = subprocess.run([sys.executable, str(VALIDATOR), str(original), str(candidate_path), str(root/'report.json'), '--dimension', 'classification', '--chunk-points', '2'], capture_output=True)
                self.assertEqual(result.returncode == 0, expected, result.stderr.decode())
                return json.loads((root/'report.json').read_text())
            self.assertTrue(run(True)['header_equal'])
            with original.open('r+b') as f:
                header = laspy.LasHeader.read_from(f)
                header.mins -= .00004
                header.maxs += .00004
                f.seek(0); header.write_to(f)
            report = run(True)
            self.assertFalse(report['header_equal'])
            self.assertTrue(report['header_compatible'])
            for kind in ('coordinate', 'classification', 'duplicate_index', 'candidate_header', 'scale', 'crs'):
                with self.subTest(kind=kind):
                    altered = laspy.read(path)
                    if kind == 'coordinate': altered.X[1] += 1
                    elif kind == 'classification': altered.classification[1] = 5
                    elif kind == 'duplicate_index': altered.OriginalPointIndex[1] = 0
                    elif kind == 'scale': altered.header.scales[0] *= 2
                    elif kind == 'crs': altered.header.add_crs(CRS.from_epsg(32632))
                    target = root/'altered.laz'; altered.write(target)
                    if kind == 'candidate_header':
                        with target.open('r+b') as f:
                            header = laspy.LasHeader.read_from(f)
                            header.maxs += .00004
                            f.seek(0); header.write_to(f)
                    run(False, target)

if __name__ == '__main__': unittest.main()
