"""Real LAS/LAZ, DTM, GIS and CSV round trips; no production data/services."""
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import struct
import unittest

import laspy
import numpy as np
from osgeo import gdal, ogr, osr

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("voxel_counts", ROOT / "src/voxel_counts.py")
voxel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(voxel)


class VoxelCountsTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.srs = osr.SpatialReference()
        self.srs.ImportFromEPSG(32632)
        header = laspy.LasHeader(point_format=0, version="1.2")
        header.scales = np.array([.001, .001, .001])
        cloud = laspy.LasData(header)
        self.xyz = np.array([
            [1000.1, 2000.1, 0], [1000.1, 2000.1, 0],
            [1000.2, 2000.2, .49], [1000.2, 2000.2, .5],
            [1000.5, 2000.5, 5], [1002.1, 2001.1, 3],
            [1001.2, 2002.2, 1], [1003.2, 2003.2, 80.25],
            [1003.1, 2003.1, -.01], [1000, 2000, .1]])
        cloud.x, cloud.y, cloud.z = self.xyz.T
        cloud.write(self.root / "source.laz")
        dtm = gdal.GetDriverByName("GTiff").Create(str(self.root / "dtm.tif"), 6, 6, 1, gdal.GDT_Float64)
        dtm.SetGeoTransform((999, 1, 0, 2005, 0, -1))
        dtm.SetProjection(self.srs.ExportToWkt())
        values = np.zeros((6, 6))
        values[2, 2] = -9999
        dtm.GetRasterBand(1).SetNoDataValue(-9999)
        dtm.GetRasterBand(1).WriteArray(values)
        dtm = None
        aoi = voxel.rectangle(1000, 2000, 1004, 2004).Difference(
            voxel.rectangle(1002, 2001, 1002.5, 2001.5))
        ds = ogr.GetDriverByName("GPKG").CreateDataSource(str(self.root / "aoi.gpkg"))
        layer = ds.CreateLayer("usable", self.srs, ogr.wkbUnknown)
        feature = ogr.Feature(layer.GetLayerDefn())
        feature.SetGeometry(aoi)
        layer.CreateFeature(feature)
        ds = None

    def tearDown(self):
        self.temporary.cleanup()

    def run_export(self, name="out", chunk=3, source="source.laz", memory=.25):
        out = self.root / name
        out.mkdir()
        paths = dict(csv=str(out / "counts.csv"), raster=str(out / "counts.tif"),
                     mask=str(out / "validity.tif"), metadata=str(out / "counts.json"))
        result = voxel.export(dict(point_cloud=str(self.root / source),
            dtm=str(self.root / "dtm.tif"), aoi=str(self.root / "aoi.gpkg"),
            temp_dir=str(out), paths=paths, chunk_points=chunk, memory_available_gib=memory))
        with open(paths["csv"]) as file:
            rows = list(csv.DictReader(file))
        return result, rows, paths

    def test_exact_counts_layers_and_roundtrip(self):
        result, rows, paths = self.run_export()
        self.assertEqual(result["point_accounting"], dict(read=10, nonfinite=0,
            excluded_domain=1, invalid_dtm=1, negative_height=1, included=7))
        observed = {(float(r["x_center"]), float(r["y_center"]), float(r["height_lower_m"])):
                    int(r["point_count"]) for r in rows}
        self.assertEqual(observed, {(1000.25, 2000.25, 0): 4,
            (1000.25, 2000.25, .5): 1, (1000.75, 2000.75, 5): 1,
            (1003.25, 2003.25, 80): 1})
        cube, mask = gdal.Open(paths["raster"]), gdal.Open(paths["mask"])
        self.assertEqual(cube.RasterCount, 161)
        self.assertEqual(cube.GetRasterBand(1).GetNoDataValue(), voxel.NODATA)
        valid = mask.ReadAsArray() == 1
        raster_total = 0
        for k in range(cube.RasterCount):
            band = cube.GetRasterBand(k+1)
            values = band.ReadAsArray()
            self.assertTrue(np.all(values[~valid] == voxel.NODATA))
            raster_total += int(values[valid].sum())
            self.assertEqual(float(band.GetMetadata()["height_lower_m"]), k*.5)
        self.assertEqual(raster_total, 7)
        for row in rows:
            x, y = float(row["x_center"]), float(row["y_center"])
            gt = cube.GetGeoTransform()
            col, r = int((x-gt[0])/.5), int((gt[3]-y)/.5)
            band = int(float(row["height_lower_m"])/.5)+1
            self.assertEqual(int(cube.GetRasterBand(band).ReadAsArray(col, r, 1, 1)[0, 0]), int(row["point_count"]))
        # Ten fine bands form [0,5); the point at exactly 5 belongs to the next bin.
        self.assertEqual(sum(int(cube.GetRasterBand(k).ReadAsArray()[valid].sum()) for k in range(1, 11)), 5)
        self.assertTrue(np.any(cube.GetRasterBand(1).ReadAsArray()[valid] == 0))
        self.assertFalse(any(Path(paths["csv"]).parent.glob("voxel_counts_*")))

    def test_batch_size_and_source_order_invariance(self):
        first, rows, paths = self.run_export("small", 1)
        cloud = laspy.read(self.root / "source.laz")
        cloud.points = cloud.points[::-1].copy()
        cloud.write(self.root / "reordered.las")
        second, rows2, paths2 = self.run_export("large", 1000, "reordered.las")
        self.assertEqual(rows, rows2)
        self.assertEqual(first["point_accounting"], second["point_accounting"])
        self.assertEqual(first["peak_batch_points"], 1)
        self.assertEqual(second["peak_batch_points"], 10)
        self.assertEqual(Path(paths["csv"]).read_bytes(), Path(paths2["csv"]).read_bytes())
        np.testing.assert_array_equal(gdal.Open(paths["raster"]).ReadAsArray(),
                                      gdal.Open(paths2["raster"]).ReadAsArray())

    def test_unknown_crs_fails(self):
        dtm = gdal.Open(str(self.root / "dtm.tif"), gdal.GA_Update)
        dtm.SetProjection("")
        dtm = None
        with self.assertRaisesRegex((ValueError, RuntimeError), "CRS|WKT"):
            self.run_export()

    def test_all_negative_keeps_valid_zero_domain(self):
        cloud = laspy.read(self.root / "source.laz")
        cloud.z = np.full(len(cloud.points), -1.)
        cloud.write(self.root / "negative.laz")
        result, rows, paths = self.run_export(source="negative.laz")
        self.assertEqual(rows, [])
        self.assertEqual(result["band_count"], 1)
        self.assertEqual(result["point_accounting"]["included"], 0)
        mask = gdal.Open(paths["mask"]).ReadAsArray()
        self.assertTrue(np.all(gdal.Open(paths["raster"]).ReadAsArray()[mask == 1] == 0))

    def test_extreme_height_fails_without_silent_cap(self):
        cloud = laspy.read(self.root / "source.laz")
        cloud.z = np.full(len(cloud.points), 40000.)
        cloud.write(self.root / "extreme.laz")
        with self.assertRaisesRegex(ValueError, "band capacity"):
            self.run_export(source="extreme.laz")

    def test_million_points_stream_in_bounded_batches(self):
        cloud = laspy.read(self.root / "source.laz")
        batch = laspy.ScaleAwarePointRecord(np.tile(cloud.points.array, 1000),
            cloud.header.point_format, cloud.header.scales, cloud.header.offsets)
        with laspy.open(self.root / "million.laz", mode="w", header=cloud.header) as writer:
            for _ in range(100):
                writer.write_points(batch)
        result, rows, _ = self.run_export(chunk=8192, source="million.laz")
        self.assertEqual(result["point_accounting"]["read"], 1000000)
        self.assertEqual(result["point_accounting"]["included"], 700000)
        self.assertEqual(result["peak_batch_points"], 8192)
        self.assertEqual(sum(int(row["point_count"]) for row in rows), 700000)
        print("Million-point streaming fixture:", json.dumps({key: result[key] for key in
            ("elapsed_seconds", "peak_batch_points", "temporary_disk_bytes")}))

    def test_counts_across_raster_blocks_and_exact_xy_edges(self):
        cloud = laspy.LasData(laspy.LasHeader(point_format=0, version="1.2"))
        cloud.x = np.array([1000.1, 1064., 1129.9])
        cloud.y = np.array([2000.1, 2064., 2129.9])
        cloud.z = np.array([0., .5, 1.])
        cloud.write(self.root / "wide.laz")
        dtm = gdal.GetDriverByName("GTiff").Create(str(self.root / "dtm.tif"), 132, 132, 1, gdal.GDT_Float64)
        dtm.SetGeoTransform((999, 1, 0, 2131, 0, -1))
        dtm.SetProjection(self.srs.ExportToWkt())
        dtm.GetRasterBand(1).Fill(0)
        dtm = None
        aoi = ogr.Open(str(self.root / "aoi.gpkg"), 1)
        layer = aoi.GetLayer(0)
        feature = layer.GetNextFeature()
        feature.SetGeometry(voxel.rectangle(1000, 2000, 1130, 2130))
        layer.SetFeature(feature)
        aoi = None
        result, rows, paths = self.run_export(chunk=1, source="wide.laz")
        self.assertEqual(result["point_accounting"]["included"], 3)
        observed = {(float(r["x_center"]), float(r["y_center"])) for r in rows}
        self.assertEqual(observed, {(1000.25, 2000.25), (1064.25, 2064.25), (1129.75, 2129.75)})
        cube = gdal.Open(paths["raster"])
        self.assertEqual(int(cube.ReadAsArray().sum()), 3)

    def test_stale_header_bounds_do_not_clip_records(self):
        _, expected, _ = self.run_export("baseline")
        # Header geometry alone changes; every stored point stays identical.
        with open(self.root / "source.laz", "r+b") as stream:
            stream.seek(179)  # Max X, LAS public header
            stream.write(struct.pack("<d", 1000.49996))
        result, actual, _ = self.run_export("stale")
        self.assertEqual(actual, expected)
        self.assertEqual(result["point_accounting"]["included"], 7)

    def test_fine_dtm_interior_hole_invalidates_whole_voxel(self):
        dtm = gdal.GetDriverByName("MEM").Create("", 5, 5, 1, gdal.GDT_Float64)
        dtm.SetGeoTransform((1000, .1, 0, 2000.5, 0, -.1))
        dtm.SetProjection(self.srs.ExportToWkt())
        values = np.zeros((5, 5)); values[2, 2] = -9999
        dtm.GetRasterBand(1).SetNoDataValue(-9999)
        dtm.GetRasterBand(1).WriteArray(values)
        grid = dict(width=1, height=1, transform=[1000,.5,0,2000.5,0,-.5],
                    crs_wkt=self.srs.ExportToWkt())
        mask, totals = voxel.make_validity_mask(str(self.root / "mask.tif"), grid,
            voxel.rectangle(1000,2000,1000.5,2000.5), voxel.RasterSampler(dtm))
        self.assertEqual(int(mask.ReadAsArray()[0,0]), 2)
        self.assertEqual(totals["invalid_dtm"], 1)

    def test_raster_directory_memory_is_admitted_before_creation(self):
        grid = dict(width=12800, height=12800)
        # This directory alone takes >30 MiB, outside GDAL's block cache.
        self.assertGreater(voxel.raster_memory_bytes(grid, 100), 30 * 1024**2)
        with self.assertRaisesRegex(ValueError, "raster directories"):
            voxel.check_raster_budget(grid, 100, 120 * 1024**2)
        voxel.check_raster_budget(grid, 100, 256 * 1024**2)

    def test_insufficient_budget_fails_before_allocating_cube(self):
        with self.assertRaisesRegex(ValueError, "memory_budget_exceeded"):
            self.run_export(memory=.01)


if __name__ == "__main__":
    unittest.main()
