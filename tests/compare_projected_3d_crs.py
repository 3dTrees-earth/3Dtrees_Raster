"""Two identical clouds differ only by whether their source CRS has a Z axis."""
import filecmp
import json
from pathlib import Path
import numpy as np
from osgeo import gdal

gdal.UseExceptions()
gdal.SetConfigOption("GDAL_PAM_ENABLED", "NO")
root = Path("/data")
for path in (root / "two").glob("foreststructure_*"):
    other = root / "three" / path.name
    if path.suffix == ".csv" and "performance" not in path.name:
        assert filecmp.cmp(path, other, shallow=False), path.name
    elif path.suffix == ".geojson":
        assert json.loads(path.read_text())["features"] == json.loads(other.read_text())["features"], path.name
    elif path.suffix == ".tif":
        left, right = gdal.Open(str(path)), gdal.Open(str(other))
        assert left.GetGeoTransform() == right.GetGeoTransform(), path.name
        assert left.RasterCount == right.RasterCount, path.name
        assert right.GetSpatialRef() and right.GetSpatialRef().GetAxesCount() == 2, path.name
        assert left.GetSpatialRef().IsSame(right.GetSpatialRef()), path.name
        for band in range(1, left.RasterCount + 1):
            np.testing.assert_array_equal(left.GetRasterBand(band).ReadAsArray(), right.GetRasterBand(band).ReadAsArray())
provenance = json.loads((root / "three/foreststructure_julia_memory_safe_run.json").read_text())
assert 'CS[Cartesian,3]' in provenance['source_crs_wkt']
assert 'CS[Cartesian,2]' in provenance['raster_crs_wkt']
