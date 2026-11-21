"""Resolve the horizontal reference for self-contained two-dimensional rasters.

The source may have a projected 3D or compound CRS. GeoTIFF publication must
not rely on PAM sidecars for that reference. This changes metadata only: source
coordinates and raster elevations/counts are never transformed.
"""
import sys

from osgeo import osr


def horizontal_crs_wkt(wkt):
    osr.UseExceptions()
    reference = osr.SpatialReference(wkt=wkt)
    if reference.GetAxesCount() <= 2:
        return wkt
    reference.DemoteTo2D()
    if reference.GetAxesCount() != 2:
        raise ValueError("cannot resolve a two-dimensional horizontal raster CRS")
    return reference.ExportToWkt(["FORMAT=WKT2_2019"])


def same_horizontal_crs(a, b):
    """Whether two CRS have the same horizontal reference (axis order ignored)."""
    osr.UseExceptions()
    first = osr.SpatialReference(wkt=horizontal_crs_wkt(osr.SpatialReference(wkt=a).ExportToWkt(["FORMAT=WKT2_2019"])))
    second = osr.SpatialReference(wkt=horizontal_crs_wkt(osr.SpatialReference(wkt=b).ExportToWkt(["FORMAT=WKT2_2019"])))
    return bool(first.IsSame(second, ["IGNORE_DATA_AXIS_TO_SRS_AXIS_MAPPING=YES",
                                      "CRITERION=EQUIVALENT_EXCEPT_AXIS_ORDER_GEOGCRS"]))


if __name__ == "__main__":
    if sys.argv[1] == "--same":
        print("true" if same_horizontal_crs(sys.argv[2], sys.argv[3]) else "false")
    else:
        print(horizontal_crs_wkt(sys.argv[1]))
