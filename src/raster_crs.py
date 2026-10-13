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


if __name__ == "__main__":
    print(horizontal_crs_wkt(sys.argv[1]))
