"""Fold CHM/DSM cells beyond the point cloud's upper X/Y bound into the last cell.

Cells own [lower, upper). Points lying exactly on the cloud's maximum X or Y open a
column or row beyond the cloud, empty but for them (3DT-3028). CHM and DSM cells
hold a maximum, so moving those points into the last cell inside the cloud is the
maximum over both cells. The raster is cropped to the cloud; only the last column
and the top row are rewritten, every other cell keeps its value.
"""
import os
import sys

import numpy as np
from osgeo import gdal

gdal.UseExceptions()


def fold(path, xmax, ymax):
    source = gdal.Open(path)
    left, size_x, _, top, _, size_y = source.GetGeoTransform()
    width, height = source.RasterXSize, source.RasterYSize
    band = source.GetRasterBand(1)
    nodata = band.GetNoDataValue()
    # Keep the columns that start before xmax and the rows that end after ymax. A bound
    # within a millimetre past a cell edge counts as on it: COPC writers may round the
    # header bounds up by a fraction of the coordinate scale.
    tolerance = 1e-3
    keep_cols = max(1, min(width, int(np.ceil((xmax - tolerance - left) / size_x))))
    drop_rows = max(0, min(height - 1, int(np.floor((top - ymax + tolerance) / -size_y))))
    if keep_cols == width and drop_rows == 0:
        return False

    def read(x, y, w, h):
        values = band.ReadAsArray(x, y, w, h)
        floats = values.astype(np.float64)
        if nodata is not None:
            floats[values == nodata] = np.nan
        return values.dtype, floats

    def restore(dtype, values):
        if nodata is not None:
            values = np.where(np.isnan(values), nodata, values)
        return values.astype(dtype)

    with np.errstate(all="ignore"):
        # Last kept column: the maximum over it and the dropped columns, row by row.
        dtype, strip = read(keep_cols - 1, 0, width - keep_cols + 1, height)
        column = np.nanmax(strip, axis=1)
        # Top kept row: the maximum over it and the dropped rows, the last column included.
        _, rows = read(0, 0, keep_cols, drop_rows + 1)
        rows[:, keep_cols - 1] = column[:drop_rows + 1]
        row = np.nanmax(rows, axis=0)
    partial = path + ".fold.tif"
    gdal.Translate(partial, source, srcWin=[0, drop_rows, keep_cols, height - drop_rows],
                   creationOptions=["TILED=YES", "COMPRESS=DEFLATE", "BIGTIFF=IF_SAFER"])
    source = None
    target = gdal.Open(partial, gdal.GA_Update)
    out = target.GetRasterBand(1)
    out.WriteArray(restore(dtype, column[drop_rows:]).reshape(-1, 1), keep_cols - 1, 0)
    out.WriteArray(restore(dtype, row).reshape(1, -1), 0, 0)
    target = None
    os.replace(partial, path)
    return True


if __name__ == "__main__":
    print("folded" if fold(sys.argv[1], float(sys.argv[2]), float(sys.argv[3])) else "unchanged")
