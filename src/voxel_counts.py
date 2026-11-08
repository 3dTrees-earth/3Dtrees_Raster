#!/usr/bin/env python3
"""Full-source voxel counts using the existing laspy/GDAL stack and disk reduction.

Read record batches (also for COPC), never a spatial query with an unbounded
number of points. SQLite holds the sparse cube; GDAL writes one 2D block/band
at a time. Neither the cloud, DTM nor dense 3D cube is materialized in memory.
"""
from __future__ import annotations

import argparse
from collections import OrderedDict
from contextlib import closing
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sqlite3
import shutil
import subprocess
import tempfile
import time

import laspy
import numpy as np
from osgeo import gdal, ogr, osr
import pyproj
from raster_crs import horizontal_crs_wkt

RESOLUTION = 0.5
BLOCK = 128
NODATA = 4294967295
MAX_BANDS = 65535  # TIFF SamplesPerPixel limit; fail explicitly, never truncate.
FIXED_MEMORY_RESERVE = 192 * 1024**2
MINIMUM_MEMORY = 256 * 1024**2
FETCH_ROWS = 65536
POINT_SCALE = 0.001
UNTWINE = os.environ.get("RASTER_UNTWINE", "/opt/untwine/bin/untwine")


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class RasterSampler:
    """Nearest raster-cell lookup with a bounded cache, including outside/NoData."""

    def __init__(self, dataset, cache_blocks=32):
        self.dataset = dataset
        self.band = dataset.GetRasterBand(1)
        self.inverse = gdal.InvGeoTransform(dataset.GetGeoTransform())
        self.nodata = self.band.GetNoDataValue()
        self.cache = OrderedDict()
        self.cache_blocks = cache_blocks

    def sample(self, x, y):
        inv = self.inverse
        col = np.floor(inv[0] + inv[1] * x + inv[2] * y).astype(np.int64)
        row = np.floor(inv[3] + inv[4] * x + inv[5] * y).astype(np.int64)
        inside = ((col >= 0) & (col < self.dataset.RasterXSize) &
                  (row >= 0) & (row < self.dataset.RasterYSize))
        values = np.full(len(x), np.nan)
        indexes = np.flatnonzero(inside)
        if not len(indexes):
            return values
        blocks_across = (self.dataset.RasterXSize + BLOCK - 1) // BLOCK
        keys = (row[indexes] // BLOCK) * blocks_across + col[indexes] // BLOCK
        order = np.argsort(keys, kind="stable")
        indexes, keys = indexes[order], keys[order]
        edges = np.r_[0, np.flatnonzero(np.diff(keys)) + 1, len(keys)]
        for start, end in zip(edges[:-1], edges[1:]):
            key = int(keys[start])
            if key not in self.cache:
                xoff = (key % blocks_across) * BLOCK
                yoff = (key // blocks_across) * BLOCK
                self.cache[key] = self.band.ReadAsArray(
                    xoff, yoff, min(BLOCK, self.dataset.RasterXSize - xoff),
                    min(BLOCK, self.dataset.RasterYSize - yoff))
                if len(self.cache) > self.cache_blocks:
                    self.cache.popitem(last=False)
            self.cache.move_to_end(key)
            idx = indexes[start:end]
            values[idx] = self.cache[key][row[idx] % BLOCK, col[idx] % BLOCK]
        if self.nodata is not None:
            values[values == self.nodata] = np.nan
        return values


def create_raster(path, grid, bands, datatype, nodata=None):
    ds = gdal.GetDriverByName("GTiff").Create(
        str(path), grid["width"], grid["height"], bands, datatype,
        options=["TILED=YES", "BLOCKXSIZE=128", "BLOCKYSIZE=128",
                 "COMPRESS=DEFLATE", "BIGTIFF=YES", "INTERLEAVE=BAND"])
    ds.SetGeoTransform(grid["transform"])
    ds.SetProjection(grid["crs_wkt"])
    if nodata is not None:
        for band in range(1, bands + 1):
            ds.GetRasterBand(band).SetNoDataValue(nodata)
    return ds


def convert_raster_format(path, resampling):
    """Rewrite a finished GeoTIFF as a COG when RASTER_FORMAT=COG."""
    if os.environ.get("RASTER_FORMAT", "GTiff") != "COG":
        return
    partial = f"{path}.partial.tif"
    converted = gdal.Translate(partial, str(path), format="COG", creationOptions=[
        "COMPRESS=DEFLATE", "BIGTIFF=IF_SAFER", f"RESAMPLING={resampling}"])
    if converted is None:
        raise ValueError(f"could not write {Path(path).name} as a COG")
    converted = None
    os.replace(partial, path)


def grid_blocks(grid):
    across = (grid["width"] + BLOCK - 1) // BLOCK
    for y in range(0, grid["height"], BLOCK):
        for x in range(0, grid["width"], BLOCK):
            yield ((y // BLOCK) * across + x // BLOCK, x, y,
                   min(BLOCK, grid["width"] - x), min(BLOCK, grid["height"] - y))


def rectangle(xmin, ymin, xmax, ymax):
    ring = ogr.Geometry(ogr.wkbLinearRing)
    for x, y in [(xmin, ymin), (xmax, ymin), (xmax, ymax),
                 (xmin, ymax), (xmin, ymin)]:
        ring.AddPoint_2D(x, y)
    polygon = ogr.Geometry(ogr.wkbPolygon)
    polygon.AddGeometry(ring)
    return polygon


def make_validity_mask(path, grid, geometry, dtm):
    """1: entire cell in AOI and valid DTM; 0: AOI/boundary; 2: invalid terrain.

    Rasterize interiors quickly, then inspect only boundary cells geometrically.
    Sample every intersected terrain cell, including fine-resolution interior
    cells. Sampling is streamed over bounded blocks, never a dense DTM window.
    Counts and zero-valued cells use exactly this same domain mask.
    """
    output = create_raster(path, grid, 1, gdal.GDT_Byte)
    output.GetRasterBand(1).SetDescription("0=excluded;1=valid;2=invalid_DTM")
    vector = ogr.GetDriverByName("Memory").CreateDataSource("")
    srs = osr.SpatialReference(wkt=grid["crs_wkt"]) if grid["crs_wkt"] else None
    layer = vector.CreateLayer("usable", srs, geom_type=ogr.wkbUnknown)
    feature = ogr.Feature(layer.GetLayerDefn())
    feature.SetGeometry(geometry)
    layer.CreateFeature(feature)
    boundary_layer = vector.CreateLayer("boundary", srs, geom_type=ogr.wkbUnknown)
    edge = ogr.Feature(boundary_layer.GetLayerDefn())
    edge.SetGeometry(geometry.Boundary())
    boundary_layer.CreateFeature(edge)
    totals = {"valid": 0, "excluded": 0, "invalid_dtm": 0}
    for _, xoff, yoff, width, height in grid_blocks(grid):
        x0 = grid["transform"][0] + xoff * RESOLUTION
        y0 = grid["transform"][3] - yoff * RESOLUTION
        scratch = gdal.GetDriverByName("MEM").Create("", width, height, 2)
        scratch.SetGeoTransform((x0, RESOLUTION, 0, y0, 0, -RESOLUTION))
        scratch.SetProjection(grid["crs_wkt"])
        gdal.RasterizeLayer(scratch, [1], layer, burn_values=[1])
        gdal.RasterizeLayer(scratch, [2], boundary_layer, burn_values=[1],
                            options=["ALL_TOUCHED=TRUE"])
        valid = scratch.GetRasterBand(1).ReadAsArray().astype(bool)
        boundary = scratch.GetRasterBand(2).ReadAsArray().astype(bool)
        for row, col in np.argwhere(boundary):
            xmin, ymax = x0 + col * RESOLUTION, y0 - row * RESOLUTION
            valid[row, col] = geometry.Contains(rectangle(
                xmin, ymax - RESOLUTION, xmin + RESOLUTION, ymax))
        xx, yy = np.meshgrid(x0 + (np.arange(width) + 0.5) * RESOLUTION,
                             y0 - (np.arange(height) + 0.5) * RESOLUTION)
        terrain_valid = np.ones(xx.size, dtype=bool)
        # Move inside each edge to avoid sampling the neighboring cell.
        transform = dtm.dataset.GetGeoTransform()
        # Spacing no larger than a DTM pixel visits every intersected pixel,
        # regardless of lattice alignment. Corners alone miss interior holes.
        half = RESOLUTION / 2
        dxs = np.linspace(-half, half, max(1, math.ceil(RESOLUTION / transform[1])) + 1)
        dys = np.linspace(-half, half, max(1, math.ceil(RESOLUTION / -transform[5])) + 1)
        for dx in dxs:
            for dy in dys:
                x = np.nextafter(xx.ravel() + dx, xx.ravel())
                y = np.nextafter(yy.ravel() + dy, yy.ravel())
                terrain_valid &= np.isfinite(dtm.sample(x, y))
        terrain_valid = terrain_valid.reshape(height, width)
        mask = np.where(valid, np.where(terrain_valid, 1, 2), 0).astype(np.uint8)
        output.GetRasterBand(1).WriteArray(mask, xoff, yoff)
        for key, value in [("excluded", 0), ("valid", 1), ("invalid_dtm", 2)]:
            totals[key] += int(np.count_nonzero(mask == value))
    output.FlushCache()
    return output, totals


def accumulate(source, connection, grid, mask, dtm, chunk_points):
    connection.execute("""CREATE TABLE counts (
        block INTEGER, height INTEGER, cell INTEGER, n INTEGER NOT NULL,
        PRIMARY KEY(block,height,cell)) WITHOUT ROWID""")
    counters = dict(read=0, nonfinite=0, excluded_domain=0, invalid_dtm=0,
                    negative_height=0, included=0)
    across = (grid["width"] + BLOCK - 1) // BLOCK
    maximum = -1
    actual_peak = 0
    for points in source.chunk_iterator(chunk_points):
        counters["read"] += len(points)
        actual_peak = max(actual_peak, len(points))
        x, y, z = np.asarray(points.x), np.asarray(points.y), np.asarray(points.z)
        finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
        counters["nonfinite"] += int(np.count_nonzero(~finite))
        x, y, z = x[finite], y[finite], z[finite]
        # Own XY edges by increasing-coordinate [lower, upper) intervals,
        # independent of GDAL's top-down row convention at exact Y boundaries.
        col_all = np.floor((x - grid["transform"][0]) / RESOLUTION)
        iy_all = np.floor((y - (grid["transform"][3] - grid["height"] * RESOLUTION)) / RESOLUTION)
        domain = mask.sample(grid["transform"][0] + (col_all + 0.5) * RESOLUTION,
                             grid["transform"][3] - grid["height"] * RESOLUTION + (iy_all + 0.5) * RESOLUTION)
        terrain = dtm.sample(x, y)
        counters["excluded_domain"] += int(np.count_nonzero((domain == 0) | np.isnan(domain)))
        invalid = (domain == 2) | ((domain == 1) & ~np.isfinite(terrain))
        counters["invalid_dtm"] += int(np.count_nonzero(invalid))
        eligible = (domain == 1) & np.isfinite(terrain)
        x, y, h = x[eligible], y[eligible], z[eligible] - terrain[eligible]
        negative = h < 0
        counters["negative_height"] += int(np.count_nonzero(negative))
        x, y, h = x[~negative], y[~negative], h[~negative]
        counters["included"] += len(h)
        if not len(h):
            continue
        if not np.all(np.isfinite(h)) or np.max(h) >= MAX_BANDS * RESOLUTION:
            raise ValueError("normalized height exceeds GeoTIFF band capacity; review source outliers (no points truncated)")
        col = np.floor((x - grid["transform"][0]) / RESOLUTION).astype(np.int64)
        row = grid["height"] - 1 - np.floor((y - (grid["transform"][3] - grid["height"] * RESOLUTION)) / RESOLUTION).astype(np.int64)
        iz = np.floor(h / RESOLUTION).astype(np.int64)
        maximum = max(maximum, int(iz.max()))
        # Pack the integer address instead of sorting Python/structured XYZ
        # tuples. A single int64 sort keeps dense batches fast and small.
        block = (row // BLOCK) * across + col // BLOCK
        keys = ((block * MAX_BANDS + iz) * (BLOCK * BLOCK) +
                (row % BLOCK) * BLOCK + col % BLOCK)
        unique, counts = np.unique(keys, return_counts=True)
        # Split the sorted addresses in NumPy and hand SQLite plain Python
        # lists; per-voxel NumPy scalar arithmetic dominated dense batches.
        connection.executemany(
            "INSERT INTO counts VALUES (?,?,?,?) ON CONFLICT(block,height,cell) DO UPDATE SET n=n+excluded.n",
            zip((unique // (MAX_BANDS * BLOCK * BLOCK)).tolist(),
                ((unique // (BLOCK * BLOCK)) % MAX_BANDS).tolist(),
                (unique % (BLOCK * BLOCK)).tolist(), counts.tolist()))
    connection.commit()
    if counters["read"] != source.header.point_count:
        raise ValueError("source header/stream point count mismatch")
    if sum(v for k, v in counters.items() if k != "read") != counters["read"]:
        raise ValueError("point accounting is not exhaustive")
    total, maximum_count, occupied = connection.execute(
        "SELECT coalesce(sum(n),0), coalesce(max(n),0), count(*) FROM counts").fetchone()
    if total != counters["included"] or maximum_count >= NODATA:
        raise ValueError("count conservation or UInt32 overflow check failed")
    return max(1, maximum + 1), counters, occupied, actual_peak


def voxel_point_header(grid, crs_wkt):
    """LAS 1.4 header for one point per occupied voxel (format 6, as COPC needs)."""
    header = laspy.LasHeader(point_format=6, version="1.4")
    header.add_extra_dims([
        laspy.ExtraBytesParams(name="point_count", type=np.uint32,
                               description="source points in the voxel"),
        laspy.ExtraBytesParams(name="height_above_dtm", type=np.float32,
                               description="voxel centre above the DTM (m)"),
    ])
    header.offsets = np.array([grid["transform"][0],
                               grid["transform"][3] - grid["height"] * RESOLUTION, 0.0])
    header.scales = np.array([POINT_SCALE] * 3)
    if crs_wkt:
        header.add_crs(pyproj.CRS.from_wkt(crs_wkt))
    return header


def write_copc(points_path, copc_path, temp_dir):
    """Index the voxel points as one COPC; untwine keeps the extra dimensions."""
    partial = Path(f"{copc_path}.partial.copc.laz")
    partial.unlink(missing_ok=True)
    environment = dict(os.environ, PROJ_DATA=str(Path(UNTWINE).parents[1] / "share" / "proj"))
    environment.pop("PROJ_LIB", None)
    completed = subprocess.run(
        [UNTWINE, "--files", str(points_path), "--output_dir", str(partial),
         "--temp_dir", str(temp_dir)],
        env=environment, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    if completed.returncode != 0 or not partial.is_file():
        raise ValueError(f"could not write the voxel COPC: {completed.stderr.strip()[-500:]}")
    os.replace(partial, copc_path)


def untwine_version():
    completed = subprocess.run([UNTWINE, "--version"], stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True)
    return completed.stdout.strip().split("(")[-1].rstrip(")") if completed.returncode == 0 else "unknown"


def check_copc(copc_path, occupied, included):
    """The COPC holds every occupied voxel once and conserves the point counts."""
    voxels = total = 0
    with laspy.open(copc_path) as reader:
        if reader.header.point_count != occupied:
            raise ValueError("voxel COPC point count differs from the occupied voxels")
        for points in reader.chunk_iterator(1_000_000):
            voxels += len(points)
            total += int(np.asarray(points["point_count"], dtype=np.int64).sum())
    if voxels != occupied or total != included:
        raise ValueError("voxel COPC lost voxels or counts")


def write_counts(connection, paths, grid, bands, mask, terrain, points_path, crs_wkt):
    raster = create_raster(paths["raster"], grid, bands, gdal.GDT_UInt32, NODATA)
    for k in range(bands):
        band = raster.GetRasterBand(k + 1)
        band.SetDescription(f"points_{k * RESOLUTION:g}_{(k + 1) * RESOLUTION:g}m_AGL")
        band.SetUnitType("count")
        band.SetMetadata({"height_lower_m": str(k * RESOLUTION),
                          "height_upper_m": str((k + 1) * RESOLUTION),
                          "height_reference": "DTM", "interval": "[lower,upper)"})
    raster.SetMetadata({"metric": "point_count", "voxel_size_m": f"{RESOLUTION:g}",
                        "height_reference": "DTM", "zero": "no recorded eligible points"})
    exported = 0
    with open(paths["csv"], "w", newline="", encoding="utf-8") as output, \
            laspy.open(points_path, mode="w", header=voxel_point_header(grid, crs_wkt)) as points:
        writer = csv.writer(output)
        writer.writerow(["x_center", "y_center", "height_lower_m", "height_upper_m", "point_count"])
        for block, xoff, yoff, width, height in grid_blocks(grid):
            valid = mask.GetRasterBand(1).ReadAsArray(xoff, yoff, width, height) == 1
            empty = np.where(valid, 0, NODATA).astype(np.uint32)
            band_index, values = 0, None
            cursor = connection.execute(
                "SELECT height,cell,n FROM counts WHERE block=? ORDER BY height,cell", (block,))
            # Bounded pieces of one block, converted to arrays: no per-voxel
            # Python loop and no unbounded fetch for dense blocks.
            while piece := cursor.fetchmany(FETCH_ROWS):
                rows = np.array(piece, dtype=np.int64)
                k, n = rows[:, 0], rows[:, 2]
                row, col = rows[:, 1] // BLOCK, rows[:, 1] % BLOCK
                if not valid[row, col].all():
                    raise ValueError("occupied voxel outside valid mask")
                # Same float64 operation order as the former scalar writer, so
                # csv.writer formats identical Python floats.
                x_center = grid["transform"][0] + (xoff + col + 0.5) * RESOLUTION
                y_center = grid["transform"][3] - (yoff + row + 0.5) * RESOLUTION
                writer.writerows(zip(x_center.tolist(), y_center.tolist(),
                    (k * RESOLUTION).tolist(), ((k + 1) * RESOLUTION).tolist(), n.tolist()))
                # One point per voxel at its band centre above the DTM cell that
                # contains the voxel column centre.
                ground = terrain.sample(x_center, y_center)
                if not np.all(np.isfinite(ground)):
                    raise ValueError("occupied voxel column without a valid DTM height")
                record = laspy.ScaleAwarePointRecord.zeros(len(n), header=points.header)
                height = (k + 0.5) * RESOLUTION
                record.x, record.y, record.z = x_center, y_center, ground + height
                record["point_count"] = n.astype(np.uint32)
                record["height_above_dtm"] = height.astype(np.float32)
                # Standard fields for viewers that ignore extra bytes (Potree 1.8):
                # the count in intensity and the 0-based height band in user data,
                # both capped at their field's maximum.
                record.intensity = np.minimum(n, np.iinfo(np.uint16).max).astype(np.uint16)
                record.user_data = np.minimum(k, np.iinfo(np.uint8).max).astype(np.uint8)
                points.write_points(record)
                exported += int(n.sum())
                edges = np.r_[0, np.flatnonzero(np.diff(k)) + 1, len(k)]
                for first, last in zip(edges[:-1], edges[1:]):
                    while band_index < k[first]:
                        raster.GetRasterBand(band_index + 1).WriteArray(
                            empty if values is None else values, xoff, yoff)
                        band_index, values = band_index + 1, None
                    if values is None:
                        values = empty.copy()
                    values[row[first:last], col[first:last]] = n[first:last]
            while band_index < bands:
                raster.GetRasterBand(band_index + 1).WriteArray(
                    empty if values is None else values, xoff, yoff)
                band_index, values = band_index + 1, None
    raster.FlushCache()
    return exported


def raster_memory_bytes(grid, bands):
    """Conservative libtiff/GDAL directory allowance, outside the block cache.

    BigTIFF needs offset and byte-count arrays per tile per band. Reserve twice
    their 16-byte payload plus per-band objects/metadata; the 192 MiB fixed reserve
    in check_raster_budget covers GDAL/SQLite caches and the bounded samplers.
    """
    tiles = math.ceil(grid["width"] / BLOCK) * math.ceil(grid["height"] / BLOCK)
    return tiles * bands * 32 + bands * 4096


def check_raster_budget(grid, bands, budget, other_raster_bytes=0):
    estimate = raster_memory_bytes(grid, bands)
    required = FIXED_MEMORY_RESERVE + estimate + other_raster_bytes
    if required > budget:
        raise ValueError(
            f"memory_budget_exceeded: raster directories require {required} bytes "
            f"including caches, available {int(budget)}; reduce AOI or increase budget")
    return dict(directory_bytes=estimate, fixed_cache_reserve_bytes=FIXED_MEMORY_RESERVE,
                required_bytes=required, available_bytes=int(budget))


def export(config):
    # One voxel edge length for XY cells and height bands, fixed per run.
    global RESOLUTION
    RESOLUTION = float(config.get("resolution", 0.5))
    if not math.isfinite(RESOLUTION) or RESOLUTION <= 0:
        raise ValueError("voxel-count resolution must be finite and positive")
    started = time.monotonic()
    gdal.UseExceptions()
    ogr.UseExceptions()
    osr.UseExceptions()
    gdal.SetCacheMax(32 * 1024 * 1024)
    dtm_ds = gdal.Open(config["dtm"])
    transform = dtm_ds.GetGeoTransform()
    if transform[2] != 0 or transform[4] != 0 or transform[1] <= 0 or transform[5] >= 0:
        raise ValueError("voxel counts require a north-up DTM")
    # Like the other products, a source without CRS gives products without CRS
    # (units taken as metres); a CRS that is present must be projected in metres.
    dtm_wkt = dtm_ds.GetProjection()
    crs = osr.SpatialReference(wkt=dtm_wkt) if dtm_wkt else None
    if crs is not None and crs.IsLocal():
        crs = None
    if crs is not None and (not crs.IsProjected() or abs(crs.GetLinearUnits() - 1) > 1e-9):
        raise ValueError("voxel counts require a projected CRS in metres when the source has a CRS")
    vector = ogr.Open(config["aoi"])
    # Keep the Feature alive while cloning its borrowed geometry reference.
    feature = vector.GetLayer(0).GetNextFeature()
    geometry = feature.GetGeometryRef().Clone()
    if geometry.IsEmpty() or not geometry.IsValid():
        raise ValueError("invalid voxel-count AOI")
    budget = float(config.get("memory_available_gib", 1)) * 1024**3
    if budget < MINIMUM_MEMORY:
        raise ValueError("memory_budget_exceeded: voxel export needs at least 256 MiB available")
    paths = config["paths"]
    with laspy.open(config["point_cloud"]) as source:
        # The product domain is the usable AOI, independent of potentially
        # stale/pre-quantization LAS header extrema. All source records stream.
        extent = geometry.GetEnvelope()
        xmin = math.floor(extent[0] / RESOLUTION) * RESOLUTION
        ymin = math.floor(extent[2] / RESOLUTION) * RESOLUTION
        xmax = math.ceil(extent[1] / RESOLUTION) * RESOLUTION
        ymax = math.ceil(extent[3] / RESOLUTION) * RESOLUTION
        if xmax <= xmin or ymax <= ymin:
            raise ValueError("AOI and source extent do not intersect")
        grid = dict(width=int(round((xmax-xmin)/RESOLUTION)),
                    height=int(round((ymax-ymin)/RESOLUTION)),
                    transform=[xmin, RESOLUTION, 0, ymax, 0, -RESOLUTION],
                    crs_wkt=horizontal_crs_wkt(dtm_wkt) if crs is not None else "", lattice_origin=[0, 0, 0])
        block_count = ((grid["width"] + BLOCK - 1) // BLOCK) * ((grid["height"] + BLOCK - 1) // BLOCK)
        if block_count > np.iinfo(np.int64).max // (MAX_BANDS * BLOCK * BLOCK):
            raise ValueError("grid exceeds integer voxel-address capacity")
        mask_bookkeeping = raster_memory_bytes(grid, 1)
        check_raster_budget(grid, 1, budget)
        # Reserve GDAL/SQLite caches, raster directories and batch temporaries.
        # Dataset 226 batches measured ~670 B of RSS per point; admit four
        # records plus 400 B of float/int64 temporaries per point.
        chunk_points = min(int(config.get("chunk_points", 2_000_000)),
                           int((budget - FIXED_MEMORY_RESERVE - mask_bookkeeping) / (4 * source.header.point_format.size + 400)))
        if chunk_points < 1:
            raise ValueError("chunk_points must be positive within the memory budget")
        terrain = RasterSampler(dtm_ds)
        mask_ds, domain = make_validity_mask(paths["mask"], grid, geometry, terrain)
        mask = RasterSampler(mask_ds)
        with tempfile.TemporaryDirectory(prefix="voxel_counts_", dir=config["temp_dir"]) as scratch:
            # closing(): the connection context manager only commits; an open
            # database on NFS keeps the scratch directory from being removed.
            with closing(sqlite3.connect(str(Path(scratch) / "counts.sqlite"))) as connection, connection:
                # Scratch database: a crash discards the whole job, so skip the
                # rollback journal and fsyncs that made dense runs I/O bound.
                connection.execute("PRAGMA journal_mode=OFF")
                connection.execute("PRAGMA synchronous=OFF")
                connection.execute("PRAGMA locking_mode=EXCLUSIVE")
                connection.execute("PRAGMA cache_size=-65536")
                connection.execute("PRAGMA temp_store=FILE")
                bands, counters, occupied, peak = accumulate(source, connection, grid, mask, terrain, chunk_points)
                raster_resources = check_raster_budget(grid, bands, budget, mask_bookkeeping)
                # Admission uses uncompressed tiled size, conservative CSV rows and
                # uncompressed voxel points (scratch LAS plus untwine's spill);
                # compression is data dependent and must not be assumed here.
                required_disk = block_count * BLOCK * BLOCK * 4 * bands + occupied * (128 + 3 * 64)
                if shutil.disk_usage(Path(paths["raster"]).parent).free < required_disk:
                    raise ValueError("disk_budget_exceeded: insufficient free space for voxel products")
                points_path = Path(scratch) / "voxel_points.las"
                point_crs = source.header.parse_crs() if crs is not None else None
                point_crs_wkt = point_crs.to_wkt() if point_crs else grid["crs_wkt"]
                if write_counts(connection, paths, grid, bands, mask_ds, terrain,
                                points_path, point_crs_wkt) != counters["included"]:
                    raise ValueError("CSV/raster export lost counts")
                # Untwine cannot index an empty cloud; without occupied voxels
                # there is no COPC and the metadata says so.
                if occupied:
                    untwine_scratch = Path(scratch) / "untwine"
                    untwine_scratch.mkdir()
                    write_copc(points_path, paths["copc"], untwine_scratch)
                    check_copc(paths["copc"], occupied, counters["included"])
                disk_bytes = sum(p.stat().st_size for p in Path(scratch).iterdir() if p.is_file())
        mask = mask_ds = None
        convert_raster_format(paths["raster"], "AVERAGE")
        convert_raster_format(paths["mask"], "NEAREST")
        source_crs = source.header.parse_crs()
        metadata = dict(
            schema_version=1, metric="point_count", units="count", voxel_size_m=[RESOLUTION]*3,
            grid=grid, band_count=bands, height_reference="metres above DTM",
            height_interval="[lower,upper)", band_bounds_formula=f"band b (1-based): [(b-1)*{RESOLUTION:g},b*{RESOLUTION:g})",
            normalization="source Z minus containing DTM cell; no height re-quantization or gap filling",
            negative_height_policy="exclude and count; never clamp to zero",
            outlier_policy="exclude nonfinite only; retain all finite nonnegative heights; fail if TIFF band capacity exceeded",
            boundary_policy=f"only whole {RESOLUTION:g} m cells contained in usable AOI; half-open XY grid",
            zero="no recorded eligible points; not proof of vegetation absence or sampling coverage",
            nodata=NODATA, validity_mask={"0": "excluded AOI/boundary", "1": "valid", "2": "invalid DTM"},
            source=dict(name=config.get("source_name", Path(config["point_cloud"]).name),
                        sha256=config.get("source_sha256") or sha256_file(config["point_cloud"]),
                        point_count=int(source.header.point_count),
                        crs_wkt=source_crs.to_wkt() if source_crs else None,
                        read_mode="all source records; no subsampling, LOD restriction or coordinate deduplication"),
            dtm=dict(name=Path(config["dtm"]).name, sha256=sha256_file(config["dtm"]),
                     source_crs_wkt=dtm_ds.GetProjection(),
                     resolution_m=[transform[1], -transform[5]]),
            point_accounting=counters, domain_cells=domain, occupied_voxels=occupied,
            chunk_points=chunk_points, peak_batch_points=peak, temporary_disk_bytes=disk_bytes,
            raster_resource_estimate=raster_resources,
            grid_extent_source="usable AOI; independent of LAS header bounds",
            elapsed_seconds=time.monotonic()-started,
            voxel_points=None if not occupied else dict(
                file=Path(paths["copc"]).name, format="COPC 1.0 (LAS 1.4, point format 6)",
                points="one per occupied voxel; same voxels and counts as the CSV",
                xy="voxel column centre", z="DTM cell containing the column centre plus the band centre",
                attributes=dict(point_count="source points in the voxel (uint32)",
                                height_above_dtm="band centre above the DTM in metres (float32)",
                                intensity="point_count capped at 65535",
                                user_data="0-based height band (height_above_dtm / voxel size, floored), capped at 255"),
                crs=("source point-cloud CRS" if point_crs else "horizontal DTM CRS") if grid["crs_wkt"]
                    else "none: the source has no CRS; coordinates as in the source, units taken as metres",
                writer=f"untwine {untwine_version()}"),
            artifacts={key: Path(value).name for key, value in paths.items()
                       if key != "copc" or occupied})
        if config.get("dtm_denoising"):
            metadata["dtm"]["denoising"] = config["dtm_denoising"]
            metadata["noise_mask_applied_to_counts"] = False
        Path(paths["metadata"]).write_text(json.dumps(metadata, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    return metadata


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    args = parser.parse_args()
    export(json.loads(args.config.read_text(encoding="utf-8")))
