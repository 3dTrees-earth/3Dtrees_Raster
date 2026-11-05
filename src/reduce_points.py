#!/usr/bin/env python3
"""Keep only the original points that decide the raster products.

Very dense clouds spend hours in repeated lidR reads, yet the results need far
fewer points. One streaming pass per mode keeps them, unchanged:

* `thin`: the lowest point per 5 cm cell, on lidR::lowest()'s grid, for the ground
  filter of clouds whose DTM windows are too large for it.
* `voxels`: for the CHM, DSM, segments and tiles, which read occupied voxels,
  instance voxels, highest points and height layers. Per voxel half-layer and
  tile one point, per instance the first point of each voxel (in source order)
  and its highest points, per voxel column the highest absolute point. Any
  subset holding these gives the same products as the whole cloud. Points too
  close to a grid edge to classify safely, or on an empty DTM cell, are all kept.

Output records stay in source order with their stored coordinates.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import laspy
import numpy as np
import polars as pl
from osgeo import gdal

from project_instance_dimension import extra_dimension, extra_params

BATCH_POINTS = 20_000_000
LOWEST_CELL_M = 0.05
EDGE_TOLERANCE = 1e-7  # in grid cells, far below any stored coordinate step
ORDER = "__order"
FILTERS = laspy.DecompressionSelection


def batches(path: Path, extra_bytes: bool):
    selection = FILTERS.XY_RETURNS_CHANNEL | FILTERS.Z
    if extra_bytes:
        selection |= FILTERS.ALL_EXTRA_BYTES
    with laspy.open(path, laz_backend=laspy.LazBackend.LazrsParallel,
                    decompression_selection=selection) as reader:
        yield from reader.chunk_iterator(BATCH_POINTS)


def table(points, first: int, dimensions: list[str], order_dimension: str | None) -> pl.DataFrame:
    """Stored coordinates and dimensions, with each record's source position."""
    order = (np.asarray(points.array[order_dimension]).astype(np.int64) if order_dimension
             else np.arange(first, first + len(points), dtype=np.int64))
    columns = {"X": points.X, "Y": points.Y, "Z": points.Z, ORDER: order}
    columns.update({name: np.asarray(points.array[name]) for name in dimensions})
    return pl.DataFrame(columns)


def rejected_ids(directory: Path | None) -> np.ndarray | None:
    if directory is None:
        return None
    parts = [np.fromfile(path, "<f8") for path in sorted(directory.glob("*.bin"))]
    return np.sort(np.concatenate(parts)).astype(np.int64) if parts else None


def accepted(order: np.ndarray, rejected: np.ndarray | None) -> np.ndarray:
    if rejected is None or not len(rejected):
        return np.ones(len(order), bool)
    position = np.minimum(np.searchsorted(rejected, order), len(rejected) - 1)
    return rejected[position] != order


def near_edge(values: np.ndarray) -> np.ndarray:
    return np.abs(values - np.rint(values)) < EDGE_TOLERANCE


def write(path: Path, source: laspy.LasHeader, selected: pl.DataFrame,
          dimensions: list[str]) -> None:
    header = laspy.LasHeader(point_format=0, version=str(source.version))
    header.scales = source.scales.copy()
    header.offsets = source.offsets.copy()
    header.generating_software = "3dtrees-foreststructure-julia-memory-safe"
    crs = source.parse_crs()
    if crs is not None:
        header.add_crs(crs, keep_compatibility=False)
    for name in dimensions:
        header.add_extra_dim(extra_params(extra_dimension(source, name)))
    selected = selected.sort(ORDER)
    points = laspy.ScaleAwarePointRecord.zeros(selected.height, header=header)
    for name in ("X", "Y", "Z"):
        setattr(points, name, selected[name].to_numpy())
    for name in dimensions:
        points.array[name] = selected[name].to_numpy()
    with laspy.open(path, mode="w", header=header, do_compress=True) as writer:
        writer.write_points(points)


def source_header(path: Path) -> laspy.LasHeader:
    with laspy.open(path) as reader:
        return reader.header


def thin(args) -> dict:
    """lidR::lowest(res = 0.05): per cell the lowest point, the first of equals."""
    header = source_header(args.input)
    keep = [args.order_dimension] if args.order_dimension else []
    rejected = rejected_ids(args.rejected_directory)
    kept, seen = None, 0
    for points in batches(args.input, bool(keep)):
        frame = table(points, seen, keep, args.order_dimension)
        seen += len(points)
        valid = accepted(frame[ORDER].to_numpy(), rejected)
        cells = frame.filter(pl.Series(valid)).with_columns(
            cell_x=pl.Series(np.floor(np.asarray(points.x)[valid] / LOWEST_CELL_M)),
            cell_y=pl.Series(np.floor(np.asarray(points.y)[valid] / LOWEST_CELL_M)))
        kept = cells if kept is None else pl.concat([kept, cells])
        kept = kept.sort(["Z", ORDER]).unique(["cell_x", "cell_y"], keep="first", maintain_order=True)
    write(args.output, header, kept, keep)
    return {"source_points": seen, "kept_points": kept.height}


def voxel_columns(frame: pl.DataFrame, x, y, height, resolution, x_edges, y_edges) -> pl.DataFrame:
    return frame.with_columns(
        height=pl.Series(height),
        voxel_x=pl.Series(np.floor(x / resolution).astype(np.int64)),
        voxel_y=pl.Series(np.floor(y / resolution).astype(np.int64)),
        voxel_z=pl.Series(np.floor(height / (resolution / 2)).astype(np.int64)),
        tile_x=pl.Series(np.searchsorted(x_edges, x, side="right")),
        tile_y=pl.Series(np.searchsorted(y_edges, y, side="right")))


def decisive(frame: pl.DataFrame, dimensions: list[str], step: float) -> pl.DataFrame:
    """The rows that fix every voxel set, apex and DSM cell of this table."""
    voxel = ["voxel_x", "voxel_y", "voxel_z"]
    ordered = frame.sort(ORDER)
    parts = [ordered.unique(voxel + ["tile_x", "tile_y"], keep="first", maintain_order=True),
             frame.sort(["Z", ORDER], descending=[True, False]).unique(
                 ["voxel_x", "voxel_y"], keep="first", maintain_order=True)]
    for name in dimensions:
        trees = ordered.filter(pl.col(name) > 0)
        parts.append(trees.unique([name] + voxel, keep="first", maintain_order=True))
        parts.append(trees.filter(pl.col("height") >= pl.col("height").max().over(name) - 2 * step))
    return pl.concat(parts).unique(ORDER, keep="first")


def extremes(frame: pl.DataFrame) -> pl.DataFrame:
    """Rows on the cloud's bounding box, so the reduced header keeps its extent."""
    rows = [frame.sort(name, descending=descending).head(1)
            for name in ("X", "Y", "Z") for descending in (False, True)]
    return pl.concat(rows)


def voxels(args) -> dict:
    header = source_header(args.input)
    scale_x, scale_y, scale_z = header.scales
    dtm = gdal.Open(str(args.dtm))
    left, cell_x, _, top, _, cell_y = dtm.GetGeoTransform()
    band = dtm.GetRasterBand(1)
    terrain = band.ReadAsArray().astype(float)
    if band.GetNoDataValue() is not None:
        terrain[terrain == band.GetNoDataValue()] = np.nan
    rows, columns = terrain.shape
    edges = json.loads(args.tile_edges.read_text())
    x_edges, y_edges = np.array(edges["x"]), np.array(edges["y"])
    resolution = args.voxel_resolution
    names = args.dimension + ([args.order_dimension] if args.order_dimension else [])
    rejected = rejected_ids(args.rejected_directory)
    kept, unsure, bounds, seen = None, [], None, 0
    for points in batches(args.input, True):
        x, y, z = np.asarray(points.x), np.asarray(points.y), np.asarray(points.z)
        frame = table(points, seen, names, args.order_dimension)
        seen += len(points)
        box = extremes(frame)
        bounds = box if bounds is None else extremes(pl.concat([bounds, box]))
        valid = accepted(frame[ORDER].to_numpy(), rejected)
        # terra::cellFromXY: the last column and row include the far raster edge.
        column_position = (x - left) / cell_x
        row_position = (y - top) / cell_y
        column = np.minimum(np.floor(column_position), columns - 1).astype(np.int64)
        row = np.minimum(np.floor(row_position), rows - 1).astype(np.int64)
        inside = (column >= 0) & (row >= 0) & (column_position <= columns) & (row_position <= rows)
        ground = np.full(len(x), np.nan)
        ground[inside] = terrain[row[inside], column[inside]]
        # lidR stores heights on the Z lattice with a zero offset.
        steps = (z - ground) / scale_z
        height = np.trunc(steps + np.copysign(0.5, steps)) * scale_z
        known = np.isfinite(height)
        edge = (near_edge(x / resolution) | near_edge(y / resolution)
                | near_edge(height / (resolution / 2)) | near_edge(column_position)
                | near_edge(row_position))
        for values, grid in ((x, x_edges), (y, y_edges)):
            nearest = np.minimum(np.searchsorted(grid, values), len(grid) - 1)
            edge |= np.abs(values - grid[nearest]) < EDGE_TOLERANCE * resolution
            edge |= np.abs(values - grid[np.maximum(nearest - 1, 0)]) < EDGE_TOLERANCE * resolution
        usable = known & (height >= -2 * scale_z) & (height <= args.maximum_height + 2 * scale_z)
        unsure.append(frame.filter(pl.Series(valid & (~known | (usable & edge)))))
        chosen = valid & usable & ~edge
        batch = voxel_columns(frame.filter(pl.Series(chosen)), x[chosen], y[chosen], height[chosen],
                              resolution, x_edges, y_edges)
        batch = decisive(batch, args.dimension, scale_z)
        kept = batch if kept is None else decisive(pl.concat([kept, batch]), args.dimension, scale_z)
    columns_out = ["X", "Y", "Z", ORDER] + names
    selected = pl.concat([kept.select(columns_out), bounds.select(columns_out)]
                         + [part.select(columns_out) for part in unsure])
    selected = selected.unique(ORDER, keep="first")
    write(args.output, header, selected, names)
    return {"source_points": seen, "kept_points": selected.height,
            "unsure_points": sum(part.height for part in unsure)}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    for name in ("thin", "voxels"):
        mode = modes.add_parser(name)
        mode.add_argument("--input", required=True, type=Path)
        mode.add_argument("--output", required=True, type=Path)
        mode.add_argument("--report", required=True, type=Path)
        mode.add_argument("--order-dimension")
        mode.add_argument("--rejected-directory", type=Path)
    mode.add_argument("--dtm", required=True, type=Path)
    mode.add_argument("--tile-edges", required=True, type=Path)
    mode.add_argument("--voxel-resolution", required=True, type=float)
    mode.add_argument("--maximum-height", required=True, type=float)
    mode.add_argument("--dimension", action="append", default=[])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = (thin if args.mode == "thin" else voxels)(args)
    args.report.write_text(json.dumps(report) + "\n")


if __name__ == "__main__":
    main()
