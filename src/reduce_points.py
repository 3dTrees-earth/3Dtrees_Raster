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

Output records keep their stored coordinates. A cloud with a source-order key
leaves either mode as a COPC, so windows read it through the spatial index and
memory is estimated per window; without a key it stays a LAZ in source order.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing
import os
from pathlib import Path
import tempfile

import laspy
import numpy as np
import polars as pl
from osgeo import gdal

from project_instance_dimension import extra_dimension, extra_params
from voxel_counts import copc_chunks, write_copc

BATCH_POINTS = 5_000_000
LOWEST_CELL_M = 0.05
EDGE_TOLERANCE = 1e-7  # in grid cells, far below any stored coordinate step
ORDER = "__order"
POSITION = "__position"
# Points per decoded run in a worker process; ten workers hold about 3 GB.
WORKER_BATCH_POINTS = 2_000_000
FILTERS = laspy.DecompressionSelection


SELECTION = FILTERS.XY_RETURNS_CHANNEL | FILTERS.Z
SELECTION_EXTRA = SELECTION | FILTERS.ALL_EXTRA_BYTES


def is_copc(path: Path) -> bool:
    with laspy.open(path) as reader:
        return any(v.user_id == "copc" and v.record_id == 1 for v in reader.header.vlrs)


def batches(path: Path, extra_bytes: bool):
    selection = SELECTION_EXTRA if extra_bytes else SELECTION
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
    header.generating_software = "3dtrees-raster"
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
    def lowest(frames):
        return pl.concat(frames).sort(["Z", ORDER]).unique(["cell_x", "cell_y"], keep="first",
                                                           maintain_order=True)

    # The lowest point per cell (first by order on ties) does not depend on
    # processing order: reduce each batch, and merge into the kept table only
    # once the pending rows outgrow it, instead of re-sorting it every batch.
    kept, pending, pending_rows, seen = None, [], 0, 0
    for points in batches(args.input, bool(keep)):
        frame = table(points, seen, keep, args.order_dimension)
        seen += len(points)
        valid = accepted(frame[ORDER].to_numpy(), rejected)
        cells = frame.filter(pl.Series(valid)).with_columns(
            cell_x=pl.Series(np.floor(np.asarray(points.x)[valid] / LOWEST_CELL_M)),
            cell_y=pl.Series(np.floor(np.asarray(points.y)[valid] / LOWEST_CELL_M)))
        pending.append(lowest([cells]))
        pending_rows += pending[-1].height
        if pending_rows > max(0 if kept is None else kept.height, BATCH_POINTS):
            kept = lowest(([] if kept is None else [kept]) + pending)
            pending, pending_rows = [], 0
    if pending:
        kept = lowest(([] if kept is None else [kept]) + pending)
    if args.order_dimension:
        # As a COPC the DTM's memory admission counts each window; a plain LAZ
        # makes every window reserve the whole thinned cloud.
        with tempfile.TemporaryDirectory(dir=args.output.parent) as scratch:
            points = Path(scratch) / "points.laz"
            write(points, header, kept, keep)
            write_copc(points, args.output, Path(scratch))
    else:
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
    """Rows on the cloud's bounding box, so the reduced header keeps its extent;
    on ties the earliest record in the file."""
    rows = [frame.sort([name, POSITION], descending=[descending, False]).head(1)
            for name in ("X", "Y", "Z") for descending in (False, True)]
    return pl.concat(rows)


class VoxelReduction:
    """Per-batch voxel selection and its running reduction. Batches can be
    folded in any order and partial reductions merged: the kept rows are first
    by source order per key, ties of the bounding-box rows go to the earliest
    record in the file (`__position`)."""

    def __init__(self, args, header):
        self.args, self.scale_z = args, header.scales[2]
        dtm = gdal.Open(str(args.dtm))
        self.left, self.cell_x, _, self.top, _, self.cell_y = dtm.GetGeoTransform()
        band = dtm.GetRasterBand(1)
        self.terrain = band.ReadAsArray().astype(float)
        if band.GetNoDataValue() is not None:
            self.terrain[self.terrain == band.GetNoDataValue()] = np.nan
        edges = json.loads(args.tile_edges.read_text())
        self.x_edges, self.y_edges = np.array(edges["x"]), np.array(edges["y"])
        self.names = args.dimension + ([args.order_dimension] if args.order_dimension else [])
        self.kept, self.unsure, self.bounds, self.seen = None, [], None, 0

    def add(self, points, rejected, position):
        args, resolution, scale_z = self.args, self.args.voxel_resolution, self.scale_z
        rows, columns = self.terrain.shape
        x, y, z = np.asarray(points.x), np.asarray(points.y), np.asarray(points.z)
        frame = table(points, self.seen, self.names, args.order_dimension).with_columns(
            pl.Series(POSITION, np.arange(position, position + len(points), dtype=np.int64)))
        self.seen += len(points)
        self.bounds = self.merge_bounds([extremes(frame)])
        valid = accepted(frame[ORDER].to_numpy(), rejected)
        # terra::cellFromXY: the last column and row include the far raster edge.
        column_position = (x - self.left) / self.cell_x
        row_position = (y - self.top) / self.cell_y
        column = np.minimum(np.floor(column_position), columns - 1).astype(np.int64)
        row = np.minimum(np.floor(row_position), rows - 1).astype(np.int64)
        inside = (column >= 0) & (row >= 0) & (column_position <= columns) & (row_position <= rows)
        ground = np.full(len(x), np.nan)
        ground[inside] = self.terrain[row[inside], column[inside]]
        # lidR stores heights on the Z lattice with a zero offset.
        steps = (z - ground) / scale_z
        height = np.trunc(steps + np.copysign(0.5, steps)) * scale_z
        known = np.isfinite(height)
        edge = (near_edge(x / resolution) | near_edge(y / resolution)
                | near_edge(height / (resolution / 2)) | near_edge(column_position)
                | near_edge(row_position))
        for values, grid in ((x, self.x_edges), (y, self.y_edges)):
            nearest = np.minimum(np.searchsorted(grid, values), len(grid) - 1)
            edge |= np.abs(values - grid[nearest]) < EDGE_TOLERANCE * resolution
            edge |= np.abs(values - grid[np.maximum(nearest - 1, 0)]) < EDGE_TOLERANCE * resolution
        usable = known & (height >= -2 * scale_z) & (height <= args.maximum_height + 2 * scale_z)
        self.unsure.append(frame.filter(pl.Series(valid & (~known | (usable & edge)))))
        chosen = valid & usable & ~edge
        batch = voxel_columns(frame.filter(pl.Series(chosen)), x[chosen], y[chosen], height[chosen],
                              resolution, self.x_edges, self.y_edges)
        self.merge_kept([decisive(batch, args.dimension, scale_z)])

    def merge_bounds(self, parts):
        return extremes(pl.concat(([] if self.bounds is None else [self.bounds]) + parts))

    def merge_kept(self, parts):
        frames = ([] if self.kept is None else [self.kept]) + parts
        self.kept = decisive(pl.concat(frames), self.args.dimension, self.scale_z)

    def partial(self):
        return self.kept, self.unsure, self.bounds, self.seen


WORKER = {}


def start_worker(args, rejected_path):
    WORKER["reader"] = laspy.CopcReader.open(args.input, decompression_selection=SELECTION_EXTRA)
    WORKER["reduction"] = VoxelReduction(args, WORKER["reader"].header)
    WORKER["rejected"] = None if rejected_path is None else np.load(rejected_path, mmap_mode="r")


def reduce_chunks(task):
    """Worker: fold a run of COPC chunks (in file order) into a partial reduction."""
    chunks, position = task
    reduction = WORKER["reduction"]
    reduction.kept, reduction.unsure, reduction.bounds, reduction.seen = None, [], None, 0
    batch, size = [], 0
    for chunk in [*chunks, None]:
        if batch and (chunk is None or size + chunk.point_count > WORKER_BATCH_POINTS):
            points = WORKER["reader"]._fetch_and_decompress_points_of_nodes(batch)
            reduction.add(points, WORKER["rejected"], position)
            position += len(points)
            batch, size = [], 0
        if chunk is not None:
            batch.append(chunk)
            size += chunk.point_count
    return reduction.partial()


def voxels(args) -> dict:
    header = source_header(args.input)
    rejected = rejected_ids(args.rejected_directory)
    reduction = VoxelReduction(args, header)
    workers = max(1, int(os.environ.get("POLARS_MAX_THREADS", "1")))
    chunks = copc_chunks(str(args.input)) if args.order_dimension and is_copc(args.input) else []
    if workers > 1 and chunks:
        # Chunk runs in worker processes; their partial reductions merge into
        # the same rows as one pass over the file.
        target = max(1, sum(c.point_count for c in chunks) // (workers * 4))
        tasks, task, size, position = [], [], 0, 0
        for chunk in chunks:
            task.append(chunk)
            size += chunk.point_count
            if size >= target:
                tasks.append((task, position))
                position += size
                task, size = [], 0
        tasks += [(task, position)] if task else []
        with tempfile.TemporaryDirectory(dir=args.output.parent) as scratch:
            rejected_path = None
            if rejected is not None:
                rejected_path = Path(scratch) / "rejected.npy"
                np.save(rejected_path, rejected)
            # One thread per worker process; the parent keeps its pools for merging.
            threads = {name: os.environ.get(name) for name in ("POLARS_MAX_THREADS", "RAYON_NUM_THREADS")}
            os.environ.update(POLARS_MAX_THREADS="1", RAYON_NUM_THREADS="1")
            try:
                with multiprocessing.get_context("spawn").Pool(
                        workers, start_worker, (args, rejected_path)) as pool:
                    for kept, unsure, bounds, seen in pool.imap_unordered(reduce_chunks, tasks):
                        if kept is not None:
                            reduction.merge_kept([kept])
                        if bounds is not None:
                            reduction.bounds = reduction.merge_bounds([bounds])
                        reduction.unsure += unsure
                        reduction.seen += seen
            finally:
                for name, value in threads.items():
                    if value is None:
                        os.environ.pop(name, None)
                    else:
                        os.environ[name] = value
    else:
        position = 0
        for points in batches(args.input, True):
            reduction.add(points, rejected, position)
            position += len(points)
    kept, unsure, bounds, seen = reduction.partial()
    names = reduction.names
    columns_out = ["X", "Y", "Z", ORDER] + names
    selected = pl.concat([kept.select(columns_out), bounds.select(columns_out)]
                         + [part.select(columns_out) for part in unsure])
    selected = selected.unique(ORDER, keep="first")
    if args.order_dimension:
        with tempfile.TemporaryDirectory(dir=args.output.parent) as scratch:
            points = Path(scratch) / "points.laz"
            write(points, header, selected, names)
            write_copc(points, args.output, Path(scratch))
    else:
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
