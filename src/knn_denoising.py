#!/usr/bin/env python3
"""kNN distance noise filter over the whole point cloud (Julian Frey's proposal).

Each point's mean distance to its k nearest neighbours, counting the point itself
as lidR::sor() does (so k = 8 averages over 7 other points); a point is noise when
that distance exceeds the cloud's mean plus `sigma` standard deviations. Rejected
source ordinals are written per fixed XY window (float64 little-endian,
ascending) so later reads load only the windows they overlap.

A COPC is processed window by window: each window's points and a buffer around
it are read through the spatial index (XYZ and the extra bytes that hold the
order key, nothing else), and the buffer grows until every point's k-th
neighbour lies inside the read, so the distances equal those of one tree over
the whole cloud. The distances are kept on disk until the cloud's mean and
standard deviation are known. A plain LAS/LAZ uses one tree over all points.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import math
import tempfile

import laspy
from laspy.copc import Bounds, load_octree_for_query
import numpy as np
from scipy.spatial import cKDTree

BYTES_PER_POINT = 96  # coordinates, ordinals, KD-tree, distances and query batches
# A window read also holds the decoded point records and the buffer's points.
WINDOW_BYTES_PER_POINT = 320
QUERY_BATCH = 2_000_000


def read_points(path: Path, dimension: str):
    with laspy.open(path) as reader:
        n = reader.header.point_count
        xyz = np.empty((n, 3))
        ids = np.empty(n)
        start = 0
        for points in reader.chunk_iterator(5_000_000):
            end = start + len(points)
            xyz[start:end, 0], xyz[start:end, 1], xyz[start:end, 2] = points.x, points.y, points.z
            ids[start:end] = points[dimension]
            start = end
    if start != n:
        raise ValueError("point count differs from the LAS header")
    if not np.all(np.isfinite(xyz)):
        raise ValueError("denoising encountered nonfinite source coordinates")
    if np.any((ids < 0) | (ids != np.floor(ids)) | (ids > 2**53)):
        raise ValueError("denoising requires exact finite nonnegative source ordinals")
    return xyz, ids


def knn_distances(xyz: np.ndarray, neighbors: int, threads: int) -> np.ndarray:
    local = xyz - xyz.min(axis=0)
    tree = cKDTree(local, leafsize=32, balanced_tree=False, compact_nodes=False)
    distances = np.empty(len(local))
    for start in range(0, len(local), QUERY_BATCH):
        d, _ = tree.query(local[start:start + QUERY_BATCH], k=neighbors, workers=threads)
        distances[start:start + QUERY_BATCH] = d[:, 1:].sum(axis=1) / (neighbors - 1)
    return distances


XYZ_AND_EXTRA_BYTES = (laspy.DecompressionSelection.xy_returns_channel()
                       | laspy.DecompressionSelection.Z | laspy.DecompressionSelection.ALL_EXTRA_BYTES)


def is_copc(path: Path) -> bool:
    with laspy.open(path) as reader:
        return any(v.user_id == "copc" and v.record_id == 1 for v in reader.header.vlrs)


def window_distances(reader, window, size, neighbors, threads, budget, header, dimension):
    """Distances and ordinals of one window's points, read with an exact buffer."""
    x0, y0 = window[0] * size, window[1] * size
    buffer = size / 4
    while True:
        lower = np.array([x0 - buffer, y0 - buffer, header.mins[2] - 1])
        upper = np.array([x0 + size + buffer, y0 + size + buffer, header.maxs[2] + 1])
        nodes = load_octree_for_query(reader.source, reader.copc_info, reader.root_page,
                                      query_bounds=Bounds(lower, upper))
        if sum(node.point_count for node in nodes) * WINDOW_BYTES_PER_POINT > budget:
            raise SystemExit(f"RASTER_FLAG memory_budget_exceeded: the kNN window at {x0:.0f} {y0:.0f} "
                             f"reads too many points for {budget / 1024**3:.1f} GiB")
        points = reader._fetch_and_decompress_points_of_nodes(nodes)
        x, y = np.asarray(points.x), np.asarray(points.y)
        inside = (x >= lower[0]) & (x <= upper[0]) & (y >= lower[1]) & (y <= upper[1])
        xyz = np.column_stack([x[inside], y[inside], np.asarray(points.z)[inside]])
        all_ids = np.asarray(points.array[dimension])[inside]
        del points, x, y
        core = ((np.floor(xyz[:, 0] / size) == window[0]) & (np.floor(xyz[:, 1] / size) == window[1]))
        if not core.any():
            return np.empty(0), np.empty(0)
        local = xyz - xyz.min(axis=0)
        tree = cKDTree(local, leafsize=32, balanced_tree=False, compact_nodes=False)
        query = local[core]
        # Batched like the whole-cloud tree: distances and indices are k wide,
        # so keep only each point's mean distance and its k-th neighbour.
        mean_distance, kth = np.empty(len(query)), np.empty(len(query))
        for start in range(0, len(query), QUERY_BATCH):
            d, _ = tree.query(query[start:start + QUERY_BATCH], k=neighbors, workers=threads)
            mean_distance[start:start + QUERY_BATCH] = d[:, 1:].sum(axis=1) / (neighbors - 1)
            kth[start:start + QUERY_BATCH] = d[:, -1]
        # Every k-th neighbour must be closer than the edge of the read, or a
        # point outside it could be nearer: then read a wider buffer.
        edge = np.minimum.reduce([xyz[core, 0] - lower[0], upper[0] - xyz[core, 0],
                                  xyz[core, 1] - lower[1], upper[1] - xyz[core, 1]])
        whole = (lower[0] <= header.mins[0] and lower[1] <= header.mins[1] and
                 upper[0] >= header.maxs[0] and upper[1] >= header.maxs[1])
        if whole or np.all(kth < edge):
            return mean_distance, all_ids[core].astype(np.float64)
        buffer *= 2


def windowed(args, n, budget) -> dict:
    """Window-by-window distances, then the global threshold over all of them."""
    with laspy.CopcReader.open(args.input, decompression_selection=XYZ_AND_EXTRA_BYTES) as reader:
        header = reader.header
        size = args.window_size
        xs = range(math.floor(header.mins[0] / size), math.floor(header.maxs[0] / size) + 1)
        ys = range(math.floor(header.mins[1] / size), math.floor(header.maxs[1] / size) + 1)
        with tempfile.TemporaryDirectory(dir=args.output_dir, prefix=".distances_") as scratch:
            windows, count, total = [], 0, 0.0
            for ix in xs:
                for iy in ys:
                    d, ids = window_distances(reader, (ix, iy), size, args.neighbors,
                                              args.threads, budget, header, args.order_dimension)
                    if not len(d):
                        continue
                    if np.any((ids < 0) | (ids != np.floor(ids)) | (ids > 2**53)):
                        raise ValueError("denoising requires exact finite nonnegative source ordinals")
                    path = Path(scratch) / f"{ix}_{iy}.npy"
                    np.save(path, np.column_stack([d, ids]))
                    windows.append(((ix, iy), path, len(d)))
                    count += len(d)
                    total += math.fsum(d)
            if count != n:
                raise ValueError(f"windows hold {count} of {n} points")
            mean = total / count
            squares = math.fsum(math.fsum((np.load(path)[:, 0] - mean) ** 2) for _, path, _ in windows)
            sd = math.sqrt(squares / (count - 1))
            threshold = mean + args.sigma * sd
            partitions, rejected_total = [], 0
            for (ix, iy), path, core in sorted(windows):
                values = np.load(path)
                rejected = np.sort(values[values[:, 0] > threshold, 1])
                filename = f"core_{ix}_{iy}.bin"
                if len(rejected):
                    rejected.astype("<f8").tofile(args.output_dir / filename)
                rejected_total += len(rejected)
                partitions.append(dict(ix=ix, iy=iy, filename=filename, core_points=int(core),
                                       rejected_points=int(len(rejected)),
                                       accepted_points=int(core - len(rejected)),
                                       insufficient_context_points=0))
    return dict(point_count=int(n), rejected_points=int(rejected_total), distance_mean_m=mean,
                distance_sd_m=sd, threshold_m=threshold, insufficient_context_points=0,
                partitions=partitions)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--order-dimension", required=True)
    parser.add_argument("--neighbors", type=int, required=True)
    parser.add_argument("--sigma", type=float, required=True)
    parser.add_argument("--window-size", type=float, required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--memory-budget-bytes", type=float, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.neighbors < 2:
        raise ValueError("--denoise-neighbors must be at least 2 (it counts the point itself)")
    with laspy.open(args.input) as reader:
        n = reader.header.point_count
    # One tree over the whole cloud when it fits; window by window otherwise.
    whole_cloud_bytes = n * BYTES_PER_POINT + 512 * 1024**2
    if n > args.neighbors and whole_cloud_bytes > args.memory_budget_bytes and is_copc(args.input):
        manifest = windowed(args, n, args.memory_budget_bytes - 512 * 1024**2)
        (args.output_dir / "manifest.json").write_text(json.dumps(manifest, allow_nan=False))
        return
    required = n * BYTES_PER_POINT + 512 * 1024**2
    if required > args.memory_budget_bytes:
        raise SystemExit(f"RASTER_FLAG memory_budget_exceeded: the kNN filter needs "
                         f"{required / 1024**3:.1f} GiB for {n} points, "
                         f"{args.memory_budget_bytes / 1024**3:.1f} GiB available")
    xyz, ids = read_points(args.input, args.order_dimension)
    if n > args.neighbors:
        distances = knn_distances(xyz, args.neighbors, args.threads)
        mean, sd = float(distances.mean()), float(distances.std(ddof=1))
        threshold = mean + args.sigma * sd
        rejected = distances > threshold
        insufficient = 0
    else:
        mean = sd = threshold = None
        rejected = np.zeros(n, bool)
        insufficient = n
    cells = np.floor(xyz[:, :2] / args.window_size).astype(np.int64)
    keys, inverse, core = np.unique(cells, axis=0, return_inverse=True, return_counts=True)
    inverse = inverse.reshape(-1)
    counts = np.bincount(inverse[rejected], minlength=len(keys))
    order = np.argsort(inverse[rejected], kind="stable")
    rejected_ids, rejected_cells = ids[rejected][order], inverse[rejected][order]
    edges = np.searchsorted(rejected_cells, np.arange(len(keys) + 1))
    partitions = []
    for i, (ix, iy) in enumerate(keys.tolist()):
        filename = f"core_{ix}_{iy}.bin"
        if counts[i]:
            np.sort(rejected_ids[edges[i]:edges[i + 1]]).astype("<f8").tofile(args.output_dir / filename)
        partitions.append(dict(ix=ix, iy=iy, filename=filename, core_points=int(core[i]),
                               rejected_points=int(counts[i]), accepted_points=int(core[i] - counts[i]),
                               insufficient_context_points=int(core[i]) if insufficient else 0))
    manifest = dict(point_count=int(n), rejected_points=int(rejected.sum()),
                    distance_mean_m=mean, distance_sd_m=sd, threshold_m=threshold,
                    insufficient_context_points=insufficient, partitions=partitions)
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, allow_nan=False))


if __name__ == "__main__":
    main()
