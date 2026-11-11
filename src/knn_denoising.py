#!/usr/bin/env python3
"""kNN distance noise filter over the whole point cloud (Julian Frey's proposal).

Each point's mean distance to its k nearest neighbours, counting the point itself
as lidR::sor() does (so k = 8 averages over 7 other points); a point is noise when
that distance exceeds the cloud's mean plus `sigma` standard deviations. Rejected
source ordinals are written per fixed XY window (float64 little-endian,
ascending) so later reads load only the windows they overlap.

A COPC is processed window by window: each window's points and a buffer around
it are read through the spatial index (XYZ and the extra bytes that hold the
order key, nothing else). A point whose k-th neighbour lies farther than the
edge of the read gets its neighbours from the box of that radius, streamed node
batch by node batch, so the distances equal those of one tree over the whole
cloud. A window too dense for the memory budget is read in quarters
(recursively). The distances are kept on disk until the cloud's mean and
standard deviation are known. A plain LAS/LAZ uses one tree over all points.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import math
import os
import tempfile

import laspy
from laspy.copc import Bounds, load_octree_for_query
import numpy as np
from scipy.spatial import cKDTree

BYTES_PER_POINT = 96  # coordinates, ordinals, KD-tree, distances and query batches
# A window read holds the decoded records of every node it touches (records plus
# float X, Y, Z and masks), and the tree, distances and copies of the points inside.
DECODED_BYTES_PER_POINT = 40
WINDOW_BYTES_PER_POINT = 320
QUERY_BATCH = 2_000_000
# A 60 m window splits down to 60 / 2^6 < 1 m before the budget is declared too small.
MAX_SPLITS = 6


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


def window_distances(reader, window, size, neighbors, threads, budget, header, dimension,
                     part=(-math.inf, math.inf, -math.inf, math.inf), depth=0):
    """Distances and ordinals of one window's points, read with a buffer.

    `part` narrows the window to [x0, x1) x [y0, y1); a read that does not fit
    the budget is split into quarters, down to MAX_SPLITS levels, so dense
    TLS windows stay within memory with the same distances.
    """
    x0, x1 = max(window[0] * size, part[0]), min((window[0] + 1) * size, part[1])
    y0, y1 = max(window[1] * size, part[2]), min((window[1] + 1) * size, part[3])
    buffer = max(x1 - x0, y1 - y0) / 4
    lower = np.array([x0 - buffer, y0 - buffer, header.mins[2] - 1])
    upper = np.array([x1 + buffer, y1 + buffer, header.maxs[2] + 1])
    nodes = load_octree_for_query(reader.source, reader.copc_info, reader.root_page,
                                  query_bounds=Bounds(lower, upper))
    # Decoded nodes cost their records plus float coordinates; the points
    # inside the read then cost WINDOW_BYTES_PER_POINT each.
    if sum(node.point_count for node in nodes) * (header.point_format.size + DECODED_BYTES_PER_POINT) > budget:
        return quarters(reader, window, size, neighbors, threads, budget, header, dimension,
                        part, depth, (x0, x1, y0, y1))
    xyz, all_ids = read_box(reader, nodes, lower, upper, dimension)
    if len(xyz) * WINDOW_BYTES_PER_POINT > budget:
        del xyz, all_ids
        return quarters(reader, window, size, neighbors, threads, budget, header, dimension,
                        part, depth, (x0, x1, y0, y1))
    core = ((np.floor(xyz[:, 0] / size) == window[0]) & (np.floor(xyz[:, 1] / size) == window[1])
            & (xyz[:, 0] >= part[0]) & (xyz[:, 0] < part[1])
            & (xyz[:, 1] >= part[2]) & (xyz[:, 1] < part[3]))
    if not core.any():
        return np.empty(0), np.empty(0)
    # One frame for every read, the cloud's minimum, as the whole-cloud tree uses.
    local = xyz - header.mins
    tree = cKDTree(local, leafsize=32, balanced_tree=False, compact_nodes=False)
    query = local[core]
    del local
    # Batched like the whole-cloud tree: distances and indices are k wide,
    # so keep only each point's mean distance and its k-th neighbour.
    mean_distance, kth = np.empty(len(query)), np.empty(len(query))
    for start in range(0, len(query), QUERY_BATCH):
        d, _ = tree.query(query[start:start + QUERY_BATCH], k=neighbors, workers=threads)
        mean_distance[start:start + QUERY_BATCH] = d[:, 1:].sum(axis=1) / (neighbors - 1)
        kth[start:start + QUERY_BATCH] = d[:, -1]
    del tree
    # A point is exact when its k-th neighbour is closer than every side of the
    # read that has cloud beyond it; otherwise a point outside could be nearer.
    sides = [(xyz[core, 0] - lower[0], lower[0] > header.mins[0]),
             (upper[0] - xyz[core, 0], upper[0] < header.maxs[0]),
             (xyz[core, 1] - lower[1], lower[1] > header.mins[1]),
             (upper[1] - xyz[core, 1], upper[1] < header.maxs[1])]
    edge = np.minimum.reduce([distance if open_ else np.full(len(distance), np.inf)
                              for distance, open_ in sides])
    open_points = np.flatnonzero(kth >= edge)
    if len(open_points):
        # The k-th distance found so far bounds the true one, so the box of
        # each open point plus that radius holds all its true neighbours.
        mean_distance[open_points] = bounded_distances(
            reader, query[open_points], kth[open_points], neighbors, threads, budget, header)
    return mean_distance, all_ids[core].astype(np.float64)


def read_box(reader, nodes, lower, upper, dimension):
    """XYZ and ordinals of the nodes' points inside the XY box."""
    points = reader._fetch_and_decompress_points_of_nodes(nodes)
    x, y = np.asarray(points.x), np.asarray(points.y)
    inside = (x >= lower[0]) & (x <= upper[0]) & (y >= lower[1]) & (y <= upper[1])
    xyz = np.column_stack([x[inside], y[inside], np.asarray(points.z)[inside]])
    ids = np.asarray(points.array[dimension])[inside] if dimension else None
    return xyz, ids


def bounded_distances(reader, query, radius, neighbors, threads, budget, header):
    """Exact mean distances of the query points (local frame) whose k nearest
    neighbours lie within `radius`: the box is streamed node batch by node
    batch, keeping each point's k smallest distances."""
    # A millimetre more than the radius absorbs rounding between the frames.
    reach = radius[:, None] + 1e-3
    lower = np.append((query[:, :2] - reach).min(axis=0) + header.mins[:2], header.mins[2] - 1)
    upper = np.append((query[:, :2] + reach).max(axis=0) + header.mins[:2], header.maxs[2] + 1)
    nodes = load_octree_for_query(reader.source, reader.copc_info, reader.root_page,
                                  query_bounds=Bounds(lower, upper))
    per_point = header.point_format.size + DECODED_BYTES_PER_POINT + WINDOW_BYTES_PER_POINT
    best = np.full((len(query), neighbors), np.inf)
    for batch in node_batches(nodes, budget // per_point):
        xyz, _ = read_box(reader, batch, lower, upper, None)
        if not len(xyz):
            continue
        tree = cKDTree(xyz - header.mins, leafsize=32, balanced_tree=False, compact_nodes=False)
        d, _ = tree.query(query, k=min(neighbors, len(xyz)), workers=threads)
        best = np.sort(np.concatenate([best, d.reshape(len(query), -1)], axis=1), axis=1)[:, :neighbors]
        del tree, xyz
    return best[:, 1:].sum(axis=1) / (neighbors - 1)


def node_batches(nodes, limit):
    """Consecutive groups of nodes holding at most `limit` points (one node at least)."""
    batch, points = [], 0
    for node in nodes:
        if batch and points + node.point_count > limit:
            yield batch
            batch, points = [], 0
        batch.append(node)
        points += node.point_count
    if batch:
        yield batch


def quarters(reader, window, size, neighbors, threads, budget, header, dimension, part, depth, bounds):
    """The distances of a part's four quarters, in a fixed order."""
    if depth == MAX_SPLITS:
        raise SystemExit(f"RASTER_FLAG memory_budget_exceeded: the kNN window at {bounds[0]:.0f} "
                         f"{bounds[2]:.0f} reads too many points for {budget / 1024**3:.1f} GiB")
    xm, ym = (bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2
    parts = [window_distances(reader, window, size, neighbors, threads, budget, header,
                              dimension, quarter, depth + 1)
             for quarter in [(part[0], xm, part[2], ym), (xm, part[1], part[2], ym),
                             (part[0], xm, ym, part[3]), (xm, part[1], ym, part[3])]]
    return np.concatenate([d for d, _ in parts]), np.concatenate([ids for _, ids in parts])


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
    # RASTER_KNN_WINDOWED=1 forces windows, so tests can compare both paths.
    windowed_wanted = (whole_cloud_bytes > args.memory_budget_bytes
                       or os.environ.get("RASTER_KNN_WINDOWED") == "1")
    if n > args.neighbors and windowed_wanted and is_copc(args.input):
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
