#!/usr/bin/env python3
"""kNN distance noise filter over the whole point cloud (Julian Frey's proposal).

Each point's mean distance to its k nearest neighbours, counting the point itself
as lidR::sor() does (so k = 8 averages over 7 other points); a point is noise when
that distance exceeds the cloud's mean plus `sigma` standard deviations. One
multithreaded KD-tree over all points gives the same decisions as lidR, without
windows. Rejected source ordinals are written per fixed XY window (float64
little-endian, ascending) so later reads load only the windows they overlap.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import laspy
import numpy as np
from scipy.spatial import cKDTree

BYTES_PER_POINT = 96  # coordinates, ordinals, KD-tree, distances and query batches
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
    required = n * BYTES_PER_POINT + 512 * 1024**2
    if required > args.memory_budget_bytes:
        raise SystemExit(f"FORESTSTRUCTURE_FLAG memory_budget_exceeded: the kNN filter needs "
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
