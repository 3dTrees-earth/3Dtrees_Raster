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
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
import json
from pathlib import Path

import math
import multiprocessing
import os
import tempfile

import laspy
from laspy.copc import Bounds, load_octree_for_query
import numpy as np
from scipy.spatial import cKDTree

BYTES_PER_POINT = 96  # coordinates, ordinals, KD-tree, distances and query batches
# A window read holds the compressed bytes and the decoded records of every node
# it touches plus float X, Y, Z, their copies inside the box and masks; then the
# tree, distances and copies of the points inside. Measured: 1887 (2.0 B TLS
# points) peaked at 59 GiB of anonymous memory with 72 and 320 for a 53 GiB budget.
DECODED_BYTES_PER_POINT = 100
WINDOW_BYTES_PER_POINT = 370
QUERY_BATCH = 2_000_000
# Nodes read at once for far neighbours: small, so the search radius shrinks early.
BOUNDED_BATCH_POINTS = 1_000_000
# A 60 m window splits down to 60 / 2^6 < 1 m before the budget is declared too small.
MAX_SPLITS = 6
WORKER_BYTES = 256 * 1024**2
# Planning stops at this many parts per worker, and splits a part only when its
# most expensive quarter costs at most this fraction of it.
MAX_PARTS_PER_WORKER = 16
SPLIT_GAIN = 0.7


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
    x0, x1, y0, y1 = part_bounds(window, size, part)
    buffer = max(x1 - x0, y1 - y0) / 4
    lower = np.array([x0 - buffer, y0 - buffer, header.mins[2] - 1])
    upper = np.array([x1 + buffer, y1 + buffer, header.maxs[2] + 1])
    nodes = load_octree_for_query(reader.source, reader.copc_info, reader.root_page,
                                  query_bounds=Bounds(lower, upper))
    # Decoding costs decoded_bytes(); the points inside the read then cost
    # WINDOW_BYTES_PER_POINT each.
    if decoded_bytes(nodes, header) > budget:
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
    ids = all_ids[core].astype(np.float64)
    del xyz, all_ids, core
    # A point is exact when its k-th neighbour is closer than every side of the
    # read that has cloud beyond it; otherwise a point outside could be nearer.
    low, high = lower - header.mins, upper - header.mins
    sides = [(query[:, 0] - low[0], lower[0] > header.mins[0]),
             (high[0] - query[:, 0], upper[0] < header.maxs[0]),
             (query[:, 1] - low[1], lower[1] > header.mins[1]),
             (high[1] - query[:, 1], upper[1] < header.maxs[1])]
    edge = np.minimum.reduce([distance if open_ else np.full(len(distance), np.inf)
                              for distance, open_ in sides])
    del sides
    open_points = np.flatnonzero(kth >= edge)
    if len(open_points):
        # The k-th distance found so far bounds the true one, so the box of
        # each open point plus that radius holds all its true neighbours.
        # (It is infinite where the read held fewer than k points.) Only the
        # window's distances and ordinals stay in memory meanwhile.
        open_query, open_kth = query[open_points], kth[open_points]
        del query, kth, edge
        mean_distance[open_points] = far_distances(
            reader, open_query, open_kth, buffer, neighbors, threads,
            budget - mean_distance.nbytes - ids.nbytes, header)
    return mean_distance, ids


def decoded_bytes(nodes, header):
    """Peak bytes of read_box() over these nodes."""
    return sum(node.byte_size + node.point_count * (header.point_format.size + DECODED_BYTES_PER_POINT)
               for node in nodes)


def read_box(reader, nodes, lower, upper, dimension):
    """XYZ and ordinals of the nodes' points inside the XY box."""
    points = reader._fetch_and_decompress_points_of_nodes(nodes)
    x = np.asarray(points.x)
    inside = (x >= lower[0]) & (x <= upper[0])
    del x
    y = np.asarray(points.y)
    inside &= (y >= lower[1]) & (y <= upper[1])
    del y
    xyz = np.empty((np.count_nonzero(inside), 3))
    for axis, name in enumerate("xyz"):
        xyz[:, axis] = np.asarray(getattr(points, name))[inside]
    ids = np.asarray(points.array[dimension])[inside] if dimension else None
    return xyz, ids


FAR_POINTS = [0]


def far_distances(reader, query, radius, buffer, neighbors, threads, budget, header):
    """Exact mean distances of points whose k-th neighbour may lie beyond the
    read. A point with fewer than k neighbours in the read has no bound yet: its
    search radius starts at twice the buffer and doubles until k neighbours lie
    within it, instead of searching the whole cloud at once."""
    FAR_POINTS[0] += len(query)
    mean = np.empty(len(query))
    reach = np.where(np.isfinite(radius), radius, 2 * buffer)
    todo = np.arange(len(query))
    while len(todo):
        found, kth = bounded_distances(reader, query[todo], reach[todo], neighbors, threads, budget, header)
        exact = kth <= reach[todo]
        mean[todo[exact]] = found[exact]
        todo = todo[~exact]
        reach[todo] *= 2
    return mean


def bounded_distances(reader, query, radius, neighbors, threads, budget, header):
    """Mean and k-th distances of the query points (local frame) over the points
    within `radius`, exact where the k-th distance is within it. The COPC nodes in the box are read nearest
    first, in batches, keeping each point's k smallest distances; a node farther
    from every point than its k-th distance so far cannot hold a neighbour and is
    skipped, so an isolated point does not read the cloud around it."""
    # A millimetre more than the radius absorbs rounding between the frames.
    reach = radius + 1e-3
    lower = np.append((query[:, :2] - reach[:, None]).min(axis=0) + header.mins[:2], header.mins[2] - 1)
    upper = np.append((query[:, :2] + reach[:, None]).max(axis=0) + header.mins[:2], header.maxs[2] + 1)
    nodes = load_octree_for_query(reader.source, reader.copc_info, reader.root_page,
                                  query_bounds=Bounds(lower, upper))
    points = query + header.mins

    def gaps(node):
        """Distance from each query point to the node's cube."""
        outside = np.maximum(node.bounds.mins - points, 0) + np.maximum(points - node.bounds.maxs, 0)
        return np.sqrt((outside ** 2).sum(axis=1))

    nearest = sorted(((gaps(node), node) for node in nodes if node.point_count > 0),
                     key=lambda pair: pair[0].min())
    best = np.full((len(query), neighbors), np.inf)
    limit = reach.copy()
    per_point = header.point_format.size + DECODED_BYTES_PER_POINT + WINDOW_BYTES_PER_POINT
    batch, cost, size = [], 0, 0

    def read(batch):
        xyz, _ = read_box(reader, batch, lower, upper, None)
        if not len(xyz):
            return best
        tree = cKDTree(xyz - header.mins, leafsize=32, balanced_tree=False, compact_nodes=False)
        d, _ = tree.query(query, k=min(neighbors, len(xyz)), workers=threads)
        return np.sort(np.concatenate([best, d.reshape(len(query), -1)], axis=1), axis=1)[:, :neighbors]

    for gap, node in nearest:
        if gap.min() > limit.max():
            break
        if not np.any(gap <= limit):
            continue
        node_cost = node.byte_size + node.point_count * per_point
        if batch and (cost + node_cost > budget or size + node.point_count > BOUNDED_BATCH_POINTS):
            best = read(batch)
            limit = np.minimum(reach, best[:, -1] + 1e-3)
            batch, cost, size = [], 0, 0
            if not np.any(gap <= limit):
                continue
        batch.append(node)
        cost += node_cost
        size += node.point_count
    if batch:
        best = read(batch)
    return best[:, 1:].sum(axis=1) / (neighbors - 1), best[:, -1]


def quarters(reader, window, size, neighbors, threads, budget, header, dimension, part, depth, bounds):
    """The distances of a part's four quarters, in a fixed order."""
    if depth == MAX_SPLITS:
        raise SystemExit(f"RASTER_FLAG memory_budget_exceeded: the kNN window at {bounds[0]:.0f} "
                         f"{bounds[2]:.0f} reads too many points for {budget / 1024**3:.1f} GiB")
    parts = [window_distances(reader, window, size, neighbors, threads, budget, header,
                              dimension, quarter, depth + 1)
             for quarter in part_quarters(part, bounds)]
    return np.concatenate([d for d, _ in parts]), np.concatenate([ids for _, ids in parts])


WORKER = {}
FULL = (-math.inf, math.inf, -math.inf, math.inf)


def part_bounds(window, size, part):
    """The XY rectangle [x0, x1) x [y0, y1) of a window part."""
    return (max(window[0] * size, part[0]), min((window[0] + 1) * size, part[1]),
            max(window[1] * size, part[2]), min((window[1] + 1) * size, part[3]))


def part_quarters(part, bounds):
    xm, ym = (bounds[0] + bounds[1]) / 2, (bounds[2] + bounds[3]) / 2
    return [(part[0], xm, part[2], ym), (xm, part[1], part[2], ym),
            (part[0], xm, ym, part[3]), (xm, part[1], ym, part[3])]


def part_cost(reader, header, window, size, part):
    """Estimated bytes of reading one part: its nodes decoded, plus the points
    inside the read in proportion to each node's overlap with it."""
    x0, x1, y0, y1 = part_bounds(window, size, part)
    buffer = max(x1 - x0, y1 - y0) / 4
    lower, upper = np.array([x0 - buffer, y0 - buffer]), np.array([x1 + buffer, y1 + buffer])
    nodes = load_octree_for_query(reader.source, reader.copc_info, reader.root_page, query_bounds=Bounds(
        np.append(lower, header.mins[2] - 1), np.append(upper, header.maxs[2] + 1)))
    inside = 0.0
    for node in nodes:
        low, high = node.bounds.mins[:2], node.bounds.maxs[:2]
        overlap = np.clip(np.minimum(high, upper) - np.maximum(low, lower), 0, None)
        inside += node.point_count * overlap.prod() / max((high - low).prod(), 1e-12)
    return decoded_bytes(nodes, header) + inside * WINDOW_BYTES_PER_POINT


def plan_parts(reader, header, size, workers, budget):
    """Windows, split into quarters while a part is estimated above its share of
    the budget or while there are too few parts to keep every worker busy. A
    split that does not make the quarters clearly cheaper (the read is
    dominated by large shallow COPC nodes, as on very dense TLS) is not taken."""
    xs = range(math.floor(header.mins[0] / size), math.floor(header.maxs[0] / size) + 1)
    ys = range(math.floor(header.mins[1] / size), math.floor(header.maxs[1] / size) + 1)
    parts = [(part_cost(reader, header, (ix, iy), size, FULL), (ix, iy), FULL, 0) for ix in xs for iy in ys]
    final = []
    while parts and len(parts) + len(final) < MAX_PARTS_PER_WORKER * workers:
        parts.sort(key=lambda item: -item[0])
        cost, window, part, depth = parts[0]
        if cost <= budget / workers and len(parts) + len(final) >= 3 * workers:
            break
        if depth == MAX_SPLITS - 1:
            final.append(parts.pop(0))
            continue
        quarters_ = [(part_cost(reader, header, window, size, q), window, q, depth + 1)
                     for q in part_quarters(part, part_bounds(window, size, part))]
        if max(q[0] for q in quarters_) > SPLIT_GAIN * cost:
            final.append(parts.pop(0))
        else:
            parts[0:1] = quarters_
    return sorted(parts + final, key=lambda item: -item[0])


def start_worker(path):
    WORKER["reader"] = laspy.CopcReader.open(path, decompression_selection=XYZ_AND_EXTRA_BYTES)


def part_distances(task):
    """Worker: one part's distances and ordinals, saved to the scratch directory."""
    window, part, depth, size, neighbors, budget, dimension, path, threads = task
    reader = WORKER["reader"]
    far = FAR_POINTS[0]
    try:
        d, ids = window_distances(reader, window, size, neighbors, threads, budget, reader.header,
                                  dimension, part, depth)
    except SystemExit as error:
        # A pool worker that exits loses its task; report the flag instead.
        raise RuntimeError(str(error)) from None
    far = FAR_POINTS[0] - far
    if not len(d):
        return window, None, 0, far
    if np.any((ids < 0) | (ids != np.floor(ids)) | (ids > 2**53)):
        raise ValueError("denoising requires exact finite nonnegative source ordinals")
    np.save(path, np.column_stack([d, ids]))
    return window, path, len(d), far


def window_values(paths):
    return np.concatenate([np.load(path) for path in paths])


def windowed(args, n, budget) -> dict:
    """Window-by-window distances in parallel worker processes, then the global
    threshold over all of them. Windows are split into parts so that every
    worker has work and a part fits its share of the budget; per-window
    sums are exact (fsum) and added in window order, so the result does not
    depend on the number of workers."""
    workers = max(1, args.threads)
    # Each worker process holds Python, laspy and the COPC hierarchy besides its part.
    budget -= workers * WORKER_BYTES
    with laspy.CopcReader.open(args.input, decompression_selection=XYZ_AND_EXTRA_BYTES) as reader:
        header = reader.header
        parts = plan_parts(reader, header, args.window_size, workers, budget)
    reservations = sorted(min(budget, cost * 1.25) for cost, *_ in parts)
    if reservations[len(reservations) // 2] > budget / 2:
        # Parts this large run one at a time anyway (very dense TLS, whose reads
        # are dominated by large shallow COPC nodes): read whole windows, so each
        # decodes those nodes once, in one worker that queries with every thread.
        parts = [(budget / 1.25, window, FULL, 0) for window in sorted({part[1] for part in parts})]
        workers = 1
    with tempfile.TemporaryDirectory(dir=args.output_dir, prefix=".distances_") as scratch:
        # Each part reserves its estimate (a quarter more for the error, at most
        # the budget) and is read within that reservation; parts start in order
        # of size while the reservations of the running parts fit the budget.
        pending = [(min(budget, cost * 1.25), (window, part, depth, args.window_size, args.neighbors,
                    min(budget, cost * 1.25), args.order_dimension, str(Path(scratch) / f"{i}.npy")))
                   for i, (cost, window, part, depth) in enumerate(parts)]
        files, count, running, far_points = {}, 0, {}, 0
        rayon = os.environ.get("RAYON_NUM_THREADS")
        os.environ["RAYON_NUM_THREADS"] = "1"
        try:
            with ProcessPoolExecutor(workers, mp_context=multiprocessing.get_context("spawn"),
                                     initializer=start_worker, initargs=(str(args.input),),
                                     max_tasks_per_child=1) as pool:
                while pending or running:
                    while pending and len(running) < workers and \
                            sum(running.values()) + pending[0][0] <= budget:
                        reserved, task = pending.pop(0)
                        # A part that leaves room for few others queries with the idle threads.
                        fit = max(1, min(workers, int(budget // max(reserved, 1))))
                        task = task + (max(1, args.threads // fit),)
                        running[pool.submit(part_distances, task)] = reserved
                    done, _ = wait(running, return_when=FIRST_COMPLETED)
                    for future in done:
                        del running[future]
                        try:
                            window, path, points, far = future.result()
                        except RuntimeError as error:
                            raise SystemExit(str(error)) from None
                        far_points += far
                        if path is not None:
                            files.setdefault(window, []).append(path)
                            count += points
        finally:
            if rayon is None:
                os.environ.pop("RAYON_NUM_THREADS")
            else:
                os.environ["RAYON_NUM_THREADS"] = rayon
        if count != n:
            raise ValueError(f"windows hold {count} of {n} points")
        print(f"KNN_WINDOWS parts={len(parts)} deepest={max(part[3] for part in parts)} "
              f"workers={workers} far_points={far_points}", flush=True)
        windows = sorted(files)
        total = 0.0
        for window in windows:
            total += math.fsum(window_values(files[window])[:, 0])
        mean = total / count
        squares = math.fsum(math.fsum((window_values(files[window])[:, 0] - mean) ** 2) for window in windows)
        sd = math.sqrt(squares / (count - 1))
        threshold = mean + args.sigma * sd
        partitions, rejected_total = [], 0
        for ix, iy in windows:
            values = window_values(files[(ix, iy)])
            rejected = np.sort(values[values[:, 0] > threshold, 1])
            filename = f"core_{ix}_{iy}.bin"
            if len(rejected):
                rejected.astype("<f8").tofile(args.output_dir / filename)
            rejected_total += len(rejected)
            partitions.append(dict(ix=ix, iy=iy, filename=filename, core_points=int(len(values)),
                                   rejected_points=int(len(rejected)),
                                   accepted_points=int(len(values) - len(rejected)),
                                   insufficient_context_points=0))
    # The parts read overlapping COPC nodes repeatedly, which leaves the input
    # in the active page cache; task admission counts that as used memory.
    with open(args.input, "rb") as source:
        os.posix_fadvise(source.fileno(), 0, 0, os.POSIX_FADV_DONTNEED)
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
    # A COPC runs window by window in parallel workers; one tree over the whole
    # cloud is single-threaded apart from the queries. A plain LAS/LAZ, or one
    # thread, uses the whole-cloud tree when it fits. RASTER_KNN_WINDOWED=1 or 0
    # forces either path, so tests can compare them.
    whole_cloud_bytes = n * BYTES_PER_POINT + 512 * 1024**2
    forced = os.environ.get("RASTER_KNN_WINDOWED")
    windowed_wanted = (whole_cloud_bytes > args.memory_budget_bytes or forced == "1"
                       or (forced != "0" and args.threads > 1))
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
