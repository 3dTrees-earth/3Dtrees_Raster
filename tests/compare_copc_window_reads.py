#!/usr/bin/env python3
"""Window reads through a COPC's hierarchy hold exactly the points of one full
read: every point of a sub-area once, with the same XYZ and ExtraBytes.

The tool reads COPC windows only through the hierarchy (load_octree_for_query
and the node decoder), never through laspy's record seek(), which returns wrong
records on LAZ with variable-size chunks such as every COPC (laspy 2.7 /
lazrs 0.8). This compares those window reads with one sequential read.

Usage: compare_copc_window_reads.py COPC [WINDOW_M] [AREA_M]
Checks the AREA_M square at the centre of the cloud's extent in WINDOW_M windows
(a corner can be empty on an irregular footprint).
"""
import sys

import laspy
from laspy.copc import Bounds, load_octree_for_query
import numpy as np

path = sys.argv[1]
size = float(sys.argv[2]) if len(sys.argv) > 2 else 20.0
area = float(sys.argv[3]) if len(sys.argv) > 3 else 60.0

with laspy.open(path) as reader:
    header = reader.header
    names = [d.name for d in header.point_format.extra_dimensions]
    centre = (np.asarray(header.mins[:2]) + np.asarray(header.maxs[:2])) / 2
    x0, y0 = np.maximum(header.mins[:2], centre - area / 2)
    x1, y1 = np.minimum(header.maxs[:2], centre + area / 2)
    full = []
    for points in reader.chunk_iterator(5_000_000):
        inside = (points.x >= x0) & (points.x < x1) & (points.y >= y0) & (points.y < y1)
        full.append(np.column_stack([np.asarray(points.x)[inside], np.asarray(points.y)[inside],
                                     np.asarray(points.z)[inside]]
                                    + [np.asarray(points[n])[inside].astype(float) for n in names]))
full = np.concatenate(full)

windows = []
with laspy.CopcReader.open(path) as reader:
    for wx in np.arange(x0, x1, size):
        for wy in np.arange(y0, y1, size):
            lower = np.array([wx, wy, header.mins[2] - 1])
            upper = np.array([min(wx + size, x1), min(wy + size, y1), header.maxs[2] + 1])
            nodes = load_octree_for_query(reader.source, reader.copc_info, reader.root_page,
                                          query_bounds=Bounds(lower, upper))
            points = reader._fetch_and_decompress_points_of_nodes(nodes)
            x, y = np.asarray(points.x), np.asarray(points.y)
            inside = (x >= lower[0]) & (x < upper[0]) & (y >= lower[1]) & (y < upper[1])
            windows.append(np.column_stack([x[inside], y[inside], np.asarray(points.z)[inside]]
                                           + [np.asarray(points[n])[inside].astype(float) for n in names]))
windows = np.concatenate(windows)


def ordered(rows):
    return rows[np.lexsort(rows.T[::-1])]


assert len(full) > 0, "no points in the checked area"
assert len(windows) == len(full), (len(windows), len(full))
assert np.array_equal(ordered(windows), ordered(full)), "window reads differ from the full read"
print(f"{len(full)} points in {area:g} m: hierarchy window reads equal one full read "
      f"({len(names)} ExtraBytes)")
