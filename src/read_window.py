#!/usr/bin/env python3
"""Copy the points of one XY window of a COPC into an uncompressed LAS,
decoding the COPC nodes in parallel processes.

lidR decodes a COPC window on one core; on dense TLS, where only one DTM window
fits in memory at a time, that read is most of the DTM time. The copy keeps the
source's point format, scales, offsets, CRS and extra bytes, and holds every
point of the touched nodes inside the window's XY box (a superset), so lidR
reading it with the catalog's selection and -inside filter gets the same points.
"""
from __future__ import annotations

import argparse
from collections import namedtuple
import copy
import multiprocessing
import os

import laspy
from laspy.copc import Bounds, load_octree_for_query
from laspy.vlrs.vlrlist import VLRList
import numpy as np

SELECTION = (laspy.DecompressionSelection.xy_returns_channel()
             | laspy.DecompressionSelection.Z | laspy.DecompressionSelection.ALL_EXTRA_BYTES)
CopcChunk = namedtuple("CopcChunk", "offset byte_size point_count")
# Points per decoding task: about 0.2 GB of records each, so ten workers stay small.
TASK_POINTS = 2_000_000
WORKER = {}


def start_worker(path):
    WORKER["reader"] = laspy.CopcReader.open(path, decompression_selection=SELECTION)


def decode(task):
    """Worker: the records of a run of COPC nodes inside the box, as bytes."""
    chunks, box = task
    points = WORKER["reader"]._fetch_and_decompress_points_of_nodes(chunks)
    x, y = np.asarray(points.x), np.asarray(points.y)
    inside = (x >= box[0]) & (x <= box[2]) & (y >= box[1]) & (y <= box[3])
    return points.array[inside].tobytes()


def tasks(chunks, box):
    task, size = [], 0
    for chunk in chunks:
        task.append(chunk)
        size += chunk.point_count
        if size >= TASK_POINTS:
            yield task, box
            task, size = [], 0
    if task:
        yield task, box


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--bounds", type=float, nargs=4, required=True, metavar=("XMIN", "YMIN", "XMAX", "YMAX"))
    parser.add_argument("--threads", type=int, required=True)
    args = parser.parse_args()
    with laspy.CopcReader.open(args.input) as reader:
        source = reader.header
        # Two coordinate steps of slack; lidR's -inside test decides the edge.
        slack = 2 * float(max(source.scales[:2]))
        box = (args.bounds[0] - slack, args.bounds[1] - slack, args.bounds[2] + slack, args.bounds[3] + slack)
        nodes = load_octree_for_query(reader.source, reader.copc_info, reader.root_page, query_bounds=Bounds(
            np.array([box[0], box[1], source.mins[2] - 1]), np.array([box[2], box[3], source.maxs[2] + 1])))
    chunks = sorted((CopcChunk(n.offset, n.byte_size, n.point_count) for n in nodes if n.point_count > 0),
                    key=lambda c: c.offset)
    header = copy.deepcopy(source)
    header.vlrs = VLRList([v for v in source.vlrs if v.user_id not in ("copc", "laszip encoded")])
    header.evlrs = VLRList()
    record = header.point_format.size
    rayon = os.environ.get("RAYON_NUM_THREADS")
    os.environ["RAYON_NUM_THREADS"] = "1"
    try:
        with laspy.open(args.output, mode="w", header=header, do_compress=False) as writer, \
                multiprocessing.get_context("spawn").Pool(args.threads, start_worker, (args.input,)) as pool:
            for data in pool.imap_unordered(decode, tasks(chunks, box)):
                if len(data) % record:
                    raise ValueError("decoded records do not match the point format")
                if data:
                    writer.write_points(laspy.PackedPointRecord(
                        np.frombuffer(data, dtype=header.point_format.dtype()), header.point_format))
    finally:
        if rayon is None:
            os.environ.pop("RAYON_NUM_THREADS")
        else:
            os.environ["RAYON_NUM_THREADS"] = rayon


if __name__ == "__main__":
    main()
