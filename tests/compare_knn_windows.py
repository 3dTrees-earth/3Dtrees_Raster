#!/usr/bin/env python3
"""Run the kNN filter with one tree and window by window under a budget that
forces quarter splits, and require identical rejected ordinals per window."""
import json
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/opt/raster")
import knn_denoising  # noqa: E402

source, root = Path(sys.argv[1]), Path(sys.argv[2])
splits = []
quarters = knn_denoising.quarters


def counted(*args, **kwargs):
    splits.append(args[9])
    return quarters(*args, **kwargs)


knn_denoising.quarters = counted
open_points = []
bounded = knn_denoising.bounded_distances


def counted_bounded(reader, query, *args):
    open_points.append(len(query))
    return bounded(reader, query, *args)


knn_denoising.bounded_distances = counted_bounded


def run(name, budget, windowed):
    out = root / name
    out.mkdir()
    os.environ["RASTER_KNN_WINDOWED"] = "1" if windowed else "0"
    sys.argv = ["knn_denoising.py", "--input", str(source), "--order-dimension", "OriginalPointIndex",
                "--neighbors", "8", "--sigma", "2", "--window-size", "60", "--threads", "2",
                "--memory-budget-bytes", str(budget), "--output-dir", str(out)]
    knn_denoising.main()
    return out, json.loads((out / "manifest.json").read_text())


whole, whole_manifest = run("knn-whole", 3 * 1024**3, False)
assert not splits and not open_points
# The fixture's three COPC nodes need about 13 MB to decode, so every read fits
# 16 MiB; ten times the cost per point inside a read caps reads at about 4,500
# points, which splits the windows and sends open points to the bounded search.
knn_denoising.WINDOW_BYTES_PER_POINT *= 10
windows, windows_manifest = run("knn-windows", 528 * 1024**2, True)
assert splits, "the budget did not split any window"
assert open_points, "no point needed a read beyond its window's buffer"
assert whole_manifest["rejected_points"] == windows_manifest["rejected_points"] > 0
cores = {(p["ix"], p["iy"]): p for p in windows_manifest["partitions"]}
for partition in whole_manifest["partitions"]:
    other = cores[(partition["ix"], partition["iy"])]
    assert partition["core_points"] == other["core_points"]
    read = [np.fromfile(d / partition["filename"], "<f8") if (d / partition["filename"]).exists()
            else np.empty(0) for d in (whole, windows)]
    assert np.array_equal(*read), partition["filename"]
print(f"kNN windows equal one tree: {whole_manifest['rejected_points']} rejected, "
      f"{len(splits)} splits down to depth {max(splits) + 1}, "
      f"{sum(open_points)} points resolved beyond the buffer")
