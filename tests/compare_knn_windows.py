#!/usr/bin/env python3
"""Run the kNN filter with one tree and window by window in parallel parts,
and require identical rejected ordinals per window."""
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

source, root = Path(sys.argv[1]), Path(sys.argv[2])


def run(name, budget, windowed, threads):
    out = root / name
    out.mkdir()
    result = subprocess.run(
        [sys.executable, "/opt/raster/knn_denoising.py", "--input", str(source),
         "--order-dimension", "OriginalPointIndex", "--neighbors", "8", "--sigma", "2",
         "--window-size", "60", "--threads", str(threads), "--memory-budget-bytes", str(budget),
         "--output-dir", str(out)],
        env=dict(os.environ, RASTER_KNN_WINDOWED="1" if windowed else "0"),
        check=True, capture_output=True, text=True)
    return out, json.loads((out / "manifest.json").read_text()), result.stdout


whole, whole_manifest, whole_log = run("knn-whole", 3 * 1024**3, False, 2)
assert "KNN_WINDOWS" not in whole_log
# The fixture's windows run in four worker processes. Its root COPC node holds
# most points, so quarters would not be cheaper and stay unsplit; splitting is
# checked on real clouds (2145 and others) outside this suite.
windows, windows_manifest, log = run("knn-windows", 2 * 1024**3, True, 4)
stats = dict(re.findall(r"(\w+)=(\d+)", log.split("KNN_WINDOWS", 1)[1].splitlines()[0]))
assert int(stats["workers"]) == 4 and int(stats["parts"]) >= len(whole_manifest["partitions"]), stats
assert int(stats["far_points"]) > 0, "no point needed a read beyond its part's buffer"
assert whole_manifest["rejected_points"] == windows_manifest["rejected_points"] > 0
cores = {(p["ix"], p["iy"]): p for p in windows_manifest["partitions"]}
for partition in whole_manifest["partitions"]:
    other = cores[(partition["ix"], partition["iy"])]
    assert partition["core_points"] == other["core_points"]
    read = [np.fromfile(d / partition["filename"], "<f8") if (d / partition["filename"]).exists()
            else np.empty(0) for d in (whole, windows)]
    assert np.array_equal(*read), partition["filename"]
print(f"kNN windows equal one tree: {whole_manifest['rejected_points']} rejected, "
      f"{stats['parts']} parts down to depth {stats['deepest']}, "
      f"{stats['far_points']} points resolved beyond the buffer")
