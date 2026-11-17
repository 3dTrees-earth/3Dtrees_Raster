#!/usr/bin/env python3
"""Crash-safety check of the kNN filter's part admission on a real COPC.

Runs the windowed filter (knn_denoising.run) with one injected fault and
compares the flagged OriginalPointIndex values with a reference noise list
(`raster_noise_points.csv` of a normal run). The faults are injected here, not
through production switches:

  plain    no fault (use it under a small container memory cap)
  low:F    admission estimates scaled by F (0.2: five times too low)
  killed   the first part's worker is killed with SIGKILL, as by the kernel

Usage: knn_crash_safety.py CASE COPC REFERENCE_CSV OUTPUT_DIR THREADS BUDGET_BYTES
Prints "identical" and exits 0 when the flagged points equal the reference.
"""
import glob
import os
import signal
import sys

import numpy as np

sys.path.insert(0, os.environ.get("RASTER_SOURCE", "/opt/raster"))
import knn_denoising  # noqa: E402

MARKER_ENV = "KNN_CRASH_SAFETY_MARKER"


def part_killed_once(task):
    """Worker: the first part to get here dies as if the kernel killed it."""
    try:
        os.close(os.open(os.environ[MARKER_ENV], os.O_CREAT | os.O_EXCL | os.O_WRONLY))
        os.kill(os.getpid(), signal.SIGKILL)
    except FileExistsError:
        pass
    return knn_denoising.part_distances(task)


def main():
    case, copc, reference, output, threads, budget = sys.argv[1:7]
    args = knn_denoising.parse_args([
        "--input", copc, "--order-dimension", "OriginalPointIndex", "--neighbors", "8",
        "--sigma", "2", "--window-size", "60", "--threads", threads,
        "--memory-budget-bytes", budget, "--output-dir", output])
    options = {}
    if case.startswith("low:"):
        options["admission_scale"] = float(case.split(":", 1)[1])
    elif case == "killed":
        os.environ[MARKER_ENV] = os.path.join(output, "killed.marker")
        options["work"] = part_killed_once
    elif case != "plain":
        raise SystemExit(f"unknown case {case}")
    # Always the windowed path, the one with part admission.
    os.environ["RASTER_KNN_WINDOWED"] = "1"
    knn_denoising.run(args, **options)
    if case == "killed" and not os.path.exists(os.environ[MARKER_ENV]):
        raise SystemExit("no worker was killed")
    files = glob.glob(os.path.join(output, "core_*.bin"))
    got = np.sort(np.concatenate([np.fromfile(f, "<f8") for f in files] or [np.empty(0)]))
    ref = np.sort(np.loadtxt(reference, delimiter=",", skiprows=1, usecols=0, ndmin=1))
    if not np.array_equal(got, ref):
        raise SystemExit(f"DIFFERENT flagged={len(got)} reference={len(ref)}")
    print(f"identical: {len(got)} flagged points")


if __name__ == "__main__":
    main()
