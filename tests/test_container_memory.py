#!/usr/bin/env python3
"""container_memory.snapshot (Python, the kNN filter) and resource_memory_snapshot
(R, the window scheduler) give the same snapshot on the same cgroup files, for
every case of tests/test_resource_memory.R; and the limit is read as specified."""
import json
import os
import subprocess
import sys
import tempfile

SOURCE = os.environ.get("RASTER_SOURCE", "/opt/raster")
sys.path.insert(0, SOURCE)
import container_memory  # noqa: E402

GIB = 1024**3
STAT_V2 = ("inactive_file", "file", "file_dirty", "file_writeback", "shmem")
STAT_V1 = tuple("total_" + name for name in ("inactive_file", "cache", "dirty", "writeback", "shmem"))


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as handle:
        handle.write(text)


def stat_lines(names, values):
    return "".join(f"{name} {int(value * GIB)}\n" for name, value in zip(names, values))


def cases():
    """(name, files) pairs; each file maps a path below the root to its text."""
    usage = f"{24 * GIB}\n"
    yield "v2", {"memory.current": usage, "memory.stat": stat_lines(STAT_V2, (18, 20, 2, 1, 3))}
    yield "v2 dirty inactive", {"memory.current": usage, "memory.stat": stat_lines(STAT_V2, (10, 20, 8, 0, 0))}
    yield "v2 nothing left", {"memory.current": usage, "memory.stat": stat_lines(STAT_V2, (3, 20, 2, 1, 3))}
    yield "v2 rss floor", {"memory.current": usage, "memory.stat": stat_lines(STAT_V2, (24, 24, 0, 0, 0))}
    for name, text in (("empty", ""), ("one field", "inactive_file 1000\n"), ("nan", "file NaN\n"),
                       ("broken", "broken\n"), ("duplicate", "file 1000\nfile 1000\n")):
        yield f"v2 stat {name}", {"memory.current": usage, "memory.stat": text}
    yield "v2 no stat", {"memory.current": usage}
    yield "no cgroup", {}
    yield "v1", {"memory/memory.usage_in_bytes": usage,
                 "memory/memory.stat": stat_lines(STAT_V1, (18, 20, 2, 1, 3))}
    yield "v1 local field only", {"memory/memory.usage_in_bytes": usage,
                                  "memory/memory.stat": "inactive_file 9999999999\n"}
    yield "v1 at root", {"memory.usage_in_bytes": usage, "memory.stat": stat_lines(STAT_V1, (18, 20, 2, 1, 3))}


R_SNAPSHOT = """
source(file.path(Sys.getenv("RASTER_SOURCE", "/opt/raster"), "resource_scheduler.R"))
args <- commandArgs(trailingOnly = TRUE)
cat(jsonlite::toJSON(lapply(args[-1], function(root)
  tryCatch(resource_memory_snapshot(root, rss = as.numeric(args[[1]])),
           error = function(error) list(error = TRUE))), auto_unbox = TRUE, digits = NA))
"""


def main():
    rss_values = (1 * GIB, 4 * GIB)
    with tempfile.TemporaryDirectory() as scratch:
        roots = []
        for index, (name, files) in enumerate(cases()):
            root = os.path.join(scratch, str(index))
            os.makedirs(root)
            for path, text in files.items():
                write(os.path.join(root, path), text)
            roots.append((name, root))
        bad = os.path.join(scratch, "bad")
        write(os.path.join(bad, "memory.current"), "not-a-counter\n")
        roots.append(("unreadable usage", bad))
        for rss in rss_values:
            r_side = json.loads(subprocess.run(
                ["Rscript", "-e", R_SNAPSHOT, str(rss)] + [root for _, root in roots],
                check=True, capture_output=True, text=True).stdout)
            for (name, root), expected in zip(roots, r_side):
                try:
                    got = container_memory.snapshot(root, rss=rss)
                except RuntimeError:
                    got = {"error": True}
                common = {key: got[key] for key in expected}
                assert common == expected, (name, rss, got, expected)
    # The limit, which only the Python side reads.
    with tempfile.TemporaryDirectory() as root:
        write(os.path.join(root, "memory.current"), f"{GIB}\n")
        for text, limit in (("max\n", None), (f"{8 * GIB}\n", 8 * GIB), (f"{2**62}\n", None), ("junk\n", None)):
            write(os.path.join(root, "memory.max"), text)
            assert container_memory.snapshot(root, rss=0)["limit_bytes"] == limit, (text, limit)
        os.remove(os.path.join(root, "memory.max"))
        assert container_memory.working_set_and_limit(root)[1] is None
    print(f"container memory: Python and R snapshots equal on {len(roots)} cgroup cases")


if __name__ == "__main__":
    main()
