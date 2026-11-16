"""The container's memory working set and limit, read from cgroup v2 or v1.

The rules are the R scheduler's (resource_memory_snapshot in
resource_scheduler.R), and tests/test_container_memory.py checks that both give
the same snapshot on the same cgroup files:

- The usage file is the first of memory.current (v2),
  memory/memory.usage_in_bytes and memory.usage_in_bytes (v1); an unreadable
  usage is an error. Without one, the process RSS stands in.
- Only demonstrably clean inactive file pages count as reclaimable:
  min(inactive_file, file) - (dirty + writeback + shmem), at least 0 and at
  most the usage (v1: the hierarchical total_* fields). A memory.stat that is
  missing, malformed, has duplicate keys or lacks a field grants nothing.
- The working set is the usage minus reclaimable pages, at least the process RSS.
- The limit is memory.max (v2) or memory.limit_in_bytes (v1); "max", a value
  of 2^60 or more, or an unreadable file means no limit (None).
"""
from __future__ import annotations

import math
import os

FIELDS_V2 = ("inactive_file", "file", "file_dirty", "file_writeback", "shmem")
FIELDS_V1 = ("total_inactive_file", "total_cache", "total_dirty", "total_writeback", "total_shmem")


def process_rss_bytes() -> float:
    try:
        with open("/proc/self/status") as handle:
            for line in handle:
                if line.startswith("VmRSS:"):
                    return float(line.split()[1]) * 1024
    except (OSError, ValueError, IndexError):
        pass
    return 0.0


def _number(text):
    try:
        value = float(text)
    except ValueError:
        return math.nan
    return value


def _reclaimable(stat_path, raw, fields):
    try:
        with open(stat_path) as handle:
            lines = [line.split() for line in handle.read().splitlines()]
    except OSError:
        return 0.0
    if not lines or any(len(tokens) != 2 for tokens in lines):
        return 0.0
    names = [tokens[0] for tokens in lines]
    if len(set(names)) != len(names) or not all(field in names for field in fields):
        return 0.0
    stats = {name: _number(value) for name, value in lines}
    values = [stats[field] for field in fields]
    if not all(math.isfinite(value) and value >= 0 for value in values):
        return 0.0
    return min(raw, max(0.0, min(values[0], values[1]) - sum(values[2:])))


def _limit(path):
    try:
        with open(path) as handle:
            text = handle.read().strip()
        if text == "max":
            return None
        value = float(text)
    except (OSError, ValueError):
        return None
    return None if not math.isfinite(value) or value >= 2**60 else value


def snapshot(root: str = "/sys/fs/cgroup", rss: float | None = None) -> dict:
    """Working set, usage, reclaimable pages and limit in bytes (see the module rules)."""
    rss = process_rss_bytes() if rss is None else float(rss)
    for usage, limit, fields, source in (
            ("memory.current", "memory.max", FIELDS_V2, "cgroup_v2"),
            ("memory/memory.usage_in_bytes", "memory/memory.limit_in_bytes", FIELDS_V1, "cgroup_v1"),
            ("memory.usage_in_bytes", "memory.limit_in_bytes", FIELDS_V1, "cgroup_v1")):
        path = os.path.join(root, usage)
        if not os.path.exists(path):
            continue
        with open(path) as handle:
            first = handle.readline().strip()
        raw = _number(first)
        if not math.isfinite(raw) or raw < 0:
            raise RuntimeError("cannot read current container memory for task admission")
        reclaimable = _reclaimable(os.path.join(os.path.dirname(path), "memory.stat"), raw, fields)
        return dict(source=source, raw_bytes=raw, reclaimable_bytes=reclaimable,
                    working_set_bytes=max(rss, raw - reclaimable), process_rss_bytes=rss,
                    limit_bytes=_limit(os.path.join(root, limit)))
    return dict(source="process_rss_fallback", raw_bytes=rss, reclaimable_bytes=0.0,
                working_set_bytes=rss, process_rss_bytes=rss, limit_bytes=None)


def working_set_and_limit(root: str = "/sys/fs/cgroup"):
    """(working set, limit) in bytes for task admission; the limit is None without one."""
    measured = snapshot(root)
    return measured["working_set_bytes"], measured["limit_bytes"]
