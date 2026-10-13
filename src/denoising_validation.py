"""Bounded validation of local-SOR accounting and rejected source identities."""

from __future__ import annotations

import argparse
import csv
import json
import sqlite3
import tempfile
from pathlib import Path


def count(value, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"invalid denoising {name}")
    if not 0 <= value <= 2**31 - 1 or int(value) != value:
        raise ValueError(f"invalid denoising {name}")
    return int(value)


def validate(metadata_path: Path, csv_path: Path, point_count: int) -> dict:
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    source_count = count(metadata.get("point_count"), "point_count")
    rejected = count(metadata.get("rejected_points"), "rejected_points")
    retained = count(metadata.get("retained_points"), "retained_points")
    if (source_count != point_count or rejected + retained != source_count
            or metadata.get("enabled") is not True
            or metadata.get("method") != "fixed_window_local_sor"
            or metadata.get("point_order_dimension") != "OriginalPointIndex"
            or metadata.get("coordinate_changes") is not False
            or metadata.get("raw_voxel_counts_filtered") is not False):
        raise ValueError("denoising source/policy accounting mismatch")
    partitions = metadata.get("partitions")
    if not isinstance(partitions, list) or not partitions:
        raise ValueError("denoising partition accounting is missing")
    if len(partitions) != count(metadata.get("partition_count"), "partition_count"):
        raise ValueError("denoising partition count mismatch")
    totals = dict(core_points=0, rejected_points=0, accepted_points=0,
                  insufficient_context_points=0)
    cells = set()
    for part in partitions:
        cell = (part["ix"], part["iy"])
        if cell in cells:
            raise ValueError("duplicate denoising scientific core")
        cells.add(cell)
        values = {key: count(part.get(key), key) for key in totals}
        if (values["accepted_points"] + values["rejected_points"] != values["core_points"]
                or values["insufficient_context_points"] > values["accepted_points"]):
            raise ValueError("denoising core counts do not reconcile")
        for key, value in values.items():
            totals[key] += value
    if (totals["core_points"] != source_count or totals["rejected_points"] != rejected
            or totals["accepted_points"] != retained
            or totals["insufficient_context_points"] != count(
                metadata.get("insufficient_context_points"), "insufficient_context_points")):
        raise ValueError("denoising partition totals do not reconcile")

    # The index stays on disk, even if most observations were rejected. A set
    # of all source IDs in Python would be unbounded for a large merged COPC.
    with tempfile.TemporaryDirectory(prefix="foreststructure_noise_validation_") as scratch:
        with sqlite3.connect(str(Path(scratch) / "ids.sqlite")) as database:
            database.execute("PRAGMA cache_size=-8192")
            database.execute("PRAGMA temp_store=FILE")
            database.execute("CREATE TABLE ids (id INTEGER PRIMARY KEY)")
            rows = 0
            batch = []
            with csv_path.open(newline="", encoding="utf-8") as stream:
                reader = csv.reader(stream)
                if next(reader, None) != ["OriginalPointIndex"]:
                    raise ValueError("invalid denoising CSV columns")
                for row in reader:
                    if len(row) != 1:
                        raise ValueError("invalid denoising CSV row")
                    ordinal = int(row[0])
                    if not 0 <= ordinal < point_count:
                        raise ValueError("denoising source ordinal outside input")
                    batch.append((ordinal,))
                    rows += 1
                    if len(batch) == 10000:
                        database.executemany("INSERT INTO ids VALUES (?)", batch)
                        batch.clear()
                database.executemany("INSERT INTO ids VALUES (?)", batch)
            if rows != rejected:
                raise ValueError("denoising CSV does not reconcile with rejected count")
    return metadata


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata", type=Path)
    parser.add_argument("csv", type=Path)
    parser.add_argument("--point-count", type=int, required=True)
    args = parser.parse_args()
    validate(args.metadata, args.csv, args.point_count)
