#!/usr/bin/env python3
"""Prove that canonical COPC conversion preserved scientific point identity.

The original LAZ and the spatially reordered COPC cannot be compared row by
row.  Instead, both sides are reduced to order-independent fingerprints of the
same logical tuples: ``(OriginalPointIndex, X, Y, Z, instance value)``.  Exact
index-range accounting separately proves that the COPC contains every source
ordinal once.  All reads are chunked, so validation does not require loading a
large point cloud into memory.
"""

from __future__ import annotations

import argparse
import json
import pathlib
from typing import Any

import laspy
import numpy as np
from pyproj import CRS


# This is the stable join key inserted before Untwine spatially reorders points.
INDEX_DIMENSION = "OriginalPointIndex"
# Validate every supported instance dimension that exists in the source cloud.
INSTANCE_DIMENSIONS = (
    "PredInstance",
    "PredInstance_SAT",
    "PredInstance_FM",
)
# Explicit masking gives Python integer additions the same wraparound as UInt64.
MASK64 = (1 << 64) - 1


def mix64(values: np.ndarray) -> np.ndarray:
    """Apply a fast deterministic UInt64 avalanche mix to an array."""

    # Copy before in-place operations so caller-owned point arrays are untouched.
    values = values.astype(np.uint64, copy=True)
    # These SplitMix64 finalizer steps spread nearby LAS integer values across
    # all 64 bits; overflow is intentional and deterministic for np.uint64.
    values ^= values >> np.uint64(30)
    values *= np.uint64(0xBF58476D1CE4E5B9)
    values ^= values >> np.uint64(27)
    values *= np.uint64(0x94D049BB133111EB)
    values ^= values >> np.uint64(31)
    return values


def tuple_hash(
    indexes: np.ndarray,
    points: laspy.ScaleAwarePointRecord,
    dimension: str,
) -> np.ndarray:
    """Hash the identity and scientific values of each point in one chunk."""

    # Seed the accumulator with OriginalPointIndex so point order becomes part
    # of logical identity even though the COPC record sequence differs.
    result = mix64(indexes.astype(np.uint64, copy=False) + np.uint64(0x9E3779B97F4A7C15))
    # Mix raw integer X/Y/Z plus one instance dimension with independent seeds;
    # raw integers avoid floating-point formatting or scale-rounding ambiguity.
    for seed, name in (
        (0x243F6A8885A308D3, "X"),
        (0x13198A2E03707344, "Y"),
        (0xA4093822299F31D0, "Z"),
        (0x082EFA98EC4E6C89, dimension),
    ):
        # Reinterpret signed or unsigned source values modulo 2^64 consistently.
        values = np.asarray(points[name]).astype(np.uint64, copy=False)
        # Fold this field into the current per-point tuple accumulator.
        result ^= mix64(values + np.uint64(seed))
        # Re-avalanche after every field so field permutations cannot cancel.
        result = mix64(result)
    # Return one 64-bit tuple hash per point for commutative aggregation.
    return result


def update_fingerprint(
    fingerprint: dict[str, int], hashes: np.ndarray
) -> None:
    """Merge chunk hashes into two independent order-neutral accumulators."""

    # Sum and XOR are commutative, so COPC physical ordering cannot affect them.
    fingerprint["sum_a"] = (
        fingerprint["sum_a"] + int(np.sum(hashes, dtype=np.uint64))
    ) & MASK64
    fingerprint["xor_a"] ^= int(np.bitwise_xor.reduce(hashes, initial=np.uint64(0)))
    # A second mixed stream makes accidental aggregate collisions substantially
    # less likely than relying on one sum/XOR pair alone.
    secondary = mix64(hashes ^ np.uint64(0x452821E638D01377))
    fingerprint["sum_b"] = (
        fingerprint["sum_b"] + int(np.sum(secondary, dtype=np.uint64))
    ) & MASK64
    fingerprint["xor_b"] ^= int(
        np.bitwise_xor.reduce(secondary, initial=np.uint64(0))
    )


def fingerprint_json(fingerprint: dict[str, int]) -> dict[str, str]:
    """Render UInt64 accumulators as fixed-width, reviewable hexadecimal."""

    return {name: f"{value:016x}" for name, value in fingerprint.items()}


def header_metadata(reader: laspy.LasReader) -> dict[str, Any]:
    """Extract the header fields whose equality is required for parity."""

    # Keep the reader header local so the report records the exact open file.
    header = reader.header
    crs = header.parse_crs()
    # Scales and offsets are included because equal displayed coordinates are
    # insufficient if later software decodes the underlying integers differently.
    return {
        "point_count": int(header.point_count),
        "crs_wkt": crs.to_wkt() if crs else None,
        "scales": header.scales.tolist(),
        "offsets": header.offsets.tolist(),
        "mins": header.mins.tolist(),
        "maxs": header.maxs.tolist(),
        "point_format": int(header.point_format.id),
        "extra_dimensions": [
            item.name for item in header.point_format.extra_dimensions
        ],
    }


def update_stored_bounds(metadata: dict[str, Any], points) -> None:
    """Accumulate actual stored XYZ extrema during the existing bounded scan."""
    if not len(points):
        return
    low = [float(np.min(points[name])) * scale + offset for name, scale, offset
           in zip(("X", "Y", "Z"), metadata["scales"], metadata["offsets"])]
    high = [float(np.max(points[name])) * scale + offset for name, scale, offset
            in zip(("X", "Y", "Z"), metadata["scales"], metadata["offsets"])]
    metadata["stored_mins"] = [min(a, b) for a, b in zip(metadata.get("stored_mins", low), low)]
    metadata["stored_maxs"] = [max(a, b) for a, b in zip(metadata.get("stored_maxs", high), high)]


def headers_compatible(source: dict[str, Any], candidate: dict[str, Any]) -> bool:
    """Require exact stored geometry, quantization and CRS, and accurate candidate bounds.

    Source writers can retain bounds from before integer quantization. Preserve
    those original bounds as evidence, but compare actual extrema from the scan.
    No numerical tolerance is applied to point coordinates or candidate bounds.
    """
    same = all(source.get(k) == candidate.get(k) for k in
               ("point_count", "scales", "offsets", "stored_mins", "stored_maxs"))
    return (same and crs_compatible(source, candidate) and
            all(candidate.get("stored_" + k) == candidate[k] for k in ("mins", "maxs")))


def crs_compatible(source: dict[str, Any], candidate: dict[str, Any]) -> bool:
    """Compare CRS meaning, preserving axis order, rather than WKT spelling."""
    left, right = source.get("crs_wkt"), candidate.get("crs_wkt")
    if left == right:
        return True
    return bool(left and right and CRS.from_wkt(left).equals(CRS.from_wkt(right)))


def source_fingerprint(
    path: pathlib.Path,
    dimensions: list[str],
    chunk_points: int,
) -> tuple[dict[str, Any], dict[str, dict[str, int]]]:
    """Fingerprint the ordered source using implicit zero-based row indices."""

    # Maintain an independent four-lane fingerprint for each scientific field.
    fingerprints = {
        dimension: {"sum_a": 0, "xor_a": 0, "sum_b": 0, "xor_b": 0}
        for dimension in dimensions
    }
    # The source's physical row order defines OriginalPointIndex.
    ordinal = 0
    # Stream rather than materialize the complete source cloud.
    with laspy.open(path) as reader:
        # Capture comparison metadata before reading any chunks.
        metadata = header_metadata(reader)
        # Fail early if the requested scientific values do not exist.
        for dimension in dimensions:
            if dimension not in reader.header.point_format.dimension_names:
                raise ValueError(f"source is missing {dimension}")
        # Iterate sequentially because the implicit ordinal depends on row order.
        for points in reader.chunk_iterator(chunk_points):
            update_stored_bounds(metadata, points)
            # Construct exactly the indices that the indexing copy should add.
            indexes = np.arange(ordinal, ordinal + len(points), dtype=np.uint64)
            # Bind XYZ and each requested instance value to those indices.
            for dimension in dimensions:
                update_fingerprint(
                    fingerprints[dimension],
                    tuple_hash(indexes, points, dimension),
                )
            # Advance by actual rows read, not by configured chunk capacity.
            ordinal += len(points)
    # Detect incomplete reads before accepting matching partial fingerprints.
    if ordinal != metadata["point_count"]:
        raise RuntimeError(
            f"read {ordinal} of {metadata['point_count']} source points"
        )
    return metadata, fingerprints


def candidate_fingerprint(
    path: pathlib.Path,
    dimensions: list[str],
    chunk_points: int,
) -> tuple[dict[str, Any], dict[str, dict[str, int]], dict[str, Any]]:
    """Fingerprint reordered COPC points and audit index completeness exactly."""

    # Use the same independent fingerprint shape as the source side.
    fingerprints = {
        dimension: {"sum_a": 0, "xor_a": 0, "sum_b": 0, "xor_b": 0}
        for dimension in dimensions
    }
    # COPC chunks arrive in spatial order; all aggregation below is order-neutral.
    with laspy.open(path) as reader:
        # Capture header evidence for exact source/candidate comparison.
        metadata = header_metadata(reader)
        # Read the actual declared dimensions once for all requirements.
        available_dimensions = set(reader.header.point_format.dimension_names)
        # Both the scientific values and the inserted order key must survive.
        for required in (*dimensions, INDEX_DIMENSION):
            if required not in available_dimensions:
                raise ValueError(f"candidate is missing {required}")
        # The candidate header defines the only valid index interval [0, n).
        point_count = metadata["point_count"]
        # One byte per expected index gives an exact missing/duplicate audit while
        # using far less memory than retaining all point records.
        seen = np.zeros(point_count, dtype=np.bool_)
        # Counters provide useful diagnostics instead of a single Boolean result.
        read_count = 0
        minimum_index = None
        maximum_index = None
        duplicate_indexes = 0
        out_of_range_indexes = 0
        # Process the spatially ordered candidate in bounded chunks.
        for points in reader.chunk_iterator(chunk_points):
            update_stored_bounds(metadata, points)
            # Decode the persisted order key as UInt64 without changing values.
            indexes = np.asarray(points[INDEX_DIMENSION], dtype=np.uint64)
            # Count every candidate record, including invalid-index records.
            read_count += len(indexes)
            # Empty chunks cannot update ranges, uniqueness, or fingerprints.
            if len(indexes) == 0:
                continue
            # Track observed extrema for a concise contiguous-range proof.
            minimum_index = min(
                int(indexes.min()),
                minimum_index if minimum_index is not None else MASK64,
            )
            maximum_index = max(
                int(indexes.max()),
                maximum_index if maximum_index is not None else 0,
            )
            # Separate invalid indices before using values to address ``seen``.
            in_range = indexes < np.uint64(point_count)
            # Retain the exact number of candidate rows outside [0, n).
            out_of_range_indexes += int(np.count_nonzero(~in_range))
            # Safe int64 indexing is valid after the UInt64 range filter.
            valid_indexes = indexes[in_range].astype(np.int64, copy=False)
            if len(valid_indexes):
                # Count duplicates within this chunk.
                unique_indexes, counts = np.unique(valid_indexes, return_counts=True)
                duplicate_indexes += int(np.sum(counts - 1))
                # Count duplicates previously observed in earlier chunks.
                duplicate_indexes += int(np.count_nonzero(seen[unique_indexes]))
                # Mark all unique valid indices from this chunk as present.
                seen[unique_indexes] = True
            # Fingerprint candidate tuples with their persisted source indices.
            for dimension in dimensions:
                update_fingerprint(
                    fingerprints[dimension],
                    tuple_hash(indexes, points, dimension),
                )
    # Assemble exact completeness evidence separately from probabilistic hashes.
    index_validation = {
        "read_count": read_count,
        "minimum": minimum_index,
        "maximum": maximum_index,
        "duplicate_indexes": duplicate_indexes,
        "out_of_range_indexes": out_of_range_indexes,
        "missing_indexes": int(np.count_nonzero(~seen)),
        # All clauses are required: correct count and extrema alone would not
        # detect a missing index replaced by a duplicate interior index.
        "unique_complete_range": (
            read_count == metadata["point_count"]
            and minimum_index == 0
            and maximum_index == metadata["point_count"] - 1
            and duplicate_indexes == 0
            and out_of_range_indexes == 0
            and not np.any(~seen)
        ),
    }
    return metadata, fingerprints, index_validation


def main() -> None:
    """Run source/candidate validation and write a durable JSON verdict."""

    # The wrapper supplies source, unpublished COPC candidate, and report paths.
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=pathlib.Path)
    parser.add_argument("candidate", type=pathlib.Path)
    parser.add_argument("report", type=pathlib.Path)
    parser.add_argument("--dimension", action="append", default=None)
    parser.add_argument("--chunk-points", type=int, default=1_000_000)
    # Convert all CLI values before any large point-cloud reads begin.
    args = parser.parse_args()
    # A non-positive chunk size would make the streaming contract invalid.
    if args.chunk_points < 1:
        raise ValueError("--chunk-points must be positive")

    # Discover the supported dimensions that actually exist unless the caller
    # explicitly repeats ``--dimension`` to select a narrower validation set.
    with laspy.open(args.source) as source_reader:
        source_dimensions = set(source_reader.header.point_format.dimension_names)
        dimensions = args.dimension or [
            name for name in INSTANCE_DIMENSIONS if name in source_dimensions
        ]
    # A geometry-only comparison cannot prove ForestStructure science parity.
    if not dimensions:
        raise ValueError("source has no supported instance dimensions to validate")

    # Compute the ordered-source evidence first.
    source_metadata, source_digests = source_fingerprint(
        args.source, dimensions, args.chunk_points
    )
    # Compute equivalent COPC evidence plus its exact index audit.
    candidate_metadata, candidate_digests, index_validation = candidate_fingerprint(
        args.candidate, dimensions, args.chunk_points
    )
    # Preserve raw header equality as evidence, separate from stored geometry.
    header_equal = all(
        source_metadata[name] == candidate_metadata[name]
        for name in ("point_count", "scales", "offsets", "mins", "maxs")
    )
    header_compatible = headers_compatible(source_metadata, candidate_metadata)
    # Report individual dimension verdicts so one corrupt field is diagnosable.
    fingerprint_equal_by_dimension = {
        dimension: source_digests[dimension] == candidate_digests[dimension]
        for dimension in dimensions
    }
    # Overall tuple parity requires every requested instance dimension.
    fingerprint_equal = all(fingerprint_equal_by_dimension.values())
    # Publication is allowed only when metadata, tuple identity, and exact index
    # uniqueness all agree.
    valid = (
        header_compatible
        and fingerprint_equal
        and index_validation["unique_complete_range"]
    )
    # Persist both the verdict and enough evidence for independent review.
    report = {
        "status": "valid" if valid else "invalid",
        "comparison": (
            "streaming dual-64-bit tuple fingerprints over "
            "(OriginalPointIndex,X,Y,Z,dimension) with exact index uniqueness"
        ),
        "dimensions": dimensions,
        "header_equal": header_equal,
        "header_compatible": header_compatible,
        "crs_equivalent": crs_compatible(source_metadata, candidate_metadata),
        "bounds_comparison": "exact extrema of stored source and candidate points; candidate header must match stored extrema",
        "fingerprint_equal": fingerprint_equal,
        "fingerprint_equal_by_dimension": fingerprint_equal_by_dimension,
        "source_fingerprints": {
            dimension: fingerprint_json(source_digests[dimension])
            for dimension in dimensions
        },
        "candidate_fingerprints": {
            dimension: fingerprint_json(candidate_digests[dimension])
            for dimension in dimensions
        },
        "index_validation": index_validation,
        "source": source_metadata,
        "candidate": candidate_metadata,
    }
    # Create the requested report parent but never the COPC output location here.
    args.report.parent.mkdir(parents=True, exist_ok=True)
    # A trailing newline keeps the report friendly to standard text tooling.
    args.report.write_text(json.dumps(report, indent=2) + "\n")
    # The shell wrapper relies on a nonzero exit to withhold invalid COPC files.
    if not valid:
        raise SystemExit(f"indexed COPC validation failed; see {args.report}")
    # Emit one concise success record; detailed evidence remains in JSON.
    print(
        "indexed COPC streaming validation passed: "
        f"points={source_metadata['point_count']} dimensions={','.join(dimensions)}"
    )


if __name__ == "__main__":
    # Avoid running validation when imported by focused unit tests.
    main()
