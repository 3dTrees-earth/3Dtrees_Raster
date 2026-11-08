#!/usr/bin/env python3
"""Build the ordered LAZ precursor used for canonical COPC conversion.

COPC physically reorders point records into a spatial hierarchy.  The original
original R workflow is order-sensitive in a few tie-breaking operations, so the
conversion must carry the source record position along with every point.  This
helper performs only that reversible preparation step: it copies every source
dimension unchanged, adds ``OriginalPointIndex``, and deliberately places the
dimensions needed by the raster tool inside lidR's first nine selectively
readable ExtraBytes slots.
"""

from __future__ import annotations

import argparse
import copy
import pathlib

import laspy
import numpy as np


# This UInt64 ExtraByte is the stable key used to restore pre-COPC record order.
DIMENSION_NAME = "OriginalPointIndex"
# These are the three supported segmentation outputs.  Moving them to the front
# keeps them addressable through lidR selectors ``1`` through ``9`` after the
# order key has occupied the first ExtraByte position.
SELECTIVE_DIMENSIONS = (
    "PredInstance",
    "PredInstance_SAT",
    "PredInstance_FM",
)


def extra_bytes_params(
    dimension: laspy.point.dims.DimensionInfo,
) -> laspy.ExtraBytesParams:
    """Clone one ExtraByte definition without changing its storage semantics."""

    # Copy every field that affects decoding; copying values without these
    # metadata fields could silently change scale, offset, or no-data meaning.
    return laspy.ExtraBytesParams(
        name=dimension.name,
        type=dimension.dtype,
        description=dimension.description,
        offsets=dimension.offsets,
        scales=dimension.scales,
        no_data=dimension.no_data,
    )


def indexed_header(source_header: laspy.LasHeader) -> laspy.LasHeader:
    """Return a source-equivalent header with deterministic ExtraByte order."""

    # Deep-copy the complete LAS header so CRS, VLRs, point format, scales, and
    # offsets remain independent from the reader-owned source object.
    output_header = copy.deepcopy(source_header)
    # Record the definitions before removing them; they are re-added below in a
    # streamable order while preserving each original definition exactly.
    original_extra_dimensions = list(source_header.point_format.extra_dimensions)
    # laspy removes ExtraBytes by name rather than by ordinal.
    original_extra_names = [item.name for item in original_extra_dimensions]
    # Remove all existing ExtraBytes so their new ordinal order is explicit.
    output_header.remove_extra_dims(original_extra_names)
    # Put the stable order key first so lidR can always request it as selector 1.
    output_header.add_extra_dim(
        laspy.ExtraBytesParams(name=DIMENSION_NAME, type=np.uint64)
    )
    # Retain the complete source definitions in a name-indexed lookup.
    dimensions_by_name = {
        item.name: item for item in original_extra_dimensions
    }
    # Put supported instance dimensions immediately after OriginalPointIndex.
    ordered_names = [
        name for name in SELECTIVE_DIMENSIONS if name in dimensions_by_name
    ]
    # Append every other source ExtraByte in its original relative order.
    ordered_names.extend(
        name for name in original_extra_names if name not in ordered_names
    )
    # Recreate ExtraBytes from their untouched definitions in the chosen order.
    for name in ordered_names:
        output_header.add_extra_dim(extra_bytes_params(dimensions_by_name[name]))
    # The caller will use this header for a streaming point-for-point copy.
    return output_header


def add_original_point_index(
    source_path: pathlib.Path,
    output_path: pathlib.Path,
    chunk_points: int,
) -> None:
    """Stream-copy ``source_path`` and attach a unique source-row ordinal."""

    # Open as a reader so even very large datasets stay memory-bounded.
    with laspy.open(source_path) as source:
        # Preserve the original header for dimension iteration and final counts.
        source_header = source.header
        # Refuse ambiguous input: an existing order key could have different
        # semantics and must never be silently overwritten.
        if DIMENSION_NAME in source_header.point_format.dimension_names:
            raise ValueError(f"{DIMENSION_NAME} already exists")
        # Construct the equivalent header with the new deterministic ExtraBytes.
        output_header = indexed_header(source_header)
        # Create only the explicit destination parent, including nested test dirs.
        output_path.parent.mkdir(parents=True, exist_ok=True)
        # ``ordinal`` is both the next source-row index and the copied-point count.
        ordinal = 0
        # Write compressed LAZ incrementally to bound peak memory.
        with laspy.open(
            output_path,
            mode="w",
            header=output_header,
            do_compress=True,
        ) as writer:
            # Preserve source order by processing chunks sequentially.
            for points in source.chunk_iterator(chunk_points):
                # Allocate records using the reordered destination header.
                indexed = laspy.ScaleAwarePointRecord.zeros(
                    len(points), header=output_header
                )
                # Copy every standard and ExtraByte dimension by name so values
                # are unchanged even though ExtraByte ordinals moved.
                for dimension in source_header.point_format.dimension_names:
                    indexed[dimension] = points[dimension]
                # Assign a globally unique, zero-based source record index to the
                # current chunk; UInt64 avoids overflow for very large clouds.
                indexed[DIMENSION_NAME] = np.arange(
                    ordinal, ordinal + len(points), dtype=np.uint64
                )
                # Append the fully populated chunk in original source order.
                writer.write_points(indexed)
                # Advance the global index by exactly the number of copied rows.
                ordinal += len(points)
        # Detect short reads or writer-loop mistakes before COPC construction.
        if ordinal != source_header.point_count:
            raise RuntimeError(
                f"indexed {ordinal} of {source_header.point_count} source points"
            )


def main() -> None:
    """Parse the small CLI surface and run the bounded indexing copy."""

    # Positional paths make the helper easy to invoke from the Docker wrapper.
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=pathlib.Path)
    parser.add_argument("output", type=pathlib.Path)
    parser.add_argument("--chunk-points", type=int, default=1_000_000)
    # Convert CLI strings to validated pathlib/integer values.
    args = parser.parse_args()
    # Reject zero/negative chunks because they cannot make streaming progress.
    if args.chunk_points < 1:
        raise ValueError("--chunk-points must be positive")
    # Perform the single intended transformation after all input validation.
    add_original_point_index(args.source, args.output, args.chunk_points)


if __name__ == "__main__":
    # Keep imports side-effect free while making direct script execution work.
    main()
