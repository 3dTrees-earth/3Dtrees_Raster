#!/usr/bin/env python3
"""Stream a compact XYZ/ExtraByte LAZ projection and verify ordered values."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import pathlib
import tempfile

import laspy
import numpy as np


def update_hash(digest: "hashlib._Hash", values: np.ndarray) -> None:
    array = np.asarray(values)
    little_endian = array.astype(array.dtype.newbyteorder("<"), copy=False)
    digest.update(little_endian.tobytes(order="C"))


def extra_dimension(header: laspy.LasHeader, name: str):
    for dimension in header.point_format.extra_dimensions:
        if dimension.name == name:
            return dimension
    raise ValueError(f"ExtraByte dimension {name!r} is absent")


def extra_params(dimension) -> laspy.ExtraBytesParams:
    kwargs: dict[str, object] = {
        "name": dimension.name,
        "type": dimension.dtype,
        "description": dimension.description or "",
    }
    if dimension.scales is not None:
        kwargs["scales"] = dimension.scales
    if dimension.offsets is not None:
        kwargs["offsets"] = dimension.offsets
    if dimension.no_data is not None:
        kwargs["no_data"] = dimension.no_data
    return laspy.ExtraBytesParams(**kwargs)


def stored_bounds(header: laspy.LasHeader, bounds: np.ndarray) -> np.ndarray:
    """Compare bounds on the stored lattice, allowing writer rounding noise."""
    return np.rint((bounds - header.offsets) / header.scales)


# Plain LAS/LAZ records are already in source order; a projection made for
# point-noise filtering numbers them so noise decisions can name each record.
ORDER_DIMENSION = "OriginalPointIndex"


def stream_projection(
    source_path: pathlib.Path,
    output_path: pathlib.Path,
    dimension_name: str | list[str],
    chunk_points: int,
    add_original_point_index: bool = False,
) -> dict[str, object]:
    source_digest = hashlib.sha256()
    projected_digest = hashlib.sha256()

    with laspy.open(source_path) as source:
        source_header = source.header
        names = [dimension_name] if isinstance(dimension_name, str) else dimension_name
        if not names or len(names) != len(set(names)):
            raise ValueError("projection requires unique dimension names")
        dimensions = [extra_dimension(source_header, name) for name in names]
        if add_original_point_index and ORDER_DIMENSION in names:
            raise ValueError(f"{ORDER_DIMENSION} already exists in the source")
        output_header = laspy.LasHeader(
            point_format=0,
            version=str(source_header.version),
        )
        output_header.scales = source_header.scales.copy()
        output_header.offsets = source_header.offsets.copy()
        output_header.generating_software = "3dtrees-foreststructure-julia-memory-safe"
        source_crs = source_header.parse_crs()
        if source_crs is not None:
            output_header.add_crs(source_crs, keep_compatibility=False)
        if add_original_point_index:
            output_header.add_extra_dim(
                laspy.ExtraBytesParams(name=ORDER_DIMENSION, type=np.uint64)
            )
        for dimension in dimensions:
            output_header.add_extra_dim(extra_params(dimension))
        output_names = ([ORDER_DIMENSION] if add_original_point_index else []) + names

        output_path.parent.mkdir(parents=True, exist_ok=True)
        descriptor, temporary_name = tempfile.mkstemp(
            prefix=f".{output_path.stem}_",
            suffix=".laz",
            dir=output_path.parent,
        )
        os.close(descriptor)
        temporary_path = pathlib.Path(temporary_name)
        try:
            with laspy.open(
                temporary_path,
                mode="w",
                header=output_header,
                do_compress=True,
            ) as writer:
                ordinal = 0
                for points in source.chunk_iterator(chunk_points):
                    projected = laspy.ScaleAwarePointRecord.zeros(
                        len(points),
                        header=output_header,
                    )
                    projected.X = points.X
                    projected.Y = points.Y
                    projected.Z = points.Z
                    if add_original_point_index:
                        projected.array[ORDER_DIMENSION] = np.arange(
                            ordinal, ordinal + len(points), dtype=np.uint64
                        )
                    ordinal += len(points)
                    for name in names:
                        projected.array[name] = points.array[name]
                    for values in (
                        points.X,
                        points.Y,
                        points.Z,
                        *(points.array[name] for name in names),
                    ):
                        update_hash(source_digest, values)
                    writer.write_points(projected)

            with laspy.open(temporary_path) as projected_source:
                expected_ordinal = 0
                for points in projected_source.chunk_iterator(chunk_points):
                    if add_original_point_index:
                        expected = np.arange(
                            expected_ordinal, expected_ordinal + len(points), dtype=np.uint64
                        )
                        if not np.array_equal(points.array[ORDER_DIMENSION], expected):
                            raise RuntimeError(f"{ORDER_DIMENSION} is not the record order")
                    expected_ordinal += len(points)
                    for values in (
                        points.X,
                        points.Y,
                        points.Z,
                        *(points.array[name] for name in names),
                    ):
                        update_hash(projected_digest, values)
                projected_header = projected_source.header

            if source_digest.digest() != projected_digest.digest():
                raise RuntimeError("ordered XYZ/instance hash differs after projection")
            if source_header.point_count != projected_header.point_count:
                raise RuntimeError("point count differs after projection")
            if not np.array_equal(source_header.scales, projected_header.scales):
                raise RuntimeError("XYZ scales differ after projection")
            if not np.array_equal(source_header.offsets, projected_header.offsets):
                raise RuntimeError("XYZ offsets differ after projection")
            if not np.array_equal(stored_bounds(source_header, source_header.mins),
                                  stored_bounds(projected_header, projected_header.mins)):
                raise RuntimeError("XYZ minimum bounds differ after projection")
            if not np.array_equal(stored_bounds(source_header, source_header.maxs),
                                  stored_bounds(projected_header, projected_header.maxs)):
                raise RuntimeError("XYZ maximum bounds differ after projection")
            if source_header.parse_crs() != projected_header.parse_crs():
                raise RuntimeError("CRS differs after projection")
            if [item.name for item in projected_header.point_format.extra_dimensions] != output_names:
                raise RuntimeError("projection contains an unexpected ExtraByte set")
            os.replace(temporary_path, output_path)
        except BaseException:
            temporary_path.unlink(missing_ok=True)
            raise

    return {
        "source": str(source_path.resolve()),
        "projection": str(output_path.resolve()),
        "dimension": names[0] if len(names) == 1 else names,
        "source_ordinal": [
            item.name for item in source_header.point_format.extra_dimensions
        ].index(names[0])
        + 1,
        "source_datatype": str(dimensions[0].dtype),
        "source_ordinals": {
            name: list(source_header.point_format.extra_dimension_names).index(name) + 1
            for name in names
        },
        "chunk_points": chunk_points,
        "point_count": source_header.point_count,
        "xyz_instance_ordered_sha256": source_digest.hexdigest(),
        "projection_selector": "xyz" + "".join(str(i + 1) for i in range(len(output_names))),
        "original_point_index_added": add_original_point_index,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=pathlib.Path)
    parser.add_argument("--output", required=True, type=pathlib.Path)
    parser.add_argument("--dimension", required=True, action="append")
    parser.add_argument("--chunk-points", type=int, default=250_000)
    parser.add_argument(
        "--add-original-point-index",
        action="store_true",
        help=f"Number the records in source order as a leading {ORDER_DIMENSION} ExtraByte",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.chunk_points < 1:
        raise ValueError("--chunk-points must be positive")
    provenance = stream_projection(
        args.source,
        args.output,
        args.dimension,
        args.chunk_points,
        args.add_original_point_index,
    )
    print(json.dumps(provenance, sort_keys=True))


if __name__ == "__main__":
    main()
