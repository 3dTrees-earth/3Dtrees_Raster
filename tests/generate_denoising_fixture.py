"""Add known noise and a supported tall crown to the small Julia fixture.

This is a synthetic test only. The source positions, labels and record order
are saved before canonical COPC conversion so later tests can verify identity.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import laspy
import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("patch", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("manifest", type=Path)
    args = parser.parse_args()
    source = laspy.read(args.patch)
    labels = ("PredInstance", "PredInstance_SAT", "PredInstance_FM")
    # A coherent crown is deliberately above the former 70 m metric ceiling.
    crown = np.array(np.meshgrid(
        np.arange(1010, 1012.01, 0.2),
        np.arange(2010, 2012.01, 0.2),
        np.arange(184, 186.01, 0.2), indexing="ij"
    )).reshape(3, -1).T
    outliers = np.array([[1017.25, 2017.25, 350], [1030.25, 2030.25, -50]])
    ordinary_count = 4 * len(source.points)
    count = ordinary_count + len(crown) + len(outliers)
    result = laspy.LasData(copy.deepcopy(source.header))
    result.points = laspy.ScaleAwarePointRecord.zeros(count, header=result.header)
    for i, (dx, dy) in enumerate(((0, 0), (40, 0), (0, 40), (40, 40))):
        block = slice(i * len(source.points), (i + 1) * len(source.points))
        result.points.array[block] = source.points.array
        result.X[block] += round(dx / source.header.scales[0])
        result.Y[block] += round(dy / source.header.scales[1])
        for dimension in labels:
            values = np.asarray(source[dimension]).copy()
            values[values > 0] += 1000 * i
            result[dimension][block] = values
    # Preserve the earlier fixture geometry: one neutral ground point extends
    # max XY by one stored tick. Exact terminal catalog boundaries are exercised
    # separately by test_segment_catalog_boundaries.R.
    edge_index = 3 * len(source.points) + np.flatnonzero(source["PredInstance"] == 0)[-1]
    result.X[edge_index] += 1
    result.Y[edge_index] += 1
    added = np.vstack((crown, outliers))
    result.x[ordinary_count:] = added[:, 0]
    result.y[ordinary_count:] = added[:, 1]
    result.z[ordinary_count:] = added[:, 2]
    result.return_number[ordinary_count:] = 1
    result.number_of_returns[ordinary_count:] = 1
    result.classification[ordinary_count:] = 1
    for dimension, label in zip(labels, (11, 101, 301)):
        result[dimension][ordinary_count:] = label
        result[dimension][-1] = 0  # The low outlier is not a tree instance.
    result.write(args.output)
    # add_original_point_index uses zero-based source-row indices.
    manifest = dict(
        fixture_version=2,
        points=count,
        high_outlier=dict(original_point_index=count - 2, xyz=outliers[0].tolist()),
        low_outlier=dict(original_point_index=count - 1, xyz=outliers[1].tolist()),
        supported_crown_indices=list(range(ordinary_count, ordinary_count + len(crown))),
        supported_instance_ids=dict(zip(labels, (11, 101, 301))),
        expected_supported_height_m=[80, 90],
        expected_high_outlier_minimum_height_m=200,
        boundary_fixture_policy="one neutral ground point extends max XY by one stored-coordinate tick; exact edges have a separate regression",
    )
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
