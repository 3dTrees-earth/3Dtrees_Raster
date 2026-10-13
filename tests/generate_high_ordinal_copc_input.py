"""Reorder the small indexed fixture without changing any point or label."""
import copy
import sys

import laspy
import numpy as np

from add_original_point_index import extra_bytes_params

source = laspy.read(sys.argv[1])
header = copy.deepcopy(source.header)
extras = list(header.point_format.extra_dimensions)
header.remove_extra_dims([d.name for d in extras])
for i in range(14):
    header.add_extra_dim(laspy.ExtraBytesParams(name=f"padding_{i}", type=np.uint8))
by_name = {d.name: d for d in extras}
ordered = ["PredInstance", "PredInstance_SAT", "PredInstance_FM", "OriginalPointIndex"]
ordered += [d.name for d in extras if d.name not in ordered]
for name in ordered:
    header.add_extra_dim(extra_bytes_params(by_name[name]))
result = laspy.LasData(header)
result.points = laspy.ScaleAwarePointRecord.zeros(len(source.points), header=header)
for name in source.point_format.dimension_names:
    result[name] = source[name]
result.write(sys.argv[2])
assert list(result.point_format.extra_dimension_names).index("PredInstance") == 14

# Header-only preflight failure fixtures must not accidentally reach terrain work.
from pathlib import Path

result.remove_extra_dims(["PredInstance", "PredInstance_FM"])
result.write(Path(sys.argv[2]).with_name("sat_only.laz"))
result.remove_extra_dim("PredInstance_SAT")
result.write(Path(sys.argv[2]).with_name("no_instances.laz"))
