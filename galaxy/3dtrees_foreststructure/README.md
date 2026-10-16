# Galaxy wrapper

`foreststructure.xml` wraps the production runner for Galaxy. It pins
`ghcr.io/3dtrees-earth/3dtrees_raster:@TOOL_VERSION@`; publish that
image tag (push the matching `v<version>` Git tag) before installing the tool.

Run the tests from the repository root:

```bash
planemo lint --fail_level warn galaxy/3dtrees_foreststructure
planemo test --docker galaxy/3dtrees_foreststructure/foreststructure.xml
```

## Test data

The fixtures are the synthetic 40 m patch from
`tests/generate_julia_memory_safe_fixture.R` (13,681 points, EPSG:32632,
distinct `PredInstance`, `PredInstance_SAT` and `PredInstance_FM` values):

- `forest_patch.laz`: the patch after `tests/add_original_point_index.py`.
- `forest_patch.copc.laz`: that file written by PDAL 2.10.2 `writers.copc`
  (`readers.las.use_eb_vlr=true`, `extra_dims=all`, `forward=all`), so it keeps
  `OriginalPointIndex` and the 0.001 m scale.
- `forest_patch_aoi.geojson`: the western half of the patch (two tiles).

With default SOR the COPC case flags two points; the expected values in the
tests come from the PR #5 runner.
