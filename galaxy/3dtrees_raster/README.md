# Galaxy wrapper

`raster.xml` is the Galaxy tool **3Dtrees: Raster** (tool id `3dtrees_raster`).
The canonical copy lives in
[3dtrees-earth/galaxytools](https://github.com/3dTrees-earth/galaxytools) under
`tools/3dtrees_raster`. This directory mirrors it, with `.shed.yml` and the test
data, so the wrapper and the runner can be tested together. Change the wrapper in
galaxytools first, then copy `raster.xml` and `.shed.yml` here unchanged.

The wrapper pins `ghcr.io/3dtrees-earth/3dtrees_raster:@TOOL_VERSION@`. Publish
that image tag (see [Releases](../../README.md#releases)) before installing a new
tool version.

Run the tests from the repository root:

```bash
planemo lint --fail_level warn galaxy/3dtrees_raster
planemo test --docker galaxy/3dtrees_raster/raster.xml
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

With default point-noise filtering the COPC case flags two points.
