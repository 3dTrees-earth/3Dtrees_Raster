# 3Dtrees: Raster

3Dtrees: Raster turns one segmented forest point cloud into raster and tile
products of forest structure:

- a digital terrain model (DTM), a canopy height model (CHM), a digital
  surface model (DSM),
- full-resolution voxel point counts as a multiband GeoTIFF, a CSV and a COPC
  point cloud with one point per occupied voxel,
- fixed-size tile metrics (vegetation density, height distribution, vertical
  complexity, canopy roughness, gaps and tree statistics) for each tree
  instance dimension,
- point-noise records and provenance for every run.

It reads every point (no subsampling) and never modifies the input file. It
runs as a Docker image, and is the raster step of the
[3Dtrees](https://github.com/3dTrees-earth) platform.

| | |
| --- | --- |
| Image | `ghcr.io/3dtrees-earth/3dtrees_raster:1.0.0` (public) |
| Licence | GPL-3.0 |

## Quick start

Put the point cloud (and optionally an area of interest) in one folder and run:

```bash
docker run --rm --network none \
  --cpus 8 --memory 30g --memory-swap 30g \
  --user "$(id -u):$(id -g)" \
  --volume "$PWD/data:/data" \
  ghcr.io/3dtrees-earth/3dtrees_raster:1.0.0 \
  --point-cloud /data/input.copc.laz \
  --instance-dimension PredInstance_SAT \
  --aoi-geojson /data/aoi.geojson \
  --output-dir /data/output \
  --temp-dir /data/work
```

- Coordinates are taken as metres. When the point cloud has a CRS, it must be
  projected in metres and every raster, vector and the voxel COPC carry it; a
  point cloud without CRS gives the same products without CRS. It needs at least
  one of the instance dimensions `PredInstance`, `PredInstance_SAT` (SegmentAnyTree),
  `PredInstance_FM` (ForestMamba) or `PredInstance_RCT` (RayCloudTools). Without
  `--instance-dimension` every one of them present is processed; with it
  (repeat the option or separate names with commas) exactly the named ones are,
  and a missing one fails the run before any analysis.
- Point-noise filtering is off by default. `--denoise knn` turns it on, on a
  canonical COPC and on plain LAS/LAZ. A plain file without `OriginalPointIndex`
  is numbered in its own record order, which costs one extra read and a scratch
  copy of XYZ and the instance dimensions.
- A COPC is recognised from its content (the `copc` info VLR), whatever its file name.
- Without `--aoi-geojson` the whole point-cloud extent is analysed.
- Without `--output-dir` and `--temp-dir`, outputs go to the working directory
  and scratch to `TMPDIR`.

The run publishes its outputs only after they validate as a complete set; a
failed run leaves no partial results.

## Parameters

### Scientific parameters

The defaults are the published 3Dtrees settings. Keep them to stay comparable
with other 3Dtrees products. Every run records the values it used under
`scientific_parameters` in its run provenance.

| Option | Default | Meaning | Allowed values |
| --- | --- | --- | --- |
| `--tile-size` | `20` | Width of the square analysis tiles in metres | > 0 |
| `--grid-search-step` | `0.5` | Step of the offset search that places the tile grid inside an AOI, in metres | > 0 and ≤ tile size |
| `--voxel-resolution` | `0.2` | Voxel edge length in metres for the tile and tree metrics | Must divide 1 m evenly, e.g. 0.1, 0.2, 0.25, 0.5, 1 |
| `--voxel-count-resolution` | `0.5` | Voxel edge length in metres of the voxel point counts, horizontal and vertical | > 0 |
| `--dtm-resolution` | `1` | DTM resolution in metres | Must divide 60 m evenly, e.g. 0.5, 1, 2, 5 |
| `--chm-resolution` | `0.5` | CHM resolution in metres | Must divide 1 m evenly |
| `--ptd-resolution` | `20` | Seed resolution in metres of the progressive TIN densification (PTD) ground classification | > 0 |
| `--vegetation-minimum-height` | `0.5` | Voxels at or below this height (m) are not vegetation for the tile metrics | ≥ 0 |
| `--gap-height-threshold` | `3` | CHM cells below this height (m) count as canopy gap | ≥ 0 |
| `--minimum-tree-voxels` | `100` | Occupied metric voxels a segment needs to count as a tree | Whole number ≥ 1 |
| `--apex-minimum-height` | `3` | A tree's highest voxel must be strictly above this height (m) | ≥ 0 |
| `--minimum-tree-thickness` | `0.5` | Smallest PCA extent (m) a tree must reach | ≥ 0 |
| `--minimum-occupied-layers` | `3` | One-metre height layers a tree must occupy | Whole number ≥ 1 |

The size rules exist because the CHM and segment passes read the cloud in
whole-metre windows and the DTM in 60 m chunks; a cell must never straddle a
window edge. The runner rejects other values before reading any points.

### Point-noise filtering

| Option | Default | Meaning |
| --- | --- | --- |
| `--denoise` | `none` | `knn` runs the kNN distance filter before the terrain model; `none` (the default) skips it |
| `--denoise-neighbors` | `8` | Nearest neighbours for each point's mean distance, counting the point itself as `lidR::sor()` does |
| `--denoise-sigma` | `2` | Points farther than the whole cloud's mean distance plus this many standard deviations are flagged |
| `--denoise-window-size` | `60` | Width in metres of the windows the noise list is stored in; does not change results |

See [Point-noise filter](#point-noise-filter)
for the method and its limits.

### Run options

| Option | Meaning |
| --- | --- |
| `--point-cloud` | Input LAS/LAZ or canonical COPC, recognised from its content (required) |
| `--instance-dimension` | Instance dimensions to process (default: every supported one present) |
| `--aoi-geojson` | Optional area of interest in the point cloud's XY coordinates; not reprojected |
| `--original-point-cloud` | Optional original LAS/LAZ of a COPC, used only for header and provenance checks |
| `--output-dir`, `--temp-dir` | Output and scratch directories |
| `--raster-format` | `COG` (default) or `GTiff` for the DTM, CHM, DSM and voxel-count rasters |
| `--version` | Print the tool version and exit |
| `--memory-budget-gib` | Optional lower internal memory budget; the container or Galaxy limit is used otherwise |
| `--target-chunk-points` | Approximate points per COPC CHM/segment read window (default 50,000,000; `0` uses fixed windows) |
| `--voxel-count-chunk-points` | Maximum records per voxel-count batch (default 2,000,000) |
| `--sensor`, `--country` | Optional labels copied into the result rows (default `NA`) |

CPU and memory come from `GALAXY_SLOTS`/`GALAXY_MEMORY_MB` or the container's
limits; see [Resources](#resources) for the operator variables. `--dataset-id` and `--export-voxel-counts` are still accepted
for old commands but are not needed.

## Outputs

Per requested instance dimension `<dimension>`, one file in each of three
folders:

| File | Contents |
| --- | --- |
| `tile_metrics/<dimension>_tile_metrics.csv` | One row of [tile metrics](#tile-metrics) per tile |
| `segment_diagnostics/<dimension>_segment_diagnostics.csv` | One row per segment with the tree-acceptance checks, the apex position (`apex_x`, `apex_y`, source CRS) and the segment's crown area (`crown_area_m2`, occupied voxel columns × voxel area) and voxel volume (`voxel_volume_m3`, occupied voxels × voxel volume), the per-segment values the tile `tree_*` aggregates are computed from |
| `tiles/<dimension>_tiles.geojson` | The tile metrics as polygons |

Shared by all dimensions, each starting with the prefix `raster`:

| File | Contents |
| --- | --- |
| `raster_dtm.tif` | DTM GeoTIFF (Cloud Optimized by default) |
| `raster_chm.tif` | CHM GeoTIFF (Cloud Optimized by default): height above the DTM of the highest voxel per cell |
| `raster_dsm.tif` | [DSM](#dsm) on the CHM grid: absolute height of the highest point per cell |
| `raster_tiles.png` | Tile overview map |
| `raster_voxel_counts.csv` | Occupied voxels: `x_center,y_center,height_lower_m,height_upper_m,point_count` |
| `raster_voxel_counts.tif` | UInt32 point counts; with voxel size *s*, band *b* holds heights `[(b-1)·s, b·s)` m above the DTM |
| `raster_voxel_counts_validity.tif` | 0 excluded (outside the AOI or on its boundary), 1 valid, 2 invalid DTM |
| `raster_voxel_counts.copc.laz` | [Voxel points](#voxel-points-copc): one COPC point per occupied voxel at DTM height plus the band centre, with `point_count` and `height_above_dtm` |
| `raster_voxel_counts.json` | Grid, CRS, source and DTM hashes, policies and point accounting |
| `raster_noise_points.csv` | With `--denoise knn` only: one rejected `OriginalPointIndex` per row (the 0-based record position for a plain LAS/LAZ without that field) |
| `raster_denoising.json` | With `--denoise knn` only: parameters and per-window point accounting |
| `raster_run.json` | Run provenance: inputs, parameters, versions, resources |
| `raster_performance.csv` | Stage timings and memory observations |
| `raster_aoi_conversion.json` | The analysed footprint and its AOI conversion checks |

Voxel counts include every recorded point, also those flagged as noise. Zero
means no recorded points, not absence of vegetation; read the validity raster
before interpreting missing CSV rows. See
[Voxel counts](#voxel-counts).

### DSM

The DSM shares the CHM's grid: resolution (`--chm-resolution`), extent, Float32
data type, nodata value (-999999), tiling, compression and COG layout.

`raster_dsm.tif` is computed natively in the CHM pass, from the same
points: denoised when the kNN filter is on, with a finite height of at least 0 m above the
DTM. Each cell holds the highest **absolute** point Z (the input's vertical
datum, not normalized) among those points. Points are assigned to cells through
the centre of their `--voxel-resolution` voxel column, exactly as the CHM
assigns its voxels. So the DSM has a value exactly where the CHM has one.
DSM − DTM differs from the CHM by the CHM's voxel rounding (at most half a
voxel) and by the DTM being sampled at the points rather than at cell centres.
The DSM is not DTM + CHM: that sum would carry the voxel rounding and the DTM
interpolation into every cell.

The run provenance records it under `surface_rasters`. Before publication the run
checks that the DSM's empty cells equal the CHM's.

### Tile metrics

Within an AOI the tile grid is shifted in `--grid-search-step` increments to fit
the most complete tiles, and only whole tiles inside the AOI are reported.
Without an AOI, tiles cover the full extent and edge tiles can be partial.
The metric columns are computed from occupied voxels of size
`--voxel-resolution`, with heights above the DTM. *H* is the shared height
domain `metric_height_upper_m` (the tallest accepted voxel, rounded up to a
whole metre); vegetation voxels are those above `--vegetation-minimum-height`.

| Column | Meaning |
| --- | --- |
| `file`, `sensor`, `country` | Input file name and the optional labels |
| `tile_id`, `tile_xmin`, `tile_ymin` | Tile number and its lower-left corner |
| `edge_tile` | `TRUE` when any of the eight neighbouring tiles is missing |
| `n_tiles_plot` | Number of tiles in the run |
| `vox_filled`, `vox_total` | Occupied vegetation voxels, and all voxels of the tile volume between the vegetation minimum height and *H* |
| `metric_height_upper_m` | Height domain *H* in metres |
| `veg_density` | `vox_filled / vox_total` |
| `zsd`, `zskew`, `zkurt`, `zq90` | Standard deviation, skewness, kurtosis and 90th percentile of vegetation voxel heights |
| `box_dim_fixed` | Box-counting (fractal) dimension of the vegetation voxels |
| `vci` | Vertical complexity index: Shannon entropy of vegetation voxels over one-metre layers, divided by log(*H*) |
| `rumple` | Rumple index of the tile CHM (surface area / ground area) |
| `gap_fraction` | Share of tile CHM cells below `--gap-height-threshold` |
| `chm_sd`, `chm_cv` | Standard deviation and coefficient of variation of the tile CHM |
| `height_max`, `height_mean` | Maximum and mean voxel height |
| `n_seg_total`, `n_trees` | Segments with their apex in the tile, and those accepted as trees |
| `tree_height_max`, `tree_height_mean`, `tree_height_gini` | Apex heights of the accepted trees |
| `tree_crownarea_mean`, `tree_crownarea_max`, `tree_crownarea_gini` | Crown area: occupied XY voxel columns × voxel size² |
| `tree_volume_mean`, `tree_volume_max`, `tree_volume_gini` | Voxel volume: occupied voxels × voxel size³ |

A segment is a tree when it has at least `--minimum-tree-voxels` voxels, its
apex is above `--apex-minimum-height`, its smallest PCA extent reaches
`--minimum-tree-thickness` and it occupies at least `--minimum-occupied-layers`
one-metre layers. The segment diagnostics list each check (`pass_vox`,
`pass_apex`, `pass_thick`, `pass_zlayer`) with the underlying values.

## Releases

Images are published to `ghcr.io/3dtrees-earth/3dtrees_raster` by
[`publish-container.yml`](.github/workflows/publish-container.yml) after the
container test suite passes:

| Trigger | Tags |
| --- | --- |
| Push to `main` | `latest`, `sha-<commit>` |
| Tag `v<X.Y.Z>` (for example a GitHub release) | `<X.Y.Z>`, `sha-<commit>` |
| Manual run with `release_tag` | that tag, `sha-<commit>` |

To release, create a GitHub release with a new tag `v<X.Y.Z>` on `main`. Do
not reuse a version for different code: run provenance records the image
version, and users identify results by it.

| Version | Notes |
| --- | --- |
| 1.0.0 | First release: tile metrics, segment diagnostics with apex position, crown area and volume, DTM, CHM and DSM (COG), voxel counts (CSV, multiband COG, COPC), optional kNN distance noise filter, PredInstance_RCT, adjustable scientific parameters. A point cloud without CRS gives products without CRS; a CRS that is present is carried into every product. Replaces the images tagged 1.0.0 on 21 Oct (`sha-3b68a08`) and 2 Nov 2026 (`sha-f3d76a1`) |

## Build from source

```bash
make build IMAGE=3dtrees-raster:local
```

## Prepare a canonical COPC

A COPC stores points in spatial order, but the ground filter and some ties
depend on the original record order. So the tool needs a COPC that carries the
source position of every point as `OriginalPointIndex`, as SmartTile supplies.
To build one from an ordered LAS/LAZ:

```bash
RASTER_IMAGE=3dtrees-raster:local \
RASTER_TEST_CPUS=10 \
RASTER_TEST_MEMORY_GIB=30 \
bash scripts/build_original_order_copc.sh \
  /data/original.laz \
  /data/original.copc.laz \
  /data/original.copc.order.json
```

The builder keeps all point values, checks every source index and publishes the
COPC only after that check passes.

## How it works

`src/raster.R` is the entry point. It validates the inputs and resources, runs
the stages in `src/analysis.R`, checks the complete set of outputs and only
then moves them into the output folder. One R process runs the stages in order.
Spatial windows go to forked R workers with one lidR thread each; streaming
passes run as Python or GDAL child processes.

| Stage | What it reads | Code |
| --- | --- | --- |
| Preflight | headers, input hash | `raster.R`, `production_contract.R` |
| AOI and tile grid | header bounds or the AOI | `aoi_conversion.R`, `analysis.R` |
| Noise filter (optional) | every point at once | `knn_denoising.py`, `denoising.R` |
| DTM | 60 m windows with a 20 m buffer | `analysis.R` |
| Voxel counts | every point, streamed | `voxel_counts.py`, `voxel_counts.R` |
| Point reduction (dense clouds) | every point, streamed | `reduce_points.py` |
| CHM and DSM | windows | `analysis.R` |
| Segments | windows, then one summary per instance | `analysis.R`, `instance_read_groups.R` |
| Tile metrics | one tile per task | `analysis.R`, `tile_scheduler.R` |
| Outputs | tables and rasters | `raster.R`, `artifact_publication.R` |

Every window restores the source record order from `OriginalPointIndex` before
any order-sensitive step. The DTM is a progressive TIN densification (PTD)
ground classification and a TIN per window, merged with GDAL. If lidR's integer
TIN rejects coordinates off the LAS scale lattice, the same window is retried
with lidR's floating-point TIN, without changing coordinates.

### Very dense point clouds

Above 10,000 points per square metre, the CHM, DSM, segments and tiles read a
reduced cloud made in one streaming pass (`src/reduce_points.py`). These
products only use occupied voxels, each instance's highest points and each voxel
column's highest absolute point, so the pass keeps the original points that fix
those and nothing else. The results are the same as from the whole cloud:
checked on 2145 (59 M points on 30 × 31 m), all outputs identical. Voxel counts
still read every point. `RASTER_REDUCTION_DENSITY` sets the threshold.

A DTM window above 300 M points would crash the ground filter, so such a cloud
first keeps its lowest point per 5 cm cell for the DTM, as `lidR::lowest()` does.
Forced onto 2145, the DTM stayed identical. `RASTER_DTM_THINNING_POINTS`
sets the limit. All other products still use the original points.

The reduced cloud is a COPC with the same `OriginalPointIndex`, built with
untwine, so windows read it through the spatial index and memory is estimated
per window. Dataset 707 (786 M points on 99 × 84 m, 10 cores) takes 12:59 min and
27 GB this way, against 1:56 h and 88 GB with thinning inside lidR.

### Point-noise filter

With `--denoise knn` the production runner removes noise points before terrain
estimation with the filter Julian Frey proposed: each point's mean distance to
its k = 8 nearest neighbours, and rejection when that distance exceeds the mean
plus two standard deviations of all points' distances. As in `lidR::sor()`, k
counts the point itself, so k = 8 averages the distances to 7 other points; the
decisions equal lidR's `sor(k = 8, m = 2)` on the whole point cloud. The filter is off by default (`--denoise none`). Rejected points are
named by `OriginalPointIndex`. A canonical COPC carries it; a plain LAS/LAZ that
lacks it is first streamed into a compact, hash-verified scratch LAZ that
numbers the records in file order (0-based), so the noise list names each
record's position in the input file. Input coordinates, classifications and
instance IDs are unchanged.

| Option | Default | Meaning |
| --- | --- | --- |
| `--denoise-neighbors` | `8` | Neighbours for each point's mean distance, the point itself included |
| `--denoise-sigma` | `2` | Standard deviations above the whole cloud's mean distance |
| `--denoise-window-size` | `60` | Width in metres of the windows the noise list is stored in |

`knn_denoising.py` reads all points once, builds one KD-tree (scipy `cKDTree`)
and queries it with all job threads. On dataset 576 (ULS, 53 M points, COPC, 10
threads) the stage takes 81 s instead of 16 min with windowed lidR, with the same
2,114,341 rejected points and byte-identical products. It needs about 96 bytes per point and fails with
`memory_budget_exceeded` when the cloud does not fit the job's memory. The
rejected ordinals are stored per fixed window anchored at `(0,0)`, so later reads
load only the windows they overlap. The mask is computed once and reused for DTM, CHM,
segment and tile metrics. The run publishes `<dataset>_noise_points.csv`, one
rejected `OriginalPointIndex` per row, and `<dataset>_denoising.json`, with the
parameters, the distance mean, standard deviation and threshold, and per-window
point accounting. Existing LAS class labels are preserved and are not rejection
criteria. No filtered point cloud is written; the source is never modified.

The mandatory voxel-count products still count **all original eligible records**,
including points flagged as noise, using the DTM estimated from accepted points.
They are raw recorded-point counts normalized against a denoised terrain model,
not counts of accepted points. A high rejected point can still increase the
voxel GeoTIFF's band count or trigger its explicit capacity limit.

A rule relative to the cloud's own spread always flags its sparsest few per
cent, including some real sparse crown, stem and far-range ground points; on
two noisy test clouds it removed 1.4 % and 4.4 % of the points. Validate it per
acquisition before relying on filtered metrics.

### Height domain

CHM, segment and tile processing keep all accepted points with a finite,
non-negative height above the DTM; there is no upper height cutoff. Tile
densities and the vertical complexity use one shared height domain,
`H = max(1, ceiling(maximum accepted voxel height))` metres over the whole cloud,
even when an AOI restricts the reported tiles. It is recorded as
`metric_height_upper_m` in each tile row, the performance CSV and the run
provenance, with the exact maximum as `observed_maximum_voxel_height_m_exact`.

### Resources

CPUs are the smallest of `GALAXY_SLOTS`, `RASTER_THREADS` and the
container quota, at most 10. Memory is the smaller of `GALAXY_MEMORY_MB` and the
container limit, with 15 % headroom; `--memory-budget-gib` can lower it. The
tool uses what the job gets. The memory request by input size belongs to the
Galaxy dispatcher (TPV).

Before each wave of windows the runner counts the points in every window from
the COPC hierarchy, without decoding them, and starts as many workers as fit
the budget. It reserves 128 MiB per task plus, per counted point, 160 bytes for
the DTM, 192 for the CHM, 256 for segments and 128 for tiles. These come from
real memory measured with one worker, plus a margin. A task that cannot fit
stops the run before it reads points, with `memory_budget_exceeded`. A plain
LAS/LAZ has no hierarchy, so every task reserves the whole file.
`RASTER_CATALOG_WORKERS` (1–10) sets the worker count explicitly.

### Galaxy

The Galaxy wrapper lives in
[bgruening/galaxytools](https://github.com/bgruening/galaxytools/tree/master/tools/3dtrees_raster)
and pins `ghcr.io/3dtrees-earth/3dtrees_raster:<version>`. Flag arguments go
to the runner; any other command, such as Galaxy's `/bin/bash tool_script.sh`,
runs as given. Stage the input as a hard link or copy with its `.laz`
extension: a symlink to an extension-less Galaxy dataset is rejected by rlas.

### Voxel counts

Voxel counts read every record, including points the noise filter flags, and
never use the reduced cloud. They decode only X, Y and Z and stream the point
cloud in batches, so memory does not grow with the cloud. A COPC's chunks are
split across one worker process per job thread; each counts its share and the
counts are summed, which gives the same result as one process. A plain LAS/LAZ
is read by one process with that many decoding threads. Each worker reserves
192 MiB, and fewer workers start when the memory budget is small. The
multiband GeoTIFF, the CSV and the COPC hold the same voxels and counts; the
metadata JSON records the grid, the CRS and the point accounting.

## Tests

| Target | What it checks |
| --- | --- |
| `make test` | Unit tests in R and Python, the image entry point, then `make test-raster` |
| `make test-raster` | A synthetic COPC end to end: products, layouts, provenance |
| `make test-voxel-counts` | Voxel counts, the AOI without occupied voxels, COPC workers against one process, and the point reduction |
| `make test-denoising` | The kNN filter against `lidR::sor()` and the uncapped heights |
| `make test-production-contract` | COPC and LAS/LAZ layouts, AOI and full extent give the same results |
| `make test-copc-all-instance-dimensions` | All instance dimensions from one COPC read |
| `make test-projected-3d-crs` | Projected 3D and compound CRSs |
| `make test-parallel-production` | One worker against 2 and 4 workers give identical results |

Set `IMAGE` to test a built image, for example `make test-raster IMAGE=3dtrees-raster:local`.

## Repository layout

| Path | Contents |
| --- | --- |
| `src/raster.R` | Entry point: arguments, validation, provenance, publication |
| `src/analysis.R` | The stages: grid, DTM, CHM, DSM, segments, tile metrics |
| `src/*.R`, `src/*.py` | Helpers named after their task |
| `tests/` | Tests, fixtures and the canonical COPC builder |
| `Dockerfile` | The published image |
