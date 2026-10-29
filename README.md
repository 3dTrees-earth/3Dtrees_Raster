# 3Dtrees: Raster

3Dtrees: Raster turns one segmented forest point cloud into raster and tile
products of forest structure:

- a digital terrain model (DTM) and a canopy height model (CHM),
- full-resolution voxel point counts as a multiband GeoTIFF and a CSV,
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

> **Naming.** The project was developed as *ForestStructure*. For compatibility
> the image keeps that name internally: the runner lives in
> `/opt/foreststructure`, operator settings use `FORESTSTRUCTURE_*` variables,
> and output files start with `foreststructure_`.

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

- The point cloud needs a projected metric CRS and at least one of the instance
  dimensions `PredInstance`, `PredInstance_SAT` (SegmentAnyTree) or
  `PredInstance_FM` (ForestMamba). Name one or more with
  `--instance-dimension` (repeat the option or separate names with commas);
  a missing dimension fails the run before any analysis.
- Point-noise filtering (`--denoise sor`, the default) works on a canonical
  COPC and on plain LAS/LAZ. A plain file without `OriginalPointIndex` is
  numbered in its own record order, which costs one extra read and a scratch
  copy of XYZ and the instance dimensions. Add `--denoise none` to skip it.
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
| `--denoise` | `sor` | `sor` runs statistical outlier removal before the terrain model; `none` disables it |
| `--denoise-neighbors` | `10` | Nearest neighbours used for each point's mean distance |
| `--denoise-sigma` | `12` | Points farther than the window mean plus this many standard deviations are flagged |
| `--denoise-window-size` | `60` | Width in metres of the fixed windows that define the statistics |
| `--denoise-buffer` | `20` | Context in metres around each window included in its statistics |

See [Point-noise classification and uncapped heights](#point-noise-classification-and-uncapped-heights)
for the method and its limits.

### Run options

| Option | Meaning |
| --- | --- |
| `--point-cloud` | Input LAS/LAZ or canonical COPC, recognised from its content (required) |
| `--instance-dimension` | Instance dimensions to process (required) |
| `--aoi-geojson` | Optional area of interest in the point cloud's XY coordinates; not reprojected |
| `--original-point-cloud` | Optional original LAS/LAZ of a COPC, used only for header and provenance checks |
| `--output-dir`, `--temp-dir` | Output and scratch directories |
| `--raster-format` | `COG` (default) or `GTiff` for the DTM, CHM and voxel-count rasters |
| `--version` | Print the tool version and exit |
| `--memory-budget-gib` | Optional lower internal memory budget; the container or Galaxy limit is used otherwise |
| `--target-chunk-points` | Approximate points per COPC CHM/segment read window (default 50,000,000; `0` uses fixed windows) |
| `--voxel-count-chunk-points` | Maximum records per voxel-count batch (default 2,000,000) |
| `--sensor`, `--country` | Optional labels copied into the result rows (default `NA`) |

CPU and memory come from `GALAXY_SLOTS`/`GALAXY_MEMORY_MB` or the container's
limits; see [Galaxy production contract](#galaxy-production-contract) for the
operator variables. `--dataset-id` and `--export-voxel-counts` are still accepted
for old commands but are not needed.

## Outputs

Per requested instance dimension `<dimension>`, one file in each of three
folders:

| File | Contents |
| --- | --- |
| `tile_metrics/<dimension>_tile_metrics.csv` | One row of [tile metrics](#tile-metrics) per tile |
| `segment_diagnostics/<dimension>_segment_diagnostics.csv` | One row per segment with the tree-acceptance checks, the apex position (`apex_x`, `apex_y`, source CRS) and the segment's crown area (`crown_area_m2`, occupied voxel columns × voxel area) and voxel volume (`voxel_volume_m3`, occupied voxels × voxel volume), the per-segment values the tile `tree_*` aggregates are computed from |
| `tiles/<dimension>_tiles.geojson` | The tile metrics as polygons |

Shared by all dimensions, each starting with the prefix `foreststructure`:

| File | Contents |
| --- | --- |
| `foreststructure_dtm.tif` | DTM GeoTIFF (Cloud Optimized by default) |
| `foreststructure_chm.tif` | CHM GeoTIFF (Cloud Optimized by default) |
| `foreststructure_tiles.png` | Tile overview map |
| `foreststructure_voxel_counts.csv` | Occupied voxels: `x_center,y_center,height_lower_m,height_upper_m,point_count` |
| `foreststructure_voxel_counts.tif` | UInt32 point counts; with voxel size *s*, band *b* holds heights `[(b-1)·s, b·s)` m above the DTM |
| `foreststructure_voxel_counts_validity.tif` | 0 excluded (outside the AOI or on its boundary), 1 valid, 2 invalid DTM |
| `foreststructure_voxel_counts.json` | Grid, CRS, source and DTM hashes, policies and point accounting |
| `foreststructure_noise_points.csv` | With SOR only: one rejected `OriginalPointIndex` per row (the 0-based record position for a plain LAS/LAZ without that field) |
| `foreststructure_denoising.json` | With SOR only: parameters and per-window point accounting |
| `foreststructure_julia_memory_safe_run.json` | Run provenance: inputs, parameters, versions, resources |
| `foreststructure_performance.csv` | Stage timings and memory observations |
| `foreststructure_aoi_conversion.json` | The analysed footprint and its AOI conversion checks |

Voxel counts include every recorded point, also those flagged as noise. Zero
means no recorded points, not absence of vegetation; read the validity raster
before interpreting missing CSV rows. See
[Full-resolution voxel point-count export](#full-resolution-voxel-point-count-export).

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
| 1.0.0 | First release: tile metrics, DTM, CHM, voxel counts, point-noise filtering and adjustable scientific parameters |

## Build from source

```bash
make build IMAGE=3dtrees_raster:local
```

`Dockerfile` is the published image. `Dockerfile.julia-memory-safe` builds
the ordered-COPC development image used by several test targets:

```bash
make build-julia-memory-safe JULIA_MEMORY_SAFE_IMAGE=3dtrees_raster:copc-local
```

## Prepare a canonical COPC

Do not use an arbitrary COPC. Build it from the ordered LAS/LAZ so it contains
the validated `OriginalPointIndex` needed to restore source record order:

```bash
FORESTSTRUCTURE_JULIA_IMAGE=3dtrees_raster:copc-local \
FORESTSTRUCTURE_TEST_CPUS=10 \
FORESTSTRUCTURE_TEST_MEMORY_GIB=30 \
bash tests/build_original_order_copc.sh \
  /data/original.laz \
  /data/original.copc.laz \
  /data/original.copc.order.json
```

The builder preserves all point values, places the supported instance
dimensions in selectively readable positions, validates every source index,
and publishes the COPC only after validation succeeds.

## Technical reference

The rest of this document explains the implementation, its validation and its
limits in detail.

### Method notes

The production workflow uses a canonical COPC for spatial streaming and source
record-order preservation. Native lidR denoising is enabled by default, and
normalized heights have no upper cutoff; these change the historical metric policy. See
[`COPC_REVIEW_GUIDE.md`](COPC_REVIEW_GUIDE.md) for the design, failure fixes,
and validation evidence.

The preserved historical scientific reference is
[`reference/Indices_Final_run.R`](reference/Indices_Final_run.R), SHA-256
`746d57b4c937001af31e4ccd1b9f14edb5cebb15d46154ad9e20d0ce39f78226`.

DTM generation first uses lidR's optimized integer TIN. If that implementation
rejects coordinates that do not match the LAS scale/offset integer grid, the
same chunk is retried with lidR 4.3.2's legacy floating-point TIN. The fallback
does not rescale coordinates or modify the LAS header.

A real replay of dataset 2153 (`8101.laz`, 40,453,840 points) reproduced the
historical integer-conversion failure with the old image and completed with the
legacy retry under 10 CPUs and a 30 GiB limit. The accepted run used 6.25 GiB
peak cgroup memory and produced all 12 expected artifacts. This proves failure
recovery, but not zero-tolerance `valid_updated` parity: dataset 2153 has no
confirmed oracle in the validation cohort, and comparison with its older
generic-COPC result found small raster and segment differences. These small
changes are expected when the fallback triangulates the original floating-point
XY coordinates instead of temporarily changing their scale for the integer TIN.
Keeping the source coordinates and LAS header unchanged is scientifically more
faithful to the measured point cloud. See
[`COPC_REVIEW_GUIDE.md`](COPC_REVIEW_GUIDE.md) and the machine-readable
validation record in `benchmarks/results/` before promoting this retry strategy.

### Resource-limited run

One conservative validated profile is 10 CPUs, a 30 GiB hard memory limit, a
25 GiB internal budget, and one LAScatalog worker:

```bash
docker run --rm --network none \
  --cpus 10 \
  --memory 30g \
  --memory-swap 30g \
  --user "$(id -u):$(id -g)" \
  --env FORESTSTRUCTURE_THREADS=10 \
  --env FORESTSTRUCTURE_CATALOG_WORKERS=1 \
  --volume /local/input:/in:ro \
  --volume /local/output:/out \
  --volume /local/work:/work \
  ghcr.io/3dtrees-earth/3dtrees_raster:1.0.0 \
  --point-cloud /in/original.copc.laz \
  --original-point-cloud /in/original.laz \
  --aoi-geojson /in/aoi.geojson \
  --instance-dimension PredInstance_SAT \
  --output-dir /out \
  --temp-dir /work \
  --memory-budget-gib 25 \
  --target-chunk-points 50000000
```

`--original-point-cloud` is optional. When supplied, it is used only for
header compatibility and provenance checks; scientific point reads still use COPC.
Extrema are compared on each header's coordinate lattice because original writers
can retain pre-quantization bounds. This lightweight check is not identity proof: the
preparation validator must still verify complete point indices and fingerprints.
Plain LAS/LAZ input is supported with and without denoising.

### Galaxy production contract

In the 3Dtrees pipeline Galaxy supplies one full-resolution COPC and optionally
an Audit AOI in that COPC's metric XY coordinate system. The COPC must retain
`OriginalPointIndex` and every supported instance dimension the job requests;
SmartTile supplies global IDs. A plain LAS/LAZ also runs, with point-noise
filtering disabled.
An AOI restricts the analysis domain. Without it, the runner reuses the existing
coordinate-lattice-snapped full-extent footprint. Neither mode subsamples points.
The existing DTM algorithm runs on SOR-accepted points. CHM, segment and tile
metrics use the same accepted-point mask and retain finite heights above 70 m.

Raster georeferencing uses the horizontal two-dimensional CRS, embedded in each
GeoTIFF. Projected 3D and compound source CRSs are retained in provenance;
coordinates and elevations are not transformed. Voxel band heights remain
metres above the DTM. Publication does not depend on GDAL PAM sidecars.

`GALAXY_SLOTS` supplies the CPU allowance; `GALAXY_MEMORY_MB` supplies the job RAM
allocation in MiB. The runner uses the smaller declared/container memory limit
and reserves 15% headroom, with a hard ceiling of 100 decimal GB (93.132 GiB).
An explicit internal budget may lower it. There is no 60 GiB internal cap when
Galaxy or the container supplies an allocation; 60 GiB remains the fallback
when neither supplies a limit. Container CPU quotas also cap threads, with an
overall maximum of 10 CPUs. `FORESTSTRUCTURE_THREADS` remains an operator cap;
without a declared CPU allowance or quota, the safe default is one thread.
Catalog concurrency defaults to the smaller of ten workers and the CPU allowance.
Memory admission may reduce it further for each wave. The job destination can set
`FORESTSTRUCTURE_CATALOG_WORKERS` explicitly; it is not a user-facing parameter.

The output directory defaults to the Galaxy job working directory and scratch
to `TMPDIR` (or R's private temp directory). A wrapper may supply `--output-dir`
and `--temp-dir` explicitly. Stage a Galaxy dataset as a hard link or copy named
`input.copc.laz` because this runner identifies COPC mode from that extension.
A symlink does not work: rlas resolves it and then rejects the extension-less
Galaxy dataset file (`File not supported`); the same applies to the AOI. Never
write beside Galaxy's original input. No dataset ID, CPU, RAM, scratch path or voxel
export toggle belongs in the Galaxy user form. Galaxy must enforce its job CPU,
RAM and scratch allocations; internal estimates are not hard OS limits.

```bash
# Within a job with GALAXY_SLOTS, GALAXY_MEMORY_MB and TMPDIR assigned by Galaxy:
Rscript /opt/foreststructure/run_julia_memory_safe.R \
  --point-cloud input.copc.laz --instance-dimension PredInstance_SAT
# Optionally add --aoi-geojson aoi.geojson.
```

The container entry point already invokes that script. Both Dockerfiles now use the production runner as their default entry point.
Arguments that start with a flag go to the runner; any other command, such as
Galaxy's `docker run IMAGE /bin/bash tool_script.sh`, runs as given.
[`galaxy/3dtrees_foreststructure`](galaxy/3dtrees_foreststructure) is the Galaxy
wrapper. It detects a COPC from its info VLR, stages the inputs as above and
pins `ghcr.io/3dtrees-earth/3dtrees_raster:<tool version>`, so that
image tag must be published before the wrapper is installed. Run its tests with
`planemo test --docker galaxy/3dtrees_foreststructure/foreststructure.xml`.
The production container and Galaxy wrapper still need to be released/installed
together; these source changes do not deploy a server tool. Validate integration with
`make test-production-contract JULIA_MEMORY_SAFE_IMAGE=<validated-runtime-image>`.
This compares input layouts and AOI/full-extent behavior using identical current
sources with denoising explicitly disabled, including 15th-field IDs, a late
order key and mandatory voxel products. `make test-denoising` exercises the default
SOR path and uncapped metric policy.

### Galaxy allocation from merged COPC size

`src/galaxy_resource_request.py` provides the dispatcher-side
`request_resources(input_bytes)` policy. Supply the recorded size of the single
merged, labelled COPC that this job receives. Do not sum original uploads or the
optional AOI, and do not apply a decompression multiplier. The function opens no
files and performs no COPC hierarchy traversal.

| Merged COPC size (decimal GB) | CPUs | RAM target (decimal GB) | Galaxy allocation (MiB) |
| --- | ---: | ---: | ---: |
| Up to 1 | 8 | 30 | 28,610 |
| Over 1, up to 5 | 10 | 30 | 28,610 |
| Over 5, up to 10 | 10 | 60 | 57,220 |
| Over 10 | 10 | 100 | 95,367 |

The MiB conversion rounds down to respect the 100 GB ceiling. These are initial
allocation tiers, not guarantees of peak memory: compression, point density and
spatial extent vary. Inside the job, the existing COPC hierarchy planner checks
each actual spatial task and reduces concurrency to fit the available memory.
A task that cannot fit fails explicitly. The tiers have not been benchmarked
as a complete policy; the parallel benchmark uses an explicit 30 **GiB** job
limit, which differs from the 30 decimal **GB** tier.

```bash
python3 src/galaxy_resource_request.py --input-bytes 13332084830
make test-galaxy-resource-request
```

[`galaxy/tpv.example.yml`](galaxy/tpv.example.yml) is an integration fragment for
Galaxy's Total Perspective Vortex dispatcher. Install the policy module on the
dispatcher's Python path and merge the rule into the site's existing TPV
configuration. Its tool pattern and `point_cloud` parameter match the wrapper. The fragment uses the parameter-specific
`helpers.get_input_size(..., estimate_uncompressed_size=False)` interface
documented by [TPV](https://total-perspective-vortex.readthedocs.io/en/latest/topics/advanced_topics.html#input-sizes).
TPV reports input size and RAM in GiB; the example converts these units explicitly.
The destination must enforce the final request and supply its actual CPU and MiB
allocation to the runtime, including any site-level reduction. The example does
not configure a runner or deploy a tool; Galaxy supports these rules through its
[dynamic job dispatcher](https://docs.galaxyproject.org/en/latest/admin/jobs.html#dynamic-destination-mapping).

### Catalog and tile concurrency

`FORESTSTRUCTURE_CATALOG_WORKERS` accepts 1–10 workers within
`FORESTSTRUCTURE_THREADS`; production automatically caps at ten when unset. The standard
CLI exposes the same upper limit through `--catalog-workers`. Tile metrics
reuse that process count. Forked workers use one native thread each to avoid
inherited OpenMP pool deadlocks; a sequential worker can use the full allowance.
Increase workers only after measuring the largest spatial reads: compressed
file size alone does not predict the memory of locally dense chunks.

On the pinned lidR 4.3.2 Linux runtime, bounded fork workers call lidR's existing
chunk engine directly. This avoids repeated future-global discovery while
retaining chunk construction, buffers, callbacks, raster writers and ordered
results. Each child handles one task and exits, releasing its working set.
Memory-admitted production runs require this pinned runtime and safe fork support;
the standard CLI retains its public-scheduler fallback. BLAS pools are limited to one thread to avoid
multiplying native threads inside workers. Scientific algorithms and source
point order are unchanged.

Before dispatch, production reads COPC hierarchy metadata to bound points in each
actual spatial read, including DTM buffers. It reserves 256 bytes per estimated
point for DTM and 512 bytes for CHM, segment and tile work, plus 128 MiB per task,
retaining another 512 MiB inside the internal budget.
Reservations across concurrent tasks must fit the available working set. Only
demonstrably clean inactive file cache is treated as reclaimable. Adaptive CHM
reserves its largest sequential child read; raster footprints remain unchanged.
An individual task that cannot fit fails before loading points. These coefficients
are conservative planning estimates, not universal allocation guarantees: keep equal hard
`--memory`/`--memory-swap` limits. A 30 GiB allocation allows a 25.5 GiB internal
budget. Provenance records requested workers, actual waves, point bounds, phase
coefficients, fixed allowances, memory snapshots and the raw container peak.
Plain LAZ uses a conservative whole-file
count fallback. COPC estimation fails explicitly above 64 MiB of hierarchy metadata
or 100,000 read windows per phase.

DTM releases the source-order column before classification and triangulation and
does not retain voxel or segment tables. Its separate coefficient also preserves
admission of the previously validated dense 1805 terrain geometry: the largest
buffered window has a metadata upper bound of 185,339,618 points, reserving
44.31 GiB at 256 bytes per point. A uniform 512-byte estimate would reserve
88.50 GiB and reject it even under its former 70 GiB allocation, although the
historical complete run peaked at 27.38 GiB. The revised estimate fits the
59.5 GiB internal budget of a 70 GiB job; a 30 GiB job still rejects that oversized
window. The historical peak supports this initial estimate for the pinned stack;
it does not establish a bound for arbitrary point distributions or native
allocators. The historical original-upload cloud lacks production instance labels.

SAT/FM fields on the same source share one read, order restoration and terrain
normalization. Their disk summaries and reductions remain independent. Different
legacy projection files retain separate reads. Terrain, canopy, tile metrics and
raw voxel counts still have distinct passes; their scientific semantics differ.
Run `make test-parallel-production` for exact current-policy parity at 1, 2 and 4 workers,
with both full extent and an AOI.

### Measured parallel qualification

The labelled COPCs contain 34.19 million points / 1.094 decimal GB (2551) and
147.96 million points / 3.910 decimal GB (1998). All candidate runs used 10 CPUs,
a 30 GiB hard limit, a 25.5 GiB internal budget and no additional swap, GPU or
network. Inputs and Audit AOIs match their saved baselines by SHA-256.

| Dataset | Worker cap | Processing minutes | Peak GiB | Average CPU cores |
| --- | ---: | ---: | ---: | ---: |
| 2551 | 1 | 10.18 | 5.60 | 1.07 |
| 2551 | 4 | 6.75 | 6.97 | 1.36 |
| 1998 | 1 | 67.90 | 3.67 | 1.13 |
| 1998 | 4 | 30.82 | 9.24 | 2.26 |
| 1998 | 10 | 18.38 | 11.41 | 3.74 |

Every scientific CSV, complete tile vector, raster cell and validity mask matches
its baseline exactly. SAT/FM callback point records fell from 295,922,732 to
147,961,366 for 1998. Shared summaries increase its sampled scratch peak from
1.93 GiB to 3.59 GiB. The ten-worker profile reduces processing time by about 40%
from the four-worker profile, while keeping its observed peak below 12 GiB.

These are single observations on a shared host. The baseline had 70 GiB allocated;
its 1998 run overlapped 1996. The four-worker candidate used the initial 512-byte
DTM coefficient, while the ten-worker candidate uses 256 bytes. This comparison
therefore measures the combined profile changes, not isolated worker scaling.
Times include analysis-container launch/teardown, but exclude preflight and the
external acceptance comparisons; pipeline-stage timings are recorded separately.
The smaller dataset has four terrain/canopy tasks and only one analysis tile,
so it offers less parallelism.

The [complete benchmark evidence](benchmarks/results/2026-10-13-parallel-memory.json)
records source/runtime identifiers, phase timings, CPU and kernel memory counters,
worker distributions, scratch observations, exact comparisons and test coverage.
The final synthetic full-extent/AOI tests also exercise automatic ten-worker
selection and high-ordinal instance fields. Galaxy allocation tiers remain initial
policy: the real benchmark's 30 GiB limit differs from the 30 decimal GB tier.

The full-resolution study benchmarks use ten workers for datasets 647/2080 and
three for dense dataset 1805, all under 10 CPUs, a 70 GiB hard limit and a
60 GiB internal budget. See the
[measured profiles and exact-parity evidence](benchmarks/results/2026-10-12-voxel-optimization.json).
These are measured
dataset profiles rather than a fixed RAM-per-compressed-GB rule.

| Dataset | Original GB | Million points | Baseline minutes | Optimized minutes | Peak GiB |
|---|---:|---:|---:|---:|---:|
| 647 | 0.445 | 44.99 | 19.72 | 3.20 | 3.77 |
| 2080 | 0.645 | 60.99 | 27.62 | 3.18 | 5.57 |
| 1805 | 7.750 | 1425.44 | 155.12 | 54.27 | 27.38 |

The baseline allowed 20 CPUs with one catalog worker. These single-run analysis
timings exclude downloads and COPC preparation; cache state and shared-host load
can affect measurements. Peaks are observations, not guaranteed minimum limits.
Every scientific raster cell, CSV count and tile vector matches the saved
baseline. Julia's separate reference dataset also retains exact metric parity.

### Adaptive COPC windows

COPC CHM and segmentation reads target approximately 50 million points per
spatial window by default. The runner reads the point count and XY extent from
the COPC header, calculates mean density, and selects:

```text
floor(sqrt(50,000,000 / mean_density)) metres
```

The result is aligned downward to a whole metre and bounded to 1–60 m. A
smaller `--chunk-size` remains an additional ceiling; values above 60 m cannot
widen the validated maximum. Use `--target-chunk-points 0` only for controlled
comparisons that need fixed metre-based behavior. Plain LAS/LAZ always retains
fixed-size behavior.

The 50-million-point value is an approximate planning target based on mean
density, not a hard per-window limit. Locally denser windows may exceed it. The
run warns without failing and records the target, density, selected size,
estimated count, observed minimum/maximum counts, and oversize count in the
performance CSV and run-provenance JSON.

CHM and segmentation use zero buffer and an absolute one-metre chunk alignment.
The runner therefore only accepts voxel and CHM resolutions that divide 1 m
evenly, so no cell straddles a read edge (the defaults are 0.2 m and 0.5 m). CHM
reads adapt, but its raster footprints retain the configured baseline size
(normally 60 m). Each read uses the unchanged normalization and voxelization;
only the maximum voxel height at each exact XY is retained across reads
(about 90,000 XY columns per 60 m parent at the default 0.2 m resolution). Raw points and
3D voxels are released before the next read. The parent is then rasterized
once with the original numeric-resolution `p2r()` and written to disk for the
VRT merge. This preserves lidR's occupied-extent edge assignment: rasterizing
each smaller read independently can move an exact-edge point by one cell.
Segmentation writes compact voxel, layer, and apex summaries to the mounted SSD
before releasing raw points and temporary tables; the existing bucket reducer
then aggregates one bounded partition at a time.

`OriginalPointIndex` makes a changed window layout scientifically deterministic:
equal-height apexes choose the lowest original index, duplicate voxels retain
their earliest original index, and voxels return to global source order before
PCA. Output artifact names and schemas are unchanged.

DTM is deliberately outside this policy. It still uses 60 m chunks and its
existing 20 m buffer because terrain interpolation needs neighboring ground
points and was not the memory-failure stage addressed by adaptive chunking.
DTM resolutions must therefore divide 60 m evenly.
For COPC, DTM still restores Julia's source-point order first; it then releases
the now-redundant `OriginalPointIndex` column before PTD/TIN. Its checks for
three non-collinear input and ground points scan fixed-size blocks instead of
materializing duplicate full-chunk XY tables. These are memory-lifetime changes
only: DTM geometry, buffer, classification, triangulation, and output remain
unchanged.

Supported instance dimensions are:

- `PredInstance`
- `PredInstance_SAT`
- `PredInstance_FM`

The production runner requires `--instance-dimension` (repeat it or separate
names with commas) and never discovers dimensions from the point cloud. The
selection is strict: any missing or unsupported field fails before analysis. SmartTile owns
global instance-ID assignment; this tool never renumbers or fabricates IDs.
Available dimensions are processed independently so their results cannot mix.

Instance fields may occur anywhere in the ExtraBytes schema, including position
15. The pinned rlas COPC reader can only load the first nine extras, even with
its all-fields selector. When a required field lies beyond nine, the production
runner reuses the bounded, hash-verified projection helper once to create a
compact scratch LAZ with XYZ, OriginalPointIndex and the selected instance fields.
It builds a LAX spatial index and uses the same adaptive spatial windows and
source-order restoration. This costs an additional full-source read, verification
read and scratch storage. Normal low-ordinal COPCs retain their direct fast path.
The mandatory voxel-count pass always reads the full original COPC, including
all duplicates, rather than a reduced or subsampled source.

### Point-noise classification and uncapped heights

The production runner defaults to `--denoise sor`, using native lidR 4.3.2
statistical outlier removal before terrain estimation. Use `--denoise none` to
disable it. Rejected points are named by `OriginalPointIndex`. A canonical COPC
carries it; a plain LAS/LAZ that lacks it is first streamed into a compact,
hash-verified scratch LAZ that numbers the records in file order (0-based),
so the noise list names each record's position in the input file. Late
ExtraByte positions use the existing projection path.
Input coordinates, classifications and instance IDs are unchanged.

SOR measures each point's mean distance to its nearest neighbours and rejects
unusually large distances relative to the local window's mean and standard
deviation. It does not impose an absolute height ceiling or trim a height
percentile. Parameters are scientific choices and are recorded in provenance:

| Option | Default | Meaning |
| --- | --- | --- |
| `--denoise-neighbors` | `10` | Neighbours used for the point-distance statistic |
| `--denoise-sigma` | `12` | Standard deviations above the window mean |
| `--denoise-window-size` | `60` | Fixed scientific core width in metres |
| `--denoise-buffer` | `20` | Context around each core in metres |

Scientific windows are anchored at `(0,0)` and assign each point to one half-open
core. Their buffered populations define the SOR statistics independently of the
adaptive execution windows and catalog worker count. Sparse windows with no more
points than the neighbour count are retained and reported. A dense window that
cannot fit the memory admission budget fails; it is not silently split into a
different statistical population. The finite context is part of the method:
this is local buffered SOR, not whole-cloud nearest-neighbour classification.

The runner computes the mask once and reuses it for DTM, CHM, segment and tile
metrics. It publishes `<dataset>_noise_points.csv`, with one rejected
`OriginalPointIndex` per row, and `<dataset>_denoising.json`, with parameters,
ownership rules and per-window point accounting. Existing LAS class labels are
preserved but are not additional rejection criteria. No filtered point cloud is
written, and the source file is never modified.

The mandatory voxel-count products still count **all original eligible records**,
including points flagged by SOR, using the DTM estimated from accepted points.
Consequently these are raw recorded-point counts normalized against a denoised
terrain model; they are not counts of accepted points. Negative-height and
invalid-terrain exclusions remain documented separately. A high rejected point
can still increase the voxel GeoTIFF's band count or trigger its explicit
capacity limit.

The conservative 12-sigma default targets extreme spatial outliers. It needs
acquisition-specific validation: sparse genuine crowns may be flagged, and
coherent noise clusters may survive. It is not a universal forest-noise rule.
A synthetic sensitivity check at sigma 6/8/10/12/15 rejected 36/4/2/2/2 points
from the earlier 56,057-point fixture. Sigma 12 retained every ordinary return;
sigma 15 missed a high outlier in the separate 568-point fixture. These checks
justify the chosen starting value, not an accuracy claim for real surveys.
For comparison, lidR's own `sor()` default is k = 10, m = 3, and the PCL tutorial
uses k = 50 with one standard deviation. These are much stricter: in close-range
(TLS, MLS, ULS) scans point spacing grows with distance from the scanner, so they
can flag real upper-crown and stand-edge points. Pass `--denoise-sigma 3` to use
lidR's stricter threshold.

CHM, segment and tile processing retain all accepted points with finite,
nonnegative normalized heights, including heights above 70 m. The vegetation
minimum height (0.5 m by default) remains part of the existing tile indices; the
raw voxel export still includes the layers below it. Point retention has no upper height
cutoff. Metric denominators and vertical-complexity normalization use one shared
height domain: `H = max(1, ceiling(maximum accepted normalized voxel height))` metres over the
full source cloud, even when an AOI restricts reported tiles. The value is recorded
as `metric_height_upper_m` in each tile row, performance CSV and provenance;
`metric_height_policy` is `observed_voxel_height_v1`. The exact maximum is gathered
during existing CHM point reads, before canopy rasterization rounds heights,
and recorded as `observed_maximum_voxel_height_m` in performance and provenance.
The companion `observed_maximum_voxel_height_m_exact` string preserves 17-digit
precision for reconstructing the domain at integer-metre boundaries.
This denominator is not a point filter.
Sparse vertical-bin accounting avoids allocating a dense array proportional to
an extreme retained height. Density and complexity values can differ from older
70 m-domain results even with `--denoise none`; cross-survey comparisons must
account for the recorded height domain.

Use `make test-denoising` for the default CLI, an explicit unfiltered control,
isolated high/low outliers, an 80–90 m supported crown, high-ordinal IDs, raw-count
semantics, the shared height domain and worker/adaptive-window invariance.
The [current synthetic acceptance](benchmarks/results/2026-10-14-denoising-height-policy.json)
retains all 1,331 supported crown points, flags only the two injected extremes
out of 56,057 records, and matches all 14 scientific artifacts exactly across
one/two workers and 60/23 m adaptive read windows. The shared height domain is
86 m after SOR and 250 m in the explicit unfiltered control. The rejected high
point is still present in the raw voxel export.
The [earlier synthetic acceptance](benchmarks/results/2026-10-13-denoising.json)
records the superseded opt-in sigma-6 policy and its 34 additional ordinary
observations flagged. Earlier real-data performance and parity benchmarks also
predate the current denoising and height-policy changes; there is no
production-scale denoising benchmark or acquisition-specific accuracy claim.
Historical Julia-oracle scripts require a pinned runtime and any source override
from their recorded scientific policy. Set `FORESTSTRUCTURE_SKIP_BUILD=1` where
the replay script supports building; do not rebuild that runtime from current
sources. These scripts keep their strict comparators and are not exact-parity
gates for the changed density/index definitions.

Both modes retain canonical source points on a terminal catalog boundary.
Earlier versions could omit such records when the source maximum lay exactly
on a CHM or segment catalog edge, changing results with execution-window size.
The dedicated boundary regression checks this correction; outputs for those
previously affected inputs can change even with denoising disabled.

### AOI handling

The AOI is interpreted in the point cloud's local XY coordinate space and is
not reprojected. The GeoJSON-to-GeoPackage round trip must preserve topology,
exact XY vertices, and exact bounds. An absolute area difference of at most
`1e-6` square coordinate units is allowed solely for floating-point summation
noise. Without an AOI, the full point-cloud extent is tiled.

### Validation tools

Reusable acceptance checks live in `validation/`; `tests/` contains regression
and integration orchestration. Validators return a nonzero exit code on an
invalid result and write machine-readable JSON where a report path is accepted.

Validate that ordered LAZ-to-COPC conversion preserved every indexed point and
supported instance value:

```bash
python3 validation/validate_indexed_copc_streaming.py \
  /data/original.laz \
  /data/original.copc.laz \
  /data/original.copc.order.json
```

Validate the actual GeoJSON-to-GeoPackage round trip with the image's pinned
R, sf, and GDAL stack:

```bash
docker run --rm --network none \
  --entrypoint Rscript \
  --volume "$PWD:/repo:ro" \
  --volume /data/aoi.geojson:/in/aoi.geojson:ro \
  ghcr.io/3dtrees-earth/3dtrees_raster:1.0.0 \
  /repo/validation/validate_aoi_roundtrip.R \
  /in/aoi.geojson
```

Validate a completed output directory before publishing it:

```bash
python3 validation/validate_output_artifacts.py \
  2001 /data/attempts/2001/output \
  --report /data/attempts/2001/output/forest_structure.validation.json
```

The report's `publish_artifacts` array is the exact file allowlist. Copy those
files plus the new validation report; do not copy other files from the attempt
directory.

GDAL may create these non-scientific PAM metadata caches while writing the DTM
and CHM:

- `<dataset>_dtm.tif.aux.xml`
- `<dataset>_dtm.tif.partial.tif.aux.xml`
- `<dataset>_chm.tif.aux.xml`
- `<dataset>_chm.tif.partial.tif.aux.xml`

`validate_output_artifacts.py` reports those exact names under
`ignored_gdal_sidecars` and excludes them from publication. No wildcard
exception is used: every other unexpected file still fails validation. This
allows completed datasets such as 2001 to be reviewed and published without
recomputing their scientific outputs.

For zero-tolerance comparison against `valid_updated`, use
`validation/compare_valid_updated_science.py` for CSV and GeoJSON outputs and
`validation/compare_valid_updated_raster.R` for DTM and CHM rasters.

### Full-resolution voxel point-count export

The production (`run_julia_memory_safe.R`) runner always exports voxel counts;
there is no opt-out. `--export-voxel-counts` remains accepted for old invocations.
The legacy `run.R` entry point retains its opt-in flag. The optional
`--voxel-count-chunk-points` (default 2,000,000) bounds the number of source records
per count batch; the memory-safe runner lowers it further to fit its remaining
memory allowance. Larger batches keep LAZ decompression parallel: on dataset 226
(451M points, SSD) 2M-point batches cut the pass from 323 s to 222 s for about
1.3 GiB more memory; 8M was no faster.
Existing index/CHM filters and output schemas are unchanged.

The export reuses the generated global DTM and usable Audit AOI. Supply the
original full-resolution LAS/LAZ or a lossless COPC with verified conversion
provenance. Every record is read, including repeated coordinates; no point
sampling or COPC level-of-detail restriction is applied. The original companion
argument is still provenance evidence, not an alternative scientific input.

The four voxel files are listed under [Outputs](#outputs). Below, *s* is the
voxel-count size `--voxel-count-resolution` (0.5 m by default); it applies to
both the XY cells and the height bands, and is recorded as `voxel_size_m` in the
metadata JSON.

The lattice is anchored at `(0,0,0)` in a verified projected metric CRS. XY and
height bins use increasing-coordinate `[lower,upper)` ownership. Raster rows
run north to south; CSV centres identify the same cells. Whole *s* cells must
be contained in the usable AOI, including exclusion boundaries. This deliberately
excludes partially covered cells instead of assigning partial observations to a
full voxel. The export does not depend on complete analysis tiles.

Normalization subtracts the containing DTM cell from source elevation, using
the existing terrain model. Unlike the legacy index path it does not fill DTM
gaps with nearest neighbors or re-quantize normalized heights to LAS scale.
A voxel column is valid only when its full XY cell has finite DTM support.
Finite negative heights are excluded and counted, never clamped into the first
layer. There is no vegetation cutoff, 70 m ceiling or automatic high-point
outlier rejection. Review anomalous heights using the accounting metadata;
exceeding TIFF's 65,535-band capacity fails rather than silently truncating.
The grid covers the usable AOI independently of LAS header extrema. Without an
explicit AOI, the standard runner derives its fallback footprint from header
extrema snapped to the stored coordinate lattice; this corrects pre-quantization
rounding but does not certify arbitrarily corrupted headers. Both entry
points stage and validate the complete voxel quartet before replacing existing
voxel products. Cells require finite support from every intersecting DTM pixel,
including interior pixels when the DTM resolution is finer than *s*.

Admission accounts for raster tile directories outside GDAL's cache as well as
bounded point batches. Extremely wide or tall cubes fail explicitly if the
available memory or conservative uncompressed output-space estimate is exceeded.
These estimates are not a substitute for a container memory limit.

A source with no eligible points emits a headers-only CSV and one zero/NoData
band for `[0, s)`. Missing or non-metric CRS fails the export explicitly, and therefore fails a
production run.

Zero means no eligible recorded points, not absence of vegetation or proof of
laser sampling. Invalid/excluded cells use UInt32 NoData `4294967295` in every
band. Sparse CSV omissions can only be interpreted with the grid and validity
mask. These counts are not GEDI PAVD/PAI. At the default 0.5 m, summing ten aligned
height bands gives a 5 m point-count interval; waveform comparison needs
additional methods.

#### Consumer read-through

In the validated dataset-647 export (EPSG:25832, default 0.5 m voxels), zero-based column 93, row 2
has centre `(474838.25, 5463373.75)`. This bounded GDAL read returns five points
in band 1, representing heights `[0,0.5)` m above the DTM:

```bash
gdal_translate -q -of XYZ -b 1 -srcwin 93 2 1 1 \
  647_voxel_counts.tif /vsistdout/
# 474838.25 5463373.75 5
```

At the same XY cell, bands 1–10 contain `[5,0,0,0,0,0,0,1,0,0]`: six points
over `[0,5)` m, exactly matching the sum of CSV rows at that centre with
`0 <= height_lower_m < 5`. The adjacent row 3, column 93 has zero in band 1
and validity 1; row 0, column 0 has NoData and validity 0. Always read the
validity mask before interpreting absent CSV rows or aggregating raster bands.
Use UInt64 or a wider accumulator for sums across many cells/bands.

#### Scaling and validation

Counting reuses the existing laspy/lazrs streaming stack. Record batches avoid
the unbounded local density of spatial reads. Counts are grouped per batch and
summed in an indexed SQLite scratch database with a fixed cache; each output
raster block/band is written separately through GDAL. No complete cloud, DTM or
dense 3D cube is loaded. Temporary disk and final output still scale with the
number of occupied voxels and the requested raster volume. SQLite sorting and
count reduction trade additional disk I/O for bounded memory. Ensure sufficient
scratch/output disk; write failures prevent promotion in both runners.

The JSON records source/read/included/rejected point accounting, peak batch
size, occupied voxels, scratch bytes and duration. Band height bounds and units
are also embedded in the TIFF, so GDAL PAM sidecars are unnecessary. The
memory-safe runner validates and promotes all four artifacts together with its
existing outputs; the external artifact validator recognizes the optional set.

Run `make test-voxel-counts` for known-count, missing-terrain, height-range,
batch/order invariance, million-point streaming, CLI publication and unchanged
legacy-output tests. `make test-copc-all-instance-dimensions` also compares
voxel CSV and rasters between the original LAZ and lossless COPC.

## Tests

Run the core container tests (CI runs `make test` before publishing an image):

```bash
make test
make test-julia-memory-safe
```

Parameters, voxel counts, point-noise filtering and the Galaxy contract have
their own targets:

```bash
make test-voxel-counts            # includes custom and invalid scientific parameters
make test-denoising
make test-production-contract
make test-galaxy-resource-request
```

Most targets take the image under test as `JULIA_MEMORY_SAFE_IMAGE=<image>`.

Verify ordered LAZ and COPC equality for all three instance dimensions:

```bash
make test-copc-all-instance-dimensions
```

Validate one real dataset against an existing oracle:

```bash
make test-valid-updated-copc-alignment \
  DATASET_ID=150 \
  COPC_LAZ=/path/to/input.copc.laz \
  ORIGINAL_LAZ=/path/to/original.laz \
  AOI_GEOJSON=/path/to/aoi.geojson \
  VALID_UPDATED_DIR=/path/to/oracle \
  INSTANCE_DIMENSION=PredInstance
```

The acceptance comparison requires zero differences in result CSVs, segment
diagnostics, tile GeoJSON, DTM cells, and CHM cells.

## Repository layout

- `src/` — analysis and container entrypoint code
- `galaxy/` — Galaxy wrapper and the TPV resource-rule example
- `validation/` — reusable COPC, AOI, artifact, and oracle validators
- `tests/` — regression and integration test orchestration
- `reference/` — unchanged scientific reference script
- `benchmarks/` — benchmark scripts and compact result data
- `Dockerfile` — published image
- `Dockerfile.julia-memory-safe` — ordered-COPC development image
- `.github/workflows/publish-container.yml` — test and publish the image
- `COPC_REVIEW_GUIDE.md` — detailed implementation explanation
