# 3Dtrees ForestStructure

ForestStructure calculates fixed-tile forest metrics from one LAS/LAZ point
cloud and an optional audited area of interest. It produces tile metrics,
segment diagnostics, vector tiles, a DTM, a CHM, and a tile overview.

The recommended workflow uses a canonical COPC for spatial streaming while
preserving the scientific behavior of the ordered source LAZ. See
[`COPC_REVIEW_GUIDE.md`](COPC_REVIEW_GUIDE.md) for the design, failure fixes,
and validation evidence.

The unchanged scientific reference is
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

## Build

Build the base image and the ordered-COPC image:

```bash
make build IMAGE=3dtrees-foreststructure:v0.1.2
make build-julia-memory-safe \
  JULIA_MEMORY_SAFE_IMAGE=3dtrees-foreststructure:copc-local
```

## Prepare a canonical COPC

Do not use an arbitrary COPC. Build it from the ordered LAS/LAZ so it contains
the validated `OriginalPointIndex` needed to restore source record order:

```bash
FORESTSTRUCTURE_JULIA_IMAGE=3dtrees-foreststructure:copc-local \
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

## Run

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
  3dtrees-foreststructure:copc-local \
  --point-cloud /in/original.copc.laz \
  --original-point-cloud /in/original.laz \
  --aoi-geojson /in/aoi.geojson \
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
Plain LAS/LAZ input is also supported directly.

### Galaxy production contract

Galaxy supplies one full-resolution COPC and optionally an Audit AOI in that
COPC's metric XY coordinate system. The COPC must retain `OriginalPointIndex`
and at least one supported instance dimension; SmartTile supplies global IDs.
An AOI restricts the analysis domain. Without it, the runner reuses the existing
coordinate-lattice-snapped full-extent footprint. Neither mode subsamples points.
DTM generation and the established scientific defaults remain unchanged.

`GALAXY_SLOTS` supplies the CPU allowance; `GALAXY_MEMORY_MB` supplies the job RAM
allocation in MiB. The runner uses the smaller declared/container memory limit
and reserves 15% headroom, subject to the existing 60 GiB default internal cap.
An explicit internal budget may lower it (maximum request 70 GiB). Container
CPU quotas also cap threads. `FORESTSTRUCTURE_THREADS` remains an operator cap;
without a declared CPU allowance or quota, the safe default is one thread.
`FORESTSTRUCTURE_CATALOG_WORKERS` defaults to one; the job destination can set a
measured higher value within the CPU allowance. It is not a user-facing parameter.

The output directory defaults to the Galaxy job working directory and scratch
to `TMPDIR` (or R's private temp directory). A wrapper may supply `--output-dir`
and `--temp-dir` explicitly. Stage a Galaxy dataset under an `input.copc.laz`
symlink because this runner identifies COPC mode from that extension. Never write
beside Galaxy's original input. No dataset ID, CPU, RAM, scratch path or voxel
export toggle belongs in the Galaxy user form. Galaxy must enforce its job CPU,
RAM and scratch allocations; internal estimates are not hard OS limits.

```bash
# Within a job with GALAXY_SLOTS, GALAXY_MEMORY_MB and TMPDIR assigned by Galaxy:
Rscript /opt/foreststructure/run_julia_memory_safe.R \
  --point-cloud input.copc.laz
# Optionally add --aoi-geojson aoi.geojson and --instance-dimension PredInstance_SAT.
```

The container entry point already invokes that script. Both Dockerfiles now use the production runner as their default entry point.
The production container and Galaxy wrapper still need to be released/installed
together; these source changes do not deploy a server tool. Validate integration with
`make test-production-contract JULIA_MEMORY_SAFE_IMAGE=<validated-baseline-image>`.
This compares candidate sources against the pinned baseline, including 15th-field
IDs, a late order key, AOI/full-extent equivalence and mandatory voxel products.

### Catalog and tile concurrency

`FORESTSTRUCTURE_CATALOG_WORKERS` accepts 1–10 workers within
`FORESTSTRUCTURE_THREADS`; the memory-safe default remains one. The standard
CLI exposes the same upper limit through `--catalog-workers`. Tile metrics
reuse that process count. Forked workers use one native thread each to avoid
inherited OpenMP pool deadlocks; a sequential worker can use the full allowance.
Increase workers only after measuring the largest spatial reads: compressed
file size alone does not predict the memory of locally dense chunks.

On the pinned lidR 4.3.2 Linux runtime, bounded fork workers call lidR's existing
chunk engine directly. This avoids repeated future-global discovery while
retaining chunk construction, buffers, callbacks, raster writers and ordered
results. Each child handles one task and exits, releasing its working set.
Other lidR versions or platforms without safe fork support retain the public
catalog scheduler. BLAS pools are limited to one thread in the image to avoid
multiplying native threads inside workers. Scientific algorithms and source
point order are unchanged.

The worker count bounds concurrency, not total memory. The internal budget
and voxel admission estimates do not reserve the sum of all concurrent DTM,
CHM and tile working sets. Keep a hard container limit with equal
`--memory`/`--memory-swap` values, leave headroom, and retain the recorded cgroup
peak. Dense inputs can need fewer workers than available CPUs.

The full-resolution study benchmarks use ten workers for datasets 647/2080 and
three for dense dataset 1805, all under 10 CPUs, a 70 GiB hard limit and a
60 GiB internal budget. See the [measured profiles and exact-parity evidence](
benchmarks/results/2026-10-12-voxel-optimization.json). These are measured
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
This alignment is compatible with the 0.2 m voxel and 0.5 m CHM grids. CHM
reads adapt, but its raster footprints retain the configured baseline size
(normally 60 m). Each read uses the unchanged normalization and voxelization;
only the maximum voxel height at each exact XY is retained across reads
(about 90,000 XY columns per 60 m parent at 0.2 m resolution). Raw points and
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

The production runner discovers available supported dimensions when no selection
is supplied. An explicit selection is strict: any missing or unsupported field
fails before analysis. No supported instance fields also fails. SmartTile owns
global instance-ID assignment; ForestStructure never renumbers or fabricates IDs.
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

## Inputs and outputs

The AOI is interpreted in the point cloud's local XY coordinate space and is
not reprojected. The GeoJSON-to-GeoPackage round trip must preserve topology,
exact XY vertices, and exact bounds. An absolute area difference of at most
`1e-6` square coordinate units is allowed solely for floating-point summation
noise. Without an AOI, the full point-cloud extent is tiled.

Production outputs use the fixed prefix `foreststructure`; no dataset ID is
required. The legacy `--dataset-id` option remains only to replay old validation
commands and is not a Galaxy input. Below, `<dataset>` denotes this output prefix.

For each available instance dimension, the run writes:

- `<dataset>_<dimension>_results.csv`
- `<dataset>_<dimension>_segment_diagnostics.csv`
- `<dataset>_<dimension>_tiles.geojson`

Shared outputs are:

- `<dataset>_dtm.tif`
- `<dataset>_chm.tif`
- `<dataset>_tiles.png`
- `<dataset>_aoi_conversion.json`
- `<dataset>_julia_memory_safe_run.json`
- `<dataset>_performance.csv`

Outputs are computed in a staging directory, validated as a complete set, and
then promoted. A failed run does not publish a partial result set.

## Validation tools

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
  3dtrees-foreststructure:copc-local \
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

## Tests

Run the core container tests:

```bash
make test
make test-julia-memory-safe
```

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
- `validation/` — reusable COPC, AOI, artifact, and oracle validators
- `tests/` — regression and integration test orchestration
- `reference/` — unchanged scientific reference script
- `benchmarks/` — benchmark scripts and compact result data
- `Dockerfile.julia-memory-safe` — ordered-COPC image
- `COPC_REVIEW_GUIDE.md` — detailed implementation explanation


## Full-resolution voxel point-count export

The production (`run_julia_memory_safe.R`) runner always exports voxel counts;
there is no opt-out. `--export-voxel-counts` remains accepted for old invocations.
The legacy `run.R` entry point retains its opt-in flag. The optional
`--voxel-count-chunk-points` (default 250,000) bounds the number of source records
per count batch; the memory-safe runner lowers it further to fit its remaining
memory allowance. Existing index/CHM filters and output schemas are unchanged.

The export reuses the generated global DTM and usable Audit AOI. Supply the
original full-resolution LAS/LAZ or a lossless COPC with verified conversion
provenance. Every record is read, including repeated coordinates; no point
sampling or COPC level-of-detail restriction is applied. The original companion
argument is still provenance evidence, not an alternative scientific input.

For dataset 150 the additional artifacts are:

| File | Contents |
| --- | --- |
| `150_voxel_counts.csv` | Occupied voxels: `x_center,y_center,height_lower_m,height_upper_m,point_count` |
| `150_voxel_counts.tif` | UInt32 counts, 0.5 m XY pixels; band b contains heights `[(b-1)*0.5,b*0.5)` metres above DTM |
| `150_voxel_counts_validity.tif` | 0: excluded AOI/boundary cell; 1: valid domain; 2: invalid DTM |
| `150_voxel_counts.json` | Grid/CRS, source and DTM hashes, policies, point accounting and resource observations |

The lattice is anchored at `(0,0,0)` in a verified projected metric CRS. XY and
height bins use increasing-coordinate `[lower,upper)` ownership. Raster rows
run north to south; CSV centres identify the same cells. Whole 0.5 m cells must
be contained in the usable AOI, including exclusion boundaries. This deliberately
excludes partially covered cells instead of assigning partial observations to a
full voxel. The export does not depend on complete 20 m Analysis Tiles.

Normalization subtracts the containing DTM cell from source elevation, using
the existing terrain model. Unlike the legacy index path it does not fill DTM
gaps with nearest neighbors or re-quantize normalized heights to LAS scale.
A voxel column is valid only when its full XY cell has finite DTM support.
Finite negative heights are excluded and counted, never clamped into the first
layer. There is no 0.5 m vegetation cutoff, 70 m ceiling or automatic high-point
outlier rejection. Review anomalous heights using the accounting metadata;
exceeding TIFF's 65,535-band capacity fails rather than silently truncating.
The grid covers the usable AOI independently of LAS header extrema. Without an
explicit AOI, the standard runner derives its fallback footprint from header
extrema snapped to the stored coordinate lattice; this corrects pre-quantization
rounding but does not certify arbitrarily corrupted headers. Both entry
points stage and validate the complete voxel quartet before replacing existing
voxel products. Cells require finite support from every intersecting DTM pixel,
including interior pixels when the DTM resolution is finer than 0.5 m.

Admission accounts for raster tile directories outside GDAL's cache as well as
bounded point batches. Extremely wide or tall cubes fail explicitly if the
available memory or conservative uncompressed output-space estimate is exceeded.
These estimates are not a substitute for a container memory limit.

A source with no eligible points emits a headers-only CSV and one zero/NoData
band for 0–0.5 m. Missing or non-metric CRS fails the export explicitly, and therefore fails a
production run.

Zero means no eligible recorded points, not absence of vegetation or proof of
laser sampling. Invalid/excluded cells use UInt32 NoData `4294967295` in every
band. Sparse CSV omissions can only be interpreted with the grid and validity
mask. These counts are not GEDI PAVD/PAI. Summing ten aligned height bands gives
a 5 m point-count interval; waveform comparison needs additional methods.

### Consumer read-through

In the validated dataset-647 export (EPSG:25832), zero-based column 93, row 2
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

### Scaling and validation

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
