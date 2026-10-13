# Ordered-COPC implementation review guide

This guide is for reviewing the COPC changes without having to reconstruct the
design from individual commits. The scientific goal is strict: processing a
canonical COPC must reproduce the ordered-LAZ/Julia result while reading point
data through COPC's spatial hierarchy.

## The problem being solved

Untwine reorganizes LAS/LAZ records into COPC spatial order. That is desirable
for bounded spatial reads, but record order is no longer the source order seen
by the original Julia/R workflow. Some lidR operations resolve ties using input
order, so feeding a generic COPC can change a small number of downstream values
even when every XYZ and instance value survived conversion.

A second, independent failure occurs in lidR's optimized DTM TIN path when an
XY coordinate cannot be represented by the LAS header's scale and offset as a
signed 32-bit integer. Dataset 2011, for example, used 0.000001 m scales around
easting 550,000 m and failed with:

```text
Internal error in C_interpolate_delaunay: xy coordinates were not converted to integer.
Scale factors are likely to be invalid.
```

## End-to-end data flow

```text
ordered LAS/LAZ
    |
    | add OriginalPointIndex = zero-based source record position
    v
indexed ordered LAZ
    |
    | Untwine spatial reordering
    v
canonical COPC
    |
    | validate header + every index exactly once + tuple fingerprints
    v
published canonical COPC
    |
    | lidR spatial/chunk/window reads from COPC
    v
bounded in-memory chunk
    |
    | sort chunk rows by OriginalPointIndex
    v
Julia-compatible point order -> unchanged scientific operations
```

The original LAS/LAZ can be mounted with `--original-point-cloud`, but it is
used only to validate point count/XYZ bounds and to preserve the historical
source filename in output metadata. Scientific point reads always use the COPC.

## Adaptive window policy

CHM and every independent `PredInstance`, `PredInstance_SAT`, and
`PredInstance_FM` pass use the same density-derived COPC window width. From the
header point count and XY area, the runner computes the mean point density and
then floors `sqrt(target_points / density)` to whole metres. The result is
bounded to 1–60 m and may be reduced further by `--chunk-size`. The default
target is 50,000,000 points; `--target-chunk-points 0` restores fixed-size mode.
Non-COPC LAS/LAZ remains fixed-size.

This is intentionally an approximate target. Spatial density is not uniform,
so an individual window can contain more than the target without invalidating
the run. Each window records its observed point count, emits a structured
warning when oversized, and contributes minimum, maximum, count, and oversize
totals to performance and provenance artifacts.

CHM and segmentation retain zero buffer. Their catalogs share the absolute
origin `(0, 0)` and whole-metre window widths, so boundaries stay aligned with
both the 0.2 m voxel lattice and the 0.5 m CHM raster. CHM produces one raster
per window and streams their VRT mosaic to the unchanged dataset-level CHM.
Segmentation writes only compact voxel, layer, and apex summaries to SSD, drops
the raw LAS and temporary tables before the next read, and reuses the existing
bucket-at-a-time disk reducer.

Changing window boundaries must not change order-dependent results. Compact
COPC summaries therefore carry `OriginalPointIndex`: apex ties use maximum
height then minimum original index, duplicate voxels retain their earliest
original index, and the merged voxel table is restored to global source order
before PCA. These internal ordinal fields never enter public scientific output
schemas.

DTM is unchanged: it continues to use 60 m chunks with a 20 m buffer. DTM
interpolation requires neighboring ground points, and this change targets the
CHM/segmentation memory peaks rather than changing terrain science.
The COPC path restores original row order before DTM and then drops the ordinal
column because PTD/TIN no longer needs it. Non-collinearity guards scan bounded
blocks rather than allocating full duplicate XY tables. Neither optimization
changes the selected points, row order, DTM parameters, or terrain algorithm.

## File-by-file review map

### `tests/add_original_point_index.py`

- Streams the ordered source instead of loading the complete cloud.
- Adds a UInt64 `OriginalPointIndex` equal to the zero-based source row.
- Moves `OriginalPointIndex`, `PredInstance`, `PredInstance_SAT`, and
  `PredInstance_FM` into lidR's first nine selectively readable ExtraBytes.
- Copies every original dimension value and ExtraByte definition unchanged.

This is the only step that defines original order. It runs before Untwine.

### `tests/build_original_order_copc.sh`

- Runs indexing, Untwine conversion, and validation in isolated containers.
- Applies the same configurable CPU and hard memory limits to every stage.
- Builds beside the destination and publishes with one rename only after all
  validation checks pass.
- Never publishes a partially built or invalid candidate.

### `validation/validate_indexed_copc_streaming.py`

- Verifies exact point-count, stored XYZ extrema, scale and offset equality,
  plus semantic CRS equivalence with axis order preserved. Equivalent WKT1 and
  WKT2 serializations do not cause a false rejection.
  Source headers may contain bounds calculated before coordinate quantization;
  the validator records this discrepancy and compares extrema from the existing
  point scan. Candidate header bounds must match its stored extrema exactly.
  The optional runtime original-companion check remains header-only and snaps
  bounds to the stored coordinate lattice to handle pre-quantization rounding.
  It cannot replace this complete identity validation; retain the source hashes
  and streaming comparison report as external provenance.
- Uses a Boolean presence array to prove that all indices from `0` through `n - 1` occur exactly
  once, with no missing, duplicate, or out-of-range values.
- Compares order-independent dual fingerprints over
  `(OriginalPointIndex, raw X, raw Y, raw Z, instance value)` for every supported
  instance dimension present in the source.
- Reads both files in bounded chunks.

The presence array gives an exact index-completeness proof and occupies one byte
per source point (about 1.33 GiB for 1.425 billion records). Fingerprints bind each index
to its geometry and scientific instance value without assuming physical order.

### `src/run_julia_memory_safe.R`

- Accepts canonical `.copc.laz` as the actual scientific input.
- Rejects COPC without `OriginalPointIndex`.
- Reads directly from COPC when the order key and requested instance dimensions
  are within lidR's first nine ExtraBytes. Otherwise it reuses the bounded,
  verified projection helper to create one compact indexed scratch LAZ, retaining
  the original order key and selected IDs. The pinned COPC reader truncates extras
  after nine even with selector `0`; a wildcard is not a valid workaround.
- Skips LAX creation for direct COPC. The high-ordinal fallback requires a LAX
  beside its scratch projection. Voxel counts always read the original source.
- Uses an optional AOI, otherwise the full-extent footprint; exports voxel counts
  unconditionally; defaults to the `foreststructure` output prefix. Galaxy owns
  CPU/memory allocation and job directories. Explicit missing instance selections
  fail, and auto discovery requires at least one supported instance field.
- Optionally validates an original LAS/LAZ companion by header count and bounds.
- Records COPC use, order restoration, selectors, and companion validation in
  the run provenance JSON.

### `src/run.R`

- `restore_original_point_order()` validates and sorts each bounded COPC read by
  `OriginalPointIndex` before DTM, CHM, segmentation, or tile calculations.
- `catalog_selection_for_dimensions()` adds only the required ExtraByte selector
  to the existing XYZ catalog/window request.
- DTM, CHM, segment, and tile functions receive the order key explicitly. There
  is no hidden global COPC mode.
- Plain LAS/LAZ passes `NULL`, so its established order and behavior are unchanged.
- `adaptive_chunk_plan()` derives one bounded whole-metre COPC window size from
  header density for CHM and segmentation only.
- Per-window observation records make the approximate target and local-density
  overshoots auditable without retaining point tables in memory.
- CHM adapts point reads, not the baseline raster footprints. Each smaller read
  is normalized and voxelized unchanged, then reduced to maximum height per
  exact voxel XY. Combining these compact columns retains the baseline occupied
  XY extent while bounding memory (about 90,000 rows for a 60 m parent). The
  original numeric-resolution `p2r()` runs once per parent, whose raster is
  written to disk and merged in the original catalog order. No fixed raster
  template, coordinate rounding, interpolation, or new rasterization is used.
  This matters because a voxel on a small raster's outer edge can otherwise be
  assigned to a different cell than in the original 60 m raster. The regression
  in `tests/test_chm_adaptive_edges.R` exercises this through real spatial reads.
- The segment disk reducer restores duplicate-voxel and PCA order globally from
  `OriginalPointIndex`, so adaptive boundaries do not change tie behavior.
- `legacy_tin()` retries only the known lidR Delaunay integer-conversion error.
  It explicitly uses lidR 4.3.2's floating-point `tDelaunay()` and
  `tInterpolate()` path with the original LAS coordinates and header.

The DTM retry does not rescale or rewrite the source LAS/LAZ or COPC. It does
not alter XY, Z, segmentation, CHM input storage, or the configured DTM
resolution. It changes only the implementation used to complete the failed
Delaunay interpolation and preserves `tin()`'s outside-hull extrapolation.

### `tests/test_copc_all_instance_dimensions.sh`

Creates one controlled fixture with deliberately different values in all three
instance dimensions. It runs ordered LAZ and canonical COPC through the same
image and requires byte-identical outputs for every dimension plus DTM and CHM.
It also proves that the fixture dimensions produce different results from one
another, preventing a vacuous routing test.

### `tests/test_valid_updated_copc_alignment.sh`

Runs a real canonical COPC under explicit CPU/memory limits, asserts COPC and
order-restoration provenance, and compares:

- result CSV values at zero tolerance;
- segment diagnostics at zero tolerance;
- output GeoJSON exactly, excluding only path metadata;
- DTM and CHM geometry, populated-cell masks, and cell values at zero tolerance.

The detailed comparisons live in
`validation/compare_valid_updated_science.py` and
`validation/compare_valid_updated_raster.R`. Both write their report before failing,
so an acceptance failure retains exact row/cell counts and sampled differences.

### Cohort acceptance scripts

- `tests/test_valid_updated_cohort.sh` proves that selected legacy `valid` or
  `pending_xy_rescale` results really differ from `valid_updated`, converts the
  ordered source to canonical COPC, and requires the new output to match the
  oracle.
- `tests/test_valid_updated_oracle_cohort.sh` applies the same conversion and
  zero-tolerance checks to datasets that have a `valid_updated` oracle but no
  legacy `valid` result, while recording conversion/runtime/resource evidence.
- `tests/compare_las_point_order.py` is a diagnostic tool used to establish the
  original root cause: it distinguishes physical sequence differences from
  changes in the unordered point-value multiset. Canonical publication uses the
  streaming validator instead because the diagnostic loads complete files.

## What did not change

- Julia's scientific parameter values.
- The released base image's R, lidR, terra, sf, or GDAL versions.
- Tree acceptance thresholds, voxel definitions, raster resolutions, or AOI
  conversion semantics.
- Original LAS/LAZ or COPC files during analysis.
- Plain LAS/LAZ behavior when no point-order dimension is supplied.

The patch changes storage access, deterministic record-order reconstruction,
input validation/provenance, and a narrowly matched DTM error recovery path.

## Acceptance evidence

The automated differential fixture covers `PredInstance`, `PredInstance_SAT`,
and `PredInstance_FM`. The ordered-COPC replay used dataset 2011:

- 23,700,806 source points;
- every `OriginalPointIndex` present exactly once;
- the historical DTM integer-conversion error was reproduced on the old image;
- the patched fallback triggered and the complete pipeline produced nine files;
- 18 segment-diagnostic rows matched the fresh ordered-LAZ control exactly;
- DTM matched across 361 cells with zero differing cells;
- CHM matched across 1,369 cells with zero differing cells;
- COPC peak cgroup memory was approximately 4.59 GiB under a 30 GiB hard limit.

Dataset 2011 had no previously promoted `valid_updated` artifact because its old
run failed at DTM. Therefore its reference was generated fresh from the ordered
original LAZ with the same patched scientific implementation. This isolates the
COPC conversion/order-restoration path while exercising the exact former failure.

### Legacy-TIN retry replay

The legacy floating-point retry introduced after that ordered-COPC validation
was replayed separately on dataset 2153 on 2026-09-02:

- Input `8101.laz` contained 40,453,840 points in 220,297,505 bytes; SHA-256
  `7bca4ef1fe055ab5f1a826050d0c85330f7f17fc4cd99580769ab0d08656a346`.
- The historical image reproduced the exact
  `C_interpolate_delaunay: xy coordinates were not converted to integer`
  failure during DTM and exited with status 1.
- The slow-retry image logged the legacy floating-point fallback, completed
  DTM, passed internal validation, and atomically promoted 12 artifacts.
- The controlled run used 10 CPUs, one LAScatalog worker, a 30 GiB hard limit,
  and a 25 GiB internal budget. Analysis took 465.808 seconds, including
  127.044 seconds for DTM, and peak cgroup memory was 6,249.4 MiB.

This replay proves that the slow retry recovers a real former failure without
rescaling XY. It is not a zero-tolerance `valid_updated` acceptance result:
dataset 2153 is absent from the confirmed oracle cohort. A diagnostic comparison
with its previously promoted generic-COPC result retained identical raster
geometry and populated masks but found 800 differing DTM cells (mean absolute
difference 0.0002 m, maximum 0.0282 m), two differing CHM cells (maximum
0.2 m), 42 differing SAT segment fields, and 63 differing FM segment fields.
Small numerical changes are expected because the retry now triangulates the
original floating-point XY coordinates rather than temporarily rescaling them
to satisfy the optimized integer TIN. Preserving the measured coordinates and
LAS header makes this fallback scientifically more faithful than the former
rescaling workaround.

That comparison also changes input order and runtime concurrency, so it cannot
attribute the differences solely to the TIN implementation. A confirmed
`valid_updated` oracle replay remains required before production promotion.

The complete structured evidence is stored in
`benchmarks/results/2026-09-02-dataset2153-slow-tin-validation.json`.

## Reviewer commands

Build the image:

```bash
make build-julia-memory-safe JULIA_MEMORY_SAFE_IMAGE=3dtrees-foreststructure:copc-local
```

Run the controlled three-dimension differential:

```bash
make test-copc-all-instance-dimensions
```

Build and validate a canonical COPC:

```bash
make build-original-order-copc \
  ORIGINAL_LAZ=/path/to/ordered.laz \
  COPC_LAZ=/path/to/canonical.copc.laz
```

Run a dataset/oracle comparison:

```bash
make test-valid-updated-copc-alignment \
  DATASET_ID=123 \
  COPC_LAZ=/path/to/canonical.copc.laz \
  ORIGINAL_LAZ=/path/to/ordered.laz \
  AOI_GEOJSON=/path/to/aoi.geojson \
  VALID_UPDATED_DIR=/path/to/oracle \
  INSTANCE_DIMENSION=PredInstance
```
