#!/usr/bin/env bash
# Treat every failed command, unset variable, or pipeline component as failure.
set -euo pipefail

# Resolve repository helpers without depending on the invocation directory.
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Pin/override the two images involved in analysis and canonical COPC creation.
image_name="${RASTER_IMAGE:-3dtrees-raster:test}"
smart_tile_image="${RASTER_SMART_TILE_IMAGE:-3dtrees-smart-tile:34f492a-v2.3.1}"
# Isolate all generated fixtures, work files, and outputs in one private root.
test_root="$(mktemp -d)"
# Preserve evidence only when explicitly requested during debugging.
cleanup() {
  if [[ "${RASTER_KEEP_TEST_ROOT:-0}" == "1" ]]; then
    echo "Preserving test root: ${test_root}"
  else
    rm -rf "${test_root}"
  fi
}
# Ensure normal exits, failures, and signals all follow the same cleanup policy.
trap cleanup EXIT
echo "COPC all-dimension test root: ${test_root}"
# Separate immutable inputs from ordered-LAZ and COPC execution outputs/work.
mkdir -p \
  "${test_root}/input" \
  "${test_root}/ordered-output" \
  "${test_root}/ordered-work" \
  "${test_root}/copc-output" \
  "${test_root}/copc-work"

# Generate a deterministic fixture in which all three instance dimensions differ.
docker run --rm --network none \
  --env RASTER_FIXTURE_EPSG=32632 \
  --user "$(id -u):$(id -g)" \
  --entrypoint Rscript \
  --volume "${repo_dir}/tests:/tests:ro" \
  --volume "${test_root}/input:/generated" \
  "${image_name}" \
  /tests/generate_test_cloud.R \
  /generated/distinct.laz /generated/aoi.geojson distinct-dimensions

# Convert that ordered LAZ to canonical COPC under a small regression-test budget.
RASTER_IMAGE="${image_name}" \
RASTER_SMART_TILE_IMAGE="${smart_tile_image}" \
RASTER_TEST_CPUS=2 \
RASTER_TEST_MEMORY_GIB=4 \
  bash "${repo_dir}/scripts/build_original_order_copc.sh" \
  "${test_root}/input/distinct.laz" \
  "${test_root}/input/distinct.copc.laz" \
  "${test_root}/input/order-validation.json"

# Run one raster input through an identical resource/science envelope.
run_raster() {
  # The first three arguments identify input, published output, and scratch work.
  local point_cloud="$1"
  local output_dir="$2"
  local work_dir="$3"
  local point_cloud_name
  # All fixture files share one read-only mount, so pass only the basename inside.
  point_cloud_name="$(basename "${point_cloud}")"
  # Remaining arguments carry COPC-only options such as the original companion.
  shift 3
  # Network isolation and hard CPU/memory limits keep the test deterministic.
  docker run --rm --network none \
    --cpus 2 \
    --memory 4g \
    --memory-swap 4g \
    --user "$(id -u):$(id -g)" \
    --env RASTER_THREADS=2 \
    --env RASTER_CATALOG_WORKERS=1 \
    --volume "${test_root}/input:/in:ro" \
    --volume "${output_dir}:/out" \
    --volume "${work_dir}:/work" \
    "${image_name}" \
    --denoise none --point-cloud "/in/${point_cloud_name}" \
    --instance-dimension PredInstance,PredInstance_SAT,PredInstance_FM \
    "$@" \
    --aoi-geojson /in/aoi.geojson \
    --dataset-id 150 \
    --output-dir /out \
    --temp-dir /work \
    --export-voxel-counts \
    --voxel-count-chunk-points 257 \
    --memory-budget-gib 3
}

# First produce the ordered-LAZ oracle using the same candidate implementation.
run_raster \
  "${test_root}/input/distinct.laz" \
  "${test_root}/ordered-output" \
  "${test_root}/ordered-work" \
  2>&1 | tee "${test_root}/ordered.log"
# Then process canonical COPC while supplying LAZ only as identity evidence.
run_raster \
  "${test_root}/input/distinct.copc.laz" \
  "${test_root}/copc-output" \
  "${test_root}/copc-work" \
  --original-point-cloud /in/distinct.laz \
  --target-chunk-points 500 \
  2>&1 | tee "${test_root}/copc.log"

# Require byte identity for every tabular/vector artifact in all three dimensions.
for dimension in PredInstance PredInstance_SAT PredInstance_FM; do
  for artifact in "tile_metrics/${dimension}_tile_metrics.csv" "segment_diagnostics/${dimension}_segment_diagnostics.csv" "tiles/${dimension}_tiles.geojson"; do
    cmp \
      "${test_root}/ordered-output/${artifact}" \
      "${test_root}/copc-output/${artifact}"
  done
done
# The observation cube is independent of both COPC order and adaptive windows.
cmp "${test_root}/ordered-output/150_voxel_counts.csv" \
    "${test_root}/copc-output/150_voxel_counts.csv"
for raster in voxel_counts.tif voxel_counts_validity.tif; do
  cmp "${test_root}/ordered-output/150_${raster}" \
      "${test_root}/copc-output/150_${raster}"
done
# The DTM must be byte-identical. The CHM grid follows the header bounds, and a
# COPC's header can reach one cell further at the far edges, where it then holds
# the points on that edge. Every shared cell must be identical.
cmp "${test_root}/ordered-output/150_dtm.tif" "${test_root}/copc-output/150_dtm.tif"
docker run --rm --network none --user "$(id -u):$(id -g)" --entrypoint python3 \
  --volume "${test_root}:/check:ro" "${image_name}" -c '
import numpy as np
from osgeo import gdal
gdal.UseExceptions()
for name in ("150_chm.tif", "150_dsm.tif"):
    a, b = (gdal.Open(f"/check/{side}-output/{name}") for side in ("ordered", "copc"))
    ga, gb = a.GetGeoTransform(), b.GetGeoTransform()
    assert ga[:3] + ga[4:] == gb[:3] + gb[4:], name
    rows = round((gb[3] - ga[3]) / ga[1])
    assert rows in (0, 1) and b.RasterXSize - a.RasterXSize in (0, 1), name
    assert b.RasterYSize - a.RasterYSize == rows, name
    x = a.GetRasterBand(1).ReadAsArray()
    y = b.GetRasterBand(1).ReadAsArray()[rows:rows + a.RasterYSize, :a.RasterXSize]
    assert np.array_equal(x, y), name
'

# Guard fixture quality: distinct dimension inputs must not collapse to the same
# results, otherwise the equality checks above would not prove dimension routing.
if cmp -s \
  "${test_root}/copc-output/tile_metrics/PredInstance_tile_metrics.csv" \
  "${test_root}/copc-output/tile_metrics/PredInstance_SAT_tile_metrics.csv"; then
  echo "PredInstance and PredInstance_SAT unexpectedly produced identical results" >&2
  exit 1
fi
# Apply the same fixture-quality guard to SAT versus foundation-model values.
if cmp -s \
  "${test_root}/copc-output/tile_metrics/PredInstance_SAT_tile_metrics.csv" \
  "${test_root}/copc-output/tile_metrics/PredInstance_FM_tile_metrics.csv"; then
  echo "PredInstance_SAT and PredInstance_FM unexpectedly produced identical results" >&2
  exit 1
fi

# Verify that provenance records COPC mode, order restoration, and all dimensions.
jq -e '
  .requested_is_copc == true and
  .point_order_dimension == "OriginalPointIndex" and
  .adaptive_chunking.enabled == true and
  .adaptive_chunking.method == "mean_density" and
  .adaptive_chunking.target_points == 500 and
  .adaptive_chunking.selected_chunk_size_m > 0 and
  .adaptive_chunking.selected_chunk_size_m < 60 and
  .adaptive_chunking.observed.chm.chunks > 1 and
  .adaptive_chunking.observed.segmentation.chunks > 1 and
  .processed_instance_dimensions == [
    "PredInstance", "PredInstance_SAT", "PredInstance_FM"
  ]
' "${test_root}/copc-output/150_run.json" >/dev/null

# Reach this line only after both differential science and provenance checks pass.
echo "COPC all-instance-dimension differential acceptance passed"
