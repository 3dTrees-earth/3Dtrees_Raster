#!/usr/bin/env bash
# Treat every failed command, unset variable, or pipeline component as failure.
set -euo pipefail

# Resolve repository helpers without depending on the invocation directory.
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Pin/override the two images involved in analysis and canonical COPC creation.
image_name="${FORESTSTRUCTURE_JULIA_IMAGE:-3dtrees-foreststructure:copc-local}"
smart_tile_image="${FORESTSTRUCTURE_SMART_TILE_IMAGE:-3dtrees-smart-tile:34f492a-v2.3.1}"
# Isolate all generated fixtures, work files, and outputs in one private root.
test_root="$(mktemp -d)"
# Preserve evidence only when explicitly requested during debugging.
cleanup() {
  if [[ "${FORESTSTRUCTURE_KEEP_TEST_ROOT:-0}" == "1" ]]; then
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
  --env FORESTSTRUCTURE_FIXTURE_EPSG=32632 \
  --user "$(id -u):$(id -g)" \
  --entrypoint Rscript \
  --volume "${repo_dir}/tests:/tests:ro" \
  --volume "${test_root}/input:/generated" \
  "${image_name}" \
  /tests/generate_julia_memory_safe_fixture.R \
  /generated/distinct.laz /generated/aoi.geojson distinct-dimensions

# Convert that ordered LAZ to canonical COPC under a small regression-test budget.
FORESTSTRUCTURE_JULIA_IMAGE="${image_name}" \
FORESTSTRUCTURE_SMART_TILE_IMAGE="${smart_tile_image}" \
FORESTSTRUCTURE_TEST_CPUS=2 \
FORESTSTRUCTURE_TEST_MEMORY_GIB=4 \
  bash "${repo_dir}/tests/build_original_order_copc.sh" \
  "${test_root}/input/distinct.laz" \
  "${test_root}/input/distinct.copc.laz" \
  "${test_root}/input/order-validation.json"

# Run one ForestStructure input through an identical resource/science envelope.
run_foreststructure() {
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
    --env FORESTSTRUCTURE_THREADS=2 \
    --env FORESTSTRUCTURE_CATALOG_WORKERS=1 \
    --volume "${test_root}/input:/in:ro" \
    --volume "${output_dir}:/out" \
    --volume "${work_dir}:/work" \
    "${image_name}" \
    --point-cloud "/in/${point_cloud_name}" \
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
run_foreststructure \
  "${test_root}/input/distinct.laz" \
  "${test_root}/ordered-output" \
  "${test_root}/ordered-work" \
  2>&1 | tee "${test_root}/ordered.log"
# Then process canonical COPC while supplying LAZ only as identity evidence.
run_foreststructure \
  "${test_root}/input/distinct.copc.laz" \
  "${test_root}/copc-output" \
  "${test_root}/copc-work" \
  --original-point-cloud /in/distinct.laz \
  --target-chunk-points 500 \
  2>&1 | tee "${test_root}/copc.log"

# Require byte identity for every tabular/vector artifact in all three dimensions.
for dimension in PredInstance PredInstance_SAT PredInstance_FM; do
  for artifact in results.csv segment_diagnostics.csv tiles.geojson; do
    cmp \
      "${test_root}/ordered-output/150_${dimension}_${artifact}" \
      "${test_root}/copc-output/150_${dimension}_${artifact}"
  done
done
# The observation cube is independent of both COPC order and adaptive windows.
cmp "${test_root}/ordered-output/150_voxel_counts.csv" \
    "${test_root}/copc-output/150_voxel_counts.csv"
for raster in voxel_counts.tif voxel_counts_validity.tif; do
  cmp "${test_root}/ordered-output/150_${raster}" \
      "${test_root}/copc-output/150_${raster}"
done
# DTM and CHM must also be byte-identical between ordered LAZ and COPC.
for raster in dtm.tif chm.tif; do
  cmp \
    "${test_root}/ordered-output/150_${raster}" \
    "${test_root}/copc-output/150_${raster}"
done

# Guard fixture quality: distinct dimension inputs must not collapse to the same
# results, otherwise the equality checks above would not prove dimension routing.
if cmp -s \
  "${test_root}/copc-output/150_PredInstance_results.csv" \
  "${test_root}/copc-output/150_PredInstance_SAT_results.csv"; then
  echo "PredInstance and PredInstance_SAT unexpectedly produced identical results" >&2
  exit 1
fi
# Apply the same fixture-quality guard to SAT versus foundation-model values.
if cmp -s \
  "${test_root}/copc-output/150_PredInstance_SAT_results.csv" \
  "${test_root}/copc-output/150_PredInstance_FM_results.csv"; then
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
' "${test_root}/copc-output/150_julia_memory_safe_run.json" >/dev/null

# Reach this line only after both differential science and provenance checks pass.
echo "COPC all-instance-dimension differential acceptance passed"
