#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
image_name="${FORESTSTRUCTURE_IMAGE:-3dtrees-foreststructure:test}"
docker_limits=(--cpus "${FORESTSTRUCTURE_CPUS:-10}" --memory "${FORESTSTRUCTURE_MEMORY:-70g}" --memory-swap "${FORESTSTRUCTURE_MEMORY:-70g}")
input_dir="$(mktemp -d)"
results_dir="$(mktemp -d)"
cleanup() {
  local status="$?"
  if [[ "${status}" != 0 || "${FORESTSTRUCTURE_KEEP_TEST_ROOT:-0}" == 1 ]]; then
    echo "Container gate inputs: ${input_dir}"
    echo "Container gate results and logs: ${results_dir} (exit ${status})"
  else
    rm -rf "${input_dir}" "${results_dir}"
  fi
}
trap cleanup EXIT

if [[ "${FORESTSTRUCTURE_SKIP_BUILD:-0}" != "1" ]]; then
  docker build --tag "${image_name}" "${repo_dir}"
fi
# Exercise the optional observational export in the standard publish gate too.
docker run --rm --network none "${docker_limits[@]}" \
  --env PYTHONDONTWRITEBYTECODE=1 --entrypoint python3 \
  --volume "${repo_dir}:/workspace:ro" "${image_name}" \
  /workspace/tests/test_voxel_counts.py
python3 "${repo_dir}/tests/test_output_artifact_validation.py"
python3 "${repo_dir}/tests/test_copc_task_counts.py"
PYTHONDONTWRITEBYTECODE=1 python3 "${repo_dir}/tests/test_galaxy_resource_request.py"
docker run --rm --network none \
  "${docker_limits[@]}" \
  --entrypoint Rscript \
  --volume "${repo_dir}:/workspace:ro" \
  --workdir /workspace/src \
  "${image_name}" \
  -e 'Sys.setenv(FORESTSTRUCTURE_SOURCE_ONLY = "1"); source("run.R"); source("../tests/test_runtime_helpers.R"); source("../tests/test_catalog_scheduler.R"); source("../tests/test_resource_scheduler.R"); source("../tests/test_resource_memory.R")'
docker run --rm --network none \
  "${docker_limits[@]}" \
  --entrypoint Rscript \
  --volume "${repo_dir}:/workspace:ro" \
  --workdir /workspace \
  "${image_name}" \
  tests/test_tile_scheduler.R
docker run --rm --network none \
  "${docker_limits[@]}" \
  --entrypoint Rscript \
  --volume "${repo_dir}:/workspace:ro" \
  --workdir /workspace/src \
  "${image_name}" \
  -e 'Sys.setenv(FORESTSTRUCTURE_SOURCE_ONLY = "1"); source("run.R"); source("../tests/test_failed_dataset_regressions.R")'
docker run --rm --network none \
  "${docker_limits[@]}" \
  --entrypoint Rscript \
  --volume "${repo_dir}:/workspace:ro" \
  --workdir /workspace/src \
  "${image_name}" \
  -e 'Sys.setenv(FORESTSTRUCTURE_SOURCE_ONLY = "1"); source("run_julia_memory_safe.R"); source("../tests/test_julia_memory_safe.R"); source("../tests/test_instance_read_groups.R"); source("../tests/test_resource_validation.R")'
docker run --rm --network none \
  "${docker_limits[@]}" \
  --entrypoint Rscript \
  --volume "${repo_dir}:/workspace:ro" \
  --workdir /workspace/src \
  "${image_name}" \
  -e 'Sys.setenv(FORESTSTRUCTURE_SOURCE_ONLY = "1"); source("run.R"); source("../tests/test_adaptive_chunks.R"); source("../tests/test_chm_adaptive_edges.R")'
docker run --rm --network none \
  "${docker_limits[@]}" \
  --entrypoint Rscript \
  --volume "${repo_dir}:/workspace:ro" \
  --workdir /workspace/src \
  "${image_name}" \
  -e 'Sys.setenv(FORESTSTRUCTURE_SOURCE_ONLY = "1"); source("run.R"); source("../tests/test_dtm_coordinate_scale.R")'
docker run --rm --network none \
  "${docker_limits[@]}" \
  --entrypoint Rscript \
  --volume "${repo_dir}:/workspace:ro" \
  --workdir /workspace/src \
  "${image_name}" \
  -e 'Sys.setenv(FORESTSTRUCTURE_SOURCE_ONLY = "1"); source("run.R"); source("../tests/test_raster_merge_semantics.R")'
docker run --rm --network none \
  "${docker_limits[@]}" \
  --entrypoint Rscript \
  --volume "${repo_dir}:/workspace:ro" \
  --workdir /workspace/src \
  "${image_name}" \
  -e 'Sys.setenv(FORESTSTRUCTURE_SOURCE_ONLY = "1"); source("run.R"); source("aoi_conversion.R"); source("../tests/test_aoi_conversion.R")'
docker run --rm --network none \
  "${docker_limits[@]}" \
  --user "$(id -u):$(id -g)" \
  --entrypoint Rscript \
  --volume "${repo_dir}/tests:/tests:ro" \
  --volume "${input_dir}:/generated" \
  "${image_name}" \
  /tests/generate_fixture.R \
  /generated/point_cloud.laz \
  /generated/point_cloud_segmented.laz \
  /generated/aoi_with_exclusion.gpkg

help_path="${results_dir}/help.txt"
docker run --rm --network none "${docker_limits[@]}" --user "$(id -u):$(id -g)" \
  --entrypoint Rscript "${image_name}" /opt/foreststructure/run.R --help > "${help_path}"
for expected_help in \
  "Optional GeoJSON or GeoPackage Audit AOI" \
  "omitted, tiles cover the complete" \
  "--dataset-id" \
  "Scientific parameters" \
  "Runtime controls" \
  "--threads" \
  "--catalog-workers" \
  "default: 2" \
  "--tile-size" \
  "default: 20" \
  "--grid-search-step" \
  "default: 0.5" \
  "--dtm-chunk-size" \
  "default: 60" \
  "--target-chunk-points" \
  "default: 50000000" \
  "--instance-dimension" \
  "--segment-diagnostics"; do
  if ! grep -Fq -- "${expected_help}" "${help_path}"; then
    echo "CLI help is missing '${expected_help}'; see ${help_path}" >&2
    exit 1
  fi
done

run_failure_case() {
  local name="$1"
  local expected_message="$2"
  local point_cloud_path="$3"
  local aoi_path="$4"
  local output_mode="$5"
  local dataset_id="$6"
  shift 6
  local output_dir="${results_dir}/failure_${name}"
  local log_path="${results_dir}/failure_${name}.log"
  mkdir -p "${output_dir}"

  if docker run --rm --network none \
    "${docker_limits[@]}" \
    --user "$(id -u):$(id -g)" \
    --volume "${repo_dir}/tests/fixtures:/fixtures:ro" \
    --volume "${input_dir}:/in:ro" \
    --volume "${output_dir}:/out:${output_mode}" \
    --entrypoint Rscript \
    "${image_name}" \
    /opt/foreststructure/run.R \
    --point-cloud "${point_cloud_path}" \
    --aoi "${aoi_path}" \
    --dataset-id "${dataset_id}" \
    --output-dir /out \
    "$@" > "${log_path}" 2>&1; then
    echo "failure case ${name} unexpectedly succeeded" >&2
    return 1
  fi
  if ! grep -Fq -- "${expected_message}" "${log_path}"; then
    echo "Failure case ${name} did not report '${expected_message}'; see ${log_path}" >&2
    tail -n 20 "${log_path}" >&2
    return 1
  fi
}

run_failure_case missing_point "exactly one existing LAS/LAZ file" \
  /in/missing.laz /fixtures/aoi.geojson rw 150
run_failure_case point_directory "exactly one existing LAS/LAZ file" \
  /in /fixtures/aoi.geojson rw 150
run_failure_case point_extension "must have a .las or .laz extension" \
  /fixtures/aoi.geojson /fixtures/aoi.geojson rw 150
run_failure_case invalid_dataset_id "--dataset-id must be a positive integer" \
  /in/point_cloud.laz /fixtures/aoi.geojson rw 0
run_failure_case invalid_tile_size "--tile-size must be greater than zero" \
  /in/point_cloud.laz /fixtures/aoi.geojson rw 150 --tile-size 0
run_failure_case invalid_threads "--threads must be zero or greater" \
  /in/point_cloud.laz /fixtures/aoi.geojson rw 150 --threads -1
run_failure_case invalid_catalog_workers "--catalog-workers must be at least one" \
  /in/point_cloud.laz /fixtures/aoi.geojson rw 150 --catalog-workers 0
run_failure_case excessive_catalog_workers "--catalog-workers must not exceed ten" \
  /in/point_cloud.laz /fixtures/aoi.geojson rw 150 --catalog-workers 11
run_failure_case invalid_target_chunk_points "--target-chunk-points must be zero or greater" \
  /in/point_cloud.laz /fixtures/aoi.geojson rw 150 --target-chunk-points -1
run_failure_case malformed_aoi "must contain only Polygon or MultiPolygon" \
  /in/point_cloud.laz /fixtures/aoi_invalid.geojson rw 150
run_failure_case nonoverlap_aoi "does not overlap the point-cloud XY extent" \
  /in/point_cloud.laz /fixtures/aoi_nonoverlap.geojson rw 150
run_failure_case unwritable_output "--output-dir must be writable" \
  /in/point_cloud.laz /fixtures/aoi.geojson ro 150

run_case() {
  local name="$1"
  local aoi_path="$2"
  local expected_tiles="$3"
  local point_cloud_path="$4"
  local expected_instance_dimension="$5"
  shift 5
  local dataset_id="150"
  local output_dir="${results_dir}/${name}"
  local log_path="${output_dir}/container.log"
  local expected_footprint_source="audit_aoi"
  local -a aoi_arguments=()
  if [[ "${aoi_path}" == "NONE" ]]; then
    expected_footprint_source="point_cloud_extent"
  else
    aoi_arguments=(--aoi "${aoi_path}")
  fi
  mkdir -p "${output_dir}"

  docker run --rm --network none \
    "${docker_limits[@]}" \
    --user "$(id -u):$(id -g)" \
    --volume "${repo_dir}/tests/fixtures:/fixtures:ro" \
    --volume "${input_dir}:/in:ro" \
    --volume "${output_dir}:/out" \
    --entrypoint Rscript \
    "${image_name}" \
    /opt/foreststructure/run.R \
    --point-cloud "${point_cloud_path}" \
    "${aoi_arguments[@]}" \
    --dataset-id "${dataset_id}" \
    --output-dir /out \
    --threads "${FORESTSTRUCTURE_CPUS:-10}" \
    "$@" 2>&1 | tee "${log_path}"

  python "${repo_dir}/tests/verify_container_output.py" \
    "${output_dir}" \
    "${expected_tiles}" \
    "${expected_instance_dimension}" \
    "${expected_footprint_source}" \
    "${dataset_id}"

  if [[ "${expected_instance_dimension}" == "NA" ]]; then
    grep -q "No configured Instance Dimensions found" "${log_path}"
  else
    grep -q "Using Instance Dimensions: ${expected_instance_dimension//|/, }" "${log_path}"
  fi
  if [[ "${expected_footprint_source}" == "point_cloud_extent" ]]; then
    grep -q "covering the complete point-cloud XY extent" "${log_path}"
  fi
}

run_case whole_cloud NONE 9 /in/point_cloud.laz NA --performance-report
run_case inclusion_only /fixtures/aoi.geojson 4 /in/point_cloud.laz NA
run_case geojson_exclusion /fixtures/aoi_with_exclusion.geojson 3 /in/point_cloud.laz NA
run_case gpkg_exclusion /in/aoi_with_exclusion.gpkg 3 /in/point_cloud.laz NA
run_case zero_tiles /fixtures/aoi_zero_tiles.geojson 0 /in/point_cloud.laz NA
run_case segmented /fixtures/aoi.geojson 4 /in/point_cloud_segmented.laz \
  'PredInstance|PredInstance_SAT'
run_case aliased /fixtures/aoi.geojson 4 /in/point_cloud_segmented.laz \
  'TreeAlias|PredInstance' \
  --instance-dimension MissingAlias \
  --instance-dimension TreeAlias \
  --instance-dimension PredInstance \
  --segment-diagnostics \
  --performance-report

docker run --rm --network none \
  "${docker_limits[@]}" \
  --user "$(id -u):$(id -g)" \
  --entrypoint Rscript \
  --volume "${repo_dir}/tests:/tests:ro" \
  --volume "${input_dir}:/inputs:ro" \
  --volume "${results_dir}:/results:ro" \
  "${image_name}" \
  /tests/verify_rasters.R /results /inputs/point_cloud.laz

# Also test the production entry point that is published by the default image.
FORESTSTRUCTURE_SKIP_BUILD=1 FORESTSTRUCTURE_JULIA_IMAGE="${image_name}" \
  bash "${repo_dir}/tests/test_julia_memory_safe.sh"
