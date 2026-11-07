#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
image_name="${RASTER_IMAGE:-3dtrees-raster:test}"
docker_limits=(--cpus "${RASTER_CPUS:-10}" --memory "${RASTER_MEMORY:-70g}" --memory-swap "${RASTER_MEMORY:-70g}")
input_dir="$(mktemp -d)"
results_dir="$(mktemp -d)"
cleanup() {
  local status="$?"
  if [[ "${status}" != 0 || "${RASTER_KEEP_TEST_ROOT:-0}" == 1 ]]; then
    echo "Container gate inputs: ${input_dir}"
    echo "Container gate results and logs: ${results_dir} (exit ${status})"
  else
    rm -rf "${input_dir}" "${results_dir}"
  fi
}
trap cleanup EXIT

if [[ "${RASTER_SKIP_BUILD:-0}" != "1" ]]; then
  docker build --tag "${image_name}" "${repo_dir}"
fi
# Exercise the optional observational export in the standard publish gate too.
docker run --rm --network none "${docker_limits[@]}" \
  --env PYTHONDONTWRITEBYTECODE=1 --entrypoint python3 \
  --volume "${repo_dir}:/workspace:ro" "${image_name}" \
  /workspace/tests/test_voxel_counts.py
python3 "${repo_dir}/tests/test_output_artifact_validation.py"
docker run --rm --network none "${docker_limits[@]}" \
  --env PYTHONDONTWRITEBYTECODE=1 --entrypoint python3 \
  --volume "${repo_dir}:/workspace:ro" "${image_name}" \
  /workspace/tests/test_denoising_validation.py
python3 "${repo_dir}/tests/test_copc_task_counts.py"
docker run --rm --network none \
  "${docker_limits[@]}" \
  --entrypoint Rscript \
  --volume "${repo_dir}:/workspace:ro" \
  --workdir /workspace/src \
  "${image_name}" \
  -e 'Sys.setenv(RASTER_SOURCE_ONLY = "1"); source("analysis.R"); source("denoising.R"); source("../tests/test_runtime_helpers.R"); source("../tests/test_catalog_scheduler.R"); source("../tests/test_resource_scheduler.R"); source("../tests/test_resource_memory.R"); source("../tests/test_denoising.R"); source("../tests/test_height_policy.R")'
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
  -e 'Sys.setenv(RASTER_SOURCE_ONLY = "1"); source("analysis.R"); source("../tests/test_failed_dataset_regressions.R")'
docker run --rm --network none \
  "${docker_limits[@]}" \
  --entrypoint Rscript \
  --volume "${repo_dir}:/workspace:ro" \
  --workdir /workspace/src \
  "${image_name}" \
  -e 'Sys.setenv(RASTER_SOURCE_ONLY = "1"); source("raster.R"); source("../tests/test_raster.R"); source("../tests/test_instance_read_groups.R"); source("../tests/test_resource_validation.R")'
docker run --rm --network none \
  "${docker_limits[@]}" \
  --entrypoint Rscript \
  --volume "${repo_dir}:/workspace:ro" \
  --workdir /workspace/src \
  "${image_name}" \
  -e 'Sys.setenv(RASTER_SOURCE_ONLY = "1"); source("analysis.R"); source("../tests/test_adaptive_chunks.R"); source("../tests/test_chm_adaptive_edges.R"); source("../tests/test_segment_catalog_boundaries.R")'
docker run --rm --network none \
  "${docker_limits[@]}" \
  --entrypoint Rscript \
  --volume "${repo_dir}:/workspace:ro" \
  --workdir /workspace/src \
  "${image_name}" \
  -e 'Sys.setenv(RASTER_SOURCE_ONLY = "1"); source("analysis.R"); source("../tests/test_dtm_coordinate_scale.R")'
docker run --rm --network none \
  "${docker_limits[@]}" \
  --entrypoint Rscript \
  --volume "${repo_dir}:/workspace:ro" \
  --workdir /workspace/src \
  "${image_name}" \
  -e 'Sys.setenv(RASTER_SOURCE_ONLY = "1"); source("analysis.R"); source("../tests/test_raster_merge_semantics.R"); source("../tests/test_surface_rasters.R")'
docker run --rm --network none \
  "${docker_limits[@]}" \
  --entrypoint Rscript \
  --volume "${repo_dir}:/workspace:ro" \
  --workdir /workspace/src \
  "${image_name}" \
  -e 'Sys.setenv(RASTER_SOURCE_ONLY = "1"); source("analysis.R"); source("aoi_conversion.R"); source("../tests/test_aoi_conversion.R")'
galaxy_command="$(docker run --rm --network none "${docker_limits[@]}" \
  --user "$(id -u):$(id -g)" "${image_name}" /bin/bash -c 'echo galaxy-command-ok')"
if [[ "${galaxy_command}" != "galaxy-command-ok" ]]; then
  echo "Image did not run an explicit command: ${galaxy_command}" >&2
  exit 1
fi
docker run --rm --network none "${docker_limits[@]}" --user "$(id -u):$(id -g)" \
  "${image_name}" --help > "${results_dir}/runner_help.txt"
if ! grep -Fq -- "--instance-dimension" "${results_dir}/runner_help.txt"; then
  echo "Flag-only run did not reach the production runner" >&2
  exit 1
fi

# Also test the production entry point that is published by the default image.
RASTER_SKIP_BUILD=1 RASTER_IMAGE="${image_name}" \
  bash "${repo_dir}/tests/test_raster.sh"
