#!/usr/bin/env bash
# Historical 70 m-policy oracle replay only: select the matching pinned image
# (and historical source override, if any); do not rebuild it from current sources.
# The current observed_voxel_height_v1 policy intentionally changes metric values/schema.
# Abort on command errors/unset variables and propagate failures through ``tee``.
set -euo pipefail

# Resolve repository paths independently of the caller's current directory.
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# Read the dataset and its three required evidence inputs from positional args.
dataset_id="${1:-1338}"
point_cloud="${2:-}"
aoi_geojson="${3:-}"
oracle_dir="${4:-}"
# Default to Julia's original instance field while allowing SAT/FM validation.
dimension="${5:-PredInstance}"
# The ordered LAS/LAZ companion is used only for identity/provenance validation.
original_point_cloud="${6:-}"
# Fail before creating work if any required source/oracle path is absent.
if [[ ! -f "${point_cloud}" || ! -f "${original_point_cloud}" ||
      ! -f "${aoi_geojson}" || ! -d "${oracle_dir}" ]]; then
  echo "usage: $0 DATASET_ID INPUT.copc.laz AOI.geojson VALID_UPDATED_DIR [DIMENSION] ORIGINAL.laz" \
    >&2
  exit 2
fi
# This regression specifically tests spatial COPC reads and order restoration.
if [[ "${point_cloud,,}" != *.copc.laz ]]; then
  echo "alignment regression requires a COPC input" >&2
  exit 2
fi

# Canonical paths make Docker mounts and provenance assertions deterministic.
point_cloud="$(realpath "${point_cloud}")"
original_point_cloud="$(realpath "${original_point_cloud}")"
aoi_geojson="$(realpath "${aoi_geojson}")"
oracle_dir="$(realpath "${oracle_dir}")"
# Allow a locally built candidate image without changing the test script.
image_name="${FORESTSTRUCTURE_JULIA_IMAGE:-3dtrees-foreststructure:copc-local}"
# Optional source override supports testing edited R code without rebuilding.
source_mount="${FORESTSTRUCTURE_SOURCE_MOUNT:-}"
# Keep the acceptance run within the production-like resource envelope.
test_cpus="${FORESTSTRUCTURE_TEST_CPUS:-10}"
container_memory_gib="${FORESTSTRUCTURE_TEST_MEMORY_GIB:-30}"
analysis_memory_gib="${FORESTSTRUCTURE_TEST_ANALYSIS_MEMORY_GIB:-25}"
# Put durable validation evidence on the large-data filesystem by default.
validation_parent="${VALID_UPDATED_VALIDATION_ROOT:-/mnt/ssds/kg281/_foreststructure_valid_updated_validation}"
mkdir -p "${validation_parent}"
# A caller-selected root makes failed-dataset reruns reproducible and inspectable.
if [[ -n "${VALID_UPDATED_RUN_ROOT:-}" ]]; then
  run_root="$(realpath -m "${VALID_UPDATED_RUN_ROOT}")"
  mkdir -p "${run_root}"
else
  run_root="$(mktemp -d "${validation_parent}/dataset${dataset_id}_XXXXXXXX")"
fi
# Keep candidate products and temporary processing state physically separate.
candidate_dir="${run_root}/candidate"
candidate_work="${run_root}/work"
mkdir -p "${candidate_dir}" "${candidate_work}"

# Use a unique name so concurrent dataset validations cannot collide.
container_name="foreststructure_valid_updated_${dataset_id}_$$_${RANDOM}"
# Always remove an interrupted candidate container while preserving its files.
cleanup_container() {
  docker rm -f "${container_name}" >/dev/null 2>&1 || true
}
# Mount edited sources read-only only when explicitly requested by the caller.
source_mount_args=()
if [[ -n "${source_mount}" ]]; then
  source_mount="$(realpath "${source_mount}")"
  source_mount_args=(--volume "${source_mount}:/opt/foreststructure:ro")
fi
# Register cleanup before starting Docker so signals cannot leave stale workers.
trap cleanup_container EXIT
# Capture Docker's status manually because the log must also stream through tee.
set +e
# Run only the canonical COPC as scientific input.  The original LAZ is mounted
# separately and passed as identity evidence; the image never computes from it.
docker run \
  --name "${container_name}" \
  --network none \
  --cpus "${test_cpus}" \
  --memory "${container_memory_gib}g" \
  --memory-swap "${container_memory_gib}g" \
  --user "$(id -u):$(id -g)" \
  --env "FORESTSTRUCTURE_THREADS=${test_cpus}" \
  --env "FORESTSTRUCTURE_CATALOG_WORKERS=${FORESTSTRUCTURE_CATALOG_WORKERS:-1}" \
  --volume "${point_cloud}:/in/input.copc.laz:ro" \
  --volume "${original_point_cloud}:/in/original.laz:ro" \
  --volume "${aoi_geojson}:/in/aoi.geojson:ro" \
  --volume "${candidate_dir}:/out" \
  --volume "${candidate_work}:/work" \
  "${source_mount_args[@]}" \
  "${image_name}" \
  --point-cloud /in/input.copc.laz \
  --original-point-cloud /in/original.laz \
  --aoi-geojson /in/aoi.geojson \
  --dataset-id "${dataset_id}" \
  --output-dir /out \
  --temp-dir /work \
  --memory-budget-gib "${analysis_memory_gib}" \
  --instance-dimension "${dimension}" \
  2>&1 | tee "${run_root}/candidate.log"
docker_status="${PIPESTATUS[0]}"
# Restore fail-fast behavior immediately after capturing the pipeline status.
set -e
# Do not compare partial artifacts from a failed ForestStructure execution.
if [[ "${docker_status}" != "0" ]]; then
  echo "candidate failed before valid_updated comparison; run_root=${run_root}" >&2
  exit "${docker_status}"
fi
# Successful runs are removed explicitly, then the failure cleanup is disarmed.
docker rm "${container_name}" >/dev/null
trap - EXIT

# Prove from machine-readable provenance that computation actually used COPC,
# restored OriginalPointIndex, and validated the mounted ordered companion.
jq -e '
  .requested_is_copc == true and
  .source_is_copc == true and
  .companion_header_identity_validated == true and
  .point_order_dimension == "OriginalPointIndex" and
  .catalog_spatial_streaming == "COPC" and
  (.requested_point_cloud | endswith("input.copc.laz")) and
  (.source_point_cloud | endswith("input.copc.laz")) and
  (.original_point_cloud | endswith("original.laz"))
' "${candidate_dir}/${dataset_id}_julia_memory_safe_run.json" >/dev/null

# Compare result rows, segment diagnostics, and GeoJSON features at zero tolerance.
/usr/bin/python3 "${repo_dir}/validation/compare_valid_updated_science.py" \
  "${dataset_id}" "${dimension}" "${oracle_dir}" "${candidate_dir}" \
  "${run_root}/science_comparison.json"

# Compare both scientific rasters cell-for-cell and require identical geometry.
for raster in dtm chm; do
  # Use the candidate image's pinned R/terra stack for a fair raster comparison.
  docker run --rm --network none \
    --cpus "${test_cpus}" \
    --memory "${container_memory_gib}g" \
    --memory-swap "${container_memory_gib}g" \
    --user "$(id -u):$(id -g)" \
    --entrypoint Rscript \
    --volume "${repo_dir}/validation:/validators:ro" \
    --volume "${oracle_dir}/${dataset_id}_${raster}.tif:/validation/expected.tif:ro" \
    --volume "${candidate_dir}/${dataset_id}_${raster}.tif:/validation/actual.tif:ro" \
    --volume "${run_root}:/report" \
    "${image_name}" /validators/compare_valid_updated_raster.R \
    "${raster}" /validation/expected.tif /validation/actual.tif \
    "/report/${raster}_comparison.json"
done

# Emit the evidence root as the single success handoff for humans and CI logs.
echo "COPC/valid_updated alignment passed: run_root=${run_root}"
