#!/usr/bin/env bash
# Historical 70 m-policy oracle replay only: select the matching pinned image
# (and historical source override, if any); do not rebuild it from current sources.
# The current observed_voxel_height_v1 policy intentionally changes metric values/schema.
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
point_cloud="${1:-${DATASET150_LAZ:-}}"
authoritative_gpkg="${2:-${DATASET150_GPKG:-}}"
julia_script="${3:-${JULIA_ORIGINAL_SCRIPT:-}}"
if [[ ! -f "${point_cloud}" || ! -f "${authoritative_gpkg}" || ! -f "${julia_script}" ]]; then
  echo "usage: $0 ORIGINAL_80_LAZ AUTHORITATIVE_80_GPKG JULIA_SCRIPT" >&2
  exit 2
fi
point_cloud="$(realpath "${point_cloud}")"
authoritative_gpkg="$(realpath "${authoritative_gpkg}")"
julia_script="$(realpath "${julia_script}")"

assert_sha256() {
  local path="$1"
  local expected="$2"
  local actual
  actual="$(sha256sum "${path}" | cut -d ' ' -f 1)"
  if [[ "${actual}" != "${expected}" ]]; then
    echo "unexpected SHA-256 for ${path}: ${actual}" >&2
    exit 2
  fi
}
assert_sha256 "${point_cloud}" d767997ce5868a98fb9dc0a688cb3d300fc8b7cee45571eadcf9b6102c9f9789
assert_sha256 "${authoritative_gpkg}" cc7befe0c55ce5e09c66c4eda062cbc947937ea005e7cc80161a028685f5c643
assert_sha256 "${julia_script}" 746d57b4c937001af31e4ccd1b9f14edb5cebb15d46154ad9e20d0ce39f78226

image_name="${FORESTSTRUCTURE_JULIA_IMAGE:-3dtrees-foreststructure:copc-local}"
catalog_workers="${FORESTSTRUCTURE_CATALOG_WORKERS:-2}"
test_cpus="${FORESTSTRUCTURE_TEST_CPUS:-15}"
validation_parent="${DATASET150_VALIDATION_ROOT:-/mnt/ssds/kg281/_foreststructure_julia_memory_safe_validation}"
mkdir -p "${validation_parent}"
run_root="$(mktemp -d "${validation_parent}/dataset150_XXXXXXXX")"
canonical_dir="${run_root}/converted"
candidate_dir="${run_root}/memory_safe"
candidate_work="${run_root}/work"
oracle_root="${run_root}/julia_original"
oracle_data="${oracle_root}/C:/Users/Julia Gäßler/Documents/Studium/Master/Masterarbeit/Indices/Data/Test2/ULS/Test"
oracle_output="${oracle_root}/C:/Users/Julia Gäßler/Documents/Studium/Master/Masterarbeit/Indices/Data/Test2"
mkdir -p "${canonical_dir}" "${candidate_dir}" "${candidate_work}" "${oracle_data}"

if [[ "${FORESTSTRUCTURE_SKIP_BUILD:-0}" != "1" ]]; then
  docker build --file "${repo_dir}/Dockerfile.julia-memory-safe" --build-arg "VCS_REF=$(git -C "${repo_dir}" rev-parse HEAD)" --build-arg "IMAGE_VERSION=julia-memory-safe-dataset150" --tag "${image_name}" "${repo_dir}"
fi

docker run --rm --network none --cpus "${test_cpus}" --memory 75g --memory-swap 75g --user "$(id -u):$(id -g)" --entrypoint Rscript --volume "${repo_dir}/tests:/tests:ro" --volume "${authoritative_gpkg}:/in/80.gpkg:ro" --volume "${canonical_dir}:/converted" "${image_name}" /tests/create_dataset150_canonical_geojson.R /in/80.gpkg /converted/80.geojson

docker run --rm --network none --cpus "${test_cpus}" --memory 75g --memory-swap 75g --user "$(id -u):$(id -g)" --entrypoint Rscript --volume "${point_cloud}:/in/80.laz:ro" --volume "${canonical_dir}:/converted" "${image_name}" /opt/foreststructure/convert_aoi.R --point-cloud /in/80.laz --aoi-geojson /converted/80.geojson --output-gpkg /converted/80.gpkg --provenance-json /converted/conversion.json

docker run --rm --network none --cpus "${test_cpus}" --memory 75g --memory-swap 75g --user "$(id -u):$(id -g)" --entrypoint bash --volume "${point_cloud}:/in/80.laz:ro" --volume "${canonical_dir}/80.gpkg:/converted/80.gpkg:ro" --volume "${julia_script}:/oracle/Indices_Final_run.R:ro" --volume "${repo_dir}/tests:/tests:ro" --volume "${oracle_root}:/validation" "${image_name}" -c '
  set -euo pipefail
  target="/validation/C:/Users/Julia Gäßler/Documents/Studium/Master/Masterarbeit/Indices/Data/Test2/ULS/Test"
  ln -s /in/80.laz "$target/80.laz"
  cp /converted/80.gpkg "$target/80.gpkg"
  cd /validation
  Rscript /tests/run_julia_original_with_dtm_capture.R \
    /oracle/Indices_Final_run.R /validation/julia_dtm.tif
' 2>&1 | tee "${run_root}/julia_original.log"

container_name="foreststructure_jms_150_$$_${RANDOM}"
cleanup_container() {
  docker rm -f "${container_name}" >/dev/null 2>&1 || true
}
trap cleanup_container EXIT
set +e
docker run --name "${container_name}" --network none --cpus "${test_cpus}" --memory 75g --memory-swap 75g --user "$(id -u):$(id -g)" --env "FORESTSTRUCTURE_THREADS=${test_cpus}" --env "FORESTSTRUCTURE_CATALOG_WORKERS=${catalog_workers}" --volume "${point_cloud}:/in/80.laz:ro" --volume "${canonical_dir}/80.geojson:/in/80.geojson:ro" --volume "${candidate_dir}:/out" --volume "${candidate_work}:/work" "${image_name}" --point-cloud /in/80.laz --aoi-geojson /in/80.geojson --dataset-id 150 --output-dir /out --temp-dir /work --memory-budget-gib 60 --sensor ULS --country Test --instance-dimension PredInstance 2>&1 | tee "${run_root}/memory_safe.log"
docker_status="${PIPESTATUS[0]}"
set -e
oom_killed="$(docker inspect "${container_name}" --format '{{.State.OOMKilled}}')"
exit_code="$(docker inspect "${container_name}" --format '{{.State.ExitCode}}')"
docker rm "${container_name}" >/dev/null
trap - EXIT
if [[ "${docker_status}" != "0" || "${exit_code}" != "0" || "${oom_killed}" != "false" ]]; then
  echo "candidate exited ${exit_code}; OOMKilled=${oom_killed}" >&2
  exit 1
fi

grep -q 'SELECTOR phase=dtm selector=xyz' "${run_root}/memory_safe.log"
grep -q 'SELECTOR phase=tile selector=xyz' "${run_root}/memory_safe.log"
grep -q 'dimension=PredInstance ordinal=1 selector=xyz1' "${run_root}/memory_safe.log"
grep -q "Using validated Julia-memory-safe catalog worker count: ${catalog_workers}" "${run_root}/memory_safe.log"
if [[ "$(grep -c "Processing .* LAScatalog chunks with ${catalog_workers} worker(s)" "${run_root}/memory_safe.log")" -lt 3 ]]; then
  echo "DTM, CHM, and segmentation did not all use ${catalog_workers} catalog worker(s)" >&2
  exit 1
fi
if grep -Eq 'selector=(xyz0|[^ ]*[rn*])' "${run_root}/memory_safe.log"; then
  echo "forbidden LAS selector found in memory-safe log" >&2
  exit 1
fi
if grep -Eqi 'coordinate.*(rescal|moved)|copc' "${run_root}/memory_safe.log"; then
  echo "coordinate rescaling or COPC evidence found in memory-safe log" >&2
  exit 1
fi
if [[ "$(jq -r '.catalog_workers' "${candidate_dir}/150_julia_memory_safe_run.json")" != "${catalog_workers}" ]]; then
  echo "candidate provenance does not record ${catalog_workers} catalog worker(s)" >&2
  exit 1
fi
shared_geometry_hash="$(jq -r '.normalized_xy_sha256' "${canonical_dir}/conversion.json")"
candidate_geometry_hash="$(jq -r '.normalized_xy_sha256' "${candidate_dir}/150_aoi_conversion.json")"
candidate_tile_count="$(jq -r '.tile_count' "${candidate_dir}/150_aoi_conversion.json")"
if [[ "${shared_geometry_hash}" != "${candidate_geometry_hash}" || "${candidate_tile_count}" != "147" ]]; then
  echo "candidate conversion differs from the Julia runtime conversion" >&2
  exit 1
fi

python3 "${repo_dir}/tests/compare_julia_original_dataset150.py" "${oracle_output}" "${candidate_dir}" "${run_root}/comparison.json"
docker run --rm --network none --user "$(id -u):$(id -g)" \
  --entrypoint Rscript \
  --volume "${repo_dir}/tests:/tests:ro" \
  --volume "${oracle_root}/julia_dtm.tif:/validation/julia_dtm.tif:ro" \
  --volume "${candidate_dir}/150_dtm.tif:/validation/candidate_dtm.tif:ro" \
  --volume "${run_root}:/report" \
  "${image_name}" /tests/compare_julia_dtm.R \
  /validation/julia_dtm.tif /validation/candidate_dtm.tif \
  /report/dtm_comparison.json
{
  echo "run_root=${run_root}"
  echo "image=${image_name}"
  echo "image_id=$(docker image inspect "${image_name}" --format '{{.Id}}')"
  echo "base_image_id=$(docker image inspect 3dtrees-foreststructure:v0.1.2 --format '{{.Id}}')"
  echo "catalog_workers=${catalog_workers}"
  echo "test_cpus=${test_cpus}"
  echo "point_cloud_sha256=$(sha256sum "${point_cloud}" | cut -d ' ' -f 1)"
  echo "canonical_geojson_sha256=$(sha256sum "${canonical_dir}/80.geojson" | cut -d ' ' -f 1)"
  echo "runtime_gpkg_sha256=$(sha256sum "${canonical_dir}/80.gpkg" | cut -d ' ' -f 1)"
  echo "julia_script_sha256=$(sha256sum "${julia_script}" | cut -d ' ' -f 1)"
  echo "oom_killed=${oom_killed}"
} > "${run_root}/validation.env"

echo "dataset-150 Julia equivalence validation passed: ${run_root}"
