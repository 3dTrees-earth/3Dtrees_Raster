#!/usr/bin/env bash
# Historical 70 m-policy oracle replay only: select the matching pinned image
# (and historical source override, if any); do not rebuild it from current sources.
# The current observed_voxel_height_v1 policy intentionally changes metric values/schema.
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ "$#" != 5 ]]; then
  echo "usage: $0 DATASET_ID ORIGINAL_LAZ AOI_GEOJSON SENSOR JULIA_SCRIPT" >&2
  exit 2
fi

dataset_id="$1"
point_cloud="$(realpath "$2")"
aoi_geojson="$(realpath "$3")"
sensor="$4"
julia_script="$(realpath "$5")"
country="Validation"

if [[ ! "${dataset_id}" =~ ^[1-9][0-9]*$ ]]; then
  echo "DATASET_ID must be a positive integer" >&2
  exit 2
fi
if [[ ! -f "${point_cloud}" || ! -f "${aoi_geojson}" || ! -f "${julia_script}" ]]; then
  echo "LAZ, AOI GeoJSON, or Julia script does not exist" >&2
  exit 2
fi
if [[ ! "${sensor}" =~ ^(TLS|MLS|ULS)$ ]]; then
  echo "SENSOR must be TLS, MLS, or ULS" >&2
  exit 2
fi
if [[ "${point_cloud,,}" == *.copc.laz ]]; then
  echo "the Julia-faithful canary requires an original LAZ, not COPC" >&2
  exit 2
fi

julia_script_sha256="$(sha256sum "${julia_script}" | cut -d ' ' -f 1)"
if [[ "${julia_script_sha256}" != "746d57b4c937001af31e4ccd1b9f14edb5cebb15d46154ad9e20d0ce39f78226" ]]; then
  echo "the Julia source-of-truth script has an unexpected SHA-256" >&2
  exit 2
fi

image_name="${FORESTSTRUCTURE_JULIA_IMAGE:-3dtrees-foreststructure:copc-local}"
validation_parent="${FORESTSTRUCTURE_CANARY_ROOT:-/mnt/ssds/kg281/_foreststructure_julia_memory_safe_validation/medium_large}"
test_cpus="${FORESTSTRUCTURE_TEST_CPUS:-15}"
catalog_workers="${FORESTSTRUCTURE_CATALOG_WORKERS:-2}"
reuse_oracle_root="${FORESTSTRUCTURE_REUSE_ORACLE_ROOT:-}"
if [[ ! "${test_cpus}" =~ ^[1-9][0-9]*$ || ! "${catalog_workers}" =~ ^[1-9][0-9]*$ ]]; then
  echo "FORESTSTRUCTURE_TEST_CPUS and FORESTSTRUCTURE_CATALOG_WORKERS must be positive integers" >&2
  exit 2
fi
mkdir -p "${validation_parent}"
run_root="$(mktemp -d "${validation_parent}/dataset${dataset_id}_XXXXXXXX")"
converted_dir="${run_root}/converted"
candidate_dir="${run_root}/memory_safe"
candidate_work="${run_root}/work"
oracle_root="${run_root}/julia_original"
oracle_reused="false"
if [[ -n "${reuse_oracle_root}" ]]; then
  oracle_root="$(realpath "${reuse_oracle_root}")"
  oracle_reused="true"
fi
oracle_output="${oracle_root}/C:/Users/Julia Gäßler/Documents/Studium/Master/Masterarbeit/Indices/Data/Test2"
oracle_data="${oracle_output}/${sensor}/${country}"
mkdir -p "${converted_dir}" "${candidate_dir}" "${candidate_work}"
if [[ "${oracle_reused}" == "false" ]]; then
  mkdir -p "${oracle_data}"
fi

point_cloud_name="$(basename "${point_cloud}")"
point_cloud_stem="${point_cloud_name%.*}"
active_container=""
current_step="preflight"
status="failed"
oracle_exit_code="not_run"
oracle_oom_killed="not_run"
candidate_exit_code="not_run"
candidate_oom_killed="not_run"

finish() {
  local shell_status=$?
  if [[ -n "${active_container}" ]]; then
    docker rm -f "${active_container}" >/dev/null 2>&1 || true
  fi
  {
    echo "status=${status}"
    echo "shell_exit_code=${shell_status}"
    echo "last_step=${current_step}"
    echo "dataset_id=${dataset_id}"
    echo "run_root=${run_root}"
    echo "image=${image_name}"
    echo "image_id=$(docker image inspect "${image_name}" --format '{{.Id}}' 2>/dev/null || echo unavailable)"
    echo "point_cloud=${point_cloud}"
    echo "point_cloud_sha256=$(sha256sum "${point_cloud}" | cut -d ' ' -f 1)"
    echo "aoi_geojson=${aoi_geojson}"
    echo "aoi_geojson_sha256=$(sha256sum "${aoi_geojson}" | cut -d ' ' -f 1)"
    echo "julia_script_sha256=${julia_script_sha256}"
    echo "sensor=${sensor}"
    echo "country=${country}"
    echo "test_cpus=${test_cpus}"
    echo "catalog_workers=${catalog_workers}"
    echo "oracle_root=${oracle_root}"
    echo "oracle_reused=${oracle_reused}"
    echo "oracle_exit_code=${oracle_exit_code}"
    echo "oracle_oom_killed=${oracle_oom_killed}"
    echo "candidate_exit_code=${candidate_exit_code}"
    echo "candidate_oom_killed=${candidate_oom_killed}"
  } > "${run_root}/validation.env"
  return "${shell_status}"
}
trap finish EXIT

echo "CANARY dataset=${dataset_id} run_root=${run_root}"

current_step="aoi_conversion"
docker run --rm --network none --cpus "${test_cpus}" \
  --memory 75g --memory-swap 75g --user "$(id -u):$(id -g)" \
  --env "FORESTSTRUCTURE_THREADS=${test_cpus}" \
  --entrypoint Rscript \
  --volume "${point_cloud}:/in/${point_cloud_name}:ro" \
  --volume "${aoi_geojson}:/in/aoi.geojson:ro" \
  --volume "${converted_dir}:/converted" \
  "${image_name}" /opt/foreststructure/convert_aoi.R \
  --point-cloud "/in/${point_cloud_name}" \
  --aoi-geojson /in/aoi.geojson \
  --output-gpkg "/converted/${point_cloud_stem}.gpkg" \
  --provenance-json /converted/conversion.json
jq -e '.tile_count > 0' "${converted_dir}/conversion.json" >/dev/null

current_step="julia_original"
if [[ "${oracle_reused}" == "true" ]]; then
  oracle_validation="$(dirname "${oracle_root}")/validation.env"
  current_sha256="$(sha256sum "${point_cloud}" | cut -d ' ' -f 1)"
  current_aoi_sha256="$(sha256sum "${aoi_geojson}" | cut -d ' ' -f 1)"
  if [[ ! -f "${oracle_root}/julia_dtm.tif" || \
        ! -f "${oracle_output}/results.csv" || \
        ! -f "${oracle_output}/segment_diagnostics.csv" || \
        ! -f "${oracle_validation}" ]]; then
    echo "reusable Julia oracle is incomplete" >&2
    exit 1
  fi
  grep -qx "point_cloud_sha256=${current_sha256}" "${oracle_validation}"
  grep -qx "aoi_geojson_sha256=${current_aoi_sha256}" "${oracle_validation}"
  grep -qx "julia_script_sha256=${julia_script_sha256}" "${oracle_validation}"
  grep -qx "oracle_exit_code=0" "${oracle_validation}"
  grep -qx "oracle_oom_killed=false" "${oracle_validation}"
  oracle_exit_code="0"
  oracle_oom_killed="false"
  echo "Reusing SHA-256-matched Julia oracle: ${oracle_root}" | tee "${run_root}/julia_original.log"
else
  active_container="foreststructure_julia_${dataset_id}_$$_${RANDOM}"
  set +e
  docker run --name "${active_container}" --network none --cpus "${test_cpus}" \
    --memory 75g --memory-swap 75g --user "$(id -u):$(id -g)" \
    --entrypoint bash \
    --env "ORACLE_SENSOR=${sensor}" \
    --env "POINT_CLOUD_NAME=${point_cloud_name}" \
    --env "POINT_CLOUD_STEM=${point_cloud_stem}" \
    --volume "${point_cloud}:/in/${point_cloud_name}:ro" \
    --volume "${converted_dir}/${point_cloud_stem}.gpkg:/converted/${point_cloud_stem}.gpkg:ro" \
    --volume "${julia_script}:/oracle/Indices_Final_run.R:ro" \
    --volume "${repo_dir}/tests:/tests:ro" \
    --volume "${oracle_root}:/validation" \
    "${image_name}" -c '
      set -euo pipefail
      target="/validation/C:/Users/Julia Gäßler/Documents/Studium/Master/Masterarbeit/Indices/Data/Test2/${ORACLE_SENSOR}/Validation"
      ln -s "/in/${POINT_CLOUD_NAME}" "${target}/${POINT_CLOUD_NAME}"
      cp "/converted/${POINT_CLOUD_STEM}.gpkg" "${target}/${POINT_CLOUD_STEM}.gpkg"
      cd /validation
      Rscript /tests/run_julia_original_with_dtm_capture.R \
        /oracle/Indices_Final_run.R /validation/julia_dtm.tif
    ' 2>&1 | tee "${run_root}/julia_original.log"
  oracle_docker_status="${PIPESTATUS[0]}"
  set -e
  oracle_oom_killed="$(docker inspect "${active_container}" --format '{{.State.OOMKilled}}')"
  oracle_exit_code="$(docker inspect "${active_container}" --format '{{.State.ExitCode}}')"
  docker rm "${active_container}" >/dev/null
  active_container=""
  if [[ "${oracle_docker_status}" != "0" || "${oracle_exit_code}" != "0" || "${oracle_oom_killed}" != "false" ]]; then
    echo "Julia original exited ${oracle_exit_code}; OOMKilled=${oracle_oom_killed}" >&2
    exit 1
  fi
fi

current_step="memory_safe_candidate"
active_container="foreststructure_candidate_${dataset_id}_$$_${RANDOM}"
set +e
docker run --name "${active_container}" --network none --cpus "${test_cpus}" \
  --memory 75g --memory-swap 75g --user "$(id -u):$(id -g)" \
  --env "FORESTSTRUCTURE_THREADS=${test_cpus}" \
  --env "FORESTSTRUCTURE_CATALOG_WORKERS=${catalog_workers}" \
  --volume "${point_cloud}:/in/${point_cloud_name}:ro" \
  --volume "${aoi_geojson}:/in/aoi.geojson:ro" \
  --volume "${candidate_dir}:/out" \
  --volume "${candidate_work}:/work" \
  "${image_name}" \
  --point-cloud "/in/${point_cloud_name}" \
  --aoi-geojson /in/aoi.geojson \
  --dataset-id "${dataset_id}" \
  --output-dir /out \
  --temp-dir /work \
  --memory-budget-gib 60 \
  --sensor "${sensor}" \
  --country "${country}" \
  --instance-dimension PredInstance \
  2>&1 | tee "${run_root}/memory_safe.log"
candidate_docker_status="${PIPESTATUS[0]}"
set -e
candidate_oom_killed="$(docker inspect "${active_container}" --format '{{.State.OOMKilled}}')"
candidate_exit_code="$(docker inspect "${active_container}" --format '{{.State.ExitCode}}')"
docker rm "${active_container}" >/dev/null
active_container=""
if [[ "${candidate_docker_status}" != "0" || "${candidate_exit_code}" != "0" || "${candidate_oom_killed}" != "false" ]]; then
  echo "memory-safe candidate exited ${candidate_exit_code}; OOMKilled=${candidate_oom_killed}" >&2
  exit 1
fi

current_step="selector_validation"
grep -q 'SELECTOR phase=dtm selector=xyz' "${run_root}/memory_safe.log"
grep -q 'SELECTOR phase=tile selector=xyz' "${run_root}/memory_safe.log"
grep -Eq 'SELECTOR phase=segment dimension=PredInstance ordinal=[1-9] selector=xyz[1-9]' "${run_root}/memory_safe.log"
if grep -E 'SELECTOR phase=(dtm|tile).*selector=.*(rn|\*)' "${run_root}/memory_safe.log"; then
  echo "forbidden LAS selector found in memory-safe log" >&2
  exit 1
fi
if grep -Eqi 'coordinate.*(rescal|moved)|copc' "${run_root}/memory_safe.log"; then
  echo "coordinate rescaling or COPC evidence found in memory-safe log" >&2
  exit 1
fi

current_step="aoi_equivalence"
shared_geometry_hash="$(jq -r '.normalized_xy_sha256' "${converted_dir}/conversion.json")"
candidate_geometry_hash="$(jq -r '.normalized_xy_sha256' "${candidate_dir}/${dataset_id}_aoi_conversion.json")"
shared_tile_count="$(jq -r '.tile_count' "${converted_dir}/conversion.json")"
candidate_tile_count="$(jq -r '.tile_count' "${candidate_dir}/${dataset_id}_aoi_conversion.json")"
if [[ "${shared_geometry_hash}" != "${candidate_geometry_hash}" || "${shared_tile_count}" != "${candidate_tile_count}" ]]; then
  echo "candidate AOI conversion differs from the Julia runtime conversion" >&2
  exit 1
fi

current_step="csv_equivalence"
python3 "${repo_dir}/tests/compare_julia_original.py" \
  "${dataset_id}" "${oracle_output}" "${candidate_dir}" \
  "${run_root}/comparison.json"

current_step="dtm_equivalence"
docker run --rm --network none --user "$(id -u):$(id -g)" \
  --entrypoint Rscript \
  --volume "${repo_dir}/tests:/tests:ro" \
  --volume "${oracle_root}/julia_dtm.tif:/validation/julia_dtm.tif:ro" \
  --volume "${candidate_dir}/${dataset_id}_dtm.tif:/validation/candidate_dtm.tif:ro" \
  --volume "${run_root}:/report" \
  "${image_name}" /tests/compare_julia_dtm.R \
  /validation/julia_dtm.tif /validation/candidate_dtm.tif \
  /report/dtm_comparison.json

current_step="artifact_validation"
jq -e \
  '.status == "valid" and .source_is_copc == false and
   .source_point_count > 0 and
   .catalog_workers == '"${catalog_workers}"' and
   (.processed_instance_dimensions == "PredInstance" or
    .processed_instance_dimensions == ["PredInstance"]) and
   .missing_instance_dimensions == []' \
  "${candidate_dir}/${dataset_id}_julia_memory_safe_run.json" >/dev/null

status="valid"
current_step="complete"
echo "Julia equivalence validation passed for dataset ${dataset_id}: ${run_root}"
