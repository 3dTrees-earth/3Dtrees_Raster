#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
image_name="${FORESTSTRUCTURE_JULIA_IMAGE:-3dtrees-foreststructure:copc-test}"
test_root="$(mktemp -d)"
trap 'rm -rf "${test_root}"' EXIT
mkdir -p "${test_root}/input" "${test_root}/output" "${test_root}/work"

if [[ "${FORESTSTRUCTURE_SKIP_BUILD:-0}" != "1" ]]; then
  docker build --file "${repo_dir}/Dockerfile.julia-memory-safe" --tag "${image_name}" "${repo_dir}"
fi
docker run --rm --network none --env FORESTSTRUCTURE_FIXTURE_EPSG=32632 --user "$(id -u):$(id -g)" --entrypoint Rscript --volume "${repo_dir}/tests:/tests:ro" --volume "${test_root}/input:/generated" "${image_name}" /tests/generate_julia_memory_safe_fixture.R /generated/high_ordinal.laz /generated/aoi.geojson

run_candidate() {
  local bucket_count="${1:-64}"
  local sensor="${2:-ULS}"
  local country="${3:-Test}"
  local memory_budget_gib="${4:-3}"
  docker run --rm --network none --cpus 2 --memory 4g --memory-swap 4g --user "$(id -u):$(id -g)" --env FORESTSTRUCTURE_THREADS=2 --env "FORESTSTRUCTURE_SEGMENT_BUCKET_COUNT=${bucket_count}" --volume "${test_root}/input/high_ordinal.laz:/in/high_ordinal.laz:ro" --volume "${test_root}/input/aoi.geojson:/in/aoi.geojson:ro" --volume "${test_root}/output:/out" --volume "${test_root}/work:/work" "${image_name}" --denoise none --point-cloud /in/high_ordinal.laz --aoi-geojson /in/aoi.geojson --dataset-id 150 --output-dir /out --temp-dir /work --memory-budget-gib "${memory_budget_gib}" --sensor "${sensor}" --country "${country}"
}

run_candidate 2>&1 | tee "${test_root}/first.log"
for dimension in PredInstance PredInstance_SAT PredInstance_FM; do
  grep -q "dimension=${dimension} .*selector=xyz1" "${test_root}/first.log"
  test -s "${test_root}/output/150_${dimension}_results.csv"
  test -s "${test_root}/output/150_${dimension}_segment_diagnostics.csv"
  test -s "${test_root}/output/150_${dimension}_tiles.geojson"
done
cmp "${test_root}/output/150_PredInstance_results.csv" "${test_root}/output/150_PredInstance_SAT_results.csv"
cmp "${test_root}/output/150_PredInstance_results.csv" "${test_root}/output/150_PredInstance_FM_results.csv"
cmp "${test_root}/output/150_PredInstance_segment_diagnostics.csv" "${test_root}/output/150_PredInstance_SAT_segment_diagnostics.csv"
cmp "${test_root}/output/150_PredInstance_segment_diagnostics.csv" "${test_root}/output/150_PredInstance_FM_segment_diagnostics.csv"
jq -e '.processed_instance_dimensions | length == 3' "${test_root}/output/150_julia_memory_safe_run.json" >/dev/null
jq -e '.projections | length == 3' "${test_root}/output/150_julia_memory_safe_run.json" >/dev/null
jq -e '.adaptive_chunking.enabled == false and .adaptive_chunking.target_points == 50000000' "${test_root}/output/150_julia_memory_safe_run.json" >/dev/null
test -z "$(find "${test_root}/work" -mindepth 1 -print -quit)"

# A second run must recompute and replace the validated artifacts, not skip.
# A high bucket count forces tiny disk partitions without changing results.
run_candidate 256 NA NA 70 2>&1 | tee "${test_root}/second.log"
grep -q 'Validated and promoted' "${test_root}/second.log"
jq -e '.segment_bucket_count == 256' "${test_root}/output/150_julia_memory_safe_run.json" >/dev/null
jq -e '.memory_budget_gib == 3.4' "${test_root}/output/150_julia_memory_safe_run.json" >/dev/null

failure_output="${test_root}/failure_output"
failure_work="${test_root}/failure_work"
mkdir -p "${failure_output}" "${failure_work}"
set +e
docker run --rm --network none --cpus 2 --memory 4g --memory-swap 4g --user "$(id -u):$(id -g)" --env FORESTSTRUCTURE_THREADS=2 --volume "${test_root}/input/high_ordinal.laz:/in/high_ordinal.laz:ro" --volume "${test_root}/input/aoi.geojson:/in/aoi.geojson:ro" --volume "${failure_output}:/out" --volume "${failure_work}:/work" "${image_name}" --denoise none --point-cloud /in/high_ordinal.laz --aoi-geojson /in/aoi.geojson --dataset-id 151 --output-dir /out --temp-dir /work --memory-budget-gib 0.0001 --instance-dimension PredInstance > "${test_root}/budget_failure.log" 2>&1
failure_status="$?"
set -e
if [[ "${failure_status}" == "0" ]]; then
  echo "impossible memory budget unexpectedly succeeded" >&2
  exit 1
fi
grep -Eq 'memory_budget_exceeded|single_segment_exceeds_memory_budget' "${test_root}/budget_failure.log"
test -z "$(find "${failure_output}" -mindepth 1 -print -quit)"
test -z "$(find "${failure_work}" -mindepth 1 -print -quit)"

echo "julia-memory-safe synthetic acceptance passed"
