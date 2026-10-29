#!/usr/bin/env bash
# Local synthetic acceptance: default conservative SOR, uncapped heights, and
# invariant science across scheduler/adaptive execution choices.
set -euo pipefail
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source_dir="${FORESTSTRUCTURE_SOURCE_DIR:-${repo}/src}"
baseline="${FORESTSTRUCTURE_BASELINE_IMAGE:-sha256:c4323f2ef0298a192048f92496230edd588d79622c68d6e4847fc065007e227e}"
untwine="${FORESTSTRUCTURE_SMART_TILE_IMAGE:-smarttile:v2.4.0-4de4393}"
baseline="$(docker image inspect --format '{{.Id}}' "${baseline}")"
untwine="$(docker image inspect --format '{{.Id}}' "${untwine}")"
root="$(mktemp -d "${TMPDIR:-/tmp}/foreststructure-denoising.XXXXXX")"
finish() {
  local status="$?"
  if [[ "${status}" != 0 || "${FORESTSTRUCTURE_KEEP_TEST_ROOT:-0}" == 1 ]]; then
    echo "Denoising acceptance evidence: ${root} (exit ${status})"
  else
    rm -rf "${root}"
  fi
}
trap finish EXIT
mkdir -p "${root}/input" "${root}/sources"
cp -a "${source_dir}/." "${root}/sources/"
limits=(--rm --network none --cpus 2 --memory 4g --memory-swap 4g --user "$(id -u):$(id -g)")
reuse="${FORESTSTRUCTURE_DENOISING_REUSE_ROOT:-}"
if [[ -n "${reuse}" ]]; then
  # Reuse only prepared synthetic inputs while iterating candidate sources.
  # Every scientific case is recomputed under the current metric policy.
  python3 - "${reuse}" <<'PY'
import json, pathlib, sys
root = pathlib.Path(sys.argv[1])
assert json.loads((root / 'input/fixture.json').read_text())['fixture_version'] == 2
PY
  cp -a "${reuse}/input/." "${root}/input/"
else
docker run "${limits[@]}" --env FORESTSTRUCTURE_FIXTURE_EPSG=32632 \
  --entrypoint Rscript -v "${repo}/tests:/tests:ro" -v "${root}/input:/data" \
  "${baseline}" /tests/generate_julia_memory_safe_fixture.R \
  /data/patch.laz /data/unused_aoi.geojson distinct-dimensions > "${root}/fixture.log" 2>&1
docker run "${limits[@]}" --entrypoint python3 -v "${repo}/tests:/tests:ro" \
  -v "${root}/input:/data" "${baseline}" /tests/generate_denoising_fixture.py \
  /data/patch.laz /data/source.laz /data/fixture.json >> "${root}/fixture.log" 2>&1
docker run "${limits[@]}" --entrypoint python3 -v "${repo}/tests:/tests:ro" \
  -v "${root}/input:/data" "${baseline}" /tests/add_original_point_index.py \
  /data/source.laz /data/indexed.laz >> "${root}/fixture.log" 2>&1
docker run "${limits[@]}" --entrypoint python3 -v "${repo}/tests:/tests:ro" \
  -v "${root}/input:/data" "${baseline}" /tests/generate_high_ordinal_copc_input.py \
  /data/indexed.laz /data/high.laz >> "${root}/fixture.log" 2>&1
docker run "${limits[@]}" --entrypoint untwine -v "${root}/input:/data" \
  "${untwine}" -i /data/high.laz -o /data/input.copc.laz >> "${root}/fixture.log" 2>&1
for dimension in PredInstance PredInstance_SAT PredInstance_FM; do
  docker run "${limits[@]}" --entrypoint python3 -v "${repo}/tests:/tests:ro" \
    -v "${root}/input:/data" "${baseline}" /tests/compare_las_point_order.py \
    /data/source.laz /data/input.copc.laz "/data/order_${dimension}.json" \
    --dimension "${dimension}" >> "${root}/fixture.log" 2>&1
done
fi
python3 - "${root}" "${baseline}" "${untwine}" <<'PY'
import hashlib, json, pathlib, sys
root = pathlib.Path(sys.argv[1])
manifest = dict(runtime_image=sys.argv[2], untwine_image=sys.argv[3],
    source_sha256=hashlib.sha256((root / 'input/input.copc.laz').read_bytes()).hexdigest(),
    source_snapshot_sha256={str(p.relative_to(root / 'sources')): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted((root / 'sources').rglob('*')) if p.is_file()},
    cpus=2, memory_bytes=4 * 1024**3, internal_memory_gib=3.4, network='none', gpu=False)
(root / 'runtime.json').write_text(json.dumps(manifest, indent=2) + '\n')
PY
run_case() {
  local name="$1" workers="$2" target="$3"; shift 3
  local case_dir="${root}/${name}"
  mkdir -p "${case_dir}/output" "${case_dir}/scratch"
  echo "Running denoising acceptance: ${name}"
  local started="${SECONDS}"
  docker run "${limits[@]}" -v "${root}/sources:/opt/foreststructure:ro" --entrypoint Rscript \
    --workdir /output --env TMPDIR=/scratch --env PYTHONDONTWRITEBYTECODE=1 \
    --env GALAXY_SLOTS=2 --env GALAXY_MEMORY_MB=4096 --env FORESTSTRUCTURE_THREADS=2 \
    --env "FORESTSTRUCTURE_CATALOG_WORKERS=${workers}" \
    -v "${root}/input:/input:ro" -v "${case_dir}/output:/output" \
    -v "${case_dir}/scratch:/scratch" "${baseline}" \
    /opt/foreststructure/run_julia_memory_safe.R --point-cloud /input/input.copc.laz \
    --instance-dimension PredInstance,PredInstance_SAT,PredInstance_FM \
    --memory-budget-gib 3.4 --target-chunk-points "${target}" "$@" \
    > "${case_dir}/run.log" 2>&1
  echo "$((SECONDS - started))" > "${case_dir}/elapsed_seconds.txt"
  python3 "${repo}/validation/validate_output_artifacts.py" foreststructure \
    "${case_dir}/output" --report "${case_dir}/artifact_validation.json"
  test -z "$(find "${case_dir}/scratch" -mindepth 1 -print -quit)"
}
# The ordinary production invocation must not run SOR without a denoise flag.
run_case none 1 50000
run_case sor-w1 1 50000 --denoise sor
run_case sor-w2 2 5000 --denoise sor
docker run "${limits[@]}" --entrypoint python3 -v "${repo}/tests:/tests:ro" \
  -v "${root}:/data" "${baseline}" /tests/compare_denoising.py /data
echo 'Denoising acceptance passed: SOR off by default, uncapped heights, raw counts and execution invariance'
