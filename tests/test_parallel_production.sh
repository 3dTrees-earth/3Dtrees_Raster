#!/usr/bin/env bash
# Differential acceptance: same scientific runtime, independent output/scratch
# directories, current sources at one reference worker and at 1/2/4 workers.
set -euo pipefail
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source_dir="${FORESTSTRUCTURE_SOURCE_DIR:-${repo}/src}"
source_dir="$(cd "${source_dir}" && pwd)"
baseline="${FORESTSTRUCTURE_BASELINE_IMAGE:-sha256:3bb86542fda4dfc189a88f6611c8da962b7f524edb4e1a473afe2684d5baf1e4}"
untwine="${FORESTSTRUCTURE_SMART_TILE_IMAGE:-smarttile:v2.4.0-4de4393}"
# Resolve mutable tags once; record and execute exactly these local image IDs.
baseline="$(docker image inspect --format '{{.Id}}' "${baseline}")"
untwine="$(docker image inspect --format '{{.Id}}' "${untwine}")"
root="$(mktemp -d "${TMPDIR:-/tmp}/foreststructure-parallel.XXXXXX")"
finish() {
  local status="$?"
  if [[ "${status}" != 0 || "${FORESTSTRUCTURE_KEEP_TEST_ROOT:-0}" == 1 ]]; then
    echo "Parallel acceptance evidence: ${root} (exit ${status})"
  else
    rm -rf "${root}"
  fi
}
trap finish EXIT
mkdir -p "${root}/input" "${root}/sources"
# Snapshot candidate sources so edits during the test cannot change later cases.
cp -a "${source_dir}/." "${root}/sources/"
limits=(--rm --network none --cpus 4 --memory 4g --memory-swap 4g --user "$(id -u):$(id -g)")
python3 - "${root}" "${source_dir}" "${baseline}" "${untwine}" <<'PY'
import hashlib, json, pathlib, sys
root, source, baseline, untwine = sys.argv[1:]
root = pathlib.Path(root)
hashes = {str(p.relative_to(root / 'sources')): hashlib.sha256(p.read_bytes()).hexdigest()
          for p in sorted((root / 'sources').rglob('*')) if p.is_file()}
(root / 'runtime.json').write_text(json.dumps(dict(
    runtime_image=baseline, reference_policy="current_source_one_worker", untwine_image=untwine, candidate_source=source,
    candidate_source_sha256=hashes, cpus=4, memory_bytes=4 * 1024**3,
    internal_memory_gib=3.4, network='none', gpu=False,
    target_chunk_points=5000, voxel_batch_points=250000,
    cases=[f'{scope}-{variant}' for scope in ('full', 'aoi')
           for variant in ('baseline-w1', 'candidate-w1', 'candidate-w2', 'candidate-w4')]
), indent=2) + '\n')
PY
docker run "${limits[@]}" --env FORESTSTRUCTURE_FIXTURE_EPSG=32632 \
  --entrypoint Rscript -v "${repo}/tests:/tests:ro" -v "${root}/input:/generated" \
  "${baseline}" /tests/generate_julia_memory_safe_fixture.R \
  /generated/patch.laz /generated/unused_aoi.geojson distinct-dimensions \
  > "${root}/fixture.log" 2>&1
# Four touching patches exercise multiple baseline raster chunks as well as
# adaptive read windows. Duplicate boundary records are deliberately retained.
docker run "${limits[@]}" --entrypoint python3 -v "${root}/input:/data" \
  "${baseline}" -c '
import copy, json, laspy, numpy as np
source = laspy.read("/data/patch.laz")
result = laspy.LasData(copy.deepcopy(source.header))
result.points = laspy.ScaleAwarePointRecord.zeros(4 * len(source.points), header=result.header)
for i, (dx, dy) in enumerate(((0, 0), (40, 0), (0, 40), (40, 40))):
    block = slice(i * len(source.points), (i + 1) * len(source.points))
    result.points.array[block] = source.points.array
    result.X[block] += round(dx / source.header.scales[0])
    result.Y[block] += round(dy / source.header.scales[1])
    for dimension in ("PredInstance", "PredInstance_SAT", "PredInstance_FM"):
        values = np.asarray(source[dimension]).copy()
        values[values > 0] += 1000 * i
        result[dimension][block] = values
result.write("/data/source.laz")
ring = [[1002, 2003], [1078, 2003], [1078, 2077], [1002, 2077], [1002, 2003]]
with open("/data/aoi.geojson", "w") as handle:
    json.dump(dict(type="Polygon", coordinates=[ring]), handle)
' >> "${root}/fixture.log" 2>&1
docker run "${limits[@]}" --entrypoint python3 -v "${repo}/tests:/tests:ro" \
  -v "${root}/input:/data" "${baseline}" \
  /tests/add_original_point_index.py /data/source.laz /data/indexed.laz \
  >> "${root}/fixture.log" 2>&1
docker run "${limits[@]}" --entrypoint python3 -v "${repo}/tests:/tests:ro" \
  -v "${root}/input:/data" "${baseline}" \
  /tests/generate_high_ordinal_copc_input.py /data/indexed.laz /data/high.laz \
  >> "${root}/fixture.log" 2>&1
docker run "${limits[@]}" --entrypoint untwine -v "${root}/input:/data" \
  "${untwine}" -i /data/high.laz -o /data/input.copc.laz \
  >> "${root}/fixture.log" 2>&1
for dimension in PredInstance PredInstance_SAT PredInstance_FM; do
  docker run "${limits[@]}" --entrypoint python3 -v "${repo}/tests:/tests:ro" \
    -v "${root}/input:/data" "${baseline}" /tests/compare_las_point_order.py \
    /data/source.laz /data/input.copc.laz "/data/order_${dimension}.json" \
    --dimension "${dimension}" >> "${root}/fixture.log" 2>&1
done
cat > "${root}/record_resources.py" <<'PY'
import json, pathlib, sys
base = pathlib.Path('/sys/fs/cgroup')
def read(name):
    path = base / name
    if not path.is_file():
        raise RuntimeError('cgroup v2 resource evidence is required: ' + name)
    return path.read_text().strip()
json.dump(dict(exit_code=int(sys.argv[1]),
    memory_peak_bytes=int(read('memory.peak')),
    memory_limit_bytes=int(read('memory.max')),
    memory_swap_limit_bytes=int(read('memory.swap.max')),
    memory_events={key: int(value) for key, value in
                   (line.split() for line in read('memory.events').splitlines())},
    cpu_stat={key: int(value) for key, value in
              (line.split() for line in read('cpu.stat').splitlines())}),
    open('/case/resources.json', 'w'), indent=2)
PY
run_case() {
  local scope="$1" kind="$2" workers="$3"
  local name="${scope}-${kind}-w${workers}"
  local case_dir="${root}/${name}"
  mkdir -p "${case_dir}/output" "${case_dir}/scratch"
  # Compare worker counts under one current scientific policy/source snapshot.
  local mount=(-v "${root}/sources:/opt/foreststructure:ro") aoi=()
  if [[ "${scope}" == aoi ]]; then aoi=(--aoi-geojson /input/aoi.geojson); fi
  local command=(docker run "${limits[@]}" "${mount[@]}" \
    --entrypoint bash --workdir /output --env TMPDIR=/scratch \
    --env PYTHONDONTWRITEBYTECODE=1 --env GALAXY_SLOTS=4 --env GALAXY_MEMORY_MB=4096 \
    --env FORESTSTRUCTURE_THREADS=4 --env "FORESTSTRUCTURE_CATALOG_WORKERS=${workers}" \
    -v "${root}/input:/input:ro" -v "${case_dir}/output:/output" \
    -v "${case_dir}/scratch:/scratch" -v "${case_dir}:/case" \
    -v "${root}/record_resources.py:/record_resources.py:ro" "${baseline}" -c \
    '"$@"; status=$?; python3 /record_resources.py "$status" || exit 99; exit "$status"' \
    bash Rscript /opt/foreststructure/run_julia_memory_safe.R \
    --denoise none --point-cloud /input/input.copc.laz --memory-budget-gib 3.4 \
    --instance-dimension PredInstance,PredInstance_SAT,PredInstance_FM \
    --target-chunk-points 5000 --voxel-count-chunk-points 250000 "${aoi[@]}")
  python3 - "${case_dir}/command.json" "${command[@]}" <<'PY'
import json, pathlib, sys
pathlib.Path(sys.argv[1]).write_text(json.dumps(sys.argv[2:], indent=2) + '\n')
PY
  echo "Running ${name}"
  local started="${SECONDS}"
  "${command[@]}" > "${case_dir}/run.log" 2>&1
  echo "$((SECONDS - started))" > "${case_dir}/elapsed_seconds.txt"
  python3 "${repo}/validation/validate_output_artifacts.py" foreststructure \
    "${case_dir}/output" --report "${case_dir}/artifact_validation.json"
  test -z "$(find "${case_dir}/scratch" -mindepth 1 -print -quit)"
}
for scope in full aoi; do
  run_case "${scope}" baseline 1
  for workers in 1 2 4; do run_case "${scope}" candidate "${workers}"; done
done
docker run "${limits[@]}" --entrypoint python3 \
  -v "${repo}/tests:/tests:ro" -v "${root}:/data" "${baseline}" \
  /tests/compare_parallel_production.py /data
echo 'Parallel production acceptance passed: current-source reference versus 1/2/4 workers, full extent and AOI'
