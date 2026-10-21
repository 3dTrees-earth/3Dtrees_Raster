#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
image="${FORESTSTRUCTURE_JULIA_IMAGE:-3dtrees-foreststructure:pr5-opt-9421302}"
untwine_image="${FORESTSTRUCTURE_SMART_TILE_IMAGE:-smarttile:v2.4.0-4de4393}"
root="$(mktemp -d)"
cleanup() {
  if [[ "${FORESTSTRUCTURE_KEEP_TEST_ROOT:-0}" == 1 ]]; then
    echo "Evidence: ${root}"
  else
    rm -rf "${root}"
  fi
}
trap cleanup EXIT
mkdir -p "${root}/input" "${root}/oracle" "${root}/full" "${root}/aoi" "${root}/narrow" "${root}/failure" "${root}/work"
# Freeze sources before any analysis so concurrent edits cannot mix policies.
mkdir -p "${root}/sources"
cp -a "${repo_dir}/src/." "${root}/sources/"
python3 - "${root}" <<'PYCODE'
import hashlib, json, pathlib, sys
root = pathlib.Path(sys.argv[1])
(root / 'source_snapshot_sha256.json').write_text(json.dumps({
    str(p.relative_to(root / 'sources')): hashlib.sha256(p.read_bytes()).hexdigest()
    for p in sorted((root / 'sources').rglob('*')) if p.is_file()
}, indent=2) + '\n')
PYCODE
limits=(--rm --network none --cpus 2 --memory 4g --memory-swap 4g --user "$(id -u):$(id -g)")
docker run "${limits[@]}" --env FORESTSTRUCTURE_FIXTURE_EPSG=32632 \
  --entrypoint Rscript -v "${repo_dir}/tests:/tests:ro" -v "${root}/input:/generated" \
  "${image}" /tests/generate_julia_memory_safe_fixture.R /generated/source.laz /generated/aoi.geojson distinct-dimensions
docker run "${limits[@]}" --entrypoint python3 -v "${repo_dir}/tests:/tests:ro" \
  -v "${root}/input:/data" "${image}" /tests/add_original_point_index.py /data/source.laz /data/low.laz
docker run "${limits[@]}" --entrypoint python3 -v "${repo_dir}/tests:/tests:ro" \
  -v "${root}/input:/data" "${image}" /tests/generate_high_ordinal_copc_input.py /data/low.laz /data/high.laz
for layout in low high; do
  docker run "${limits[@]}" --entrypoint untwine -v "${root}/input:/data" \
    "${untwine_image}" -i "/data/${layout}.laz" -o "/data/${layout}.copc.laz"
done
run() {
  local name="$1" cloud="$2"; shift 2
  local source_mount=(-v "${root}/sources:/opt/foreststructure:ro")
  # Reference and candidates share the current scientific height policy.
  local paths=(--output-dir /out --temp-dir /work)
  if [[ "${name}" == full ]]; then paths=(); fi
  docker run "${limits[@]}" "${source_mount[@]}" \
    --workdir /out --env TMPDIR=/work \
    --env GALAXY_SLOTS=2 --env GALAXY_MEMORY_MB=4096 --env FORESTSTRUCTURE_THREADS=2 \
    --env FORESTSTRUCTURE_CATALOG_WORKERS=1 \
    -v "${root}/input:/in:ro" -v "${root}/${name}:/out" -v "${root}/work:/work" \
    "${image}" --denoise none --point-cloud "/in/${cloud}" --original-point-cloud /in/source.laz \
    "${paths[@]}" "$@" > "${root}/${name}.log" 2>&1
}
run oracle low.copc.laz --dataset-id 150 --instance-dimension PredInstance,PredInstance_SAT,PredInstance_FM --aoi-geojson /in/aoi.geojson --export-voxel-counts --memory-budget-gib 3
run full high.copc.laz --target-chunk-points 500 --instance-dimension PredInstance,PredInstance_SAT,PredInstance_FM
run aoi low.copc.laz --aoi-geojson /in/aoi.geojson --instance-dimension PredInstance,PredInstance_SAT,PredInstance_FM
python3 "${repo_dir}/tests/compare_production_contract.py" "${root}"
jq -e '.voxel_counts.enabled and .footprint_source == "point_cloud_extent" and .scheduler_resources.threads == 2 and .memory_budget_gib == 3.4 and .copc_projection.dimension == ["OriginalPointIndex", "PredInstance", "PredInstance_SAT", "PredInstance_FM"] and .catalog_spatial_streaming == "indexed LAZ projection"' "${root}/full/foreststructure_julia_memory_safe_run.json"
python3 "${repo_dir}/validation/validate_output_artifacts.py" foreststructure "${root}/full"
cat > "${root}/input/narrow.geojson" <<'JSON'
{"type":"Polygon","coordinates":[[[1001,2001],[1005,2001],[1005,2005],[1001,2005],[1001,2001]]]}
JSON
run narrow high.copc.laz --aoi-geojson /in/narrow.geojson --instance-dimension PredInstance,PredInstance_SAT,PredInstance_FM
jq -e '.features | length == 0' "${root}/narrow/tiles/PredInstance_tiles.geojson"
jq -e '.point_accounting.included > 0 and .point_accounting.excluded_domain > 0' "${root}/narrow/foreststructure_voxel_counts.json"
if run failure high.copc.laz --instance-dimension missing; then
  echo 'Missing dimension unexpectedly succeeded' >&2; exit 1
fi
for cloud in sat_only no_instances; do
  request=(--instance-dimension PredInstance_SAT)
  if [[ "${cloud}" == sat_only ]]; then request=(--instance-dimension PredInstance_FM); fi
  if docker run "${limits[@]}" -v "${root}/sources:/opt/foreststructure:ro" \
    -v "${root}/input:/in:ro" -v "${root}/failure:/out" -v "${root}/work:/work" \
    "${image}" --denoise none --point-cloud "/in/${cloud}.laz" --output-dir /out --temp-dir /work \
    "${request[@]}" > "${root}/${cloud}.log" 2>&1; then
    echo "${cloud} unexpectedly succeeded" >&2; exit 1
  fi
done
grep -q 'missing requested instance dimension' "${root}/sat_only.log"
grep -q 'missing requested instance dimension' "${root}/no_instances.log"
test -z "$(find "${root}/failure" -mindepth 1 -print -quit)"
test -z "$(find "${root}/work" -mindepth 1 -print -quit)"
echo 'Production CLI, full-extent/AOI, high-ordinal COPC and exact science parity passed'
