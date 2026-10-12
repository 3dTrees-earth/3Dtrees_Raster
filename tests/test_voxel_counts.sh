#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
image="${FORESTSTRUCTURE_JULIA_IMAGE:-3dtrees-foreststructure:copc-local}"
root="$(mktemp -d)"
trap 'rm -rf "${root}"' EXIT
mkdir -p "${root}/input" "${root}/baseline" "${root}/counts" "${root}/narrow" "${root}/work"
docker run --rm --network none --cpus 2 --memory 4g --user "$(id -u):$(id -g)" \
  --env FORESTSTRUCTURE_FIXTURE_EPSG=32632 \
  --entrypoint Rscript --volume "${repo_dir}/tests:/tests:ro" \
  --volume "${root}/input:/generated" "${image}" \
  /tests/generate_julia_memory_safe_fixture.R /generated/source.laz /generated/aoi.geojson distinct-dimensions
run() {
  local output="$1"
  shift
  docker run --rm --network none --cpus 2 --memory 4g --memory-swap 4g \
    --user "$(id -u):$(id -g)" --env FORESTSTRUCTURE_THREADS=2 \
    --volume "${root}/input:/in:ro" --volume "${root}/${output}:/out" \
    --volume "${root}/work:/work" "${image}" \
    --point-cloud /in/source.laz --aoi-geojson /in/aoi.geojson --dataset-id 150 \
    --output-dir /out --temp-dir /work --memory-budget-gib 3 "$@"
}
run baseline
run counts --export-voxel-counts --voxel-count-chunk-points 257
for dimension in PredInstance PredInstance_SAT PredInstance_FM; do
  for suffix in results.csv segment_diagnostics.csv tiles.geojson; do
    cmp "${root}/baseline/150_${dimension}_${suffix}" "${root}/counts/150_${dimension}_${suffix}"
  done
done
for suffix in dtm.tif chm.tif; do
  cmp "${root}/baseline/150_${suffix}" "${root}/counts/150_${suffix}"
done
python3 "${repo_dir}/validation/validate_output_artifacts.py" 150 "${root}/counts"
jq -e '.voxel_counts.enabled == true' "${root}/counts/150_julia_memory_safe_run.json"
jq -e '.point_accounting.included > 0 and .peak_batch_points <= 257' "${root}/counts/150_voxel_counts.json"
cat > "${root}/input/narrow.geojson" <<'JSON'
{"type":"Polygon","coordinates":[[[1001,2001],[1005,2001],[1005,2005],[1001,2005],[1001,2001]]]}
JSON
run narrow --export-voxel-counts --aoi-geojson /in/narrow.geojson
jq -e '.features | length == 0' "${root}/narrow/150_PredInstance_tiles.geojson"
jq -e '.point_accounting.included > 0' "${root}/narrow/150_voxel_counts.json"
python3 "${repo_dir}/validation/validate_output_artifacts.py" 150 "${root}/narrow"
test -z "$(find "${root}/work" -mindepth 1 -print -quit)"
echo "voxel-count CLI, publication and unchanged-science acceptance passed"
