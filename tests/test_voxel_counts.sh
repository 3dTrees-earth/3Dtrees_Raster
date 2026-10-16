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
    --denoise none --point-cloud /in/source.laz --aoi-geojson /in/aoi.geojson --dataset-id 150 \
    --instance-dimension PredInstance,PredInstance_SAT,PredInstance_FM \
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
# Advanced scientific parameters reach the analysis, the voxel export and provenance.
mkdir -p "${root}/custom" "${root}/invalid"
run custom --export-voxel-counts --voxel-count-resolution 1 --voxel-resolution 0.25 \
  --tile-size 10 --chm-resolution 1 --minimum-tree-voxels 50
jq -e '.voxel_size_m == [1, 1, 1]' "${root}/custom/150_voxel_counts.json"
jq -e '.scientific_parameters | .tile_size == 10 and .voxel_resolution == 0.25
  and .voxel_count_resolution == 1 and .chm_resolution == 1 and .minimum_tree_voxels == 50
  and .dtm_resolution == 1' "${root}/custom/150_julia_memory_safe_run.json"
python3 "${repo_dir}/validation/validate_output_artifacts.py" 150 "${root}/custom"
if run invalid --export-voxel-counts --tile-size 0 2>"${root}/invalid.log"; then
  echo "--tile-size 0 was accepted" >&2
  exit 1
fi
grep -q -- "--tile-size must be greater than zero" "${root}/invalid.log"
if run invalid --export-voxel-counts --voxel-resolution 0.3 2>"${root}/invalid.log"; then
  echo "--voxel-resolution 0.3 was accepted" >&2
  exit 1
fi
grep -q -- "--voxel-resolution must divide 1 m evenly" "${root}/invalid.log"
test -z "$(find "${root}/work" -mindepth 1 -print -quit)"
echo "voxel-count CLI, publication and unchanged-science acceptance passed"
