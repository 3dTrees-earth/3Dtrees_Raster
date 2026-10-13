#!/usr/bin/env bash
set -euo pipefail
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
image="${FORESTSTRUCTURE_JULIA_IMAGE:-3dtrees-foreststructure:copc-local}"
root="$(mktemp -d)"
trap 'rm -rf "${root}"' EXIT
mkdir -p "${root}/input" "${root}/two" "${root}/three" "${root}/work"
limits=(--rm --network none --cpus 2 --memory 4g --memory-swap 4g --user "$(id -u):$(id -g)")
docker run "${limits[@]}" --env FORESTSTRUCTURE_FIXTURE_EPSG=32632 \
  --entrypoint Rscript -v "${repo}/tests:/tests:ro" -v "${root}/input:/generated" \
  "${image}" /tests/generate_julia_memory_safe_fixture.R /generated/two.laz /generated/aoi.geojson distinct-dimensions
docker run "${limits[@]}" --entrypoint python3 -v "${root}/input:/data" "${image}" -c \
  'import laspy; p=laspy.read("/data/two.laz"); p.header.version=laspy.header.Version(1,4); p.header.add_crs(p.header.parse_crs().to_3d(),keep_compatibility=False); p.write("/data/three.laz")'
for mode in two three; do
  docker run "${limits[@]}" --env GALAXY_SLOTS=2 --env GALAXY_MEMORY_MB=4096 \
    --env TMPDIR=/work --workdir /out -v "${root}/input/${mode}.laz:/in/input.laz:ro" \
    -v "${root}/input/aoi.geojson:/in/aoi.geojson:ro" -v "${root}/${mode}:/out" \
    -v "${root}/work:/work" "${image}" --point-cloud /in/input.laz --aoi-geojson /in/aoi.geojson
done
docker run "${limits[@]}" --entrypoint python3 -v "${root}:/data:ro" \
  -v "${repo}/tests:/tests:ro" "${image}" /tests/compare_projected_3d_crs.py
echo 'Projected 3D CLI, embedded raster CRS and unchanged scientific values passed'
