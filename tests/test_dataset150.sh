#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
point_cloud="${1:-${DATASET150_LAZ:-}}"
if [[ -z "${point_cloud}" || ! -f "${point_cloud}" ]]; then
  echo "usage: DATASET150_LAZ=/path/to/80.laz make test-dataset150" >&2
  exit 2
fi
point_cloud="$(realpath "${point_cloud}")"
expected_laz_sha256="d767997ce5868a98fb9dc0a688cb3d300fc8b7cee45571eadcf9b6102c9f9789"
actual_laz_sha256="$(sha256sum "${point_cloud}" | cut -d ' ' -f 1)"
if [[ "$(basename "${point_cloud}")" != "80.laz" ||
      "${actual_laz_sha256}" != "${expected_laz_sha256}" ]]; then
  echo "dataset 150 regression requires the original 80.laz, not a COPC or another LAZ" >&2
  exit 2
fi

image_name="${FORESTSTRUCTURE_IMAGE:-3dtrees-foreststructure:test}"
fixture_dir="${repo_dir}/tests/fixtures/dataset150"
threads="${FORESTSTRUCTURE_CPUS:-10}"
memory="${FORESTSTRUCTURE_MEMORY:-100g}"
temporary_output="false"
if [[ -n "${DATASET150_RESULTS_DIR:-}" ]]; then
  output_dir="${DATASET150_RESULTS_DIR}"
  mkdir -p "${output_dir}"
else
  output_dir="$(mktemp -d)"
  temporary_output="true"
fi
if [[ "${temporary_output}" == "true" ]]; then
  trap 'rm -rf "${output_dir}"' EXIT
fi

if [[ "${FORESTSTRUCTURE_SKIP_BUILD:-0}" != "1" ]]; then
  docker build --tag "${image_name}" "${repo_dir}"
fi

docker run --rm --network none \
  --cpus "${threads}" \
  --memory "${memory}" \
  --memory-swap "${memory}" \
  --user "$(id -u):$(id -g)" \
  --volume "${point_cloud}:/in/$(basename "${point_cloud}"):ro" \
  --volume "${fixture_dir}:/fixtures:ro" \
  --volume "${output_dir}:/out" \
  --entrypoint Rscript \
  "${image_name}" \
  /opt/foreststructure/run.R \
  --point-cloud "/in/$(basename "${point_cloud}")" \
  --aoi /fixtures/150_upstream_tiles.geojson \
  --dataset-id 150 \
  --output-dir /out \
  --threads "${threads}" \
  --catalog-workers 2 \
  --dtm-chunk-size 60 \
  --chunk-size 60 \
  --dtm-buffer 20 \
  --performance-report

python "${repo_dir}/tests/compare_dataset150.py" \
  "${output_dir}" \
  "${fixture_dir}" \
  "$(basename "${point_cloud}")"

for artifact in \
  150_PredInstance_results.csv \
  150_PredInstance_segment_diagnostics.csv \
  150_dtm.tif \
  150_chm.tif \
  150_PredInstance_tiles.geojson \
  150_tiles.png \
  150_performance.csv; do
  test -s "${output_dir}/${artifact}"
done
test ! -e "${output_dir}/chm"
if find "${output_dir}" -maxdepth 1 -type d -name '.*_chunks' -print -quit |
  grep -q .; then
  echo "temporary raster chunk directories were not cleaned" >&2
  exit 1
fi

docker run --rm --network none \
  --entrypoint Rscript \
  --volume "${output_dir}:/out:ro" \
  "${image_name}" \
  -e '
    library(terra)
    dtm <- rast("/out/150_dtm.tif")
    chm <- rast("/out/150_chm.tif")
    stopifnot(!inMemory(dtm), !inMemory(chm))
    stopifnot(all(abs(res(dtm) - c(1, 1)) < 1e-9))
    stopifnot(all(abs(res(chm) - c(0.5, 0.5)) < 1e-9))
    stopifnot(!nzchar(crs(dtm)), !nzchar(crs(chm)))
  '

echo "dataset 150 artifacts verified in ${output_dir}"
