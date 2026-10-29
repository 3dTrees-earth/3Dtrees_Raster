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
  docker run --rm --network none --cpus 2 --memory 4g --memory-swap 4g --user "$(id -u):$(id -g)" --env FORESTSTRUCTURE_THREADS=2 --env "FORESTSTRUCTURE_SEGMENT_BUCKET_COUNT=${bucket_count}" --volume "${test_root}/input/high_ordinal.laz:/in/high_ordinal.laz:ro" --volume "${test_root}/input/aoi.geojson:/in/aoi.geojson:ro" --volume "${test_root}/output:/out" --volume "${test_root}/work:/work" "${image_name}" --denoise none --point-cloud /in/high_ordinal.laz --aoi-geojson /in/aoi.geojson --dataset-id 150 --output-dir /out --temp-dir /work --memory-budget-gib "${memory_budget_gib}" --sensor "${sensor}" --country "${country}" --instance-dimension PredInstance,PredInstance_SAT,PredInstance_FM
}

test "$(docker run --rm --network none "${image_name}" --version)" = "1.0.0"

run_candidate 2>&1 | tee "${test_root}/first.log"
for dimension in PredInstance PredInstance_SAT PredInstance_FM; do
  grep -q "dimension=${dimension} .*selector=xyz1" "${test_root}/first.log"
  test -s "${test_root}/output/tile_metrics/${dimension}_tile_metrics.csv"
  test -s "${test_root}/output/segment_diagnostics/${dimension}_segment_diagnostics.csv"
  test -s "${test_root}/output/tiles/${dimension}_tiles.geojson"
done
cmp "${test_root}/output/tile_metrics/PredInstance_tile_metrics.csv" "${test_root}/output/tile_metrics/PredInstance_SAT_tile_metrics.csv"
cmp "${test_root}/output/tile_metrics/PredInstance_tile_metrics.csv" "${test_root}/output/tile_metrics/PredInstance_FM_tile_metrics.csv"
cmp "${test_root}/output/segment_diagnostics/PredInstance_segment_diagnostics.csv" "${test_root}/output/segment_diagnostics/PredInstance_SAT_segment_diagnostics.csv"
cmp "${test_root}/output/segment_diagnostics/PredInstance_segment_diagnostics.csv" "${test_root}/output/segment_diagnostics/PredInstance_FM_segment_diagnostics.csv"
jq -e '.processed_instance_dimensions | length == 3' "${test_root}/output/150_julia_memory_safe_run.json" >/dev/null
jq -e '.raster_format == "COG"' "${test_root}/output/150_julia_memory_safe_run.json" >/dev/null
raster_layouts() {
  docker run --rm --network none --entrypoint python3 --volume "$1:/check:ro" "${image_name}" -c '
import sys
from osgeo import gdal
for name in sys.argv[1:]:
    print(gdal.Info("/check/" + name, format="json").get("metadata", {}).get("IMAGE_STRUCTURE", {}).get("LAYOUT", "GTiff"))
' "${@:2}"
}
test "$(raster_layouts "${test_root}/output" 150_dtm.tif 150_chm.tif 150_dsm.tif 150_chm_filled.tif | sort -u)" = "COG"
# DSM and gap-filled CHM share the CHM grid (also checked before promotion).
jq -e '.surface_rasters.chm_filled.filled_cells >= 0 and .surface_rasters.chm_filled.empty_cells >= .surface_rasters.chm_filled.filled_cells' "${test_root}/output/150_julia_memory_safe_run.json" >/dev/null
docker run --rm --network none --entrypoint python3 --volume "${test_root}/output:/check:ro" "${image_name}" -c '
from osgeo import gdal
import numpy as np
gdal.UseExceptions()
read = lambda name: gdal.Open("/check/" + name)
chm, dsm, filled = read("150_chm.tif"), read("150_dsm.tif"), read("150_chm_filled.tif")
for other in (dsm, filled):
    assert other.GetGeoTransform() == chm.GetGeoTransform()
    assert (other.RasterXSize, other.RasterYSize) == (chm.RasterXSize, chm.RasterYSize)
    assert other.GetRasterBand(1).GetNoDataValue() == chm.GetRasterBand(1).GetNoDataValue()
    assert other.GetRasterBand(1).DataType == chm.GetRasterBand(1).DataType
nodata = chm.GetRasterBand(1).GetNoDataValue()
c, d, f = (r.GetRasterBand(1).ReadAsArray() for r in (chm, dsm, filled))
assert np.array_equal(c == nodata, d == nodata)
assert np.array_equal(f[c != nodata], c[c != nodata])
'
jq -e '.projections | length == 3' "${test_root}/output/150_julia_memory_safe_run.json" >/dev/null
jq -e '.adaptive_chunking.enabled == false and .adaptive_chunking.target_points == 50000000' "${test_root}/output/150_julia_memory_safe_run.json" >/dev/null
test -z "$(find "${test_root}/work" -mindepth 1 -print -quit)"

# A second run must recompute and replace the validated artifacts, not skip.
# A high bucket count forces tiny disk partitions without changing results.
run_candidate 256 NA NA 70 2>&1 | tee "${test_root}/second.log"
grep -q 'Validated and promoted' "${test_root}/second.log"
jq -e '.segment_bucket_count == 256' "${test_root}/output/150_julia_memory_safe_run.json" >/dev/null
jq -e '.memory_budget_gib == 3.4' "${test_root}/output/150_julia_memory_safe_run.json" >/dev/null

# Default point-noise filtering on plain LAZ numbers records in file order.
sor_output="${test_root}/sor_output"
sor_work="${test_root}/sor_work"
mkdir -p "${sor_output}" "${sor_work}"
docker run --rm --network none --cpus 2 --memory 4g --memory-swap 4g --user "$(id -u):$(id -g)" --env FORESTSTRUCTURE_THREADS=2 --volume "${test_root}/input/high_ordinal.laz:/in/high_ordinal.laz:ro" --volume "${test_root}/input/aoi.geojson:/in/aoi.geojson:ro" --volume "${sor_output}:/out" --volume "${sor_work}:/work" "${image_name}" --point-cloud /in/high_ordinal.laz --aoi-geojson /in/aoi.geojson --dataset-id 152 --output-dir /out --temp-dir /work --memory-budget-gib 3 --raster-format GTiff --denoise sor --instance-dimension PredInstance,PredInstance_SAT 2>&1 | tee "${test_root}/sor.log"
test "$(head -n 1 "${sor_output}/152_noise_points.csv")" = "OriginalPointIndex"
jq -e '.enabled == true' "${sor_output}/152_denoising.json" >/dev/null
test "$(raster_layouts "${sor_output}" 152_dtm.tif 152_chm.tif 152_dsm.tif 152_chm_filled.tif 152_voxel_counts.tif 152_voxel_counts_validity.tif | sort -u)" = "GTiff"
jq -e '.requested_is_copc == false and .copc_projection.original_point_index_added == true' "${sor_output}/152_julia_memory_safe_run.json" >/dev/null
for dimension in PredInstance PredInstance_SAT; do
  test -s "${sor_output}/tile_metrics/${dimension}_tile_metrics.csv"
  test -s "${sor_output}/segment_diagnostics/${dimension}_segment_diagnostics.csv"
  test -s "${sor_output}/tiles/${dimension}_tiles.geojson"
done
test -z "$(find "${sor_output}" -maxdepth 1 -name '152_PredInstance*' -print -quit)"
cmp "${sor_output}/tile_metrics/PredInstance_tile_metrics.csv" "${sor_output}/tile_metrics/PredInstance_SAT_tile_metrics.csv"
test -z "$(find "${sor_work}" -mindepth 1 -print -quit)"

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
