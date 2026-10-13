#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -lt 2 || "$#" -gt 3 ]]; then
  echo "usage: $0 IMAGE WORK_DIRECTORY [THREADS]" >&2
  exit 2
fi

image_name="$1"
work_dir="$2"
threads="${3:-1}"
dataset_id="999999"
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fixture_dir="${work_dir}/fixtures"
image_slug="${image_name//[\/:]/_}_threads${threads}"
run_root="${work_dir}/runs/${image_slug}"
mkdir -p "${fixture_dir}" "${run_root}"

prepare_fixture() {
  local name="$1"
  local point_count="$2"
  local aoi_width="$3"
  if [[ -f "${fixture_dir}/${name}.laz" && -f "${fixture_dir}/${name}.geojson" ]]; then
    return
  fi
  docker run --rm --network none \
    --user "$(id -u):$(id -g)" \
    --entrypoint Rscript \
    --volume "${repo_dir}/benchmarks:/benchmarks:ro" \
    --volume "${fixture_dir}:/fixtures" \
    "${image_name}" \
    /benchmarks/generate_fixture.R \
    "/fixtures/${name}.laz" \
    "/fixtures/${name}.geojson" \
    "${point_count}" \
    "${aoi_width}"
}

run_case() {
  local name="$1"
  local output_dir="${run_root}/${name}"
  mkdir -p "${output_dir}"
  docker run --rm --network none \
    --user "$(id -u):$(id -g)" \
    --volume "${fixture_dir}:/in:ro" \
    --volume "${output_dir}:/out" \
    --entrypoint Rscript \
    "${image_name}" \
    /opt/foreststructure/run.R \
    --point-cloud "/in/${name}.laz" \
    --aoi "/in/${name}.geojson" \
    --dataset-id "${dataset_id}" \
    --output-dir /out \
    --threads "${threads}" \
    --instance-dimension PredInstance \
    --segment-diagnostics \
    --performance-report \
    2>&1 | tee "${output_dir}/container.log"
  sed -n '1,2p' "${output_dir}/${dataset_id}_performance.csv"
}

prepare_fixture small 50000 60
prepare_fixture large 500000 120
run_case small
run_case large
