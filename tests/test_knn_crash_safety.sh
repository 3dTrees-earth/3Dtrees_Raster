#!/usr/bin/env bash
# kNN filter crash safety on a real COPC: every case must flag exactly the points
# of the reference noise list (raster_noise_points.csv of a normal run).
#   COPC=... REFERENCE=... RASTER_IMAGE=... bash tests/test_knn_crash_safety.sh
# Cases: a 4 GB container; admission estimates 5x and 20x too low; a killed worker.
set -euo pipefail
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
image="${RASTER_IMAGE:-3dtrees-raster:test}"
threads="${RASTER_KNN_SAFETY_THREADS:-10}"
root="$(mktemp -d)"
trap 'rm -rf "${root}"' EXIT
run_case() {  # name case memory_gib
  local name=$1 case=$2 memory=$3
  mkdir -p "${root}/${name}"
  local budget=$(( memory * 1024 * 1024 * 1024 * 85 / 100 ))
  docker run --rm --network none --user "$(id -u):$(id -g)" --cpus "${threads}" \
    --memory "${memory}g" --memory-swap "${memory}g" --env OMP_NUM_THREADS=1 \
    --entrypoint python3 --volume "${repo}/tests:/tests:ro" \
    --volume "${COPC}:/in/cloud.copc.laz:ro" --volume "${REFERENCE}:/in/reference.csv:ro" \
    --volume "${root}/${name}:/out" "${image}" \
    /tests/knn_crash_safety.py "${case}" /in/cloud.copc.laz /in/reference.csv /out "${threads}" "${budget}" \
    > "${root}/${name}.log" 2>&1 || { tail -20 "${root}/${name}.log"; echo "kNN crash safety ${name}: FAILED"; exit 1; }
  echo "kNN crash safety ${name}: $(tail -1 "${root}/${name}.log")$(grep -o 'KNN_RETRY[^$]*' "${root}/${name}.log" | head -1 | sed 's/^/; /')"
}
run_case cap4g plain 4
run_case low5x low:0.2 8
run_case low20x low:0.05 8
run_case killed killed 8
echo "kNN crash safety passed"
