#!/usr/bin/env bash

# Stop on the first failed command, reject unset variables, and propagate
# failures from any command inside a pipeline.
set -euo pipefail

# Resolve the repository root from this script's location so the script works
# independently of the caller's current working directory.
repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Read the ordered LAS/LAZ source path from the first positional argument.
source_laz="${1:-}"
# Read the destination COPC path from the second positional argument.
output_copc="${2:-}"
# Use the optional third argument for the validation report, or place the
# report beside the destination COPC when no explicit path is supplied.
report="${3:-${output_copc}.order.json}"

# Reject a missing source file or destination before creating temporary data.
if [[ ! -f "${source_laz}" || -z "${output_copc}" ]]; then
  # Show the exact required and optional arguments to help the caller recover.
  echo "usage: $0 ORIGINAL.laz OUTPUT.copc.laz [ORDER_REPORT.json]" >&2
  # Exit with the conventional code for invalid command-line usage.
  exit 2
fi

# Canonicalize the existing source path to make Docker mounts unambiguous.
source_laz="$(realpath "${source_laz}")"
# Canonicalize the not-yet-created destination path without requiring it to
# exist already.
output_copc="$(realpath -m "${output_copc}")"
# Canonicalize the not-yet-created report path for the same reason.
report="$(realpath -m "${report}")"
# Keep the output parent separately because temporary files must live on the
# same filesystem as the final atomic move.
output_parent="$(dirname "${output_copc}")"
# Create both destination directories before allocating temporary work.
mkdir -p "${output_parent}" "$(dirname "${report}")"
# Create a unique, private work directory beside the final COPC.
work_dir="$(mktemp -d "${output_parent}/.indexed_copc_XXXXXXXX")"
# Name the intermediate ordered LAZ that will receive OriginalPointIndex.
indexed_laz="${work_dir}/indexed_original.laz"
# Name the intermediate COPC that is withheld until all validation passes.
temporary_copc="${work_dir}/indexed.copc.laz"

# Remove every intermediate artifact on normal exit, validation failure, or
# interruption; the final COPC is outside this directory after publication.
cleanup() {
  # Delete only the unique work directory created above.
  rm -rf "${work_dir}"
}
# Register cleanup immediately so later failures cannot leave large LAS files.
trap cleanup EXIT

# Select the ForestStructure image that supplies the indexing and validation
# Python dependencies, while allowing callers to override the default tag.
forest_image="${FORESTSTRUCTURE_JULIA_IMAGE:-3dtrees-foreststructure:copc-local}"
# Select the Smart Tile image that supplies Untwine for COPC construction.
smart_tile_image="${FORESTSTRUCTURE_SMART_TILE_IMAGE:-3dtrees-smart-tile:34f492a-v2.3.1}"
# Limit each conversion stage to the requested CPU allocation.
test_cpus="${FORESTSTRUCTURE_TEST_CPUS:-10}"
# Apply the same hard memory and swap limit to every conversion stage.
container_memory_gib="${FORESTSTRUCTURE_TEST_MEMORY_GIB:-30}"
# Preserve host ownership on files written by the otherwise isolated containers.
host_uid="$(id -u)"
# Preserve the host group for the same output-permission reason.
host_gid="$(id -g)"

# `docker run` starts the isolated indexing process.
# `--rm` removes the container as soon as indexing exits.
# `--network none` prevents the conversion from contacting external services.
# `--cpus` enforces the configured CPU ceiling.
# `--memory` enforces the configured resident-memory ceiling.
# `--memory-swap` prevents swap from extending that memory ceiling.
# `--user` makes generated files belong to the invoking host user and group.
# The first `--volume` exposes the Python helper code read-only.
# The second `--volume` exposes the ordered source cloud read-only.
# The third `--volume` is the only writable mount and contains temporary data.
# `--entrypoint python3` runs the indexing helper instead of the image default.
# `forest_image` supplies laspy and the other pinned Python dependencies.
# `add_original_point_index.py` copies points in source order and appends the key.
# The final two arguments are the container input and intermediate output paths.
docker run --rm --network none \
  --cpus "${test_cpus}" \
  --memory "${container_memory_gib}g" \
  --memory-swap "${container_memory_gib}g" \
  --user "${host_uid}:${host_gid}" \
  --volume "${repo_dir}/tests:/tests:ro" \
  --volume "${source_laz}:/in/original.laz:ro" \
  --volume "${work_dir}:/work" \
  --entrypoint python3 \
  "${forest_image}" \
  /tests/add_original_point_index.py \
  /in/original.laz /work/indexed_original.laz

# `docker run` starts the isolated COPC construction process.
# `--rm` removes the container as soon as Untwine exits.
# `--network none` prevents COPC construction from contacting external services.
# `--cpus`, `--memory`, and `--memory-swap` enforce the same resource ceilings.
# `--user` preserves host ownership on the generated COPC.
# `--volume` exposes only the private work directory to this container.
# `--entrypoint untwine` selects the COPC builder instead of the image default.
# `smart_tile_image` supplies the pinned Untwine implementation.
# `-i` selects the indexed ordered LAZ produced by the first container.
# `-o` writes the spatially organized COPC to its unpublished temporary path.
docker run --rm --network none \
  --cpus "${test_cpus}" \
  --memory "${container_memory_gib}g" \
  --memory-swap "${container_memory_gib}g" \
  --user "${host_uid}:${host_gid}" \
  --volume "${work_dir}:/work" \
  --entrypoint untwine \
  "${smart_tile_image}" \
  -i /work/indexed_original.laz -o /work/indexed.copc.laz

# `docker run` starts the isolated validation process.
# `--rm` removes the validator container when the comparison exits.
# `--network none` guarantees validation uses only the mounted evidence.
# `--cpus`, `--memory`, and `--memory-swap` enforce the same resource ceilings.
# `--user` preserves host ownership on the generated JSON report.
# The first `--volume` exposes the validator source code read-only.
# The second `--volume` exposes the original ordered cloud read-only.
# The third `--volume` exposes the unpublished COPC candidate read-only.
# The fourth `--volume` exposes the report directory as the sole write target.
# `--entrypoint python3` runs the validator instead of the image default.
# `forest_image` supplies the pinned LAS/COPC Python dependencies.
# `validate_indexed_copc_streaming.py` performs the streaming comparison.
# The next two arguments identify the original and candidate container paths.
# The report argument retains the caller-selected filename inside `/report`.
# Redirecting stdout keeps successful batch logs concise; failures still abort.
docker run --rm --network none \
  --cpus "${test_cpus}" \
  --memory "${container_memory_gib}g" \
  --memory-swap "${container_memory_gib}g" \
  --user "${host_uid}:${host_gid}" \
  --volume "${repo_dir}/validation:/validators:ro" \
  --volume "${source_laz}:/in/original.laz:ro" \
  --volume "${temporary_copc}:/in/indexed.copc.laz:ro" \
  --volume "$(dirname "${report}"):/report" \
  --entrypoint python3 \
  "${forest_image}" \
  /validators/validate_indexed_copc_streaming.py \
  /in/original.laz /in/indexed.copc.laz \
  "/report/$(basename "${report}")" \
  >/dev/null

# `jq -e` exits nonzero unless the entire validation expression is true.
# `.status` requires the validator's overall verdict to be valid.
# `.header_compatible` requires exact stored bounds, quantization and CRS.
# `.fingerprint_equal` requires matching indexed XYZ/scientific tuples.
# `.unique_complete_range` rejects missing, duplicate, or out-of-range indices.
# `extra_dimensions[0]` keeps OriginalPointIndex in the first streamable slot.
# The nested expression checks every requested scientific dimension in turn.
# `index($dimension)` obtains each candidate ExtraByte's zero-based position.
# `// 99` turns a missing dimension into a guaranteed validation failure.
# `< 9` limits each dimension to lidR's selectively streamable ExtraByte range.
# `all` requires that positional constraint for every requested dimension.
# `.source.scales` preserves the exact XYZ quantization scales.
# `.source.offsets` preserves the exact XYZ coordinate offsets.
# The report path is passed explicitly and successful JSON output is suppressed.
jq -e '
  .status == "valid" and
  .header_compatible == true and
  .fingerprint_equal == true and
  .index_validation.unique_complete_range == true and
  .candidate.extra_dimensions[0] == "OriginalPointIndex" and
  (. as $report |
    [$report.dimensions[] as $dimension |
      (($report.candidate.extra_dimensions | index($dimension)) // 99) < 9] |
    all) and
  .source.scales == .candidate.scales and
  .source.offsets == .candidate.offsets
' "${report}" >/dev/null

# Publish only the fully validated temporary COPC; because it was built beside
# the destination, this rename is atomic on the target filesystem.
mv "${temporary_copc}" "${output_copc}"
# Print the final path as a concise success record for callers and batch logs.
echo "original-order COPC built: ${output_copc}"
