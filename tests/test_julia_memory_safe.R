stopifnot(identical(JULIA_MEMORY_SAFE_SEGMENT_SELECTION, "xyz"))
stopifnot(identical(validated_catalog_worker_count("1", 20L), 1L))
stopifnot(identical(validated_catalog_worker_count("2", 20L), 2L))
invalid_catalog_workers <- vapply(c("0", "11", "1.5", "not-a-number"), function(value) {
  tryCatch(
    {
      validated_catalog_worker_count(value, 20L)
      FALSE
    },
    error = function(error) grepl(
      "FORESTSTRUCTURE_CATALOG_WORKERS must be an integer from 1 to 10",
      conditionMessage(error),
      fixed = TRUE
    )
  )
}, logical(1))
stopifnot(all(invalid_catalog_workers))
worker_budget_rejected <- tryCatch(
  {
    validated_catalog_worker_count("2", 1L)
    FALSE
  },
  error = function(error) grepl(
    "exceeds FORESTSTRUCTURE_THREADS",
    conditionMessage(error),
    fixed = TRUE
  )
)
stopifnot(worker_budget_rejected)
stopifnot(
  segment_selection_from_names(
    c("PredInstance", "unrelated"),
    "PredInstance",
    JULIA_MEMORY_SAFE_SEGMENT_SELECTION
  ) == "xyz1"
)
beyond_nine <- tryCatch({
  segment_selection_from_names(c(paste0("extra_", 1:14), "PredInstance"),
    "PredInstance", JULIA_MEMORY_SAFE_SEGMENT_SELECTION)
  FALSE
}, error = function(e) grepl("verified streaming projection", conditionMessage(e)))
stopifnot(beyond_nine)

dtm_body <- paste(deparse(body(dtm_chunk)), collapse = "\n")
stopifnot(grepl("last_returns = FALSE", dtm_body, fixed = TRUE))
tile_body <- paste(deparse(body(calculate_tile_metrics)), collapse = "\n")
# Tile reads must stay spatial/selective while requesting OriginalPointIndex in
# COPC mode; checking a literal ``select = "xyz"`` would enforce the old LAZ-only
# implementation and prevent deterministic order restoration.
stopifnot(grepl(
  "catalog_selection_for_dimensions",
  tile_body,
  fixed = TRUE
))
stopifnot(grepl("restore_original_point_order", tile_body, fixed = TRUE))
segment_body <- paste(deparse(body(accumulate_segments)), collapse = "\n")
stopifnot(!grepl('select = "*"', segment_body, fixed = TRUE))
stopifnot(!grepl('readLAS(point_cloud)', segment_body, fixed = TRUE))

unlimited <- list(memory_budget_gib = Inf)
stopifnot(isTRUE(assert_memory_budget(.Machine$double.xmax, unlimited, "test")))
impossible <- tryCatch(
  {
    assert_memory_budget(1, list(memory_budget_gib = 1e-12), "test")
    FALSE
  },
  error = function(error) grepl(
    "FORESTSTRUCTURE_FLAG memory_budget_exceeded",
    conditionMessage(error),
    fixed = TRUE
  )
)
stopifnot(impossible)

store_directory <- tempfile("segment_adapter_test_")
dir.create(store_directory)
piece <- function(voxels, apex_z) {
  list(
    voxels = data.table(
      instance_dimension = "PredInstance",
      instance_id = 1,
      voxel_x = voxels,
      voxel_y = voxels,
      voxel_z = voxels
    ),
    layers = data.table(
      instance_dimension = "PredInstance",
      instance_id = 1,
      occupied_layer = voxels
    ),
    apex = data.table(
      instance_dimension = "PredInstance",
      instance_id = 1,
      apex_x = 10,
      apex_y = 20,
      apex_z = apex_z
    )
  )
}
first_path <- file.path(store_directory, "first.rds")
second_path <- file.path(store_directory, "second.rds")
saveRDS(piece(1:3, 3), first_path, compress = FALSE)
saveRDS(piece(3:5, 5), second_path, compress = FALSE)
store <- list(
  chunks = list(
    setNames(first_path, "1"),
    setNames(second_path, "1")
  ),
  bucket_count = 1L
)
adapter_parameters <- list(
  voxel_resolution = 0.2,
  minimum_tree_voxels = 1L,
  apex_minimum_height = 0,
  minimum_tree_thickness = 0,
  minimum_occupied_layers = 1L,
  memory_budget_gib = Inf
)
disk_result <- finalize_segment_store(store, adapter_parameters, "disk")
memory_result <- finalize_segment_store(store, adapter_parameters, "memory")
stopifnot(identical(disk_result, memory_result))
unlink(store_directory, recursive = TRUE, force = TRUE)
