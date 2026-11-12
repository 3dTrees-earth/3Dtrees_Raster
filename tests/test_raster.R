stopifnot(identical(SEGMENT_SELECTION, "xyz"))
stopifnot(identical(validated_catalog_worker_count("1", 20L), 1L))
stopifnot(identical(validated_catalog_worker_count("2", 20L), 2L))
invalid_catalog_workers <- vapply(c("0", "11", "1.5", "not-a-number"), function(value) {
  tryCatch(
    {
      validated_catalog_worker_count(value, 20L)
      FALSE
    },
    error = function(error) grepl(
      "RASTER_CATALOG_WORKERS must be an integer from 1 to 10",
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
    "exceeds RASTER_THREADS",
    conditionMessage(error),
    fixed = TRUE
  )
)
stopifnot(worker_budget_rejected)
stopifnot(
  segment_selection_from_names(
    c("PredInstance", "unrelated"),
    "PredInstance",
    SEGMENT_SELECTION
  ) == "xyz1"
)
beyond_nine <- tryCatch({
  segment_selection_from_names(c(paste0("extra_", 1:14), "PredInstance"),
    "PredInstance", SEGMENT_SELECTION)
  FALSE
}, error = function(e) grepl("verified streaming projection", conditionMessage(e)))
stopifnot(beyond_nine)

dtm_body <- paste(deparse(body(dtm_chunk)), collapse = "\n")
stopifnot(grepl("last_returns = FALSE", dtm_body, fixed = TRUE))
# Single tiles and tile blocks both read through read_tile_points().
tile_body <- paste(deparse(body(read_tile_points)), collapse = "\n")
stopifnot(grepl("read_tile_points", paste(deparse(body(calculate_tile_metrics)), collapse = "\n"),
                fixed = TRUE))
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
    "RASTER_FLAG memory_budget_exceeded",
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

# A block read split into tiles gives each tile the rows LASlib's -inside
# filter would keep (xmin <= X < xmax, ymin <= Y < ymax), in block order.
block_points <- LAS(data.table::data.table(
  X = c(0, 5, 10, 19.999, 20, 25, 40, 10, 30),
  Y = c(0, 5, 10, 19.999, 0, 25, 5, 20, 39.999),
  Z = 1
))
block_windows <- list(c(0, 0, 20, 20), c(20, 0, 40, 20), c(0, 20, 20, 40), c(20, 20, 40, 40))
block_rows <- tile_point_rows(block_points, block_windows)
for (index in seq_along(block_windows)) {
  window <- block_windows[[index]]
  expected <- which(block_points$X >= window[[1]] & block_points$X < window[[3]] &
                    block_points$Y >= window[[2]] & block_points$Y < window[[4]])
  stopifnot(identical(as.integer(block_rows[[index]]), as.integer(expected)))
}
stopifnot(identical(unname(lengths(block_rows)), c(4L, 1L, 1L, 2L)))

# Overfull CHM reads are replaced by their quarters, which tile the read
# exactly with shared edges.
split_parent <- lidR:::LAScluster(list(x = 30, y = 30), 60, 60, 0, lidR:::LIDRRECTANGLE,
  "cloud.laz", "parent", sf::st_crs(NA), list())
split_plan <- list()
split_plan[[cluster_key(split_parent)]] <- list(c(0, 0, 30, 30), c(30, 0, 30, 30),
                                                c(0, 30, 30, 30), c(30, 30, 30, 30))
split_children <- adaptive_chm_read_clusters(split_parent, 60, NULL, split_plan)
stopifnot(length(split_children) == 4L)
stopifnot(identical(vapply(split_children, cluster_key, character(1)), c(
  "0.000000 0.000000 30.000000 30.000000", "30.000000 0.000000 60.000000 30.000000",
  "0.000000 30.000000 30.000000 60.000000", "30.000000 30.000000 60.000000 60.000000")))
stopifnot(length(adaptive_chm_read_clusters(split_parent, 60, NULL, list())) == 1L)
