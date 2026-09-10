# Behavioral checks for the COPC adaptive chunk planner.  The planner is a pure
# seam: callers provide public-header facts and receive the runtime chunk policy
# that is also written to performance/provenance artifacts.

plan <- adaptive_chunk_plan(
  point_count = 400000000,
  xy_bounds = c(xmin = 0, ymin = 0, xmax = 50, ymax = 50),
  configured_chunk_size = 60,
  target_chunk_points = 50000000,
  is_copc = TRUE
)
stopifnot(
  isTRUE(plan$enabled),
  identical(plan$method, "mean_density"),
  identical(plan$target_points, 50000000),
  identical(plan$chunk_size_m, 17),
  identical(plan$alignment_m, 1),
  identical(plan$estimated_points_per_full_chunk, 46240000)
)

disabled <- adaptive_chunk_plan(
  400000000,
  c(xmin = 0, ymin = 0, xmax = 50, ymax = 50),
  60,
  0,
  TRUE
)
stopifnot(!disabled$enabled, disabled$chunk_size_m == 60)

plain_laz <- adaptive_chunk_plan(
  400000000,
  c(xmin = 0, ymin = 0, xmax = 50, ymax = 50),
  60,
  50000000,
  FALSE
)
stopifnot(!plain_laz$enabled, plain_laz$chunk_size_m == 60)

ceiling_limited <- adaptive_chunk_plan(
  1000000,
  c(xmin = 0, ymin = 0, xmax = 1000, ymax = 1000),
  120,
  50000000,
  TRUE
)
stopifnot(ceiling_limited$enabled, ceiling_limited$chunk_size_m == 60)

degenerate_rejected <- tryCatch(
  {
    adaptive_chunk_plan(
      100,
      c(xmin = 1, ymin = 0, xmax = 1, ymax = 10),
      60,
      50000000,
      TRUE
    )
    FALSE
  },
  error = function(error) grepl("positive XY area", conditionMessage(error))
)
stopifnot(degenerate_rejected)

# The target is deliberately approximate: a locally dense window is recorded
# and warned about, but it remains a successful observation rather than an
# analysis failure.
observation_directory <- tempfile("adaptive_observations_")
oversize_log <- capture.output(
  record_chunk_observation(
    observation_directory,
    "segment",
    50000001,
    50000000
  ),
  type = "message"
)
observation_summary <- summarize_chunk_observations(observation_directory)
stopifnot(
  any(grepl("adaptive_chunk_oversize", oversize_log, fixed = TRUE)),
  observation_summary$count == 1L,
  observation_summary$maximum_points == 50000001,
  observation_summary$oversized_count == 1L
)
unlink(observation_directory, recursive = TRUE, force = TRUE)

# Equal-height apexes from different spatial windows must follow the original
# source order rather than whichever worker result happens to be bound first.
segment_directory <- tempfile("adaptive_segment_tie_")
dir.create(segment_directory)
segment_piece <- function(apex_x, apex_order, voxel_x, voxel_first_order) {
  list(
    voxels = data.table(
      instance_dimension = "PredInstance",
      instance_id = 1,
      voxel_x = voxel_x,
      voxel_y = 1,
      voxel_z = 1,
      voxel_first_order = voxel_first_order
    ),
    layers = data.table(
      instance_dimension = "PredInstance",
      instance_id = 1,
      occupied_layer = 1
    ),
    apex = data.table(
      instance_dimension = "PredInstance",
      instance_id = 1,
      apex_x = apex_x,
      apex_y = 1,
      apex_z = 10,
      apex_order = apex_order
    )
  )
}
late_path <- file.path(segment_directory, "late.rds")
early_path <- file.path(segment_directory, "early.rds")
saveRDS(segment_piece(10, 100, 2, 100), late_path, compress = FALSE)
saveRDS(segment_piece(20, 50, 1, 50), early_path, compress = FALSE)
tie_store <- list(
  chunks = list(setNames(late_path, "1"), setNames(early_path, "1")),
  bucket_count = 1L
)
tie_parameters <- list(
  voxel_resolution = 0.2,
  minimum_tree_voxels = 1L,
  apex_minimum_height = 0,
  minimum_tree_thickness = 0,
  minimum_occupied_layers = 1L,
  memory_budget_gib = Inf
)
tie_result <- finalize_segment_store(tie_store, tie_parameters, "disk")
stopifnot(tie_result$apex_x[[1]] == 20)
unlink(segment_directory, recursive = TRUE, force = TRUE)

# One instance can cross many windows and place the same voxel in two pieces.
# The reducer must keep the globally earliest occurrence and feed PCA an order
# independent of worker completion or RDS enumeration order.
cross_window_piece <- function(voxels, orders) {
  list(
    voxels = data.table(
      instance_dimension = "PredInstance",
      instance_id = 9,
      voxel_x = voxels[, 1],
      voxel_y = voxels[, 2],
      voxel_z = voxels[, 3],
      voxel_first_order = orders
    ),
    layers = data.table(
      instance_dimension = "PredInstance",
      instance_id = 9,
      occupied_layer = unique(voxels[, 3])
    ),
    apex = data.table(
      instance_dimension = "PredInstance",
      instance_id = 9,
      apex_x = 0,
      apex_y = 0,
      apex_z = 8,
      apex_order = min(orders)
    )
  )
}
piece_a <- cross_window_piece(
  rbind(c(4, 1, 3), c(1, 1, 1), c(2, 4, 2)),
  c(40, 20, 30)
)
piece_b <- cross_window_piece(
  rbind(c(3, 2, 5), c(1, 1, 1)),
  c(10, 5)
)
merged_forward <- merge_segment_pieces(list(piece_a, piece_b))
merged_reverse <- merge_segment_pieces(list(piece_b, piece_a))
stopifnot(
  identical(merged_forward$voxels, merged_reverse$voxels),
  nrow(merged_forward$voxels) == 4L,
  identical(merged_forward$voxels$voxel_first_order, c(5, 10, 30, 40))
)
pca_forward <- finalize_segments(merged_forward, tie_parameters)
pca_reverse <- finalize_segments(merged_reverse, tie_parameters)
stopifnot(identical(pca_forward, pca_reverse))

# Both adaptive scientific stages must share an absolute one-metre-compatible
# origin so a smaller window cannot split a voxel or CHM cell differently.
chm_catalog_body <- paste(deparse(body(write_global_chm)), collapse = "\n")
segment_catalog_body <- paste(deparse(body(accumulate_segments)), collapse = "\n")
stopifnot(
  grepl("opt_chunk_alignment", chm_catalog_body, fixed = TRUE),
  grepl("opt_chunk_alignment", segment_catalog_body, fixed = TRUE)
)

# Adaptive window selection must not replace the baseline's numeric-resolution
# CHM rasterization with a fixed per-window template.  The latter changes only
# a few NA/value cells at exact window edges, but valid_updated parity is an
# exact contract rather than a tolerance-based comparison.
chm_chunk_body <- paste(deparse(body(chm_chunk)), collapse = "\n")
stopifnot(
  grepl("res = chm_resolution", chm_chunk_body, fixed = TRUE),
  !grepl("chunk_template", chm_chunk_body, fixed = TRUE)
)

# Sorting a very dense window must not materialize a second complete LAS table.
# The order vector is unavoidable, but data.table's in-place reordering avoids
# the multi-gigabyte row-subset copy that previously pushed dataset 157 over its
# 64 GiB internal memory gate after otherwise successful processing.
order_restore_body <- paste(
  deparse(body(restore_original_point_order)),
  collapse = "\n"
)
stopifnot(
  grepl("setorderv", order_restore_body, fixed = TRUE),
  !grepl("[order(point_order)", order_restore_body, fixed = TRUE),
  grepl("is.unsorted", order_restore_body, fixed = TRUE),
  !grepl("anyDuplicated", order_restore_body, fixed = TRUE)
)

# Allocation-safe validation must still reject a duplicate source ordinal.
duplicate_order_las <- lidR::LAS(data.frame(
  X = c(0, 1),
  Y = c(0, 1),
  Z = c(1, 2),
  OriginalPointIndex = c(7, 7)
))
duplicate_order_rejected <- tryCatch(
  {
    restore_original_point_order(
      duplicate_order_las,
      "OriginalPointIndex",
      "adaptive test"
    )
    FALSE
  },
  error = function(error) grepl("duplicate", conditionMessage(error))
)
stopifnot(duplicate_order_rejected)

# DTM geometry guards must use bounded block scans rather than materializing a
# second unique XY table for a dense catalog chunk.  Exercise points spanning
# scan blocks, exact collinearity, ground-only filtering, and too-few points.
stopifnot(
  !xy_has_nondegenerate_points(c(0, 1), c(0, 1), block_size = 1L),
  !xy_has_nondegenerate_points(c(0, 1, 2), c(0, 1, 2), block_size = 1L),
  xy_has_nondegenerate_points(c(0, 1, 2), c(0, 1, 0), block_size = 1L),
  !xy_has_nondegenerate_points(
    c(0, 1, 2, 3), c(0, 1, 0, 2), c(1L, 1L, 1L, 2L), 2L,
    block_size = 2L
  ),
  xy_has_nondegenerate_points(
    c(0, 1, 2, 3), c(0, 1, 0, 2), c(2L, 2L, 2L, 1L), 2L,
    block_size = 2L
  )
)
geometry_guard_body <- paste(
  deparse(body(xy_has_nondegenerate_points)),
  collapse = "\n"
)
stopifnot(
  !grepl("unique", geometry_guard_body, fixed = TRUE),
  !grepl("as.matrix", geometry_guard_body, fixed = TRUE)
)

# Once DTM order is restored, the order key itself is scientifically unused by
# PTD/TIN and must be released before those high-water-memory operations.
dtm_body <- paste(deparse(body(dtm_chunk)), collapse = "\n")
stopifnot(grepl("release_point_order_dimension", dtm_body, fixed = TRUE))

# One LAScatalog worker must not fan the very large source-order radix sort out
# to several data.table buffers.  Scientific operations and their parameters
# are unchanged; this is the memory-bound execution policy for dense windows.
stopifnot(identical(DATA_TABLE_THREADS, 1L))
