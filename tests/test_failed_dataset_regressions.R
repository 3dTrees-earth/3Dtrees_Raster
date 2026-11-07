stopifnot(segment_selection_from_names(
  c(paste0("extra_", 1:9), "PredInstance"),
  "PredInstance"
) == "xyzrn0")

stopifnot(identical(SEGMENT_CATALOG_BASE_SELECTION, "xyzrn"))
stopifnot(segment_selection_from_names(
  c("PredInstance", "unrelated"),
  "PredInstance"
) == "xyzrn1")
stopifnot(segment_selection_from_names(
  c("unrelated"),
  "PredInstance"
) == "xyzrn")

stopifnot(identical(DTM_CATALOG_SELECTION, "xyz"))
fine_scale_las <- LAS(data.frame(
  X = c(226990.0000001, 226991.0000002, 226990.5000003),
  Y = c(378990.0000001, 378990.5000002, 378991.0000003),
  Z = c(357, 358, 359),
  Classification = rep(2L, 3L)
))
fine_scale_las@header@PHB[["X scale factor"]] <- 1e-7
fine_scale_las@header@PHB[["Y scale factor"]] <- 1e-7
# The explicit fallback must execute lidR's legacy floating-point triangulation
# directly and leave the source LAS quantization metadata unchanged.
fine_scale_header_before <- fine_scale_las@header@PHB
legacy_surface <- rasterize_terrain(
  fine_scale_las,
  res = 1,
  algorithm = legacy_tin()
)
stopifnot(inherits(legacy_surface, "SpatRaster"))
stopifnot(identical(fine_scale_las@header@PHB, fine_scale_header_before))

# With unequal XY scales, lidR selects its own legacy branch before attempting
# integer triangulation. Its raster must match the explicit retry implementation.
automatic_slow_las <- fine_scale_las
automatic_slow_las@header@PHB[["Y scale factor"]] <- 1e-6
automatic_slow_surface <- suppressMessages(rasterize_terrain(
  automatic_slow_las,
  res = 0.25,
  algorithm = tin()
))
explicit_slow_surface <- rasterize_terrain(
  fine_scale_las,
  res = 0.25,
  algorithm = legacy_tin()
)
stopifnot(identical(
  terra::values(explicit_slow_surface, mat = FALSE),
  terra::values(automatic_slow_surface, mat = FALSE)
))
dtm_body <- paste(deparse(body(write_global_dtm)), collapse = "\n")
stopifnot(grepl(
  "catalog_selection_for_dimensions",
  dtm_body,
  fixed = TRUE
))
stopifnot(grepl("DTM_CATALOG_SELECTION", dtm_body, fixed = TRUE))
stopifnot(grepl(
  "complete_chunk_raster_paths(results, \"DTM\")",
  dtm_body,
  fixed = TRUE
))
stopifnot(!grepl("sort(", dtm_body, fixed = TRUE))
stopifnot(!grepl(
  "populated_or_all_empty_chunk_raster_paths(results, \"DTM\")",
  dtm_body,
  fixed = TRUE
))

chm_body <- paste(deparse(body(write_global_chm)), collapse = "\n")
stopifnot(grepl(
  paste0(
    "build_chunk_virtual_raster\\(\\s*results,\\s*work_directory,",
    "\\s*\"CHM\",\\s*source_crs\\s*\\)"
  ),
  chm_body,
  perl = TRUE
))
stopifnot(grepl(
  "stream_virtual_raster(covering, output_path, \"CHM\", source_crs)",
  chm_body,
  fixed = TRUE
))

segment_body <- paste(deparse(body(segment_chunk)), collapse = "\n")
stopifnot(grepl("lapply(selected", segment_body, fixed = TRUE))
stopifnot(!grepl("rbindlist(lapply(selected", segment_body, fixed = TRUE))

point_bounds <- c(
  xmin = -1.1641532182693481e-10,
  ymin = 4294967.296,
  xmax = 999.9989999999525,
  ymax = 4295967.2949999999
)
covering_extent <- terra::ext(0, 1000, 4294967, 4295968)
stopifnot(raster_extent_covers_point_bounds(
  covering_extent,
  point_bounds,
  1
))
stopifnot(!raster_extent_covers_point_bounds(
  terra::ext(0.01, 1000, 4294967, 4295968),
  point_bounds,
  1
))

stopifnot(assert_lidr_point_count_supported(MAX_LIDR_POINT_COUNT) ==
  MAX_LIDR_POINT_COUNT)
unsupported <- tryCatch(
  {
    assert_lidr_point_count_supported(MAX_LIDR_POINT_COUNT + 1)
    FALSE
  },
  error = function(error) grepl(
    "RASTER_FLAG unsupported_lidr_point_count",
    conditionMessage(error),
    fixed = TRUE
  )
)
stopifnot(unsupported)

complete <- list("one", "two")
stopifnot(identical(assert_catalog_completed(complete, "test"), complete))
incomplete <- tryCatch(
  {
    assert_catalog_completed(list("one", NULL), "test")
    FALSE
  },
  error = function(error) grepl("incomplete chunks", conditionMessage(error))
)
stopifnot(incomplete)

collapsed <- c(xmin = 10, ymin = 20, xmax = 10, ymax = 20)
collapsed <- ensure_nonempty_grid_bounds(collapsed, 1)
stopifnot(collapsed[["xmax"]] > collapsed[["xmin"]])
stopifnot(collapsed[["ymax"]] > collapsed[["ymin"]])
