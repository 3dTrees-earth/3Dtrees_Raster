# Run after sourcing run.R. Exercise real catalog reads, normalization, disk
# summaries and reduction, including overlapping IDs and source-order apex ties.
local({
  root <- tempfile("instance_read_groups_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  path <- file.path(root, "instances.las")
  other_path <- file.path(root, "other_instances.las")
  alias <- file.path(root, "alias.las")
  dimensions <- c("PredInstance_SAT", "PredInstance_FM")
  cloud <- LAS(data.frame(
    X = 1000 + c(0.1, 9.9, 10, 10.1, 19.9, 2, 5, 12, 15, 9.99, 10.01, 18),
    Y = 2000 + c(0.1, 2, 2, 5, 9.9, 12, 15, 12, 15, 18, 18, 19.9),
    Z = 100 + c(1, 8, 8, 4, 2, 6, 7, 5, 3, 0, -1, 71)
  ))
  sf::st_crs(cloud) <- 32632
  cloud <- add_lasattribute(cloud, c(30, 20, 10, 40:48),
                           "OriginalPointIndex", "source order")
  cloud <- add_lasattribute(cloud, c(rep(1, 5), rep(2, 4), 0, 1, 1),
                           dimensions[[1]], "SAT ID")
  cloud <- add_lasattribute(cloud, c(rep(2, 5), rep(1, 4), 0, 2, 2),
                           dimensions[[2]], "FM ID")
  writeLAS(cloud, path)
  cloud@data[, PredInstance_FM := PredInstance_FM + 100]
  writeLAS(cloud, other_path)
  stopifnot(file.symlink(path, alias))

  grouped <- instance_read_groups(dimensions, path,
                                  list(PredInstance_FM = alias))
  separate <- instance_read_groups(dimensions, path, sequential = TRUE)
  different <- instance_read_groups(dimensions, path,
                                    list(PredInstance_FM = other_path))
  stopifnot(length(grouped) == 1L,
            identical(grouped[[1]]$dimensions, dimensions),
            length(separate) == 2L, length(different) == 2L,
            identical(different[[2]]$source, other_path))

  dtm <- terra::rast(xmin = 999, xmax = 1021, ymin = 1999, ymax = 2021, res = 1)
  terra::crs(dtm) <- "EPSG:32632"
  terra::values(dtm) <- 100
  dtm_path <- file.path(root, "dtm.tif")
  terra::writeRaster(dtm, dtm_path)
  parameters <- list(voxel_resolution = 0.2, minimum_tree_voxels = 1L,
    apex_minimum_height = 0, minimum_tree_thickness = 0,
    minimum_occupied_layers = 1L, memory_budget_gib = Inf)
  counts <- new.env(parent = emptyenv())
  counts$read <- counts$normalized <- 0L
  runtime <- environment(segment_chunk)
  original_read <- get("readLAS", envir = runtime)
  original_normalize <- get("normalize_height", envir = runtime)
  had_read <- exists("readLAS", envir = runtime, inherits = FALSE)
  had_normalize <- exists("normalize_height", envir = runtime, inherits = FALSE)
  on.exit({
    if (had_read) assign("readLAS", original_read, envir = runtime) else
      rm("readLAS", envir = runtime)
    if (had_normalize) assign("normalize_height", original_normalize, envir = runtime) else
      rm("normalize_height", envir = runtime)
  }, add = TRUE)
  assign("readLAS", function(...) {
    value <- original_read(...)
    counts$read <- counts$read + nrow(value@data)
    value
  }, envir = runtime)
  assign("normalize_height", function(las, ...) {
    counts$normalized <- counts$normalized + nrow(las@data)
    original_normalize(las, ...)
  }, envir = runtime)

  run_groups <- function(groups, label) {
    counts$read <- counts$normalized <- 0L
    result <- list()
    for (index in seq_along(groups)) {
      group <- groups[[index]]
      directory <- file.path(root, paste0(label, "_", index))
      store <- accumulate_segments(group$source, dtm_path, 10, 0.2, 70,
        group$dimensions, 1L, directory, bucket_count = 2L,
        base_selection = JULIA_MEMORY_SAFE_SEGMENT_SELECTION,
        point_order_dimension = "OriginalPointIndex")
      for (dimension in group$dimensions) {
        single <- segment_store_for_dimension(store, dimension)
        # Every selected RDS contains only this dimension, including ID bucket
        # collisions between SAT and FM. Ownership selection does not read RDS.
        for (piece_path in unlist(single$chunks, use.names = FALSE)) {
          piece <- readRDS(piece_path)
          stopifnot(all(vapply(piece, function(table)
            all(table$instance_dimension == dimension), logical(1))))
        }
        result[dimension] <- list(finalize_segment_store(single, parameters))
      }
    }
    list(segments = result, read = counts$read, normalized = counts$normalized)
  }
  fused <- run_groups(grouped, "fused")
  independent <- run_groups(separate, "independent")
  stopifnot(identical(fused$segments, independent$segments),
            fused$read > 0, independent$read == 2 * fused$read,
            independent$normalized == 2 * fused$normalized)
  sat <- fused$segments[[dimensions[[1]]]]
  fm <- fused$segments[[dimensions[[2]]]]
  # The equal-height apex at X=1010 has the earlier persisted ordinal even
  # though it was written later and falls into a different catalog chunk.
  stopifnot(sat$apex_x[sat$instance_id == 1] == 1010,
            fm$apex_x[fm$instance_id == 2] == 1010,
            identical(sat$instance_id, c(1, 2)),
            identical(fm$instance_id, c(1, 2)))
  alternate <- run_groups(different, "different")
  stopifnot(identical(alternate$segments[[dimensions[[1]]]], sat),
            all(alternate$segments[[dimensions[[2]]]]$instance_id >= 100))
  # A present dimension with no positive IDs has a valid empty store. Its empty
  # bucket metadata must not remove its output name or borrow the other IDs.
  cloud@data[, PredInstance_FM := 0]
  empty_fm_path <- file.path(root, "empty_fm.las")
  writeLAS(cloud, empty_fm_path)
  empty_fm <- run_groups(instance_read_groups(dimensions, empty_fm_path), "empty_fm")
  stopifnot(identical(names(empty_fm$segments), dimensions),
    identical(empty_fm$segments[[dimensions[[1]]]], sat),
    is.null(empty_fm$segments[[dimensions[[2]]]]))
  cloud@data[, PredInstance_SAT := 0]
  empty_all_path <- file.path(root, "empty_all.las")
  writeLAS(cloud, empty_all_path)
  empty_all <- run_groups(instance_read_groups(dimensions, empty_all_path), "empty_all")
  stopifnot(identical(names(empty_all$segments), dimensions),
            all(vapply(empty_all$segments, is.null, logical(1))))
  message("instance read groups: exact segment parity; shared reads and normalization halved")
})
