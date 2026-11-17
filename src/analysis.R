suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(future)
  library(lidR)
  library(sf)
  library(terra)
})

# Runtime invariants chosen for the validated production image. Catalog
# workers use one lidR thread; segment aggregation uses a small, fixed
# data.table pool and deterministic disk buckets to bound peak memory.  A single
# table thread prevents dense OriginalPointIndex radix sorts from allocating one
# large scratch buffer per thread inside the single LAScatalog worker.
DATA_TABLE_THREADS <- 1L
SEGMENT_BUCKET_COUNT <- 64L
# The image's rlas reports larger counts as 64-bit numbers, and lidR only reads
# spatial windows, so a file may hold up to 2^53 points.
MAX_LIDR_POINT_COUNT <- 2^53
DTM_CATALOG_SELECTION <- "xyz"
SEGMENT_CATALOG_BASE_SELECTION <- "xyzrn"
SEGMENT_SELECTION <- "xyz"

# COPC stores points in spatial hierarchy order, while the original forest-structure script
# implicitly resolves several equal-value ties by input record order.  Sort each
# independently read chunk by the persisted source ordinal before applying any
# order-sensitive DTM, CHM, segmentation, or tile operation.  Plain LAS/LAZ runs
# pass NULL and retain their already-original order without extra work.
restore_original_point_order <- function(las, point_order_dimension = NULL,
                                         phase = "analysis") {
  # Preserve the historical path when no key was requested, and avoid touching
  # genuinely empty chunks returned by LAScatalog edge queries.
  if (is.null(point_order_dimension) || is.empty(las)) return(las)
  # A missing key means selective reading or COPC construction violated the
  # canonical input contract; continuing would produce plausible but wrong ties.
  if (!point_order_dimension %in% names(las@data)) {
    stop(sprintf(
      "%s did not load required point-order dimension %s",
      phase,
      point_order_dimension
    ))
  }
  # Missing ordinals make deterministic reconstruction impossible.  Do not bind
  # the complete column separately: an extra reference can prevent in-place work
  # on very dense windows.
  if (anyNA(las@data[[point_order_dimension]])) {
    stop(sprintf(
      "%s received missing or duplicate %s values",
      phase,
      point_order_dimension
    ))
  }
  # Reorder the complete point table in place so all dimensions remain aligned.
  # A row subset would allocate a second full LAS table; on dense windows that
  # copy alone can consume several GiB and breach the internal memory budget.
  data.table::setorderv(las@data, point_order_dimension)
  # Once sorted, strict monotonicity detects duplicate ordinals in one linear C
  # scan without the large hash table allocated by anyDuplicated().
  if (is.unsorted(las@data[[point_order_dimension]], strictly = TRUE)) {
    stop(sprintf(
      "%s received missing or duplicate %s values",
      phase,
      point_order_dimension
    ))
  }
  # Return the same LAS object shape expected by existing scientific functions.
  las
}

parameter_or <- function(parameters, name, default) {
  value <- parameters[[name]]
  if (is.null(value)) default else value
}

validated_catalog_worker_count <- function(value, thread_budget) {
  workers <- suppressWarnings(as.integer(value))
  if (length(workers) != 1L || is.na(workers) ||
      suppressWarnings(as.numeric(value)) != workers || workers < 1L || workers > 10L) {
    stop("RASTER_CATALOG_WORKERS must be an integer from 1 to 10")
  }
  if (workers > thread_budget) {
    stop("RASTER_CATALOG_WORKERS exceeds RASTER_THREADS")
  }
  workers
}

current_rss_bytes <- function() {
  status <- tryCatch(
    readLines("/proc/self/status", warn = FALSE),
    error = function(error) character()
  )
  line <- grep("^VmRSS:", status, value = TRUE)
  if (length(line) != 1L) return(0)
  kib <- suppressWarnings(as.numeric(sub(
    "^VmRSS:[[:space:]]*([0-9]+).*$",
    "\\1",
    line
  )))
  if (!is.finite(kib)) 0 else kib * 1024
}

assert_memory_budget <- function(required_bytes, parameters, operation,
                                 single_segment = FALSE) {
  budget_gib <- parameter_or(parameters, "memory_budget_gib", Inf)
  if (!is.finite(budget_gib)) return(invisible(TRUE))
  rss <- current_rss_bytes()
  projected <- rss + max(0, required_bytes)
  if (projected <= budget_gib * 1024^3) return(invisible(TRUE))
  flag <- if (single_segment) {
    "single_segment_exceeds_memory_budget"
  } else {
    "memory_budget_exceeded"
  }
  stop(sprintf(
    paste0(
      "RASTER_FLAG %s: %s requires an estimated %.3f GiB ",
      "at %.3f GiB RSS with a %.3f GiB internal budget"
    ),
    flag,
    operation,
    required_bytes / 1024^3,
    rss / 1024^3,
    budget_gib
  ))
}

# A read with selector "0" (every ExtraByte: how fields after the ninth are
# reached, patches/rlas-1.8.6-all-extra-bytes.patch) also brings ExtraBytes the
# phase does not use. Drop them right after the read, so the phase's memory and
# its selective-read check see only the requested fields. Standard LAS fields
# are never dropped here.
keep_requested_extra_bytes <- function(las, keep) {
  descriptions <- las@header@VLR$Extra_Bytes[["Extra Bytes Description"]]
  if (is.null(descriptions) || !length(descriptions)) return(las)
  extra <- vapply(descriptions, function(description) as.character(description$name), character(1))
  for (name in intersect(setdiff(names(las@data), keep), extra)) {
    las <- lidR::remove_lasattribute(las, name)
  }
  las
}

assert_loaded_dimensions <- function(las, allowed, phase) {
  present <- names(las@data)
  unexpected <- setdiff(present, allowed)
  missing <- setdiff(c("X", "Y", "Z"), present)
  if (length(unexpected) > 0L || length(missing) > 0L) {
    stop(sprintf(
      paste0(
        "RASTER_FLAG selective_read_violation: %s loaded [%s]; ",
        "allowed [%s]"
      ),
      phase,
      paste(present, collapse = ","),
      paste(allowed, collapse = ",")
    ))
  }
  invisible(TRUE)
}

configure_catalog_worker <- function() {
  lidR::set_lidr_threads(1L)
  data.table::setDTthreads(DATA_TABLE_THREADS)
  invisible(NULL)
}

read_las_uint_le <- function(connection, bytes) {
  raw_value <- readBin(connection, what = "raw", n = bytes)
  if (length(raw_value) != bytes) stop("LAS/LAZ header is truncated")
  sum(as.numeric(raw_value) * 256^(seq_len(bytes) - 1L))
}

point_count_from_las_header <- function(path) {
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  signature <- rawToChar(readBin(connection, what = "raw", n = 4L))
  if (!identical(signature, "LASF")) {
    stop("point cloud does not have a LAS/LAZ header")
  }

  seek(connection, where = 24L, origin = "start")
  version <- as.integer(readBin(connection, what = "raw", n = 2L))
  if (length(version) != 2L) stop("LAS/LAZ header is truncated")

  seek(connection, where = 107L, origin = "start")
  legacy_count <- read_las_uint_le(connection, 4L)
  is_las_14 <- version[[1]] > 1L ||
    (version[[1]] == 1L && version[[2]] >= 4L)
  if (!is_las_14) return(legacy_count)

  seek(connection, where = 247L, origin = "start")
  extended_count <- read_las_uint_le(connection, 8L)
  if (extended_count > 0) extended_count else legacy_count
}

assert_lidr_point_count_supported <- function(point_count) {
  if (!is.finite(point_count) || point_count < 0) {
    stop("LAS/LAZ header contains an invalid point count")
  }
  if (point_count > MAX_LIDR_POINT_COUNT) {
    stop(sprintf(
      paste(
        "RASTER_FLAG unsupported_lidr_point_count:",
        "single point-cloud file contains %.0f points; lidR/R supports at most %.0f"
      ),
      point_count,
      MAX_LIDR_POINT_COUNT
    ))
  }
  invisible(point_count)
}

# Translate a point-count target into the square metre-based window understood
# by LAScatalog.  COPC makes those spatial reads efficient; plain LAS/LAZ keeps
# its established fixed-size behavior so this runtime control cannot silently
# change the legacy scientific path.
adaptive_chunk_plan <- function(point_count, xy_bounds, configured_chunk_size,
                                target_chunk_points, is_copc,
                                alignment_m = 1) {
  if (!is.finite(configured_chunk_size) || configured_chunk_size <= 0) {
    stop("configured chunk size must be greater than zero")
  }
  if (!is.finite(target_chunk_points) || target_chunk_points < 0) {
    stop("target chunk points must be zero or greater")
  }
  disabled <- !isTRUE(is_copc) || target_chunk_points == 0
  if (disabled) {
    return(list(
      enabled = FALSE,
      method = "fixed_metres",
      target_points = target_chunk_points,
      alignment_m = alignment_m,
      chunk_size_m = configured_chunk_size,
      estimated_points_per_full_chunk = NA_real_
    ))
  }
  required_names <- c("xmin", "ymin", "xmax", "ymax")
  if (!all(required_names %in% names(xy_bounds))) {
    stop("adaptive chunk planning requires named XY bounds")
  }
  width <- as.numeric(xy_bounds[["xmax"]] - xy_bounds[["xmin"]])
  height <- as.numeric(xy_bounds[["ymax"]] - xy_bounds[["ymin"]])
  area <- width * height
  if (!is.finite(area) || area <= 0) {
    stop("adaptive chunk planning requires a positive XY area")
  }
  if (!is.finite(point_count) || point_count < 0) {
    stop("adaptive chunk planning requires a non-negative point count")
  }
  if (!is.finite(alignment_m) || alignment_m <= 0) {
    stop("adaptive chunk alignment must be greater than zero")
  }
  density <- point_count / area
  ideal_size <- if (density == 0) {
    configured_chunk_size
  } else {
    sqrt(target_chunk_points / density)
  }
  # Sixty metres is the scientifically validated upper bound.  A smaller
  # --chunk-size remains an additional operator-supplied ceiling, while a
  # larger value cannot silently widen the established window geometry.
  maximum_size <- min(configured_chunk_size, 60)
  aligned_size <- floor(min(ideal_size, maximum_size) / alignment_m) *
    alignment_m
  chunk_size <- max(alignment_m, aligned_size)
  list(
    enabled = TRUE,
    method = "mean_density",
    target_points = target_chunk_points,
    alignment_m = alignment_m,
    point_count = point_count,
    width_m = width,
    height_m = height,
    area_m2 = area,
    mean_density_points_m2 = density,
    ideal_chunk_size_m = ideal_size,
    configured_maximum_chunk_size_m = configured_chunk_size,
    effective_maximum_chunk_size_m = maximum_size,
    chunk_size_m = chunk_size,
    estimated_points_per_full_chunk = density * chunk_size^2
  )
}

script_argument <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
script_directory <- if (length(script_argument)) {
  dirname(normalizePath(sub("^--file=", "", script_argument[[1]])))
} else {
  getwd()
}
source(file.path(script_directory, "tile_scheduler.R"))
source(file.path(script_directory, "resource_scheduler.R"))
source(file.path(script_directory, "catalog_scheduler.R"))
source(file.path(script_directory, "voxel_counts.R"))
source(file.path(script_directory, "artifact_publication.R"))
source(file.path(script_directory, "denoising.R"))

# CLI and artifact naming ---------------------------------------------------

artifact_paths <- function(output_dir, dataset_id) {
  prefix <- file.path(output_dir, as.character(dataset_id))
  list(
    prefix = prefix,
    dtm = paste0(prefix, "_dtm.tif"),
    chm = paste0(prefix, "_chm.tif"),
    dsm = paste0(prefix, "_dsm.tif"),
    tiles_png = paste0(prefix, "_tiles.png"),
    performance = paste0(prefix, "_performance.csv")
  )
}

dimension_artifact_paths <- function(paths, instance_dimension = NULL) {
  suffix <- if (is.null(instance_dimension)) {
    ""
  } else {
    safe <- gsub("[^A-Za-z0-9_-]", "_", instance_dimension)
    if (!nzchar(safe)) stop("Instance Dimension cannot form a safe artifact name")
    paste0("_", safe)
  }
  list(
    results = paste0(paths$prefix, suffix, "_results.csv"),
    segment_diagnostics = paste0(
      paths$prefix, suffix, "_segment_diagnostics.csv"
    ),
    tiles_geojson = paste0(paths$prefix, suffix, "_tiles.geojson")
  )
}

# The stages, one file each.
source(file.path(script_directory, "analysis_grid.R"))
source(file.path(script_directory, "analysis_rasters.R"))
source(file.path(script_directory, "analysis_reduction.R"))
source(file.path(script_directory, "analysis_segments.R"))
source(file.path(script_directory, "analysis_tiles.R"))
source(file.path(script_directory, "analysis_run.R"))
