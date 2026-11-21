# Point reduction for dense clouds, chunk observations, the CHM and DSM chunk
# passes and the split of dense reads. Sourced by analysis.R; uses its
# constants and helpers.

# Point reduction -----------------------------------------------------------

# A DTM window holding more points than this is thinned for the ground filter,
# which crashes on about 640 M points and only needs the low ones.
dtm_thinning_points <- function() {
  as.numeric(Sys.getenv("RASTER_DTM_THINNING_POINTS", "3e8"))
}

# Below this mean density voxels hold too few points for reduction to pay off.
reduction_density_points_m2 <- function() {
  as.numeric(Sys.getenv("RASTER_REDUCTION_DENSITY", "1e4"))
}

dense_point_cloud <- function(point_cloud) {
  bounds <- point_cloud_xy_bounds(point_cloud)
  area <- (bounds[["xmax"]] - bounds[["xmin"]]) * (bounds[["ymax"]] - bounds[["ymin"]])
  area > 0 && point_count_from_las_header(point_cloud) / area >= reduction_density_points_m2()
}

run_point_reduction <- function(mode, source, output, parameters, ...) {
  script <- file.path(script_directory, "reduce_points.py")
  report <- paste0(output, ".json")
  denoising <- parameters$denoising
  arguments <- c(mode, "--input", shQuote(source), "--output", shQuote(output),
    "--report", shQuote(report),
    if (!is.null(parameters$point_order_dimension)) {
      c("--order-dimension", shQuote(parameters$point_order_dimension))
    },
    if (!is.null(denoising)) c("--rejected-directory", shQuote(denoising$directory)),
    ...)
  environment <- sprintf("RAYON_NUM_THREADS=%d POLARS_MAX_THREADS=%d",
                         parameters$threads, parameters$threads)
  status <- system2("env", c(strsplit(environment, " ")[[1L]], "python3", shQuote(script), arguments))
  if (!identical(status, 0L) || !file.exists(output)) {
    stop(sprintf("point reduction (%s) failed", mode))
  }
  if (!grepl("\\.copc\\.laz$", output)) rlas::writelax(output)
  result <- jsonlite::read_json(report)
  message(sprintf("Point reduction (%s): kept %.0f of %.0f points",
                  mode, result$kept_points, result$source_points))
  output
}

# The largest DTM window, from the COPC hierarchy or the point count of the file.
dtm_needs_thinning <- function(parameters) {
  if (is.null(parameters$memory_policy)) return(FALSE)
  catalog <- readLAScatalog(parameters$point_cloud)
  opt_chunk_size(catalog) <- parameters$dtm_chunk_size
  opt_chunk_buffer(catalog) <- parameters$dtm_buffer
  opt_progress(catalog) <- FALSE
  settings <- lidR:::engine_parse_options(list(raster_alignment = parameters$dtm_resolution))
  chunks <- lidR:::engine_chunks(catalog, settings$raster_alignment)
  counts <- task_memory_estimates(lapply(chunks, catalog_task_window),
                                  parameters$memory_policy, "DTM")$counts
  length(counts) > 0L && max(counts) > dtm_thinning_points()
}

# CHM, DSM, segments and tiles read only occupied voxels, so they get the same
# results from the original points reduce_points.py keeps.
reduced_analysis_cloud <- function(parameters, tiles, dtm_path, instance_dimensions,
                                   work_directory) {
  bounds <- lapply(seq_len(nrow(tiles)), function(index) sf::st_bbox(tiles[index, ]))
  # calculate_tile_metrics() hands these bounds to LASlib with ten decimals.
  edges <- function(names) {
    sort(unique(as.numeric(sprintf("%.10f", unlist(lapply(bounds, `[`, names))))))
  }
  edges_path <- file.path(work_directory, "tile_edges.json")
  jsonlite::write_json(list(x = edges(c("xmin", "xmax")), y = edges(c("ymin", "ymax"))),
                       edges_path, digits = NA)
  output <- if (is.null(parameters$point_order_dimension)) "analysis_points.laz" else
    "analysis_points.copc.laz"
  run_point_reduction(
    "voxels", parameters$point_cloud, file.path(work_directory, output),
    parameters, "--dtm", shQuote(dtm_path), "--tile-edges", shQuote(edges_path),
    "--voxel-resolution", parameters$voxel_resolution,
    "--maximum-height", if (is.finite(parameters$maximum_height)) parameters$maximum_height else "inf",
    unlist(lapply(instance_dimensions, function(name) c("--dimension", shQuote(name)))))
}

build_dtm_overlap_mosaic <- function(chunk_paths, work_directory,
                                     fallback_crs = NA_character_) {
  reference_crs <- chunk_reference_crs(chunk_paths, "DTM")
  if ((is.na(reference_crs) || !nzchar(reference_crs)) &&
      !is.na(fallback_crs) && nzchar(fallback_crs)) {
    reference_crs <- fallback_crs
  }
  mosaic_path <- file.path(work_directory, "dtm_mosaic.tif")
  mosaic <- terra::mosaic(
    terra::sprc(lapply(chunk_paths, terra::rast)),
    filename = mosaic_path,
    overwrite = TRUE,
    wopt = list(
      datatype = "FLT8S",
      gdal = c("TILED=YES", "COMPRESS=LZW", "BIGTIFF=IF_SAFER")
    )
  )
  if (terra::inMemory(mosaic) || !file.exists(mosaic_path)) {
    stop("global DTM overlap mosaic was not written disk-to-disk")
  }

  mosaic_vrt_path <- file.path(work_directory, "dtm_mosaic.vrt")
  arguments <- c("-overwrite")
  if (!is.na(reference_crs) && nzchar(reference_crs)) {
    arguments <- c(arguments, "-a_srs", shQuote(reference_crs))
  }
  arguments <- c(
    arguments,
    shQuote(mosaic_vrt_path),
    shQuote(mosaic_path)
  )
  status <- system2(
    "gdalbuildvrt",
    arguments
  )
  if (!identical(status, 0L) || !file.exists(mosaic_vrt_path)) {
    stop("could not create the global DTM publication VRT")
  }
  result <- terra::rast(mosaic_vrt_path)
  if (!is.na(reference_crs) && nzchar(reference_crs)) {
    actual_crs <- terra::crs(result)
    if (is.na(actual_crs) || !nzchar(actual_crs)) {
      stop("global DTM publication VRT did not persist its reference CRS")
    }
  }
  result
}

# Persist tiny per-window counters independently from scientific artifacts so
# worker processes can report observed point pressure without sharing memory.
record_chunk_observation <- function(directory, phase, point_count,
                                     target_chunk_points) {
  if (is.null(directory)) return(invisible(NULL))
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  oversized <- is.finite(target_chunk_points) && target_chunk_points > 0 &&
    point_count > target_chunk_points
  message(sprintf(
    "ADAPTIVE_CHUNK_OBSERVED phase=%s points=%.0f target_points=%.0f oversized=%s",
    phase,
    point_count,
    target_chunk_points,
    tolower(as.character(oversized))
  ))
  if (oversized) {
    message(sprintf(
      paste0(
        "RASTER_FLAG adaptive_chunk_oversize: %s chunk contained ",
        "%.0f points above the approximate %.0f-point target"
      ),
      phase,
      point_count,
      target_chunk_points
    ))
  }
  path <- tempfile(
    # Include the worker PID so the optional two-worker compatibility mode
    # cannot collide even if forked workers inherit identical RNG state.
    pattern = sprintf("%s_chunk_%d_", phase, Sys.getpid()),
    tmpdir = directory,
    fileext = ".rds"
  )
  saveRDS(
    list(point_count = point_count, oversized = oversized, maximum_voxel_height = NA_real_),
    path,
    compress = FALSE
  )
  invisible(path)
}

record_chunk_height <- function(path, heights) {
  if (is.null(path)) return(invisible(NULL))
  observation <- readRDS(path)
  observation$maximum_voxel_height <- max(heights)
  saveRDS(observation, path, compress = FALSE)
  invisible(path)
}

maximum_observed_height <- function(observations) {
  heights <- vapply(observations, function(value)
    parameter_or(value, "maximum_voxel_height", NA_real_), numeric(1))
  heights <- heights[is.finite(heights)]
  if (length(heights)) max(heights) else NA_real_
}

summarize_chunk_observations <- function(directory) {
  paths <- if (is.null(directory) || !dir.exists(directory)) {
    character()
  } else {
    list.files(directory, pattern = "\\.rds$", full.names = TRUE)
  }
  if (length(paths) == 0L) {
    return(list(count = 0L, minimum_points = NA_real_, maximum_points = NA_real_,
                oversized_count = 0L, maximum_voxel_height = NA_real_))
  }
  observations <- lapply(paths, readRDS)
  points <- vapply(observations, `[[`, numeric(1), "point_count")
  oversized <- vapply(observations, `[[`, logical(1), "oversized")
  list(
    count = length(paths),
    minimum_points = min(points),
    maximum_points = max(points),
    oversized_count = sum(oversized),
    maximum_voxel_height = maximum_observed_height(observations)
  )
}

combine_chunk_observations <- function(summaries) {
  summaries <- Filter(function(summary) summary$count > 0L, summaries)
  if (length(summaries) == 0L) return(summarize_chunk_observations(NULL))
  list(
    count = sum(vapply(summaries, `[[`, integer(1), "count")),
    minimum_points = min(vapply(summaries, `[[`, numeric(1), "minimum_points")),
    maximum_points = max(vapply(summaries, `[[`, numeric(1), "maximum_points")),
    oversized_count = sum(vapply(summaries, `[[`, integer(1), "oversized_count")),
    maximum_voxel_height = maximum_observed_height(summaries)
  )
}

# DSM: the highest absolute point Z (`Zref`, kept by normalize_height) per
# voxel column, keyed by the exact column centres voxelize_points() gives the
# CHM. Rasterized on the CHM chunk raster as template, every DSM cell is built
# from the same points as its CHM cell and the two share grid and empty cells.
dsm_column_maxima <- function(data, voxel_resolution) {
  data[, .(Z = max(Zref)), by = .(
    X = lidR:::f_grid(X, voxel_resolution, 0),
    Y = lidR:::f_grid(Y, voxel_resolution, 0)
  )]
}

dsm_chunk_raster <- function(chm_raster, maxima, header, crs) {
  if (is.null(maxima) || nrow(maxima) == 0L) {
    dsm <- terra::rast(chm_raster)
    terra::values(dsm) <- NA_real_
    return(dsm)
  }
  compact <- LAS(maxima, header = header, crs = crs, check = FALSE)
  rasterize_canopy(compact, res = chm_raster, algorithm = p2r())
}

# The DSM chunk sits beside the CHM chunk lidR writes from the returned raster.
dsm_chunk_path <- function(chm_chunk_path) {
  file.path(dirname(chm_chunk_path), sub("^chm_", "dsm_", basename(chm_chunk_path)))
}

write_dsm_chunk <- function(chunk, chm_raster, maxima = NULL, header = NULL) {
  if (is.null(chm_raster)) return(invisible(NULL))
  dsm <- dsm_chunk_raster(chm_raster, maxima, header, sf::st_crs(chunk))
  terra::writeRaster(dsm, dsm_chunk_path(paste0(chunk@save, ".tif")), overwrite = TRUE)
  invisible(NULL)
}

# Process one CHM chunk after restoring the source record sequence when COPC is
# active.  This keeps order-dependent voxel/canopy ties reproducible.
# With `dsm = TRUE` the same points also give the chunk's DSM, written beside it.
# Ground points the interpolated DTM lies slightly above (paths, ditches) get height 0,
# so bare ground reads as 0 m canopy instead of an empty CHM and DSM cell (3DT-3028).
# Points more than half a voxel below the DTM are still dropped.
ground_to_zero <- function(las, voxel_resolution) {
  below <- is.finite(las@data$Z) & las@data$Z < 0 & las@data$Z >= -voxel_resolution / 2
  if (any(below)) las@data$Z[below] <- 0
  las
}

chm_chunk <- function(chunk, dtm_path, voxel_resolution, maximum_height,
                      chm_resolution, point_order_dimension = NULL,
                      observation_directory = NULL,
                      target_chunk_points = 0, denoising = NULL,
                      boundary_read_adjustment = NULL, dsm = FALSE) {
  configure_catalog_worker()
  finish <- function(chm_raster, maxima = NULL, header = NULL) {
    if (isTRUE(dsm)) write_dsm_chunk(chunk, chm_raster, maxima, header)
    chm_raster
  }
  las <- keep_requested_extra_bytes(
    readLAS(adjust_catalog_boundary_read(chunk, boundary_read_adjustment)), point_order_dimension)
  observation_path <- record_chunk_observation(
    observation_directory,
    "chm",
    nrow(las@data),
    target_chunk_points
  )
  if (is.empty(las)) return(finish(empty_chunk_raster(chunk, chm_resolution)))
  # Sort before normalization, filtering, voxelization, or canopy rasterization.
  las <- restore_original_point_order(las, point_order_dimension, "CHM")
  las <- apply_denoising_mask(las, denoising)
  if (is.empty(las)) return(finish(empty_chunk_raster(chunk, chm_resolution)))
  las <- ground_to_zero(normalize_height(las, terra::rast(dtm_path)), voxel_resolution)
  las <- filter_poi(las, is.finite(Z) & Z >= 0 & Z <= maximum_height)
  if (is.empty(las)) return(finish(empty_chunk_raster(chunk, chm_resolution)))
  dsm_maxima <- if (isTRUE(dsm)) dsm_column_maxima(las@data, voxel_resolution)
  voxels <- voxelize_points(las, res = voxel_resolution)
  if (is.empty(voxels)) return(finish(empty_chunk_raster(chunk, chm_resolution)))
  record_chunk_height(observation_path, voxels@data$Z)
  # Preserve the baseline lidR rasterization semantics.  Passing the numeric
  # resolution lets lidR derive the local extent from the chunk's populated
  # voxels; a fixed core template changes NA/value ownership for points lying
  # on window edges and therefore does not reproduce valid_updated exactly.
  chunk_raster <- rasterize_canopy(
    voxels,
    res = chm_resolution,
    algorithm = p2r()
  )
  header <- voxels@header
  # A worker needs only the raster after this point. Drop the raw points and
  # voxel table before LAScatalog schedules the next spatial window.
  rm(las, voxels)
  gc(verbose = FALSE)
  finish(chunk_raster, dsm_maxima, header)
}

# Reduce only vertical voxel columns: p2r() uses their maximum, so lower voxels
# cannot affect either its result or its XY extent. Never round XY here: the
# exact voxel coordinates are needed for lidR's raster-edge ownership rules.
chm_column_maxima <- function(points) {
  points[, .(Z = max(Z)), by = .(X, Y)]
}

# Keep the baseline raster footprint while reading smaller COPC windows. A
# numeric-resolution p2r() call derives its extent from occupied voxels; doing
# that separately per adaptive read can move a boundary point by one cell.
# Instead retain at most one maximum per 0.2 m XY column of the baseline 60 m
# parent (about 90,000 rows), then rasterize that parent exactly once.
# One half-open read rectangle inside `parent`, as adaptive children are built.
rectangle_cluster <- function(parent, rectangle) {
  child <- lidR:::LAScluster(
    list(x = rectangle[[1]] + rectangle[[3]] / 2, y = rectangle[[2]] + rectangle[[4]] / 2),
    rectangle[[3]], rectangle[[4]], 0, lidR:::LIDRRECTANGLE,
    parent@files, parent@name, sf::st_crs(parent), parent@index
  )
  child@select <- parent@select
  child
}

cluster_key <- function(cluster) {
  bounds <- sf::st_bbox(cluster)
  sprintf("%.6f %.6f %.6f %.6f", bounds$xmin, bounds$ymin, bounds$xmax, bounds$ymax)
}

# A CHM or segment read whose COPC count is far above the chunk target, a dense
# patch in a cloud whose window size came from the mean density, is split into
# quarters.
read_split_points <- function() {
  as.numeric(Sys.getenv("RASTER_SPLIT_POINTS", "5e7"))
}

# Rectangles (xmin, ymin, width, height) that replace each overfull read, keyed
# by that read. Quarters keep the half-open read boundaries, so every point
# still belongs to exactly one read: CHM and DSM keep the maximum per cell, and
# segment summaries are merged by source order, whatever the windows.
plan_dense_reads <- function(reads, boundary_read_adjustment, memory_policy, phase,
                             limit = read_split_points(), max_depth = 4L) {
  if (is.null(memory_policy) || !length(reads)) return(list())
  candidates <- lapply(reads, function(child) {
    bounds <- sf::st_bbox(child)
    list(root = cluster_key(child), parent = child,
         rectangle = c(bounds$xmin, bounds$ymin, bounds$xmax - bounds$xmin,
                       bounds$ymax - bounds$ymin))
  })
  plan <- list()
  depth <- 0L
  while (length(candidates) && depth < max_depth) {
    windows <- lapply(candidates, function(candidate) catalog_task_window(
      adjust_catalog_boundary_read(rectangle_cluster(candidate$parent, candidate$rectangle),
                                   boundary_read_adjustment)))
    counts <- task_memory_estimates(windows, memory_policy, phase)$counts
    following <- list()
    for (i in which(counts > limit)) {
      candidate <- candidates[[i]]
      r <- unname(candidate$rectangle)
      w <- r[[3]] / 2
      h <- r[[4]] / 2
      quarters <- list(c(r[[1]], r[[2]], w, h), c(r[[1]] + w, r[[2]], w, h),
                       c(r[[1]], r[[2]] + h, w, h), c(r[[1]] + w, r[[2]] + h, w, h))
      kept <- Filter(function(existing) !identical(existing, r), plan[[candidate$root]])
      plan[[candidate$root]] <- c(kept, quarters)
      following <- c(following, lapply(quarters, function(quarter)
        list(root = candidate$root, parent = candidate$parent, rectangle = quarter)))
    }
    candidates <- following
    depth <- depth + 1L
  }
  if (length(plan)) {
    message(sprintf("%s: %d dense read(s) split into %d smaller reads above %.0f points",
                    phase, length(plan), sum(lengths(plan)), limit))
  }
  plan
}

# The chunks apply_catalog_chunks() will schedule with these options.
scheduled_chunks <- function(catalog, .options) {
  settings <- lidR:::engine_parse_options(.options)
  lidR:::engine_chunks(catalog, if (settings$check_alignment) settings$raster_alignment else FALSE)
}

# One segment chunk read as several dense-patch pieces, merged as one result.
combine_segment_parts <- function(parts) {
  parts <- Filter(function(part) !is.null(part) && length(part$paths), parts)
  if (!length(parts)) return(list(paths = NULL))
  dimensions <- unique(unlist(lapply(parts, function(part) names(part$dimension_paths))))
  list(
    paths = do.call(c, lapply(parts, `[[`, "paths")),
    dimension_paths = setNames(lapply(dimensions, function(dimension)
      do.call(c, lapply(parts, function(part) part$dimension_paths[[dimension]]))), dimensions)
  )
}

adaptive_chm_read_clusters <- function(chunk, read_size, source_bounds = NULL,
                                       dense_plan = NULL) {
  bounds <- sf::st_bbox(chunk)
  # LAScatalog probes callbacks with a zero-area dummy chunk before scheduling.
  if (bounds$xmax <= bounds$xmin || bounds$ymax <= bounds$ymin) return(list())
  x_starts <- seq(bounds$xmin, bounds$xmax, by = read_size)
  y_starts <- seq(bounds$ymin, bounds$ymax, by = read_size)
  # Header bounds safely exclude empty edge queries without moving the parent
  # grid. Keep a child starting exactly at the source maximum: it owns that edge.
  if (!is.null(source_bounds)) {
    x_starts <- x_starts[x_starts <= source_bounds[["xmax"]] &
                           x_starts + read_size > source_bounds[["xmin"]]]
    y_starts <- y_starts[y_starts <= source_bounds[["ymax"]] &
                           y_starts + read_size > source_bounds[["ymin"]]]
  }
  children <- list()
  # Parent chunks remain in catalog order; children stay inside their parent.
  # Reuse lidR's half-open spatial-filter constructor, with no buffer, so each
  # point belongs to precisely the same baseline parent as before adaptation.
  for (y in y_starts[y_starts < bounds$ymax]) {
    for (x in x_starts[x_starts < bounds$xmax]) {
      width <- min(read_size, bounds$xmax - x)
      height <- min(read_size, bounds$ymax - y)
      child <- lidR:::LAScluster(
        list(x = x + width / 2, y = y + height / 2),
        width, height, 0, lidR:::LIDRRECTANGLE,
        chunk@files, chunk@name, sf::st_crs(chunk), chunk@index
      )
      child@select <- chunk@select
      rectangles <- if (length(dense_plan)) dense_plan[[cluster_key(child)]]
      if (is.null(rectangles)) {
        children[[length(children) + 1L]] <- child
      } else {
        for (rectangle in rectangles) {
          children[[length(children) + 1L]] <- rectangle_cluster(child, rectangle)
        }
      }
    }
  }
  children
}

adaptive_chm_chunk <- function(chunk, dtm_path, voxel_resolution,
                               maximum_height, chm_resolution,
                               point_order_dimension = NULL,
                               observation_directory = NULL,
                               target_chunk_points = 0, read_size,
                               source_bounds = NULL, denoising = NULL,
                               boundary_read_adjustment = NULL, dsm = FALSE,
                               dense_plan = NULL) {
  # Bound both lidR and data.table scratch allocations inside this worker.
  configure_catalog_worker()
  bounds <- sf::st_bbox(chunk)
  if (bounds$xmax <= bounds$xmin || bounds$ymax <= bounds$ymin) return(NULL)
  children <- adaptive_chm_read_clusters(chunk, read_size, source_bounds, dense_plan)
  maxima <- data.table(X = numeric(), Y = numeric(), Z = numeric())
  dsm_maxima <- data.table(X = numeric(), Y = numeric(), Z = numeric())
  header <- NULL
  for (child in children) {
    las <- keep_requested_extra_bytes(
      readLAS(adjust_catalog_boundary_read(child, boundary_read_adjustment)), point_order_dimension)
    observation_path <- record_chunk_observation(observation_directory, "chm", nrow(las@data),
                             target_chunk_points)
    if (is.empty(las)) next
    las <- restore_original_point_order(las, point_order_dimension, "CHM")
    las <- apply_denoising_mask(las, denoising)
    if (is.empty(las)) next
    las <- ground_to_zero(normalize_height(las, terra::rast(dtm_path)), voxel_resolution)
    las <- filter_poi(las, is.finite(Z) & Z >= 0 & Z <= maximum_height)
    if (is.empty(las)) next
    if (isTRUE(dsm)) {
      dsm_maxima <- rbindlist(list(dsm_maxima, dsm_column_maxima(las@data, voxel_resolution)))
      dsm_maxima <- dsm_maxima[, .(Z = max(Z)), by = .(X, Y)]
    }
    voxels <- voxelize_points(las, res = voxel_resolution)
    if (is.empty(voxels)) next
    record_chunk_height(observation_path, voxels@data$Z)
    # Preserve the LAS header to avoid re-quantizing floating voxel centers
    # when reconstructing the compact cloud for the unchanged p2r() call.
    header <- voxels@header
    maxima <- chm_column_maxima(rbindlist(list(
      maxima, chm_column_maxima(voxels@data)
    )))
    # Release every dense table before the next spatial query. Only the tiny
    # two-dimensional maximum table survives across reads, never 3D voxels.
    rm(las, voxels)
    gc(verbose = FALSE)
  }
  if (nrow(maxima) == 0L) {
    chm_raster <- empty_chunk_raster(chunk, chm_resolution)
    if (isTRUE(dsm)) write_dsm_chunk(chunk, chm_raster)
    return(chm_raster)
  }
  compact <- LAS(maxima, header = header, crs = sf::st_crs(chunk), check = FALSE)
  chm_raster <- rasterize_canopy(compact, res = chm_resolution, algorithm = p2r())
  if (isTRUE(dsm)) write_dsm_chunk(chunk, chm_raster, dsm_maxima, header)
  chm_raster
}

# Build the disk-backed CHM. Raster chunks retain the configured baseline size;
# only raw point reads adapt. Completed parent rasters still stream through VRT.
write_global_chm <- function(point_cloud, dtm_path, output_path, work_directory,
                             chunk_size, voxel_resolution, maximum_height,
                             chm_resolution, workers, source_crs = NULL,
                             point_order_dimension = NULL,
                             target_chunk_points = 0,
                             raster_chunk_size = chunk_size, memory_policy = NULL,
                             denoising = NULL, dsm_path = NULL) {
  prepare_raster_work_directory(work_directory, "CHM")
  on.exit(unlink(work_directory, recursive = TRUE, force = TRUE), add = TRUE)
  observation_directory <- file.path(work_directory, "observations")
  dir.create(observation_directory)

  catalog <- readLAScatalog(point_cloud)
  if (!is.null(point_order_dimension)) {
    catalog <- catalog_with_source_boundaries(catalog, point_cloud)
  }
  opt_chunk_size(catalog) <- raster_chunk_size
  opt_chunk_alignment(catalog) <- c(0, 0)
  opt_chunk_buffer(catalog) <- 0
  # Keep the COPC read spatial/selective while adding only the required key.
  opt_select(catalog) <- catalog_selection_for_dimensions(
    point_cloud,
    point_order_dimension,
    "xyz"
  )
  opt_progress(catalog) <- FALSE
  opt_output_files(catalog) <- file.path(work_directory, "chm_{ID}")
  if (is.null(source_crs)) source_crs <- spatial_reference_wkt(catalog)

  # Fixed-size runs keep the original callback unchanged. Adaptive runs share
  # its normalization, voxelization and p2r semantics but bound the raw reads.
  source_bounds <- if (is.null(point_order_dimension)) {
    point_cloud_xy_bounds(point_cloud)
  } else {
    point_cloud_xyz_bounds(point_cloud)[c("xmin", "ymin", "xmax", "ymax")]
  }
  boundary_read_adjustment <- if (is.null(point_order_dimension)) NULL else
    catalog_boundary_read_adjustment(point_cloud)
  dense_plan <- if (is.null(memory_policy)) list() else plan_dense_reads(
    unlist(lapply(scheduled_chunks(catalog, list(automerge = FALSE, drop_null = FALSE,
                                                 raster_alignment = chm_resolution)),
                  adaptive_chm_read_clusters, chunk_size, source_bounds), recursive = FALSE),
    boundary_read_adjustment, memory_policy, "CHM")
  has_dense_reads <- function(chunk) {
    any(vapply(adaptive_chm_read_clusters(chunk, chunk_size, source_bounds),
               function(child) cluster_key(child) %in% names(dense_plan), logical(1)))
  }
  callback <- if (chunk_size < raster_chunk_size || length(dense_plan)) {
    function(chunk, ...) {
      if (chunk_size < raster_chunk_size || has_dense_reads(chunk)) {
        adaptive_chm_chunk(chunk, ..., read_size = chunk_size, source_bounds = source_bounds,
                           dense_plan = dense_plan)
      } else {
        chm_chunk(chunk, ...)
      }
    }
  } else {
    chm_chunk
  }
  # A parent processes these child reads sequentially. Admission must reserve
  # their largest live window, rather than all points in the 60 m parent.
  task_windows_fn <- if (chunk_size < raster_chunk_size || length(dense_plan)) {
    function(chunk) lapply(
      adaptive_chm_read_clusters(chunk, chunk_size, source_bounds, dense_plan),
      function(child) catalog_task_window(adjust_catalog_boundary_read(child, boundary_read_adjustment))
    )
  } else if (!is.null(boundary_read_adjustment)) {
    function(chunk) catalog_task_window(adjust_catalog_boundary_read(chunk, boundary_read_adjustment))
  } else {
    NULL
  }
  results <- apply_catalog_chunks(
      catalog, callback, workers, "CHM",
      memory_policy = memory_policy,
      task_windows_fn = task_windows_fn,
      dtm_path = dtm_path,
      voxel_resolution = voxel_resolution,
      maximum_height = maximum_height,
      chm_resolution = chm_resolution,
      # Each CHM worker sorts only its bounded chunk, never the complete cloud.
      point_order_dimension = point_order_dimension,
      observation_directory = observation_directory,
      target_chunk_points = target_chunk_points,
      denoising = denoising,
      boundary_read_adjustment = boundary_read_adjustment,
      dsm = !is.null(dsm_path),
      .options = list(
        automerge = FALSE,
        drop_null = FALSE,
        raster_alignment = chm_resolution
      )
  )
  result <- build_chunk_virtual_raster(
    results,
    work_directory,
    "CHM",
    source_crs
  )
  covering <- cover_virtual_raster(
    result,
    point_cloud,
    chm_resolution,
    work_directory,
    "CHM",
    voxel_resolution = if (is.null(point_order_dimension)) NULL else voxel_resolution
  )
  stream_virtual_raster(covering, output_path, "CHM", source_crs)
  if (!is.null(dsm_path)) {
    # The DSM mosaic repeats the CHM's VRT, covering-extent and streaming steps
    # on the DSM chunks, so both rasters share one grid and nodata policy.
    dsm_results <- lapply(results, function(path) if (!is.null(path)) dsm_chunk_path(path))
    dsm_covering <- cover_virtual_raster(
      build_chunk_virtual_raster(dsm_results, work_directory, "DSM", source_crs),
      point_cloud,
      chm_resolution,
      work_directory,
      "DSM",
      voxel_resolution = if (is.null(point_order_dimension)) NULL else voxel_resolution
    )
    stream_virtual_raster(dsm_covering, dsm_path, "DSM", source_crs)
  }
  summarize_chunk_observations(observation_directory)
}
