# The disk-backed DTM pipeline: catalog chunks, ground classification, chunk
# rasters, their CRS and the global DTM. Sourced by analysis.R; uses its
# constants and helpers.

# Disk-backed DTM and CHM pipeline -----------------------------------------

# Determine whether an XY sequence contains a triangle without building the
# full unique-coordinate table used by the former guard.  Dense DTM windows can
# contain hundreds of millions of rows, so even a temporary duplicate XY table
# can exceed the process budget before PTD or TIN begins.  Fixed-size blocks
# bound every temporary vector while preserving the same yes/no geometry guard.
xy_has_nondegenerate_points <- function(x, y, classes = NULL,
                                        eligible_class = NULL,
                                        block_size = 1000000L) {
  if (length(x) != length(y)) stop("XY vectors must have equal length")
  if (!is.null(classes) && length(classes) != length(x)) {
    stop("classification and XY vectors must have equal length")
  }
  if (!is.null(classes) && length(eligible_class) != 1L) {
    stop("an eligible class is required for a classified XY scan")
  }
  block_size <- suppressWarnings(as.integer(block_size))
  if (length(block_size) != 1L || is.na(block_size) || block_size < 1L) {
    stop("XY geometry scan block size must be a positive integer")
  }
  point_count <- length(x)
  if (point_count < 3L) return(FALSE)

  first_point <- NULL
  second_point <- NULL
  for (start in seq.int(1L, point_count, by = block_size)) {
    finish <- min(point_count, as.double(start) + block_size - 1)
    rows <- seq.int(start, finish)
    block_x <- x[rows]
    block_y <- y[rows]
    eligible <- is.finite(block_x) & is.finite(block_y)
    if (!is.null(classes)) {
      eligible <- eligible & !is.na(classes[rows]) &
        classes[rows] == eligible_class
    }
    candidates <- which(eligible)
    if (length(candidates) == 0L) next

    if (is.null(first_point)) {
      first_row <- candidates[[1L]]
      first_point <- c(block_x[[first_row]], block_y[[first_row]])
    }
    if (is.null(second_point)) {
      different <- candidates[
        block_x[candidates] != first_point[[1L]] |
          block_y[candidates] != first_point[[2L]]
      ]
      if (length(different) > 0L) {
        second_row <- different[[1L]]
        second_point <- c(block_x[[second_row]], block_y[[second_row]])
      }
    }
    if (is.null(second_point)) next

    delta_x <- second_point[[1L]] - first_point[[1L]]
    delta_y <- second_point[[2L]] - first_point[[2L]]
    candidate_dx <- block_x[candidates] - first_point[[1L]]
    candidate_dy <- block_y[candidates] - first_point[[2L]]
    cross_product <- delta_x * candidate_dy - delta_y * candidate_dx
    magnitude <- pmax(
      1,
      abs(delta_x * candidate_dy) + abs(delta_y * candidate_dx)
    )
    if (any(abs(cross_product) > 1e-12 * magnitude)) return(TRUE)
  }
  FALSE
}

# The source ordinal is required to establish the original row order, but
# PTD and TIN do not read it afterward.  Releasing this UInt64 column before
# those dense DTM operations shortens its lifetime without changing DTM chunks,
# buffers, classification, triangulation, coordinates, or row order.
release_point_order_dimension <- function(las, point_order_dimension = NULL) {
  if (is.null(point_order_dimension) ||
      !point_order_dimension %in% names(las@data)) return(las)
  las@data[, (point_order_dimension) := NULL]
  gc(verbose = FALSE)
  las
}

# Process one DTM catalog chunk.  The optional order field is selected by the
# catalog caller and restored immediately after reading.
# lidR decodes a COPC window on one core. A task given more threads (dense
# windows run one at a time) copies its window to an uncompressed LAS with that
# many processes; lidR then reads the same points, with the same selection and
# filter, from the copy.
read_catalog_chunk <- function(chunk) {
  threads <- suppressWarnings(as.integer(Sys.getenv("RASTER_TASK_THREADS", "1")))
  if (is.na(threads) || threads < 2L || length(chunk@files) != 1L || !nzchar(chunk@save) ||
      !is_copc_point_cloud(chunk@files)) {
    return(readLAS(chunk))
  }
  copy <- paste0(chunk@save, "_points.las")
  on.exit(unlink(copy), add = TRUE)
  status <- system2("python3", c(shQuote(file.path(script_directory, "read_window.py")),
    "--input", shQuote(chunk@files), "--output", shQuote(copy),
    "--bounds", format(catalog_task_window(chunk), digits = 17, scientific = FALSE),
    "--threads", threads))
  if (!identical(status, 0L)) stop("parallel COPC window read failed")
  chunk@files <- copy
  readLAS(chunk)
}

dtm_chunk <- function(chunk, dtm_resolution, ptd_resolution,
                      point_order_dimension = NULL, denoising = NULL) {
  configure_catalog_worker()
  las <- keep_requested_extra_bytes(read_catalog_chunk(chunk), point_order_dimension)
  if (is.empty(las)) return(empty_chunk_raster(chunk, dtm_resolution))
  # Treat a missing order key as a hard error before classification or TIN work.
  assert_loaded_dimensions(
    las,
    c("X", "Y", "Z", "buffer", point_order_dimension),
    "DTM"
  )
  # Match the record sequence of the original LAZ-based workflow.
  las <- restore_original_point_order(las, point_order_dimension, "DTM")
  # Rejected points often define the window's edge; keep the unfiltered
  # window's grid so the DTM still covers the full extent (NA left to TIN
  # extrapolation, as without the filter).
  layout <- if (!is.null(denoising)) lidR:::raster_layout(las, dtm_resolution, format = "terra")
  las <- apply_denoising_mask(las, denoising)
  # The order is now encoded by the rows themselves; do not retain a redundant
  # eight-byte ordinal per point through the memory-intensive PTD/TIN stages.
  las <- release_point_order_dimension(las, point_order_dimension)
  if (!las_has_nondegenerate_xy(las)) {
    message(
      "DTM chunk has fewer than three non-collinear input points; ",
      "leaving it empty for buffered neighbouring chunks"
    )
    return(empty_chunk_raster(chunk, dtm_resolution))
  }
  # The original DTM does not use return-number semantics. Being explicit avoids
  # lidR changing behaviour when ReturnNumber happens to be present.
  las <- classify_ground(
    las,
    ptd(res = ptd_resolution),
    last_returns = FALSE
  )
  surface <- classified_ground_surface(las, if (is.null(layout)) dtm_resolution else layout)
  if (is.null(surface)) empty_chunk_raster(chunk, dtm_resolution) else surface
}

# Build the same floating-point Delaunay algorithm that lidR 4.3.2 uses when
# its optimized integer path rejects a LAS scale, offset, or coordinate. Calling
# the legacy primitives explicitly makes the retry independent of lidR's sampled
# fast-path precheck and preserves every source coordinate and header value.
legacy_tin <- function(..., extrapolate = knnidw(3, 1, 50)) {
  # Resolve the extrapolation algorithm once, matching lidR::tin() semantics.
  extrapolate <- lazyeval::uq(extrapolate)
  # Define the DTM plugin callback that receives one LAS chunk and raster grid.
  algorithm <- function(las, where) {
    # A numeric matrix forces tDelaunay() onto its legacy floating-point branch.
    points <- as.matrix(las@data)
    # tInterpolate() expects candidate raster locations as a numeric XY matrix.
    query <- as.matrix(where)
    # lidR 4.3.2 drops matrix dimensions inside tInterpolate() for one query.
    single_query <- nrow(query) == 1L
    # Duplicate that one location solely to retain matrix shape in the legacy API.
    interpolation_query <- if (single_query) rbind(query, query) else query
    # Build the legacy triangulation without changing LAS coordinates or headers.
    triangulation <- lidR:::tDelaunay(points, trim = 0)
    # Interpolate grid elevations with the configured lidR thread allowance.
    elevations <- lidR:::tInterpolate(
      triangulation,
      points,
      interpolation_query,
      threads = lidR::get_lidr_threads()
    )
    # Discard the duplicate result introduced only for the one-cell API guard.
    if (single_query) elevations <- elevations[[1L]]
    # Find grid cells outside the convex hull, which Delaunay leaves as NA.
    outside <- is.na(elevations)
    # Preserve lidR::tin()'s configured nearest-neighbour extrapolation behavior.
    if (any(outside)) {
      # Tell the extrapolation plugin that it runs in spatial-interpolation mode.
      lidR.context <- "spatial_interpolation"
      # Pass only unresolved XY cells to avoid repeating successful interpolation.
      outside_coordinates <- data.frame(
        X = where$X[outside],
        Y = where$Y[outside]
      )
      # Fill outside-hull cells with the same default used by lidR::tin().
      elevations[outside] <- extrapolate(las, outside_coordinates)
    }
    # Return one elevation per requested raster location to the DTM plugin.
    elevations
  }
  # Mark the callback as a lidR DTM algorithm with threaded interpolation support.
  lidR:::plugin_dtm(algorithm, omp = TRUE)
}

las_has_nondegenerate_xy <- function(las) {
  xy_has_nondegenerate_points(las@data[["X"]], las@data[["Y"]])
}

chunk_has_valid_extent <- function(chunk) {
  chunk_bounds <- st_bbox(chunk)
  all(is.finite(chunk_bounds)) &&
    chunk_bounds[["xmax"]] > chunk_bounds[["xmin"]] &&
    chunk_bounds[["ymax"]] > chunk_bounds[["ymin"]]
}

ensure_nonempty_grid_bounds <- function(bounds, resolution) {
  if (bounds[["xmax"]] <= bounds[["xmin"]]) {
    bounds[["xmax"]] <- bounds[["xmin"]] + resolution
  }
  if (bounds[["ymax"]] <= bounds[["ymin"]]) {
    bounds[["ymax"]] <- bounds[["ymin"]] + resolution
  }
  bounds
}

empty_chunk_raster <- function(chunk, resolution) {
  if (!chunk_has_valid_extent(chunk)) return(NULL)
  chunk_bounds <- st_bbox(chunk)
  bounds <- ensure_nonempty_grid_bounds(
    covering_grid_bounds(chunk_bounds, resolution),
    resolution
  )
  raster <- terra::rast(
    xmin = bounds[["xmin"]],
    xmax = bounds[["xmax"]],
    ymin = bounds[["ymin"]],
    ymax = bounds[["ymax"]],
    resolution = resolution
  )
  terra::values(raster) <- NA_real_
  raster
}

# `dtm_resolution` is a cell size or, with the noise filter, the unfiltered
# window's raster template.
classified_ground_surface <- function(las, dtm_resolution) {
  classifications <- las[["Classification"]]
  if (is.null(classifications)) {
    message(
      "DTM chunk contains no PTD ground points; leaving it empty for ",
      "buffered neighbouring chunks"
    )
    return(NULL)
  }
  if (!xy_has_nondegenerate_points(
    las@data[["X"]],
    las@data[["Y"]],
    classifications,
    2L
  )) {
    message(
      "DTM chunk has fewer than three non-collinear PTD ground points; ",
      "leaving it empty for buffered neighbouring chunks"
    )
    return(NULL)
  }
  tryCatch(
    rasterize_terrain(las, res = dtm_resolution, algorithm = tin()),
    error = function(error) {
      # Retry only the known lidR integer-conversion failure.  Unrelated TIN
      # errors must retain their existing behavior and diagnostic path.
      if (grepl(
        "xy coordinates were not converted to integer",
        conditionMessage(error),
        fixed = TRUE
      )) {
        message(
          paste(
            "Retrying DTM TIN with lidR's legacy floating-point implementation;",
            "coordinates and header are unchanged"
          )
        )
        # Repeat the same DTM resolution with the original LAS and slow TIN path.
        return(rasterize_terrain(
          las,
          res = dtm_resolution,
          algorithm = legacy_tin()
        ))
      }
      if (grepl(
        paste0(
          "cannot triangulate less than 3 points|No ground points found|",
          "Invalid grid dimensions"
        ),
        conditionMessage(error),
        ignore.case = TRUE
      )) {
        message(
          "DTM chunk has insufficient ground geometry for TIN; ",
          "leaving it empty for buffered neighbouring chunks"
        )
        return(NULL)
      }
      stop(error)
    }
  )
}

assert_catalog_completed <- function(results, phase) {
  if (length(results) == 0L || any(vapply(results, function(result) {
    is.null(result) || inherits(result, "try-error")
  }, logical(1)))) {
    stop(sprintf(
      "%s LAScatalog processing failed or returned incomplete chunks",
      phase
    ))
  }
  invisible(results)
}

with_catalog_workers <- function(workers, phase, expression) {
  workers <- max(1L, as.integer(workers))
  previous_lidr_threads <- lidR::get_lidr_threads()
  on.exit(
    lidR::set_lidr_threads(previous_lidr_threads),
    add = TRUE
  )
  message(sprintf(
    "Processing %s LAScatalog chunks with %d worker(s)",
    phase,
    workers
  ))
  if (workers == 1L) return(force(expression))

  previous_plan <- future::plan()
  on.exit(future::plan(previous_plan), add = TRUE)
  future::plan(future::multisession, workers = workers)
  force(expression)
}

prepare_raster_work_directory <- function(path, label) {
  if (dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
  if (!dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop(sprintf("could not create temporary %s chunk directory", label))
  }
  invisible(path)
}

virtual_raster_source <- function(result, label) {
  if (!inherits(result, "SpatRaster") || terra::inMemory(result)) {
    stop(sprintf(
      "global %s processing did not produce a disk-backed raster mosaic",
      label
    ))
  }
  sources <- terra::sources(result)
  if (length(sources) != 1L ||
      !grepl("\\.vrt$", sources[[1]], ignore.case = TRUE)) {
    stop(sprintf(
      "global %s processing did not produce one virtual raster mosaic",
      label
    ))
  }
  sources[[1]]
}

complete_chunk_raster_paths <- function(results, label) {
  assert_catalog_completed(results, label)
  paths <- unname(unlist(results, use.names = FALSE))
  if (length(paths) == 0L || any(!file.exists(paths))) {
    stop(sprintf(
      "%s LAScatalog processing did not write every chunk raster",
      label
    ))
  }

  paths
}

populated_or_all_empty_chunk_raster_paths <- function(paths, label) {
  has_data <- vapply(paths, function(path) {
    count <- suppressWarnings(terra::global(terra::rast(path), "notNA"))
    is.finite(count[[1L, 1L]]) && count[[1L, 1L]] > 0
  }, logical(1))
  if (any(has_data)) {
    return(paths[has_data])
  }

  message(sprintf(
    "global %s processing produced only empty chunks; publishing an all-NoData raster",
    label
  ))
  paths
}

chunk_reference_crs <- function(paths, label) {
  candidates <- vapply(paths, function(path) {
    terra::crs(terra::rast(path))
  }, character(1))
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  if (length(candidates) > 1L) {
    stop(sprintf("%s chunk rasters contain conflicting CRS metadata", label))
  }
  if (length(candidates) == 1L) candidates[[1L]] else NA_character_
}

spatial_reference_wkt <- function(value) {
  reference <- sf::st_crs(value)
  if (is.na(reference) || is.null(reference$wkt) ||
      is.na(reference$wkt) || !nzchar(reference$wkt)) {
    return(NA_character_)
  }
  reference$wkt
}

usable_crs_wkt <- function(value) {
  length(value) == 1L && !is.null(value) && !is.na(value) && nzchar(value)
}

# A geocentric CRS (EPSG:4978) on a cloud near 0,0,0 is Standardization 1.2.1's placeholder
# for uploads in local coordinates (3DT-3026). Real Earth-centred coordinates lie about
# 6,370 km from the centre; nothing that close to the origin can be on the Earth's surface.
LOCAL_PLACEHOLDER_MAX_CENTRE_DISTANCE_M <- 1e6
is_geocentric_crs <- function(wkt) {
  usable_crs_wkt(wkt) &&
    grepl("^\\s*(GEOCCS|GEODCRS|GEODETICCRS)\\b", wkt) &&
    grepl("Cartesian", wkt, fixed = TRUE)
}

# Whether two CRS share the same horizontal reference (raster_crs.py, GDAL).
same_horizontal_crs <- function(a, b) {
  result <- system2("python3", c(
    shQuote(file.path(script_directory, "raster_crs.py")), "--same",
    shQuote(a), shQuote(b)), stdout = TRUE)
  status <- attr(result, "status")
  if ((!is.null(status) && status != 0L) || !length(result)) {
    stop("could not compare the point-cloud CRS with the supplied CRS")
  }
  identical(trimws(result[[length(result)]]), "true")
}

# The output CRS: the point cloud's header, else its original companion, else the
# CRS the caller supplied from the dataset record (3DT-3029). A placeholder header
# means local coordinates: the products get no CRS and nothing is supplied to them.
resolve_analysis_crs <- function(point_cloud_crs, original_companion_crs = NA_character_,
                                 supplied_crs = NA_character_,
                                 same_crs = same_horizontal_crs,
                                 centre_distance_m = 0) {
  header <- if (usable_crs_wkt(point_cloud_crs)) {
    list(wkt = point_cloud_crs, source = "point_cloud")
  } else if (usable_crs_wkt(original_companion_crs)) {
    list(wkt = original_companion_crs, source = "original_companion")
  } else NULL
  if (!is.null(header) && is_geocentric_crs(header$wkt)) {
    if (!is.finite(centre_distance_m) || centre_distance_m > LOCAL_PLACEHOLDER_MAX_CENTRE_DISTANCE_M) {
      stop(paste(
        "the point cloud is in Earth-centred (ECEF, geocentric) coordinates;",
        "reproject it to a projected CRS in metres before the raster tool"
      ))
    }
    if (usable_crs_wkt(supplied_crs)) {
      message("Point cloud is in local coordinates (geocentric placeholder CRS); ignoring --crs")
    }
    return(list(wkt = NA_character_, source = "local"))
  }
  if (!is.null(header)) {
    if (usable_crs_wkt(supplied_crs) && !same_crs(header$wkt, supplied_crs)) {
      stop(paste(
        "the supplied --crs differs from the point cloud's CRS; the header CRS",
        "is kept only when they agree"
      ))
    }
    return(header)
  }
  if (usable_crs_wkt(supplied_crs)) {
    return(list(wkt = supplied_crs, source = "supplied"))
  }
  list(wkt = NA_character_, source = "unavailable")
}

# ``--crs`` as WKT: an EPSG code (EPSG:25832) or any CRS PROJ reads, projected in metres.
supplied_crs_wkt <- function(value) {
  if (is.null(value) || !nzchar(value)) return(NA_character_)
  reference <- sf::st_crs(value)
  if (is.na(reference)) stop(sprintf("--crs %s is not a CRS PROJ can read", value))
  if (isTRUE(reference$IsGeographic) || is_geocentric_crs(reference$wkt)) {
    stop(sprintf("--crs %s must be a projected CRS in metres", value))
  }
  units <- reference$units_gdal
  if (!is.null(units) && !is.na(units) && !units %in% c("metre", "meter", "Meter")) {
    stop(sprintf("--crs %s must be a projected CRS in metres", value))
  }
  reference$wkt
}

build_chunk_virtual_raster <- function(results, work_directory, label,
                                       fallback_crs = NA_character_) {
  all_paths <- complete_chunk_raster_paths(results, label)
  reference_crs <- chunk_reference_crs(all_paths, label)
  if ((is.na(reference_crs) || !nzchar(reference_crs)) &&
      !is.na(fallback_crs) && nzchar(fallback_crs)) {
    reference_crs <- fallback_crs
  }
  paths <- populated_or_all_empty_chunk_raster_paths(all_paths, label)
  vrt_path <- file.path(
    work_directory,
    sprintf("%s_chunks.vrt", tolower(label))
  )
  input_list_path <- file.path(
    work_directory,
    sprintf("%s_chunks.txt", tolower(label))
  )
  unlink(vrt_path, force = TRUE)
  writeLines(paths, input_list_path, useBytes = TRUE)
  arguments <- c("-overwrite")
  if (!is.na(reference_crs) && nzchar(reference_crs)) {
    arguments <- c(arguments, "-a_srs", shQuote(reference_crs))
  }
  arguments <- c(
    arguments,
    "-input_file_list", shQuote(input_list_path),
    shQuote(vrt_path)
  )
  status <- system2("gdalbuildvrt", arguments)
  if (!identical(status, 0L) || !file.exists(vrt_path)) {
    stop(sprintf("could not build the global %s chunk VRT", label))
  }
  result <- terra::rast(vrt_path)
  if (!is.na(reference_crs) && nzchar(reference_crs)) {
    actual_crs <- terra::crs(result)
    if (is.na(actual_crs) || !nzchar(actual_crs)) {
      stop(sprintf(
        "global %s chunk VRT did not persist its reference CRS",
        label
      ))
    }
  }
  result
}

gdal_raster_has_crs <- function(path) {
  info <- system2(
    "gdalinfo",
    c("--config", "GDAL_PAM_ENABLED", "NO", "-json", shQuote(path)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(info, "status")
  if (!is.null(status) && status != 0L) {
    stop(sprintf("could not inspect the global raster CRS for %s", path))
  }
  any(grepl('"coordinateSystem"', info, fixed = TRUE))
}

align_raster_file_crs <- function(path, source_crs, label) {
  if (is.null(source_crs)) return(invisible(path))

  source_has_crs <- !is.na(source_crs) && nzchar(source_crs)
  definition <- if (source_has_crs) source_crs else "None"
  if (source_has_crs) {
    horizontal <- system2("python3", c(
      shQuote(file.path(script_directory, "raster_crs.py")), shQuote(source_crs)),
      stdout = TRUE)
    status <- attr(horizontal, "status")
    if ((!is.null(status) && status != 0L) || !length(horizontal)) {
      stop(sprintf("could not resolve the horizontal %s raster CRS", label))
    }
    definition <- paste(horizontal, collapse = "\n")
  }
  status <- system2(
    "gdal_edit.py",
    c("-a_srs", shQuote(definition), shQuote(path))
  )
  if (!identical(status, 0L)) {
    stop(sprintf(
      "could not align the global %s raster CRS with the point-cloud source",
      label
    ))
  }

  output_has_crs <- gdal_raster_has_crs(path)
  if (!identical(output_has_crs, source_has_crs)) {
    stop(sprintf(
      "global %s raster CRS does not match the point-cloud source state",
      label
    ))
  }
  invisible(path)
}

stream_virtual_raster <- function(result, output_path, label,
                                  source_crs = NULL) {
  source_path <- virtual_raster_source(result, label)
  partial_path <- paste0(output_path, ".partial.tif")
  unlink(partial_path, force = TRUE)
  on.exit(unlink(partial_path, force = TRUE), add = TRUE)
  status <- system2(
    "gdal_translate",
    c(
      "-of", "GTiff",
      "-co", "TILED=YES",
      "-co", "COMPRESS=DEFLATE",
      "-co", "BIGTIFF=IF_SAFER",
      shQuote(source_path),
      shQuote(partial_path)
    )
  )
  if (!identical(status, 0L) || !file.exists(partial_path)) {
    stop(sprintf(
      "could not stream the global %s virtual mosaic to GeoTIFF",
      label
    ))
  }
  if (!file.rename(partial_path, output_path)) {
    stop(sprintf("could not atomically publish the global %s GeoTIFF", label))
  }
  align_raster_file_crs(output_path, source_crs, label)
  convert_raster_format(output_path, label)
}

# RASTER_FORMAT=COG rewrites a finished GeoTIFF as a Cloud
# Optimized GeoTIFF; the CRS is aligned first because editing a COG in place
# would break its layout.
convert_raster_format <- function(path, label, resampling = "AVERAGE") {
  if (!identical(Sys.getenv("RASTER_FORMAT", "GTiff"), "COG")) {
    return(invisible(path))
  }
  partial_path <- paste0(path, ".partial.tif")
  unlink(partial_path, force = TRUE)
  on.exit(unlink(partial_path, force = TRUE), add = TRUE)
  status <- system2(
    "gdal_translate",
    c(
      "-of", "COG",
      "-co", "COMPRESS=DEFLATE",
      "-co", "BIGTIFF=IF_SAFER",
      "-co", paste0("RESAMPLING=", resampling),
      shQuote(path),
      shQuote(partial_path)
    )
  )
  if (!identical(status, 0L) || !file.exists(partial_path) ||
      !file.rename(partial_path, path)) {
    stop(sprintf("could not write the global %s raster as a COG", label))
  }
  invisible(path)
}

align_grid_coordinate <- function(value, resolution, direction) {
  scaled <- value / resolution
  nearest <- round(scaled)
  if (abs(scaled - nearest) < 1e-9) scaled <- nearest
  if (identical(direction, "lower")) {
    floor(scaled) * resolution
  } else {
    ceiling(scaled) * resolution
  }
}

covering_grid_bounds <- function(bounds, resolution) {
  c(
    xmin = align_grid_coordinate(bounds[["xmin"]], resolution, "lower"),
    ymin = align_grid_coordinate(bounds[["ymin"]], resolution, "lower"),
    xmax = align_grid_coordinate(bounds[["xmax"]], resolution, "upper"),
    ymax = align_grid_coordinate(bounds[["ymax"]], resolution, "upper")
  )
}

lidr_catalog_grid_bounds <- function(point_cloud, resolution) {
  required_version <- "4.3.2"
  actual_version <- as.character(utils::packageVersion("lidR"))
  if (!identical(actual_version, required_version)) {
    stop(sprintf(
      paste0(
        "the published DTM layout requires pinned lidR %s; ",
        "found %s"
      ),
      required_version,
      actual_version
    ))
  }

  # The original rasterize_terrain() call delegates numeric-resolution
  # geometry to this pinned lidR layout. In particular, it retains an extra
  # cell when a LAS bound lies exactly on the grid alignment. Reusing the same
  # implementation avoids approximating that half-cell rounding convention.
  layout <- lidR:::raster_layout(
    readLAScatalog(point_cloud),
    resolution,
    start = c(0, 0),
    buffer = 0,
    format = "template"
  )
  bounds <- c(
    xmin = as.numeric(layout$xmin),
    ymin = as.numeric(layout$ymin),
    xmax = as.numeric(layout$xmax),
    ymax = as.numeric(layout$ymax)
  )
  if (any(!is.finite(bounds)) ||
      bounds[["xmax"]] <= bounds[["xmin"]] ||
      bounds[["ymax"]] <= bounds[["ymin"]]) {
    stop("pinned lidR returned an invalid DTM layout")
  }
  bounds
}

raster_extent_covers_point_bounds <- function(raster_extent, point_bounds,
                                               resolution) {
  extent_values <- c(
    xmin = as.numeric(raster_extent$xmin),
    ymin = as.numeric(raster_extent$ymin),
    xmax = as.numeric(raster_extent$xmax),
    ymax = as.numeric(raster_extent$ymax)
  )
  coordinate_scale <- max(abs(c(extent_values, point_bounds)), 1)
  tolerance <- max(
    resolution * 1e-9,
    .Machine$double.eps * coordinate_scale * 16
  )
  extent_values[["xmin"]] <= point_bounds[["xmin"]] + tolerance &&
    extent_values[["xmax"]] >= point_bounds[["xmax"]] - tolerance &&
    extent_values[["ymin"]] <= point_bounds[["ymin"]] + tolerance &&
    extent_values[["ymax"]] >= point_bounds[["ymax"]] - tolerance
}

cover_virtual_raster <- function(result, point_cloud, resolution,
                                 work_directory, label,
                                 extent_mode = c("point_bounds", "lidr"),
                                 voxel_resolution = NULL) {
  extent_mode <- match.arg(extent_mode)
  source_path <- virtual_raster_source(result, label)
  point_bounds <- if (is.null(voxel_resolution)) {
    point_cloud_xy_bounds(point_cloud)
  } else {
    point_cloud_xyz_bounds(point_cloud)[c("xmin", "ymin", "xmax", "ymax")]
  }
  target_bounds <- if (identical(extent_mode, "lidr")) {
    lidr_catalog_grid_bounds(point_cloud, resolution)
  } else {
    covering_grid_bounds(point_bounds, resolution)
  }
  if (!is.null(voxel_resolution)) {
    # Reuse pinned voxelize_points' exact XY center arithmetic. Extend only to
    # cells needed for source extrema, never to empty catalog-chunk rasters.
    # Strict upper coverage retains a voxel center exactly on a raster edge.
    centers <- lidR:::f_grid(point_bounds, voxel_resolution, 0)
    target_bounds[1:2] <- pmin(target_bounds[1:2], floor(centers[1:2] / resolution) * resolution)
    target_bounds[3:4] <- pmax(target_bounds[3:4], (floor(centers[3:4] / resolution) + 1) * resolution)
  }
  if (identical(extent_mode, "lidr") &&
      !raster_extent_covers_point_bounds(
        terra::ext(result),
        target_bounds,
        resolution
      )) {
    stop(sprintf(
      paste0(
        "global %s mosaic does not contain the complete pinned-lidR ",
        "publication extent"
      ),
      label
    ))
  }
  covering_path <- file.path(
    work_directory,
    sprintf("%s_covering.vrt", tolower(label))
  )
  number <- function(value) {
    format(value, scientific = FALSE, trim = TRUE, digits = 15)
  }
  status <- system2(
    "gdalbuildvrt",
    c(
      "-overwrite",
      "-resolution", "user",
      "-tr", number(resolution), number(resolution),
      "-te",
      number(target_bounds[["xmin"]]),
      number(target_bounds[["ymin"]]),
      number(target_bounds[["xmax"]]),
      number(target_bounds[["ymax"]]),
      "-vrtnodata", "-999999",
      shQuote(covering_path),
      shQuote(source_path)
    )
  )
  if (!identical(status, 0L) || !file.exists(covering_path)) {
    stop(sprintf(
      "could not create the point-cloud-covering global %s VRT",
      label
    ))
  }
  covering <- terra::rast(covering_path)
  if (terra::inMemory(covering)) {
    stop(sprintf("global %s covering VRT materialized in memory", label))
  }
  covering_extent <- terra::ext(covering)
  if (!raster_extent_covers_point_bounds(
    covering_extent,
    point_bounds,
    resolution
  )) {
    stop(sprintf(
      "global %s covering VRT does not contain the point-cloud XY extent",
      label
    ))
  }
  covering
}

# Build the disk-backed global DTM while carrying the COPC order key through
# lidR's selective catalog reads.
write_global_dtm <- function(point_cloud, output_path, work_directory,
                             chunk_size, dtm_buffer, dtm_resolution,
                             ptd_resolution, workers, source_crs = NULL,
                             point_order_dimension = NULL, memory_policy = NULL,
                             denoising = NULL, extent_cloud = point_cloud) {
  prepare_raster_work_directory(work_directory, "DTM")
  on.exit(unlink(work_directory, recursive = TRUE, force = TRUE), add = TRUE)

  catalog <- readLAScatalog(point_cloud)
  # A thinned cloud keeps the source's chunk layout and raster extent.
  source_extent <- readLAScatalog(extent_cloud)@data
  for (name in c("Min.X", "Min.Y", "Max.X", "Max.Y")) catalog@data[[name]] <- source_extent[[name]]
  opt_chunk_size(catalog) <- chunk_size
  opt_chunk_buffer(catalog) <- dtm_buffer
  # Plain LAZ requests xyz; canonical COPC adds the OriginalPointIndex ordinal.
  opt_select(catalog) <- catalog_selection_for_dimensions(
    point_cloud,
    point_order_dimension,
    DTM_CATALOG_SELECTION
  )
  opt_progress(catalog) <- FALSE
  opt_output_files(catalog) <- file.path(work_directory, "dtm_{ID}")
  catalog@output_options$drivers$SpatRaster$param$datatype <- "FLT8S"
  if (is.null(source_crs)) source_crs <- spatial_reference_wkt(catalog)

  results <- apply_catalog_chunks(
      catalog, dtm_chunk, workers, "DTM",
      memory_policy = memory_policy,
      dtm_resolution = dtm_resolution,
      ptd_resolution = ptd_resolution,
      # Every worker restores source record order before ground/TIN operations.
      point_order_dimension = point_order_dimension,
      denoising = denoising,
      .options = list(
        automerge = FALSE,
        drop_null = FALSE,
        raster_alignment = dtm_resolution
      )
  )
  chunk_paths <- complete_chunk_raster_paths(results, "DTM")
  mosaic <- build_dtm_overlap_mosaic(
    chunk_paths,
    work_directory,
    fallback_crs = source_crs
  )
  covering <- cover_virtual_raster(
    mosaic,
    extent_cloud,
    dtm_resolution,
    work_directory,
    "DTM",
    extent_mode = "lidr"
  )
  stream_virtual_raster(covering, output_path, "DTM", source_crs)
}
