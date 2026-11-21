# The complete set of incoming outputs and their checks before promotion:
# expected files, the performance record's resources, the surface rasters.
# Sourced by raster.R.

expected_artifacts <- function(directory, dataset_id, dimensions, denoising = FALSE,
                               debug = TRUE) {
  prefix <- file.path(directory, as.character(dataset_id))
  common <- c(
    paste0(prefix, "_dtm.tif"),
    paste0(prefix, "_chm.tif"),
    paste0(prefix, "_dsm.tif"),
    if (debug) paste0(prefix, "_performance.csv"),
    paste0(prefix, "_aoi_conversion.json"),
    paste0(prefix, "_run.json")
  )
  # Structural metrics come per instance dimension; without one there are none (3DT-3007).
  suffixes <- if (length(dimensions)) paste0("_", dimensions) else character()
  dimension_files <- unlist(lapply(suffixes, function(suffix) c(
    paste0(prefix, suffix, "_results.csv"),
    paste0(prefix, suffix, "_segment_diagnostics.csv"),
    paste0(prefix, suffix, "_tiles.geojson")
  )))
  c(common, if (length(dimensions)) paste0(prefix, "_tiles.png"), dimension_files,
    if (denoising) paste0(prefix, c("_noise_points.csv", "_denoising.json")))
}

validate_performance_resources <- function(performance, memory_budget_gib,
                                           source_point_count,
                                           memory_limit_gib = memory_budget_gib) {
  # A debug report; it never fails the run (Kilian, 21 Nov).
  if (nrow(performance) != 1L ||
      !isTRUE(performance$point_count[[1]] == source_point_count) ||
      !is.finite(performance$peak_rss_mib[[1]]) ||
      !is.finite(performance$cgroup_memory_peak_mib[[1]]) ||
      !is.finite(performance$temporary_disk_peak_mib[[1]])) {
    message("Warning: the performance report is incomplete or does not match the source")
  }
  invisible(TRUE)
}

# Galaxy enforces the memory grant itself, so memory peaks only warn, in every run.
# The kernel's cgroup peak also counts reclaimable disk cache (dataset 97: 77.4 GB
# peak, 0.7 GB main-process RSS, 40.9 GB grant).
warn_memory_peaks <- function(performance, memory_budget_gib,
                              memory_limit_gib = memory_budget_gib) {
  if (isTRUE(performance$peak_rss_mib[[1]] >= memory_budget_gib * 1024)) {
    message(sprintf("Warning: peak RSS %.0f MiB reached the %.0f MiB internal budget",
      performance$peak_rss_mib[[1]], memory_budget_gib * 1024))
  }
  if (isTRUE(performance$cgroup_memory_peak_mib[[1]] > memory_limit_gib * 1024)) {
    message(sprintf(paste(
      "Warning: cgroup memory peak %.0f MiB exceeds the %.0f MiB limit (it includes",
      "reclaimable disk cache; main-process peak RSS %.0f MiB)"),
      performance$cgroup_memory_peak_mib[[1]], memory_limit_gib * 1024,
      performance$peak_rss_mib[[1]]))
  }
  invisible(TRUE)
}

# The DSM must share the CHM's grid, CRS and empty cells.
validate_surface_rasters <- function(chm, dsm) {
  if (terra::inMemory(dsm) ||
      !isTRUE(terra::compareGeom(chm, dsm, crs = TRUE, stopOnError = FALSE)) ||
      !identical(terra::crs(chm), terra::crs(dsm))) {
    stop("incoming DSM grid differs from the CHM")
  }
  mismatched <- terra::global(xor(is.na(dsm), is.na(chm)), "sum", na.rm = TRUE)[[1]]
  if (mismatched != 0) {
    stop("incoming DSM empty cells differ from the CHM")
  }
  invisible(TRUE)
}

validate_incoming_artifacts <- function(directory, dataset_id, dimensions,
                                        memory_budget_gib, source_point_count,
                                        point_cloud_name, sensor, country, voxel_counts = FALSE,
                                        memory_limit_gib = memory_budget_gib,
                                        denoising = FALSE, dtm_resolution = 1,
                                        chm_resolution = 0.5,
                                        voxel_count_resolution = 0.5,
                                        debug = TRUE) {
  matches_label <- function(values, expected) {
    if (identical(expected, "NA")) {
      return(all(is.na(values) | as.character(values) == "NA"))
    }
    all(!is.na(values) & as.character(values) == expected)
  }
  paths <- c(expected_artifacts(directory, dataset_id, dimensions, denoising, debug),
             if (voxel_counts) validate_voxel_count_artifacts(directory, dataset_id, voxel_count_resolution))
  absent <- paths[!file.exists(paths) | file.info(paths)$size <= 0]
  if (length(absent) > 0L) {
    stop(sprintf(
      "incoming artifact validation failed; missing/empty: %s",
      paste(basename(absent), collapse = ", ")
    ))
  }
  if (denoising) {
    status <- system2("python3", c(
      shQuote(file.path(script_directory, "denoising_validation.py")),
      shQuote(file.path(directory, paste0(dataset_id, "_denoising.json"))),
      shQuote(file.path(directory, paste0(dataset_id, "_noise_points.csv"))),
      "--point-count", format(source_point_count, scientific = FALSE, trim = TRUE)
    ))
    if (status != 0L) stop("incoming denoising accounting validation failed")
  }
  # Written by every run; published (and checked) only with --debug.
  performance <- fread(
    file.path(directory, paste0(dataset_id, "_performance.csv"))
  )
  if (debug) {
    validate_performance_resources(performance, memory_budget_gib,
                                   source_point_count, memory_limit_gib)
  }
  dtm <- terra::rast(file.path(directory, paste0(dataset_id, "_dtm.tif")))
  chm <- terra::rast(file.path(directory, paste0(dataset_id, "_chm.tif")))
  if (terra::inMemory(dtm) || terra::inMemory(chm) ||
      any(abs(terra::res(dtm) - dtm_resolution) > 1e-9) ||
      any(abs(terra::res(chm) - chm_resolution) > 1e-9) ||
      !identical(terra::crs(dtm), terra::crs(chm))) {
    stop("incoming raster resolution, storage mode, or CRS validation failed")
  }
  validate_surface_rasters(
    chm,
    terra::rast(file.path(directory, paste0(dataset_id, "_dsm.tif")))
  )
  for (dimension in dimensions) {
    results <- fread(file.path(
      directory,
      sprintf("%s_%s_results.csv", dataset_id, dimension)
    ))
    diagnostics <- fread(file.path(
      directory,
      sprintf("%s_%s_segment_diagnostics.csv", dataset_id, dimension)
    ))
    required_results <- names(empty_result_table())
    required_diagnostics <- names(empty_segment_diagnostics())
    if (!identical(names(results), required_results) ||
        !identical(names(diagnostics), required_diagnostics)) {
      stop(sprintf(
        "incoming CSV schema differs from the published schema for %s",
        dimension
      ))
    }
    if (nrow(results) != performance$tile_count[[1]] ||
        any(results$file != point_cloud_name) ||
        !matches_label(results$sensor, sensor) ||
        !matches_label(results$country, country) ||
        any(diagnostics$file != point_cloud_name) ||
        !matches_label(diagnostics$sensor, sensor) ||
        !matches_label(diagnostics$country, country)) {
      stop(sprintf(
        "incoming row count or metadata validation failed for %s",
        dimension
      ))
    }
    tile_features <- st_read(
      file.path(
        directory,
        sprintf("%s_%s_tiles.geojson", dataset_id, dimension)
      ),
      quiet = TRUE,
      stringsAsFactors = FALSE
    )
    if (nrow(tile_features) != nrow(results)) {
      stop(sprintf("incoming tile GeoJSON row count differs for %s", dimension))
    }
  }
  invisible(paths)
}
