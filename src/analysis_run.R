# Orchestration (run_analysis) and the optional performance report. Sourced by
# analysis.R; uses its constants and helpers.

# Orchestration and optional performance report ----------------------------

elapsed_seconds <- function(started_at) {
  unname(proc.time()[["elapsed"]] - started_at)
}

peak_rss_mib <- function() {
  status <- tryCatch(readLines("/proc/self/status", warn = FALSE),
    error = function(error) character())
  line <- grep("^VmHWM:", status, value = TRUE)
  if (length(line) == 0) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", line[[1]])) / 1024
}

cgroup_memory_peak_mib <- function() {
  candidates <- c(
    "/sys/fs/cgroup/memory.peak",
    "/sys/fs/cgroup/memory/memory.max_usage_in_bytes"
  )
  path <- candidates[file.exists(candidates)][1]
  if (is.na(path)) return(NA_real_)
  value <- suppressWarnings(as.numeric(readLines(path, n = 1L, warn = FALSE)))
  if (!is.finite(value)) NA_real_ else value / 1024^2
}

temporary_directory_bytes <- function(path) {
  if (!dir.exists(path)) return(0)
  paths <- list.files(
    path,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    no.. = TRUE
  )
  if (length(paths) == 0L) return(0)
  # The staged source cloud is a symlink and is not temporary disk usage.
  paths <- paths[!nzchar(Sys.readlink(paths))]
  sizes <- file.info(paths)$size
  sum(sizes[is.finite(sizes)], na.rm = TRUE)
}

write_performance_report <- function(parameters, footprint_source, point_count, tile_count,
                                     instance_dimensions, catalog_workers,
                                     timings, path) {
  report <- data.table(
    point_cloud = basename(parameters$point_cloud),
    footprint_source = footprint_source,
    point_count = point_count,
    tile_count = tile_count,
    instance_dimension = if (length(instance_dimensions) == 0) {
      NA_character_
    } else {
      paste(instance_dimensions, collapse = "|")
    },
    threads_requested = parameters$threads,
    threads_effective = get_lidr_threads(),
    catalog_workers = catalog_workers,
    peak_rss_mib = round(peak_rss_mib(), 3),
    cgroup_memory_peak_mib = round(cgroup_memory_peak_mib(), 3),
    temporary_disk_peak_mib = round(timings$temporary_disk_peak_mib, 3),
    grid_seconds = round(timings$grid, 6),
    dtm_seconds = round(timings$dtm, 6),
    chm_seconds = round(timings$chm, 6),
    segment_seconds = round(timings$segments, 6),
    tile_seconds = round(timings$tiles, 6),
    output_seconds = round(timings$outputs, 6),
    total_seconds = round(timings$total, 6),
    tile_size = parameters$tile_size,
    grid_search_step = parameters$grid_search_step,
    ptd_resolution = parameters$ptd_resolution,
    dtm_resolution = parameters$dtm_resolution,
    maximum_height = parameters$maximum_height,
    metric_height_policy = if (is.finite(parameters$maximum_height)) "explicit_legacy_ceiling" else "observed_voxel_height_v1",
    metric_height_upper_m = parameters$metric_height_upper_m,
    observed_maximum_voxel_height_m = parameters$observed_maximum_voxel_height_m,
    observed_maximum_voxel_height_m_exact = sprintf("%.17g", parameters$observed_maximum_voxel_height_m),
    voxel_resolution = parameters$voxel_resolution,
    vegetation_minimum_height = parameters$vegetation_minimum_height,
    chm_resolution = parameters$chm_resolution,
    gap_height_threshold = parameters$gap_height_threshold,
    minimum_tree_voxels = parameters$minimum_tree_voxels,
    apex_minimum_height = parameters$apex_minimum_height,
    minimum_tree_thickness = parameters$minimum_tree_thickness,
    minimum_occupied_layers = parameters$minimum_occupied_layers,
    dtm_chunk_size = parameters$dtm_chunk_size,
    chunk_size = parameters$chunk_size,
    target_chunk_points = parameters$adaptive_chunking$target_points,
    adaptive_chunking_enabled = parameters$adaptive_chunking$enabled,
    adaptive_chunking_method = parameters$adaptive_chunking$method,
    adaptive_chunk_size_m = parameters$adaptive_chunking$chunk_size_m,
    adaptive_mean_density_points_m2 = parameter_or(
      parameters$adaptive_chunking,
      "mean_density_points_m2",
      NA_real_
    ),
    adaptive_estimated_points_per_full_chunk = parameter_or(
      parameters$adaptive_chunking,
      "estimated_points_per_full_chunk",
      NA_real_
    ),
    adaptive_chm_observed_chunks = timings$chm_observations$count,
    adaptive_chm_observed_minimum_points =
      timings$chm_observations$minimum_points,
    adaptive_chm_observed_maximum_points =
      timings$chm_observations$maximum_points,
    adaptive_chm_oversized_chunks = timings$chm_observations$oversized_count,
    adaptive_segment_observed_chunks = timings$segment_observations$count,
    adaptive_segment_observed_minimum_points =
      timings$segment_observations$minimum_points,
    adaptive_segment_observed_maximum_points =
      timings$segment_observations$maximum_points,
    adaptive_segment_oversized_chunks =
      timings$segment_observations$oversized_count,
    dtm_buffer = parameters$dtm_buffer,
    dtm_storage_mode = "disk_backed_vrt",
    chm_storage_mode = "disk_backed_vrt",
    execution_mode = if (identical(parameters$segment_catalog_base_selection,
                                  SEGMENT_SELECTION)) {
      "streaming"
    } else {
      "default"
    },
    internal_memory_budget_gib = parameters$memory_budget_gib,
    segment_catalog_base_selection =
      parameters$segment_catalog_base_selection,
    instance_dimension_candidates = paste(parameters$instance_dimension, collapse = "|")
  )
  fwrite(report, path, na = "NA")
  message(sprintf(
    paste0(
      "PERF total=%.3fs grid=%.3fs dtm=%.3fs chm=%.3fs ",
      "segments=%.3fs tiles=%.3fs outputs=%.3fs peak_rss=%.1fMiB ",
      "cgroup_peak=%.1fMiB"
    ),
    timings$total, timings$grid, timings$dtm, timings$chm, timings$segments,
    timings$tiles, timings$outputs, report$peak_rss_mib[[1]],
    report$cgroup_memory_peak_mib[[1]]
  ))
  invisible(path)
}

run_analysis <- function(parameters) {
  parameters$sensor <- parameter_or(parameters, "sensor", NA_character_)
  parameters$country <- parameter_or(parameters, "country", NA_character_)
  parameters$memory_budget_gib <- parameter_or(
    parameters, "memory_budget_gib", Inf
  )
  parameters$segment_catalog_base_selection <- parameter_or(
    parameters,
    "segment_catalog_base_selection",
    SEGMENT_CATALOG_BASE_SELECTION
  )
  parameters$segment_bucket_count <- parameter_or(
    parameters,
    "segment_bucket_count",
    SEGMENT_BUCKET_COUNT
  )
  # NULL retains original LAZ behavior; COPC runners provide OriginalPointIndex.
  parameters$point_order_dimension <- parameter_or(
    parameters, "point_order_dimension", NULL
  )
  parameters$point_cloud_display_name <- parameter_or(
    parameters,
    "point_cloud_display_name",
    basename(parameters$point_cloud)
  )
  work_root <- parameter_or(parameters, "temp_dir", parameters$output_dir)
  if (!dir.exists(work_root) || file.access(work_root, mode = 2) != 0) {
    stop("analysis temporary directory must exist and be writable")
  }
  run_work_directory <- file.path(
    work_root,
    sprintf(".%s_raster_%d", parameters$dataset_id, Sys.getpid())
  )
  dir.create(run_work_directory, recursive = TRUE, showWarnings = FALSE)
  temporary_disk_peak_bytes <- temporary_directory_bytes(run_work_directory)
  on.exit(unlink(run_work_directory, recursive = TRUE, force = TRUE), add = TRUE)
  total_started <- proc.time()[["elapsed"]]
  point_count <- point_count_from_las_header(parameters$point_cloud)
  assert_lidr_point_count_supported(point_count)
  parameters$adaptive_chunking <- adaptive_chunk_plan(
    point_count,
    point_cloud_xy_bounds(parameters$point_cloud),
    parameters$chunk_size,
    parameter_or(parameters, "target_chunk_points", 0),
    parameter_or(parameters, "adaptive_spatial_input",
      grepl("\\.copc\\.laz$", parameters$point_cloud, ignore.case = TRUE))
  )
  parameters$effective_chunk_size <- parameters$adaptive_chunking$chunk_size_m
  message(sprintf(
    paste0(
      "ADAPTIVE_CHUNK enabled=%s method=%s target_points=%.0f density_points_m2=%s ",
      "maximum_chunk_size_m=%.3f selected_chunk_size_m=%.3f ",
      "estimated_points_per_full_chunk=%s"
    ),
    tolower(as.character(parameters$adaptive_chunking$enabled)),
    parameters$adaptive_chunking$method,
    parameters$adaptive_chunking$target_points,
    if (is.finite(parameter_or(
      parameters$adaptive_chunking,
      "mean_density_points_m2",
      NA_real_
    ))) {
      sprintf("%.6f", parameters$adaptive_chunking$mean_density_points_m2)
    } else {
      "NA"
    },
    parameters$chunk_size,
    parameters$effective_chunk_size,
    if (is.finite(parameters$adaptive_chunking$estimated_points_per_full_chunk)) {
      sprintf("%.0f", parameters$adaptive_chunking$estimated_points_per_full_chunk)
    } else {
      "NA"
    }
  ))
  source_crs <- parameter_or(
    parameters,
    "source_crs",
    spatial_reference_wkt(readLAScatalog(parameters$point_cloud))
  )
  paths <- artifact_paths(parameters$output_dir, parameters$dataset_id)
  segment_work_directory <- file.path(run_work_directory, "segment_chunks")
  on.exit(
    unlink(segment_work_directory, recursive = TRUE, force = TRUE),
    add = TRUE
  )
  if (parameters$threads > 0) {
    set_lidr_threads(parameters$threads)
    message(sprintf("Using %d lidR threads", get_lidr_threads()))
  } else {
    message(sprintf("Using lidR default thread count: %d", get_lidr_threads()))
  }
  thread_budget <- max(
    1L,
    as.integer(
      if (parameters$threads > 0) parameters$threads else get_lidr_threads()
    )
  )
  catalog_workers <- min(
    thread_budget,
    max(1L, as.integer(parameters$catalog_workers))
  )
  message(sprintf(
    "Using %d memory-bounded LAScatalog worker(s) within the %d-thread allowance",
    catalog_workers,
    thread_budget
  ))
  message("RASTER_STAGE stage=grid status=started")
  grid_started <- proc.time()[["elapsed"]]
  aoi <- read_analysis_footprint(parameters$point_cloud, parameters$aoi)
  if (identical(aoi$source, "point_cloud_extent")) {
    message("No Audit AOI supplied; covering the complete point-cloud XY extent")
  } else {
    message(sprintf("Using Audit AOI: %s", basename(parameters$aoi)))
  }
  tiles <- build_analysis_tiles(
    aoi,
    parameters$tile_size,
    parameters$grid_search_step,
    thread_budget
  )
  grid_seconds <- elapsed_seconds(grid_started)
  message(sprintf("Optimized grid contains %d complete Analysis Tiles", nrow(tiles)))
  message("RASTER_STAGE stage=grid status=completed")

  message("RASTER_STAGE stage=dtm status=started")
  dtm_started <- proc.time()[["elapsed"]]
  dtm_path <- paths$dtm
  reduce_points <- isTRUE(parameter_or(parameters, "reduce_points", FALSE))
  dtm_source <- if (reduce_points && dtm_needs_thinning(parameters)) {
    run_point_reduction("thin", parameters$point_cloud,
      file.path(run_work_directory, if (is.null(parameters$point_order_dimension)) "dtm_points.laz"
                else "dtm_points.copc.laz"), parameters)
  } else {
    parameters$point_cloud
  }
  dtm_policy <- parameters$memory_policy
  if (!is.null(dtm_policy)) dtm_policy$source <- dtm_source
  write_global_dtm(
    dtm_source,
    dtm_path,
    file.path(run_work_directory, "dtm_chunks"),
    parameters$dtm_chunk_size,
    parameters$dtm_buffer,
    parameters$dtm_resolution,
    parameters$ptd_resolution,
    catalog_workers,
    source_crs,
    # Propagate the key through DTM catalog selection and per-chunk restoration.
    parameters$point_order_dimension,
    memory_policy = dtm_policy,
    denoising = parameters$denoising,
    extent_cloud = parameters$point_cloud
  )
  dtm <- terra::rast(dtm_path)
  if (terra::inMemory(dtm)) {
    stop("published DTM unexpectedly materialized in memory")
  }
  dtm_seconds <- elapsed_seconds(dtm_started)
  temporary_disk_peak_bytes <- max(
    temporary_disk_peak_bytes,
    temporary_directory_bytes(run_work_directory)
  )
  message("RASTER_STAGE stage=dtm status=completed")

  if (isTRUE(parameter_or(parameters, "export_voxel_counts", FALSE))) {
    message("RASTER_STAGE stage=voxel_counts status=started")
    voxel_export <- write_voxel_counts(parameters, dtm_path, aoi, run_work_directory,
                                       workers = thread_budget)
    temporary_disk_peak_bytes <- max(temporary_disk_peak_bytes,
      temporary_directory_bytes(run_work_directory) + voxel_export$temporary_disk_bytes)
    message("RASTER_STAGE stage=voxel_counts status=completed")
  }

  instance_dimensions <- select_instance_dimensions(
    parameters$point_cloud,
    parameters$instance_dimension
  )
  analysis_cloud <- parameters$point_cloud
  analysis_policy <- parameters$memory_policy
  if (reduce_points && dense_point_cloud(parameters$point_cloud)) {
    message("RASTER_STAGE stage=point_reduction status=started")
    analysis_cloud <- reduced_analysis_cloud(parameters, tiles, dtm_path,
                                             instance_dimensions, run_work_directory)
    if (!is.null(analysis_policy)) analysis_policy$source <- analysis_cloud
    message("RASTER_STAGE stage=point_reduction status=completed")
  }

  message("RASTER_STAGE stage=chm status=started")
  chm_started <- proc.time()[["elapsed"]]
  # The DSM is a product of the production runner only.
  surface_rasters <- isTRUE(parameter_or(parameters, "surface_rasters", FALSE))
  chm_path <- paths$chm
  chm_observations <- write_global_chm(
    analysis_cloud,
    dtm_path,
    chm_path,
    file.path(run_work_directory, "chm_chunks"),
    parameters$effective_chunk_size,
    parameters$voxel_resolution,
    parameters$maximum_height,
    parameters$chm_resolution,
    catalog_workers,
    source_crs,
    # Propagate the same key through CHM selection and per-chunk restoration.
    parameters$point_order_dimension,
    if (parameters$adaptive_chunking$enabled) {
      parameters$adaptive_chunking$target_points
    } else {
      0
    },
    raster_chunk_size = parameters$chunk_size,
    memory_policy = analysis_policy,
    denoising = parameters$denoising,
    dsm_path = if (surface_rasters) paths$dsm
  )
  chm_seconds <- elapsed_seconds(chm_started)
  temporary_disk_peak_bytes <- max(
    temporary_disk_peak_bytes,
    temporary_directory_bytes(run_work_directory)
  )
  message("RASTER_STAGE stage=chm status=completed")

  parameters$observed_maximum_voxel_height_m <- chm_observations$maximum_voxel_height
  parameters$metric_height_upper_m <- metric_height_reference(
    parameters$observed_maximum_voxel_height_m,
    parameters$maximum_height, parameters$vegetation_minimum_height
  )
  message(sprintf("Metric height reference: %.10g m; normalized height ceiling: %s",
    parameters$metric_height_upper_m,
    if (is.finite(parameters$maximum_height)) as.character(parameters$maximum_height) else "none"))

  message("RASTER_STAGE stage=segmentation status=started")
  segment_started <- proc.time()[["elapsed"]]
  segments_by_dimension <- list()
  segment_observation_summaries <- list()
  if (length(instance_dimensions) == 0) {
    message("No configured Instance Dimensions found; tree and segment metrics will be NA")
  } else {
    message(sprintf(
      "Using Instance Dimensions: %s",
      paste(instance_dimensions, collapse = ", ")
    ))
    read_directory <- file.path(segment_work_directory, "read")
    message(sprintf(
      "Reading segments for %s from %s",
      paste(instance_dimensions, collapse = ", "), basename(analysis_cloud)
    ))
    # One read of the cloud supplies every instance dimension: the order is
    # restored and heights are normalized once.
    accumulated <- accumulate_segments(
      analysis_cloud,
      dtm_path,
      parameters$effective_chunk_size,
      parameters$voxel_resolution,
      parameters$maximum_height,
      instance_dimensions,
      catalog_workers,
      read_directory,
      bucket_count = parameters$segment_bucket_count,
      base_selection = parameters$segment_catalog_base_selection,
      point_order_dimension = parameters$point_order_dimension,
      target_chunk_points = if (parameters$adaptive_chunking$enabled) {
        parameters$adaptive_chunking$target_points
      } else {
        0
      },
      memory_policy = analysis_policy,
      denoising = parameters$denoising
    )
    segment_observation_summaries[[1]] <- accumulated$observations
    temporary_disk_peak_bytes <- max(
      temporary_disk_peak_bytes,
      temporary_directory_bytes(run_work_directory)
    )
    finalize_started <- proc.time()[["elapsed"]]
    segments_by_dimension <- finalize_dimensions(
      accumulated, instance_dimensions, parameters, tiles
    )
    for (instance_dimension in instance_dimensions) {
      dimension_segments <- segments_by_dimension[[instance_dimension]]
      if (!is.null(dimension_segments)) {
        message(sprintf(
          "Segment pass %s found %d segments (%d accepted trees)",
          instance_dimension, nrow(dimension_segments),
          sum(dimension_segments$is_tree, na.rm = TRUE)
        ))
      }
    }
    message(sprintf("Finalized %d segment pass(es) in %.1f s",
      length(instance_dimensions), elapsed_seconds(finalize_started)))
    rm(accumulated)
    unlink(read_directory, recursive = TRUE, force = TRUE)
    gc(verbose = FALSE)
    unlink(segment_work_directory, recursive = TRUE, force = TRUE)
  }

  segment_seconds <- elapsed_seconds(segment_started)
  segment_observations <- combine_chunk_observations(
    segment_observation_summaries
  )
  message("RASTER_STAGE stage=segmentation status=completed")

  # Without an instance dimension (none given or found), only the general products are made:
  # no tile metrics, tiles or segment diagnostics (3DT-3007).
  structural_metrics <- isTRUE(parameter_or(parameters, "structural_metrics", TRUE))
  tile_seconds <- 0
  if (structural_metrics) {
    edge_flags <- compute_edge_flags(tiles, parameters$tile_size)
    message("RASTER_STAGE stage=tiles status=started")
    tile_started <- proc.time()[["elapsed"]]
    tile_indices <- seq_len(nrow(tiles))
    # Match the decimal bounds passed to the tile reader's -inside filter.
    tile_windows <- lapply(tile_indices, function(index)
      as.numeric(sprintf("%.10f", st_bbox(tiles[index, ]))))
    # Neighbouring tiles share one read: every per-tile COPC query decoded the
    # coarse octree levels again, which dominated large areas.
    analysis_bounds <- point_cloud_xy_bounds(analysis_cloud)
    analysis_density <- point_count_from_las_header(analysis_cloud) /
      max(1, (analysis_bounds[["xmax"]] - analysis_bounds[["xmin"]]) *
             (analysis_bounds[["ymax"]] - analysis_bounds[["ymin"]]))
    tile_blocks <- plan_tile_blocks(tile_windows, parameters$tile_size,
      analysis_density, analysis_policy, workers = min(thread_budget, catalog_workers))
    block_windows <- lapply(tile_blocks, function(block) tile_block_window(tile_windows[block]))
    tile_plan <- tile_worker_plan(thread_budget, length(tile_blocks), catalog_workers)
    tile_workers <- tile_plan$process_workers
    tile_threads_per_worker <- tile_plan$threads_per_worker
    calculate_block <- function(block_index) {
      if (tile_workers > 1L) {
        set_lidr_threads(tile_threads_per_worker)
      }
      worker_dtm <- rast(dtm_path)
      block <- tile_blocks[[block_index]]
      block_las <- if (length(block) > 1L) {
        window <- block_windows[[block_index]]
        read_tile_points(analysis_cloud,
          c(xmin = window[[1]], ymin = window[[2]], xmax = window[[3]], ymax = window[[4]]),
          parameters)
      }
      block_rows <- if (!is.null(block_las)) {
        if (is.empty(block_las)) lapply(block, function(index) integer())
        else tile_point_rows(block_las, tile_windows[block])
      }
      lapply(seq_along(block), function(position) {
        index <- block[[position]]
        metrics <- calculate_tile_metrics(
          analysis_cloud,
          tiles[index, ],
          worker_dtm,
          parameters,
          block_las = block_las,
          block_rows = block_rows[[position]]
        )
        build_tile_row(
          parameters$point_cloud_display_name,
          tiles[index, ],
          edge_flags[[index]],
          nrow(tiles),
          metrics,
          empty_tree_metrics(),
          parameters$sensor,
          parameters$country
        )
      })
    }
    if (length(tile_indices) > 0L) {
      message(sprintf(
        "Calculating %d tile metric rows in %d block(s) with %d process worker(s), %d lidR thread(s) each",
        length(tile_indices), length(tile_blocks),
        tile_workers, tile_threads_per_worker
      ))
    }
    base_result <- collect_tile_block_rows(tile_indices, tile_blocks, calculate_block,
      tile_workers, task_windows = block_windows, memory_policy = analysis_policy)
    tile_seconds <- elapsed_seconds(tile_started)
    message("RASTER_STAGE stage=tiles status=completed")
  }

  message("RASTER_STAGE stage=output status=started")
  output_started <- proc.time()[["elapsed"]]
  output_paths <- list()
  if (structural_metrics) {
    output_dimensions <- if (length(instance_dimensions) == 0) {
      list(NULL)
    } else {
      as.list(instance_dimensions)
    }
    output_paths <- lapply(output_dimensions, function(instance_dimension) {
      target_dimension <- instance_dimension
      dimension_paths <- dimension_artifact_paths(paths, instance_dimension)
      dimension_segments <- if (is.null(instance_dimension)) {
        NULL
      } else {
        segments_by_dimension[[target_dimension]]
      }
      result <- result_with_tree_metrics(base_result, dimension_segments)
      write_upstream_csv(result, dimension_paths$results)
      write_segment_diagnostics(
        dimension_segments,
        parameters$point_cloud_display_name,
        dimension_paths$segment_diagnostics,
        parameters$sensor,
        parameters$country
      )
      write_tile_geojson(
        tiles,
        result,
        source_crs,
        dimension_paths$tiles_geojson
      )
      message(sprintf(
        "Wrote %d Analysis Tile rows for Instance Dimension %s to %s",
        nrow(result),
        if (is.null(instance_dimension)) "NA" else instance_dimension,
        dimension_paths$results
      ))
      dimension_paths
    })
    write_layout_png(
      aoi,
      tiles,
      edge_flags,
      parameters$point_cloud_display_name,
      parameters$tile_size,
      paths$tiles_png
    )
  }
  output_seconds <- elapsed_seconds(output_started)
  if (parameters$performance_report) {
    write_performance_report(
      parameters,
      aoi$source,
      point_count,
      nrow(tiles),
      instance_dimensions,
      catalog_workers,
      list(
        grid = grid_seconds,
        dtm = dtm_seconds,
        chm = chm_seconds,
        segments = segment_seconds,
        tiles = tile_seconds,
        outputs = output_seconds,
        temporary_disk_peak_mib = temporary_disk_peak_bytes / 1024^2,
        chm_observations = chm_observations,
        segment_observations = segment_observations,
        total = elapsed_seconds(total_started)
      ),
      paths$performance
    )
  }
  message("RASTER_STAGE stage=output status=completed")
  invisible(output_paths)
}
