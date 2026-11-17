# Memory-bounded global segment aggregation and finalization. Sourced by
# analysis.R; uses its constants and helpers.

# Memory-bounded global segment aggregation --------------------------------

available_extra_dimensions <- function(point_cloud) {
  header <- readLASheader(point_cloud)
  extra_bytes <- header@VLR$Extra_Bytes[["Extra Bytes Description"]]
  if (is.null(extra_bytes) || length(extra_bytes) == 0) return(character())
  unique(vapply(extra_bytes, function(description) {
    as.character(description$name)
  }, character(1)))
}

select_instance_dimensions <- function(point_cloud, candidates) {
  available <- available_extra_dimensions(point_cloud)
  unique(candidates[candidates %in% available])
}

segment_catalog_selection <- function(point_cloud, instance_dimensions,
                                      base_selection = SEGMENT_CATALOG_BASE_SELECTION) {
  header <- readLASheader(point_cloud)
  extra_bytes <- header@VLR$Extra_Bytes[["Extra Bytes Description"]]
  if (is.null(extra_bytes) || length(extra_bytes) == 0) {
    return(base_selection)
  }
  available <- vapply(extra_bytes, function(description) {
    as.character(description$name)
  }, character(1))
  segment_selection_from_names(available, instance_dimensions, base_selection)
}

# Convert named ExtraByte requirements into lidR's compact LAS selection syntax.
# lidR addresses positions 1-9 individually; selector 0 means all ExtraBytes.
catalog_selection_for_dimensions <- function(point_cloud, dimensions,
                                             base_selection = "xyz") {
  # NULL becomes character(0); remove absent placeholders passed by LAZ mode.
  dimensions <- dimensions[!is.na(dimensions) & nzchar(dimensions)]
  # Avoid reading any ExtraByte when the caller needs only standard XYZ fields.
  if (length(dimensions) == 0L) return(base_selection)
  # Inspect VLR definitions without loading point records.
  header <- readLASheader(point_cloud)
  extra_bytes <- header@VLR$Extra_Bytes[["Extra Bytes Description"]]
  if (is.null(extra_bytes) || length(extra_bytes) == 0L) {
    stop("requested catalog ExtraByte dimensions are absent")
  }
  # Preserve declared ordinal order because selectors are position-based.
  available <- vapply(extra_bytes, function(description) {
    as.character(description$name)
  }, character(1))
  # A missing required field would make later order restoration/science unsafe.
  missing <- setdiff(dimensions, available)
  if (length(missing) > 0L) {
    stop(sprintf(
      "catalog is missing required ExtraByte dimension(s): %s",
      paste(missing, collapse = ", ")
    ))
  }
  # Translate names to one-based lidR ExtraByte selector positions.
  positions <- unique(match(dimensions, available))
  # Position 0 requests all ExtraBytes when any requirement lies beyond 9.
  if (any(positions > 9L)) return(paste0(base_selection, "0"))
  # Concatenate individual positions for the minimal selective read.
  paste0(base_selection, paste(positions, collapse = ""))
}

segment_selection_from_names <- function(available, instance_dimensions,
                                         base_selection = SEGMENT_CATALOG_BASE_SELECTION) {
  positions <- unique(match(instance_dimensions, available))
  positions <- positions[!is.na(positions)]
  if (length(positions) == 0) return(base_selection)

  # lidR exposes individual ExtraBytes selectors only for positions 1-9.
  # Position 0 requests all ExtraBytes (every one with the patched rlas) while
  # still excluding unrelated standard LAS attributes, so dimensions beyond the
  # ninth field are read from the source itself, COPC included.
  if (any(positions > 9L)) return(paste0(base_selection, "0"))
  paste0(base_selection, paste(positions, collapse = ""))
}

segment_bucket <- function(instance_id, bucket_count) {
  as.integer((as.numeric(instance_id) %% bucket_count) + 1L)
}

write_segment_chunk_store <- function(voxels, layers, apex, work_directory,
                                      bucket_count) {
  voxels[, bucket := segment_bucket(instance_id, bucket_count)]
  layers[, bucket := segment_bucket(instance_id, bucket_count)]
  apex[, bucket := segment_bucket(instance_id, bucket_count)]
  paths <- character()
  voxel_columns <- c(
    "instance_dimension", "instance_id", "voxel_x", "voxel_y", "voxel_z",
    intersect("voxel_first_order", names(voxels))
  )
  apex_columns <- c(
    "instance_dimension", "instance_id", "apex_x", "apex_y", "apex_z",
    intersect("apex_order", names(apex))
  )
  for (bucket_index in sort(unique(voxels$bucket))) {
    path <- tempfile(
      pattern = sprintf("segment_%d_%03d_", Sys.getpid(), bucket_index),
      tmpdir = work_directory,
      fileext = ".rds"
    )
    saveRDS(
      list(
        voxels = voxels[bucket == bucket_index, ..voxel_columns],
        layers = layers[bucket == bucket_index, .(
          instance_dimension, instance_id, occupied_layer
        )],
        apex = apex[bucket == bucket_index, ..apex_columns]
      ),
      path,
      compress = FALSE
    )
    paths[[as.character(bucket_index)]] <- path
  }
  paths
}

# Read and aggregate one bounded segmentation chunk.  COPC order restoration is
# intentionally performed before instance filtering and apex tie-breaking.
segment_chunk <- function(chunk, dtm_path, voxel_resolution, maximum_height,
                          instance_dimensions, work_directory, bucket_count,
                          strict_dimensions = FALSE,
                          point_order_dimension = NULL,
                          observation_directory = NULL,
                          target_chunk_points = 0, denoising = NULL,
                          boundary_read_adjustment = NULL) {
  configure_catalog_worker()
  las <- keep_requested_extra_bytes(
    readLAS(adjust_catalog_boundary_read(chunk, boundary_read_adjustment)),
    c(instance_dimensions, point_order_dimension))
  record_chunk_observation(
    observation_directory,
    "segment",
    nrow(las@data),
    target_chunk_points
  )
  if (is.empty(las)) {
    if (!chunk_has_valid_extent(chunk)) return(NULL)
    return(list(paths = NULL))
  }
  # Recover original LAZ record order inside this spatially streamed chunk.
  las <- restore_original_point_order(
    las,
    point_order_dimension,
    "segment"
  )
  if (isTRUE(strict_dimensions)) {
    assert_loaded_dimensions(
      las,
      c(
        "X", "Y", "Z", "buffer", instance_dimensions,
        # Strict COPC mode proves that selector construction loaded the key too.
        point_order_dimension
      ),
      "segment"
    )
  }
  selected <- instance_dimensions[instance_dimensions %in% names(las@data)]
  if (length(selected) == 0) return(list(paths = NULL))
  las <- apply_denoising_mask(las, denoising)
  if (is.empty(las)) return(list(paths = NULL))
  dtm <- terra::rast(dtm_path)
  las <- normalize_height(las, dtm)
  las_data <- as.data.table(las@data)
  dimension_paths <- lapply(selected, function(instance_dimension) {
    source <- data.table(
      X = las_data$X,
      Y = las_data$Y,
      Z = las_data$Z,
      instance_id = las_data[[instance_dimension]],
      source_order = if (is.null(point_order_dimension)) {
        seq_len(nrow(las_data))
      } else {
        las_data[[point_order_dimension]]
      }
    )
    source <- source[
      is.finite(instance_id) & instance_id > 0 &
        is.finite(Z) & Z >= 0 & Z <= maximum_height
    ]
    if (nrow(source) == 0) return(character())

    segment_keys <- "instance_id"
    apex <- source[
      source[, .I[which.max(Z)], by = segment_keys]$V1,
      .(
        instance_id,
        apex_x = X,
        apex_y = Y,
        apex_z = Z,
        apex_order = source_order
      )
    ]
    source[, `:=`(
      voxel_x = floor(X / voxel_resolution),
      voxel_y = floor(Y / voxel_resolution),
      voxel_z = floor(Z / voxel_resolution)
    )]

    voxels <- unique(source[, .(
      instance_id,
      voxel_x,
      voxel_y,
      voxel_z,
      voxel_first_order = source_order
    )], by = c("instance_id", "voxel_x", "voxel_y", "voxel_z"))
    layers <- unique(source[, .(
      instance_id,
      occupied_layer = floor(Z)
    )])
    voxels[, instance_dimension := instance_dimension]
    setcolorder(voxels, c("instance_dimension", setdiff(names(voxels), "instance_dimension")))
    layers[, instance_dimension := instance_dimension]
    setcolorder(layers, c("instance_dimension", setdiff(names(layers), "instance_dimension")))
    apex[, instance_dimension := instance_dimension]
    setcolorder(apex, c("instance_dimension", setdiff(names(apex), "instance_dimension")))

    # Plain LAS/LAZ retains the legacy compact-summary schema and reducer order.
    # Only indexed COPC needs global source ordinals to neutralize spatial
    # reordering and adaptive window-boundary changes.
    if (is.null(point_order_dimension)) {
      voxels[, voxel_first_order := NULL]
      apex[, apex_order := NULL]
    }

    paths <- write_segment_chunk_store(
      voxels,
      layers,
      apex,
      work_directory,
      bucket_count
    )
    rm(source, voxels, layers, apex)
    gc(verbose = FALSE)
    paths
  })
  names(dimension_paths) <- selected
  paths <- do.call(c, dimension_paths)
  # Preserve numeric bucket names for the existing disk reducer.
  names(paths) <- unlist(lapply(dimension_paths, names), use.names = FALSE)
  # The compact RDS summaries now own everything needed by the disk reducer.
  # Release the normalized LAS and its table view before another window loads.
  rm(las, las_data, dtm)
  gc(verbose = FALSE)
  list(paths = paths, dimension_paths = dimension_paths)
}

# Configure a spatially streamed segment catalog.  The order key travels beside
# requested instance dimensions but is not included in scientific aggregation.
accumulate_segments <- function(point_cloud, dtm_path, chunk_size,
                                voxel_resolution, maximum_height,
                                instance_dimensions, workers, work_directory,
                                bucket_count = SEGMENT_BUCKET_COUNT,
                                base_selection = SEGMENT_CATALOG_BASE_SELECTION,
                                point_order_dimension = NULL,
                                target_chunk_points = 0, memory_policy = NULL,
                                denoising = NULL) {
  prepare_raster_work_directory(work_directory, "segment")
  observation_directory <- file.path(work_directory, "observations")
  dir.create(observation_directory)
  catalog <- readLAScatalog(point_cloud)
  if (!is.null(point_order_dimension)) {
    catalog <- catalog_with_source_boundaries(catalog, point_cloud)
  }
  opt_chunk_size(catalog) <- chunk_size
  opt_chunk_alignment(catalog) <- c(0, 0)
  opt_chunk_buffer(catalog) <- 0
  # Request both the segmentation values and source-row order in one COPC scan.
  selection <- segment_catalog_selection(
    point_cloud,
    c(instance_dimensions, point_order_dimension),
    base_selection
  )
  opt_select(catalog) <- selection
  opt_progress(catalog) <- FALSE
  # The planner adds the columns this selection brings (phase_bytes_per_point).
  if (!is.null(memory_policy)) memory_policy$selections$segment <- selection
  boundary_read_adjustment <- if (is.null(point_order_dimension)) NULL else
    catalog_boundary_read_adjustment(point_cloud)
  segment_options <- list(automerge = FALSE, drop_null = FALSE)
  dense_plan <- if (is.null(memory_policy)) list() else plan_dense_reads(
    scheduled_chunks(catalog, segment_options), boundary_read_adjustment, memory_policy, "segment")
  task_windows_fn <- if (length(dense_plan)) {
    function(chunk) {
      rectangles <- dense_plan[[cluster_key(chunk)]]
      reads <- if (is.null(rectangles)) list(chunk) else
        lapply(rectangles, rectangle_cluster, parent = chunk)
      lapply(reads, function(read)
        catalog_task_window(adjust_catalog_boundary_read(read, boundary_read_adjustment)))
    }
  } else if (!is.null(boundary_read_adjustment)) {
    function(chunk)
      catalog_task_window(adjust_catalog_boundary_read(chunk, boundary_read_adjustment))
  } else {
    NULL
  }
  callback <- function(chunk, ...) {
    rectangles <- if (length(dense_plan)) dense_plan[[cluster_key(chunk)]]
    if (is.null(rectangles)) return(segment_chunk(chunk, ...))
    combine_segment_parts(lapply(rectangles, function(rectangle)
      segment_chunk(rectangle_cluster(chunk, rectangle), ...)))
  }
  message(sprintf(
    "Reading segment dimensions with LAS selection %s",
    selection
  ))
  chunks <- apply_catalog_chunks(
      catalog, callback, workers, "segment",
      memory_policy = memory_policy,
      task_windows_fn = task_windows_fn,
      dtm_path = dtm_path,
      voxel_resolution = voxel_resolution,
      maximum_height = maximum_height,
      instance_dimensions = instance_dimensions,
      work_directory = work_directory,
      bucket_count = bucket_count,
      strict_dimensions = identical(
        base_selection,
        SEGMENT_SELECTION
      ),
      # The worker uses this only to sort its in-memory bounded chunk.
      point_order_dimension = point_order_dimension,
      observation_directory = observation_directory,
      target_chunk_points = target_chunk_points,
      denoising = denoising,
      boundary_read_adjustment = boundary_read_adjustment,
      .options = segment_options
  )
  assert_catalog_completed(chunks, "segment")
  dimension_chunks <- setNames(lapply(instance_dimensions, function(dimension) {
    paths <- lapply(chunks, function(chunk) chunk$dimension_paths[[dimension]])
    Filter(function(paths) length(paths) > 0L, paths)
  }), instance_dimensions)
  chunks <- lapply(chunks, `[[`, "paths")
  chunks <- Filter(function(paths) length(paths) > 0L, chunks)
  observations <- summarize_chunk_observations(observation_directory)
  if (length(chunks) == 0) {
    return(list(chunks = list(), dimension_chunks = dimension_chunks,
                bucket_count = bucket_count,
                observations = observations))
  }
  list(
    chunks = chunks,
    dimension_chunks = dimension_chunks,
    bucket_count = bucket_count,
    observations = observations
  )
}

finalize_segments <- function(accumulated, parameters) {
  voxels <- accumulated$voxels
  layers <- accumulated$layers
  apex <- accumulated$apex

  segment_keys <- c("instance_dimension", "instance_id")
  voxel_metrics <- voxels[, .(n_vox = .N), by = segment_keys]
  voxel_metrics[, voxel_volume := n_vox * parameters$voxel_resolution^3]
  crown_metrics <- voxels[, .(
    crown_area = uniqueN(data.table(voxel_x, voxel_y)) *
      parameters$voxel_resolution^2
  ), by = segment_keys]
  layer_metrics <- layers[, .(
    n_occupied_layers = uniqueN(occupied_layer)
  ), by = segment_keys]
  extent_metrics <- voxels[, {
    if (.N < 3) {
      .(pca_extent_1 = NA_real_, pca_extent_2 = NA_real_, pca_extent_3 = NA_real_)
    } else {
      assert_memory_budget(.N * 3 * 8 * 6, parameters, "segment PCA", TRUE)
      coordinates <- as.matrix(.SD) * parameters$voxel_resolution
      extents <- tryCatch({
        rotated <- prcomp(coordinates, center = TRUE)$x
        apply(rotated, 2, function(values) diff(range(values)))
      }, error = function(error) rep(NA_real_, 3))
      .(
        pca_extent_1 = extents[[1]],
        pca_extent_2 = extents[[2]],
        pca_extent_3 = extents[[3]]
      )
    }
  }, by = segment_keys, .SDcols = c("voxel_x", "voxel_y", "voxel_z")]

  segments <- Reduce(
    function(left, right) merge(left, right, by = segment_keys, all = TRUE),
    list(apex, voxel_metrics, crown_metrics, layer_metrics, extent_metrics)
  )
  segments[, `:=`(
    pass_voxels = n_vox >= parameters$minimum_tree_voxels,
    pass_apex = apex_z > parameters$apex_minimum_height,
    pass_thickness = !is.na(pca_extent_3) &
      pca_extent_3 >= parameters$minimum_tree_thickness,
    pass_occupied_layers = n_occupied_layers >= parameters$minimum_occupied_layers
  )]
  segments[, is_tree := pass_voxels & pass_apex &
    pass_thickness & pass_occupied_layers]
  segments[, fail_reason := paste0(
    ifelse(!pass_voxels, "voxels ", ""),
    ifelse(!pass_apex, "apex ", ""),
    ifelse(!pass_thickness, "thickness ", ""),
    ifelse(!pass_occupied_layers, "layers ", "")
  )]
  segments
}

# Merge disk pieces in source-order semantics rather than worker completion
# order.  Internal ordinal fields are removed before scientific output schemas
# are assembled, so callers see the historical columns unchanged.
merge_segment_pieces <- function(pieces) {
  voxels <- rbindlist(lapply(pieces, function(piece) piece[["voxels"]]),
                      use.names = TRUE, fill = TRUE)
  voxel_keys <- c(
    "instance_dimension", "instance_id", "voxel_x", "voxel_y", "voxel_z"
  )
  if ("voxel_first_order" %in% names(voxels)) {
    setorderv(voxels, c(voxel_keys, "voxel_first_order"))
    voxels <- unique(voxels, by = voxel_keys)
    setorderv(voxels, c("instance_dimension", "instance_id", "voxel_first_order"))
  } else {
    voxels <- unique(voxels, by = voxel_keys)
  }

  layers <- unique(rbindlist(
    lapply(pieces, function(piece) piece[["layers"]]),
    use.names = TRUE,
    fill = TRUE
  ))
  apex <- rbindlist(lapply(pieces, function(piece) piece[["apex"]]),
                    use.names = TRUE, fill = TRUE)
  apex_keys <- c("instance_dimension", "instance_id")
  if ("apex_order" %in% names(apex)) {
    setorderv(apex, c(apex_keys, "apex_z", "apex_order"), c(1L, 1L, -1L, 1L))
    apex <- unique(apex, by = apex_keys)
    apex[, apex_order := NULL]
  } else {
    apex <- apex[
      apex[, .I[which.max(apex_z)], by = apex_keys]$V1
    ]
  }
  list(voxels = voxels, layers = layers, apex = apex)
}

segment_partition_estimate <- function(paths) {
  sizes <- file.info(paths)$size
  sizes[!is.finite(sizes)] <- 0
  # Uncompressed RDS expands into list/data.table objects and is copied by
  # rbindlist(), unique(), and PCA preparation. Splitting more often is safe.
  sum(sizes) * 8
}

load_segment_partition <- function(paths, parameters) {
  assert_memory_budget(
    segment_partition_estimate(paths),
    parameters,
    "segment partition aggregation"
  )
  pieces <- lapply(paths, readRDS)
  merged <- merge_segment_pieces(pieces)
  voxels <- merged$voxels
  layers <- merged$layers
  apex <- merged$apex
  result <- finalize_segments(
    list(voxels = voxels, layers = layers, apex = apex),
    parameters
  )
  rm(pieces, voxels, layers, apex)
  gc(verbose = FALSE)
  result
}

split_segment_partition <- function(paths, work_directory, depth, parameters) {
  children <- list(character(), character())
  keys <- character()
  divisor <- SEGMENT_BUCKET_COUNT * 2^depth
  for (path in paths) {
    assert_memory_budget(
      file.info(path)$size * 8,
      parameters,
      "one segment chunk summary"
    )
    piece <- readRDS(path)
    piece_keys <- unique(rbindlist(lapply(piece, function(values) {
      values[, .(instance_dimension, instance_id)]
    }), use.names = TRUE))
    keys <- unique(c(
      keys,
      paste(piece_keys$instance_dimension, piece_keys$instance_id)
    ))
    for (branch in 0:1) {
      subset_piece <- lapply(piece, function(values) {
        values[(floor(abs(as.numeric(instance_id)) / divisor) %% 2) == branch]
      })
      if (!any(vapply(subset_piece, nrow, integer(1)) > 0L)) next
      child_path <- tempfile(
        pattern = sprintf("segment_split_%02d_%d_", depth, branch),
        tmpdir = work_directory,
        fileext = ".rds"
      )
      saveRDS(subset_piece, child_path, compress = FALSE)
      children[[branch + 1L]] <- c(children[[branch + 1L]], child_path)
    }
    rm(piece, piece_keys)
    gc(verbose = FALSE)
  }
  list(paths = Filter(length, children), unique_keys = length(keys))
}

finalize_segment_partition <- function(paths, parameters, depth = 0L) {
  budget_gib <- parameter_or(parameters, "memory_budget_gib", Inf)
  estimate <- segment_partition_estimate(paths)
  if (!is.finite(budget_gib) ||
      current_rss_bytes() + estimate <= budget_gib * 1024^3) {
    return(load_segment_partition(paths, parameters))
  }

  split <- split_segment_partition(
    paths,
    dirname(paths[[1]]),
    depth,
    parameters
  )
  child_paths <- unlist(split$paths, use.names = FALSE)
  on.exit(unlink(child_paths, force = TRUE), add = TRUE)
  if (split$unique_keys <= 1L || depth >= 52L) {
    assert_memory_budget(
      estimate,
      parameters,
      "one complete segment",
      single_segment = TRUE
    )
  }
  message(sprintf(
    "Splitting oversized segment partition at depth %d into %d child partition(s)",
    depth + 1L,
    length(split$paths)
  ))
  results <- lapply(split$paths, function(child) {
    finalize_segment_partition(child, parameters, depth + 1L)
  })
  rbindlist(results, use.names = TRUE)
}

finalize_segment_store_in_memory <- function(store, parameters) {
  paths <- unlist(store$chunks, use.names = FALSE)
  if (length(paths) == 0L) return(NULL)
  pieces <- lapply(paths, readRDS)
  merged <- merge_segment_pieces(pieces)
  voxels <- merged$voxels
  layers <- merged$layers
  apex <- merged$apex
  segments <- finalize_segments(
    list(voxels = voxels, layers = layers, apex = apex),
    parameters
  )
  setkeyv(segments, c("instance_dimension", "instance_id"))
  segments
}

# The most workers whose largest bucket estimates fit the budget together.
segment_bucket_workers <- function(estimates, threads, available) {
  largest <- sort(estimates, decreasing = TRUE)
  fits <- cumsum(largest) <= available
  max(1L, min(as.integer(threads), length(largest), sum(fits)))
}

finalize_segment_store <- function(store, parameters,
                                   adapter = c("disk", "memory")) {
  adapter <- match.arg(adapter)
  if (identical(adapter, "memory")) {
    return(finalize_segment_store_in_memory(store, parameters))
  }
  previous_threads <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(previous_threads), add = TRUE)
  data.table::setDTthreads(DATA_TABLE_THREADS)
  buckets <- lapply(seq_len(store$bucket_count), function(bucket_index) {
    key <- as.character(bucket_index)
    unlist(lapply(store$chunks, function(chunk_paths) {
      unname(chunk_paths[names(chunk_paths) == key])
    }), use.names = FALSE)
  })
  buckets <- Filter(length, buckets)
  # Buckets hold disjoint segments: finalize them in forked workers, as many
  # as their largest estimates allow within the budget; results keep bucket
  # order, so the table is the same as from one process.
  estimates <- vapply(buckets, segment_partition_estimate, numeric(1))
  available <- parameter_or(parameters, "memory_budget_gib", Inf) * 1024^3 - current_rss_bytes()
  workers <- segment_bucket_workers(estimates, parameter_or(parameters, "threads", 1L), available)
  worker_parameters <- parameters
  if (is.finite(available)) {
    worker_parameters$memory_budget_gib <- (current_rss_bytes() + available / workers) / 1024^3
  }
  finalized <- run_bounded_tasks(seq_along(buckets), function(index) {
    result <- finalize_segment_partition(buckets[[index]], worker_parameters)
    gc(verbose = FALSE)
    result
  }, workers)
  if (any(vapply(finalized, inherits, logical(1), "try-error"))) {
    stop("segment finalization worker failed")
  }
  finalized <- Filter(Negate(is.null), finalized)
  if (length(finalized) == 0) return(NULL)
  segments <- rbindlist(finalized, use.names = TRUE)
  setkeyv(segments, c("instance_dimension", "instance_id"))
  segments
}

# A chunk writes each dimension to separate RDS files already. Carry their path
# ownership through catalog collection so reducing SAT never loads FM summaries.
# This metadata is bounded by chunk/bucket counts, not by point or voxel counts.
segment_store_for_dimension <- function(store, instance_dimension) {
  chunks <- store$dimension_chunks[[instance_dimension]]
  if (is.null(chunks)) {
    stop(sprintf("segment store has no path ownership for %s", instance_dimension))
  }
  list(chunks = chunks, bucket_count = store$bucket_count)
}

# The instance dimensions of one read have disjoint stores: finalize them in
# forked workers, each with its share of the threads and of the memory left,
# and return them by name in read order, so the tables equal sequential ones.
# Each worker removes its dimension's summaries once its table is built.
finalize_dimensions <- function(accumulated, dimensions, parameters, tiles) {
  finalize_one <- function(dimension, dimension_parameters) {
    store <- segment_store_for_dimension(accumulated, dimension)
    segments <- assign_segments_to_tiles(
      finalize_segment_store(store, dimension_parameters), tiles
    )
    unlink(unlist(store$chunks, use.names = FALSE), force = TRUE)
    gc(verbose = FALSE)
    list(segments = segments)
  }
  workers <- max(1L, min(length(dimensions),
    as.integer(parameter_or(parameters, "threads", 1L))))
  dimension_parameters <- parameters
  if (workers > 1L) {
    dimension_parameters$threads <- max(1L,
      as.integer(parameter_or(parameters, "threads", 1L)) %/% workers)
    budget_gib <- parameter_or(parameters, "memory_budget_gib", Inf)
    if (is.finite(budget_gib)) {
      rss <- current_rss_bytes()
      dimension_parameters$memory_budget_gib <-
        (rss + (budget_gib * 1024^3 - rss) / workers) / 1024^3
    }
  }
  # Forked workers inherit no multi-threaded OpenMP pool (see tile_worker_plan).
  previous_threads <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(previous_threads), add = TRUE)
  data.table::setDTthreads(DATA_TABLE_THREADS)
  results <- run_bounded_tasks(seq_along(dimensions), function(index) {
    finalize_one(dimensions[[index]], dimension_parameters)
  }, workers)
  for (index in seq_along(dimensions)) {
    result <- results[[index]]
    if (inherits(result, "try-error")) {
      stop(conditionMessage(attr(result, "condition")), call. = FALSE)
    }
    if (!is.list(result) || !"segments" %in% names(result)) {
      stop(sprintf("segment finalization of %s did not return (worker killed?)",
                   dimensions[[index]]))
    }
  }
  setNames(lapply(results, `[[`, "segments"), dimensions)
}
