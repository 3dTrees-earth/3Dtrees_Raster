# Scientific tile metrics and output adapters. Sourced by analysis.R; uses its
# constants and helpers.

# Scientific tile metrics and output adapters ------------------------------

assign_segments_to_tiles <- function(segments, tiles) {
  if (is.null(segments)) return(NULL)
  segments[, tile_id := NA_integer_]
  if (nrow(tiles) == 0) return(segments)

  bounds <- t(vapply(seq_len(nrow(tiles)), function(index) {
    as.numeric(st_bbox(tiles[index, ]))
  }, numeric(4)))
  colnames(bounds) <- c("xmin", "ymin", "xmax", "ymax")
  segments[, tile_id := mapply(function(apex_x, apex_y) {
    matches <- which(
      apex_x >= bounds[, "xmin"] & apex_x < bounds[, "xmax"] &
        apex_y >= bounds[, "ymin"] & apex_y < bounds[, "ymax"]
    )
    if (length(matches) == 0) NA_integer_ else tiles$tile_id[[matches[[1]]]]
  }, apex_x, apex_y)]
  segments
}

gini <- function(values) {
  values <- values[!is.na(values)]
  count <- length(values)
  if (count < 2) return(NA_real_)
  if (any(values < 0)) values <- values - min(values)
  values <- sort(values)
  total <- sum(values)
  if (total == 0) return(0)
  coefficient <- (2 * sum(seq_len(count) * values) / (count * total)) -
    (count + 1) / count
  coefficient * count / (count - 1)
}

empty_tree_metrics <- function() {
  list(
    n_seg_total = NA_integer_,
    n_trees = NA_integer_,
    tree_height_max = NA_real_,
    tree_height_mean = NA_real_,
    tree_height_gini = NA_real_,
    tree_crownarea_mean = NA_real_,
    tree_crownarea_max = NA_real_,
    tree_crownarea_gini = NA_real_,
    tree_volume_mean = NA_real_,
    tree_volume_max = NA_real_,
    tree_volume_gini = NA_real_
  )
}

tree_metrics_for_tile <- function(segments, tile_id) {
  if (is.null(segments)) return(empty_tree_metrics())
  target_tile_id <- tile_id
  in_tile <- segments[segments$tile_id == target_tile_id]
  trees <- in_tile[is_tree == TRUE]
  metrics <- empty_tree_metrics()
  metrics$n_seg_total <- nrow(in_tile)
  metrics$n_trees <- nrow(trees)
  if (nrow(trees) == 0) return(metrics)

  metrics$tree_height_max <- safe_round(max(trees$apex_z))
  metrics$tree_height_mean <- safe_round(mean(trees$apex_z))
  metrics$tree_height_gini <- safe_round(gini(trees$apex_z))
  metrics$tree_crownarea_mean <- safe_round(mean(trees$crown_area))
  metrics$tree_crownarea_max <- safe_round(max(trees$crown_area))
  metrics$tree_crownarea_gini <- safe_round(gini(trees$crown_area))
  metrics$tree_volume_mean <- safe_round(mean(trees$voxel_volume))
  metrics$tree_volume_max <- safe_round(max(trees$voxel_volume))
  metrics$tree_volume_gini <- safe_round(gini(trees$voxel_volume))
  metrics
}

result_with_tree_metrics <- function(base_result, segments) {
  result <- copy(base_result)
  tree_columns <- names(empty_tree_metrics())
  if (nrow(result) == 0) return(result)
  for (index in seq_len(nrow(result))) {
    metrics <- tree_metrics_for_tile(segments, result$tile_id[[index]])
    for (name in tree_columns) {
      set(result, i = index, j = name, value = metrics[[name]])
    }
  }
  result
}

empty_segment_diagnostics <- function() {
  data.table(
    file = character(),
    sensor = character(),
    country = character(),
    tile_id = integer(),
    PredInstance = numeric(),
    n_vox = integer(),
    apex_z = numeric(),
    pc_ext1 = numeric(),
    pc_ext2 = numeric(),
    pc_ext3 = numeric(),
    n_zlayer = integer(),
    pass_vox = logical(),
    pass_apex = logical(),
    pass_thick = logical(),
    pass_zlayer = logical(),
    is_tree = logical(),
    apex_in_tile = logical(),
    apex_x = numeric(),
    apex_y = numeric(),
    crown_area_m2 = numeric(),
    voxel_volume_m3 = numeric()
  )
}

write_upstream_csv <- function(values, path) {
  utils::write.csv(
    as.data.frame(values),
    path,
    row.names = FALSE,
    na = "NA"
  )
  invisible(path)
}

write_segment_diagnostics <- function(segments, point_cloud, path,
                                      sensor = NA_character_,
                                      country = NA_character_) {
  diagnostics <- empty_segment_diagnostics()
  if (!is.null(segments) && nrow(segments) > 0) {
    diagnostics <- segments[, .(
      file = basename(point_cloud),
      sensor = sensor,
      country = country,
      tile_id,
      PredInstance = instance_id,
      n_vox,
      apex_z = round(apex_z, 4),
      pc_ext1 = round(pca_extent_1, 4),
      pc_ext2 = round(pca_extent_2, 4),
      pc_ext3 = round(pca_extent_3, 4),
      n_zlayer = n_occupied_layers,
      pass_vox = pass_voxels,
      pass_apex,
      pass_thick = pass_thickness,
      pass_zlayer = pass_occupied_layers,
      is_tree,
      apex_in_tile = !is.na(tile_id),
      apex_x = round(apex_x, 4),
      apex_y = round(apex_y, 4),
      crown_area_m2 = round(crown_area, 4),
      voxel_volume_m3 = round(voxel_volume, 4)
    )]
  }
  write_upstream_csv(diagnostics, path)
}

empty_metric_values <- function(voxel_total, metric_height_upper_m) {
  list(
    vox_filled = 0L,
    vox_total = voxel_total,
    metric_height_upper_m = metric_height_upper_m,
    veg_density = NA_real_,
    zsd = NA_real_,
    zskew = NA_real_,
    zkurt = NA_real_,
    zq90 = NA_real_,
    box_dim_fixed = NA_real_,
    vci = NA_real_,
    rumple = NA_real_,
    gap_fraction = NA_real_,
    chm_sd = NA_real_,
    chm_cv = NA_real_,
    height_max = NA_real_,
    height_mean = NA_real_
  )
}

safe_round <- function(value, digits = 4) {
  if (length(value) == 0 || !is.finite(value[[1]])) return(NA_real_)
  round(as.numeric(value[[1]]), digits)
}

# A shared finite reference keeps density and VCI comparable between tiles.
# This is a metric denominator, never a point filter. Reuse the exact scalar
# maxima from the CHM pass before lidR rounds its raster to three decimals.
# No extra cloud read, raster scan or dense vertical histogram is required.
metric_height_reference <- function(observed, maximum_height = Inf,
                                    vegetation_minimum_height = 0.5) {
  if (is.finite(maximum_height)) return(maximum_height)
  if (!is.finite(observed)) {
    stop("Cannot derive the metric height reference without finite accepted voxel heights")
  }
  reference <- max(1, ceiling(observed), floor(vegetation_minimum_height) + 1)
  if (!is.finite(reference) || reference > 2^53 - 1) {
    stop("Observed canopy height exceeds exact numeric layer-index support")
  }
  reference
}

# Preserve the four used lidR::stdmetrics_z formulas without its unused dense
# entropy histogram (whose memory grows with the largest retained height).
tile_height_statistics <- function(z) {
  n <- length(z)
  centered <- z - mean(z)
  sum_squared <- sum(centered^2)
  list(
    zsd = stats::sd(z),
    zskew = (sum(centered^3) / n) / (sum_squared / n)^(3 / 2),
    zkurt = n * sum(centered^4) / sum_squared^2,
    # Reuse lidR's constructed probability, including its floating rounding.
    zq90 = unname(stats::quantile(z, seq(0.05, 0.95, 0.05)[[18L]]))
  )
}

vertical_complexity <- function(z, height_upper_m) {
  # The original script used cut(..., seq(0, upper, 1), include.lowest=TRUE), i.e.
  # right-closed one-metre layers. Count occupied layers only; empty layers
  # contribute zero entropy, but remain in the normalization denominator.
  layer_count <- floor(height_upper_m)
  if (!is.finite(layer_count) || layer_count > 2^53 - 1 || layer_count < 1) {
    stop("Metric height reference cannot represent one-metre layers exactly")
  }
  eligible <- z[is.finite(z) & z >= 0 & z <= layer_count]
  if (!length(eligible) || layer_count == 1) return(NA_real_)
  counts <- data.table(layer = pmax(1, ceiling(eligible)))[, .N, by = layer]$N
  proportions <- counts / sum(counts)
  -sum(proportions * log(proportions)) / log(layer_count)
}

# Calculate tile metrics from a spatial window.  COPC still
# supplies only the requested window; order restoration occurs after that read.
# Points of one window in source order, without points the noise filter rejected.
read_tile_points <- function(point_cloud, bounds, parameters) {
  las <- readLAS(
    point_cloud,
    # Add OriginalPointIndex to XYZ only when the active input is canonical COPC.
    select = catalog_selection_for_dimensions(
      point_cloud,
      parameters$point_order_dimension,
      "xyz"
    ),
    filter = sprintf(
      "-inside %.10f %.10f %.10f %.10f",
      bounds[["xmin"]], bounds[["ymin"]], bounds[["xmax"]], bounds[["ymax"]]
    )
  )
  if (is.empty(las)) return(las)
  las <- keep_requested_extra_bytes(las, parameters$point_order_dimension)
  # Fail closed if selective COPC reading omitted the required order key.
  assert_loaded_dimensions(
    las,
    c("X", "Y", "Z", parameters$point_order_dimension),
    "tile"
  )
  # Match the point sequence used by the original tile calculations.
  las <- restore_original_point_order(
    las,
    parameters$point_order_dimension,
    "tile"
  )
  apply_denoising_mask(las, parameters$denoising)
}

# Row numbers of each tile's points in a block read, in one pass over the
# points. LASlib's -inside keeps xmin <= X < xmax and ymin <= Y < ymax at the
# same ten-decimal bounds; grid tiles do not overlap, so each point has at most
# one tile. Rows stay in the block's restored source order.
tile_point_rows <- function(block_las, windows) {
  bounds <- t(vapply(windows, function(window) as.numeric(sprintf("%.10f", window)), numeric(4)))
  xs <- sort(unique(bounds[, 1]))
  ys <- sort(unique(bounds[, 2]))
  cell <- function(ix, iy) ix + iy * (length(xs) + 1)
  tile <- match(cell(findInterval(block_las@data$X, xs), findInterval(block_las@data$Y, ys)),
                cell(match(bounds[, 1], xs), match(bounds[, 2], ys)))
  inside <- which(!is.na(tile))
  inside <- inside[block_las@data$X[inside] < bounds[tile[inside], 3] &
                   block_las@data$Y[inside] < bounds[tile[inside], 4]]
  split(inside, factor(tile[inside], levels = seq_along(windows)))
}

# The same subset filter_poi() would return for these rows.
tile_points_from_block <- function(block_las, rows) {
  if (length(rows) == npoints(block_las)) block_las else block_las[rows]
}

calculate_tile_metrics <- function(point_cloud, tile, dtm, parameters, block_las = NULL,
                                   block_rows = NULL) {
  bounds <- st_bbox(tile)
  height_reference <- parameter_or(parameters, "metric_height_upper_m", parameters$maximum_height)
  voxel_total <- round(
    (parameters$tile_size / parameters$voxel_resolution)^2 *
      ((height_reference - parameters$vegetation_minimum_height) /
         parameters$voxel_resolution)
  )
  if (!is.finite(voxel_total) || voxel_total <= 0 || voxel_total > 2^53 - 1) {
    stop("Metric voxel volume exceeds exact numeric support or is not positive")
  }
  empty <- empty_metric_values(voxel_total, height_reference)

  las <- if (is.null(block_las)) {
    read_tile_points(point_cloud, bounds, parameters)
  } else {
    tile_points_from_block(block_las, block_rows)
  }
  if (is.empty(las)) return(empty)

  las <- normalize_height(las, dtm)
  las <- filter_poi(las, is.finite(Z) & Z >= 0 & Z <= parameters$maximum_height)
  if (is.empty(las)) return(empty)

  voxels <- voxelize_points(las, res = parameters$voxel_resolution)
  if (is.empty(voxels)) return(empty)

  height_max <- max(voxels@data$Z, na.rm = TRUE)
  if (!is.finite(parameters$maximum_height) && height_max > height_reference) {
    stop("Observed tile voxel height exceeds the shared metric height reference")
  }
  height_mean <- mean(voxels@data$Z, na.rm = TRUE)
  chm <- rasterize_canopy(
    voxels,
    res = parameters$chm_resolution,
    algorithm = p2r()
  )
  chm_values <- terra::values(chm, mat = FALSE)
  observed <- chm_values[is.finite(chm_values)]
  chm_mean <- if (length(observed)) mean(observed) else NA_real_
  chm_sd <- if (length(observed) > 1) stats::sd(observed) else NA_real_
  gap_fraction <- if (length(observed)) {
    mean(observed < parameters$gap_height_threshold)
  } else {
    NA_real_
  }
  rumple <- tryCatch(rumple_index(chm), error = function(error) NA_real_)

  vegetation <- filter_poi(voxels, Z > parameters$vegetation_minimum_height)
  if (is.empty(vegetation)) {
    empty$height_max <- safe_round(height_max)
    empty$height_mean <- safe_round(height_mean)
    empty$rumple <- safe_round(rumple)
    empty$gap_fraction <- safe_round(gap_fraction)
    empty$chm_sd <- safe_round(chm_sd)
    empty$chm_cv <- safe_round(chm_sd / chm_mean)
    return(empty)
  }

  vegetation_data <- as.data.table(vegetation@data)
  filled <- nrow(vegetation_data)
  z_metrics <- tile_height_statistics(vegetation_data$Z)

  voxel_indices <- unique(vegetation_data[, .(
    vx = floor(X / parameters$voxel_resolution),
    vy = floor(Y / parameters$voxel_resolution),
    vz = floor(Z / parameters$voxel_resolution)
  )])
  box_widths <- parameters$voxel_resolution * 2^(0:30)
  box_widths <- box_widths[box_widths < parameters$tile_size]
  scale_factors <- unique(round(box_widths / parameters$voxel_resolution))
  counts <- vapply(scale_factors, function(scale_factor) {
    uniqueN(voxel_indices[, .(
      vx = vx %/% scale_factor,
      vy = vy %/% scale_factor,
      vz = vz %/% scale_factor
    )])
  }, numeric(1))
  usable_counts <- is.finite(counts) & counts > 0
  box_dimension <- if (sum(usable_counts) >= 2) {
    unname(coef(lm(
      log(counts[usable_counts]) ~ log(1 / box_widths[usable_counts])
    ))[[2]])
  } else {
    NA_real_
  }

  vci <- vertical_complexity(vegetation_data$Z, height_reference)

  list(
    vox_filled = filled,
    vox_total = voxel_total,
    metric_height_upper_m = height_reference,
    veg_density = safe_round(filled / voxel_total),
    zsd = safe_round(z_metrics$zsd),
    zskew = safe_round(z_metrics$zskew),
    zkurt = safe_round(z_metrics$zkurt),
    zq90 = safe_round(z_metrics$zq90),
    box_dim_fixed = safe_round(box_dimension),
    vci = safe_round(vci),
    rumple = safe_round(rumple),
    gap_fraction = safe_round(gap_fraction),
    chm_sd = safe_round(chm_sd),
    chm_cv = safe_round(chm_sd / chm_mean),
    height_max = safe_round(height_max),
    height_mean = safe_round(height_mean)
  )
}

build_tile_row <- function(point_cloud_name, tile, edge_tile, tile_count,
                           metrics, tree_metrics, sensor = NA_character_,
                           country = NA_character_) {
  bounds <- st_bbox(tile)
  as.data.table(c(list(
    file = point_cloud_name,
    sensor = sensor,
    country = country,
    tile_id = tile$tile_id[[1]],
    tile_xmin = as.numeric(bounds[["xmin"]]),
    tile_ymin = as.numeric(bounds[["ymin"]]),
    edge_tile = edge_tile,
    n_tiles_plot = as.integer(tile_count)
  ), metrics, tree_metrics))
}

empty_result_table <- function() {
  data.table(
    file = character(),
    sensor = character(),
    country = character(),
    tile_id = integer(),
    tile_xmin = numeric(),
    tile_ymin = numeric(),
    edge_tile = logical(),
    n_tiles_plot = integer(),
    vox_filled = integer(),
    vox_total = numeric(),
    metric_height_upper_m = numeric(),
    veg_density = numeric(),
    zsd = numeric(),
    zskew = numeric(),
    zkurt = numeric(),
    zq90 = numeric(),
    box_dim_fixed = numeric(),
    vci = numeric(),
    rumple = numeric(),
    gap_fraction = numeric(),
    chm_sd = numeric(),
    chm_cv = numeric(),
    height_max = numeric(),
    height_mean = numeric(),
    n_seg_total = integer(),
    n_trees = integer(),
    tree_height_max = numeric(),
    tree_height_mean = numeric(),
    tree_height_gini = numeric(),
    tree_crownarea_mean = numeric(),
    tree_crownarea_max = numeric(),
    tree_crownarea_gini = numeric(),
    tree_volume_mean = numeric(),
    tree_volume_max = numeric(),
    tree_volume_gini = numeric()
  )
}

apply_output_crs <- function(geometry, source_crs) {
  if (is.null(source_crs) || is.na(source_crs) || !nzchar(source_crs)) {
    st_crs(geometry) <- NA
  } else {
    st_crs(geometry) <- st_crs(source_crs)
  }
  geometry
}

write_tile_geojson <- function(tiles, result, source_crs, path) {
  if (nrow(tiles) == 0) {
    writeLines('{"type":"FeatureCollection","features":[]}', path, useBytes = TRUE)
    return(invisible(path))
  }
  output <- apply_output_crs(tiles, source_crs)
  attributes <- as.data.frame(result)
  for (name in setdiff(names(attributes), "tile_id")) {
    output[[name]] <- attributes[[name]]
  }
  st_write(output, path, driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)
  invisible(path)
}

write_layout_png <- function(aoi, tiles, edge_flags, point_cloud_name, tile_size, path) {
  png(path, width = 1600, height = 1200, res = 150)
  on.exit(dev.off(), add = TRUE)
  par(mar = c(5, 5, 5, 2) + 0.1)
  plot_bounds <- st_bbox(aoi$inclusion)
  if (nrow(tiles) > 0) {
    tile_bounds <- st_bbox(tiles)
    plot_bounds[c("xmin", "ymin")] <- pmin(
      plot_bounds[c("xmin", "ymin")],
      tile_bounds[c("xmin", "ymin")]
    )
    plot_bounds[c("xmax", "ymax")] <- pmax(
      plot_bounds[c("xmax", "ymax")],
      tile_bounds[c("xmax", "ymax")]
    )
  }
  plot(
    aoi$inclusion,
    col = adjustcolor("#2e8b57", alpha.f = 0.12),
    border = "#17633b",
    lwd = 3,
    axes = TRUE,
    asp = 1,
    xlim = as.numeric(plot_bounds[c("xmin", "xmax")]),
    ylim = as.numeric(plot_bounds[c("ymin", "ymax")]),
    xlab = "Point-cloud-local X",
    ylab = "Point-cloud-local Y",
    main = sprintf("%s — %d valid Analysis Tiles", point_cloud_name, nrow(tiles))
  )
  if (length(aoi$exclusion) > 0) {
    plot(
      aoi$exclusion,
      add = TRUE,
      col = adjustcolor("#dc2626", alpha.f = 0.25),
      border = "#991b1b",
      lwd = 2
    )
  }
  if (nrow(tiles) > 0) {
    tile_colors <- ifelse(edge_flags, "#f59e0b", "#2563eb")
    for (index in seq_len(nrow(tiles))) {
      plot(st_geometry(tiles[index, ]), add = TRUE, border = tile_colors[[index]], lwd = 2)
    }
    centers <- st_coordinates(st_centroid(st_geometry(tiles)))
    text(centers[, 1], centers[, 2], labels = tiles$tile_id, cex = 0.85)
  } else {
    text(mean(par("usr")[1:2]), mean(par("usr")[3:4]), "0 valid tiles", cex = 1.4)
  }

  bounds <- st_bbox(aoi$inclusion)
  span_x <- as.numeric(bounds[["xmax"]] - bounds[["xmin"]])
  span_y <- as.numeric(bounds[["ymax"]] - bounds[["ymin"]])
  arrow_x <- as.numeric(bounds[["xmax"]] - span_x * 0.08)
  arrow_y <- as.numeric(bounds[["ymax"]] - span_y * 0.18)
  arrows(arrow_x, arrow_y, arrow_x, arrow_y + span_y * 0.1, length = 0.12, lwd = 2)
  text(arrow_x, arrow_y + span_y * 0.125, "N", font = 2)

  scale_x <- as.numeric(bounds[["xmin"]] + span_x * 0.06)
  scale_y <- as.numeric(bounds[["ymin"]] + span_y * 0.06)
  segments(scale_x, scale_y, scale_x + tile_size, scale_y, lwd = 4)
  segments(scale_x, scale_y - span_y * 0.01, scale_x, scale_y + span_y * 0.01, lwd = 2)
  segments(
    scale_x + tile_size,
    scale_y - span_y * 0.01,
    scale_x + tile_size,
    scale_y + span_y * 0.01,
    lwd = 2
  )
  text(scale_x + tile_size / 2, scale_y + span_y * 0.035, sprintf("%g m", tile_size))
  boundary_label <- if (identical(aoi$source, "point_cloud_extent")) {
    "Point-cloud extent"
  } else {
    "Audit AOI"
  }
  legend_labels <- c(boundary_label, "Edge tile", "Interior tile")
  legend_colors <- c("#17633b", "#f59e0b", "#2563eb")
  legend_widths <- c(3, 2, 2)
  if (length(aoi$exclusion) > 0) {
    legend_labels <- append(legend_labels, "Exclusion", after = 1)
    legend_colors <- append(legend_colors, "#991b1b", after = 1)
    legend_widths <- append(legend_widths, 2, after = 1)
  }
  legend(
    "bottomright",
    legend = legend_labels,
    col = legend_colors,
    lwd = legend_widths,
    bg = "white"
  )
  invisible(path)
}
