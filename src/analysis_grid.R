# The audit footprint and the deterministic tile grid. Sourced by analysis.R;
# uses its constants and helpers.

# Audit footprint and deterministic tile grid ------------------------------

normalize_aoi_geometry <- function(features, label, allow_empty = FALSE) {
  if (is.null(features) || nrow(features) == 0) {
    if (allow_empty) return(st_sfc(crs = NA_crs_))
    stop(sprintf("Audit AOI contains no %s geometry", label))
  }
  geometry_types <- unique(as.character(st_geometry_type(features, by_geometry = TRUE)))
  if (!all(geometry_types %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop(sprintf("Audit AOI %s must contain only Polygon or MultiPolygon geometry", label))
  }

  # Audit AOIs are already expressed in the point cloud's local XY space.
  # GeoJSON readers may infer WGS84 from the container format; that metadata is
  # not authoritative for the local audit handoff and must not trigger a transform.
  st_crs(features) <- NA
  features <- st_make_valid(features)
  geometry <- st_union(st_geometry(features))
  if (length(geometry) == 0 || any(st_is_empty(geometry))) {
    if (allow_empty) return(st_sfc(crs = NA_crs_))
    stop(sprintf("Audit AOI %s is empty after geometry validation", label))
  }
  geometry
}

read_geojson_aoi <- function(path) {
  features <- st_read(path, quiet = TRUE, stringsAsFactors = FALSE)
  if (nrow(features) == 0) stop("Audit AOI contains no geometry")
  if (!"role" %in% names(features)) {
    return(list(inclusion = features, exclusion = NULL))
  }

  roles <- tolower(trimws(as.character(features$role)))
  aliases <- c(
    include = "include",
    inclusion = "include",
    exclude = "exclude",
    exclusion = "exclude"
  )
  unknown <- unique(roles[!roles %in% names(aliases)])
  if (length(unknown) > 0) {
    stop(sprintf("Audit AOI contains unsupported feature role(s): %s", paste(unknown, collapse = ", ")))
  }
  roles <- unname(aliases[roles])
  list(
    inclusion = features[roles == "include", , drop = FALSE],
    exclusion = features[roles == "exclude", , drop = FALSE]
  )
}

read_geopackage_aoi <- function(path) {
  layers <- st_layers(path)$name
  if (!"include" %in% layers) {
    stop("GeoPackage Audit AOI must contain an 'include' layer")
  }
  list(
    inclusion = st_read(path, layer = "include", quiet = TRUE, stringsAsFactors = FALSE),
    exclusion = if ("exclude" %in% layers) {
      st_read(path, layer = "exclude", quiet = TRUE, stringsAsFactors = FALSE)
    } else {
      NULL
    }
  )
}

read_audit_aoi <- function(path) {
  raw <- if (grepl("\\.gpkg$", path, ignore.case = TRUE)) {
    read_geopackage_aoi(path)
  } else {
    read_geojson_aoi(path)
  }
  inclusion <- normalize_aoi_geometry(raw$inclusion, "inclusion")
  exclusion <- normalize_aoi_geometry(raw$exclusion, "exclusion", allow_empty = TRUE)
  usable <- if (length(exclusion) > 0) {
    st_difference(inclusion, exclusion)
  } else {
    inclusion
  }
  if (length(usable) == 0 || all(st_is_empty(usable))) {
    stop("Audit AOI exclusions remove the complete inclusion geometry")
  }
  list(
    inclusion = inclusion,
    exclusion = exclusion,
    usable = usable,
    source = "audit_aoi"
  )
}

# Snap header bounds to their stored coordinate lattice. Some writers retain
# pre-quantization extrema, while COPC reports the actual stored coordinates.
# This is only an inexpensive header compatibility check, not identity proof.
point_cloud_xyz_bounds <- function(path) {
  # lidR exposes LAS/COPC public-header fields through the same header object.
  header <- readLASheader(path)
  # Return a named vector so comparison order is explicit and reviewable.
  bounds <- c(
    xmin = as.numeric(header@PHB[["Min X"]]),
    ymin = as.numeric(header@PHB[["Min Y"]]),
    zmin = as.numeric(header@PHB[["Min Z"]]),
    xmax = as.numeric(header@PHB[["Max X"]]),
    ymax = as.numeric(header@PHB[["Max Y"]]),
    zmax = as.numeric(header@PHB[["Max Z"]])
  )
  scales <- rep(vapply(c("X", "Y", "Z"), function(axis) {
    as.numeric(header@PHB[[paste(axis, "scale factor")]])
  }, numeric(1)), 2L)
  offsets <- rep(vapply(c("X", "Y", "Z"), function(axis) {
    as.numeric(header@PHB[[paste(axis, "offset")]])
  }, numeric(1)), 2L)
  if (any(!is.finite(c(bounds, scales, offsets))) || any(scales <= 0)) {
    stop("invalid LAS coordinate lattice in companion header")
  }
  round((bounds - offsets) / scales) * scales + offsets
}

point_cloud_xy_bounds <- function(path) {
  header <- readLASheader(path)
  bounds <- c(
    xmin = as.numeric(header@PHB[["Min X"]]),
    ymin = as.numeric(header@PHB[["Min Y"]]),
    xmax = as.numeric(header@PHB[["Max X"]]),
    ymax = as.numeric(header@PHB[["Max Y"]])
  )
  if (any(!is.finite(bounds)) ||
      bounds[["xmax"]] <= bounds[["xmin"]] ||
      bounds[["ymax"]] <= bounds[["ymin"]]) {
    stop("point-cloud header does not contain a usable two-dimensional XY extent")
  }
  bounds
}

# The pinned catalog index requires positive-area overlap. A source point at
# the maximum XY exactly on a chunk boundary otherwise has no scheduled owner.
# Enlarge only in-memory file lookup bounds; source coordinates, scientific
# chunk origins and the half-open read filters are unchanged.
catalog_with_source_boundaries <- function(catalog, point_cloud) {
  bounds <- point_cloud_xyz_bounds(point_cloud)[c("xmin", "ymin", "xmax", "ymax")]
  epsilon <- pmax(1e-9, abs(bounds[3:4]) * .Machine$double.eps * 16)
  catalog@data$Min.X <- bounds[[1L]]
  catalog@data$Min.Y <- bounds[[2L]]
  catalog@data$Max.X <- bounds[[3L]] + epsilon[[1L]]
  catalog@data$Max.Y <- bounds[[4L]] + epsilon[[2L]]
  catalog
}

# LASlib checks a file's raw header bounds before opening its spatial index.
# When a rounded-down header maximum precedes the last stored coordinate, an
# exact terminal query is skipped. Widen that query within the final lattice
# interval only, so no preceding source coordinate can be read twice.
catalog_boundary_read_adjustment <- function(point_cloud) {
  raw <- point_cloud_xy_bounds(point_cloud)
  stored <- point_cloud_xyz_bounds(point_cloud)
  header <- readLASheader(point_cloud)
  axes <- lapply(c("x", "y"), function(axis) {
    maximum <- paste0(axis, "max")
    if (raw[[maximum]] >= stored[[maximum]]) return(NULL)
    scale <- as.numeric(header@PHB[[paste(toupper(axis), "scale factor")]])
    previous <- stored[[maximum]] - scale
    lower <- previous + (raw[[maximum]] - previous) / 2
    if (!is.finite(lower) || !(lower > previous && lower < raw[[maximum]])) {
      stop("unsupported source lattice: cannot safely query rounded header maximum")
    }
    list(raw_maximum = raw[[maximum]], stored_maximum = stored[[maximum]], lower = lower)
  })
  if (all(vapply(axes, is.null, logical(1)))) NULL else axes
}

adjust_catalog_boundary_read <- function(chunk, adjustment = NULL) {
  if (is.null(adjustment)) return(chunk)
  fields <- strsplit(trimws(chunk@filter), "[[:space:]]+")[[1L]]
  start <- match("-inside", fields)
  if (is.na(start) || start + 4L > length(fields)) stop("boundary adjustment requires a rectangular catalog read")
  bounds <- suppressWarnings(as.numeric(fields[start + seq_len(4L)]))
  if (any(!is.finite(bounds))) stop("invalid catalog boundary query")
  for (axis in seq_len(2L)) {
    value <- adjustment[[axis]]
    if (!is.null(value) && bounds[[axis]] >= value$raw_maximum &&
        bounds[[axis]] <= value$stored_maximum) {
      fields[[start + axis]] <- format(value$lower, scientific = FALSE, trim = TRUE, digits = 17)
    }
  }
  chunk@filter <- paste(fields, collapse = " ")
  chunk
}

read_point_cloud_extent <- function(path) {
  # Use stored-coordinate extrema for the fallback AOI too: a pre-quantization
  # header must not invalidate a complete boundary voxel.
  bounds <- point_cloud_xyz_bounds(path)[c("xmin", "ymin", "xmax", "ymax")]
  ring <- matrix(c(
    bounds[["xmin"]], bounds[["ymin"]],
    bounds[["xmax"]], bounds[["ymin"]],
    bounds[["xmax"]], bounds[["ymax"]],
    bounds[["xmin"]], bounds[["ymax"]],
    bounds[["xmin"]], bounds[["ymin"]]
  ), ncol = 2, byrow = TRUE)
  extent <- st_sfc(st_polygon(list(ring)), crs = NA_crs_)
  list(
    inclusion = extent,
    exclusion = st_sfc(crs = NA_crs_),
    usable = extent,
    source = "point_cloud_extent"
  )
}

read_analysis_footprint <- function(point_cloud, aoi = NULL) {
  if (is.null(aoi)) return(read_point_cloud_extent(point_cloud))
  audit <- read_audit_aoi(aoi)
  point_cloud_extent <- read_point_cloud_extent(point_cloud)
  overlap <- suppressWarnings(st_intersection(
    audit$usable,
    point_cloud_extent$usable
  ))
  overlap_area <- if (length(overlap) == 0L) 0 else sum(as.numeric(st_area(overlap)))
  if (!is.finite(overlap_area) || overlap_area <= 0) {
    stop(
      "Audit AOI does not overlap the point-cloud XY extent; ",
      "both inputs must use the same local coordinate frame"
    )
  }
  audit
}

make_grid_at <- function(geometry, tile_size, offset_x, offset_y,
                         reference_bounds = st_bbox(geometry)) {
  st_make_grid(
    geometry,
    cellsize = tile_size,
    offset = c(
      reference_bounds[["xmin"]] - offset_x,
      reference_bounds[["ymin"]] - offset_y
    )
  )
}

count_complete_tiles <- function(grid, geometry) {
  sum(lengths(st_within(grid, geometry)) > 0)
}

build_optimized_tiles <- function(geometry, tile_size, search_step, threads = 1L) {
  offsets <- seq(0, tile_size - search_step, by = search_step)
  reference_bounds <- st_bbox(geometry)
  components <- st_cast(geometry, "POLYGON", warn = FALSE)
  candidates <- data.frame(
    offset_x = rep(offsets, each = length(offsets)),
    offset_y = rep(offsets, times = length(offsets))
  )
  maximum_possible_tiles <- floor(
    as.numeric(st_area(geometry)) / tile_size^2
  )
  count_candidate <- function(index) {
    sum(vapply(seq_along(components), function(component_index) {
      component <- components[component_index]
      count_complete_tiles(
        make_grid_at(
          component,
          tile_size,
          candidates$offset_x[[index]],
          candidates$offset_y[[index]],
          reference_bounds
        ),
        component
      )
    }, integer(1)))
  }

  first_count <- count_candidate(1L)
  if (first_count >= maximum_possible_tiles || nrow(candidates) == 1L) {
    best_index <- 1L
  } else {
    remaining <- seq.int(2L, nrow(candidates))
    workers <- min(max(1L, as.integer(threads)), length(remaining))
    message(sprintf(
      "Optimizing %d remaining grid offsets with %d worker(s)",
      length(remaining),
      workers
    ))
    remaining_counts <- if (workers == 1L) {
      vapply(remaining, count_candidate, integer(1))
    } else {
      unlist(
        parallel::mclapply(
          remaining,
          count_candidate,
          mc.cores = workers,
          mc.preschedule = TRUE,
          mc.set.seed = FALSE
        ),
        use.names = FALSE
      )
    }
    counts <- c(first_count, remaining_counts)
    best_index <- which.max(counts)
  }

  best_x <- candidates$offset_x[[best_index]]
  best_y <- candidates$offset_y[[best_index]]
  grid <- make_grid_at(geometry, tile_size, best_x, best_y)
  complete <- lengths(st_within(grid, geometry)) > 0
  if (!any(complete)) {
    return(st_sf(
      tile_id = integer(),
      geometry = st_sfc(crs = st_crs(geometry))
    ))
  }

  tiles <- st_sf(geometry = grid[complete])
  bounds <- t(vapply(seq_len(nrow(tiles)), function(index) {
    as.numeric(st_bbox(tiles[index, ]))
  }, numeric(4)))
  order_index <- order(bounds[, 2], bounds[, 1])
  tiles <- tiles[order_index, , drop = FALSE]
  tiles$tile_id <- seq_len(nrow(tiles))
  tiles[, c("tile_id", "geometry")]
}

build_covering_tiles <- function(geometry, tile_size) {
  bounds <- st_bbox(geometry)
  span_x <- as.numeric(bounds[["xmax"]] - bounds[["xmin"]])
  span_y <- as.numeric(bounds[["ymax"]] - bounds[["ymin"]])
  columns <- as.integer(ceiling(span_x / tile_size))
  rows <- as.integer(ceiling(span_y / tile_size))
  covered_width <- columns * tile_size
  covered_height <- rows * tile_size
  origin_x <- as.numeric(bounds[["xmin"]] - (covered_width - span_x) / 2)
  origin_y <- as.numeric(bounds[["ymin"]] - (covered_height - span_y) / 2)

  polygons <- vector("list", columns * rows)
  index <- 1L
  for (row in seq_len(rows)) {
    ymin <- origin_y + (row - 1) * tile_size
    ymax <- ymin + tile_size
    for (column in seq_len(columns)) {
      xmin <- origin_x + (column - 1) * tile_size
      xmax <- xmin + tile_size
      ring <- matrix(c(
        xmin, ymin,
        xmax, ymin,
        xmax, ymax,
        xmin, ymax,
        xmin, ymin
      ), ncol = 2, byrow = TRUE)
      polygons[[index]] <- st_polygon(list(ring))
      index <- index + 1L
    }
  }
  st_sf(
    tile_id = seq_along(polygons),
    geometry = st_sfc(polygons, crs = st_crs(geometry))
  )
}

build_analysis_tiles <- function(footprint, tile_size, search_step, threads = 1L) {
  if (identical(footprint$source, "point_cloud_extent")) {
    return(build_covering_tiles(footprint$usable, tile_size))
  }
  build_optimized_tiles(footprint$usable, tile_size, search_step, threads)
}

compute_edge_flags <- function(tiles, tile_size) {
  if (nrow(tiles) == 0) return(logical())
  centroids <- st_coordinates(st_centroid(st_geometry(tiles)))
  index_x <- as.integer(round(centroids[, 1] / tile_size))
  index_y <- as.integer(round(centroids[, 2] / tile_size))
  keys <- paste(index_x, index_y, sep = "_")

  vapply(seq_len(nrow(tiles)), function(index) {
    neighbours <- expand.grid(dx = -1:1, dy = -1:1)
    neighbours <- neighbours[neighbours$dx != 0 | neighbours$dy != 0, ]
    neighbour_keys <- paste(
      index_x[[index]] + neighbours$dx,
      index_y[[index]] + neighbours$dy,
      sep = "_"
    )
    !all(neighbour_keys %in% keys)
  }, logical(1))
}
