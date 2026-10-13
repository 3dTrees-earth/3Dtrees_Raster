suppressPackageStartupMessages({
  library(jsonlite)
  library(sf)
})

sha256_file <- function(path) {
  output <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop(sprintf("cannot hash %s: %s", path, paste(output, collapse = " ")))
  }
  sub("[[:space:]].*$", "", output[[1]])
}

polygonal_union <- function(features, label, allow_empty = FALSE) {
  if (is.null(features) || nrow(features) == 0L) {
    if (allow_empty) return(st_sfc(crs = NA_crs_))
    stop(sprintf("AOI contains no %s geometry", label))
  }
  st_crs(features) <- NA
  features <- st_zm(features, drop = TRUE, what = "ZM")
  features <- st_make_valid(features)
  geometry <- suppressWarnings(st_collection_extract(
    st_geometry(features),
    "POLYGON"
  ))
  if (length(geometry) == 0L || all(st_is_empty(geometry))) {
    if (allow_empty) return(st_sfc(crs = NA_crs_))
    stop(sprintf("AOI %s is empty or non-polygonal after validation", label))
  }
  st_union(geometry)
}

effective_geojson_geometry <- function(path) {
  features <- st_read(path, quiet = TRUE, stringsAsFactors = FALSE)
  if (nrow(features) == 0L) stop("AOI GeoJSON contains no features")
  if (!"role" %in% names(features)) {
    return(polygonal_union(features, "inclusion"))
  }

  roles <- tolower(trimws(as.character(features$role)))
  role_aliases <- c(
    include = "include",
    inclusion = "include",
    exclude = "exclude",
    exclusion = "exclude"
  )
  unknown <- unique(roles[!roles %in% names(role_aliases)])
  if (length(unknown) > 0L) {
    stop(sprintf(
      "AOI contains unsupported role(s): %s",
      paste(unknown, collapse = ", ")
    ))
  }
  roles <- unname(role_aliases[roles])
  inclusion <- polygonal_union(
    features[roles == "include", , drop = FALSE],
    "inclusion"
  )
  exclusion <- polygonal_union(
    features[roles == "exclude", , drop = FALSE],
    "exclusion",
    allow_empty = TRUE
  )
  effective <- if (length(exclusion) == 0L) {
    inclusion
  } else {
    suppressWarnings(st_difference(inclusion, exclusion))
  }
  effective <- st_make_valid(effective)
  if (length(effective) == 0L || all(st_is_empty(effective)) ||
      sum(as.numeric(st_area(effective))) <= 0) {
    stop("AOI exclusions remove the complete positive-area footprint")
  }
  st_union(effective)
}

geometry_hash <- function(geometry) {
  path <- tempfile(fileext = ".wkb")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  connection <- file(path, open = "wb")
  binaries <- st_as_binary(geometry, EWKB = FALSE)
  for (binary in binaries) writeBin(binary, connection)
  close(connection)
  sha256_file(path)
}

geometry_summary <- function(geometry) {
  bounds <- st_bbox(geometry)
  list(
    normalized_xy_sha256 = geometry_hash(geometry),
    area = unname(sum(as.numeric(st_area(geometry)))),
    bounds = unname(as.numeric(bounds[c("xmin", "ymin", "xmax", "ymax")]))
  )
}

canonical_ring_signature <- function(ring) {
  ring <- as.matrix(ring)
  if (ncol(ring) < 2L) stop("AOI polygon ring has fewer than two coordinates")
  ring <- unname(ring[, 1:2, drop = FALSE])
  if (nrow(ring) < 4L ||
      !identical(as.numeric(ring[1L, ]), as.numeric(ring[nrow(ring), ]))) {
    stop("AOI polygon contains an invalid or unclosed ring")
  }
  # The closing coordinate repeats the chosen ring start. Remove it before
  # sorting so rotating the start vertex does not change the signature.
  vertices <- ring[-nrow(ring), , drop = FALSE]
  vertices <- vertices[order(vertices[, 1L], vertices[, 2L]), , drop = FALSE]
  paste(
    apply(vertices, 1L, function(vertex) {
      paste(sprintf("%a", as.numeric(vertex)), collapse = ",")
    }),
    collapse = ";"
  )
}

canonical_polygon_coordinate_signature <- function(geometry) {
  polygons <- suppressWarnings(st_cast(geometry, "POLYGON"))
  signatures <- vapply(seq_along(polygons), function(index) {
    rings <- unclass(polygons[[index]])
    exterior <- canonical_ring_signature(rings[[1L]])
    holes <- if (length(rings) > 1L) {
      sort(vapply(rings[-1L], canonical_ring_signature, character(1)))
    } else {
      character()
    }
    paste(c(exterior, holes), collapse = "|ring|")
  }, character(1))
  sort(signatures)
}

aoi_areas_equal <- function(expected_area, actual_area) {
  is.finite(expected_area) &&
    is.finite(actual_area) &&
    abs(expected_area - actual_area) <= 1e-6
}

assert_geometry_identical <- function(expected, actual, label) {
  # GDAL may rotate a polygon ring's starting vertex while writing GeoPackage.
  # That is the same polygon, but st_equals_exact() and raw WKB are sensitive to
  # coordinate order. Require topological equality plus an exact, order-neutral
  # inventory of every non-closing XY vertex. This accepts rotation, direction,
  # ring order and polygon order, while rejecting added, removed or moved
  # vertices even when an added collinear vertex is topologically harmless.
  relation <- st_equals(expected, actual, sparse = FALSE)
  expected_summary <- geometry_summary(expected)
  actual_summary <- geometry_summary(actual)
  # sf/GDAL may rotate rings while writing. GEOS then accumulates the same
  # shoelace terms in a different order. Permit at most 1e-6 square coordinate
  # units of absolute area noise. Topology, exact vertices and exact bounds are
  # checked independently below, so this cannot admit a moved coordinate.
  area_equal <- aoi_areas_equal(
    expected_summary$area,
    actual_summary$area
  )
  if (!isTRUE(relation[[1, 1]]) ||
      !isTRUE(area_equal) ||
      !identical(expected_summary$bounds, actual_summary$bounds) ||
      !identical(
        canonical_polygon_coordinate_signature(expected),
        canonical_polygon_coordinate_signature(actual)
      )) {
    stop(sprintf("%s changed normalized geometry, area, or bounds", label))
  }
  invisible(TRUE)
}

assert_point_cloud_overlap <- function(geometry, point_bounds) {
  point_extent <- st_sfc(st_polygon(list(matrix(c(
    point_bounds[["xmin"]], point_bounds[["ymin"]],
    point_bounds[["xmax"]], point_bounds[["ymin"]],
    point_bounds[["xmax"]], point_bounds[["ymax"]],
    point_bounds[["xmin"]], point_bounds[["ymax"]],
    point_bounds[["xmin"]], point_bounds[["ymin"]]
  ), ncol = 2, byrow = TRUE))), crs = NA_crs_)
  overlap <- suppressWarnings(st_intersection(geometry, point_extent))
  area <- if (length(overlap) == 0L) 0 else sum(as.numeric(st_area(overlap)))
  if (!is.finite(area) || area <= 0) {
    stop("converted AOI does not overlap the point-cloud-local XY bounds")
  }
  invisible(area)
}

tile_geometry_hash <- function(tiles) {
  if (nrow(tiles) == 0L) return(sha256_file("/dev/null"))
  geometry_hash(st_geometry(tiles))
}

convert_geojson_to_gpkg <- function(source_geojson, target_gpkg, point_bounds,
                                    tile_size = 20, search_step = 0.5,
                                    tile_threads = 1L) {
  if (!file.exists(source_geojson) ||
      !grepl("\\.(geojson|json)$", source_geojson, ignore.case = TRUE)) {
    stop("--aoi-geojson must be one existing GeoJSON file")
  }
  if (!dir.exists(dirname(target_gpkg)) ||
      file.access(dirname(target_gpkg), mode = 2) != 0) {
    stop("generated GeoPackage directory must exist and be writable")
  }

  source_geometry <- effective_geojson_geometry(source_geojson)
  assert_point_cloud_overlap(source_geometry, point_bounds)
  source_tiles <- build_optimized_tiles(
    source_geometry,
    tile_size,
    search_step,
    tile_threads
  )
  temporary_gpkg <- tempfile(
    pattern = paste0(".", tools::file_path_sans_ext(basename(target_gpkg)), "_"),
    tmpdir = dirname(target_gpkg),
    fileext = ".gpkg"
  )
  on.exit(unlink(temporary_gpkg, force = TRUE), add = TRUE)
  st_write(
    st_sf(role = "include", geometry = source_geometry),
    temporary_gpkg,
    layer = "include",
    driver = "GPKG",
    delete_dsn = TRUE,
    quiet = TRUE
  )
  converted <- st_read(
    temporary_gpkg,
    layer = "include",
    quiet = TRUE,
    stringsAsFactors = FALSE
  )
  st_crs(converted) <- NA
  converted_geometry <- polygonal_union(converted, "converted inclusion")
  assert_geometry_identical(
    source_geometry,
    converted_geometry,
    "GeoJSON to GeoPackage conversion"
  )
  converted_tiles <- build_optimized_tiles(
    converted_geometry,
    tile_size,
    search_step,
    tile_threads
  )
  if (!identical(source_tiles$tile_id, converted_tiles$tile_id) ||
      !identical(tile_geometry_hash(source_tiles), tile_geometry_hash(converted_tiles))) {
    stop("GeoJSON to GeoPackage conversion changed the optimized tile grid")
  }
  unlink(target_gpkg, force = TRUE)
  if (!file.rename(temporary_gpkg, target_gpkg)) {
    stop("cannot atomically promote generated GeoPackage")
  }

  source_summary <- geometry_summary(source_geometry)
  list(
    source_geojson = normalizePath(source_geojson, mustWork = TRUE),
    generated_gpkg = normalizePath(target_gpkg, mustWork = TRUE),
    source_geojson_sha256 = sha256_file(source_geojson),
    generated_gpkg_sha256 = sha256_file(target_gpkg),
    normalized_xy_sha256 = source_summary$normalized_xy_sha256,
    area = source_summary$area,
    bounds = source_summary$bounds,
    tile_count = nrow(source_tiles),
    tile_geometry_sha256 = tile_geometry_hash(source_tiles)
  )
}
