local({
  work_directory <- tempfile("foreststructure_vrt_test_")
  dir.create(work_directory)
  on.exit(unlink(work_directory, recursive = TRUE, force = TRUE), add = TRUE)

  write_test_raster <- function(name, values, crs = "", datatype = NULL) {
    path <- file.path(work_directory, paste0(name, ".tif"))
    raster <- terra::rast(xmin = 0, xmax = 2, ymin = 0, ymax = 1, res = 1)
    terra::values(raster) <- values
    if (nzchar(crs)) terra::crs(raster) <- crs
    if (is.null(datatype)) {
      terra::writeRaster(raster, path, overwrite = TRUE)
    } else {
      terra::writeRaster(
        raster,
        path,
        overwrite = TRUE,
        datatype = datatype
      )
    }
    path
  }

  gdal_raster_has_crs <- function(path) {
    info <- system2(
      "gdalinfo",
      c("-json", shQuote(path)),
      stdout = TRUE,
      stderr = TRUE
    )
    status <- attr(info, "status")
    if (!is.null(status) && status != 0L) {
      stop(sprintf("gdalinfo failed for %s", path))
    }
    any(grepl('"coordinateSystem"', info, fixed = TRUE))
  }

  first_path <- write_test_raster("first", c(1, 1))
  second_path <- write_test_raster("second", c(3, 3))
  empty_no_crs_path <- write_test_raster(
    "empty_no_crs",
    c(NA_real_, NA_real_)
  )
  empty_path <- write_test_raster(
    "empty",
    c(NA_real_, NA_real_),
    "EPSG:4326"
  )
  empty_raster <- terra::rast(empty_path)

  result <- build_chunk_virtual_raster(
    list(first_path, second_path, empty_path),
    work_directory,
    "TEST"
  )

  # The CHM path preserves lidR 4.3.2's disk-backed GDAL VRT semantics. Later
  # populated chunk rasters win in overlaps; CHM values must never be averaged.
  stopifnot(!terra::inMemory(result))
  stopifnot(identical(as.numeric(terra::values(result)), c(3, 3)))
  stopifnot(!any(as.numeric(terra::values(result)) == 2))
  stopifnot(terra::same.crs(result, empty_raster))

  empty_result <- build_chunk_virtual_raster(
    list(empty_no_crs_path, empty_path),
    work_directory,
    "EMPTY"
  )
  stopifnot(!terra::inMemory(empty_result))
  stopifnot(all(is.na(as.numeric(terra::values(empty_result)))))

  fallback_result <- build_chunk_virtual_raster(
    list(first_path, second_path),
    work_directory,
    "FALLBACK",
    terra::crs(empty_raster)
  )
  stopifnot(!terra::inMemory(fallback_result))
  stopifnot(terra::same.crs(fallback_result, empty_raster))

  # A sparse first DTM chunk can be all-NoData and therefore CRS-less, while a
  # later populated chunk carries the source CRS. The overlap mosaic must scan
  # every chunk instead of inheriting only the first chunk's missing metadata.
  crsless_dtm_path <- file.path(work_directory, "a_empty_dtm.tif")
  crsless_dtm <- terra::rast(
    xmin = 1886455,
    xmax = 1886457,
    ymin = 5771842,
    ymax = 5771843,
    res = 1
  )
  terra::values(crsless_dtm) <- c(NA_real_, NA_real_)
  terra::crs(crsless_dtm) <- ""
  terra::writeRaster(crsless_dtm, crsless_dtm_path, overwrite = TRUE)
  stopifnot(!nzchar(terra::crs(terra::rast(crsless_dtm_path))))
  populated_crs_path <- file.path(work_directory, "b_populated_dtm.tif")
  populated_crs <- crsless_dtm
  terra::values(populated_crs) <- c(4, 5)
  terra::crs(populated_crs) <- "EPSG:2193"
  terra::writeRaster(populated_crs, populated_crs_path, overwrite = TRUE)
  dtm_work_directory <- file.path(work_directory, "dtm")
  dir.create(dtm_work_directory)
  dtm_result <- build_dtm_overlap_mosaic(
    c(crsless_dtm_path, populated_crs_path),
    dtm_work_directory
  )
  stopifnot(identical(as.numeric(terra::values(dtm_result)), c(4, 5)))
  stopifnot(terra::same.crs(dtm_result, terra::rast(populated_crs_path)))

  fallback_dtm_work_directory <- file.path(work_directory, "fallback_dtm")
  dir.create(fallback_dtm_work_directory)
  fallback_dtm_result <- build_dtm_overlap_mosaic(
    crsless_dtm_path,
    fallback_dtm_work_directory,
    terra::crs(populated_crs)
  )
  stopifnot(all(is.na(as.numeric(terra::values(fallback_dtm_result)))))
  stopifnot(
    terra::same.crs(fallback_dtm_result, terra::rast(populated_crs_path))
  )

  # LAScatalog returns chunk rasters in numeric catalog order. Keep that order:
  # overlapping FLT8S means are not associative at floating-point boundaries,
  # and lexically sorting dtm_1...dtm_16 changes the scientific result.
  order_sensitive_values <- rep(NA_real_, 16L)
  order_sensitive_values[c(5L, 6L, 9L, 10L)] <- c(
    5.4698000000000002,
    5.319,
    5.319,
    5.3012000000000006
  )
  order_sensitive_paths <- vapply(seq_len(16L), function(index) {
    write_test_raster(
      paste0("dtm_", index),
      rep(order_sensitive_values[[index]], 2L),
      datatype = "FLT8S"
    )
  }, character(1))
  catalog_results <- as.list(order_sensitive_paths)
  numeric_paths <- complete_chunk_raster_paths(catalog_results, "DTM-ORDER")
  stopifnot(identical(numeric_paths, unname(order_sensitive_paths)))

  numeric_order_directory <- file.path(work_directory, "dtm_numeric_order")
  lexical_order_directory <- file.path(work_directory, "dtm_lexical_order")
  dir.create(numeric_order_directory)
  dir.create(lexical_order_directory)
  numeric_order_result <- build_dtm_overlap_mosaic(
    numeric_paths,
    numeric_order_directory
  )
  lexical_order_result <- build_dtm_overlap_mosaic(
    sort(numeric_paths),
    lexical_order_directory
  )
  numeric_order_value <- as.numeric(terra::values(numeric_order_result))[[1L]]
  lexical_order_value <- as.numeric(terra::values(lexical_order_result))[[1L]]
  stopifnot(identical(numeric_order_value, 5.3522500000000006))
  stopifnot(!identical(lexical_order_value, numeric_order_value))

  dtm_writer <- paste(deparse(body(write_global_dtm)), collapse = "\n")
  stopifnot(grepl(
    'complete_chunk_raster_paths(results, "DTM")',
    dtm_writer,
    fixed = TRUE
  ))
  stopifnot(!grepl("sort(", dtm_writer, fixed = TRUE))
  stopifnot(grepl('extent_mode = "lidr"', dtm_writer, fixed = TRUE))

  point_cloud_path <- file.path(work_directory, "fallback_bounds.las")
  point_cloud <- LAS(data.frame(
    X = c(0.25, 1.75),
    Y = c(0.25, 0.75),
    Z = c(1, 2)
  ))
  writeLAS(point_cloud, point_cloud_path)
  covering_fallback <- cover_virtual_raster(
    fallback_result,
    point_cloud_path,
    1,
    work_directory,
    "FALLBACK"
  )
  streamed_fallback_path <- file.path(work_directory, "fallback_streamed.tif")
  stream_virtual_raster(
    covering_fallback,
    streamed_fallback_path,
    "FALLBACK"
  )
  streamed_fallback <- terra::rast(streamed_fallback_path)
  stopifnot(terra::same.crs(streamed_fallback, empty_raster))

  # The original script delegates DTM geometry to lidR's raster layout. When a LAS bound is
  # exactly grid-aligned, lidR 4.3.2 retains an additional border cell instead
  # of applying ordinary floor/ceiling point-bound coverage.
  exact_boundary_cloud_path <- file.path(
    work_directory,
    "exact_boundary.las"
  )
  exact_boundary_cloud <- LAS(data.frame(
    X = c(0, 2, 2, 0),
    Y = c(0, 0, 2, 2),
    Z = c(1, 1.1, 1.2, 1.3)
  ))
  writeLAS(exact_boundary_cloud, exact_boundary_cloud_path)
  lidr_layout <- lidR:::raster_layout(
    readLAScatalog(exact_boundary_cloud_path),
    1,
    start = c(0, 0),
    buffer = 0,
    format = "template"
  )
  exact_boundary_source_path <- file.path(
    work_directory,
    "exact_boundary_source.tif"
  )
  exact_boundary_source <- terra::rast(
    xmin = lidr_layout$xmin,
    xmax = lidr_layout$xmax,
    ymin = lidr_layout$ymin,
    ymax = lidr_layout$ymax,
    resolution = 1
  )
  terra::values(exact_boundary_source) <- seq_len(
    terra::ncell(exact_boundary_source)
  )
  terra::writeRaster(
    exact_boundary_source,
    exact_boundary_source_path,
    overwrite = TRUE
  )
  exact_boundary_result <- build_chunk_virtual_raster(
    list(exact_boundary_source_path),
    work_directory,
    "LIDR-LAYOUT-SOURCE"
  )
  lidr_covering <- cover_virtual_raster(
    exact_boundary_result,
    exact_boundary_cloud_path,
    1,
    work_directory,
    "LIDR-LAYOUT",
    extent_mode = "lidr"
  )
  stopifnot(terra::compareGeom(
    lidr_covering,
    exact_boundary_source,
    crs = FALSE,
    ext = TRUE,
    rowcol = TRUE,
    res = TRUE,
    stopOnError = FALSE
  ))

  # Lock down the two production grid shapes that exposed the publication
  # regression. These are header-layout tests only: they do not need the large
  # point clouds or run PTD/TIN, so the failure signal stays fast and exact.
  write_layout_cloud <- function(name, xmax, ymax) {
    path <- file.path(work_directory, paste0(name, ".las"))
    cloud <- LAS(data.frame(
      X = c(0, xmax),
      Y = c(0, ymax),
      Z = c(1, 2)
    ))
    writeLAS(cloud, path)
    path
  }

  layout_970_path <- write_layout_cloud("layout_970", 238.99, 239.99)
  stopifnot(identical(
    unname(lidr_catalog_grid_bounds(layout_970_path, 1)),
    c(-1, -1, 239, 240)
  ))

  layout_2062_path <- write_layout_cloud("layout_2062", 306.99, 197.99)
  stopifnot(identical(
    unname(lidr_catalog_grid_bounds(layout_2062_path, 1)),
    c(-1, -1, 307, 198)
  ))

  # A non-aligned cloud must not receive an unconditional border. Its pinned
  # lidR layout is the ordinary covering grid.
  unaligned_layout_path <- file.path(work_directory, "layout_unaligned.las")
  writeLAS(LAS(data.frame(
    X = c(0.25, 2.75),
    Y = c(0.25, 3.75),
    Z = c(1, 2)
  )), unaligned_layout_path)
  stopifnot(identical(
    lidr_catalog_grid_bounds(unaligned_layout_path, 1),
    covering_grid_bounds(
      point_cloud_xy_bounds(unaligned_layout_path),
      1
    )
  ))

  # The additional lidR border is empty. Publication must preserve the inner
  # cells exactly and retain NA in the outer cells; it must not interpolate or
  # recompute any value.
  border_source_path <- file.path(work_directory, "layout_border_source.tif")
  border_source <- terra::rast(
    xmin = lidr_layout$xmin,
    xmax = lidr_layout$xmax,
    ymin = lidr_layout$ymin,
    ymax = lidr_layout$ymax,
    resolution = 1
  )
  border_xy <- terra::xyFromCell(border_source, seq_len(terra::ncell(border_source)))
  border_values <- rep(NA_real_, terra::ncell(border_source))
  inner <- border_xy[, "x"] >= 0 & border_xy[, "y"] >= 0
  border_values[inner] <- seq_len(sum(inner))
  terra::values(border_source) <- border_values
  terra::writeRaster(border_source, border_source_path, overwrite = TRUE)
  border_result <- build_chunk_virtual_raster(
    list(border_source_path),
    work_directory,
    "LIDR-EMPTY-BORDER"
  )
  border_covering <- cover_virtual_raster(
    border_result,
    exact_boundary_cloud_path,
    1,
    work_directory,
    "LIDR-EMPTY-BORDER",
    extent_mode = "lidr"
  )
  stopifnot(terra::compareGeom(
    border_covering,
    border_source,
    crs = FALSE,
    ext = TRUE,
    rowcol = TRUE,
    res = TRUE,
    stopOnError = FALSE
  ))
  border_actual <- as.numeric(terra::values(border_covering))
  border_expected <- as.numeric(terra::values(border_source))
  stopifnot(identical(is.na(border_actual), is.na(border_expected)))
  stopifnot(identical(
    border_actual[!is.na(border_actual)],
    border_expected[!is.na(border_expected)]
  ))

  # terra/GDAL can infer WGS84 when a CRS-less local raster happens to fit
  # longitude/latitude bounds. The LAS header is authoritative: explicitly
  # publishing a missing source CRS must remove that inferred metadata.
  local_crsless_path <- write_test_raster(
    "local_crsless",
    c(7, 8)
  )
  stopifnot(gdal_raster_has_crs(local_crsless_path))
  local_result <- build_chunk_virtual_raster(
    list(local_crsless_path),
    work_directory,
    "LOCAL-CRSLESS"
  )
  local_streamed_path <- file.path(
    work_directory,
    "local_crsless_streamed.tif"
  )
  stream_virtual_raster(
    local_result,
    local_streamed_path,
    "LOCAL-CRSLESS",
    source_crs = NA_character_
  )
  stopifnot(!gdal_raster_has_crs(local_streamed_path))

  geometry <- sf::st_sfc(sf::st_point(c(1, 1)), crs = 4326)
  stopifnot(is.na(sf::st_crs(apply_output_crs(
    geometry,
    NA_character_
  ))))
  known_crs <- sf::st_crs(32632)
  stopifnot(sf::st_crs(apply_output_crs(
    geometry,
    known_crs$wkt
  )) == known_crs)

  tile_polygon <- sf::st_polygon(list(rbind(
    c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0)
  )))
  tile_output_path <- file.path(work_directory, "crsless_tiles.geojson")
  write_tile_geojson(
    sf::st_sf(tile_id = 1L, geometry = sf::st_sfc(tile_polygon)),
    data.frame(tile_id = 1L),
    NA_character_,
    tile_output_path
  )
  tile_json <- paste(readLines(tile_output_path, warn = FALSE), collapse = "")
  stopifnot(!grepl('"crs"', tile_json, fixed = TRUE))

  conflicting_path <- write_test_raster(
    "conflicting",
    c(NA_real_, NA_real_),
    "EPSG:3857"
  )
  conflicting <- tryCatch(
    {
      build_chunk_virtual_raster(
        list(empty_path, conflicting_path),
        work_directory,
        "CONFLICTING"
      )
      FALSE
    },
    error = function(error) grepl(
      "conflicting CRS metadata",
      conditionMessage(error),
      fixed = TRUE
    )
  )
  stopifnot(conflicting)
})
