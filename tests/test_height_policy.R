# Exact Julia statistics for ordinary heights, including one-metre boundaries.
local({
  legacy_vci <- function(z, height) {
    layers <- cut(z, seq(0, height, 1), include.lowest = TRUE)
    counts <- table(factor(layers, levels = levels(layers)))
    p <- counts / sum(counts)
    -sum(p[p > 0] * log(p[p > 0])) / log(length(counts))
  }
  for (z in list(c(0.5, 1, 1.01, 2, 2, 69.9, 70), rep(5, 7),
                 c(0.7, 1.1, 12.3, 23.4, 42, 69.7), 4)) {
    actual <- tile_height_statistics(z)
    expected <- lidR::stdmetrics_z(z)
    for (name in names(actual)) {
      stopifnot(isTRUE(all.equal(actual[[name]], expected[[name]], tolerance = 0)))
    }
    stopifnot(isTRUE(all.equal(vertical_complexity(z, 70), legacy_vci(z, 70), tolerance = 1e-15)))
  }
  stopifnot(is.na(vertical_complexity(c(0.1, 0.8), 1)))
  # A huge supported range must not allocate a vector for every empty metre.
  z <- c(0.7, 1, 80, 80, 1e12)
  expected <- -sum(c(2, 2, 1) / 5 * log(c(2, 2, 1) / 5)) / log(1e12)
  stopifnot(isTRUE(all.equal(vertical_complexity(z, 1e12), expected, tolerance = 1e-15)),
    is.finite(tile_height_statistics(z)$zq90))
  invalid <- tryCatch(vertical_complexity(1, 2^53), error = identity)
  stopifnot(inherits(invalid, "error"))

  root <- tempfile("height_reference_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  stopifnot(metric_height_reference(80.1) == 81,
    metric_height_reference(80.1, 70) == 70,
    metric_height_reference(0.1) == 1)
  invalid <- tryCatch(metric_height_reference(NA_real_), error = identity)
  stopifnot(inherits(invalid, "error"))

  # lidR's exact vertical grid can put an integer-valued source height just
  # above that integer. lidR rounds the CHM to three decimals, so deriving the
  # reference from the raster would silently omit the final VCI layer.
  cloud_path <- file.path(root, "precision.las")
  cloud <- LAS(data.frame(X = c(1, 2, 3), Y = c(1, 2, 3), Z = c(1, 25, 2)))
  sf::st_crs(cloud) <- 32632
  writeLAS(cloud, cloud_path)
  heights <- voxelize_points(cloud, 0.2)@data$Z
  stopifnot(max(heights) > 25)
  dtm <- terra::rast(xmin = 0, xmax = 5, ymin = 0, ymax = 5, res = 1, crs = "EPSG:32632")
  terra::values(dtm) <- 0
  dtm_path <- file.path(root, "dtm.tif")
  terra::writeRaster(dtm, dtm_path, datatype = "FLT8S")
  path <- file.path(root, "precision_chm.tif")
  observations <- write_global_chm(cloud_path, dtm_path, path, file.path(root, "chunks"),
    10, 0.2, Inf, 0.5, 1L)
  stopifnot(as.numeric(terra::global(terra::rast(path), "max", na.rm = TRUE)[1, 1]) == 25,
    observations$maximum_voxel_height == max(heights),
    metric_height_reference(observations$maximum_voxel_height) == 26,
    all(heights <= metric_height_reference(observations$maximum_voxel_height)),
    isTRUE(all.equal(vertical_complexity(heights, 26), log(3) / log(26), tolerance = 1e-15)))
  adaptive <- write_global_chm(cloud_path, dtm_path, file.path(root, "adaptive.tif"),
    file.path(root, "adaptive_chunks"), 1, 0.2, Inf, 0.5, 1L, raster_chunk_size = 10)
  stopifnot(adaptive$maximum_voxel_height == observations$maximum_voxel_height)

  # A 17-digit string survives fwrite/fread and JSON without making a 26 m
  # reference look inconsistent with an apparently exact 25 m maximum.
  exact_path <- file.path(root, "precision.csv")
  fwrite(data.table(observed_maximum_voxel_height_m_exact = sprintf("%.17g", max(heights))), exact_path)
  exact <- fread(exact_path, colClasses = c(observed_maximum_voxel_height_m_exact = "character"))
  serialized <- jsonlite::fromJSON(jsonlite::toJSON(as.list(exact[1]), auto_unbox = TRUE))
  stopifnot(is.character(serialized$observed_maximum_voxel_height_m_exact),
    as.numeric(serialized$observed_maximum_voxel_height_m_exact) == max(heights))

  # Exercise argparse's Python bridge: an R Inf literal as the parser default
  # crashes before argument validation, even though direct helper calls pass.
  cloud <- file.path(root, "parser_only.laz")
  file.create(cloud)
  code <- paste(
    'Sys.setenv(FORESTSTRUCTURE_SOURCE_ONLY="1"); source("run.R");',
    'stopifnot(identical(parse_parameters()$maximum_height, Inf))'
  )
  output <- system2(file.path(R.home("bin"), "Rscript"), c(
    "-e", shQuote(code), "--point-cloud", shQuote(cloud),
    "--dataset-id", "150", "--output-dir", shQuote(root)
  ), stdout = TRUE, stderr = TRUE)
  stopifnot(is.null(attr(output, "status")))
})
