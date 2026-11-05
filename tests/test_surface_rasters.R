# DSM: source run.R with FORESTSTRUCTURE_SOURCE_ONLY=1 first.

# Native DSM: the highest absolute point Z per CHM cell, from the points the
# CHM uses, keyed by the voxel-column centres the CHM is rasterized from.
local({
  set.seed(7)
  count <- 4000L
  points <- data.frame(
    X = runif(count, 0, 10), Y = runif(count, 0, 10), Z = 100 + runif(count, 0, 25)
  )
  # An empty strip makes some cells nodata in both rasters.
  points <- points[points$X < 4 | points$X > 5.5, ]
  las <- suppressMessages(lidR::LAS(points))
  sf::st_crs(las) <- 32632
  dtm <- terra::rast(xmin = -5, xmax = 15, ymin = -5, ymax = 15, resolution = 1,
                     crs = "EPSG:32632")
  terra::values(dtm) <- 100
  normalized <- lidR::normalize_height(las, dtm)
  voxels <- lidR::voxelize_points(normalized, res = 0.2)
  chm <- lidR::rasterize_canopy(voxels, res = 0.5, algorithm = lidR::p2r())
  dsm <- dsm_chunk_raster(chm, dsm_column_maxima(normalized@data, 0.2),
                          voxels@header, sf::st_crs(las))
  stopifnot(terra::compareGeom(chm, dsm),
            identical(is.na(terra::values(dsm)), is.na(terra::values(chm))),
            any(is.na(terra::values(chm))))
  # Independent maximum per cell of each point's voxel-column centre.
  cells <- terra::cellFromXY(chm, cbind(
    lidR:::f_grid(normalized@data$X, 0.2, 0), lidR:::f_grid(normalized@data$Y, 0.2, 0)))
  expected <- rep(NA_real_, terra::ncell(chm))
  maxima <- tapply(normalized@data$Zref, cells, max)
  expected[as.integer(names(maxima))] <- as.numeric(maxima)
  # lidR's raster holds single precision, as the published Float32 GeoTIFF does.
  actual <- as.vector(terra::values(dsm))
  stopifnot(identical(is.na(actual), is.na(expected)),
            max(abs(actual - expected), na.rm = TRUE) < 1e-4)
  # DSM minus the flat DTM differs from the voxelized CHM by at most half a voxel.
  difference <- terra::values(dsm) - 100 - terra::values(chm)
  stopifnot(max(abs(difference), na.rm = TRUE) <= 0.1 + 1e-9)
  # A chunk without points still writes an all-nodata DSM on the CHM grid.
  empty <- dsm_chunk_raster(chm, NULL, NULL, sf::st_crs(las))
  stopifnot(terra::compareGeom(chm, empty), all(is.na(terra::values(empty))))
})
message("DSM checks passed")
