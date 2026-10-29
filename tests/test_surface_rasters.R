# DSM and gap-filled CHM: source run.R with FORESTSTRUCTURE_SOURCE_ONLY=1 first.

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

# Gap filling: only interior empty cells are filled, each with the mean of its
# valid 8 neighbours, repeated until none is left.
local({
  chm <- terra::rast(nrows = 10, ncols = 10, xmin = 0, xmax = 10, ymin = 0, ymax = 10,
                     crs = "EPSG:32632")
  values <- matrix(5, 10, 10)
  values[3:5, 3:5] <- NA          # interior 3x3 hole: ring first, centre second
  values[2, 2] <- 9               # a neighbour of the hole with a distinct value
  values[1:2, 8] <- NA            # gap open to the grid border
  values[8, 3] <- NA              # interior single cell
  values[7, 2] <- 1
  values[7, 9] <- NA              # interior to the grid, but next to the footprint edge
  values[, 10] <- NA              # outside the footprint
  terra::values(chm) <- as.vector(t(values))
  # Footprint: x from 0 to 9, so column 10 is outside.
  footprint <- sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(9, 0), c(9, 10), c(0, 10), c(0, 0)))))
  result <- fill_chm_interior_gaps(chm, footprint)
  filled <- matrix(terra::values(result$raster), 10, 10, byrow = TRUE)
  stopifnot(result$filled_cells == 10, result$iterations == 2L,
            result$empty_cells == sum(is.na(values)))
  # Cells with points keep their values; open gaps and outside cells stay empty.
  stopifnot(identical(filled[!is.na(values)], values[!is.na(values)]),
            all(is.na(filled[1:2, 8])), is.na(filled[7, 9]), all(is.na(filled[, 10])),
            !anyNA(filled[3:5, 3:5]), !is.na(filled[8, 3]))
  # Single cell: the mean of its eight neighbours (one of them 1, seven 5).
  stopifnot(isTRUE(all.equal(filled[8, 3], (1 + 7 * 5) / 8)))
  # Corner of the hole: five valid neighbours, one of them 9.
  stopifnot(isTRUE(all.equal(filled[3, 3], (9 + 4 * 5) / 5)))
  # Centre: filled in the second pass from its eight filled neighbours.
  ring <- filled[3:5, 3:5][-5]
  stopifnot(isTRUE(all.equal(filled[4, 4], mean(ring))))

  # Without interior gaps nothing changes and no pass runs.
  terra::values(chm) <- 5
  untouched <- fill_chm_interior_gaps(chm, footprint)
  stopifnot(untouched$filled_cells == 0, untouched$iterations == 0L,
            all(terra::values(untouched$raster) == 5))
})

message("DSM and gap-filled CHM checks passed")
