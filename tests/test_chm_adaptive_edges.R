# Exercise the real catalog/read/voxelize/raster/mosaic path. A sparse adaptive
# read starts at Y=7.5 after voxelization, but the baseline parent starts at 0.
# Rasterizing that read independently moves the 19.2 m maximum up one cell.
local({
  root <- tempfile("chm_adaptive_edges_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  cloud <- LAS(data.frame(
    X = 1260 + c(2.7, 2.7, 2.7, 0.1, 59.99, 60, 60.01, 67, 67.01),
    Y = 2100 + c(7.49, 7.7, 8.7, 0.1, 59.99, 60, 60.01, 67, 67.01),
    Z = c(19.15, 14.55, 9.15, 1, 2, 3, 4, 5, 6)
  ))
  sf::st_crs(cloud) <- 32632
  path <- file.path(root, "edges.las")
  writeLAS(cloud, path)
  dtm <- terra::rast(xmin = 1259, xmax = 1330, ymin = 2099, ymax = 2170, res = 1)
  terra::crs(dtm) <- "EPSG:32632"
  terra::values(dtm) <- 0
  dtm_path <- file.path(root, "dtm.tif")
  terra::writeRaster(dtm, dtm_path)
  render <- function(size) {
    output <- file.path(root, sprintf("chm_%d.tif", size))
    write_global_chm(path, dtm_path, output,
                     file.path(root, sprintf("work_%d", size)),
                     size, 0.2, 70, 0.5, 1L, raster_chunk_size = 60)
    terra::rast(output)
  }
  baseline <- render(60)
  adaptive <- render(7)
  stopifnot(terra::compareGeom(baseline, adaptive))
  # Ignore band labels derived from the two different output filenames, but
  # compare every scientific value and NA exactly, without numeric tolerance.
  a <- as.numeric(terra::values(baseline))
  b <- as.numeric(terra::values(adaptive))
  if (!identical(a, b)) {
    cells <- which(xor(is.na(a), is.na(b)) | (!is.na(a) & !is.na(b) & a != b))
    print(cbind(terra::xyFromCell(baseline, cells), baseline = a[cells],
                adaptive = b[cells]))
  }
  stopifnot(identical(a, b))
})
