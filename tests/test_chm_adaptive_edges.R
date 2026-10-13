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
  catalog <- readLAScatalog(path)
  opt_chunk_size(catalog) <- 60
  opt_chunk_alignment(catalog) <- c(0, 0)
  opt_chunk_buffer(catalog) <- 0
  opt_select(catalog) <- "xyz"
  opt_progress(catalog) <- FALSE
  parent <- lidR:::engine_chunks(catalog, FALSE)[[1]]
  children <- adaptive_chm_read_clusters(parent, 7)
  child_bounds <- function(values) do.call(rbind, lapply(values, function(child)
    as.numeric(sf::st_bbox(child))))
  actual_bounds <- child_bounds(children)
  stopifnot(length(children) == 81L,
    identical(actual_bounds[1, ], c(1260, 2100, 1267, 2107)),
    identical(actual_bounds[2, ], c(1267, 2100, 1274, 2107)),
    identical(actual_bounds[10, ], c(1260, 2107, 1267, 2114)),
    identical(actual_bounds[81, ], c(1316, 2156, 1320, 2160)),
    all(vapply(children, function(child)
      identical(child@select, parent@select) &&
      identical(child@files, parent@files) &&
      identical(child@index, parent@index), logical(1))))
  # Entirely disjoint children are omitted, but a child starting exactly at
  # the source maximum remains: its half-open read owns that boundary point.
  bounded <- adaptive_chm_read_clusters(parent, 7,
    c(xmin = 1267, ymin = 2107, xmax = 1274, ymax = 2114))
  stopifnot(identical(child_bounds(bounded), rbind(
    c(1267, 2107, 1274, 2114), c(1274, 2107, 1281, 2114),
    c(1267, 2114, 1274, 2121), c(1274, 2114, 1281, 2121))))
  stopifnot(length(adaptive_chm_read_clusters(parent, 7,
    c(xmin = 2000, ymin = 3000, xmax = 2010, ymax = 3010))) == 0L)
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
