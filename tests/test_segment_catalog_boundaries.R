# Canonical input must retain points exactly at maximum XY, independently of
# adaptive execution windows. The old positive-overlap catalog index lost them.
local({
  root <- tempfile("segment_catalog_boundaries_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  xyz <- expand.grid(X = c(0, 30, 60), Y = c(0, 30, 60))
  xyz$Z <- as.numeric(seq_len(nrow(xyz)))
  xyz <- rbind(xyz, data.frame(X = 59.999, Y = 30, Z = 1.5))
  cloud <- lidR::LAS(xyz)
  sf::st_crs(cloud) <- 32632
  cloud <- lidR::add_lasattribute(cloud, as.numeric(seq_len(nrow(xyz)) - 1L),
                                  "OriginalPointIndex", "source order")
  cloud <- lidR::add_lasattribute(cloud, rep(1, nrow(xyz)), "PredInstance_FM", "tree ID")
  path <- file.path(root, "boundary.las")
  lidR::writeLAS(cloud, path)
  dtm <- terra::rast(xmin = -1, xmax = 62, ymin = -1, ymax = 62, res = 1,
                    crs = "EPSG:32632")
  terra::values(dtm) <- 0
  dtm_path <- file.path(root, "dtm.tif")
  terra::writeRaster(dtm, dtm_path)
  parameters <- list(voxel_resolution = 0.2, minimum_tree_voxels = 1L,
    apex_minimum_height = 0, minimum_tree_thickness = 0,
    minimum_occupied_layers = 1L, memory_budget_gib = Inf)
  segments <- function(size, workers) {
    store <- accumulate_segments(path, dtm_path, size, 0.2, 70, "PredInstance_FM",
      workers, file.path(root, paste0("segments-", size)), bucket_count = 2L,
      base_selection = JULIA_MEMORY_SAFE_SEGMENT_SELECTION,
      point_order_dimension = "OriginalPointIndex")
    records <- lapply(list.files(file.path(root, paste0("segments-", size), "observations"),
                                full.names = TRUE), readRDS)
    stopifnot(sum(vapply(records, `[[`, numeric(1), "point_count")) == nrow(xyz))
    finalize_segment_store(store, parameters)
  }
  complete <- segments(60, 1L)
  adaptive <- segments(23, 2L)
  stopifnot(identical(complete, adaptive), complete$n_vox == nrow(xyz),
            complete$apex_x == 60, complete$apex_y == 60, complete$apex_z == 9)
  canopy <- function(size, workers, source_path = path, label = as.character(size)) {
    output <- file.path(root, paste0("chm-", label, ".tif"))
    write_global_chm(source_path, dtm_path, output, file.path(root, paste0("chm-work-", label)),
      size, 0.2, 70, 0.5, workers, point_order_dimension = "OriginalPointIndex",
      raster_chunk_size = 60)
    terra::rast(output)
  }
  full_chm <- canopy(60, 1L)
  adaptive_chm <- canopy(23, 2L)
  stopifnot(terra::compareGeom(full_chm, adaptive_chm),
            identical(as.numeric(terra::values(full_chm)), as.numeric(terra::values(adaptive_chm))),
            max(terra::values(full_chm), na.rm = TRUE) > 8.9,
            identical(unname(as.vector(terra::ext(full_chm))), c(-0.5, 60.5, -0.5, 60.5)))
  # COPC writers can store tiny extrema rounding differences. Source identity
  # remains on its .001 m lattice, and adaptive pruning must use that lattice.
  rounded_path <- file.path(root, "rounded-header.las")
  file.copy(path, rounded_path)
  connection <- file(rounded_path, "r+b")
  seek(connection, 179, "start", rw = "write")
  writeBin(59.99995, connection, size = 8L, endian = "little")
  close(connection)
  adjustment <- catalog_boundary_read_adjustment(rounded_path)
  stopifnot(adjustment[[1]]$lower > 59.999, adjustment[[1]]$lower < 59.99995)
  catalog <- catalog_with_source_boundaries(lidR::readLAScatalog(rounded_path), rounded_path)
  lidR::opt_chunk_size(catalog) <- 60
  lidR::opt_chunk_buffer(catalog) <- 0
  lidR::opt_chunk_alignment(catalog) <- c(0, 0)
  lidR::opt_progress(catalog) <- FALSE
  chunks <- lidR:::engine_chunks(catalog, FALSE)
  terminal <- Filter(function(chunk) sf::st_bbox(chunk)$xmin == 60 && sf::st_bbox(chunk)$ymin == 0, chunks)[[1L]]
  actual <- adjust_catalog_boundary_read(terminal, adjustment)
  query <- catalog_task_window(actual)
  stopifnot(query[[1]] > 59.999, query[[1]] < 59.99995,
            identical(sf::st_bbox(actual), sf::st_bbox(terminal)))
  terminal_points <- lidR::readLAS(actual)
  stopifnot(lidR::npoints(terminal_points) == 2, all(terminal_points$X == 60))
  rounded_chm <- canopy(23, 2L, rounded_path, "rounded-header")
  stopifnot(terra::compareGeom(full_chm, rounded_chm),
    identical(as.numeric(terra::values(full_chm)), as.numeric(terra::values(rounded_chm))))
})
cat("canonical segment and CHM terminal-boundary checks passed\n")
