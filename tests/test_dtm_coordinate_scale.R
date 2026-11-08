ground <- expand.grid(
  X = 292142 + seq(0, 0.009, by = 0.001),
  Y = 6456591 + seq(0, 0.009, by = 0.001)
)
ground$Z <- 0
ground$Classification <- 2L

las <- LAS(ground)
las <- las_reoffset(las, xoffset = 292142, yoffset = 6456591, zoffset = 0)
las <- las_rescale(las, xscale = 0.00025, yscale = 0.00025, zscale = 0.001)

# The original implementation passes the classified LAS directly to TIN.
# Preserve that first attempt, but retry the specific lidR integer-conversion
# failure with the legacy floating-point TIN and the unchanged LAS coordinates.
stopifnot(exists("legacy_tin", inherits = FALSE))
surface_body <- paste(deparse(body(classified_ground_surface)), collapse = "\n")
stopifnot(grepl("legacy_tin", surface_body, fixed = TRUE))
stopifnot(!grepl("tin_compatible_xy_scale", surface_body, fixed = TRUE))

coordinates_before <- las@data[, c("X", "Y")]
scale_before <- unlist(las@header@PHB[c("X scale factor", "Y scale factor")])
offset_before <- unlist(las@header@PHB[c("X offset", "Y offset")])
dtm <- classified_ground_surface(las, 1)
stopifnot(inherits(dtm, "SpatRaster"))
stopifnot(identical(las@data[, c("X", "Y")], coordinates_before))
stopifnot(identical(
  unlist(las@header@PHB[c("X scale factor", "Y scale factor")]),
  scale_before
))
stopifnot(identical(
  unlist(las@header@PHB[c("X offset", "Y offset")]),
  offset_before
))

fine_ground <- expand.grid(
  X = 226990 + seq(0, 0.09, by = 0.01),
  Y = 378990 + seq(0, 0.09, by = 0.01)
)
fine_ground$Z <- 357
fine_ground$Classification <- 2L
fine_las <- LAS(fine_ground)
fine_las@header@PHB[["X scale factor"]] <- 1e-7
fine_las@header@PHB[["Y scale factor"]] <- 1e-7
fine_coordinates_before <- fine_las@data[, c("X", "Y")]
fine_scale_before <- unlist(
  fine_las@header@PHB[c("X scale factor", "Y scale factor")]
)
fine_dtm <- classified_ground_surface(fine_las, 1)
stopifnot(inherits(fine_dtm, "SpatRaster"))
stopifnot(identical(fine_las@data[, c("X", "Y")], fine_coordinates_before))
stopifnot(identical(
  unlist(fine_las@header@PHB[c("X scale factor", "Y scale factor")]),
  fine_scale_before
))

# Exercise the exact error-dispatch seam deterministically. The first mocked
# rasterization represents lidR's optimized integer TIN failure; the second
# call must receive the explicit legacy algorithm without changing the LAS.
original_rasterize_terrain <- rasterize_terrain
original_legacy_tin <- legacy_tin
on.exit({
  rasterize_terrain <- original_rasterize_terrain
  legacy_tin <- original_legacy_tin
}, add = TRUE)
rasterize_calls <- 0L
legacy_retry_called <- FALSE
legacy_tin <- function(...) {
  legacy_retry_called <<- TRUE
  original_legacy_tin(...)
}
rasterize_terrain <- function(las, res, algorithm) {
  rasterize_calls <<- rasterize_calls + 1L
  if (rasterize_calls == 1L) {
    stop(paste(
      "Internal error in C_interpolate_delaunay:",
      "xy coordinates were not converted to integer"
    ))
  }
  force(algorithm)
  "legacy-retry-result"
}
retry_coordinates_before <- fine_las@data[, c("X", "Y")]
retry_header_before <- fine_las@header@PHB
retry_result <- classified_ground_surface(fine_las, 1)
stopifnot(identical(retry_result, "legacy-retry-result"))
stopifnot(rasterize_calls == 2L)
stopifnot(legacy_retry_called)
stopifnot(identical(fine_las@data[, c("X", "Y")], retry_coordinates_before))
stopifnot(identical(fine_las@header@PHB, retry_header_before))
