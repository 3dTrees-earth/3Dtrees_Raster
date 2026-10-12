# Source run_julia_memory_safe.R in source-only mode before this test.
local({
  root <- tempfile("companion_bounds_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  source <- file.path(root, "source.las")
  candidate <- file.path(root, "candidate.las")
  cloud <- lidR::LAS(data.frame(X = c(1000,1000.5),
                                Y = c(2000,2000.5), Z = c(0,1)))
  lidR::writeLAS(cloud, source)
  file.copy(source, candidate)
  patch_max_x <- function(path, value) {
    connection <- file(path, "r+b")
    on.exit(close(connection), add = TRUE)
    seek(connection, 179, "start", rw = "write")
    writeBin(value, connection, size = 8L, endian = "little")
  }
  patch_max_x(source, 1000.49996)
  stopifnot(isTRUE(validate_original_companion(candidate, source)))
  fallback <- read_analysis_footprint(source)
  reference <- read_analysis_footprint(candidate)
  stopifnot(identical(sf::st_bbox(fallback$usable), sf::st_bbox(reference$usable)),
            as.numeric(sf::st_area(fallback$usable)) == .25)
  patch_max_x(source, 1000.49)
  error <- tryCatch(validate_original_companion(candidate, source), error = identity)
  stopifnot(inherits(error, "error"))
})
