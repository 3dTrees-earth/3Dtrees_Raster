Sys.setenv(RASTER_SOURCE_ONLY = "1")

runtime_arguments <- commandArgs(trailingOnly = TRUE)
arguments <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", arguments, value = TRUE)
script_directory <- dirname(normalizePath(sub("^--file=", "", script_argument[[1]])))
run_script <- normalizePath(file.path(script_directory, "..", "src", "analysis.R"))
commandArgs <- function(trailingOnly = FALSE) {
  if (trailingOnly) character() else paste0("--file=", run_script)
}
source(run_script)
rm(commandArgs)

primary <- resolve_analysis_crs("PRIMARY_WKT", "COMPANION_WKT")
stopifnot(identical(primary$wkt, "PRIMARY_WKT"))
stopifnot(identical(primary$source, "point_cloud"))

fallback <- resolve_analysis_crs(NA_character_, "COMPANION_WKT")
stopifnot(identical(fallback$wkt, "COMPANION_WKT"))
stopifnot(identical(fallback$source, "original_companion"))

missing <- resolve_analysis_crs(NA_character_, NA_character_)
stopifnot(is.na(missing$wkt))
stopifnot(identical(missing$source, "unavailable"))

if (length(runtime_arguments) == 2L) {
  point_cloud_crs <- spatial_reference_wkt(readLAScatalog(runtime_arguments[[1]]))
  original_crs <- spatial_reference_wkt(readLAScatalog(runtime_arguments[[2]]))
  real_headers <- resolve_analysis_crs(point_cloud_crs, original_crs)
  stopifnot(identical(real_headers$source, "original_companion"))
  stopifnot(identical(sf::st_crs(real_headers$wkt)$epsg, 4978L))
}

message("original companion CRS fallback acceptance passed")
