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

# 3DT-3029: a known CRS is used only when no header carries one.
utm <- supplied_crs_wkt("EPSG:32610")
supplied <- resolve_analysis_crs(NA_character_, NA_character_, utm)
stopifnot(identical(supplied$source, "supplied"), identical(supplied$wkt, utm))
agreeing <- resolve_analysis_crs(utm, NA_character_, supplied_crs_wkt("EPSG:32610"))
stopifnot(identical(agreeing$source, "point_cloud"))
conflict <- tryCatch(
  resolve_analysis_crs(utm, NA_character_, supplied_crs_wkt("EPSG:25832")),
  error = function(error) conditionMessage(error)
)
stopifnot(is.character(conflict), grepl("differs", conflict))
stopifnot(inherits(try(supplied_crs_wkt("EPSG:4326"), silent = TRUE), "try-error"))
stopifnot(inherits(try(supplied_crs_wkt("EPSG:4978"), silent = TRUE), "try-error"))
stopifnot(is.na(supplied_crs_wkt(NULL)))

# 3DT-3026: Standardization's geocentric placeholder means local coordinates.
geocentric <- sf::st_crs(4978)$wkt
stopifnot(is_geocentric_crs(geocentric), !is_geocentric_crs(utm))
local <- resolve_analysis_crs(geocentric, NA_character_, utm)
stopifnot(is.na(local$wkt), identical(local$source, "local"))

if (length(runtime_arguments) == 2L) {
  point_cloud_crs <- spatial_reference_wkt(readLAScatalog(runtime_arguments[[1]]))
  original_crs <- spatial_reference_wkt(readLAScatalog(runtime_arguments[[2]]))
  real_headers <- resolve_analysis_crs(point_cloud_crs, original_crs)
  # The fixture's companion carries the geocentric placeholder: local coordinates.
  stopifnot(identical(real_headers$source, "local"))
}

message("original companion CRS fallback acceptance passed")
