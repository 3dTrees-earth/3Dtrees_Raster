# Optional observational export. The existing global DTM and usable AOI are
# shared; Julia's index filters and order-sensitive calculations stay separate.
voxel_count_paths <- function(directory, dataset_id) {
  prefix <- file.path(directory, paste0(dataset_id, "_voxel_counts"))
  list(csv = paste0(prefix, ".csv"), raster = paste0(prefix, ".tif"),
       metadata = paste0(prefix, ".json"), mask = paste0(prefix, "_validity.tif"))
}

add_voxel_count_arguments <- function(parser, mandatory = FALSE) {
  parser$add_argument("--export-voxel-counts", action = "store_true", default = mandatory,
                      help = if (mandatory) "Already enabled; accepted for compatibility" else
                        "Export full-source voxel point counts as CSV and GeoTIFF")
  parser$add_argument("--voxel-count-chunk-points", type = "integer", default = 2000000L,
                      help = "Maximum source records per voxel-count batch (default: 2000000)")
}

validate_voxel_count_parameters <- function(parameters) {
  n <- parameters$voxel_count_chunk_points
  if (length(n) != 1L || !is.finite(n) || n < 1L) {
    stop("--voxel-count-chunk-points must be positive")
  }
}

write_voxel_counts <- function(parameters, dtm_path, aoi, work_directory) {
  staging <- tempfile(".voxel_incoming_", tmpdir = parameters$output_dir)
  dir.create(staging)
  on.exit(unlink(staging, recursive = TRUE, force = TRUE), add = TRUE)
  paths <- voxel_count_paths(staging, parameters$dataset_id)
  chunk_points <- parameter_or(parameters, "voxel_count_chunk_points", 2000000L)
  if (!is.finite(chunk_points) || chunk_points < 1L) {
    stop("--voxel-count-chunk-points must be positive")
  }
  # The helper runs beside the R process, so reserve the parent's live memory.
  budget <- parameter_or(parameters, "memory_budget_gib", Inf)
  available <- if (is.finite(budget)) budget - current_rss_bytes() / 1024^3 else 1
  assert_memory_budget(256 * 1024^2, parameters, "voxel point-count export")
  usable_path <- file.path(work_directory, "voxel_usable.gpkg")
  sf::st_write(sf::st_sf(geometry = aoi$usable), usable_path, quiet = TRUE)
  config <- list(
    point_cloud = parameter_or(parameters, "voxel_point_cloud", parameters$point_cloud),
    dtm = dtm_path, aoi = usable_path,
    temp_dir = work_directory, paths = paths, chunk_points = chunk_points,
    source_name = basename(parameter_or(parameters, "voxel_point_cloud", parameters$point_cloud)),
    source_sha256 = parameter_or(parameters, "point_cloud_sha256", NULL),
    dtm_denoising = if (is.null(parameters$denoising)) NULL else parameters$denoising$summary,
    memory_available_gib = available,
    resolution = parameter_or(parameters, "voxel_count_resolution", 0.5)
  )
  config_path <- file.path(work_directory, "voxel_config.json")
  jsonlite::write_json(config, config_path, auto_unbox = TRUE, null = "null")
  status <- system2("python3", c(shQuote(file.path(script_directory, "voxel_counts.py")),
                                 shQuote(config_path)))
  if (!identical(status, 0L)) stop("voxel point-count export failed")
  artifacts <- validate_voxel_count_artifacts(staging, parameters$dataset_id,
    parameter_or(parameters, "voxel_count_resolution", 0.5))
  metadata <- jsonlite::read_json(paths$metadata, simplifyVector = TRUE)
  promote_artifacts(artifacts, parameters$output_dir, parameters$dataset_id,
                    replace_dataset = FALSE)
  metadata
}

validate_voxel_count_artifacts <- function(directory, dataset_id, resolution = 0.5) {
  paths <- voxel_count_paths(directory, dataset_id)
  if (any(!file.exists(unlist(paths))) || any(file.info(unlist(paths))$size <= 0)) {
    stop("voxel-count artifacts missing or empty")
  }
  metadata <- jsonlite::read_json(paths$metadata, simplifyVector = TRUE)
  cube <- terra::rast(paths$raster)
  mask <- terra::rast(paths$mask)
  if (!nzchar(terra::crs(cube)) || !nzchar(terra::crs(mask))) {
    stop("voxel-count GeoTIFF CRS is missing; it must persist without sidecars")
  }
  if (terra::nlyr(cube) != metadata$band_count ||
      any(abs(terra::res(cube) - resolution) > 1e-9) ||
      any(abs(unlist(metadata$voxel_size_m) - resolution) > 1e-9) ||
      !terra::compareGeom(cube[[1]], mask, stopOnError = FALSE) ||
      !nzchar(terra::crs(cube)) ||
      metadata$point_accounting$read != metadata$source$point_count ||
      metadata$point_accounting$included < 0 ||
      sum(unlist(metadata$point_accounting)[names(metadata$point_accounting) != "read"]) !=
        metadata$point_accounting$read ||
      terra::ncol(cube) != metadata$grid$width ||
      terra::nrow(cube) != metadata$grid$height ||
      !identical(names(data.table::fread(paths$csv, nrows = 0)),
                 c("x_center", "y_center", "height_lower_m", "height_upper_m", "point_count"))) {
    stop("voxel-count grid or provenance validation failed")
  }
  invisible(unlist(paths))
}
