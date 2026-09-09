suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(future)
  library(lidR)
  library(sf)
  library(terra)
})

# Runtime invariants chosen for the validated production image. Catalog
# workers use one lidR thread; segment aggregation uses a small, fixed
# data.table pool and deterministic disk buckets to bound peak memory.
DATA_TABLE_THREADS <- 5L
SEGMENT_BUCKET_COUNT <- 64L
MAX_LIDR_POINT_COUNT <- .Machine$integer.max
DTM_CATALOG_SELECTION <- "xyz"
SEGMENT_CATALOG_BASE_SELECTION <- "xyzrn"
JULIA_MEMORY_SAFE_SEGMENT_SELECTION <- "xyz"

# COPC stores points in spatial hierarchy order, while the original Julia script
# implicitly resolves several equal-value ties by input record order.  Sort each
# independently read chunk by the persisted source ordinal before applying any
# order-sensitive DTM, CHM, segmentation, or tile operation.  Plain LAS/LAZ runs
# pass NULL and retain their already-original order without extra work.
restore_original_point_order <- function(las, point_order_dimension = NULL,
                                         phase = "analysis") {
  # Preserve the historical path when no key was requested, and avoid touching
  # genuinely empty chunks returned by LAScatalog edge queries.
  if (is.null(point_order_dimension) || is.empty(las)) return(las)
  # A missing key means selective reading or COPC construction violated the
  # canonical input contract; continuing would produce plausible but wrong ties.
  if (!point_order_dimension %in% names(las@data)) {
    stop(sprintf(
      "%s did not load required point-order dimension %s",
      phase,
      point_order_dimension
    ))
  }
  # Missing ordinals make deterministic reconstruction impossible.  Do not bind
  # the complete column separately: an extra reference can prevent in-place work
  # on very dense windows.
  if (anyNA(las@data[[point_order_dimension]])) {
    stop(sprintf(
      "%s received missing or duplicate %s values",
      phase,
      point_order_dimension
    ))
  }
  # Reorder the complete point table in place so all dimensions remain aligned.
  # A row subset would allocate a second full LAS table; on dense windows that
  # copy alone can consume several GiB and breach the internal memory budget.
  data.table::setorderv(las@data, point_order_dimension)
  # Once sorted, strict monotonicity detects duplicate ordinals in one linear C
  # scan without the large hash table allocated by anyDuplicated().
  if (is.unsorted(las@data[[point_order_dimension]], strictly = TRUE)) {
    stop(sprintf(
      "%s received missing or duplicate %s values",
      phase,
      point_order_dimension
    ))
  }
  # Return the same LAS object shape expected by existing scientific functions.
  las
}

parameter_or <- function(parameters, name, default) {
  value <- parameters[[name]]
  if (is.null(value)) default else value
}

validated_catalog_worker_count <- function(value, thread_budget) {
  workers <- suppressWarnings(as.integer(value))
  if (length(workers) != 1L || is.na(workers) ||
      !workers %in% c(1L, 2L)) {
    stop("FORESTSTRUCTURE_CATALOG_WORKERS must be 1 or 2")
  }
  if (workers > thread_budget) {
    stop("FORESTSTRUCTURE_CATALOG_WORKERS exceeds FORESTSTRUCTURE_THREADS")
  }
  workers
}

current_rss_bytes <- function() {
  status <- tryCatch(
    readLines("/proc/self/status", warn = FALSE),
    error = function(error) character()
  )
  line <- grep("^VmRSS:", status, value = TRUE)
  if (length(line) != 1L) return(0)
  kib <- suppressWarnings(as.numeric(sub(
    "^VmRSS:[[:space:]]*([0-9]+).*$",
    "\\1",
    line
  )))
  if (!is.finite(kib)) 0 else kib * 1024
}

assert_memory_budget <- function(required_bytes, parameters, operation,
                                 single_segment = FALSE) {
  budget_gib <- parameter_or(parameters, "memory_budget_gib", Inf)
  if (!is.finite(budget_gib)) return(invisible(TRUE))
  rss <- current_rss_bytes()
  projected <- rss + max(0, required_bytes)
  if (projected <= budget_gib * 1024^3) return(invisible(TRUE))
  flag <- if (single_segment) {
    "single_segment_exceeds_memory_budget"
  } else {
    "memory_budget_exceeded"
  }
  stop(sprintf(
    paste0(
      "FORESTSTRUCTURE_FLAG %s: %s requires an estimated %.3f GiB ",
      "at %.3f GiB RSS with a %.3f GiB internal budget"
    ),
    flag,
    operation,
    required_bytes / 1024^3,
    rss / 1024^3,
    budget_gib
  ))
}

assert_loaded_dimensions <- function(las, allowed, phase) {
  present <- names(las@data)
  unexpected <- setdiff(present, allowed)
  missing <- setdiff(c("X", "Y", "Z"), present)
  if (length(unexpected) > 0L || length(missing) > 0L) {
    stop(sprintf(
      paste0(
        "FORESTSTRUCTURE_FLAG selective_read_violation: %s loaded [%s]; ",
        "allowed [%s]"
      ),
      phase,
      paste(present, collapse = ","),
      paste(allowed, collapse = ",")
    ))
  }
  invisible(TRUE)
}

configure_catalog_worker <- function() {
  lidR::set_lidr_threads(1L)
  data.table::setDTthreads(DATA_TABLE_THREADS)
  invisible(NULL)
}

read_las_uint_le <- function(connection, bytes) {
  raw_value <- readBin(connection, what = "raw", n = bytes)
  if (length(raw_value) != bytes) stop("LAS/LAZ header is truncated")
  sum(as.numeric(raw_value) * 256^(seq_len(bytes) - 1L))
}

point_count_from_las_header <- function(path) {
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  signature <- rawToChar(readBin(connection, what = "raw", n = 4L))
  if (!identical(signature, "LASF")) {
    stop("point cloud does not have a LAS/LAZ header")
  }

  seek(connection, where = 24L, origin = "start")
  version <- as.integer(readBin(connection, what = "raw", n = 2L))
  if (length(version) != 2L) stop("LAS/LAZ header is truncated")

  seek(connection, where = 107L, origin = "start")
  legacy_count <- read_las_uint_le(connection, 4L)
  is_las_14 <- version[[1]] > 1L ||
    (version[[1]] == 1L && version[[2]] >= 4L)
  if (!is_las_14) return(legacy_count)

  seek(connection, where = 247L, origin = "start")
  extended_count <- read_las_uint_le(connection, 8L)
  if (extended_count > 0) extended_count else legacy_count
}

assert_lidr_point_count_supported <- function(point_count) {
  if (!is.finite(point_count) || point_count < 0) {
    stop("LAS/LAZ header contains an invalid point count")
  }
  if (point_count > MAX_LIDR_POINT_COUNT) {
    stop(sprintf(
      paste(
        "FORESTSTRUCTURE_FLAG unsupported_lidr_point_count:",
        "single point-cloud file contains %.0f points; lidR/R supports at most %.0f"
      ),
      point_count,
      MAX_LIDR_POINT_COUNT
    ))
  }
  invisible(point_count)
}

# Translate a point-count target into the square metre-based window understood
# by LAScatalog.  COPC makes those spatial reads efficient; plain LAS/LAZ keeps
# its established fixed-size behavior so this runtime control cannot silently
# change the legacy scientific path.
adaptive_chunk_plan <- function(point_count, xy_bounds, configured_chunk_size,
                                target_chunk_points, is_copc,
                                alignment_m = 1) {
  if (!is.finite(configured_chunk_size) || configured_chunk_size <= 0) {
    stop("configured chunk size must be greater than zero")
  }
  if (!is.finite(target_chunk_points) || target_chunk_points < 0) {
    stop("target chunk points must be zero or greater")
  }
  disabled <- !isTRUE(is_copc) || target_chunk_points == 0
  if (disabled) {
    return(list(
      enabled = FALSE,
      method = "fixed_metres",
      target_points = target_chunk_points,
      alignment_m = alignment_m,
      chunk_size_m = configured_chunk_size,
      estimated_points_per_full_chunk = NA_real_
    ))
  }
  required_names <- c("xmin", "ymin", "xmax", "ymax")
  if (!all(required_names %in% names(xy_bounds))) {
    stop("adaptive chunk planning requires named XY bounds")
  }
  width <- as.numeric(xy_bounds[["xmax"]] - xy_bounds[["xmin"]])
  height <- as.numeric(xy_bounds[["ymax"]] - xy_bounds[["ymin"]])
  area <- width * height
  if (!is.finite(area) || area <= 0) {
    stop("adaptive chunk planning requires a positive XY area")
  }
  if (!is.finite(point_count) || point_count < 0) {
    stop("adaptive chunk planning requires a non-negative point count")
  }
  if (!is.finite(alignment_m) || alignment_m <= 0) {
    stop("adaptive chunk alignment must be greater than zero")
  }
  density <- point_count / area
  ideal_size <- if (density == 0) {
    configured_chunk_size
  } else {
    sqrt(target_chunk_points / density)
  }
  # Sixty metres is the scientifically validated upper bound.  A smaller
  # --chunk-size remains an additional operator-supplied ceiling, while a
  # larger value cannot silently widen the established window geometry.
  maximum_size <- min(configured_chunk_size, 60)
  aligned_size <- floor(min(ideal_size, maximum_size) / alignment_m) *
    alignment_m
  chunk_size <- max(alignment_m, aligned_size)
  list(
    enabled = TRUE,
    method = "mean_density",
    target_points = target_chunk_points,
    alignment_m = alignment_m,
    point_count = point_count,
    width_m = width,
    height_m = height,
    area_m2 = area,
    mean_density_points_m2 = density,
    ideal_chunk_size_m = ideal_size,
    configured_maximum_chunk_size_m = configured_chunk_size,
    effective_maximum_chunk_size_m = maximum_size,
    chunk_size_m = chunk_size,
    estimated_points_per_full_chunk = density * chunk_size^2
  )
}

script_argument <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
script_directory <- if (length(script_argument)) {
  dirname(normalizePath(sub("^--file=", "", script_argument[[1]])))
} else {
  getwd()
}
source(file.path(script_directory, "tile_scheduler.R"))

# CLI and artifact naming ---------------------------------------------------

parse_parameters <- function() {
  parser <- ArgumentParser(
    description = "Compute forest-structure indices for optimized Analysis Tiles"
  )
  inputs <- parser$add_argument_group("Inputs")
  scientific <- parser$add_argument_group("Scientific parameters")
  runtime <- parser$add_argument_group("Runtime controls")
  segments <- parser$add_argument_group("Optional segment controls")

  inputs$add_argument(
    "--point-cloud", "--point_cloud",
    required = TRUE,
    help = "Path to exactly one LAS or LAZ point cloud"
  )
  inputs$add_argument(
    "--dataset-id", "--dataset_id",
    required = TRUE,
    type = "integer",
    help = "Positive 3Dtrees dataset identifier used to prefix every output"
  )
  inputs$add_argument(
    "--aoi",
    required = FALSE,
    default = NULL,
    help = paste(
      "Optional GeoJSON or GeoPackage Audit AOI; when omitted, tiles cover",
      "the complete point-cloud XY extent"
    )
  )
  inputs$add_argument(
    "--output-dir", "--output_dir",
    required = TRUE,
    help = "Existing writable output directory"
  )
  scientific$add_argument(
    "--tile-size", "--tile_size",
    type = "double",
    default = 20,
    help = "Fixed Analysis Tile width in metres (default: 20)"
  )
  scientific$add_argument(
    "--grid-search-step", "--grid_search_step",
    type = "double",
    default = 0.5,
    help = "Offset-search resolution in metres; placement is always optimized (default: 0.5)"
  )
  scientific$add_argument(
    "--ptd-resolution", "--ptd_resolution",
    type = "double",
    default = 20,
    help = "PTD seed resolution in metres (default: 20)"
  )
  scientific$add_argument(
    "--dtm-resolution", "--dtm_resolution",
    type = "double",
    default = 1,
    help = "Global DTM resolution in metres (default: 1)"
  )
  scientific$add_argument(
    "--maximum-height", "--maximum_height",
    type = "double",
    default = 70,
    help = "Maximum normalized height retained in metres (default: 70)"
  )
  scientific$add_argument(
    "--voxel-resolution", "--voxel_resolution",
    type = "double",
    default = 0.2,
    help = "Voxel edge length in metres (default: 0.2)"
  )
  scientific$add_argument(
    "--vegetation-minimum-height", "--vegetation_minimum_height",
    type = "double",
    default = 0.5,
    help = "Minimum vegetation height in metres (default: 0.5)"
  )
  scientific$add_argument(
    "--chm-resolution", "--chm_resolution",
    type = "double",
    default = 0.5,
    help = "Full-point-cloud CHM resolution in metres (default: 0.5)"
  )
  scientific$add_argument(
    "--gap-height-threshold", "--gap_height_threshold",
    type = "double",
    default = 3,
    help = "CHM gap-height threshold in metres (default: 3)"
  )
  runtime$add_argument(
    "--dtm-chunk-size", "--dtm_chunk_size",
    type = "double",
    default = 60,
    help = "LAScatalog DTM chunk width in metres (default: 60)"
  )
  runtime$add_argument(
    "--chunk-size", "--chunk_size",
    type = "double",
    default = 60,
    help = paste(
      "Maximum LAScatalog CHM/segment chunk width in metres; adaptive COPC",
      "planning may select a smaller whole-metre width (default: 60)"
    )
  )
  runtime$add_argument(
    "--target-chunk-points", "--target_chunk_points",
    type = "double",
    default = 50000000,
    help = paste(
      "Approximate COPC CHM/segment points per spatial chunk; 0 disables",
      "adaptive sizing (default: 50000000)"
    )
  )
  runtime$add_argument(
    "--dtm-buffer", "--dtm_buffer",
    type = "double",
    default = 20,
    help = "DTM chunk buffer in metres (default: 20)"
  )
  runtime$add_argument(
    "--threads",
    type = "integer",
    default = 0,
    help = "lidR thread count; 0 preserves the package/container default (default: 0)"
  )
  runtime$add_argument(
    "--catalog-workers", "--catalog_workers",
    type = "integer",
    default = 2,
    help = paste(
      "Maximum concurrent LAScatalog processes; kept separate from the lidR",
      "thread allowance to bound memory (default: 2)"
    )
  )
  runtime$add_argument(
    "--performance-report", "--performance_report",
    action = "store_true",
    default = FALSE,
    help = "Write the dataset-prefixed performance CSV with timings and peak RSS"
  )
  segments$add_argument(
    "--instance-dimension", "--instance_dimension",
    action = "append",
    default = NULL,
    help = paste(
      "Instance-ID dimensions to process; repeat the option or provide",
      "comma-separated names. Missing dimensions are skipped. Defaults to",
      "PredInstance, PredInstance_SAT, PredInstance_FM"
    )
  )
  segments$add_argument(
    "--segment-diagnostics", "--segment_diagnostics",
    action = "store_true",
    default = FALSE,
    help = "Deprecated compatibility flag; segment diagnostics are always written"
  )
  segments$add_argument(
    "--minimum-tree-voxels", "--minimum_tree_voxels",
    type = "integer",
    default = 100,
    help = "Minimum occupied voxels for an accepted tree (default: 100)"
  )
  segments$add_argument(
    "--apex-minimum-height", "--apex_minimum_height",
    type = "double",
    default = 3,
    help = "Strict minimum tree-apex height in metres (default: 3)"
  )
  segments$add_argument(
    "--minimum-tree-thickness", "--minimum_tree_thickness",
    type = "double",
    default = 0.5,
    help = "Minimum smallest PCA extent in metres (default: 0.5)"
  )
  segments$add_argument(
    "--minimum-occupied-layers", "--minimum_occupied_layers",
    type = "integer",
    default = 3,
    help = "Minimum occupied one-metre height layers (default: 3)"
  )
  args <- parser$parse_args()

  raw_instance_dimensions <- unlist(
    args$instance_dimension,
    recursive = TRUE,
    use.names = FALSE
  )
  if (is.null(raw_instance_dimensions) || length(raw_instance_dimensions) == 0 ||
      all(is.na(raw_instance_dimensions))) {
    args$instance_dimension <- c(
      "PredInstance", "PredInstance_SAT", "PredInstance_FM"
    )
  } else {
    candidates <- trimws(unlist(strsplit(
      as.character(raw_instance_dimensions),
      ",",
      fixed = TRUE
    )))
    args$instance_dimension <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
    if (length(args$instance_dimension) == 0) {
      stop("--instance-dimension must contain at least one non-empty name")
    }
  }

  if (!file.exists(args$point_cloud) || dir.exists(args$point_cloud)) {
    stop("--point-cloud must be exactly one existing LAS/LAZ file")
  }
  if (!is.finite(args$dataset_id) || args$dataset_id < 1) {
    stop("--dataset-id must be a positive integer")
  }
  if (!grepl("\\.la[sz]$", args$point_cloud, ignore.case = TRUE)) {
    stop("--point-cloud must have a .las or .laz extension")
  }
  if (!is.null(args$aoi)) {
    if (!file.exists(args$aoi) || dir.exists(args$aoi)) {
      stop("--aoi must be exactly one existing GeoJSON or GeoPackage file")
    }
    if (!grepl("\\.(geojson|json|gpkg)$", args$aoi, ignore.case = TRUE)) {
      stop("--aoi must have a .geojson, .json, or .gpkg extension")
    }
  }
  if (!dir.exists(args$output_dir)) {
    stop("--output-dir must be an existing directory")
  }
  if (file.access(args$output_dir, mode = 2) != 0) {
    stop("--output-dir must be writable")
  }

  positive <- c(
    "tile_size", "grid_search_step", "ptd_resolution", "dtm_resolution",
    "maximum_height", "voxel_resolution", "chm_resolution", "dtm_chunk_size",
    "chunk_size", "dtm_buffer"
  )
  for (name in positive) {
    if (!is.finite(args[[name]]) || args[[name]] <= 0) {
      stop(sprintf("--%s must be greater than zero", gsub("_", "-", name)))
    }
  }
  if (!is.finite(args$target_chunk_points) || args$target_chunk_points < 0) {
    stop("--target-chunk-points must be zero or greater")
  }
  if (args$grid_search_step > args$tile_size) {
    stop("--grid-search-step must not exceed --tile-size")
  }
  if (!is.finite(args$vegetation_minimum_height) || args$vegetation_minimum_height < 0) {
    stop("--vegetation-minimum-height must be zero or greater")
  }
  if (!is.finite(args$gap_height_threshold) || args$gap_height_threshold < 0) {
    stop("--gap-height-threshold must be zero or greater")
  }
  if (!is.finite(args$minimum_tree_voxels) || args$minimum_tree_voxels < 1) {
    stop("--minimum-tree-voxels must be at least one")
  }
  if (!is.finite(args$apex_minimum_height) || args$apex_minimum_height < 0) {
    stop("--apex-minimum-height must be zero or greater")
  }
  if (!is.finite(args$minimum_tree_thickness) || args$minimum_tree_thickness < 0) {
    stop("--minimum-tree-thickness must be zero or greater")
  }
  if (!is.finite(args$minimum_occupied_layers) || args$minimum_occupied_layers < 1) {
    stop("--minimum-occupied-layers must be at least one")
  }
  if (!is.finite(args$threads) || args$threads < 0) {
    stop("--threads must be zero or greater")
  }
  if (!is.finite(args$catalog_workers) || args$catalog_workers < 1) {
    stop("--catalog-workers must be at least one")
  }
  args
}

artifact_paths <- function(output_dir, dataset_id) {
  prefix <- file.path(output_dir, as.character(dataset_id))
  list(
    prefix = prefix,
    dtm = paste0(prefix, "_dtm.tif"),
    chm = paste0(prefix, "_chm.tif"),
    tiles_png = paste0(prefix, "_tiles.png"),
    performance = paste0(prefix, "_performance.csv")
  )
}

dimension_artifact_paths <- function(paths, instance_dimension = NULL) {
  suffix <- if (is.null(instance_dimension)) {
    ""
  } else {
    safe <- gsub("[^A-Za-z0-9_-]", "_", instance_dimension)
    if (!nzchar(safe)) stop("Instance Dimension cannot form a safe artifact name")
    paste0("_", safe)
  }
  list(
    results = paste0(paths$prefix, suffix, "_results.csv"),
    segment_diagnostics = paste0(
      paths$prefix, suffix, "_segment_diagnostics.csv"
    ),
    tiles_geojson = paste0(paths$prefix, suffix, "_tiles.geojson")
  )
}

# Audit footprint and deterministic tile grid ------------------------------

normalize_aoi_geometry <- function(features, label, allow_empty = FALSE) {
  if (is.null(features) || nrow(features) == 0) {
    if (allow_empty) return(st_sfc(crs = NA_crs_))
    stop(sprintf("Audit AOI contains no %s geometry", label))
  }
  geometry_types <- unique(as.character(st_geometry_type(features, by_geometry = TRUE)))
  if (!all(geometry_types %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop(sprintf("Audit AOI %s must contain only Polygon or MultiPolygon geometry", label))
  }

  # Audit AOIs are already expressed in the point cloud's local XY space.
  # GeoJSON readers may infer WGS84 from the container format; that metadata is
  # not authoritative for the local audit handoff and must not trigger a transform.
  st_crs(features) <- NA
  features <- st_make_valid(features)
  geometry <- st_union(st_geometry(features))
  if (length(geometry) == 0 || any(st_is_empty(geometry))) {
    if (allow_empty) return(st_sfc(crs = NA_crs_))
    stop(sprintf("Audit AOI %s is empty after geometry validation", label))
  }
  geometry
}

read_geojson_aoi <- function(path) {
  features <- st_read(path, quiet = TRUE, stringsAsFactors = FALSE)
  if (nrow(features) == 0) stop("Audit AOI contains no geometry")
  if (!"role" %in% names(features)) {
    return(list(inclusion = features, exclusion = NULL))
  }

  roles <- tolower(trimws(as.character(features$role)))
  aliases <- c(
    include = "include",
    inclusion = "include",
    exclude = "exclude",
    exclusion = "exclude"
  )
  unknown <- unique(roles[!roles %in% names(aliases)])
  if (length(unknown) > 0) {
    stop(sprintf("Audit AOI contains unsupported feature role(s): %s", paste(unknown, collapse = ", ")))
  }
  roles <- unname(aliases[roles])
  list(
    inclusion = features[roles == "include", , drop = FALSE],
    exclusion = features[roles == "exclude", , drop = FALSE]
  )
}

read_geopackage_aoi <- function(path) {
  layers <- st_layers(path)$name
  if (!"include" %in% layers) {
    stop("GeoPackage Audit AOI must contain an 'include' layer")
  }
  list(
    inclusion = st_read(path, layer = "include", quiet = TRUE, stringsAsFactors = FALSE),
    exclusion = if ("exclude" %in% layers) {
      st_read(path, layer = "exclude", quiet = TRUE, stringsAsFactors = FALSE)
    } else {
      NULL
    }
  )
}

read_audit_aoi <- function(path) {
  raw <- if (grepl("\\.gpkg$", path, ignore.case = TRUE)) {
    read_geopackage_aoi(path)
  } else {
    read_geojson_aoi(path)
  }
  inclusion <- normalize_aoi_geometry(raw$inclusion, "inclusion")
  exclusion <- normalize_aoi_geometry(raw$exclusion, "exclusion", allow_empty = TRUE)
  usable <- if (length(exclusion) > 0) {
    st_difference(inclusion, exclusion)
  } else {
    inclusion
  }
  if (length(usable) == 0 || all(st_is_empty(usable))) {
    stop("Audit AOI exclusions remove the complete inclusion geometry")
  }
  list(
    inclusion = inclusion,
    exclusion = exclusion,
    usable = usable,
    source = "audit_aoi"
  )
}

point_cloud_xy_bounds <- function(path) {
  header <- readLASheader(path)
  bounds <- c(
    xmin = as.numeric(header@PHB[["Min X"]]),
    ymin = as.numeric(header@PHB[["Min Y"]]),
    xmax = as.numeric(header@PHB[["Max X"]]),
    ymax = as.numeric(header@PHB[["Max Y"]])
  )
  if (any(!is.finite(bounds)) ||
      bounds[["xmax"]] <= bounds[["xmin"]] ||
      bounds[["ymax"]] <= bounds[["ymin"]]) {
    stop("point-cloud header does not contain a usable two-dimensional XY extent")
  }
  bounds
}

read_point_cloud_extent <- function(path) {
  bounds <- point_cloud_xy_bounds(path)
  ring <- matrix(c(
    bounds[["xmin"]], bounds[["ymin"]],
    bounds[["xmax"]], bounds[["ymin"]],
    bounds[["xmax"]], bounds[["ymax"]],
    bounds[["xmin"]], bounds[["ymax"]],
    bounds[["xmin"]], bounds[["ymin"]]
  ), ncol = 2, byrow = TRUE)
  extent <- st_sfc(st_polygon(list(ring)), crs = NA_crs_)
  list(
    inclusion = extent,
    exclusion = st_sfc(crs = NA_crs_),
    usable = extent,
    source = "point_cloud_extent"
  )
}

read_analysis_footprint <- function(point_cloud, aoi = NULL) {
  if (is.null(aoi)) return(read_point_cloud_extent(point_cloud))
  audit <- read_audit_aoi(aoi)
  point_cloud_extent <- read_point_cloud_extent(point_cloud)
  overlap <- suppressWarnings(st_intersection(
    audit$usable,
    point_cloud_extent$usable
  ))
  overlap_area <- if (length(overlap) == 0L) 0 else sum(as.numeric(st_area(overlap)))
  if (!is.finite(overlap_area) || overlap_area <= 0) {
    stop(
      "Audit AOI does not overlap the point-cloud XY extent; ",
      "both inputs must use the same local coordinate frame"
    )
  }
  audit
}

make_grid_at <- function(geometry, tile_size, offset_x, offset_y,
                         reference_bounds = st_bbox(geometry)) {
  st_make_grid(
    geometry,
    cellsize = tile_size,
    offset = c(
      reference_bounds[["xmin"]] - offset_x,
      reference_bounds[["ymin"]] - offset_y
    )
  )
}

count_complete_tiles <- function(grid, geometry) {
  sum(lengths(st_within(grid, geometry)) > 0)
}

build_optimized_tiles <- function(geometry, tile_size, search_step, threads = 1L) {
  offsets <- seq(0, tile_size - search_step, by = search_step)
  reference_bounds <- st_bbox(geometry)
  components <- st_cast(geometry, "POLYGON", warn = FALSE)
  candidates <- data.frame(
    offset_x = rep(offsets, each = length(offsets)),
    offset_y = rep(offsets, times = length(offsets))
  )
  maximum_possible_tiles <- floor(
    as.numeric(st_area(geometry)) / tile_size^2
  )
  count_candidate <- function(index) {
    sum(vapply(seq_along(components), function(component_index) {
      component <- components[component_index]
      count_complete_tiles(
        make_grid_at(
          component,
          tile_size,
          candidates$offset_x[[index]],
          candidates$offset_y[[index]],
          reference_bounds
        ),
        component
      )
    }, integer(1)))
  }

  first_count <- count_candidate(1L)
  if (first_count >= maximum_possible_tiles || nrow(candidates) == 1L) {
    best_index <- 1L
  } else {
    remaining <- seq.int(2L, nrow(candidates))
    workers <- min(max(1L, as.integer(threads)), length(remaining))
    message(sprintf(
      "Optimizing %d remaining grid offsets with %d worker(s)",
      length(remaining),
      workers
    ))
    remaining_counts <- if (workers == 1L) {
      vapply(remaining, count_candidate, integer(1))
    } else {
      unlist(
        parallel::mclapply(
          remaining,
          count_candidate,
          mc.cores = workers,
          mc.preschedule = TRUE,
          mc.set.seed = FALSE
        ),
        use.names = FALSE
      )
    }
    counts <- c(first_count, remaining_counts)
    best_index <- which.max(counts)
  }

  best_x <- candidates$offset_x[[best_index]]
  best_y <- candidates$offset_y[[best_index]]
  grid <- make_grid_at(geometry, tile_size, best_x, best_y)
  complete <- lengths(st_within(grid, geometry)) > 0
  if (!any(complete)) {
    return(st_sf(
      tile_id = integer(),
      geometry = st_sfc(crs = st_crs(geometry))
    ))
  }

  tiles <- st_sf(geometry = grid[complete])
  bounds <- t(vapply(seq_len(nrow(tiles)), function(index) {
    as.numeric(st_bbox(tiles[index, ]))
  }, numeric(4)))
  order_index <- order(bounds[, 2], bounds[, 1])
  tiles <- tiles[order_index, , drop = FALSE]
  tiles$tile_id <- seq_len(nrow(tiles))
  tiles[, c("tile_id", "geometry")]
}

build_covering_tiles <- function(geometry, tile_size) {
  bounds <- st_bbox(geometry)
  span_x <- as.numeric(bounds[["xmax"]] - bounds[["xmin"]])
  span_y <- as.numeric(bounds[["ymax"]] - bounds[["ymin"]])
  columns <- as.integer(ceiling(span_x / tile_size))
  rows <- as.integer(ceiling(span_y / tile_size))
  covered_width <- columns * tile_size
  covered_height <- rows * tile_size
  origin_x <- as.numeric(bounds[["xmin"]] - (covered_width - span_x) / 2)
  origin_y <- as.numeric(bounds[["ymin"]] - (covered_height - span_y) / 2)

  polygons <- vector("list", columns * rows)
  index <- 1L
  for (row in seq_len(rows)) {
    ymin <- origin_y + (row - 1) * tile_size
    ymax <- ymin + tile_size
    for (column in seq_len(columns)) {
      xmin <- origin_x + (column - 1) * tile_size
      xmax <- xmin + tile_size
      ring <- matrix(c(
        xmin, ymin,
        xmax, ymin,
        xmax, ymax,
        xmin, ymax,
        xmin, ymin
      ), ncol = 2, byrow = TRUE)
      polygons[[index]] <- st_polygon(list(ring))
      index <- index + 1L
    }
  }
  st_sf(
    tile_id = seq_along(polygons),
    geometry = st_sfc(polygons, crs = st_crs(geometry))
  )
}

build_analysis_tiles <- function(footprint, tile_size, search_step, threads = 1L) {
  if (identical(footprint$source, "point_cloud_extent")) {
    return(build_covering_tiles(footprint$usable, tile_size))
  }
  build_optimized_tiles(footprint$usable, tile_size, search_step, threads)
}

compute_edge_flags <- function(tiles, tile_size) {
  if (nrow(tiles) == 0) return(logical())
  centroids <- st_coordinates(st_centroid(st_geometry(tiles)))
  index_x <- as.integer(round(centroids[, 1] / tile_size))
  index_y <- as.integer(round(centroids[, 2] / tile_size))
  keys <- paste(index_x, index_y, sep = "_")

  vapply(seq_len(nrow(tiles)), function(index) {
    neighbours <- expand.grid(dx = -1:1, dy = -1:1)
    neighbours <- neighbours[neighbours$dx != 0 | neighbours$dy != 0, ]
    neighbour_keys <- paste(
      index_x[[index]] + neighbours$dx,
      index_y[[index]] + neighbours$dy,
      sep = "_"
    )
    !all(neighbour_keys %in% keys)
  }, logical(1))
}

# Disk-backed DTM and CHM pipeline -----------------------------------------

# Process one DTM catalog chunk.  The optional order field is selected by the
# catalog caller and restored immediately after reading.
dtm_chunk <- function(chunk, dtm_resolution, ptd_resolution,
                      point_order_dimension = NULL) {
  lidR::set_lidr_threads(1L)
  las <- readLAS(chunk)
  if (is.empty(las)) return(empty_chunk_raster(chunk, dtm_resolution))
  # Treat a missing order key as a hard error before classification or TIN work.
  assert_loaded_dimensions(
    las,
    c("X", "Y", "Z", "buffer", point_order_dimension),
    "DTM"
  )
  # Match the record sequence seen by Julia's original LAZ-based workflow.
  las <- restore_original_point_order(las, point_order_dimension, "DTM")
  if (!las_has_nondegenerate_xy(las)) {
    message(
      "DTM chunk has fewer than three non-collinear input points; ",
      "leaving it empty for buffered neighbouring chunks"
    )
    return(empty_chunk_raster(chunk, dtm_resolution))
  }
  # Julia's DTM does not use return-number semantics. Being explicit avoids
  # lidR changing behaviour when ReturnNumber happens to be present.
  las <- classify_ground(
    las,
    ptd(res = ptd_resolution),
    last_returns = FALSE
  )
  surface <- classified_ground_surface(las, dtm_resolution)
  if (is.null(surface)) empty_chunk_raster(chunk, dtm_resolution) else surface
}

# Build the same floating-point Delaunay algorithm that lidR 4.3.2 uses when
# its optimized integer path rejects a LAS scale, offset, or coordinate. Calling
# the legacy primitives explicitly makes the retry independent of lidR's sampled
# fast-path precheck and preserves every source coordinate and header value.
legacy_tin <- function(..., extrapolate = knnidw(3, 1, 50)) {
  # Resolve the extrapolation algorithm once, matching lidR::tin() semantics.
  extrapolate <- lazyeval::uq(extrapolate)
  # Define the DTM plugin callback that receives one LAS chunk and raster grid.
  algorithm <- function(las, where) {
    # A numeric matrix forces tDelaunay() onto its legacy floating-point branch.
    points <- as.matrix(las@data)
    # tInterpolate() expects candidate raster locations as a numeric XY matrix.
    query <- as.matrix(where)
    # lidR 4.3.2 drops matrix dimensions inside tInterpolate() for one query.
    single_query <- nrow(query) == 1L
    # Duplicate that one location solely to retain matrix shape in the legacy API.
    interpolation_query <- if (single_query) rbind(query, query) else query
    # Build the legacy triangulation without changing LAS coordinates or headers.
    triangulation <- lidR:::tDelaunay(points, trim = 0)
    # Interpolate grid elevations with the configured lidR thread allowance.
    elevations <- lidR:::tInterpolate(
      triangulation,
      points,
      interpolation_query,
      threads = lidR::get_lidr_threads()
    )
    # Discard the duplicate result introduced only for the one-cell API guard.
    if (single_query) elevations <- elevations[[1L]]
    # Find grid cells outside the convex hull, which Delaunay leaves as NA.
    outside <- is.na(elevations)
    # Preserve lidR::tin()'s configured nearest-neighbour extrapolation behavior.
    if (any(outside)) {
      # Tell the extrapolation plugin that it runs in spatial-interpolation mode.
      lidR.context <- "spatial_interpolation"
      # Pass only unresolved XY cells to avoid repeating successful interpolation.
      outside_coordinates <- data.frame(
        X = where$X[outside],
        Y = where$Y[outside]
      )
      # Fill outside-hull cells with the same default used by lidR::tin().
      elevations[outside] <- extrapolate(las, outside_coordinates)
    }
    # Return one elevation per requested raster location to the DTM plugin.
    elevations
  }
  # Mark the callback as a lidR DTM algorithm with threaded interpolation support.
  lidR:::plugin_dtm(algorithm, omp = TRUE)
}

las_has_nondegenerate_xy <- function(las) {
  xy <- unique(as.data.table(las@data)[, .(X, Y)])
  if (nrow(xy) < 3L) return(FALSE)
  qr(scale(as.matrix(xy), center = TRUE, scale = FALSE))$rank >= 2L
}

chunk_has_valid_extent <- function(chunk) {
  chunk_bounds <- st_bbox(chunk)
  all(is.finite(chunk_bounds)) &&
    chunk_bounds[["xmax"]] > chunk_bounds[["xmin"]] &&
    chunk_bounds[["ymax"]] > chunk_bounds[["ymin"]]
}

ensure_nonempty_grid_bounds <- function(bounds, resolution) {
  if (bounds[["xmax"]] <= bounds[["xmin"]]) {
    bounds[["xmax"]] <- bounds[["xmin"]] + resolution
  }
  if (bounds[["ymax"]] <= bounds[["ymin"]]) {
    bounds[["ymax"]] <- bounds[["ymin"]] + resolution
  }
  bounds
}

empty_chunk_raster <- function(chunk, resolution) {
  if (!chunk_has_valid_extent(chunk)) return(NULL)
  chunk_bounds <- st_bbox(chunk)
  bounds <- ensure_nonempty_grid_bounds(
    covering_grid_bounds(chunk_bounds, resolution),
    resolution
  )
  raster <- terra::rast(
    xmin = bounds[["xmin"]],
    xmax = bounds[["xmax"]],
    ymin = bounds[["ymin"]],
    ymax = bounds[["ymax"]],
    resolution = resolution
  )
  terra::values(raster) <- NA_real_
  raster
}

classified_ground_surface <- function(las, dtm_resolution) {
  classifications <- las[["Classification"]]
  if (is.null(classifications) || !any(classifications == 2L, na.rm = TRUE)) {
    message(
      "DTM chunk contains no PTD ground points; leaving it empty for ",
      "buffered neighbouring chunks"
    )
    return(NULL)
  }
  ground_xy <- unique(
    as.data.table(las@data)[Classification == 2L, .(X, Y)]
  )
  ground_rank <- if (nrow(ground_xy) < 3L) {
    0L
  } else {
    qr(scale(as.matrix(ground_xy), center = TRUE, scale = FALSE))$rank
  }
  if (nrow(ground_xy) < 3L || ground_rank < 2L) {
    message(
      "DTM chunk has fewer than three non-collinear PTD ground points; ",
      "leaving it empty for buffered neighbouring chunks"
    )
    return(NULL)
  }
  tryCatch(
    rasterize_terrain(las, res = dtm_resolution, algorithm = tin()),
    error = function(error) {
      # Retry only the known lidR integer-conversion failure.  Unrelated TIN
      # errors must retain their existing behavior and diagnostic path.
      if (grepl(
        "xy coordinates were not converted to integer",
        conditionMessage(error),
        fixed = TRUE
      )) {
        message(
          paste(
            "Retrying DTM TIN with lidR's legacy floating-point implementation;",
            "coordinates and header are unchanged"
          )
        )
        # Repeat the same DTM resolution with the original LAS and slow TIN path.
        return(rasterize_terrain(
          las,
          res = dtm_resolution,
          algorithm = legacy_tin()
        ))
      }
      if (grepl(
        paste0(
          "cannot triangulate less than 3 points|No ground points found|",
          "Invalid grid dimensions"
        ),
        conditionMessage(error),
        ignore.case = TRUE
      )) {
        message(
          "DTM chunk has insufficient ground geometry for TIN; ",
          "leaving it empty for buffered neighbouring chunks"
        )
        return(NULL)
      }
      stop(error)
    }
  )
}

assert_catalog_completed <- function(results, phase) {
  if (length(results) == 0L || any(vapply(results, is.null, logical(1)))) {
    stop(sprintf(
      "%s LAScatalog processing failed or returned incomplete chunks",
      phase
    ))
  }
  invisible(results)
}

with_catalog_workers <- function(workers, phase, expression) {
  workers <- max(1L, as.integer(workers))
  previous_lidr_threads <- lidR::get_lidr_threads()
  on.exit(
    lidR::set_lidr_threads(previous_lidr_threads),
    add = TRUE
  )
  message(sprintf(
    "Processing %s LAScatalog chunks with %d worker(s)",
    phase,
    workers
  ))
  if (workers == 1L) return(force(expression))

  previous_plan <- future::plan()
  on.exit(future::plan(previous_plan), add = TRUE)
  future::plan(future::multisession, workers = workers)
  force(expression)
}

prepare_raster_work_directory <- function(path, label) {
  if (dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
  if (!dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop(sprintf("could not create temporary %s chunk directory", label))
  }
  invisible(path)
}

virtual_raster_source <- function(result, label) {
  if (!inherits(result, "SpatRaster") || terra::inMemory(result)) {
    stop(sprintf(
      "global %s processing did not produce a disk-backed raster mosaic",
      label
    ))
  }
  sources <- terra::sources(result)
  if (length(sources) != 1L ||
      !grepl("\\.vrt$", sources[[1]], ignore.case = TRUE)) {
    stop(sprintf(
      "global %s processing did not produce one virtual raster mosaic",
      label
    ))
  }
  sources[[1]]
}

complete_chunk_raster_paths <- function(results, label) {
  assert_catalog_completed(results, label)
  paths <- unname(unlist(results, use.names = FALSE))
  if (length(paths) == 0L || any(!file.exists(paths))) {
    stop(sprintf(
      "%s LAScatalog processing did not write every chunk raster",
      label
    ))
  }

  paths
}

populated_or_all_empty_chunk_raster_paths <- function(paths, label) {
  has_data <- vapply(paths, function(path) {
    count <- suppressWarnings(terra::global(terra::rast(path), "notNA"))
    is.finite(count[[1L, 1L]]) && count[[1L, 1L]] > 0
  }, logical(1))
  if (any(has_data)) {
    return(paths[has_data])
  }

  message(sprintf(
    "global %s processing produced only empty chunks; publishing an all-NoData raster",
    label
  ))
  paths
}

chunk_reference_crs <- function(paths, label) {
  candidates <- vapply(paths, function(path) {
    terra::crs(terra::rast(path))
  }, character(1))
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  if (length(candidates) > 1L) {
    stop(sprintf("%s chunk rasters contain conflicting CRS metadata", label))
  }
  if (length(candidates) == 1L) candidates[[1L]] else NA_character_
}

spatial_reference_wkt <- function(value) {
  reference <- sf::st_crs(value)
  if (is.na(reference) || is.null(reference$wkt) ||
      is.na(reference$wkt) || !nzchar(reference$wkt)) {
    return(NA_character_)
  }
  reference$wkt
}

usable_crs_wkt <- function(value) {
  length(value) == 1L && !is.null(value) && !is.na(value) && nzchar(value)
}

resolve_analysis_crs <- function(point_cloud_crs, original_companion_crs = NA_character_) {
  if (usable_crs_wkt(point_cloud_crs)) {
    return(list(wkt = point_cloud_crs, source = "point_cloud"))
  }
  if (usable_crs_wkt(original_companion_crs)) {
    return(list(wkt = original_companion_crs, source = "original_companion"))
  }
  list(wkt = NA_character_, source = "unavailable")
}

build_chunk_virtual_raster <- function(results, work_directory, label,
                                       fallback_crs = NA_character_) {
  all_paths <- complete_chunk_raster_paths(results, label)
  reference_crs <- chunk_reference_crs(all_paths, label)
  if ((is.na(reference_crs) || !nzchar(reference_crs)) &&
      !is.na(fallback_crs) && nzchar(fallback_crs)) {
    reference_crs <- fallback_crs
  }
  paths <- populated_or_all_empty_chunk_raster_paths(all_paths, label)
  vrt_path <- file.path(
    work_directory,
    sprintf("%s_chunks.vrt", tolower(label))
  )
  input_list_path <- file.path(
    work_directory,
    sprintf("%s_chunks.txt", tolower(label))
  )
  unlink(vrt_path, force = TRUE)
  writeLines(paths, input_list_path, useBytes = TRUE)
  arguments <- c("-overwrite")
  if (!is.na(reference_crs) && nzchar(reference_crs)) {
    arguments <- c(arguments, "-a_srs", shQuote(reference_crs))
  }
  arguments <- c(
    arguments,
    "-input_file_list", shQuote(input_list_path),
    shQuote(vrt_path)
  )
  status <- system2("gdalbuildvrt", arguments)
  if (!identical(status, 0L) || !file.exists(vrt_path)) {
    stop(sprintf("could not build the global %s chunk VRT", label))
  }
  result <- terra::rast(vrt_path)
  if (!is.na(reference_crs) && nzchar(reference_crs)) {
    actual_crs <- terra::crs(result)
    if (is.na(actual_crs) || !nzchar(actual_crs)) {
      stop(sprintf(
        "global %s chunk VRT did not persist its reference CRS",
        label
      ))
    }
  }
  result
}

gdal_raster_has_crs <- function(path) {
  info <- system2(
    "gdalinfo",
    c("-json", shQuote(path)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(info, "status")
  if (!is.null(status) && status != 0L) {
    stop(sprintf("could not inspect the global raster CRS for %s", path))
  }
  any(grepl('"coordinateSystem"', info, fixed = TRUE))
}

align_raster_file_crs <- function(path, source_crs, label) {
  if (is.null(source_crs)) return(invisible(path))

  source_has_crs <- !is.na(source_crs) && nzchar(source_crs)
  definition <- if (source_has_crs) source_crs else "None"
  status <- system2(
    "gdal_edit.py",
    c("-a_srs", shQuote(definition), shQuote(path))
  )
  if (!identical(status, 0L)) {
    stop(sprintf(
      "could not align the global %s raster CRS with the point-cloud source",
      label
    ))
  }

  output_has_crs <- gdal_raster_has_crs(path)
  if (!identical(output_has_crs, source_has_crs)) {
    stop(sprintf(
      "global %s raster CRS does not match the point-cloud source state",
      label
    ))
  }
  invisible(path)
}

stream_virtual_raster <- function(result, output_path, label,
                                  source_crs = NULL) {
  source_path <- virtual_raster_source(result, label)
  partial_path <- paste0(output_path, ".partial.tif")
  unlink(partial_path, force = TRUE)
  on.exit(unlink(partial_path, force = TRUE), add = TRUE)
  status <- system2(
    "gdal_translate",
    c(
      "-of", "GTiff",
      "-co", "TILED=YES",
      "-co", "COMPRESS=LZW",
      "-co", "BIGTIFF=IF_SAFER",
      shQuote(source_path),
      shQuote(partial_path)
    )
  )
  if (!identical(status, 0L) || !file.exists(partial_path)) {
    stop(sprintf(
      "could not stream the global %s virtual mosaic to GeoTIFF",
      label
    ))
  }
  if (!file.rename(partial_path, output_path)) {
    stop(sprintf("could not atomically publish the global %s GeoTIFF", label))
  }
  align_raster_file_crs(output_path, source_crs, label)
  invisible(output_path)
}

align_grid_coordinate <- function(value, resolution, direction) {
  scaled <- value / resolution
  nearest <- round(scaled)
  if (abs(scaled - nearest) < 1e-9) scaled <- nearest
  if (identical(direction, "lower")) {
    floor(scaled) * resolution
  } else {
    ceiling(scaled) * resolution
  }
}

covering_grid_bounds <- function(bounds, resolution) {
  c(
    xmin = align_grid_coordinate(bounds[["xmin"]], resolution, "lower"),
    ymin = align_grid_coordinate(bounds[["ymin"]], resolution, "lower"),
    xmax = align_grid_coordinate(bounds[["xmax"]], resolution, "upper"),
    ymax = align_grid_coordinate(bounds[["ymax"]], resolution, "upper")
  )
}

lidr_catalog_grid_bounds <- function(point_cloud, resolution) {
  required_version <- "4.3.2"
  actual_version <- as.character(utils::packageVersion("lidR"))
  if (!identical(actual_version, required_version)) {
    stop(sprintf(
      paste0(
        "Julia-faithful DTM layout requires pinned lidR %s; ",
        "found %s"
      ),
      required_version,
      actual_version
    ))
  }

  # Julia's original rasterize_terrain() call delegates numeric-resolution
  # geometry to this pinned lidR layout. In particular, it retains an extra
  # cell when a LAS bound lies exactly on the grid alignment. Reusing the same
  # implementation avoids approximating that half-cell rounding convention.
  layout <- lidR:::raster_layout(
    readLAScatalog(point_cloud),
    resolution,
    start = c(0, 0),
    buffer = 0,
    format = "template"
  )
  bounds <- c(
    xmin = as.numeric(layout$xmin),
    ymin = as.numeric(layout$ymin),
    xmax = as.numeric(layout$xmax),
    ymax = as.numeric(layout$ymax)
  )
  if (any(!is.finite(bounds)) ||
      bounds[["xmax"]] <= bounds[["xmin"]] ||
      bounds[["ymax"]] <= bounds[["ymin"]]) {
    stop("pinned lidR returned an invalid Julia-faithful DTM layout")
  }
  bounds
}

raster_extent_covers_point_bounds <- function(raster_extent, point_bounds,
                                               resolution) {
  extent_values <- c(
    xmin = as.numeric(raster_extent$xmin),
    ymin = as.numeric(raster_extent$ymin),
    xmax = as.numeric(raster_extent$xmax),
    ymax = as.numeric(raster_extent$ymax)
  )
  coordinate_scale <- max(abs(c(extent_values, point_bounds)), 1)
  tolerance <- max(
    resolution * 1e-9,
    .Machine$double.eps * coordinate_scale * 16
  )
  extent_values[["xmin"]] <= point_bounds[["xmin"]] + tolerance &&
    extent_values[["xmax"]] >= point_bounds[["xmax"]] - tolerance &&
    extent_values[["ymin"]] <= point_bounds[["ymin"]] + tolerance &&
    extent_values[["ymax"]] >= point_bounds[["ymax"]] - tolerance
}

cover_virtual_raster <- function(result, point_cloud, resolution,
                                 work_directory, label,
                                 extent_mode = c("point_bounds", "lidr")) {
  extent_mode <- match.arg(extent_mode)
  source_path <- virtual_raster_source(result, label)
  point_bounds <- point_cloud_xy_bounds(point_cloud)
  target_bounds <- if (identical(extent_mode, "lidr")) {
    lidr_catalog_grid_bounds(point_cloud, resolution)
  } else {
    covering_grid_bounds(point_bounds, resolution)
  }
  if (identical(extent_mode, "lidr") &&
      !raster_extent_covers_point_bounds(
        terra::ext(result),
        target_bounds,
        resolution
      )) {
    stop(sprintf(
      paste0(
        "global %s mosaic does not contain the complete pinned-lidR ",
        "publication extent"
      ),
      label
    ))
  }
  covering_path <- file.path(
    work_directory,
    sprintf("%s_covering.vrt", tolower(label))
  )
  number <- function(value) {
    format(value, scientific = FALSE, trim = TRUE, digits = 15)
  }
  status <- system2(
    "gdalbuildvrt",
    c(
      "-overwrite",
      "-resolution", "user",
      "-tr", number(resolution), number(resolution),
      "-te",
      number(target_bounds[["xmin"]]),
      number(target_bounds[["ymin"]]),
      number(target_bounds[["xmax"]]),
      number(target_bounds[["ymax"]]),
      "-vrtnodata", "-999999",
      shQuote(covering_path),
      shQuote(source_path)
    )
  )
  if (!identical(status, 0L) || !file.exists(covering_path)) {
    stop(sprintf(
      "could not create the point-cloud-covering global %s VRT",
      label
    ))
  }
  covering <- terra::rast(covering_path)
  if (terra::inMemory(covering)) {
    stop(sprintf("global %s covering VRT materialized in memory", label))
  }
  covering_extent <- terra::ext(covering)
  if (!raster_extent_covers_point_bounds(
    covering_extent,
    point_bounds,
    resolution
  )) {
    stop(sprintf(
      "global %s covering VRT does not contain the point-cloud XY extent",
      label
    ))
  }
  covering
}

# Build the disk-backed global DTM while carrying the COPC order key through
# lidR's selective catalog reads.
write_global_dtm <- function(point_cloud, output_path, work_directory,
                             chunk_size, dtm_buffer, dtm_resolution,
                             ptd_resolution, workers, source_crs = NULL,
                             point_order_dimension = NULL) {
  prepare_raster_work_directory(work_directory, "DTM")
  on.exit(unlink(work_directory, recursive = TRUE, force = TRUE), add = TRUE)

  catalog <- readLAScatalog(point_cloud)
  opt_chunk_size(catalog) <- chunk_size
  opt_chunk_buffer(catalog) <- dtm_buffer
  # Plain LAZ requests xyz; canonical COPC adds the OriginalPointIndex ordinal.
  opt_select(catalog) <- catalog_selection_for_dimensions(
    point_cloud,
    point_order_dimension,
    DTM_CATALOG_SELECTION
  )
  opt_progress(catalog) <- FALSE
  opt_output_files(catalog) <- file.path(work_directory, "dtm_{ID}")
  catalog@output_options$drivers$SpatRaster$param$datatype <- "FLT8S"
  if (is.null(source_crs)) source_crs <- spatial_reference_wkt(catalog)

  results <- with_catalog_workers(
    workers,
    "DTM",
    catalog_apply(
      catalog,
      dtm_chunk,
      dtm_resolution = dtm_resolution,
      ptd_resolution = ptd_resolution,
      # Every worker restores source record order before ground/TIN operations.
      point_order_dimension = point_order_dimension,
      .options = list(
        automerge = FALSE,
        drop_null = FALSE,
        raster_alignment = dtm_resolution
      )
    )
  )
  chunk_paths <- complete_chunk_raster_paths(results, "DTM")
  mosaic <- build_dtm_overlap_mosaic(
    chunk_paths,
    work_directory,
    fallback_crs = source_crs
  )
  covering <- cover_virtual_raster(
    mosaic,
    point_cloud,
    dtm_resolution,
    work_directory,
    "DTM",
    extent_mode = "lidr"
  )
  stream_virtual_raster(covering, output_path, "DTM", source_crs)
}

build_dtm_overlap_mosaic <- function(chunk_paths, work_directory,
                                     fallback_crs = NA_character_) {
  reference_crs <- chunk_reference_crs(chunk_paths, "DTM")
  if ((is.na(reference_crs) || !nzchar(reference_crs)) &&
      !is.na(fallback_crs) && nzchar(fallback_crs)) {
    reference_crs <- fallback_crs
  }
  mosaic_path <- file.path(work_directory, "dtm_mosaic.tif")
  mosaic <- terra::mosaic(
    terra::sprc(lapply(chunk_paths, terra::rast)),
    filename = mosaic_path,
    overwrite = TRUE,
    wopt = list(
      datatype = "FLT8S",
      gdal = c("TILED=YES", "COMPRESS=LZW", "BIGTIFF=IF_SAFER")
    )
  )
  if (terra::inMemory(mosaic) || !file.exists(mosaic_path)) {
    stop("global DTM overlap mosaic was not written disk-to-disk")
  }

  mosaic_vrt_path <- file.path(work_directory, "dtm_mosaic.vrt")
  arguments <- c("-overwrite")
  if (!is.na(reference_crs) && nzchar(reference_crs)) {
    arguments <- c(arguments, "-a_srs", shQuote(reference_crs))
  }
  arguments <- c(
    arguments,
    shQuote(mosaic_vrt_path),
    shQuote(mosaic_path)
  )
  status <- system2(
    "gdalbuildvrt",
    arguments
  )
  if (!identical(status, 0L) || !file.exists(mosaic_vrt_path)) {
    stop("could not create the global DTM publication VRT")
  }
  result <- terra::rast(mosaic_vrt_path)
  if (!is.na(reference_crs) && nzchar(reference_crs)) {
    actual_crs <- terra::crs(result)
    if (is.na(actual_crs) || !nzchar(actual_crs)) {
      stop("global DTM publication VRT did not persist its reference CRS")
    }
  }
  result
}

# Persist tiny per-window counters independently from scientific artifacts so
# worker processes can report observed point pressure without sharing memory.
record_chunk_observation <- function(directory, phase, point_count,
                                     target_chunk_points) {
  if (is.null(directory)) return(invisible(NULL))
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  oversized <- is.finite(target_chunk_points) && target_chunk_points > 0 &&
    point_count > target_chunk_points
  message(sprintf(
    "ADAPTIVE_CHUNK_OBSERVED phase=%s points=%.0f target_points=%.0f oversized=%s",
    phase,
    point_count,
    target_chunk_points,
    tolower(as.character(oversized))
  ))
  if (oversized) {
    message(sprintf(
      paste0(
        "FORESTSTRUCTURE_FLAG adaptive_chunk_oversize: %s chunk contained ",
        "%.0f points above the approximate %.0f-point target"
      ),
      phase,
      point_count,
      target_chunk_points
    ))
  }
  path <- tempfile(
    # Include the worker PID so the optional two-worker compatibility mode
    # cannot collide even if forked workers inherit identical RNG state.
    pattern = sprintf("%s_chunk_%d_", phase, Sys.getpid()),
    tmpdir = directory,
    fileext = ".rds"
  )
  saveRDS(
    list(point_count = point_count, oversized = oversized),
    path,
    compress = FALSE
  )
  invisible(path)
}

summarize_chunk_observations <- function(directory) {
  paths <- if (is.null(directory) || !dir.exists(directory)) {
    character()
  } else {
    list.files(directory, pattern = "\\.rds$", full.names = TRUE)
  }
  if (length(paths) == 0L) {
    return(list(count = 0L, minimum_points = NA_real_, maximum_points = NA_real_,
                oversized_count = 0L))
  }
  observations <- lapply(paths, readRDS)
  points <- vapply(observations, `[[`, numeric(1), "point_count")
  oversized <- vapply(observations, `[[`, logical(1), "oversized")
  list(
    count = length(paths),
    minimum_points = min(points),
    maximum_points = max(points),
    oversized_count = sum(oversized)
  )
}

combine_chunk_observations <- function(summaries) {
  summaries <- Filter(function(summary) summary$count > 0L, summaries)
  if (length(summaries) == 0L) return(summarize_chunk_observations(NULL))
  list(
    count = sum(vapply(summaries, `[[`, integer(1), "count")),
    minimum_points = min(vapply(summaries, `[[`, numeric(1), "minimum_points")),
    maximum_points = max(vapply(summaries, `[[`, numeric(1), "maximum_points")),
    oversized_count = sum(vapply(summaries, `[[`, integer(1), "oversized_count"))
  )
}

# Process one CHM chunk after restoring the source record sequence when COPC is
# active.  This retains Julia parity for order-dependent voxel/canopy ties.
chm_chunk <- function(chunk, dtm_path, voxel_resolution, maximum_height,
                      chm_resolution, point_order_dimension = NULL,
                      observation_directory = NULL,
                      target_chunk_points = 0) {
  lidR::set_lidr_threads(1L)
  las <- readLAS(chunk)
  record_chunk_observation(
    observation_directory,
    "chm",
    nrow(las@data),
    target_chunk_points
  )
  if (is.empty(las)) return(empty_chunk_raster(chunk, chm_resolution))
  # Sort before normalization, filtering, voxelization, or canopy rasterization.
  las <- restore_original_point_order(las, point_order_dimension, "CHM")
  las <- normalize_height(las, terra::rast(dtm_path))
  las <- filter_poi(las, is.finite(Z) & Z >= 0 & Z <= maximum_height)
  if (is.empty(las)) return(empty_chunk_raster(chunk, chm_resolution))
  voxels <- voxelize_points(las, res = voxel_resolution)
  if (is.empty(voxels)) return(empty_chunk_raster(chunk, chm_resolution))
  # A fixed core template prevents lidR from deriving a data-dependent local
  # raster extent.  Every window therefore owns exactly its aligned grid cells,
  # including sparse edge windows and points lying on chunk boundaries.
  chunk_template <- empty_chunk_raster(chunk, chm_resolution)
  template_crs <- spatial_reference_wkt(las)
  if (!is.na(template_crs) && nzchar(template_crs)) {
    terra::crs(chunk_template) <- template_crs
  }
  chunk_raster <- rasterize_canopy(
    voxels,
    res = chunk_template,
    algorithm = p2r()
  )
  # A worker needs only the raster after this point. Drop the raw points and
  # voxel table before LAScatalog schedules the next spatial window.
  rm(las, voxels, chunk_template)
  gc(verbose = FALSE)
  chunk_raster
}

# Build the disk-backed CHM and ensure catalog workers receive the order key.
write_global_chm <- function(point_cloud, dtm_path, output_path, work_directory,
                             chunk_size, voxel_resolution, maximum_height,
                             chm_resolution, workers, source_crs = NULL,
                             point_order_dimension = NULL,
                             target_chunk_points = 0) {
  prepare_raster_work_directory(work_directory, "CHM")
  on.exit(unlink(work_directory, recursive = TRUE, force = TRUE), add = TRUE)
  observation_directory <- file.path(work_directory, "observations")
  dir.create(observation_directory)

  catalog <- readLAScatalog(point_cloud)
  opt_chunk_size(catalog) <- chunk_size
  opt_chunk_alignment(catalog) <- c(0, 0)
  opt_chunk_buffer(catalog) <- 0
  # Keep the COPC read spatial/selective while adding only the required key.
  opt_select(catalog) <- catalog_selection_for_dimensions(
    point_cloud,
    point_order_dimension,
    "xyz"
  )
  opt_progress(catalog) <- FALSE
  opt_output_files(catalog) <- file.path(work_directory, "chm_{ID}")
  if (is.null(source_crs)) source_crs <- spatial_reference_wkt(catalog)

  results <- with_catalog_workers(
    workers,
    "CHM",
    catalog_apply(
      catalog,
      chm_chunk,
      dtm_path = dtm_path,
      voxel_resolution = voxel_resolution,
      maximum_height = maximum_height,
      chm_resolution = chm_resolution,
      # Each CHM worker sorts only its bounded chunk, never the complete cloud.
      point_order_dimension = point_order_dimension,
      observation_directory = observation_directory,
      target_chunk_points = target_chunk_points,
      .options = list(
        automerge = FALSE,
        drop_null = FALSE,
        raster_alignment = chm_resolution
      )
    )
  )
  result <- build_chunk_virtual_raster(
    results,
    work_directory,
    "CHM",
    source_crs
  )
  covering <- cover_virtual_raster(
    result,
    point_cloud,
    chm_resolution,
    work_directory,
    "CHM"
  )
  stream_virtual_raster(covering, output_path, "CHM", source_crs)
  summarize_chunk_observations(observation_directory)
}

# Memory-bounded global segment aggregation --------------------------------

available_extra_dimensions <- function(point_cloud) {
  header <- readLASheader(point_cloud)
  extra_bytes <- header@VLR$Extra_Bytes[["Extra Bytes Description"]]
  if (is.null(extra_bytes) || length(extra_bytes) == 0) return(character())
  unique(vapply(extra_bytes, function(description) {
    as.character(description$name)
  }, character(1)))
}

select_instance_dimensions <- function(point_cloud, candidates) {
  available <- available_extra_dimensions(point_cloud)
  unique(candidates[candidates %in% available])
}

segment_catalog_selection <- function(point_cloud, instance_dimensions,
                                      base_selection = SEGMENT_CATALOG_BASE_SELECTION) {
  header <- readLASheader(point_cloud)
  extra_bytes <- header@VLR$Extra_Bytes[["Extra Bytes Description"]]
  if (is.null(extra_bytes) || length(extra_bytes) == 0) {
    return(base_selection)
  }
  available <- vapply(extra_bytes, function(description) {
    as.character(description$name)
  }, character(1))
  segment_selection_from_names(available, instance_dimensions, base_selection)
}

# Convert named ExtraByte requirements into lidR's compact LAS selection syntax.
# lidR addresses positions 1-9 individually; selector 0 means all ExtraBytes.
catalog_selection_for_dimensions <- function(point_cloud, dimensions,
                                             base_selection = "xyz") {
  # NULL becomes character(0); remove absent placeholders passed by LAZ mode.
  dimensions <- dimensions[!is.na(dimensions) & nzchar(dimensions)]
  # Avoid reading any ExtraByte when the caller needs only standard XYZ fields.
  if (length(dimensions) == 0L) return(base_selection)
  # Inspect VLR definitions without loading point records.
  header <- readLASheader(point_cloud)
  extra_bytes <- header@VLR$Extra_Bytes[["Extra Bytes Description"]]
  if (is.null(extra_bytes) || length(extra_bytes) == 0L) {
    stop("requested catalog ExtraByte dimensions are absent")
  }
  # Preserve declared ordinal order because selectors are position-based.
  available <- vapply(extra_bytes, function(description) {
    as.character(description$name)
  }, character(1))
  # A missing required field would make later order restoration/science unsafe.
  missing <- setdiff(dimensions, available)
  if (length(missing) > 0L) {
    stop(sprintf(
      "catalog is missing required ExtraByte dimension(s): %s",
      paste(missing, collapse = ", ")
    ))
  }
  # Translate names to one-based lidR ExtraByte selector positions.
  positions <- unique(match(dimensions, available))
  # Position 0 requests all ExtraBytes when any requirement lies beyond 9.
  if (any(positions > 9L)) return(paste0(base_selection, "0"))
  # Concatenate individual positions for the minimal selective read.
  paste0(base_selection, paste(positions, collapse = ""))
}

segment_selection_from_names <- function(available, instance_dimensions,
                                         base_selection = SEGMENT_CATALOG_BASE_SELECTION) {
  positions <- unique(match(instance_dimensions, available))
  positions <- positions[!is.na(positions)]
  if (length(positions) == 0) return(base_selection)

  # lidR exposes individual ExtraBytes selectors only for positions 1-9.
  # Position 0 requests all ExtraBytes while still excluding unrelated
  # standard LAS attributes. This fallback preserves dimensions beyond the
  # ninth ExtraBytes field without returning to the all-attribute wildcard.
  if (any(positions > 9L)) {
    if (identical(base_selection, JULIA_MEMORY_SAFE_SEGMENT_SELECTION)) {
      stop(
        "Julia-memory-safe selection requires a verified streaming ",
        "projection for ExtraByte positions beyond nine"
      )
    }
    return(paste0(base_selection, "0"))
  }
  paste0(base_selection, paste(positions, collapse = ""))
}

segment_bucket <- function(instance_id, bucket_count) {
  as.integer((as.numeric(instance_id) %% bucket_count) + 1L)
}

write_segment_chunk_store <- function(voxels, layers, apex, work_directory,
                                      bucket_count) {
  voxels[, bucket := segment_bucket(instance_id, bucket_count)]
  layers[, bucket := segment_bucket(instance_id, bucket_count)]
  apex[, bucket := segment_bucket(instance_id, bucket_count)]
  paths <- character()
  voxel_columns <- c(
    "instance_dimension", "instance_id", "voxel_x", "voxel_y", "voxel_z",
    intersect("voxel_first_order", names(voxels))
  )
  apex_columns <- c(
    "instance_dimension", "instance_id", "apex_x", "apex_y", "apex_z",
    intersect("apex_order", names(apex))
  )
  for (bucket_index in sort(unique(voxels$bucket))) {
    path <- tempfile(
      pattern = sprintf("segment_%d_%03d_", Sys.getpid(), bucket_index),
      tmpdir = work_directory,
      fileext = ".rds"
    )
    saveRDS(
      list(
        voxels = voxels[bucket == bucket_index, ..voxel_columns],
        layers = layers[bucket == bucket_index, .(
          instance_dimension, instance_id, occupied_layer
        )],
        apex = apex[bucket == bucket_index, ..apex_columns]
      ),
      path,
      compress = FALSE
    )
    paths[[as.character(bucket_index)]] <- path
  }
  paths
}

# Read and aggregate one bounded segmentation chunk.  COPC order restoration is
# intentionally performed before instance filtering and apex tie-breaking.
segment_chunk <- function(chunk, dtm_path, voxel_resolution, maximum_height,
                          instance_dimensions, work_directory, bucket_count,
                          strict_dimensions = FALSE,
                          point_order_dimension = NULL,
                          observation_directory = NULL,
                          target_chunk_points = 0) {
  configure_catalog_worker()
  las <- readLAS(chunk)
  record_chunk_observation(
    observation_directory,
    "segment",
    nrow(las@data),
    target_chunk_points
  )
  if (is.empty(las)) {
    if (!chunk_has_valid_extent(chunk)) return(NULL)
    return(list(paths = NULL))
  }
  # Recover original LAZ record order inside this spatially streamed chunk.
  las <- restore_original_point_order(
    las,
    point_order_dimension,
    "segment"
  )
  if (isTRUE(strict_dimensions)) {
    assert_loaded_dimensions(
      las,
      c(
        "X", "Y", "Z", "buffer", instance_dimensions,
        # Strict COPC mode proves that selector construction loaded the key too.
        point_order_dimension
      ),
      "segment"
    )
  }
  selected <- instance_dimensions[instance_dimensions %in% names(las@data)]
  if (length(selected) == 0) return(list(paths = NULL))
  dtm <- terra::rast(dtm_path)
  las <- normalize_height(las, dtm)
  las_data <- as.data.table(las@data)
  dimension_paths <- lapply(selected, function(instance_dimension) {
    source <- data.table(
      X = las_data$X,
      Y = las_data$Y,
      Z = las_data$Z,
      instance_id = las_data[[instance_dimension]],
      source_order = if (is.null(point_order_dimension)) {
        seq_len(nrow(las_data))
      } else {
        las_data[[point_order_dimension]]
      }
    )
    source <- source[
      is.finite(instance_id) & instance_id > 0 &
        is.finite(Z) & Z >= 0 & Z <= maximum_height
    ]
    if (nrow(source) == 0) return(character())

    segment_keys <- "instance_id"
    apex <- source[
      source[, .I[which.max(Z)], by = segment_keys]$V1,
      .(
        instance_id,
        apex_x = X,
        apex_y = Y,
        apex_z = Z,
        apex_order = source_order
      )
    ]
    source[, `:=`(
      voxel_x = floor(X / voxel_resolution),
      voxel_y = floor(Y / voxel_resolution),
      voxel_z = floor(Z / voxel_resolution)
    )]

    voxels <- unique(source[, .(
      instance_id,
      voxel_x,
      voxel_y,
      voxel_z,
      voxel_first_order = source_order
    )], by = c("instance_id", "voxel_x", "voxel_y", "voxel_z"))
    layers <- unique(source[, .(
      instance_id,
      occupied_layer = floor(Z)
    )])
    voxels[, instance_dimension := instance_dimension]
    setcolorder(voxels, c("instance_dimension", setdiff(names(voxels), "instance_dimension")))
    layers[, instance_dimension := instance_dimension]
    setcolorder(layers, c("instance_dimension", setdiff(names(layers), "instance_dimension")))
    apex[, instance_dimension := instance_dimension]
    setcolorder(apex, c("instance_dimension", setdiff(names(apex), "instance_dimension")))

    # Plain LAS/LAZ retains the legacy compact-summary schema and reducer order.
    # Only indexed COPC needs global source ordinals to neutralize spatial
    # reordering and adaptive window-boundary changes.
    if (is.null(point_order_dimension)) {
      voxels[, voxel_first_order := NULL]
      apex[, apex_order := NULL]
    }

    paths <- write_segment_chunk_store(
      voxels,
      layers,
      apex,
      work_directory,
      bucket_count
    )
    rm(source, voxels, layers, apex)
    gc(verbose = FALSE)
    paths
  })
  paths <- do.call(c, dimension_paths)
  # The compact RDS summaries now own everything needed by the disk reducer.
  # Release the normalized LAS and its table view before another window loads.
  rm(las, las_data, dtm, dimension_paths)
  gc(verbose = FALSE)
  list(paths = paths)
}

# Configure a spatially streamed segment catalog.  The order key travels beside
# requested instance dimensions but is not included in scientific aggregation.
accumulate_segments <- function(point_cloud, dtm_path, chunk_size,
                                voxel_resolution, maximum_height,
                                instance_dimensions, workers, work_directory,
                                bucket_count = SEGMENT_BUCKET_COUNT,
                                base_selection = SEGMENT_CATALOG_BASE_SELECTION,
                                point_order_dimension = NULL,
                                target_chunk_points = 0) {
  prepare_raster_work_directory(work_directory, "segment")
  observation_directory <- file.path(work_directory, "observations")
  dir.create(observation_directory)
  catalog <- readLAScatalog(point_cloud)
  opt_chunk_size(catalog) <- chunk_size
  opt_chunk_alignment(catalog) <- c(0, 0)
  opt_chunk_buffer(catalog) <- 0
  # Request both the segmentation values and source-row order in one COPC scan.
  selection <- segment_catalog_selection(
    point_cloud,
    c(instance_dimensions, point_order_dimension),
    base_selection
  )
  opt_select(catalog) <- selection
  opt_progress(catalog) <- FALSE
  message(sprintf(
    "Reading segment dimensions with LAS selection %s",
    selection
  ))
  chunks <- with_catalog_workers(
    workers,
    "segment",
    catalog_apply(
      catalog,
      segment_chunk,
      dtm_path = dtm_path,
      voxel_resolution = voxel_resolution,
      maximum_height = maximum_height,
      instance_dimensions = instance_dimensions,
      work_directory = work_directory,
      bucket_count = bucket_count,
      strict_dimensions = identical(
        base_selection,
        JULIA_MEMORY_SAFE_SEGMENT_SELECTION
      ),
      # The worker uses this only to sort its in-memory bounded chunk.
      point_order_dimension = point_order_dimension,
      observation_directory = observation_directory,
      target_chunk_points = target_chunk_points,
      .options = list(automerge = FALSE, drop_null = FALSE)
    )
  )
  assert_catalog_completed(chunks, "segment")
  chunks <- lapply(chunks, `[[`, "paths")
  chunks <- Filter(function(paths) length(paths) > 0L, chunks)
  observations <- summarize_chunk_observations(observation_directory)
  if (length(chunks) == 0) {
    return(list(chunks = list(), bucket_count = bucket_count,
                observations = observations))
  }
  list(
    chunks = chunks,
    bucket_count = bucket_count,
    observations = observations
  )
}

finalize_segments <- function(accumulated, parameters) {
  voxels <- accumulated$voxels
  layers <- accumulated$layers
  apex <- accumulated$apex

  segment_keys <- c("instance_dimension", "instance_id")
  voxel_metrics <- voxels[, .(n_vox = .N), by = segment_keys]
  voxel_metrics[, voxel_volume := n_vox * parameters$voxel_resolution^3]
  crown_metrics <- voxels[, .(
    crown_area = uniqueN(data.table(voxel_x, voxel_y)) *
      parameters$voxel_resolution^2
  ), by = segment_keys]
  layer_metrics <- layers[, .(
    n_occupied_layers = uniqueN(occupied_layer)
  ), by = segment_keys]
  extent_metrics <- voxels[, {
    if (.N < 3) {
      .(pca_extent_1 = NA_real_, pca_extent_2 = NA_real_, pca_extent_3 = NA_real_)
    } else {
      assert_memory_budget(.N * 3 * 8 * 6, parameters, "segment PCA", TRUE)
      coordinates <- as.matrix(.SD) * parameters$voxel_resolution
      extents <- tryCatch({
        rotated <- prcomp(coordinates, center = TRUE)$x
        apply(rotated, 2, function(values) diff(range(values)))
      }, error = function(error) rep(NA_real_, 3))
      .(
        pca_extent_1 = extents[[1]],
        pca_extent_2 = extents[[2]],
        pca_extent_3 = extents[[3]]
      )
    }
  }, by = segment_keys, .SDcols = c("voxel_x", "voxel_y", "voxel_z")]

  segments <- Reduce(
    function(left, right) merge(left, right, by = segment_keys, all = TRUE),
    list(apex, voxel_metrics, crown_metrics, layer_metrics, extent_metrics)
  )
  segments[, `:=`(
    pass_voxels = n_vox >= parameters$minimum_tree_voxels,
    pass_apex = apex_z > parameters$apex_minimum_height,
    pass_thickness = !is.na(pca_extent_3) &
      pca_extent_3 >= parameters$minimum_tree_thickness,
    pass_occupied_layers = n_occupied_layers >= parameters$minimum_occupied_layers
  )]
  segments[, is_tree := pass_voxels & pass_apex &
    pass_thickness & pass_occupied_layers]
  segments[, fail_reason := paste0(
    ifelse(!pass_voxels, "voxels ", ""),
    ifelse(!pass_apex, "apex ", ""),
    ifelse(!pass_thickness, "thickness ", ""),
    ifelse(!pass_occupied_layers, "layers ", "")
  )]
  segments
}

# Merge disk pieces in source-order semantics rather than worker completion
# order.  Internal ordinal fields are removed before scientific output schemas
# are assembled, so callers see the historical columns unchanged.
merge_segment_pieces <- function(pieces) {
  voxels <- rbindlist(lapply(pieces, function(piece) piece[["voxels"]]),
                      use.names = TRUE, fill = TRUE)
  voxel_keys <- c(
    "instance_dimension", "instance_id", "voxel_x", "voxel_y", "voxel_z"
  )
  if ("voxel_first_order" %in% names(voxels)) {
    setorderv(voxels, c(voxel_keys, "voxel_first_order"))
    voxels <- unique(voxels, by = voxel_keys)
    setorderv(voxels, c("instance_dimension", "instance_id", "voxel_first_order"))
  } else {
    voxels <- unique(voxels, by = voxel_keys)
  }

  layers <- unique(rbindlist(
    lapply(pieces, function(piece) piece[["layers"]]),
    use.names = TRUE,
    fill = TRUE
  ))
  apex <- rbindlist(lapply(pieces, function(piece) piece[["apex"]]),
                    use.names = TRUE, fill = TRUE)
  apex_keys <- c("instance_dimension", "instance_id")
  if ("apex_order" %in% names(apex)) {
    setorderv(apex, c(apex_keys, "apex_z", "apex_order"), c(1L, 1L, -1L, 1L))
    apex <- unique(apex, by = apex_keys)
    apex[, apex_order := NULL]
  } else {
    apex <- apex[
      apex[, .I[which.max(apex_z)], by = apex_keys]$V1
    ]
  }
  list(voxels = voxels, layers = layers, apex = apex)
}

segment_partition_estimate <- function(paths) {
  sizes <- file.info(paths)$size
  sizes[!is.finite(sizes)] <- 0
  # Uncompressed RDS expands into list/data.table objects and is copied by
  # rbindlist(), unique(), and PCA preparation. Splitting more often is safe.
  sum(sizes) * 8
}

load_segment_partition <- function(paths, parameters) {
  assert_memory_budget(
    segment_partition_estimate(paths),
    parameters,
    "segment partition aggregation"
  )
  pieces <- lapply(paths, readRDS)
  merged <- merge_segment_pieces(pieces)
  voxels <- merged$voxels
  layers <- merged$layers
  apex <- merged$apex
  result <- finalize_segments(
    list(voxels = voxels, layers = layers, apex = apex),
    parameters
  )
  rm(pieces, voxels, layers, apex)
  gc(verbose = FALSE)
  result
}

split_segment_partition <- function(paths, work_directory, depth, parameters) {
  children <- list(character(), character())
  keys <- character()
  divisor <- SEGMENT_BUCKET_COUNT * 2^depth
  for (path in paths) {
    assert_memory_budget(
      file.info(path)$size * 8,
      parameters,
      "one segment chunk summary"
    )
    piece <- readRDS(path)
    piece_keys <- unique(rbindlist(lapply(piece, function(values) {
      values[, .(instance_dimension, instance_id)]
    }), use.names = TRUE))
    keys <- unique(c(
      keys,
      paste(piece_keys$instance_dimension, piece_keys$instance_id)
    ))
    for (branch in 0:1) {
      subset_piece <- lapply(piece, function(values) {
        values[(floor(abs(as.numeric(instance_id)) / divisor) %% 2) == branch]
      })
      if (!any(vapply(subset_piece, nrow, integer(1)) > 0L)) next
      child_path <- tempfile(
        pattern = sprintf("segment_split_%02d_%d_", depth, branch),
        tmpdir = work_directory,
        fileext = ".rds"
      )
      saveRDS(subset_piece, child_path, compress = FALSE)
      children[[branch + 1L]] <- c(children[[branch + 1L]], child_path)
    }
    rm(piece, piece_keys)
    gc(verbose = FALSE)
  }
  list(paths = Filter(length, children), unique_keys = length(keys))
}

finalize_segment_partition <- function(paths, parameters, depth = 0L) {
  budget_gib <- parameter_or(parameters, "memory_budget_gib", Inf)
  estimate <- segment_partition_estimate(paths)
  if (!is.finite(budget_gib) ||
      current_rss_bytes() + estimate <= budget_gib * 1024^3) {
    return(load_segment_partition(paths, parameters))
  }

  split <- split_segment_partition(
    paths,
    dirname(paths[[1]]),
    depth,
    parameters
  )
  child_paths <- unlist(split$paths, use.names = FALSE)
  on.exit(unlink(child_paths, force = TRUE), add = TRUE)
  if (split$unique_keys <= 1L || depth >= 52L) {
    assert_memory_budget(
      estimate,
      parameters,
      "one complete segment",
      single_segment = TRUE
    )
  }
  message(sprintf(
    "Splitting oversized segment partition at depth %d into %d child partition(s)",
    depth + 1L,
    length(split$paths)
  ))
  results <- lapply(split$paths, function(child) {
    finalize_segment_partition(child, parameters, depth + 1L)
  })
  rbindlist(results, use.names = TRUE)
}

finalize_segment_store_in_memory <- function(store, parameters) {
  paths <- unlist(store$chunks, use.names = FALSE)
  if (length(paths) == 0L) return(NULL)
  pieces <- lapply(paths, readRDS)
  merged <- merge_segment_pieces(pieces)
  voxels <- merged$voxels
  layers <- merged$layers
  apex <- merged$apex
  segments <- finalize_segments(
    list(voxels = voxels, layers = layers, apex = apex),
    parameters
  )
  setkeyv(segments, c("instance_dimension", "instance_id"))
  segments
}

finalize_segment_store <- function(store, parameters,
                                   adapter = c("disk", "memory")) {
  adapter <- match.arg(adapter)
  if (identical(adapter, "memory")) {
    return(finalize_segment_store_in_memory(store, parameters))
  }
  previous_threads <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(previous_threads), add = TRUE)
  data.table::setDTthreads(DATA_TABLE_THREADS)
  finalized <- vector("list", store$bucket_count)
  for (bucket_index in seq_len(store$bucket_count)) {
    key <- as.character(bucket_index)
    paths <- unlist(lapply(store$chunks, function(chunk_paths) {
      unname(chunk_paths[names(chunk_paths) == key])
    }), use.names = FALSE)
    if (length(paths) == 0) next

    finalized[[bucket_index]] <- finalize_segment_partition(paths, parameters)
    gc(verbose = FALSE)
  }
  finalized <- Filter(Negate(is.null), finalized)
  if (length(finalized) == 0) return(NULL)
  segments <- rbindlist(finalized, use.names = TRUE)
  setkeyv(segments, c("instance_dimension", "instance_id"))
  segments
}

# Scientific tile metrics and output adapters ------------------------------

assign_segments_to_tiles <- function(segments, tiles) {
  if (is.null(segments)) return(NULL)
  segments[, tile_id := NA_integer_]
  if (nrow(tiles) == 0) return(segments)

  bounds <- t(vapply(seq_len(nrow(tiles)), function(index) {
    as.numeric(st_bbox(tiles[index, ]))
  }, numeric(4)))
  colnames(bounds) <- c("xmin", "ymin", "xmax", "ymax")
  segments[, tile_id := mapply(function(apex_x, apex_y) {
    matches <- which(
      apex_x >= bounds[, "xmin"] & apex_x < bounds[, "xmax"] &
        apex_y >= bounds[, "ymin"] & apex_y < bounds[, "ymax"]
    )
    if (length(matches) == 0) NA_integer_ else tiles$tile_id[[matches[[1]]]]
  }, apex_x, apex_y)]
  segments
}

gini <- function(values) {
  values <- values[!is.na(values)]
  count <- length(values)
  if (count < 2) return(NA_real_)
  if (any(values < 0)) values <- values - min(values)
  values <- sort(values)
  total <- sum(values)
  if (total == 0) return(0)
  coefficient <- (2 * sum(seq_len(count) * values) / (count * total)) -
    (count + 1) / count
  coefficient * count / (count - 1)
}

empty_tree_metrics <- function() {
  list(
    n_seg_total = NA_integer_,
    n_trees = NA_integer_,
    tree_height_max = NA_real_,
    tree_height_mean = NA_real_,
    tree_height_gini = NA_real_,
    tree_crownarea_mean = NA_real_,
    tree_crownarea_max = NA_real_,
    tree_crownarea_gini = NA_real_,
    tree_volume_mean = NA_real_,
    tree_volume_max = NA_real_,
    tree_volume_gini = NA_real_
  )
}

tree_metrics_for_tile <- function(segments, tile_id) {
  if (is.null(segments)) return(empty_tree_metrics())
  target_tile_id <- tile_id
  in_tile <- segments[segments$tile_id == target_tile_id]
  trees <- in_tile[is_tree == TRUE]
  metrics <- empty_tree_metrics()
  metrics$n_seg_total <- nrow(in_tile)
  metrics$n_trees <- nrow(trees)
  if (nrow(trees) == 0) return(metrics)

  metrics$tree_height_max <- safe_round(max(trees$apex_z))
  metrics$tree_height_mean <- safe_round(mean(trees$apex_z))
  metrics$tree_height_gini <- safe_round(gini(trees$apex_z))
  metrics$tree_crownarea_mean <- safe_round(mean(trees$crown_area))
  metrics$tree_crownarea_max <- safe_round(max(trees$crown_area))
  metrics$tree_crownarea_gini <- safe_round(gini(trees$crown_area))
  metrics$tree_volume_mean <- safe_round(mean(trees$voxel_volume))
  metrics$tree_volume_max <- safe_round(max(trees$voxel_volume))
  metrics$tree_volume_gini <- safe_round(gini(trees$voxel_volume))
  metrics
}

result_with_tree_metrics <- function(base_result, segments) {
  result <- copy(base_result)
  tree_columns <- names(empty_tree_metrics())
  if (nrow(result) == 0) return(result)
  for (index in seq_len(nrow(result))) {
    metrics <- tree_metrics_for_tile(segments, result$tile_id[[index]])
    for (name in tree_columns) {
      set(result, i = index, j = name, value = metrics[[name]])
    }
  }
  result
}

empty_segment_diagnostics <- function() {
  data.table(
    file = character(),
    sensor = character(),
    country = character(),
    tile_id = integer(),
    PredInstance = numeric(),
    n_vox = integer(),
    apex_z = numeric(),
    pc_ext1 = numeric(),
    pc_ext2 = numeric(),
    pc_ext3 = numeric(),
    n_zlayer = integer(),
    pass_vox = logical(),
    pass_apex = logical(),
    pass_thick = logical(),
    pass_zlayer = logical(),
    is_tree = logical(),
    apex_in_tile = logical()
  )
}

write_upstream_csv <- function(values, path) {
  utils::write.csv(
    as.data.frame(values),
    path,
    row.names = FALSE,
    na = "NA"
  )
  invisible(path)
}

write_segment_diagnostics <- function(segments, point_cloud, path,
                                      sensor = NA_character_,
                                      country = NA_character_) {
  diagnostics <- empty_segment_diagnostics()
  if (!is.null(segments) && nrow(segments) > 0) {
    diagnostics <- segments[, .(
      file = basename(point_cloud),
      sensor = sensor,
      country = country,
      tile_id,
      PredInstance = instance_id,
      n_vox,
      apex_z = round(apex_z, 4),
      pc_ext1 = round(pca_extent_1, 4),
      pc_ext2 = round(pca_extent_2, 4),
      pc_ext3 = round(pca_extent_3, 4),
      n_zlayer = n_occupied_layers,
      pass_vox = pass_voxels,
      pass_apex,
      pass_thick = pass_thickness,
      pass_zlayer = pass_occupied_layers,
      is_tree,
      apex_in_tile = !is.na(tile_id)
    )]
  }
  write_upstream_csv(diagnostics, path)
}

empty_metric_values <- function(voxel_total) {
  list(
    vox_filled = 0L,
    vox_total = voxel_total,
    veg_density = NA_real_,
    zsd = NA_real_,
    zskew = NA_real_,
    zkurt = NA_real_,
    zq90 = NA_real_,
    box_dim_fixed = NA_real_,
    vci = NA_real_,
    rumple = NA_real_,
    gap_fraction = NA_real_,
    chm_sd = NA_real_,
    chm_cv = NA_real_,
    height_max = NA_real_,
    height_mean = NA_real_
  )
}

safe_round <- function(value, digits = 4) {
  if (length(value) == 0 || !is.finite(value[[1]])) return(NA_real_)
  round(as.numeric(value[[1]]), digits)
}

# Calculate Julia-compatible tile metrics from a spatial window.  COPC still
# supplies only the requested window; order restoration occurs after that read.
calculate_tile_metrics <- function(point_cloud, tile, dtm, parameters) {
  bounds <- st_bbox(tile)
  voxel_total <- round(
    (parameters$tile_size / parameters$voxel_resolution)^2 *
      ((parameters$maximum_height - parameters$vegetation_minimum_height) /
         parameters$voxel_resolution)
  )
  empty <- empty_metric_values(voxel_total)

  las <- readLAS(
    point_cloud,
    # Add OriginalPointIndex to XYZ only when the active input is canonical COPC.
    select = catalog_selection_for_dimensions(
      point_cloud,
      parameters$point_order_dimension,
      "xyz"
    ),
    filter = sprintf(
      "-inside %.10f %.10f %.10f %.10f",
      bounds[["xmin"]], bounds[["ymin"]], bounds[["xmax"]], bounds[["ymax"]]
    )
  )
  if (is.empty(las)) return(empty)
  # Fail closed if selective COPC reading omitted the required order key.
  assert_loaded_dimensions(
    las,
    c("X", "Y", "Z", parameters$point_order_dimension),
    "tile"
  )
  # Match the point sequence used by the original Julia tile calculations.
  las <- restore_original_point_order(
    las,
    parameters$point_order_dimension,
    "tile"
  )

  las <- normalize_height(las, dtm)
  las <- filter_poi(las, is.finite(Z) & Z >= 0 & Z <= parameters$maximum_height)
  if (is.empty(las)) return(empty)

  voxels <- voxelize_points(las, res = parameters$voxel_resolution)
  if (is.empty(voxels)) return(empty)

  height_max <- max(voxels@data$Z, na.rm = TRUE)
  height_mean <- mean(voxels@data$Z, na.rm = TRUE)
  chm <- rasterize_canopy(
    voxels,
    res = parameters$chm_resolution,
    algorithm = p2r()
  )
  chm_values <- terra::values(chm, mat = FALSE)
  observed <- chm_values[is.finite(chm_values)]
  chm_mean <- if (length(observed)) mean(observed) else NA_real_
  chm_sd <- if (length(observed) > 1) stats::sd(observed) else NA_real_
  gap_fraction <- if (length(observed)) {
    mean(observed < parameters$gap_height_threshold)
  } else {
    NA_real_
  }
  rumple <- tryCatch(rumple_index(chm), error = function(error) NA_real_)

  vegetation <- filter_poi(voxels, Z > parameters$vegetation_minimum_height)
  if (is.empty(vegetation)) {
    empty$height_max <- safe_round(height_max)
    empty$height_mean <- safe_round(height_mean)
    empty$rumple <- safe_round(rumple)
    empty$gap_fraction <- safe_round(gap_fraction)
    empty$chm_sd <- safe_round(chm_sd)
    empty$chm_cv <- safe_round(chm_sd / chm_mean)
    return(empty)
  }

  vegetation_data <- as.data.table(vegetation@data)
  filled <- nrow(vegetation_data)
  z_metrics <- cloud_metrics(vegetation, .stdmetrics_z)

  voxel_indices <- unique(vegetation_data[, .(
    vx = floor(X / parameters$voxel_resolution),
    vy = floor(Y / parameters$voxel_resolution),
    vz = floor(Z / parameters$voxel_resolution)
  )])
  box_widths <- parameters$voxel_resolution * 2^(0:30)
  box_widths <- box_widths[box_widths < parameters$tile_size]
  scale_factors <- unique(round(box_widths / parameters$voxel_resolution))
  counts <- vapply(scale_factors, function(scale_factor) {
    uniqueN(voxel_indices[, .(
      vx = vx %/% scale_factor,
      vy = vy %/% scale_factor,
      vz = vz %/% scale_factor
    )])
  }, numeric(1))
  usable_counts <- is.finite(counts) & counts > 0
  box_dimension <- if (sum(usable_counts) >= 2) {
    unname(coef(lm(
      log(counts[usable_counts]) ~ log(1 / box_widths[usable_counts])
    ))[[2]])
  } else {
    NA_real_
  }

  layer_breaks <- seq(0, parameters$maximum_height, by = 1)
  layers <- cut(vegetation_data$Z, breaks = layer_breaks, include.lowest = TRUE)
  layer_counts <- table(factor(layers, levels = levels(layers)))
  proportions <- layer_counts / sum(layer_counts)
  nonzero <- proportions[proportions > 0]
  vci <- if (length(layer_counts) > 1) {
    -sum(nonzero * log(nonzero)) / log(length(layer_counts))
  } else {
    NA_real_
  }

  list(
    vox_filled = filled,
    vox_total = voxel_total,
    veg_density = safe_round(filled / voxel_total),
    zsd = safe_round(z_metrics$zsd),
    zskew = safe_round(z_metrics$zskew),
    zkurt = safe_round(z_metrics$zkurt),
    zq90 = safe_round(z_metrics$zq90),
    box_dim_fixed = safe_round(box_dimension),
    vci = safe_round(vci),
    rumple = safe_round(rumple),
    gap_fraction = safe_round(gap_fraction),
    chm_sd = safe_round(chm_sd),
    chm_cv = safe_round(chm_sd / chm_mean),
    height_max = safe_round(height_max),
    height_mean = safe_round(height_mean)
  )
}

build_tile_row <- function(point_cloud_name, tile, edge_tile, tile_count,
                           metrics, tree_metrics, sensor = NA_character_,
                           country = NA_character_) {
  bounds <- st_bbox(tile)
  as.data.table(c(list(
    file = point_cloud_name,
    sensor = sensor,
    country = country,
    tile_id = tile$tile_id[[1]],
    tile_xmin = as.numeric(bounds[["xmin"]]),
    tile_ymin = as.numeric(bounds[["ymin"]]),
    edge_tile = edge_tile,
    n_tiles_plot = as.integer(tile_count)
  ), metrics, tree_metrics))
}

empty_result_table <- function() {
  data.table(
    file = character(),
    sensor = character(),
    country = character(),
    tile_id = integer(),
    tile_xmin = numeric(),
    tile_ymin = numeric(),
    edge_tile = logical(),
    n_tiles_plot = integer(),
    vox_filled = integer(),
    vox_total = integer(),
    veg_density = numeric(),
    zsd = numeric(),
    zskew = numeric(),
    zkurt = numeric(),
    zq90 = numeric(),
    box_dim_fixed = numeric(),
    vci = numeric(),
    rumple = numeric(),
    gap_fraction = numeric(),
    chm_sd = numeric(),
    chm_cv = numeric(),
    height_max = numeric(),
    height_mean = numeric(),
    n_seg_total = integer(),
    n_trees = integer(),
    tree_height_max = numeric(),
    tree_height_mean = numeric(),
    tree_height_gini = numeric(),
    tree_crownarea_mean = numeric(),
    tree_crownarea_max = numeric(),
    tree_crownarea_gini = numeric(),
    tree_volume_mean = numeric(),
    tree_volume_max = numeric(),
    tree_volume_gini = numeric()
  )
}

apply_output_crs <- function(geometry, source_crs) {
  if (is.null(source_crs) || is.na(source_crs) || !nzchar(source_crs)) {
    st_crs(geometry) <- NA
  } else {
    st_crs(geometry) <- st_crs(source_crs)
  }
  geometry
}

write_tile_geojson <- function(tiles, result, source_crs, path) {
  if (nrow(tiles) == 0) {
    writeLines('{"type":"FeatureCollection","features":[]}', path, useBytes = TRUE)
    return(invisible(path))
  }
  output <- apply_output_crs(tiles, source_crs)
  attributes <- as.data.frame(result)
  for (name in setdiff(names(attributes), "tile_id")) {
    output[[name]] <- attributes[[name]]
  }
  st_write(output, path, driver = "GeoJSON", delete_dsn = TRUE, quiet = TRUE)
  invisible(path)
}

write_layout_png <- function(aoi, tiles, edge_flags, point_cloud_name, tile_size, path) {
  png(path, width = 1600, height = 1200, res = 150)
  on.exit(dev.off(), add = TRUE)
  par(mar = c(5, 5, 5, 2) + 0.1)
  plot_bounds <- st_bbox(aoi$inclusion)
  if (nrow(tiles) > 0) {
    tile_bounds <- st_bbox(tiles)
    plot_bounds[c("xmin", "ymin")] <- pmin(
      plot_bounds[c("xmin", "ymin")],
      tile_bounds[c("xmin", "ymin")]
    )
    plot_bounds[c("xmax", "ymax")] <- pmax(
      plot_bounds[c("xmax", "ymax")],
      tile_bounds[c("xmax", "ymax")]
    )
  }
  plot(
    aoi$inclusion,
    col = adjustcolor("#2e8b57", alpha.f = 0.12),
    border = "#17633b",
    lwd = 3,
    axes = TRUE,
    asp = 1,
    xlim = as.numeric(plot_bounds[c("xmin", "xmax")]),
    ylim = as.numeric(plot_bounds[c("ymin", "ymax")]),
    xlab = "Point-cloud-local X",
    ylab = "Point-cloud-local Y",
    main = sprintf("%s — %d valid Analysis Tiles", point_cloud_name, nrow(tiles))
  )
  if (length(aoi$exclusion) > 0) {
    plot(
      aoi$exclusion,
      add = TRUE,
      col = adjustcolor("#dc2626", alpha.f = 0.25),
      border = "#991b1b",
      lwd = 2
    )
  }
  if (nrow(tiles) > 0) {
    tile_colors <- ifelse(edge_flags, "#f59e0b", "#2563eb")
    for (index in seq_len(nrow(tiles))) {
      plot(st_geometry(tiles[index, ]), add = TRUE, border = tile_colors[[index]], lwd = 2)
    }
    centers <- st_coordinates(st_centroid(st_geometry(tiles)))
    text(centers[, 1], centers[, 2], labels = tiles$tile_id, cex = 0.85)
  } else {
    text(mean(par("usr")[1:2]), mean(par("usr")[3:4]), "0 valid tiles", cex = 1.4)
  }

  bounds <- st_bbox(aoi$inclusion)
  span_x <- as.numeric(bounds[["xmax"]] - bounds[["xmin"]])
  span_y <- as.numeric(bounds[["ymax"]] - bounds[["ymin"]])
  arrow_x <- as.numeric(bounds[["xmax"]] - span_x * 0.08)
  arrow_y <- as.numeric(bounds[["ymax"]] - span_y * 0.18)
  arrows(arrow_x, arrow_y, arrow_x, arrow_y + span_y * 0.1, length = 0.12, lwd = 2)
  text(arrow_x, arrow_y + span_y * 0.125, "N", font = 2)

  scale_x <- as.numeric(bounds[["xmin"]] + span_x * 0.06)
  scale_y <- as.numeric(bounds[["ymin"]] + span_y * 0.06)
  segments(scale_x, scale_y, scale_x + tile_size, scale_y, lwd = 4)
  segments(scale_x, scale_y - span_y * 0.01, scale_x, scale_y + span_y * 0.01, lwd = 2)
  segments(
    scale_x + tile_size,
    scale_y - span_y * 0.01,
    scale_x + tile_size,
    scale_y + span_y * 0.01,
    lwd = 2
  )
  text(scale_x + tile_size / 2, scale_y + span_y * 0.035, sprintf("%g m", tile_size))
  boundary_label <- if (identical(aoi$source, "point_cloud_extent")) {
    "Point-cloud extent"
  } else {
    "Audit AOI"
  }
  legend_labels <- c(boundary_label, "Edge tile", "Interior tile")
  legend_colors <- c("#17633b", "#f59e0b", "#2563eb")
  legend_widths <- c(3, 2, 2)
  if (length(aoi$exclusion) > 0) {
    legend_labels <- append(legend_labels, "Exclusion", after = 1)
    legend_colors <- append(legend_colors, "#991b1b", after = 1)
    legend_widths <- append(legend_widths, 2, after = 1)
  }
  legend(
    "bottomright",
    legend = legend_labels,
    col = legend_colors,
    lwd = legend_widths,
    bg = "white"
  )
  invisible(path)
}

# Orchestration and optional performance report ----------------------------

elapsed_seconds <- function(started_at) {
  unname(proc.time()[["elapsed"]] - started_at)
}

peak_rss_mib <- function() {
  status <- tryCatch(readLines("/proc/self/status", warn = FALSE),
    error = function(error) character())
  line <- grep("^VmHWM:", status, value = TRUE)
  if (length(line) == 0) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", line[[1]])) / 1024
}

cgroup_memory_peak_mib <- function() {
  candidates <- c(
    "/sys/fs/cgroup/memory.peak",
    "/sys/fs/cgroup/memory/memory.max_usage_in_bytes"
  )
  path <- candidates[file.exists(candidates)][1]
  if (is.na(path)) return(NA_real_)
  value <- suppressWarnings(as.numeric(readLines(path, n = 1L, warn = FALSE)))
  if (!is.finite(value)) NA_real_ else value / 1024^2
}

temporary_directory_bytes <- function(path) {
  if (!dir.exists(path)) return(0)
  paths <- list.files(
    path,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    no.. = TRUE
  )
  if (length(paths) == 0L) return(0)
  # The staged source cloud is a symlink and is not temporary disk usage.
  paths <- paths[!nzchar(Sys.readlink(paths))]
  sizes <- file.info(paths)$size
  sum(sizes[is.finite(sizes)], na.rm = TRUE)
}

write_performance_report <- function(parameters, footprint_source, point_count, tile_count,
                                     instance_dimensions, catalog_workers,
                                     timings, path) {
  report <- data.table(
    point_cloud = basename(parameters$point_cloud),
    footprint_source = footprint_source,
    point_count = point_count,
    tile_count = tile_count,
    instance_dimension = if (length(instance_dimensions) == 0) {
      NA_character_
    } else {
      paste(instance_dimensions, collapse = "|")
    },
    threads_requested = parameters$threads,
    threads_effective = get_lidr_threads(),
    catalog_workers = catalog_workers,
    peak_rss_mib = round(peak_rss_mib(), 3),
    cgroup_memory_peak_mib = round(cgroup_memory_peak_mib(), 3),
    temporary_disk_peak_mib = round(timings$temporary_disk_peak_mib, 3),
    grid_seconds = round(timings$grid, 6),
    dtm_seconds = round(timings$dtm, 6),
    chm_seconds = round(timings$chm, 6),
    segment_seconds = round(timings$segments, 6),
    tile_seconds = round(timings$tiles, 6),
    output_seconds = round(timings$outputs, 6),
    total_seconds = round(timings$total, 6),
    tile_size = parameters$tile_size,
    grid_search_step = parameters$grid_search_step,
    ptd_resolution = parameters$ptd_resolution,
    dtm_resolution = parameters$dtm_resolution,
    maximum_height = parameters$maximum_height,
    voxel_resolution = parameters$voxel_resolution,
    vegetation_minimum_height = parameters$vegetation_minimum_height,
    chm_resolution = parameters$chm_resolution,
    gap_height_threshold = parameters$gap_height_threshold,
    minimum_tree_voxels = parameters$minimum_tree_voxels,
    apex_minimum_height = parameters$apex_minimum_height,
    minimum_tree_thickness = parameters$minimum_tree_thickness,
    minimum_occupied_layers = parameters$minimum_occupied_layers,
    dtm_chunk_size = parameters$dtm_chunk_size,
    chunk_size = parameters$chunk_size,
    target_chunk_points = parameters$adaptive_chunking$target_points,
    adaptive_chunking_enabled = parameters$adaptive_chunking$enabled,
    adaptive_chunking_method = parameters$adaptive_chunking$method,
    adaptive_chunk_size_m = parameters$adaptive_chunking$chunk_size_m,
    adaptive_mean_density_points_m2 = parameter_or(
      parameters$adaptive_chunking,
      "mean_density_points_m2",
      NA_real_
    ),
    adaptive_estimated_points_per_full_chunk = parameter_or(
      parameters$adaptive_chunking,
      "estimated_points_per_full_chunk",
      NA_real_
    ),
    adaptive_chm_observed_chunks = timings$chm_observations$count,
    adaptive_chm_observed_minimum_points =
      timings$chm_observations$minimum_points,
    adaptive_chm_observed_maximum_points =
      timings$chm_observations$maximum_points,
    adaptive_chm_oversized_chunks = timings$chm_observations$oversized_count,
    adaptive_segment_observed_chunks = timings$segment_observations$count,
    adaptive_segment_observed_minimum_points =
      timings$segment_observations$minimum_points,
    adaptive_segment_observed_maximum_points =
      timings$segment_observations$maximum_points,
    adaptive_segment_oversized_chunks =
      timings$segment_observations$oversized_count,
    dtm_buffer = parameters$dtm_buffer,
    dtm_storage_mode = "disk_backed_vrt",
    chm_storage_mode = "disk_backed_vrt",
    execution_mode = if (isTRUE(parameters$sequential_instance_dimensions)) {
      "julia_memory_safe"
    } else {
      "default"
    },
    internal_memory_budget_gib = parameters$memory_budget_gib,
    segment_catalog_base_selection =
      parameters$segment_catalog_base_selection,
    instance_dimension_candidates = paste(parameters$instance_dimension, collapse = "|")
  )
  fwrite(report, path, na = "NA")
  message(sprintf(
    paste0(
      "PERF total=%.3fs grid=%.3fs dtm=%.3fs chm=%.3fs ",
      "segments=%.3fs tiles=%.3fs outputs=%.3fs peak_rss=%.1fMiB ",
      "cgroup_peak=%.1fMiB"
    ),
    timings$total, timings$grid, timings$dtm, timings$chm, timings$segments,
    timings$tiles, timings$outputs, report$peak_rss_mib[[1]],
    report$cgroup_memory_peak_mib[[1]]
  ))
  invisible(path)
}

run_analysis <- function(parameters) {
  parameters$sensor <- parameter_or(parameters, "sensor", NA_character_)
  parameters$country <- parameter_or(parameters, "country", NA_character_)
  parameters$memory_budget_gib <- parameter_or(
    parameters, "memory_budget_gib", Inf
  )
  parameters$sequential_instance_dimensions <- parameter_or(
    parameters, "sequential_instance_dimensions", FALSE
  )
  parameters$segment_catalog_base_selection <- parameter_or(
    parameters,
    "segment_catalog_base_selection",
    SEGMENT_CATALOG_BASE_SELECTION
  )
  parameters$segment_bucket_count <- parameter_or(
    parameters,
    "segment_bucket_count",
    SEGMENT_BUCKET_COUNT
  )
  parameters$dimension_point_clouds <- parameter_or(
    parameters, "dimension_point_clouds", list()
  )
  # NULL retains original LAZ behavior; COPC runners provide OriginalPointIndex.
  parameters$point_order_dimension <- parameter_or(
    parameters, "point_order_dimension", NULL
  )
  parameters$point_cloud_display_name <- parameter_or(
    parameters,
    "point_cloud_display_name",
    basename(parameters$point_cloud)
  )
  work_root <- parameter_or(parameters, "temp_dir", parameters$output_dir)
  if (!dir.exists(work_root) || file.access(work_root, mode = 2) != 0) {
    stop("analysis temporary directory must exist and be writable")
  }
  run_work_directory <- file.path(
    work_root,
    sprintf(".%s_forest_structure_%d", parameters$dataset_id, Sys.getpid())
  )
  dir.create(run_work_directory, recursive = TRUE, showWarnings = FALSE)
  temporary_disk_peak_bytes <- temporary_directory_bytes(run_work_directory)
  on.exit(unlink(run_work_directory, recursive = TRUE, force = TRUE), add = TRUE)
  total_started <- proc.time()[["elapsed"]]
  point_count <- point_count_from_las_header(parameters$point_cloud)
  assert_lidr_point_count_supported(point_count)
  parameters$adaptive_chunking <- adaptive_chunk_plan(
    point_count,
    point_cloud_xy_bounds(parameters$point_cloud),
    parameters$chunk_size,
    parameter_or(parameters, "target_chunk_points", 0),
    grepl("\\.copc\\.laz$", parameters$point_cloud, ignore.case = TRUE)
  )
  parameters$effective_chunk_size <- parameters$adaptive_chunking$chunk_size_m
  message(sprintf(
    paste0(
      "ADAPTIVE_CHUNK enabled=%s method=%s target_points=%.0f density_points_m2=%s ",
      "maximum_chunk_size_m=%.3f selected_chunk_size_m=%.3f ",
      "estimated_points_per_full_chunk=%s"
    ),
    tolower(as.character(parameters$adaptive_chunking$enabled)),
    parameters$adaptive_chunking$method,
    parameters$adaptive_chunking$target_points,
    if (is.finite(parameter_or(
      parameters$adaptive_chunking,
      "mean_density_points_m2",
      NA_real_
    ))) {
      sprintf("%.6f", parameters$adaptive_chunking$mean_density_points_m2)
    } else {
      "NA"
    },
    parameters$chunk_size,
    parameters$effective_chunk_size,
    if (is.finite(parameters$adaptive_chunking$estimated_points_per_full_chunk)) {
      sprintf("%.0f", parameters$adaptive_chunking$estimated_points_per_full_chunk)
    } else {
      "NA"
    }
  ))
  source_crs <- parameter_or(
    parameters,
    "source_crs",
    spatial_reference_wkt(readLAScatalog(parameters$point_cloud))
  )
  paths <- artifact_paths(parameters$output_dir, parameters$dataset_id)
  segment_work_directory <- file.path(run_work_directory, "segment_chunks")
  on.exit(
    unlink(segment_work_directory, recursive = TRUE, force = TRUE),
    add = TRUE
  )
  if (parameters$threads > 0) {
    set_lidr_threads(parameters$threads)
    message(sprintf("Using %d lidR threads", get_lidr_threads()))
  } else {
    message(sprintf("Using lidR default thread count: %d", get_lidr_threads()))
  }
  thread_budget <- max(
    1L,
    as.integer(
      if (parameters$threads > 0) parameters$threads else get_lidr_threads()
    )
  )
  catalog_workers <- min(
    thread_budget,
    max(1L, as.integer(parameters$catalog_workers))
  )
  message(sprintf(
    "Using %d memory-bounded LAScatalog worker(s) within the %d-thread allowance",
    catalog_workers,
    thread_budget
  ))
  message("FORESTSTRUCTURE_STAGE stage=grid status=started")
  grid_started <- proc.time()[["elapsed"]]
  aoi <- read_analysis_footprint(parameters$point_cloud, parameters$aoi)
  if (identical(aoi$source, "point_cloud_extent")) {
    message("No Audit AOI supplied; covering the complete point-cloud XY extent")
  } else {
    message(sprintf("Using Audit AOI: %s", basename(parameters$aoi)))
  }
  tiles <- build_analysis_tiles(
    aoi,
    parameters$tile_size,
    parameters$grid_search_step,
    thread_budget
  )
  grid_seconds <- elapsed_seconds(grid_started)
  message(sprintf("Optimized grid contains %d complete Analysis Tiles", nrow(tiles)))
  message("FORESTSTRUCTURE_STAGE stage=grid status=completed")

  message("FORESTSTRUCTURE_STAGE stage=dtm status=started")
  dtm_started <- proc.time()[["elapsed"]]
  dtm_path <- paths$dtm
  write_global_dtm(
    parameters$point_cloud,
    dtm_path,
    file.path(run_work_directory, "dtm_chunks"),
    parameters$dtm_chunk_size,
    parameters$dtm_buffer,
    parameters$dtm_resolution,
    parameters$ptd_resolution,
    catalog_workers,
    source_crs,
    # Propagate the key through DTM catalog selection and per-chunk restoration.
    parameters$point_order_dimension
  )
  dtm <- terra::rast(dtm_path)
  if (terra::inMemory(dtm)) {
    stop("published DTM unexpectedly materialized in memory")
  }
  dtm_seconds <- elapsed_seconds(dtm_started)
  temporary_disk_peak_bytes <- max(
    temporary_disk_peak_bytes,
    temporary_directory_bytes(run_work_directory)
  )
  message("FORESTSTRUCTURE_STAGE stage=dtm status=completed")

  message("FORESTSTRUCTURE_STAGE stage=chm status=started")
  chm_started <- proc.time()[["elapsed"]]
  chm_path <- paths$chm
  chm_observations <- write_global_chm(
    parameters$point_cloud,
    dtm_path,
    chm_path,
    file.path(run_work_directory, "chm_chunks"),
    parameters$effective_chunk_size,
    parameters$voxel_resolution,
    parameters$maximum_height,
    parameters$chm_resolution,
    catalog_workers,
    source_crs,
    # Propagate the same key through CHM selection and per-chunk restoration.
    parameters$point_order_dimension,
    if (parameters$adaptive_chunking$enabled) {
      parameters$adaptive_chunking$target_points
    } else {
      0
    }
  )
  chm_seconds <- elapsed_seconds(chm_started)
  temporary_disk_peak_bytes <- max(
    temporary_disk_peak_bytes,
    temporary_directory_bytes(run_work_directory)
  )
  message("FORESTSTRUCTURE_STAGE stage=chm status=completed")

  message("FORESTSTRUCTURE_STAGE stage=segmentation status=started")
  segment_started <- proc.time()[["elapsed"]]
  instance_dimensions <- select_instance_dimensions(
    parameters$point_cloud,
    parameters$instance_dimension
  )
  segments_by_dimension <- list()
  segment_observation_summaries <- list()
  if (length(instance_dimensions) == 0) {
    message("No configured Instance Dimensions found; tree and segment metrics will be NA")
  } else {
    message(sprintf(
      "Using Instance Dimensions: %s",
      paste(instance_dimensions, collapse = ", ")
    ))
    if (isTRUE(parameters$sequential_instance_dimensions)) {
      for (instance_dimension in instance_dimensions) {
        dimension_source <- parameters$dimension_point_clouds[[instance_dimension]]
        if (is.null(dimension_source)) dimension_source <- parameters$point_cloud
        dimension_directory <- file.path(
          segment_work_directory,
          gsub("[^A-Za-z0-9_-]", "_", instance_dimension)
        )
        message(sprintf(
          "Starting independent segment pass for %s from %s",
          instance_dimension,
          basename(dimension_source)
        ))
        accumulated <- accumulate_segments(
          dimension_source,
          dtm_path,
          parameters$effective_chunk_size,
          parameters$voxel_resolution,
          parameters$maximum_height,
          instance_dimension,
          catalog_workers,
          dimension_directory,
          bucket_count = parameters$segment_bucket_count,
          base_selection = parameters$segment_catalog_base_selection,
          # Restore source order separately in each dimension's COPC pass.
          point_order_dimension = parameters$point_order_dimension,
          target_chunk_points = if (parameters$adaptive_chunking$enabled) {
            parameters$adaptive_chunking$target_points
          } else {
            0
          }
        )
        segment_observation_summaries[[instance_dimension]] <-
          accumulated$observations
        dimension_segments <- if (is.null(accumulated)) {
          NULL
        } else {
          assign_segments_to_tiles(
            finalize_segment_store(accumulated, parameters),
            tiles
          )
        }
        segments_by_dimension[instance_dimension] <- list(dimension_segments)
        if (!is.null(dimension_segments)) {
          message(sprintf(
            "Segment pass %s found %d segments (%d accepted trees)",
            instance_dimension,
            nrow(dimension_segments),
            sum(dimension_segments$is_tree, na.rm = TRUE)
          ))
        }
        temporary_disk_peak_bytes <- max(
          temporary_disk_peak_bytes,
          temporary_directory_bytes(run_work_directory)
        )
        rm(accumulated)
        unlink(dimension_directory, recursive = TRUE, force = TRUE)
        gc(verbose = FALSE)
      }
    } else {
      accumulated <- accumulate_segments(
        parameters$point_cloud,
        dtm_path,
        parameters$effective_chunk_size,
        parameters$voxel_resolution,
        parameters$maximum_height,
        instance_dimensions,
        catalog_workers,
        segment_work_directory,
        bucket_count = parameters$segment_bucket_count,
        base_selection = parameters$segment_catalog_base_selection,
        # Restore source order in the combined legacy segmentation pass too.
        point_order_dimension = parameters$point_order_dimension,
        target_chunk_points = if (parameters$adaptive_chunking$enabled) {
          parameters$adaptive_chunking$target_points
        } else {
          0
        }
      )
      segment_observation_summaries[["combined"]] <- accumulated$observations
      segments <- if (is.null(accumulated)) {
        NULL
      } else {
        assign_segments_to_tiles(
          finalize_segment_store(accumulated, parameters),
          tiles
        )
      }
      for (instance_dimension in instance_dimensions) {
        target_dimension <- instance_dimension
        segments_by_dimension[instance_dimension] <- list(
          if (is.null(segments)) NULL else
            segments[instance_dimension == target_dimension]
        )
      }
      if (!is.null(segments)) {
        message(sprintf(
          "Global segment pass found %d segments (%d accepted trees) across %d dimension(s)",
          nrow(segments),
          sum(segments$is_tree, na.rm = TRUE),
          length(instance_dimensions)
        ))
      }
      temporary_disk_peak_bytes <- max(
        temporary_disk_peak_bytes,
        temporary_directory_bytes(run_work_directory)
      )
      rm(accumulated)
    }
    unlink(segment_work_directory, recursive = TRUE, force = TRUE)
  }

  segment_seconds <- elapsed_seconds(segment_started)
  segment_observations <- combine_chunk_observations(
    segment_observation_summaries
  )
  message("FORESTSTRUCTURE_STAGE stage=segmentation status=completed")

  edge_flags <- compute_edge_flags(tiles, parameters$tile_size)
  message("FORESTSTRUCTURE_STAGE stage=tiles status=started")
  tile_started <- proc.time()[["elapsed"]]
  tile_indices <- seq_len(nrow(tiles))
  tile_plan <- tile_worker_plan(thread_budget, length(tile_indices))
  tile_workers <- tile_plan$process_workers
  tile_threads_per_worker <- tile_plan$threads_per_worker
  calculate_row <- function(index) {
    if (tile_workers > 1L) {
      set_lidr_threads(tile_threads_per_worker)
    }
    worker_dtm <- rast(dtm_path)
    metrics <- calculate_tile_metrics(
      parameters$point_cloud,
      tiles[index, ],
      worker_dtm,
      parameters
    )
    build_tile_row(
      parameters$point_cloud_display_name,
      tiles[index, ],
      edge_flags[[index]],
      nrow(tiles),
      metrics,
      empty_tree_metrics(),
      parameters$sensor,
      parameters$country
    )
  }
  if (length(tile_indices) > 0L) {
    message(sprintf(
      "Calculating %d tile metric rows sequentially with %d lidR thread(s)",
      length(tile_indices),
      tile_threads_per_worker
    ))
  }
  base_result <- collect_tile_rows(tile_indices, calculate_row, tile_workers)
  tile_seconds <- elapsed_seconds(tile_started)
  message("FORESTSTRUCTURE_STAGE stage=tiles status=completed")

  message("FORESTSTRUCTURE_STAGE stage=output status=started")
  output_started <- proc.time()[["elapsed"]]
  output_dimensions <- if (length(instance_dimensions) == 0) {
    list(NULL)
  } else {
    as.list(instance_dimensions)
  }
  output_paths <- lapply(output_dimensions, function(instance_dimension) {
    target_dimension <- instance_dimension
    dimension_paths <- dimension_artifact_paths(paths, instance_dimension)
    dimension_segments <- if (is.null(instance_dimension)) {
      NULL
    } else {
      segments_by_dimension[[target_dimension]]
    }
    result <- result_with_tree_metrics(base_result, dimension_segments)
    write_upstream_csv(result, dimension_paths$results)
    write_segment_diagnostics(
      dimension_segments,
      parameters$point_cloud_display_name,
      dimension_paths$segment_diagnostics,
      parameters$sensor,
      parameters$country
    )
    write_tile_geojson(
      tiles,
      result,
      source_crs,
      dimension_paths$tiles_geojson
    )
    message(sprintf(
      "Wrote %d Analysis Tile rows for Instance Dimension %s to %s",
      nrow(result),
      if (is.null(instance_dimension)) "NA" else instance_dimension,
      dimension_paths$results
    ))
    dimension_paths
  })
  write_layout_png(
    aoi,
    tiles,
    edge_flags,
    parameters$point_cloud_display_name,
    parameters$tile_size,
    paths$tiles_png
  )
  output_seconds <- elapsed_seconds(output_started)
  if (parameters$performance_report) {
    write_performance_report(
      parameters,
      aoi$source,
      point_count,
      nrow(tiles),
      instance_dimensions,
      catalog_workers,
      list(
        grid = grid_seconds,
        dtm = dtm_seconds,
        chm = chm_seconds,
        segments = segment_seconds,
        tiles = tile_seconds,
        outputs = output_seconds,
        temporary_disk_peak_mib = temporary_disk_peak_bytes / 1024^2,
        chm_observations = chm_observations,
        segment_observations = segment_observations,
        total = elapsed_seconds(total_started)
      ),
      paths$performance
    )
  }
  message("FORESTSTRUCTURE_STAGE stage=output status=completed")
  invisible(output_paths)
}

main <- function() {
  run_analysis(parse_parameters())
}

if (!identical(Sys.getenv("FORESTSTRUCTURE_SOURCE_ONLY"), "1")) {
  tryCatch(
    main(),
    error = function(error) {
      message("forest-structure analysis failed: ", conditionMessage(error))
      quit(status = 1)
    }
  )
}
