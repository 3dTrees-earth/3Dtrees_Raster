suppressPackageStartupMessages({
  library(argparse)
  library(jsonlite)
  library(lidR)
})

# Hard ceiling validated for this execution mode; the caller may request less.
MAX_MEMORY_BUDGET_GIB <- 70
# Canonical COPC files carry this UInt64 source-row key as their first ExtraByte.
POINT_ORDER_DIMENSION <- "OriginalPointIndex"

# Identify COPC from the standard double extension while accepting case variants.
is_copc_point_cloud <- function(path) {
  grepl("\\.copc\\.laz$", path, ignore.case = TRUE)
}

# Read exact header bounds without loading any point records.  These values form
# the inexpensive companion-file identity check below.
point_cloud_xyz_bounds <- function(path) {
  # lidR exposes LAS/COPC public-header fields through the same header object.
  header <- readLASheader(path)
  # Return a named vector so comparison order is explicit and reviewable.
  c(
    xmin = as.numeric(header@PHB[["Min X"]]),
    ymin = as.numeric(header@PHB[["Min Y"]]),
    zmin = as.numeric(header@PHB[["Min Z"]]),
    xmax = as.numeric(header@PHB[["Max X"]]),
    ymax = as.numeric(header@PHB[["Max Y"]]),
    zmax = as.numeric(header@PHB[["Max Z"]])
  )
}

# The original LAZ is optional at runtime and is never used for computation.
# When supplied, prove that it describes the same cloud as the requested COPC so
# provenance cannot accidentally name an unrelated ordered source.
validate_original_companion <- function(requested_point_cloud, point_cloud) {
  # Compare header counts without paying the cost of two full point scans.
  requested_count <- point_count_from_las_header(requested_point_cloud)
  source_count <- point_count_from_las_header(point_cloud)
  # Compare all six XYZ extrema at a tolerance far below source quantization.
  requested_bounds <- point_cloud_xyz_bounds(requested_point_cloud)
  source_bounds <- point_cloud_xyz_bounds(point_cloud)
  # Either disagreement means the companion cannot safely describe this COPC.
  if (requested_count != source_count ||
      any(abs(requested_bounds - source_bounds) > 1e-9)) {
    stop(paste(
      "--original-point-cloud does not match the COPC header point count",
      "and XYZ bounds"
    ))
  }
  # Return invisibly because success is a guard, not a scientific value.
  invisible(TRUE)
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
# Remember whether this file was intentionally sourced for focused tests before
# temporarily suppressing the entry points of the shared implementation files.
julia_memory_safe_source_only <- identical(
  Sys.getenv("FORESTSTRUCTURE_SOURCE_ONLY"),
  "1"
)
Sys.setenv(FORESTSTRUCTURE_SOURCE_ONLY = "1")
source(file.path(script_directory, "run.R"))
source(file.path(script_directory, "aoi_conversion.R"))

parse_julia_memory_safe_parameters <- function() {
  parser <- ArgumentParser(
    description = paste(
      "Julia-faithful ForestStructure analysis with selective ExtraBytes",
      "and disk-backed memory guards"
    )
  )
  # ``--point-cloud`` is the sole scientific input and may now be canonical COPC.
  parser$add_argument("--point-cloud", required = TRUE)
  # ``--original-point-cloud`` is optional identity/provenance evidence only.
  parser$add_argument("--original-point-cloud", default = NULL)
  parser$add_argument("--aoi-geojson", required = TRUE)
  parser$add_argument("--dataset-id", required = TRUE, type = "integer")
  parser$add_argument("--output-dir", required = TRUE)
  parser$add_argument("--temp-dir", required = TRUE)
  parser$add_argument("--memory-budget-gib", type = "double", default = 60)
  parser$add_argument(
    "--target-chunk-points",
    type = "double",
    default = 50000000,
    help = paste(
      "Approximate COPC CHM/segment points per spatial chunk; 0 disables",
      "adaptive sizing (default: 50000000)"
    )
  )
  parser$add_argument("--sensor", default = "NA")
  parser$add_argument("--country", default = "NA")
  parser$add_argument(
    "--instance-dimension",
    action = "append",
    default = NULL
  )
  arguments <- parser$parse_args()
  raw_instance_dimensions <- unlist(
    arguments$instance_dimension,
    recursive = TRUE,
    use.names = FALSE
  )
  if (is.null(raw_instance_dimensions) ||
      length(raw_instance_dimensions) == 0L ||
      all(is.na(raw_instance_dimensions))) {
    arguments$instance_dimension <- c(
      "PredInstance",
      "PredInstance_SAT",
      "PredInstance_FM"
    )
  } else {
    raw_instance_dimensions <- raw_instance_dimensions[
      !is.na(raw_instance_dimensions)
    ]
    arguments$instance_dimension <- unique(trimws(unlist(strsplit(
      as.character(raw_instance_dimensions),
      ",",
      fixed = TRUE
    ))))
    arguments$instance_dimension <- arguments$instance_dimension[
      nzchar(arguments$instance_dimension)
    ]
  }
  if (!file.exists(arguments$point_cloud) ||
      !grepl("\\.la[sz]$", arguments$point_cloud, ignore.case = TRUE)) {
    stop("--point-cloud must be exactly one existing LAS/LAZ file")
  }
  # COPC accepts an optional non-COPC companion; computation still reads COPC.
  if (is_copc_point_cloud(arguments$point_cloud)) {
    # Reject a missing, wrongly typed, or second-COPC companion before staging.
    if (!is.null(arguments$original_point_cloud) &&
        (!file.exists(arguments$original_point_cloud) ||
        !grepl(
          "\\.la[sz]$",
          arguments$original_point_cloud,
          ignore.case = TRUE
        ) ||
        is_copc_point_cloud(arguments$original_point_cloud))) {
      stop(paste(
        "--original-point-cloud must be an existing non-COPC LAS/LAZ",
        "when supplied for COPC identity validation"
      ))
    }
  # A companion has no role when the requested source is already ordered LAS/LAZ.
  } else if (!is.null(arguments$original_point_cloud)) {
    stop("--original-point-cloud is only valid when --point-cloud is a COPC")
  }
  if (!file.exists(arguments$aoi_geojson) ||
      !grepl("\\.(geojson|json)$", arguments$aoi_geojson, ignore.case = TRUE)) {
    stop("--aoi-geojson must be exactly one existing GeoJSON file")
  }
  if (!is.finite(arguments$dataset_id) || arguments$dataset_id < 1) {
    stop("--dataset-id must be a positive integer")
  }
  for (directory in c(arguments$output_dir, arguments$temp_dir)) {
    if (!dir.exists(directory) || file.access(directory, mode = 2) != 0) {
      stop("output and temporary directories must exist and be writable")
    }
  }
  if (!is.finite(arguments$memory_budget_gib) ||
      arguments$memory_budget_gib <= 0 ||
      arguments$memory_budget_gib > MAX_MEMORY_BUDGET_GIB) {
    stop(sprintf(
      "--memory-budget-gib must be greater than zero and at most %s",
      MAX_MEMORY_BUDGET_GIB
    ))
  }
  if (!is.finite(arguments$target_chunk_points) ||
      arguments$target_chunk_points < 0) {
    stop("--target-chunk-points must be zero or greater")
  }
  arguments
}

extra_byte_definitions <- function(point_cloud) {
  header <- readLASheader(point_cloud)
  descriptions <- header@VLR$Extra_Bytes[["Extra Bytes Description"]]
  if (is.null(descriptions)) return(list())
  lapply(seq_along(descriptions), function(index) {
    description <- descriptions[[index]]
    list(
      name = as.character(description$name),
      ordinal = index,
      data_type = as.integer(description$data_type)
    )
  })
}

run_projection <- function(source, output, dimension) {
  output_lines <- system2(
    "python3",
    c(
      file.path(script_directory, "project_instance_dimension.py"),
      "--source", source,
      "--output", output,
      "--dimension", dimension
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output_lines, "status")
  if (!is.null(status) && status != 0L) {
    stop(sprintf(
      "streaming projection failed for %s: %s",
      dimension,
      paste(output_lines, collapse = " ")
    ))
  }
  fromJSON(tail(output_lines, 1L), simplifyVector = FALSE)
}

expected_artifacts <- function(directory, dataset_id, dimensions) {
  prefix <- file.path(directory, as.character(dataset_id))
  common <- c(
    paste0(prefix, "_dtm.tif"),
    paste0(prefix, "_chm.tif"),
    paste0(prefix, "_tiles.png"),
    paste0(prefix, "_performance.csv"),
    paste0(prefix, "_aoi_conversion.json"),
    paste0(prefix, "_julia_memory_safe_run.json")
  )
  suffixes <- if (length(dimensions) == 0L) "" else paste0("_", dimensions)
  dimension_files <- unlist(lapply(suffixes, function(suffix) c(
    paste0(prefix, suffix, "_results.csv"),
    paste0(prefix, suffix, "_segment_diagnostics.csv"),
    paste0(prefix, suffix, "_tiles.geojson")
  )))
  c(common, dimension_files)
}

validate_incoming_artifacts <- function(directory, dataset_id, dimensions,
                                        memory_budget_gib, source_point_count,
                                        point_cloud_name, sensor, country) {
  matches_label <- function(values, expected) {
    if (identical(expected, "NA")) {
      return(all(is.na(values) | as.character(values) == "NA"))
    }
    all(!is.na(values) & as.character(values) == expected)
  }
  paths <- expected_artifacts(directory, dataset_id, dimensions)
  absent <- paths[!file.exists(paths) | file.info(paths)$size <= 0]
  if (length(absent) > 0L) {
    stop(sprintf(
      "incoming artifact validation failed; missing/empty: %s",
      paste(basename(absent), collapse = ", ")
    ))
  }
  performance <- fread(
    file.path(directory, paste0(dataset_id, "_performance.csv"))
  )
  if (nrow(performance) != 1L ||
      performance$point_count[[1]] != source_point_count ||
      !is.finite(performance$peak_rss_mib[[1]]) ||
      performance$peak_rss_mib[[1]] >= memory_budget_gib * 1024 ||
      !is.finite(performance$cgroup_memory_peak_mib[[1]]) ||
      performance$cgroup_memory_peak_mib[[1]] >= memory_budget_gib * 1024 ||
      !is.finite(performance$temporary_disk_peak_mib[[1]])) {
    stop("incoming artifact validation failed: performance/resource checks failed")
  }
  dtm <- terra::rast(file.path(directory, paste0(dataset_id, "_dtm.tif")))
  chm <- terra::rast(file.path(directory, paste0(dataset_id, "_chm.tif")))
  if (terra::inMemory(dtm) || terra::inMemory(chm) ||
      any(abs(terra::res(dtm) - c(1, 1)) > 1e-9) ||
      any(abs(terra::res(chm) - c(0.5, 0.5)) > 1e-9) ||
      !identical(terra::crs(dtm), terra::crs(chm))) {
    stop("incoming raster resolution, storage mode, or CRS validation failed")
  }
  for (dimension in dimensions) {
    results <- fread(file.path(
      directory,
      sprintf("%s_%s_results.csv", dataset_id, dimension)
    ))
    diagnostics <- fread(file.path(
      directory,
      sprintf("%s_%s_segment_diagnostics.csv", dataset_id, dimension)
    ))
    required_results <- names(empty_result_table())
    required_diagnostics <- names(empty_segment_diagnostics())
    if (!identical(names(results), required_results) ||
        !identical(names(diagnostics), required_diagnostics)) {
      stop(sprintf(
        "incoming CSV schema differs from Julia for %s",
        dimension
      ))
    }
    if (nrow(results) != performance$tile_count[[1]] ||
        any(results$file != point_cloud_name) ||
        !matches_label(results$sensor, sensor) ||
        !matches_label(results$country, country) ||
        any(diagnostics$file != point_cloud_name) ||
        !matches_label(diagnostics$sensor, sensor) ||
        !matches_label(diagnostics$country, country)) {
      stop(sprintf(
        "incoming row count or metadata validation failed for %s",
        dimension
      ))
    }
    tile_features <- st_read(
      file.path(
        directory,
        sprintf("%s_%s_tiles.geojson", dataset_id, dimension)
      ),
      quiet = TRUE,
      stringsAsFactors = FALSE
    )
    if (nrow(tile_features) != nrow(results)) {
      stop(sprintf("incoming tile GeoJSON row count differs for %s", dimension))
    }
  }
  invisible(paths)
}

promote_artifacts <- function(incoming_artifacts, output_directory, dataset_id) {
  # Promote only paths returned by validate_incoming_artifacts(). GDAL may
  # leave non-scientific PAM sidecars beside TIFFs; those files are cleaned up
  # with the staging directory and never enter the published artifact set.
  incoming <- normalizePath(incoming_artifacts, mustWork = TRUE)
  prefix_pattern <- sprintf("^%s(_|\\.)", dataset_id)
  existing <- list.files(
    output_directory,
    pattern = prefix_pattern,
    full.names = TRUE,
    recursive = FALSE,
    all.files = FALSE
  )
  backup_directory <- tempfile(
    pattern = sprintf(".backup_%s_", dataset_id),
    tmpdir = output_directory
  )
  dir.create(backup_directory)
  promoted <- character()
  success <- FALSE
  on.exit({
    if (!success) {
      unlink(promoted, force = TRUE)
      backup <- list.files(backup_directory, full.names = TRUE)
      for (path in backup) {
        file.rename(path, file.path(output_directory, basename(path)))
      }
    }
    unlink(backup_directory, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  for (path in existing) {
    if (!file.rename(path, file.path(backup_directory, basename(path)))) {
      stop(sprintf("cannot stage existing artifact %s for replacement", path))
    }
  }
  for (path in incoming) {
    target <- file.path(output_directory, basename(path))
    if (!file.rename(path, target)) {
      stop(sprintf("cannot atomically promote %s", basename(path)))
    }
    promoted <- c(promoted, target)
  }
  success <- TRUE
  invisible(promoted)
}

package_versions <- function() {
  packages <- c(
    "lidR", "data.table", "sf", "terra", "rlas", "future", "argparse",
    "jsonlite"
  )
  as.list(setNames(
    vapply(packages, function(package) {
      as.character(packageVersion(package))
    }, character(1)),
    packages
  ))
}

python_package_versions <- function() {
  output <- system2(
    "python3",
    c(
      "-c",
      shQuote(
        paste(
          "import importlib.metadata as m;",
          "print(m.version('laspy'));",
          "print(m.version('lazrs'))"
        )
      )
    ),
    stdout = TRUE
  )
  list(laspy = output[[1]], lazrs = output[[2]])
}

julia_memory_safe_main <- function() {
  # Parse and validate every external path/resource setting before creating work.
  arguments <- parse_julia_memory_safe_parameters()
  message("FORESTSTRUCTURE_STAGE stage=preflight status=started")
  # Retain the user-requested path separately from later local staging paths.
  requested_point_cloud <- normalizePath(arguments$point_cloud, mustWork = TRUE)
  # Normalize optional companion evidence once; NULL remains an explicit absence.
  original_point_cloud <- if (is.null(arguments$original_point_cloud)) {
    NULL
  } else {
    normalizePath(arguments$original_point_cloud, mustWork = TRUE)
  }
  # This single flag controls order restoration, spatial-index handling, and
  # provenance so COPC behavior cannot diverge across independent checks.
  requested_is_copc <- is_copc_point_cloud(requested_point_cloud)
  # Scientific computation always reads the requested cloud (COPC when supplied).
  point_cloud <- requested_point_cloud
  # Preserve Julia's historical source filename in result rows when a companion
  # is available, even though point data are spatially streamed from COPC.
  point_cloud_display_name <- if (is.null(original_point_cloud)) {
    basename(point_cloud)
  } else {
    basename(original_point_cloud)
  }
  # Non-COPC records are already ordered; COPC must be sorted per chunk by key.
  point_order_dimension <- if (requested_is_copc) {
    POINT_ORDER_DIMENSION
  } else {
    NULL
  }
  # Enforce the canonical-COPC contract before running expensive analysis.
  if (requested_is_copc) {
    # Read declared ExtraBytes from the COPC header without materializing points.
    definitions <- extra_byte_definitions(requested_point_cloud)
    # Reduce definitions to names for the mandatory order-key membership check.
    available <- vapply(definitions, function(definition) {
      definition$name
    }, character(1))
    # A generic COPC cannot reproduce order-sensitive Julia outputs safely.
    if (!POINT_ORDER_DIMENSION %in% available) {
      stop(paste(
        "COPC scientific input requires the OriginalPointIndex ExtraByte;",
        "rebuild COPC from the ordered source before ForestStructure"
      ))
    }
    # Validate optional provenance evidence, but do not substitute it for COPC.
    if (!is.null(original_point_cloud)) {
      validate_original_companion(requested_point_cloud, original_point_cloud)
    }
    # Make the selected streaming/order mode explicit in operational logs.
    message(sprintf(
      "Using COPC spatial streaming with %s order restoration from %s",
      POINT_ORDER_DIMENSION,
      basename(point_cloud)
    ))
  }
  point_cloud_crs <- spatial_reference_wkt(readLAScatalog(point_cloud))
  original_companion_crs <- if (is.null(original_point_cloud)) {
    NA_character_
  } else {
    spatial_reference_wkt(readLAScatalog(original_point_cloud))
  }
  resolved_crs <- resolve_analysis_crs(
    point_cloud_crs,
    original_companion_crs
  )
  if (identical(resolved_crs$source, "original_companion")) {
    message(paste(
      "COPC CRS is unavailable; using the validated original companion",
      "header CRS for output metadata"
    ))
  }
  aoi_geojson <- normalizePath(arguments$aoi_geojson, mustWork = TRUE)
  input_sha256 <- sha256_file(point_cloud)
  job_directory <- tempfile(
    pattern = sprintf("foreststructure_%s_", arguments$dataset_id),
    tmpdir = arguments$temp_dir
  )
  incoming_directory <- tempfile(
    pattern = sprintf(".incoming_%s_", arguments$dataset_id),
    tmpdir = arguments$output_dir
  )
  dir.create(job_directory, recursive = TRUE)
  dir.create(incoming_directory)
  on.exit({
    unlink(job_directory, recursive = TRUE, force = TRUE)
    unlink(incoming_directory, recursive = TRUE, force = TRUE)
  }, add = TRUE)

  # Stage by symlink so all processing stays on the requested COPC without an
  # expensive duplicate copy into the per-job work directory.
  staged_point_cloud <- file.path(job_directory, basename(point_cloud))
  if (!file.symlink(point_cloud, staged_point_cloud)) {
    stop("cannot stage original point cloud by local-SSD symbolic link")
  }
  # Plain LAS/LAZ benefits from a sidecar LAX; COPC already contains its spatial
  # hierarchy and must not be rewritten or supplemented with a LAX index.
  if (!requested_is_copc) {
    lax_path <- sub("\\.la[sz]$", ".lax", staged_point_cloud, ignore.case = TRUE)
    tryCatch(
      rlas::writelax(staged_point_cloud),
      error = function(error) message(
        "Spatial-index creation failed; continuing with selective scans: ",
        conditionMessage(error)
      )
    )
    if (file.exists(lax_path)) {
      message(sprintf("Built local spatial index %s", basename(lax_path)))
    }
  }

  message("FORESTSTRUCTURE_STAGE stage=aoi_conversion status=started")
  generated_gpkg <- file.path(
    job_directory,
    paste0(tools::file_path_sans_ext(basename(point_cloud)), ".gpkg")
  )
  aoi_provenance <- convert_geojson_to_gpkg(
    aoi_geojson,
    generated_gpkg,
    point_cloud_xy_bounds(staged_point_cloud),
    tile_threads = 1L
  )
  write_json(
    aoi_provenance,
    file.path(
      incoming_directory,
      paste0(arguments$dataset_id, "_aoi_conversion.json")
    ),
    auto_unbox = TRUE,
    pretty = TRUE
  )
  message("FORESTSTRUCTURE_STAGE stage=aoi_conversion status=completed")

  message("FORESTSTRUCTURE_STAGE stage=dimension_projection status=started")
  definitions <- extra_byte_definitions(staged_point_cloud)
  available <- vapply(definitions, function(definition) {
    definition$name
  }, character(1))
  dimensions <- arguments$instance_dimension[
    arguments$instance_dimension %in% available
  ]
  missing <- setdiff(arguments$instance_dimension, dimensions)
  if (length(missing) > 0L) {
    message(sprintf(
      "Skipping missing Instance Dimensions: %s",
      paste(missing, collapse = ", ")
    ))
  }
  # lidR can selectively request only ExtraBytes 1-9 from COPC.  Requiring the
  # key and requested instance fields in that window keeps reads spatial and
  # avoids silently falling back to repeated full-cloud scans.
  if (requested_is_copc) {
    selective_dimensions <- unique(c(point_order_dimension, dimensions))
    selective_ordinals <- match(selective_dimensions, available)
    unsupported <- selective_dimensions[selective_ordinals > 9L]
    if (length(unsupported) > 0L) {
      stop(paste(
        "COPC selective streaming requires OriginalPointIndex and requested",
        "instance dimensions among the first nine ExtraBytes; rebuild the",
        "COPC with build_original_order_copc.sh. Unsupported:",
        paste(unsupported, collapse = ", ")
      ))
    }
  }

  # Track the actual source and selector used for each instance dimension; these
  # maps feed both execution and the durable provenance report.
  dimension_sources <- list()
  selectors <- list()
  projections <- list()
  # Resolve each requested segmentation field independently.
  for (dimension in dimensions) {
    # ExtraByte ordinal maps directly to lidR selector digits 1-9.
    ordinal <- match(dimension, available)
    definition <- definitions[[ordinal]]
    # Canonical COPC guarantees streamable placement.  Plain LAZ may still need
    # a one-dimension projection when its field is beyond lidR's selector window.
    if (ordinal <= 9L || requested_is_copc) {
      dimension_sources[[dimension]] <- staged_point_cloud
      selector <- if (ordinal <= 9L) paste0("xyz", ordinal) else "xyz0"
    } else {
      projection_path <- file.path(
        job_directory,
        sprintf("%s_%s_projected.laz", arguments$dataset_id, dimension)
      )
      projections[[dimension]] <- run_projection(
        staged_point_cloud,
        projection_path,
        dimension
      )
      dimension_sources[[dimension]] <- projection_path
      selector <- "xyz1"
    }
    selectors[[dimension]] <- list(
      requested_dimension = dimension,
      source_ordinal = ordinal,
      source_data_type = definition$data_type,
      effective_selector = selector
    )
    message(sprintf(
      "SELECTOR phase=segment dimension=%s ordinal=%d selector=%s",
      dimension,
      ordinal,
      selector
    ))
  }
  # Resolve the order-key selector once for DTM, CHM, and tile reads.
  point_order_ordinal <- if (is.null(point_order_dimension)) {
    NA_integer_
  } else {
    match(point_order_dimension, available)
  }
  # ``xyz`` is sufficient for already ordered LAZ; canonical COPC normally uses
  # ``xyz1`` because OriginalPointIndex is deliberately the first ExtraByte.
  point_order_selector <- if (is.na(point_order_ordinal)) {
    "xyz"
  } else if (point_order_ordinal <= 9L) {
    paste0("xyz", point_order_ordinal)
  } else {
    "xyz0"
  }
  message(sprintf("SELECTOR phase=dtm selector=%s", point_order_selector))
  message(sprintf("SELECTOR phase=tile selector=%s", point_order_selector))
  message("FORESTSTRUCTURE_STAGE stage=dimension_projection status=completed")

  thread_count <- suppressWarnings(as.integer(Sys.getenv(
    "FORESTSTRUCTURE_THREADS",
    "20"
  )))
  if (!is.finite(thread_count) || thread_count < 1L) thread_count <- 20L
  catalog_worker_count <- validated_catalog_worker_count(
    Sys.getenv("FORESTSTRUCTURE_CATALOG_WORKERS", "1"),
    thread_count
  )
  message(sprintf(
    "Using validated Julia-memory-safe catalog worker count: %d",
    catalog_worker_count
  ))
  segment_bucket_count <- suppressWarnings(as.integer(Sys.getenv(
    "FORESTSTRUCTURE_SEGMENT_BUCKET_COUNT",
    as.character(SEGMENT_BUCKET_COUNT)
  )))
  if (!is.finite(segment_bucket_count) ||
      segment_bucket_count < 1L ||
      segment_bucket_count > 65536L) {
    stop("FORESTSTRUCTURE_SEGMENT_BUCKET_COUNT must be an integer from 1 to 65536")
  }
  # Preserve Julia's scientific constants exactly.  Only storage, selection,
  # resource controls, and record-order restoration differ in this runner.
  parameters <- list(
    point_cloud = staged_point_cloud,
    dataset_id = arguments$dataset_id,
    aoi = generated_gpkg,
    output_dir = incoming_directory,
    temp_dir = job_directory,
    tile_size = 20,
    grid_search_step = 0.5,
    ptd_resolution = 20,
    dtm_resolution = 1,
    maximum_height = 70,
    voxel_resolution = 0.2,
    vegetation_minimum_height = 0.5,
    chm_resolution = 0.5,
    gap_height_threshold = 3,
    dtm_chunk_size = 60,
    chunk_size = 60,
    target_chunk_points = arguments$target_chunk_points,
    dtm_buffer = 20,
    threads = thread_count,
    catalog_workers = catalog_worker_count,
    performance_report = TRUE,
    instance_dimension = arguments$instance_dimension,
    segment_diagnostics = TRUE,
    minimum_tree_voxels = 100L,
    apex_minimum_height = 3,
    minimum_tree_thickness = 0.5,
    minimum_occupied_layers = 3L,
    sensor = arguments$sensor,
    country = arguments$country,
    memory_budget_gib = arguments$memory_budget_gib,
    sequential_instance_dimensions = TRUE,
    segment_catalog_base_selection = JULIA_MEMORY_SAFE_SEGMENT_SELECTION,
    segment_bucket_count = segment_bucket_count,
    # Per-dimension sources support legacy LAZ projection without changing COPC.
    dimension_point_clouds = dimension_sources,
    # Result metadata should identify the original ordered filename when known.
    point_cloud_display_name = point_cloud_display_name,
    # NULL disables sorting for LAZ; the key enables deterministic COPC sorting.
    point_order_dimension = point_order_dimension,
    # Header-only fallback preserves output CRS without reading companion points.
    source_crs = resolved_crs$wkt
  )
  message("FORESTSTRUCTURE_STAGE stage=analysis status=started")
  run_analysis(parameters)
  message("FORESTSTRUCTURE_STAGE stage=analysis status=completed")
  performance <- fread(file.path(
    incoming_directory,
    paste0(arguments$dataset_id, "_performance.csv")
  ))
  temporary_disk_peak_mib <- performance$temporary_disk_peak_mib[[1]] +
    temporary_directory_bytes(job_directory) / 1024^2

  # Persist enough provenance to prove that the run used COPC spatial streaming,
  # restored original order, and optionally validated the original companion.
  run_provenance <- list(
    status = "candidate",
    mode = "julia-memory-safe",
    rerun_policy = "always_recompute_then_replace_after_validation",
    source_point_cloud = normalizePath(point_cloud, mustWork = TRUE),
    source_point_cloud_sha256 = input_sha256,
    source_point_count = point_count_from_las_header(point_cloud),
    # ``source_*`` records the file actually read for scientific computation.
    source_is_copc = requested_is_copc,
    # ``requested_*`` records the caller-facing input before local symlinking.
    requested_point_cloud = requested_point_cloud,
    requested_is_copc = requested_is_copc,
    # This is true only when an explicit companion passed count/bounds checks.
    companion_header_identity_validated = (
      requested_is_copc && !is.null(original_point_cloud)
    ),
    original_point_cloud = original_point_cloud,
    output_crs_source = resolved_crs$source,
    # These fields make order restoration and streaming mode machine-auditable.
    point_order_dimension = point_order_dimension,
    catalog_spatial_streaming = if (requested_is_copc) "COPC" else "LAS/LAZ",
    requested_instance_dimensions = arguments$instance_dimension,
    processed_instance_dimensions = dimensions,
    missing_instance_dimensions = missing,
    extra_byte_definitions = definitions,
    selectors = selectors,
    projections = projections,
    catalog_workers = catalog_worker_count,
    data_table_threads = DATA_TABLE_THREADS,
    segment_bucket_count = segment_bucket_count,
    memory_budget_gib = arguments$memory_budget_gib,
    parent_peak_rss_mib = performance$peak_rss_mib[[1]],
    cgroup_memory_peak_mib = performance$cgroup_memory_peak_mib[[1]],
    temporary_disk_peak_mib = temporary_disk_peak_mib,
    adaptive_chunking = list(
      enabled = performance$adaptive_chunking_enabled[[1]],
      method = performance$adaptive_chunking_method[[1]],
      target_points = performance$target_chunk_points[[1]],
      configured_maximum_chunk_size_m = performance$chunk_size[[1]],
      selected_chunk_size_m = performance$adaptive_chunk_size_m[[1]],
      mean_density_points_m2 = performance$adaptive_mean_density_points_m2[[1]],
      estimated_points_per_full_chunk =
        performance$adaptive_estimated_points_per_full_chunk[[1]],
      observed = list(
        chm = list(
          chunks = performance$adaptive_chm_observed_chunks[[1]],
          minimum_points =
            performance$adaptive_chm_observed_minimum_points[[1]],
          maximum_points =
            performance$adaptive_chm_observed_maximum_points[[1]],
          oversized_chunks =
            performance$adaptive_chm_oversized_chunks[[1]]
        ),
        segmentation = list(
          chunks = performance$adaptive_segment_observed_chunks[[1]],
          minimum_points =
            performance$adaptive_segment_observed_minimum_points[[1]],
          maximum_points =
            performance$adaptive_segment_observed_maximum_points[[1]],
          oversized_chunks =
            performance$adaptive_segment_oversized_chunks[[1]]
        )
      )
    ),
    scientific_parameters = list(
      tile_size = 20,
      grid_search_step = 0.5,
      ptd_resolution = 20,
      dtm_resolution = 1,
      dtm_chunk_size = 60,
      dtm_buffer = 20,
      segment_chunk_size = performance$adaptive_chunk_size_m[[1]],
      voxel_resolution = 0.2,
      maximum_height = 70,
      vegetation_minimum_height = 0.5,
      chm_resolution = 0.5,
      gap_height_threshold = 3,
      minimum_tree_voxels = 100,
      apex_minimum_height = 3,
      minimum_tree_thickness = 0.5,
      minimum_occupied_layers = 3,
      classify_ground_last_returns = FALSE
    ),
    package_versions = package_versions(),
    python_package_versions = python_package_versions()
  )
  write_json(
    run_provenance,
    file.path(
      incoming_directory,
      paste0(arguments$dataset_id, "_julia_memory_safe_run.json")
    ),
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )
  message("FORESTSTRUCTURE_STAGE stage=internal_validation status=started")
  validated_artifacts <- validate_incoming_artifacts(
    incoming_directory,
    arguments$dataset_id,
    dimensions,
    arguments$memory_budget_gib,
    point_count_from_las_header(point_cloud),
    point_cloud_display_name,
    arguments$sensor,
    arguments$country
  )
  message("FORESTSTRUCTURE_STAGE stage=internal_validation status=completed")
  run_provenance$status <- "valid"
  write_json(
    run_provenance,
    file.path(
      incoming_directory,
      paste0(arguments$dataset_id, "_julia_memory_safe_run.json")
    ),
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )
  message("FORESTSTRUCTURE_STAGE stage=artifact_promotion status=started")
  promoted <- promote_artifacts(
    validated_artifacts,
    arguments$output_dir,
    arguments$dataset_id
  )
  message("FORESTSTRUCTURE_STAGE stage=artifact_promotion status=completed")
  message(sprintf(
    "Validated and promoted %d freshly recomputed artifact(s) for dataset %s",
    length(promoted),
    arguments$dataset_id
  ))
  invisible(promoted)
}

if (!julia_memory_safe_source_only) {
  tryCatch(
    julia_memory_safe_main(),
    error = function(error) {
      message("julia-memory-safe analysis failed: ", conditionMessage(error))
      quit(status = 1)
    }
  )
}
