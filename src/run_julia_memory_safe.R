suppressPackageStartupMessages({
  library(argparse)
  library(jsonlite)
  library(lidR)
})

# Canonical COPC files carry this UInt64 source-row key as their first ExtraByte.
POINT_ORDER_DIMENSION <- "OriginalPointIndex"

# Identify COPC from the standard double extension while accepting case variants.
is_copc_point_cloud <- function(path) {
  grepl("\\.copc\\.laz$", path, ignore.case = TRUE)
}

# The original LAZ is optional at runtime and is never used for computation.
# Reject incompatible companions without re-scanning the source. Lossless
# preparation validation (indices and fingerprints) remains the identity proof.
validate_original_companion <- function(requested_point_cloud, point_cloud) {
  # Compare header counts without paying the cost of two full point scans.
  requested_count <- point_count_from_las_header(requested_point_cloud)
  source_count <- point_count_from_las_header(point_cloud)
  # Compare the six extrema on their stored coordinate lattices.
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
source(file.path(script_directory, "production_contract.R"))

parse_julia_memory_safe_parameters <- function() {
  parser <- ArgumentParser(
    description = paste(
      "ForestStructure analysis with uncapped heights, selective ExtraBytes",
      "and disk-backed memory guards"
    )
  )
  add_voxel_count_arguments(parser, mandatory = TRUE)
  add_denoising_arguments(parser)
  # ``--point-cloud`` is the sole scientific input and may now be canonical COPC.
  parser$add_argument("--point-cloud", required = TRUE)
  # ``--original-point-cloud`` is optional identity/provenance evidence only.
  parser$add_argument("--original-point-cloud", default = NULL)
  parser$add_argument("--aoi-geojson", default = NULL,
    help = "Optional Audit AOI in source XY coordinates; otherwise use the full extent")
  # Retain old batch invocations without exposing database identifiers to Galaxy.
  parser$add_argument("--dataset-id", type = "integer", default = NULL,
    help = "Deprecated output-prefix compatibility; not needed for production jobs")
  parser$add_argument("--output-dir", default = getwd())
  parser$add_argument("--temp-dir", default = Sys.getenv("TMPDIR", tempdir()))
  parser$add_argument("--memory-budget-gib", type = "double", default = NULL)
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
    required = TRUE,
    help = paste(
      "Required; repeat or comma-separate PredInstance, PredInstance_SAT or",
      "PredInstance_FM. Every selected dimension must exist in the point cloud"
    )
  )
  arguments <- parser$parse_args()
  validate_voxel_count_parameters(arguments)
  validate_denoising_options(arguments)
  raw_instance_dimensions <- unlist(
    arguments$instance_dimension,
    recursive = TRUE,
    use.names = FALSE
  )
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
  if (!length(arguments$instance_dimension)) {
    stop("--instance-dimension must name at least one instance dimension")
  }
  if (!file.exists(arguments$point_cloud) ||
      !grepl("\\.la[sz]$", arguments$point_cloud, ignore.case = TRUE)) {
    stop("--point-cloud must be exactly one existing LAS/LAZ file")
  }
  # COPC accepts an optional non-COPC companion; computation still reads COPC.
  if (arguments$denoise != "none" && !is_copc_point_cloud(arguments$point_cloud)) {
    stop("--denoise requires a canonical COPC with OriginalPointIndex")
  }
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
  if (!is.null(arguments$aoi_geojson) && (!file.exists(arguments$aoi_geojson) ||
      !grepl("\\.(geojson|json)$", arguments$aoi_geojson, ignore.case = TRUE))) {
    stop("--aoi-geojson must be exactly one existing GeoJSON file")
  }
  if (!is.null(arguments$dataset_id) && (!is.finite(arguments$dataset_id) || arguments$dataset_id < 1)) {
    stop("--dataset-id must be a positive integer")
  }
  if (is.null(arguments$dataset_id)) arguments$dataset_id <- "foreststructure"
  for (directory in c(arguments$output_dir, arguments$temp_dir)) {
    if (!dir.exists(directory) || file.access(directory, mode = 2) != 0) {
      stop("output and temporary directories must exist and be writable")
    }
  }
  arguments$resources <- production_resources(arguments$memory_budget_gib)
  arguments$memory_budget_gib <- arguments$resources$memory_budget_gib
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
      shQuote(file.path(script_directory, "project_instance_dimension.py")),
      "--source", shQuote(source),
      "--output", shQuote(output),
      unlist(lapply(dimension, function(name) c("--dimension", shQuote(name))))
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output_lines, "status")
  if (!is.null(status) && status != 0L) {
    stop(sprintf(
      "streaming projection failed for %s: %s",
      paste(dimension, collapse = ", "),
      paste(output_lines, collapse = " ")
    ))
  }
  fromJSON(tail(output_lines, 1L), simplifyVector = FALSE)
}

expected_artifacts <- function(directory, dataset_id, dimensions, voxel_counts = FALSE,
                               denoising = FALSE) {
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
  c(common, dimension_files, if (voxel_counts) unlist(voxel_count_paths(directory, dataset_id)),
    if (denoising) paste0(prefix, c("_noise_points.csv", "_denoising.json")))
}

validate_performance_resources <- function(performance, memory_budget_gib,
                                           source_point_count,
                                           memory_limit_gib = memory_budget_gib) {
  if (nrow(performance) != 1L ||
      performance$point_count[[1]] != source_point_count ||
      !is.finite(performance$peak_rss_mib[[1]]) ||
      performance$peak_rss_mib[[1]] >= memory_budget_gib * 1024 ||
      !is.finite(performance$cgroup_memory_peak_mib[[1]]) ||
      # Kernel usage includes reclaimable disk cache. Compare that peak with the
      # allocation; the smaller internal budget governs working-set admission.
      performance$cgroup_memory_peak_mib[[1]] > memory_limit_gib * 1024 ||
      !is.finite(performance$temporary_disk_peak_mib[[1]])) {
    stop("incoming artifact validation failed: performance/resource checks failed")
  }
  invisible(TRUE)
}

validate_incoming_artifacts <- function(directory, dataset_id, dimensions,
                                        memory_budget_gib, source_point_count,
                                        point_cloud_name, sensor, country, voxel_counts = FALSE,
                                        memory_limit_gib = memory_budget_gib,
                                        denoising = FALSE) {
  matches_label <- function(values, expected) {
    if (identical(expected, "NA")) {
      return(all(is.na(values) | as.character(values) == "NA"))
    }
    all(!is.na(values) & as.character(values) == expected)
  }
  paths <- expected_artifacts(directory, dataset_id, dimensions, voxel_counts, denoising)
  if (voxel_counts) validate_voxel_count_artifacts(directory, dataset_id)
  absent <- paths[!file.exists(paths) | file.info(paths)$size <= 0]
  if (length(absent) > 0L) {
    stop(sprintf(
      "incoming artifact validation failed; missing/empty: %s",
      paste(basename(absent), collapse = ", ")
    ))
  }
  if (denoising) {
    status <- system2("python3", c(
      shQuote(file.path(script_directory, "denoising_validation.py")),
      shQuote(file.path(directory, paste0(dataset_id, "_denoising.json"))),
      shQuote(file.path(directory, paste0(dataset_id, "_noise_points.csv"))),
      "--point-count", format(source_point_count, scientific = FALSE, trim = TRUE)
    ))
    if (status != 0L) stop("incoming denoising accounting validation failed")
  }
  performance <- fread(
    file.path(directory, paste0(dataset_id, "_performance.csv"))
  )
  validate_performance_resources(performance, memory_budget_gib,
                                 source_point_count, memory_limit_gib)
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
  aoi_geojson <- if (is.null(arguments$aoi_geojson)) NULL else
    normalizePath(arguments$aoi_geojson, mustWork = TRUE)
  definitions <- extra_byte_definitions(point_cloud)
  source_definitions <- definitions
  available <- vapply(definitions, function(definition) definition$name, character(1))
  dimensions <- production_instance_dimensions(available, arguments$instance_dimension)
  arguments$instance_dimension <- dimensions
  missing <- character()
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
  copc_projection <- NULL
  required_dimensions <- unique(c(point_order_dimension, dimensions))
  if (requested_is_copc && any(match(required_dimensions, available) > 9L)) {
    # Pinned rlas truncates COPC ExtraBytes after field nine even with selector 0.
    # Reuse the verified bounded projection once, retaining the source-order key
    # and all selected IDs, then spatially index that compact scratch LAZ.
    projected <- file.path(job_directory, "analysis_projection.laz")
    copc_projection <- run_projection(staged_point_cloud, projected, required_dimensions)
    staged_point_cloud <- projected
    rlas::writelax(staged_point_cloud)
    if (!file.exists(sub("\\.laz$", ".lax", staged_point_cloud))) {
      stop("failed to spatially index high-ordinal COPC projection")
    }
    definitions <- extra_byte_definitions(staged_point_cloud)
    available <- vapply(definitions, function(definition) definition$name, character(1))
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

  thread_count <- arguments$resources$threads
  worker_request <- Sys.getenv("FORESTSTRUCTURE_CATALOG_WORKERS", "")
  catalog_worker_count <- validated_catalog_worker_count(
    if (nzchar(worker_request)) worker_request else as.character(min(10L, thread_count)),
    thread_count
  )
  message(sprintf(
    "Using validated Julia-memory-safe catalog worker count: %d",
    catalog_worker_count
  ))

  message("FORESTSTRUCTURE_STAGE stage=aoi_conversion status=started")
  generated_gpkg <- NULL
  if (is.null(aoi_geojson)) {
    footprint <- read_analysis_footprint(staged_point_cloud)
    aoi_provenance <- c(list(source = "point_cloud_extent", source_geojson = NULL),
                        geometry_summary(footprint$usable))
  } else {
    generated_gpkg <- file.path(job_directory, "analysis_aoi.gpkg")
    aoi_provenance <- convert_geojson_to_gpkg(
      aoi_geojson, generated_gpkg, point_cloud_xy_bounds(staged_point_cloud),
      tile_threads = thread_count
    )
    aoi_provenance$source <- "audit_aoi"
  }
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
  # High-ordinal COPC fields are now in the verified indexed scratch projection.

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
    # Plain LAZ retains its existing verified streaming projection fallback.
    if (ordinal <= 9L || requested_is_copc) {
      dimension_sources[[dimension]] <- staged_point_cloud
      selector <- catalog_selection_for_dimensions(staged_point_cloud,
        c(point_order_dimension, dimension), "xyz")
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

  segment_bucket_count <- suppressWarnings(as.integer(Sys.getenv(
    "FORESTSTRUCTURE_SEGMENT_BUCKET_COUNT",
    as.character(SEGMENT_BUCKET_COUNT)
  )))
  if (!is.finite(segment_bucket_count) ||
      segment_bucket_count < 1L ||
      segment_bucket_count > 65536L) {
    stop("FORESTSTRUCTURE_SEGMENT_BUCKET_COUNT must be an integer from 1 to 65536")
  }
  # Preserve Julia's other scientific constants; accepted normalized heights
  # have no ceiling. Metric denominators use the exact CHM-pass voxel maxima.
  parameters <- list(
    point_cloud = staged_point_cloud,
    point_cloud_sha256 = input_sha256,
    voxel_point_cloud = point_cloud,
    adaptive_spatial_input = requested_is_copc,
    export_voxel_counts = arguments$export_voxel_counts,
    voxel_count_chunk_points = arguments$voxel_count_chunk_points,
    dataset_id = arguments$dataset_id,
    aoi = generated_gpkg,
    output_dir = incoming_directory,
    temp_dir = job_directory,
    tile_size = 20,
    grid_search_step = 0.5,
    ptd_resolution = 20,
    dtm_resolution = 1,
    maximum_height = Inf,
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
    # Original COPC hierarchy also bounds projected reads: projection preserves
    # coordinates and records. Never estimate those from the flat LAZ header.
    memory_policy = resource_policy(point_cloud, arguments$memory_budget_gib, job_directory),
    sequential_instance_dimensions = FALSE,
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
  denoise_started <- proc.time()[["elapsed"]]
  parameters$denoising <- prepare_denoising_mask(
    staged_point_cloud, point_order_dimension, arguments,
    file.path(job_directory, "denoising"), catalog_worker_count,
    parameters$memory_policy
  )
  denoise_seconds <- if (is.null(parameters$denoising)) 0 else elapsed_seconds(denoise_started)
  denoise_provenance <- list(enabled = FALSE, method = "none")
  if (!is.null(parameters$denoising)) {
    parameters$denoising$summary$source_sha256 <- input_sha256
    parameters$denoising$summary$source_name <- basename(requested_point_cloud)
    parameters$denoising$summary$elapsed_seconds <- denoise_seconds
    parameters$denoising$summary$voxel_policy <- paste(
      "original source records; noise mask is not applied to voxel counts;",
      "normalization uses the DTM fitted after denoising"
    )
    denoise_provenance <- parameters$denoising$summary
    assert_memory_budget(parameters$denoising$lookup_reserve_bytes,
                         parameters, "denoising quality artifact export")
    write_denoising_artifacts(parameters$denoising,
      file.path(incoming_directory, as.character(arguments$dataset_id)))
    # A small adaptive read can intersect a large classification partition.
    # Reserve the largest complete mask partition for every downstream worker;
    # lookup releases each partition before opening the next one.
    parameters$memory_policy$task_fixed_bytes <-
      parameters$memory_policy$task_fixed_bytes + parameters$denoising$lookup_reserve_bytes
  }
  message("FORESTSTRUCTURE_STAGE stage=analysis status=started")
  run_analysis(parameters)
  message("FORESTSTRUCTURE_STAGE stage=analysis status=completed")
  performance <- fread(file.path(
    incoming_directory,
    paste0(arguments$dataset_id, "_performance.csv")
  ), colClasses = c(observed_maximum_voxel_height_m_exact = "character"))
  if (!is.null(parameters$denoising)) {
    performance[, denoise_seconds := denoise_seconds]
    performance[, total_seconds := total_seconds + denoise_seconds]
    fwrite(performance, file.path(incoming_directory,
      paste0(arguments$dataset_id, "_performance.csv")), na = "NA")
  }
  temporary_disk_peak_mib <- performance$temporary_disk_peak_mib[[1]] +
    temporary_directory_bytes(job_directory) / 1024^2

  # Persist enough provenance to prove that the run used COPC spatial streaming,
  # restored original order, and optionally validated the original companion.
  run_provenance <- list(
    status = "candidate",
    mode = if (is.null(parameters$denoising)) "julia-memory-safe" else "julia-memory-safe-local-sor",
    denoising = denoise_provenance,
    voxel_counts = list(enabled = arguments$export_voxel_counts,
      metadata = if (arguments$export_voxel_counts) paste0(arguments$dataset_id, "_voxel_counts.json") else NULL),
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
    source_crs_wkt = resolved_crs$wkt,
    raster_crs_wkt = terra::crs(terra::rast(file.path(
      incoming_directory, paste0(arguments$dataset_id, "_dtm.tif")))),
    raster_crs_policy = "horizontal CRS embedded in GeoTIFF; no coordinate or height transformation",
    # These fields make order restoration and streaming mode machine-auditable.
    point_order_dimension = point_order_dimension,
    catalog_spatial_streaming = if (!is.null(copc_projection)) "indexed LAZ projection" else
      if (requested_is_copc) "COPC" else "LAS/LAZ",
    copc_projection = copc_projection,
    requested_instance_dimensions = arguments$instance_dimension,
    processed_instance_dimensions = dimensions,
    instance_read_groups = instance_read_groups(
      dimensions, staged_point_cloud, dimension_sources,
      parameters$sequential_instance_dimensions
    ),
    missing_instance_dimensions = missing,
    extra_byte_definitions = definitions,
    source_extra_byte_definitions = source_definitions,
    selectors = selectors,
    projections = projections,
    catalog_workers = catalog_worker_count,
    catalog_worker_policy = if (nzchar(worker_request)) "explicit_cap" else "auto_up_to_ten",
    resource_admission = resource_plan_reports(parameters$memory_policy),
    data_table_threads = DATA_TABLE_THREADS,
    scheduler_resources = arguments$resources,
    footprint_source = aoi_provenance$source,
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
      # Read windows adapt; baseline raster footprints do not. Record both so
      # the scientific geometry cannot be mistaken for the memory policy.
      chm_raster_chunk_size_m = performance$chunk_size[[1]],
      chm_reduction = "maximum_voxel_height_per_exact_xy_within_baseline_chunk",
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
      maximum_height = NULL,
      normalized_height_policy = "finite_nonnegative_without_upper_cutoff",
      metric_height_policy = performance$metric_height_policy[[1]],
      metric_height_upper_m = performance$metric_height_upper_m[[1]],
      observed_maximum_voxel_height_m = performance$observed_maximum_voxel_height_m[[1]],
      observed_maximum_voxel_height_m_exact = performance$observed_maximum_voxel_height_m_exact[[1]],
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
    arguments$country,
    voxel_counts = arguments$export_voxel_counts,
    memory_limit_gib = if (is.null(arguments$resources$memory_limit_gib)) {
      arguments$memory_budget_gib
    } else arguments$resources$memory_limit_gib,
    denoising = !is.null(parameters$denoising)
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
