# Conservative local SOR decisions are computed once in fixed scientific windows.
# Adaptive analysis windows and worker counts never define the SOR population.
# Only source ordinals are persisted; original coordinates and classes stay intact.
denoising_defaults <- function() {
  list(denoise = "sor", denoise_neighbors = 10L, denoise_sigma = 12,
       denoise_window_size = 60, denoise_buffer = 20)
}

add_denoising_arguments <- function(parser) {
  defaults <- denoising_defaults()
  parser$add_argument("--denoise", choices = c("none", "sor"), default = defaults$denoise,
    help = "Conservative local lidR SOR before terrain/metrics (default: sor); raw voxel counts remain unfiltered")
  parser$add_argument("--denoise-neighbors", type = "integer", default = defaults$denoise_neighbors,
    help = sprintf("SOR nearest-neighbor count (default: %s)", defaults$denoise_neighbors))
  parser$add_argument("--denoise-sigma", type = "double", default = defaults$denoise_sigma,
    help = sprintf("SOR distance mean plus this many standard deviations (default: %s)", defaults$denoise_sigma))
  parser$add_argument("--denoise-window-size", type = "double", default = defaults$denoise_window_size,
    help = sprintf("Fixed scientific SOR core width in metres (default: %s)", defaults$denoise_window_size))
  parser$add_argument("--denoise-buffer", type = "double", default = defaults$denoise_buffer,
    help = sprintf("Fixed scientific SOR context beyond each core in metres (default: %s)", defaults$denoise_buffer))
}

validate_denoising_options <- function(arguments) {
  defaults <- denoising_defaults()
  settings <- lapply(names(defaults), function(name) parameter_or(arguments, name, defaults[[name]]))
  names(settings) <- names(defaults)
  if (length(settings$denoise) != 1L || is.na(settings$denoise) ||
      !settings$denoise %in% c("none", "sor")) stop("--denoise must be none or sor")
  for (name in names(settings)[-1L]) {
    value <- settings[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0) {
      stop(sprintf("--%s must be finite and positive", gsub("_", "-", name)))
    }
  }
  if (settings$denoise_neighbors != floor(settings$denoise_neighbors) ||
      settings$denoise_neighbors >= .Machine$integer.max) {
    stop("--denoise-neighbors must be a positive supported integer")
  }
  settings
}

denoising_partition_name <- function(ix, iy) {
  paste0("core_", format(ix, scientific = FALSE, trim = TRUE), "_",
         format(iy, scientific = FALSE, trim = TRUE), ".rds")
}

denoising_core_membership <- function(x, y, ix, iy, size) {
  # Use the same global index arithmetic at classification and lookup time,
  # including when a user chooses a noninteger scientific window width.
  floor(x / size) == ix & floor(y / size) == iy
}

denoising_ids <- function(las, dimension) {
  if (is.null(dimension) || length(dimension) != 1L ||
      !dimension %in% names(las@data)) stop("denoising requires the source point-order dimension")
  ids <- las@data[[dimension]]
  if (!is.numeric(ids) || any(!is.finite(ids) | ids < 0 | ids != floor(ids) | ids > 2^53)) {
    stop("denoising requires exact finite nonnegative source ordinals")
  }
  ids
}

# Isolate lidR's temporary Classification column from any existing class labels.
# Calling multiple classify_noise passes directly on source classes would reset
# existing class-18 labels. This module records only its independent SOR decision.
denoising_sor_decisions <- function(las, neighbors, sigma) {
  count <- lidR::npoints(las)
  if (count <= neighbors) {
    return(list(rejected = rep(FALSE, count), insufficient_context = TRUE))
  }
  classified <- las
  classified@data <- data.table::copy(las@data)
  if ("Classification" %in% names(classified@data)) {
    data.table::set(classified@data, j = "Classification", value = NULL)
  }
  classified <- lidR::classify_noise(classified, lidR::sor(k = neighbors, m = sigma))
  list(rejected = classified@data[["Classification"]] == 18L,
       insufficient_context = FALSE)
}

denoising_chunk <- function(chunk, point_order_dimension, settings, directory) {
  configure_catalog_worker()
  if (!chunk_has_valid_extent(chunk)) return(NULL)
  bounds <- as.numeric(chunk@bbox)
  ix <- round(bounds[[1L]] / settings$denoise_window_size)
  iy <- round(bounds[[2L]] / settings$denoise_window_size)
  filename <- denoising_partition_name(ix, iy)
  result <- list(ix = ix, iy = iy, filename = filename, core_points = 0,
       context_points = 0, rejected_points = 0, accepted_points = 0,
       insufficient_context_points = 0)
  # Read the catalog's exact buffered filter directly. readLAS(LAScluster)
  # discards an entire read when its floating bbox test sees only buffer points;
  # that test can disagree with global index ownership at decimal boundaries.
  las <- lidR::readLAS(chunk@files, select = chunk@select, filter = chunk@filter)
  if (lidR::is.empty(las)) return(result)
  result$context_points <- lidR::npoints(las)
  las <- restore_original_point_order(las, point_order_dimension, "denoising")
  ids <- denoising_ids(las, point_order_dimension)
  if (any(!is.finite(las@data$X) | !is.finite(las@data$Y) | !is.finite(las@data$Z))) {
    stop("denoising encountered nonfinite source coordinates")
  }
  core <- denoising_core_membership(las@data$X, las@data$Y, ix, iy,
                                   settings$denoise_window_size)
  if (!any(core)) return(result)
  decisions <- denoising_sor_decisions(las, settings$denoise_neighbors, settings$denoise_sigma)
  rejected <- ids[core & decisions$rejected]
  if (length(rejected)) saveRDS(rejected, file.path(directory, filename), compress = FALSE)
  list(ix = ix, iy = iy, filename = filename, core_points = sum(core),
       context_points = lidR::npoints(las), rejected_points = length(rejected),
       accepted_points = sum(core) - length(rejected),
       insufficient_context_points = if (decisions$insufficient_context) sum(core) else 0)
}

prepare_denoising_mask <- function(point_cloud, point_order_dimension, options,
                                   work_directory, workers, memory_policy = NULL) {
  settings <- validate_denoising_options(options)
  if (settings$denoise == "none") return(NULL)
  if (is.null(point_order_dimension) || length(point_order_dimension) != 1L ||
      !nzchar(point_order_dimension)) stop("denoising requires OriginalPointIndex")
  if (utils::packageVersion("lidR") != "4.3.2") {
    stop("denoising requires the validated lidR 4.3.2 runtime")
  }
  input_points <- point_count_from_las_header(point_cloud)
  assert_lidr_point_count_supported(input_points)
  message("FORESTSTRUCTURE_STAGE stage=denoise status=started")
  if (dir.exists(work_directory)) stop("denoising work directory must be new")
  dir.create(work_directory, recursive = TRUE)
  if (!dir.exists(work_directory)) stop("cannot create denoising work directory")
  catalog <- lidR::readLAScatalog(point_cloud)
  bounds <- point_cloud_xyz_bounds(point_cloud)[c("xmin", "ymin", "xmax", "ymax")]
  grid_shape <- floor(bounds[3:4] / settings$denoise_window_size) -
    floor(bounds[1:2] / settings$denoise_window_size) + 1
  # The shared hierarchy-admission helper accepts at most 100,000 windows.
  # Reject an extreme footprint before lidR allocates its rectangular grid.
  if (any(!is.finite(grid_shape)) || any(grid_shape < 1) || prod(grid_shape) > 100000) {
    stop("denoising fixed scientific grid exceeds the supported 100000-window planning limit")
  }
  # lidR's catalog index uses strict rectangle overlap and would omit the core
  # starting exactly at a file's maximum XY. Expand only catalog lookup bounds
  # by a representable epsilon, so such boundary points get their own half-open
  # core. Scientific grid coordinates and actual point coordinates do not move.
  catalog <- catalog_with_source_boundaries(catalog, point_cloud)
  lidR::opt_chunk_size(catalog) <- settings$denoise_window_size
  lidR::opt_chunk_buffer(catalog) <- settings$denoise_buffer
  lidR::opt_chunk_alignment(catalog) <- c(0, 0)
  lidR::opt_select(catalog) <- catalog_selection_for_dimensions(point_cloud, point_order_dimension, "xyz")
  lidR::opt_progress(catalog) <- FALSE
  pieces <- apply_catalog_chunks(catalog, denoising_chunk, workers, "denoise",
    point_order_dimension = point_order_dimension, settings = settings,
    directory = work_directory, memory_policy = memory_policy,
    .options = list(automerge = FALSE, drop_null = FALSE, need_buffer = TRUE))
  manifest <- data.table::rbindlist(pieces, fill = TRUE)
  if (!nrow(manifest) || sum(manifest$core_points) != input_points ||
      anyDuplicated(manifest[, .(ix, iy)])) {
    stop(sprintf("denoising core ownership does not cover the source exactly once: expected %.0f, observed %.0f",
                 input_points, sum(manifest$core_points)))
  }
  data.table::setorderv(manifest, c("iy", "ix"))
  rejected_points <- sum(manifest$rejected_points)
  maximum <- max(manifest$rejected_points)
  summary <- list(schema_version = 1L, enabled = TRUE, algorithm = "lidR::sor",
    method = "fixed_window_local_sor", point_order_dimension = point_order_dimension,
    lidR_version = as.character(utils::packageVersion("lidR")),
    k = settings$denoise_neighbors, m = settings$denoise_sigma,
    window_size_m = settings$denoise_window_size, buffer_m = settings$denoise_buffer,
    origin = c(0, 0), boundary_ownership = "half-open: lower inclusive, upper exclusive",
    statistics_population = "all source points in each fixed buffered scientific window",
    existing_classification_policy = "preserved; not used as a rejection criterion",
    insufficient_context_policy = "retain points when buffered population is no larger than k",
    point_count = input_points, rejected_points = rejected_points,
    retained_points = input_points - rejected_points,
    insufficient_context_points = sum(manifest$insufficient_context_points),
    partition_count = nrow(manifest),
    raw_voxel_counts_filtered = FALSE,
    coordinate_changes = FALSE)
  jsonlite::write_json(summary, file.path(work_directory, "summary.json"),
                       auto_unbox = TRUE, pretty = TRUE, digits = NA)
  message(sprintf("FORESTSTRUCTURE_STAGE stage=denoise status=completed points=%.0f rejected=%.0f",
                  input_points, rejected_points))
  list(directory = work_directory, core_size = settings$denoise_window_size,
       origin = c(0, 0), point_order_dimension = point_order_dimension,
       partitions = manifest, summary = summary, max_partition_points = maximum,
       # Include decompression, hash/match work and transient R vector copies.
       # Callers add this bound per simultaneous downstream worker reservation.
       lookup_reserve_bytes = 16 * 1024^2 + 128 * maximum)
}

read_denoising_partition <- function(mask, row) {
  ids <- readRDS(file.path(mask$directory, row$filename))
  if (!is.numeric(ids) || length(ids) != row$rejected_points ||
      any(!is.finite(ids) | ids < 0 | ids != floor(ids) | ids > 2^53) ||
      is.unsorted(ids, strictly = TRUE)) stop("invalid denoising mask partition")
  ids
}

apply_denoising_mask <- function(las, mask) {
  if (is.null(mask) || lidR::is.empty(las)) return(las)
  ids <- denoising_ids(las, mask$point_order_dimension)
  xmin <- floor(min(las@data$X) / mask$core_size)
  xmax <- floor(max(las@data$X) / mask$core_size)
  ymin <- floor(min(las@data$Y) / mask$core_size)
  ymax <- floor(max(las@data$Y) / mask$core_size)
  rows <- which(mask$partitions$ix >= xmin & mask$partitions$ix <= xmax &
                mask$partitions$iy >= ymin & mask$partitions$iy <= ymax &
                mask$partitions$rejected_points > 0)
  if (!length(rows)) return(las)
  keep <- rep(TRUE, length(ids))
  for (i in rows) {
    rejected <- read_denoising_partition(mask, mask$partitions[i, ])
    keep[ids %in% rejected] <- FALSE
    rm(rejected)
    # A partition must die before loading the next; admission uses their max,
    # not their sum, even for small adaptive windows touching several cores.
    gc(verbose = FALSE)
  }
  lidR::filter_poi(las, keep)
}

# A bounded, inspectable quality product; each ID occurs once and identifies a
# source observation, without moving/deleting a record in the immutable COPC.
write_denoising_artifacts <- function(mask, prefix) {
  if (is.null(mask)) return(character())
  csv <- paste0(prefix, "_noise_points.csv")
  metadata <- paste0(prefix, "_denoising.json")
  writeLines("OriginalPointIndex", csv)
  for (i in which(mask$partitions$rejected_points > 0)) {
    ids <- read_denoising_partition(mask, mask$partitions[i, ])
    data.table::fwrite(data.table::data.table(OriginalPointIndex = ids), csv,
                      append = TRUE, col.names = FALSE, scipen = 999)
    rm(ids)
    gc(verbose = FALSE)
  }
  jsonlite::write_json(c(mask$summary, list(partitions = as.data.frame(mask$partitions))),
                      metadata, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  c(csv, metadata)
}
