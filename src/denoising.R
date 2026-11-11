# kNN distance outlier filter (Julian Frey's proposal): a point is noise when its
# mean distance to its k nearest neighbours exceeds the mean plus `sigma`
# standard deviations of that distance over the whole point cloud. knn_denoising.py
# computes it in one multithreaded pass; the rejected source ordinals are kept
# per fixed XY window so later reads load only the windows they overlap.
# Coordinates and classes stay untouched.
denoising_defaults <- function() {
  list(denoise = "knn", denoise_neighbors = 8L, denoise_sigma = 2,
       denoise_window_size = 60)
}

add_denoising_arguments <- function(parser) {
  defaults <- denoising_defaults()
  parser$add_argument("--denoise", choices = c("none", "knn"), default = defaults$denoise,
    help = "kNN distance outlier filter before terrain/metrics (default: knn; none switches it off); raw voxel counts remain unfiltered")
  parser$add_argument("--denoise-neighbors", type = "integer", default = defaults$denoise_neighbors,
    help = sprintf("Nearest neighbours per point, counting the point itself as lidR::sor() does (default: %s)", defaults$denoise_neighbors))
  parser$add_argument("--denoise-sigma", type = "double", default = defaults$denoise_sigma,
    help = sprintf("Reject above the global mean distance plus this many standard deviations (default: %s)", defaults$denoise_sigma))
  parser$add_argument("--denoise-window-size", type = "double", default = defaults$denoise_window_size,
    help = sprintf("Width in metres of the windows the noise list is stored in; does not change results (default: %s)", defaults$denoise_window_size))
}

validate_denoising_options <- function(arguments) {
  defaults <- denoising_defaults()
  settings <- lapply(names(defaults), function(name) parameter_or(arguments, name, defaults[[name]]))
  names(settings) <- names(defaults)
  if (length(settings$denoise) != 1L || is.na(settings$denoise) ||
      !settings$denoise %in% c("none", "knn")) stop("--denoise must be none or knn")
  for (name in names(settings)[-1L]) {
    value <- settings[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0) {
      stop(sprintf("--%s must be finite and positive", gsub("_", "-", name)))
    }
  }
  if (settings$denoise_neighbors != floor(settings$denoise_neighbors) ||
      settings$denoise_neighbors < 2 || settings$denoise_neighbors >= .Machine$integer.max) {
    stop("--denoise-neighbors must be an integer of at least 2 (it counts the point itself)")
  }
  settings
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

prepare_denoising_mask <- function(point_cloud, point_order_dimension, options,
                                   work_directory, workers, memory_policy = NULL) {
  settings <- validate_denoising_options(options)
  if (settings$denoise == "none") return(NULL)
  if (is.null(point_order_dimension) || length(point_order_dimension) != 1L ||
      !nzchar(point_order_dimension)) stop("denoising requires OriginalPointIndex")
  input_points <- point_count_from_las_header(point_cloud)
  assert_lidr_point_count_supported(input_points)
  message("RASTER_STAGE stage=denoise status=started")
  if (dir.exists(work_directory)) stop("denoising work directory must be new")
  dir.create(work_directory, recursive = TRUE)
  if (!dir.exists(work_directory)) stop("cannot create denoising work directory")
  bounds <- point_cloud_xyz_bounds(point_cloud)[c("xmin", "ymin", "xmax", "ymax")]
  grid_shape <- floor(bounds[3:4] / settings$denoise_window_size) -
    floor(bounds[1:2] / settings$denoise_window_size) + 1
  if (any(!is.finite(grid_shape)) || any(grid_shape < 1) || prod(grid_shape) > 100000) {
    stop("denoising fixed scientific grid exceeds the supported 100000-window planning limit")
  }
  budget <- if (is.null(memory_policy)) 8 * 1024^3 else
    memory_policy$budget_bytes - memory_policy$reserve_bytes - current_rss_bytes()
  status <- system2("python3", c(shQuote(file.path(script_directory, "knn_denoising.py")),
    "--input", shQuote(point_cloud), "--order-dimension", shQuote(point_order_dimension),
    "--neighbors", settings$denoise_neighbors, "--sigma", format(settings$denoise_sigma, digits = 17),
    "--window-size", format(settings$denoise_window_size, digits = 17),
    "--threads", workers, "--memory-budget-bytes", format(budget, scientific = FALSE),
    "--output-dir", shQuote(work_directory)))
  if (!identical(status, 0L)) stop("kNN denoising failed")
  result <- jsonlite::read_json(file.path(work_directory, "manifest.json"), simplifyVector = TRUE)
  manifest <- data.table::as.data.table(result$partitions)
  if (!nrow(manifest) || result$point_count != input_points ||
      sum(manifest$core_points) != input_points || anyDuplicated(manifest[, .(ix, iy)])) {
    stop("denoising core ownership does not cover the source exactly once")
  }
  data.table::setorderv(manifest, c("iy", "ix"))
  rejected_points <- sum(manifest$rejected_points)
  maximum <- max(manifest$rejected_points)
  summary <- list(schema_version = 3L, enabled = TRUE,
    algorithm = "mean distance to the k nearest neighbours, the point itself counted as in lidR::sor (scipy cKDTree)",
    method = "global_knn_distance", point_order_dimension = point_order_dimension,
    k = settings$denoise_neighbors, sigma = settings$denoise_sigma,
    distance_mean_m = result$distance_mean_m, distance_sd_m = result$distance_sd_m,
    threshold_m = result$threshold_m,
    window_size_m = settings$denoise_window_size,
    origin = c(0, 0), boundary_ownership = "half-open: lower inclusive, upper exclusive",
    statistics_population = "all source points (whole point cloud)",
    neighbour_population = "all source points (whole point cloud)",
    existing_classification_policy = "preserved; not used as a rejection criterion",
    insufficient_context_policy = "retain all points when the cloud holds no more than k points",
    point_count = input_points, rejected_points = rejected_points,
    retained_points = input_points - rejected_points,
    insufficient_context_points = sum(manifest$insufficient_context_points),
    partition_count = nrow(manifest),
    raw_voxel_counts_filtered = FALSE,
    coordinate_changes = FALSE)
  jsonlite::write_json(summary, file.path(work_directory, "summary.json"),
                       auto_unbox = TRUE, pretty = TRUE, digits = NA)
  message(sprintf("RASTER_STAGE stage=denoise status=completed points=%.0f rejected=%.0f",
                  input_points, rejected_points))
  list(directory = work_directory, core_size = settings$denoise_window_size,
       origin = c(0, 0), point_order_dimension = point_order_dimension,
       partitions = manifest, summary = summary, max_partition_points = maximum,
       # Include decompression, hash/match work and transient R vector copies.
       # Callers add this bound per simultaneous downstream worker reservation.
       lookup_reserve_bytes = 16 * 1024^2 + 128 * maximum)
}

read_denoising_partition <- function(mask, row) {
  path <- file.path(mask$directory, row$filename)
  ids <- readBin(path, "double", n = file.size(path) / 8, size = 8, endian = "little")
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
