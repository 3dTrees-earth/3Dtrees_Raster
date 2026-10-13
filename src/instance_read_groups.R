# Fuse only dimensions that really share an input file. Legacy high-ordinal LAZ
# projections each contain one dimension and must retain their independent scan.
instance_read_groups <- function(instance_dimensions, point_cloud,
                                 dimension_point_clouds = list(),
                                 sequential = FALSE) {
  groups <- list()
  canonical_sources <- character()
  for (dimension in instance_dimensions) {
    source <- dimension_point_clouds[[dimension]]
    if (is.null(source)) source <- point_cloud
    canonical_source <- normalizePath(source, mustWork = TRUE)
    index <- if (isTRUE(sequential)) NA_integer_ else
      match(canonical_source, canonical_sources)
    if (is.na(index)) {
      index <- length(groups) + 1L
      groups[[index]] <- list(source = source, dimensions = character())
      canonical_sources[[index]] <- canonical_source
    }
    groups[[index]]$dimensions <- c(groups[[index]]$dimensions, dimension)
  }
  groups
}

# A chunk writes each dimension to separate RDS files already. Carry their path
# ownership through catalog collection so reducing SAT never loads FM summaries.
# This metadata is bounded by chunk/bucket counts, not by point or voxel counts.
segment_store_for_dimension <- function(store, instance_dimension) {
  chunks <- store$dimension_chunks[[instance_dimension]]
  if (is.null(chunks)) {
    stop(sprintf("segment store has no path ownership for %s", instance_dimension))
  }
  list(chunks = chunks, bucket_count = store$bucket_count)
}
