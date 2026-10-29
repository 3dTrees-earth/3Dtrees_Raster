promote_artifacts <- function(incoming_artifacts, output_directory, dataset_id, replace_dataset = TRUE) {
  # Promote only paths returned by validate_incoming_artifacts(). GDAL may
  # leave non-scientific PAM sidecars beside TIFFs; those files are cleaned up
  # with the staging directory and never enter the published artifact set.
  incoming <- normalizePath(incoming_artifacts, mustWork = TRUE)
  prefix_pattern <- sprintf("^%s(_|\\.)", dataset_id)
  existing <- if (replace_dataset) list.files(
    output_directory,
    pattern = prefix_pattern,
    full.names = TRUE,
    recursive = FALSE,
    all.files = FALSE
  ) else {
    targets <- file.path(output_directory, basename(incoming))
    targets[file.exists(targets)]
  }
  # An output directory can also contain inputs or operator notes. Replace
  # only recognized generated products, never every dataset-prefixed file.
  if (replace_dataset) {
    generated_suffix <- paste0(
      "(_(dtm|chm|chm_filled|dsm|voxel_counts(_validity)?)\\.tif(\\.aux\\.xml)?|",
      "_(results|segment_diagnostics|performance|voxel_counts|noise_points)\\.csv|",
      "_tiles\\.(geojson|png)|",
      "_(aoi_conversion|julia_memory_safe_run|voxel_counts|denoising)\\.json)$"
    )
    existing <- existing[grepl(generated_suffix, basename(existing))]
  }
  backup_directory <- tempfile(
    pattern = sprintf(".backup_%s_", dataset_id),
    tmpdir = output_directory
  )
  dir.create(backup_directory)
  promoted <- character()
  success <- FALSE
  on.exit({
    restored <- TRUE
    if (!success) {
      unlink(promoted, force = TRUE)
      backup <- list.files(backup_directory, full.names = TRUE)
      for (path in backup) {
        restored <- file.rename(path, file.path(output_directory, basename(path))) && restored
      }
    }
    if (success || restored) {
      unlink(backup_directory, recursive = TRUE, force = TRUE)
    } else {
      warning(sprintf("artifact rollback incomplete; originals retained in %s",
                      backup_directory))
    }
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

