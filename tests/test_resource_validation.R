# Run after sourcing raster.R with SOURCE_ONLY enabled. The
# publication gate distinguishes parent working memory from whole-cgroup cache.
local({
  performance <- data.table(point_count = 12, peak_rss_mib = 2 * 1024,
    cgroup_memory_peak_mib = 3.75 * 1024, temporary_disk_peak_mib = 100)
  stopifnot(isTRUE(validate_performance_resources(performance, 3.4, 12, 4)))
  # The report never fails a run; an incomplete one is only a warning.
  check_failure <- function(record, limit = 4, points = 12) {
    stopifnot(isTRUE(validate_performance_resources(record, 3.4, points, limit)))
  }
  # A cgroup peak above the allocation is disk cache and only warns (dataset 97);
  # working memory above the budget still fails.
  boundary <- copy(performance)
  boundary$cgroup_memory_peak_mib <- 4 * 1024
  stopifnot(isTRUE(validate_performance_resources(boundary, 3.4, 12, 4)))
  boundary$cgroup_memory_peak_mib <- 19 * 1024
  stopifnot(isTRUE(validate_performance_resources(boundary, 3.4, 12, 4)))
  # Galaxy enforces the grant; the debug report only warns on working memory too.
  boundary <- copy(performance)
  boundary$peak_rss_mib <- 3.4 * 1024
  stopifnot(isTRUE(validate_performance_resources(boundary, 3.4, 12, 4)))
  # Missing resource evidence must not be treated as zero usage.
  for (column in c("peak_rss_mib", "cgroup_memory_peak_mib", "temporary_disk_peak_mib")) {
    invalid <- copy(performance)
    invalid[[column]] <- NA_real_
    check_failure(invalid)
  }
  check_failure(performance, points = 13)
  check_failure(rbind(performance, performance))
  # Without a known allocation the budget is the reference for the warning.
  stopifnot(isTRUE(validate_performance_resources(performance, 3.4, 12)))
})
cat("publication resource boundaries passed\n")
