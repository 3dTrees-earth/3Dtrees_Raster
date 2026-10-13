# Run after sourcing run_julia_memory_safe.R with SOURCE_ONLY enabled. The
# publication gate distinguishes parent working memory from whole-cgroup cache.
local({
  performance <- data.table(point_count = 12, peak_rss_mib = 2 * 1024,
    cgroup_memory_peak_mib = 3.75 * 1024, temporary_disk_peak_mib = 100)
  stopifnot(isTRUE(validate_performance_resources(performance, 3.4, 12, 4)))
  check_failure <- function(record, limit = 4, points = 12) {
    error <- tryCatch(validate_performance_resources(record, 3.4, points, limit),
                      error = identity)
    stopifnot(inherits(error, "error"), grepl("performance/resource checks failed",
                                             conditionMessage(error), fixed = TRUE))
  }
  # The allocation itself is allowed; any measured excess still fails.
  boundary <- copy(performance)
  boundary$cgroup_memory_peak_mib <- 4 * 1024
  stopifnot(isTRUE(validate_performance_resources(boundary, 3.4, 12, 4)))
  boundary$cgroup_memory_peak_mib <- 4 * 1024 + 0.001
  check_failure(boundary)
  boundary <- copy(performance)
  boundary$peak_rss_mib <- 3.4 * 1024
  check_failure(boundary)
  # Missing resource evidence must not be treated as zero usage.
  for (column in c("peak_rss_mib", "cgroup_memory_peak_mib", "temporary_disk_peak_mib")) {
    invalid <- copy(performance)
    invalid[[column]] <- NA_real_
    check_failure(invalid)
  }
  check_failure(performance, points = 13)
  check_failure(rbind(performance, performance))
  # Without a known allocation, the existing conservative budget remains the
  # fallback for the whole-cgroup peak too.
  error <- tryCatch(validate_performance_resources(performance, 3.4, 12), error = identity)
  stopifnot(inherits(error, "error"))
})
cat("publication resource boundaries passed\n")
