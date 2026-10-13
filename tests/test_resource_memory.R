# Unit fixtures for aggregate admission accounting; no live cgroup is modified.
local({
  root <- tempfile("resource-memory-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  gib <- 1024^3
  stat <- function(inactive = 18, file = 20, dirty = 2, writeback = 1, shmem = 3) {
    values <- c(inactive_file = inactive, file = file, file_dirty = dirty,
                file_writeback = writeback, shmem = shmem) * gib
    writeLines(paste(names(values), format(values, scientific = FALSE, trim = TRUE)),
               file.path(root, "memory.stat"))
  }
  writeLines(as.character(24 * gib), file.path(root, "memory.current"))
  stat()
  measured <- resource_memory_snapshot(root, rss = gib)
  stopifnot(identical(measured$source, "cgroup_v2"),
            measured$raw_bytes == 24 * gib,
            measured$reclaimable_bytes == 12 * gib,
            measured$working_set_bytes == 12 * gib)
  # A naïve min(inactive_file, file - dirty) would discard dirty inactive pages.
  stat(inactive = 10, file = 20, dirty = 8, writeback = 0, shmem = 0)
  measured <- resource_memory_snapshot(root, rss = gib)
  stopifnot(measured$reclaimable_bytes == 2 * gib)
  # Dirty, writeback and tmpfs can consume the entire inactive allowance.
  stat(inactive = 3, file = 20, dirty = 2, writeback = 1, shmem = 3)
  stopifnot(resource_memory_snapshot(root, rss = gib)$reclaimable_bytes == 0)
  # The RSS floor protects resident parent mappings and counter skew.
  stat(inactive = 24, file = 24, dirty = 0, writeback = 0, shmem = 0)
  stopifnot(resource_memory_snapshot(root, rss = 4 * gib)$working_set_bytes == 4 * gib)
  # Missing or unusable statistics grant no speculative cache deduction.
  for (lines in list(character(), "inactive_file 1000", "file NaN", "broken",
                     c("file 1000", "file 1000"))) {
    writeLines(lines, file.path(root, "memory.stat"))
    stopifnot(resource_memory_snapshot(root, rss = gib)$working_set_bytes == 24 * gib)
  }
  unlink(file.path(root, "memory.stat"))
  stopifnot(resource_memory_snapshot(root, rss = gib)$working_set_bytes == 24 * gib)
  writeLines("not-a-counter", file.path(root, "memory.current"))
  error <- tryCatch(resource_memory_snapshot(root, rss = gib), error = identity)
  stopifnot(inherits(error, "error"))
  unlink(file.path(root, "memory.current"))
  fallback <- resource_memory_snapshot(root, rss = gib)
  stopifnot(identical(fallback$source, "process_rss_fallback"), fallback$raw_bytes == gib)
  # cgroup v1 requires hierarchical total_* statistics to match parent usage.
  dir.create(file.path(root, "memory"))
  writeLines(as.character(24 * gib), file.path(root, "memory/memory.usage_in_bytes"))
  values <- c(total_inactive_file = 18, total_cache = 20, total_dirty = 2,
              total_writeback = 1, total_shmem = 3) * gib
  writeLines(paste(names(values), format(values, scientific = FALSE, trim = TRUE)),
             file.path(root, "memory/memory.stat"))
  measured <- resource_memory_snapshot(root, rss = gib)
  stopifnot(identical(measured$source, "cgroup_v1"), measured$working_set_bytes == 12 * gib)
  writeLines("inactive_file 9999999999", file.path(root, "memory/memory.stat"))
  stopifnot(resource_memory_snapshot(root, rss = gib)$reclaimable_bytes == 0)
})
cat("resource working-set accounting checks passed\n")
