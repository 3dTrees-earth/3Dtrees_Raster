# Scientific callbacks and task order remain owned by the existing schedulers.
# This module only reserves aggregate memory before dispatching bounded waves.
# Count only demonstrably clean inactive file pages as reclaimable. Inactive
# lists may still contain dirty/writeback pages; subtracting their global totals
# from the inactive count is conservative even when some belong to active lists.
# Keep anonymous, tmpfs/shared, dirty, writeback and kernel memory in admission.
# Counters are not an atomic snapshot, so leave the scheduler/container headroom.
resource_memory_snapshot <- function(root = "/sys/fs/cgroup", rss = current_rss_bytes()) {
  candidates <- c(file.path(root, "memory.current"),
                  file.path(root, "memory/memory.usage_in_bytes"),
                  file.path(root, "memory.usage_in_bytes"))
  paths <- candidates[file.exists(candidates)]
  if (!length(paths)) {
    return(list(source = "process_rss_fallback", raw_bytes = rss,
                reclaimable_bytes = 0, working_set_bytes = rss,
                process_rss_bytes = rss))
  }
  path <- paths[[1L]]
  version_two <- identical(basename(path), "memory.current")
  raw <- suppressWarnings(as.numeric(readLines(path, n = 1L, warn = FALSE)))
  if (length(raw) != 1L || !is.finite(raw) || raw < 0) {
    stop("cannot read current container memory for task admission")
  }
  reclaimable <- 0
  stats_path <- file.path(dirname(path), "memory.stat")
  if (file.exists(stats_path)) {
    lines <- strsplit(trimws(readLines(stats_path, warn = FALSE)), "[[:space:]]+")
    if (length(lines) && all(lengths(lines) == 2L)) {
      names <- vapply(lines, `[[`, character(1), 1L)
      values <- suppressWarnings(as.numeric(vapply(lines, `[[`, character(1), 2L)))
      stats <- setNames(values, names)
      fields <- if (version_two) {
        c("inactive_file", "file", "file_dirty", "file_writeback", "shmem")
      } else {
        # v1 usage includes descendants; never subtract local-only counters.
        c("total_inactive_file", "total_cache", "total_dirty", "total_writeback", "total_shmem")
      }
      if (!anyDuplicated(names) && all(fields %in% names) &&
          all(is.finite(stats[fields]) & stats[fields] >= 0)) {
        selected <- unname(stats[fields])
        reclaimable <- min(raw, max(0, min(selected[1:2]) - sum(selected[3:5])))
      }
    }
  }
  list(source = if (version_two) "cgroup_v2" else "cgroup_v1",
       raw_bytes = raw, reclaimable_bytes = reclaimable,
       working_set_bytes = max(rss, raw - reclaimable), process_rss_bytes = rss)
}

resource_memory_current <- function(root = "/sys/fs/cgroup") {
  resource_memory_snapshot(root)$working_set_bytes
}

resource_policy <- function(source, budget_gib, work_directory,
                            helper = file.path(script_directory, "copc_task_counts.py")) {
  if (!is.finite(budget_gib) || budget_gib <= 0) {
    stop("memory-aware scheduling requires a finite positive memory budget")
  }
  directory <- file.path(work_directory, "resource_plans")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  list(source = source, budget_bytes = budget_gib * 1024^3,
       reserve_bytes = 512 * 1024^2, directory = directory, helper = helper,
       # DTM releases its source ordinal before PTD/TIN and retains no voxel
       # or segment tables. Its separate coefficient preserves admission of
       # known dense terrain windows; see the historical 1805 regression.
       # Estimates include temporary allocations but are not allocator bounds.
       bytes_per_point = c(DTM = 256, CHM = 512, segment = 512, tiles = 512,
                           denoise = 768),
       task_fixed_bytes = 128 * 1024^2)
}

task_memory_estimates <- function(windows, policy, phase) {
  if (!length(windows)) return(list(bytes = numeric(), counts = numeric(), mode = "empty"))
  prefix <- tempfile(paste0(tolower(phase), "_"), tmpdir = policy$directory)
  input <- paste0(prefix, "_windows.json")
  output <- paste0(prefix, "_counts.json")
  # Adaptive CHM tasks read children sequentially while retaining only a compact
  # XY-maxima table. Reserve the largest child, not the sum of their raw points.
  grouped <- lapply(windows, function(window) if (is.list(window)) window else list(window))
  flat_windows <- unlist(grouped, recursive = FALSE)
  if (!length(flat_windows)) {
    return(list(bytes = rep(policy$task_fixed_bytes, length(windows)),
                counts = rep(0, length(windows)), mode = "empty", prefix = prefix))
  }
  jsonlite::write_json(unname(flat_windows), input, auto_unbox = FALSE, digits = NA)
  status <- system2("python3", c(shQuote(policy$helper), "--input", shQuote(policy$source),
    "--windows", shQuote(input), "--output", shQuote(output), "--allow-header-fallback"))
  if (status != 0L || !file.exists(output)) stop("could not estimate spatial task memory")
  estimate <- jsonlite::read_json(output, simplifyVector = TRUE)
  flat_counts <- as.numeric(estimate$point_count_upper_bounds)
  if (length(flat_counts) != length(flat_windows) || any(!is.finite(flat_counts) | flat_counts < 0)) {
    stop("invalid spatial point-count estimates")
  }
  ends <- cumsum(lengths(grouped))
  starts <- c(0L, head(ends, -1L))
  counts <- vapply(seq_along(grouped), function(i) {
    if (ends[[i]] == starts[[i]]) 0 else max(flat_counts[seq.int(starts[[i]] + 1L, ends[[i]])])
  }, numeric(1))
  coefficient <- unname(policy$bytes_per_point[[phase]])
  if (is.null(coefficient) || !is.finite(coefficient)) stop("unknown memory planning phase")
  list(bytes = policy$task_fixed_bytes + counts * coefficient, counts = counts,
       mode = estimate$mode, prefix = prefix)
}

admitted_task_wave <- function(costs, remaining, workers, available_bytes) {
  if (!length(remaining)) return(integer())
  if (any(!is.finite(costs) | costs < 0) || !is.finite(available_bytes)) {
    stop("invalid aggregate task memory estimate")
  }
  candidates <- head(remaining, max(1L, as.integer(workers)))
  fits <- cumsum(costs[candidates]) <= available_bytes
  if (!any(fits)) {
    stop(sprintf(paste0("FORESTSTRUCTURE_FLAG memory_budget_exceeded: one spatial task ",
      "requires an estimated %.3f GiB with %.3f GiB available; geometry was not changed"),
      costs[[candidates[[1]]]] / 1024^3, max(0, available_bytes) / 1024^3))
  }
  candidates[fits]
}

run_memory_admitted_tasks <- function(indices, calculate, workers, windows,
                                      policy = NULL, phase = "tiles",
                                      memory_current = resource_memory_current,
                                      runner = run_bounded_tasks) {
  if (is.null(policy)) return(runner(indices, calculate, workers))
  if (length(indices) != length(windows)) stop("task windows do not match scheduled tasks")
  if (!length(indices)) return(list())
  estimates <- task_memory_estimates(windows, policy, phase)
  message(sprintf("RESOURCE_PLAN phase=%s requested_workers=%d tasks=%d max_task_gib=%.3f mode=%s",
    phase, workers, length(indices), max(estimates$bytes) / 1024^3, estimates$mode))
  remaining <- seq_along(indices)
  results <- vector("list", length(indices))
  waves <- list()
  while (length(remaining)) {
    gc(verbose = FALSE)
    current <- memory_current()
    snapshot <- resource_memory_snapshot()
    available <- policy$budget_bytes - policy$reserve_bytes - current
    wave <- admitted_task_wave(estimates$bytes, remaining, workers, available)
    reserved <- sum(estimates$bytes[wave])
    message(sprintf("RESOURCE_WAVE phase=%s workers=%d reserved_gib=%.3f current_gib=%.3f budget_gib=%.3f",
      phase, length(wave), reserved / 1024^3, current / 1024^3, policy$budget_bytes / 1024^3))
    started <- proc.time()[["elapsed"]]
    batch <- runner(indices[wave], calculate, length(wave))
    if (length(batch) != length(wave) || any(vapply(batch, inherits, logical(1), "try-error"))) {
      stop(sprintf("%s memory-admitted worker failed or returned incomplete results", phase))
    }
    results[wave] <- batch
    waves[[length(waves) + 1L]] <- list(tasks = unname(indices[wave]), workers = length(wave),
      reserved_bytes = reserved, memory_before_bytes = current,
      memory_before = snapshot, memory_after = resource_memory_snapshot(),
      memory_after_bytes = memory_current(), elapsed_seconds = elapsed_seconds(started))
    remaining <- remaining[-seq_along(wave)]
  }
  jsonlite::write_json(list(phase = phase, requested_workers = workers,
    budget_bytes = policy$budget_bytes, reserve_bytes = policy$reserve_bytes,
    bytes_per_point = unname(policy$bytes_per_point[[phase]]),
    task_fixed_bytes = policy$task_fixed_bytes,
    estimate_mode = estimates$mode, point_count_upper_bounds = estimates$counts,
    task_reservations_bytes = estimates$bytes, waves = waves),
    paste0(estimates$prefix, "_admission.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
  results
}

resource_plan_reports <- function(policy) {
  if (is.null(policy)) return(list())
  paths <- sort(list.files(policy$directory, pattern = "_admission[.]json$", full.names = TRUE))
  lapply(paths, jsonlite::read_json, simplifyVector = FALSE)
}
