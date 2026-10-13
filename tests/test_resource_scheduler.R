# Exercise admission around the existing fork runner, including shrinking
# headroom, empty callback results, and fail-before-dispatch behavior.
local({
  root <- tempfile("resource-scheduler-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  path <- file.path(root, "points.las")
  las <- lidR::LAS(data.frame(X = as.numeric(seq_len(100)), Y = rep(1, 100), Z = rep(2, 100)))
  lidR::writeLAS(las, path)
  policy <- resource_policy(path, 1, root)
  policy$reserve_bytes <- 0
  windows <- rep(list(c(0, 0, 101, 2)), 5)
  estimates <- task_memory_estimates(windows, policy, "tiles")
  stopifnot(identical(estimates$mode, "header_fallback"), all(estimates$counts == 100))
  cost <- estimates$bytes[[1]]
  # First wave fits two, subsequent waves only one. The wrapper must recheck
  # aggregate memory rather than independently giving each worker the budget.
  dispatched <- list()
  headroom <- function() policy$budget_bytes - if (!length(dispatched)) 2.1 * cost else 1.1 * cost
  runner <- function(indices, calculate, workers) {
    stopifnot(sum(estimates$bytes[indices]) <= policy$budget_bytes - headroom(),
              workers == length(indices))
    dispatched[[length(dispatched) + 1L]] <<- indices
    run_bounded_tasks(indices, calculate, workers)
  }
  calculate <- function(i) if (i == 3L) NULL else i * 10L
  result <- run_memory_admitted_tasks(seq_len(5), calculate, 4L, windows, policy,
    memory_current = headroom, runner = runner)
  stopifnot(identical(result, list(10L, 20L, NULL, 40L, 50L)),
    identical(lengths(dispatched), c(2L, 1L, 1L, 1L)))
  report <- resource_plan_reports(policy)
  stopifnot(length(report) == 1L, length(report[[1]]$waves) == 4L)
  fails <- function(expression, pattern) {
    error <- tryCatch({ force(expression); NULL }, error = identity)
    stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error), fixed = TRUE))
  }
  fails(run_memory_admitted_tasks(1L, function(i) stop("must not run"), 4L,
    windows[1L], policy, memory_current = function() policy$budget_bytes),
    "memory_budget_exceeded")
  dispatched <- list()
  fails(run_memory_admitted_tasks(seq_len(5), function(i) {
    if (i == 1L) stop("worker fixture failure") else i
  }, 4L, windows, policy, memory_current = headroom, runner = runner), "worker failed")
  stopifnot(length(dispatched) == 1L)
  # Catalog task bounds must include the unchanged scientific DTM buffer.
  catalog <- lidR::readLAScatalog(path)
  lidR::opt_progress(catalog) <- FALSE
  lidR::opt_chunk_size(catalog) <- 60
  lidR::opt_chunk_buffer(catalog) <- 20
  chunks <- lidR:::engine_chunks(catalog, FALSE)
  stopifnot(length(chunks) >= 1L)
  for (chunk in chunks) {
    bounds <- as.numeric(chunk@bbox)
    buffered <- as.numeric(chunk@bbbox)
    stopifnot(all(buffered == bounds + c(-20, -20, 20, 20)))
    estimated <- catalog_task_window(chunk)
    stopifnot(all(estimated[1:2] <= buffered[1:2]), all(estimated[3:4] >= buffered[3:4]))
  }
  chunk <- chunks[[1]]
  chunk@filter <- "-inside -500 -600 700 800"
  stopifnot(identical(catalog_task_window(chunk), c(-500, -600, 700, 800)))
  grouped <- task_memory_estimates(list(windows[1:3], list(), windows[[1]]), policy, "CHM")
  stopifnot(identical(grouped$counts, c(100, 0, 100)),
            grouped$bytes[[1]] == estimates$bytes[[1]])
})
cat("aggregate resource admission checks passed\n")
