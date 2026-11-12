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
  # Two tasks fit at first, later only one. The pool must count the running
  # tasks' estimates rather than giving each worker the budget.
  headroom <- function() policy$budget_bytes - 2.1 * cost
  calculate <- function(i) if (i == 3L) NULL else i * 10L
  result <- run_memory_admitted_tasks(seq_len(5), calculate, 4L, windows, policy,
    memory_current = headroom)
  stopifnot(identical(result, list(10L, 20L, NULL, 40L, 50L)))
  report <- resource_plan_reports(policy)
  admissions <- report[[1]]$admissions
  stopifnot(length(report) == 1L, length(admissions) == 5L,
            max(vapply(admissions, `[[`, integer(1), "workers")) == 2L,
            all(vapply(admissions, function(a) a$reserved_bytes <= 2.1 * cost, logical(1))),
            report[[1]]$bytes_per_point == 128,
            report[[1]]$task_fixed_bytes == 128 * 1024^2)
  fails <- function(expression, pattern) {
    error <- tryCatch({ force(expression); NULL }, error = identity)
    stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error), fixed = TRUE))
  }
  fails(run_memory_admitted_tasks(1L, function(i) stop("must not run"), 4L,
    windows[1L], policy, memory_current = function() policy$budget_bytes),
    "memory_budget_exceeded")
  fails(run_memory_admitted_tasks(seq_len(5), function(i) {
    if (i == 1L) stop("worker fixture failure") else i
  }, 4L, windows, policy, memory_current = headroom), "worker failed")
  # A worker that dies returns nothing, which is not a NULL result.
  fails(run_memory_admitted_tasks(seq_len(2), function(i) {
    if (i == 2L) tools::pskill(Sys.getpid(), tools::SIGKILL) else i
  }, 2L, windows[1:2], policy, memory_current = headroom), "worker failed")
  # One task at a time runs in the parent process.
  stopifnot(identical(run_memory_admitted_tasks(seq_len(3), function(i) i, 1L,
    windows[1:3], policy, memory_current = headroom), list(1L, 2L, 3L)))
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
            grouped$bytes[[1]] == policy$task_fixed_bytes + 100 * 192)
  # The previously validated 1805 COPC's largest buffered DTM window has this
  # hierarchy upper bound and a 27.38 GiB whole-run peak, which its estimate covers.
  dense_policy <- resource_policy(path, 59.5, root)
  dense_count <- 185339618
  dense_cost <- dense_policy$task_fixed_bytes +
    dense_count * unname(dense_policy$bytes_per_point[["DTM"]])
  retained_parent <- 512 * 1024^2
  available <- function(budget_gib) budget_gib * 1024^3 -
    dense_policy$reserve_bytes - retained_parent
  stopifnot(dense_cost / 1024^3 > 27.38,
            2 * dense_cost <= available(59.5), 3 * dense_cost > available(59.5),
            dense_cost > available(25.5),
            identical(unname(dense_policy$bytes_per_point[c("CHM", "segment", "tiles")]),
                      c(192, 256, 128)))
  dtm_estimates <- task_memory_estimates(windows, dense_policy, "DTM")
  stopifnot(all(dtm_estimates$bytes == dense_policy$task_fixed_bytes + 100 * 160))
})
cat("aggregate resource admission checks passed\n")
