# Stable result order and one task per child bound each worker's object lifetime.
run_bounded_tasks <- function(indices, calculate, workers,
                              parallel_runner = parallel::mclapply) {
  if (length(indices) > 0L && workers > 1L) {
    parallel_runner(indices, calculate, mc.cores = workers,
                    mc.preschedule = FALSE, mc.set.seed = FALSE)
  } else {
    lapply(indices, calculate)
  }
}

tile_worker_plan <- function(requested_threads, tile_count, requested_workers = 1L) {
  requested_threads <- max(1L, as.integer(requested_threads))
  workers <- max(1L, min(as.integer(requested_workers), requested_threads,
                         as.integer(tile_count)))
  list(
    process_workers = workers,
    # Native OpenMP pools inherited across fork can deadlock when restarted
    # with multiple threads. Forked tasks must stay single-threaded, as catalog
    # callbacks already do; a sequential process may use the full allowance.
    threads_per_worker = if (workers > 1L) 1L else requested_threads
  )
}

validate_tile_rows <- function(tile_indices, rows) {
  valid_rows <- sum(vapply(
    rows,
    function(row) {
      !is.null(row) && !inherits(row, "try-error") && NROW(row) == 1L
    },
    logical(1)
  ))
  if (length(rows) != length(tile_indices) || valid_rows != length(tile_indices)) {
    stop(sprintf(
      "tile worker results incomplete: expected %d rows, received %d valid rows",
      length(tile_indices),
      valid_rows
    ))
  }
  invisible(rows)
}

collect_tile_rows <- function(
    tile_indices,
    calculate_row,
    tile_workers,
    parallel_runner = parallel::mclapply) {
  rows <- run_bounded_tasks(tile_indices, calculate_row, tile_workers,
                            parallel_runner)
  validate_tile_rows(tile_indices, rows)
  if (length(rows)) {
    data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  } else if (exists("empty_result_table", mode = "function", inherits = TRUE)) {
    empty_result_table()
  } else {
    data.table::data.table()
  }
}
