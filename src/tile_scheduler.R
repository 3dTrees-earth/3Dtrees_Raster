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
    parallel_runner = parallel::mclapply,
    task_windows = NULL, memory_policy = NULL) {
  runner <- function(indices, calculate, workers) {
    run_bounded_tasks(indices, calculate, workers, parallel_runner)
  }
  rows <- if (is.null(memory_policy)) {
    runner(tile_indices, calculate_row, tile_workers)
  } else {
    run_memory_admitted_tasks(tile_indices, calculate_row, tile_workers,
      task_windows, memory_policy, "tiles", runner = runner)
  }
  validate_tile_rows(tile_indices, rows)
  if (length(rows)) {
    data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  } else if (exists("empty_result_table", mode = "function", inherits = TRUE)) {
    empty_result_table()
  } else {
    data.table::data.table()
  }
}

# Tiles read together aim at this many points per block.
tile_block_points <- function() {
  as.numeric(Sys.getenv("RASTER_TILE_BLOCK_POINTS", "2e7"))
}

# Group neighbouring tiles into square blocks of about tile_block_points()
# points at the mean density, but small enough that every worker gets at least
# two blocks. A block whose COPC upper-bound count exceeds four times the target
# (a dense patch) goes back to one read per tile.
plan_tile_blocks <- function(tile_windows, tile_size, density, memory_policy = NULL,
                             workers = 1L) {
  if (!length(tile_windows)) return(list())
  side <- if (is.finite(density) && density > 0) {
    max(1L, as.integer(floor(sqrt(tile_block_points() / (density * tile_size^2)))))
  } else {
    1L
  }
  xmin <- vapply(tile_windows, `[[`, numeric(1), 1L)
  ymin <- vapply(tile_windows, `[[`, numeric(1), 2L)
  tile_column <- round((xmin - min(xmin)) / tile_size)
  tile_row <- round((ymin - min(ymin)) / tile_size)
  block_count <- function(side) {
    length(unique(paste(floor(tile_row / side), floor(tile_column / side))))
  }
  while (side > 1L && block_count(side) < 2L * workers) side <- side - 1L
  if (side == 1L) return(as.list(seq_along(tile_windows)))
  column <- floor(tile_column / side)
  row <- floor(tile_row / side)
  key <- paste(row, column)
  blocks <- unname(split(seq_along(tile_windows), factor(key, levels = unique(key))))
  if (!is.null(memory_policy)) {
    windows <- lapply(blocks, function(block) tile_block_window(tile_windows[block]))
    counts <- task_memory_estimates(windows, memory_policy, "tiles")$counts
    dense <- counts > 4 * tile_block_points()
    blocks <- unlist(lapply(seq_along(blocks), function(i) {
      if (dense[[i]]) as.list(blocks[[i]]) else blocks[i]
    }), recursive = FALSE)
  }
  blocks
}

tile_block_window <- function(windows) {
  c(min(vapply(windows, `[[`, numeric(1), 1L)), min(vapply(windows, `[[`, numeric(1), 2L)),
    max(vapply(windows, `[[`, numeric(1), 3L)), max(vapply(windows, `[[`, numeric(1), 4L)))
}

# Run one task per block and return the rows in tile order. A failed block
# leaves its tiles without rows, so validation reports them.
collect_tile_block_rows <- function(
    tile_indices,
    blocks,
    calculate_block,
    workers,
    parallel_runner = parallel::mclapply,
    task_windows = NULL, memory_policy = NULL) {
  runner <- function(indices, calculate, workers) {
    run_bounded_tasks(indices, calculate, workers, parallel_runner)
  }
  block_indices <- seq_along(blocks)
  results <- if (is.null(memory_policy)) {
    runner(block_indices, calculate_block, workers)
  } else {
    run_memory_admitted_tasks(block_indices, calculate_block, workers,
      task_windows, memory_policy, "tiles", runner = runner)
  }
  rows <- vector("list", length(tile_indices))
  for (i in seq_along(blocks)) {
    result <- results[[i]]
    if (is.list(result) && !inherits(result, "try-error") &&
        length(result) == length(blocks[[i]])) {
      rows[match(blocks[[i]], tile_indices)] <- result
    }
  }
  validate_tile_rows(tile_indices, rows)
  if (length(rows)) {
    data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
  } else if (exists("empty_result_table", mode = "function", inherits = TRUE)) {
    empty_result_table()
  } else {
    data.table::data.table()
  }
}
