source(file.path("src", "tile_scheduler.R"))

calculate_row <- function(index) {
  data.frame(tile_id = index)
}

error <- tryCatch(
  {
    validate_tile_rows(1:3, lapply(1:2, calculate_row))
    NULL
  },
  error = identity
)

if (is.null(error) ||
    !grepl("tile worker", conditionMessage(error), fixed = TRUE)) {
  stop("incomplete parallel tile results were not rejected")
}

plan <- tile_worker_plan(requested_threads = 10L, tile_count = 34L)
stopifnot(
  identical(plan$process_workers, 1L),
  identical(plan$threads_per_worker, 10L)
)

plan <- tile_worker_plan(10L, 34L, 10L)
stopifnot(plan$process_workers == 10L, plan$threads_per_worker == 1L)
plan <- tile_worker_plan(10L, 2L, 10L)
stopifnot(plan$process_workers == 2L, plan$threads_per_worker == 5L)
serial <- collect_tile_rows(1:6, calculate_row, 1L)
parallel <- collect_tile_rows(1:6, calculate_row, 2L)
stopifnot(identical(serial, parallel))
