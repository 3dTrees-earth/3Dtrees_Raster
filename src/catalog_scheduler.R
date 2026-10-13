# Keep lidR's chunk construction, callbacks, raster writers and ordering. Only
# dispatch changes: scanning future globals for every tiny chunk dominated the
# sparse-cloud profile. Disk paths / compact R summaries are safe to send back
# from the same bounded fork runner already used for tile metrics.
catalog_task_window <- function(chunk) {
  bounds <- as.numeric(chunk@bbbox)
  # Include expansions caused by lidR's numeric-to-text LASlib filter bounds.
  fields <- strsplit(trimws(chunk@filter), "[[:space:]]+")[[1L]]
  start <- match("-inside", fields)
  if (is.na(start) || start + 4L > length(fields)) {
    stop("memory admission requires rectangular catalog filters")
  }
  actual <- suppressWarnings(as.numeric(fields[start + seq_len(4L)]))
  if (any(!is.finite(actual))) stop("invalid catalog read bounds")
  c(pmin(bounds[1:2], actual[1:2]), pmax(bounds[3:4], actual[3:4]))
}

apply_catalog_chunks <- function(catalog, callback, workers, phase, ...,
                                  .options = NULL, memory_policy = NULL,
                                  task_windows_fn = NULL) {
  workers <- max(1L, as.integer(workers))
  # The internal engine is pinned and tested at this version. Preserve lidR's
  # public path on other versions or platforms where safe forking is disabled.
  if (utils::packageVersion("lidR") != "4.3.2" ||
      (workers > 1L && !parallelly::supportsMulticore())) {
    if (!is.null(memory_policy)) {
      stop("memory-admitted catalog scheduling requires pinned lidR 4.3.2 and fork support")
    }
    return(with_catalog_workers(workers, phase,
      lidR::catalog_apply(catalog, callback, ..., .options = .options)))
  }
  previous_threads <- lidR::get_lidr_threads()
  on.exit(lidR::set_lidr_threads(previous_threads), add = TRUE)
  settings <- lidR:::engine_parse_options(.options)
  if (settings$autoread || settings$autocrop ||
      !identical(settings$automerge, FALSE) || settings$drop_null) {
    stop("direct catalog scheduler requires explicit callbacks and unmerged results")
  }
  lidR:::assert_fun_is_null_with_empty_cluster(catalog, callback, ...)
  lidR:::assert_processing_constraints_are_repected(
    catalog, settings$need_buffer, settings$need_output_file)
  alignment <- if (settings$check_alignment) settings$raster_alignment else FALSE
  chunks <- lidR:::engine_chunks(catalog, alignment)
  message(sprintf("Processing %s LAScatalog chunks with %d worker(s) [direct bounded scheduler]",
                  phase, workers))
  process <- function(index) {
    previous_options <- options(lidR.no.future = TRUE, lidR.progress = FALSE,
                                lidR.verbose = FALSE)
    on.exit(options(previous_options), add = TRUE)
    result <- lidR:::engine_apply(
      chunks[index], callback, catalog@processing_options, catalog@output_options,
      .GLOBALS = settings$globals, ...)[[1L]]
    if (inherits(result, c("SpatRaster", "SpatVector"))) {
      stop("catalog callbacks must return disk paths or serializable summaries")
    }
    result
  }
  windows <- if (is.null(memory_policy)) {
    NULL
  } else if (is.null(task_windows_fn)) {
    lapply(chunks, catalog_task_window)
  } else {
    lapply(chunks, task_windows_fn)
  }
  results <- run_memory_admitted_tasks(seq_along(chunks), process, workers,
    windows, memory_policy, phase)
  if (length(results) != length(chunks)) stop("catalog worker results incomplete")
  assert_catalog_completed(results, phase)
  results
}
