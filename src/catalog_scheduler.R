# Keep lidR's chunk construction, callbacks, raster writers and ordering. Only
# dispatch changes: scanning future globals for every tiny chunk dominated the
# sparse-cloud profile. Disk paths / compact R summaries are safe to send back
# from the same bounded fork runner already used for tile metrics.
apply_catalog_chunks <- function(catalog, callback, workers, phase, ...,
                                  .options = NULL) {
  workers <- max(1L, as.integer(workers))
  # The internal engine is pinned and tested at this version. Preserve lidR's
  # public path on other versions or platforms where safe forking is disabled.
  if (utils::packageVersion("lidR") != "4.3.2" ||
      (workers > 1L && !parallelly::supportsMulticore())) {
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
  results <- run_bounded_tasks(seq_along(chunks), process, workers)
  if (length(results) != length(chunks)) stop("catalog worker results incomplete")
  assert_catalog_completed(results, phase)
  results
}
