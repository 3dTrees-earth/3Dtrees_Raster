# The direct scheduler must preserve lidR's chunk layout, order, results and
# caller state, and reject failed workers instead of returning partial output.
local({
  root <- tempfile("catalog_scheduler_")
  dir.create(root)
  on.exit(unlink(root, recursive=TRUE, force=TRUE), add=TRUE)
  path <- file.path(root, "source.las")
  xyz <- expand.grid(X=c(.25,1.25,2.25), Y=c(.25,1.25,2.25))
  xyz$Z <- as.double(seq_len(nrow(xyz)))
  lidR::writeLAS(lidR::LAS(xyz), path)
  catalog <- lidR::readLAScatalog(path)
  lidR::opt_chunk_size(catalog) <- 1
  lidR::opt_chunk_buffer(catalog) <- 0
  lidR::opt_progress(catalog) <- FALSE
  summarize <- function(chunk) {
    configure_catalog_worker()
    las <- lidR::readLAS(chunk)
    if (lidR::is.empty(las)) return(NULL)
    list(n=nrow(las@data), z=sum(las@data$Z))
  }
  options <- list(automerge=FALSE, drop_null=FALSE)
  expected <- with_catalog_workers(1L,"reference",
    lidR::catalog_apply(catalog, summarize, .options=options))
  previous_threads <- lidR::get_lidr_threads()
  on.exit(lidR::set_lidr_threads(previous_threads), add=TRUE)
  lidR::set_lidr_threads(2L)
  previous_options <- getOption("lidR.no.future")
  previous_verbose <- options(lidR.verbose=TRUE)
  on.exit(options(previous_verbose), add=TRUE)
  one <- apply_catalog_chunks(catalog,summarize,1L,"test",.options=options)
  two <- apply_catalog_chunks(catalog,summarize,2L,"test",.options=options)
  stopifnot(identical(expected,one),identical(expected,two),
            lidR::get_lidr_threads()==2L,
            identical(getOption("lidR.no.future"),previous_options),
            isTRUE(getOption("lidR.verbose")))
  fail <- function(chunk) {
    las <- lidR::readLAS(chunk)
    if (lidR::is.empty(las)) return(NULL)
    stop("injected worker error")
  }
  for (workers in c(1L,2L)) {
    error <- suppressWarnings(tryCatch(
      apply_catalog_chunks(catalog,fail,workers,"test failure",.options=options),
      error=identity))
    stopifnot(inherits(error,"error"),
              grepl(if (workers == 1L) "injected worker error" else
                      "failed or returned incomplete",conditionMessage(error)),
              isTRUE(getOption("lidR.verbose")),
              lidR::get_lidr_threads()==2L)
  }
})
