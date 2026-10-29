# Source run.R (FORESTSTRUCTURE_SOURCE_ONLY=1) and denoising.R before this file.
local({
  directory <- tempfile("denoising-test-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  settings <- list(denoise = "none", denoise_neighbors = 10L, denoise_sigma = 12,
                   denoise_window_size = 60, denoise_buffer = 20)
  parser <- argparse::ArgumentParser(add_help = FALSE)
  add_denoising_arguments(parser)
  stopifnot(identical(validate_denoising_options(list()), settings),
            isTRUE(all.equal(validate_denoising_options(parser$parse_args(character())), settings)))
  rejects <- function(expression) inherits(tryCatch(force(expression), error = identity), "error")
  stopifnot(is.null(prepare_denoising_mask("unused", NULL, list(denoise = "none"), directory, 1L)),
            is.null(prepare_denoising_mask("unused", NULL, list(), directory, 1L)),
            rejects(prepare_denoising_mask("unused", NULL, list(denoise = "sor"), directory, 1L)))
  for (field in names(settings)[-1L]) {
    for (value in list(NA_real_, Inf, 0, -1, "bad", numeric())) {
      invalid <- settings
      invalid[[field]] <- value
      stopifnot(rejects(validate_denoising_options(invalid)))
    }
  }
  invalid <- settings
  invalid$denoise_neighbors <- 1.5
  stopifnot(rejects(validate_denoising_options(invalid)),
    rejects(validate_denoising_options(list(denoise = "ivf"))),
    rejects(prepare_denoising_mask("unused", NULL, modifyList(settings, list(denoise = "sor")),
                                   directory, 1L)))
  stopifnot(identical(denoising_core_membership(c(-1, 0, 59.99, 60),
    c(1, 0, 59.99, 1), 0, 0, 60), c(FALSE, TRUE, TRUE, FALSE)))

  # Dense ordinary structure, a supported tall crown, isolated high/low noise,
  # exact scientific boundaries, negative grid coordinates, and sparse context.
  ordinary <- expand.grid(X = seq(1, 21, by = 1), Y = seq(1, 21, by = 1))
  ordinary$Z <- 20 + sin(ordinary$X) * 0.05
  crown <- expand.grid(X = seq(10, 12, by = 0.2), Y = seq(10, 12, by = 0.2))
  crown$Z <- 80 + sin(crown$X) * 0.05
  boundary <- data.frame(X = c(-60, 0, 60, 120), Y = c(0, 0, 60, 120), Z = 20)
  noise <- data.frame(X = c(10, 11), Y = c(10, 11), Z = c(500, -500))
  points <- rbind(ordinary, crown, boundary, noise)
  source <- lidR::LAS(data.frame(points, Classification = rep(18L, nrow(points))))
  source <- lidR::add_lasattribute(source, as.numeric(seq_len(nrow(points)) - 1L),
                                  "OriginalPointIndex", "source order")
  source <- lidR::add_lasattribute(source, rep(17, nrow(points)), "PredInstance", "tree ID")
  before <- data.table::copy(source@data)
  direct <- denoising_sor_decisions(source, settings$denoise_neighbors, settings$denoise_sigma)
  stopifnot(identical(source@data, before),
            all(tail(direct$rejected, 2)),
            sum(direct$rejected) == 2,
            !any(direct$rejected[nrow(ordinary) + seq_len(nrow(crown))]))
  tiny <- lidR::filter_poi(source, OriginalPointIndex < 5)
  sparse <- denoising_sor_decisions(tiny, settings$denoise_neighbors, settings$denoise_sigma)
  stopifnot(sparse$insufficient_context, !any(sparse$rejected))

  path <- file.path(directory, "source.las")
  lidR::writeLAS(source, path)
  # SOR with omitted options must invoke the same conservative classifier as
  # explicit options; only the two injected extremes may be removed here.
  first <- prepare_denoising_mask(path, "OriginalPointIndex", list(denoise = "sor"),
                                  file.path(directory, "one"), 1L)
  second <- prepare_denoising_mask(path, "OriginalPointIndex",
                                   modifyList(settings, list(denoise = "sor")),
                                   file.path(directory, "two"), 2L)
  stopifnot(identical(first$summary, second$summary),
            identical(first$partitions, second$partitions),
            first$summary$point_count == nrow(points),
            first$summary$rejected_points == 2,
            first$summary$retained_points + first$summary$rejected_points == nrow(points),
            first$summary$insufficient_context_points > 0,
            first$lookup_reserve_bytes >= 128 * first$max_partition_points)
  for (i in which(first$partitions$rejected_points > 0)) {
    stopifnot(identical(read_denoising_partition(first, first$partitions[i, ]),
                        read_denoising_partition(second, second$partitions[i, ])))
  }
  filtered <- apply_denoising_mask(source, first)
  high_low <- tail(source@data$OriginalPointIndex, 2)
  stopifnot(!any(high_low %in% filtered@data$OriginalPointIndex),
            identical(filtered@data$OriginalPointIndex, head(source@data$OriginalPointIndex, -2L)),
            all(filtered@data$PredInstance == 17),
            all(filtered@data$Classification == 18L),
            identical(source@data, before),
            all((nrow(ordinary) + seq_len(nrow(crown)) - 1L) %in%
                  filtered@data$OriginalPointIndex))
  # Downstream execution reads can vary freely after the fixed mask is saved.
  source_ids <- source@data$OriginalPointIndex
  pieces <- lapply(split(source_ids, floor(source@data$X / 7)), function(ids) {
    subset <- lidR::filter_poi(source, OriginalPointIndex %in% ids)
    apply_denoising_mask(subset, first)@data$OriginalPointIndex
  })
  stopifnot(identical(sort(unlist(pieces, use.names = FALSE)),
                      sort(filtered@data$OriginalPointIndex)))
  shuffled <- lidR::filter_poi(source, rep(TRUE, lidR::npoints(source)))
  shuffled@data <- shuffled@data[rev(seq_len(nrow(shuffled@data)))]
  stopifnot(identical(sort(apply_denoising_mask(shuffled, first)@data$OriginalPointIndex),
                      sort(filtered@data$OriginalPointIndex)))
  products <- write_denoising_artifacts(first, file.path(directory, "quality"))
  csv <- data.table::fread(products[[1]])
  metadata <- jsonlite::read_json(products[[2]], simplifyVector = TRUE)
  stopifnot(identical(names(csv), "OriginalPointIndex"),
            nrow(csv) == first$summary$rejected_points,
            metadata$point_count == nrow(points), !anyDuplicated(csv$OriginalPointIndex),
            !any(csv$OriginalPointIndex %in% filtered@data$OriginalPointIndex))
  decimal <- lidR::LAS(expand.grid(X = c(-0.6, 0, 0.3, 0.6), Y = c(0, 0.3, 0.6), Z = 20))
  decimal <- lidR::add_lasattribute(decimal, as.numeric(seq_len(lidR::npoints(decimal)) - 1),
                                    "OriginalPointIndex", "source order")
  decimal_path <- file.path(directory, "decimal.las")
  lidR::writeLAS(decimal, decimal_path)
  decimal_settings <- modifyList(settings, list(denoise = "sor"))
  decimal_settings$denoise_window_size <- 0.3
  decimal_settings$denoise_buffer <- 0.1
  decimal_mask <- prepare_denoising_mask(decimal_path, "OriginalPointIndex", decimal_settings,
                                         file.path(directory, "decimal-mask"), 1L)
  stopifnot(decimal_mask$summary$point_count == lidR::npoints(decimal),
            sum(decimal_mask$partitions$core_points) == lidR::npoints(decimal))
  too_many_windows <- modifyList(settings, list(denoise = "sor"))
  too_many_windows$denoise_window_size <- 0.001
  stopifnot(rejects(prepare_denoising_mask(decimal_path, "OriginalPointIndex", too_many_windows,
                                          file.path(directory, "wide-mask"), 1L)))
  broken <- first
  bad_row <- which(broken$partitions$rejected_points > 0)[[1L]]
  saveRDS(c(1, 1), file.path(broken$directory, broken$partitions$filename[[bad_row]]))
  stopifnot(rejects(apply_denoising_mask(source, broken)))
})
cat("fixed-window local SOR and reusable mask checks passed\n")
