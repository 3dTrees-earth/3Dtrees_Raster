local({
  root <- tempfile("artifact_promotion_")
  incoming <- file.path(root, "incoming")
  output <- file.path(root, "output")
  dir.create(incoming, recursive = TRUE)
  dir.create(output)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)

  expected <- file.path(incoming, c("2001_dtm.tif", "2001_chm.tif"))
  stopifnot(all(file.create(expected)))
  sidecars <- file.path(incoming, c(
    "2001_dtm.tif.aux.xml",
    "2001_chm.tif.partial.tif.aux.xml"
  ))
  stopifnot(all(file.create(sidecars)))
  stale_sidecar <- file.path(output, "2001_dtm.tif.aux.xml")
  stopifnot(file.create(stale_sidecar))

  source_path <- file.path(output, "2001_original.laz")
  writeLines("preserve input", source_path)
  promoted <- promote_artifacts(expected, output, 2001)
  stopifnot(readLines(source_path) == "preserve input")

  stopifnot(identical(sort(basename(promoted)), sort(basename(expected))))
  stopifnot(all(file.exists(file.path(output, basename(expected)))))
  stopifnot(!any(file.exists(file.path(output, basename(sidecars)))))
  stopifnot(all(file.exists(sidecars)))
})

# A failed optional export must preserve every previously published artifact.
local({
  root <- tempfile("voxel_publication_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  old <- voxel_count_paths(root, 7)
  for (path in old) writeLines("previous valid product", path)
  original <- tools::md5sum(unlist(old))
  writer <- write_voxel_counts
  environment(writer) <- new.env(parent = environment(write_voxel_counts))
  environment(writer)$system2 <- function(command, args) {
    config <- jsonlite::read_json(file.path(root, "voxel_config.json"))
    writeLines("partial new mask", config$paths$mask)
    1L
  }
  parameters <- list(output_dir = root, dataset_id = 7, point_cloud = "unused")
  aoi <- list(usable = sf::st_as_sfc("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))", crs = 32632))
  error <- tryCatch(writer(parameters, "unused", aoi, root), error = identity)
  stopifnot(inherits(error, "error"),
            grepl("voxel point-count export failed", conditionMessage(error)),
            identical(original, tools::md5sum(unlist(old))),
            !length(list.files(root, pattern = "voxel_incoming", all.files = TRUE)))
  # Replacing only the optional quartet must preserve the existing DTM.
  incoming <- file.path(root, "incoming"); dir.create(incoming)
  writeLines("new counts", file.path(incoming, "7_voxel_counts.csv"))
  writeLines("terrain", file.path(root, "7_dtm.tif"))
  promote_artifacts(file.path(incoming, "7_voxel_counts.csv"), root, 7,
                    replace_dataset = FALSE)
  stopifnot(readLines(file.path(root, "7_dtm.tif")) == "terrain")
})

# The standard CLI also preserves the DTM and the entire previous generation
# when an error occurs after terrain generation (not just the voxel quartet).
local({
  root <- tempfile("standard_publication_")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  old <- file.path(root, c("7_dtm.tif", "7_voxel_counts.csv"))
  for (path in old) writeLines("previous generation", path)
  hashes <- tools::md5sum(old)
  runner <- run_standard_analysis
  environment(runner) <- new.env(parent = environment(run_standard_analysis))
  environment(runner)$run_analysis <- function(parameters) {
    writeLines("new terrain", artifact_paths(parameters$output_dir, 7)$dtm)
    stop("normalized height exceeds GeoTIFF band capacity")
  }
  error <- tryCatch(runner(list(output_dir=root, dataset_id=7)), error=identity)
  stopifnot(inherits(error, "error"), identical(hashes, tools::md5sum(old)),
            !length(list.files(root, pattern="incoming", all.files=TRUE)))
})

local({
  root <- tempfile("promotion_rollback_")
  incoming <- file.path(root, "incoming")
  dir.create(incoming, recursive=TRUE)
  on.exit(unlink(root, recursive=TRUE, force=TRUE), add=TRUE)
  names <- c("7_dtm.tif", "7_chm.tif")
  for (name in names) {
    writeLines("old", file.path(root, name))
    writeLines("new", file.path(incoming, name))
  }
  original <- tools::md5sum(file.path(root, names))
  promote <- promote_artifacts
  environment(promote) <- new.env(parent=environment(promote_artifacts))
  environment(promote)$file.rename <- function(from, to) {
    if (identical(from, file.path(incoming, names[[2]]))) return(FALSE)
    base::file.rename(from, to)
  }
  error <- tryCatch(promote(file.path(incoming,names), root,7),error=identity)
  stopifnot(inherits(error,"error"),
            identical(original,tools::md5sum(file.path(root,names))))
})
