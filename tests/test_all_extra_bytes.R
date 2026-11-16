# The patched rlas reads every ExtraByte with selector "0"
# (patches/rlas-1.8.6-all-extra-bytes.patch). Unpatched it returned only the
# first nine, so an instance ID stored as the eleventh field (3462's
# PredInstance_RCT) was lost and the whole cloud had to be copied.
suppressPackageStartupMessages(library(lidR))

set.seed(7)
n <- 500L
las <- LAS(data.table::data.table(
  X = runif(n, 0, 10), Y = runif(n, 0, 10), Z = runif(n, 0, 5)
), header = LASheader(data.table::data.table(X = 0, Y = 0, Z = 0)))
names <- c("OriginalPointIndex", paste0("field_", 2:10), "PredInstance_RCT")
for (i in seq_along(names)) {
  las <- add_lasattribute(las, as.integer(seq_len(n) * i), names[[i]], names[[i]])
}
path <- tempfile(fileext = ".las")
writeLAS(las, path)

all_fields <- readLAS(path, select = "xyz0")
stopifnot(all(names %in% names(all_fields@data)))
stopifnot(identical(all_fields@data$PredInstance_RCT, as.integer(seq_len(n) * 11L)))
stopifnot(identical(all_fields@data$field_10, as.integer(seq_len(n) * 10L)))

# Digit selectors still read only the named fields.
first <- readLAS(path, select = "xyz1")
stopifnot("OriginalPointIndex" %in% names(first@data))
stopifnot(!any(c("field_2", "PredInstance_RCT") %in% names(first@data)))

# An ExtraByte unselect keeps the nine-field behaviour, as before the patch.
unselected <- readLAS(path, select = "xyz0 -2")
stopifnot(!"field_2" %in% names(unselected@data))

unlink(path)
cat("rlas selector 0 reads all", length(names), "ExtraBytes\n")
