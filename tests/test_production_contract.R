source("/opt/foreststructure/production_contract.R")
source("/opt/foreststructure/aoi_conversion.R")
odd_path <- tempfile("input with spaces and ' quote ")
writeLines("fixture", odd_path)
stopifnot(grepl("^[0-9a-f]{64}$", sha256_file(odd_path)))
unlink(odd_path)
fails <- function(expr, text) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(text, conditionMessage(error), fixed = TRUE))
}
available <- c("OriginalPointIndex", "PredInstance_SAT")
stopifnot(identical(production_instance_dimensions(available), "PredInstance_SAT"))
fails(production_instance_dimensions(available, "PredInstance"), "missing requested")
fails(production_instance_dimensions(character()), "no supported")
fails(production_instance_dimensions(available, "junk"), "unsupported")
fails(production_instance_dimensions(available, ""), "empty")
root <- tempfile(); dir.create(root)
writeLines("4294967296", file.path(root, "memory.max"))
writeLines("200000 100000", file.path(root, "cpu.max"))
fake_env <- function(name, default) {
  values <- c(GALAXY_SLOTS = "10", GALAXY_MEMORY_MB = "8192", FORESTSTRUCTURE_THREADS = "20")
  if (name %in% names(values)) unname(values[[name]]) else default
}
resources <- production_resources(env = fake_env, cgroup_root = root)
stopifnot(resources$threads == 2L, resources$memory_budget_gib == 3.4,
          resources$memory_limit_gib == 4)
stopifnot(production_resources(2, fake_env, root)$memory_budget_gib == 2)
fails(production_resources(-1, fake_env, root), "must be positive")
bad_env <- function(name, default) if (name == "GALAXY_SLOTS") "2.5" else default
fails(production_resources(env = bad_env, cgroup_root = root), "integral")
unlink(root, recursive = TRUE)
cat("production contract unit checks passed\n")
