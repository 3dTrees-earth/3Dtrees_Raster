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
stopifnot(identical(production_instance_dimensions(available, "PredInstance_SAT"), "PredInstance_SAT"))
fails(production_instance_dimensions(available), "is required")
fails(production_instance_dimensions(available, NULL), "is required")
fails(production_instance_dimensions(available, "PredInstance"), "missing requested")
fails(production_instance_dimensions(character(), "PredInstance_SAT"), "missing requested")
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
writeLines("max", file.path(root, "memory.max"))
writeLines("max 100000", file.path(root, "cpu.max"))
large_env <- function(name, default) {
  values <- c(GALAXY_SLOTS = "32", GALAXY_MEMORY_MB = "95367")
  if (name %in% names(values)) unname(values[[name]]) else default
}
large <- production_resources(env = large_env, cgroup_root = root)
stopifnot(large$threads == 10L,
          large$memory_budget_gib == 95367 / 1024 * 0.85,
          large$memory_budget_gib > 70,
          production_resources(65, large_env, root)$memory_budget_gib == 65)
unlimited_env <- function(name, default) default
unallocated <- production_resources(env = unlimited_env, cgroup_root = root)
stopifnot(unallocated$threads == 1L, unallocated$memory_budget_gib == 60,
          is.null(unallocated$memory_limit_gib))
oversized_env <- function(name, default) {
  if (name == "GALAXY_MEMORY_MB") "204800" else default
}
capped <- production_resources(env = oversized_env, cgroup_root = root)
stopifnot(capped$memory_limit_gib * 1024^3 == 100e9,
          capped$memory_budget_gib == 100e9 / 1024^3 * 0.85)
fails(production_resources(100, large_env, root), "must not exceed 100 GB")
dir.create(file.path(root, "memory"))
v1_limit <- file.path(root, "memory/memory.limit_in_bytes")
writeLines("9223372036854771712", v1_limit)
v1_unlimited <- production_resources(env = unlimited_env, cgroup_root = root)
stopifnot(v1_unlimited$memory_budget_gib == 60,
          is.null(v1_unlimited$memory_limit_gib))
writeLines("21474836480", v1_limit)
v1_limited <- production_resources(env = large_env, cgroup_root = root)
stopifnot(v1_limited$memory_limit_gib == 20,
          v1_limited$memory_budget_gib == 17)
unlink(root, recursive = TRUE)
cat("production contract unit checks passed\n")
