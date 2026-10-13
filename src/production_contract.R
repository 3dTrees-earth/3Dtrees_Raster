# Resolve scheduler allocations once, before any point-cloud work. Galaxy owns
# job limits; the analysis leaves headroom inside the smaller declared/kernel cap.
positive_resource <- function(value, label, integer = FALSE) {
  number <- suppressWarnings(as.numeric(value))
  if (length(number) != 1L || !is.finite(number) || number <= 0 ||
      (integer && number != floor(number))) {
    stop(paste(label, "must be positive", if (integer) "and integral"))
  }
  number
}

production_resources <- function(memory_budget = NULL, env = Sys.getenv,
                                 cgroup_root = "/sys/fs/cgroup") {
  declared <- function(name, integer = FALSE) {
    value <- env(name, "")
    if (!nzchar(value)) return(Inf)
    positive_resource(value, name, integer)
  }
  read_limit <- function(name) {
    path <- file.path(cgroup_root, name)
    if (!file.exists(path)) return(Inf)
    value <- suppressWarnings(as.numeric(readLines(path, warn = FALSE, n = 1L)))
    # cgroup v1 exposes its page-rounded signed-long maximum for unlimited
    # memory, unlike v2's literal "max". It is not a scheduler allocation.
    if (grepl("memory.limit_in_bytes$", name) && length(value) == 1L &&
        is.finite(value) && value >= 2^60) return(Inf)
    if (length(value) == 1L && is.finite(value) && value > 0) value else Inf
  }
  threads <- min(declared("GALAXY_SLOTS", TRUE), declared("FORESTSTRUCTURE_THREADS", TRUE))
  cpu_path <- file.path(cgroup_root, "cpu.max")
  if (file.exists(cpu_path)) {
    fields <- strsplit(readLines(cpu_path, warn = FALSE, n = 1L), " +")[[1]]
    cpu <- suppressWarnings(as.numeric(fields))
    if (length(cpu) == 2L && all(is.finite(cpu)) && all(cpu > 0)) {
      threads <- min(threads, max(1, floor(cpu[[1]] / cpu[[2]])))
    }
  }
  if (!is.finite(threads)) threads <- 1
  threads <- min(10, threads)
  memory_limit <- min(declared("GALAXY_MEMORY_MB") / 1024,
                      read_limit("memory.max") / 1024^3,
                      read_limit("memory/memory.limit_in_bytes") / 1024^3)
  # Galaxy chooses a size-dependent allocation. Use its available headroom
  # rather than retaining the former 60 GiB cap on larger job destinations.
  # The deployment ceiling is 100 decimal GB; CLI values remain binary GiB.
  maximum_memory_gib <- 100e9 / 1024^3
  requested <- if (is.null(memory_budget)) {
    if (is.finite(memory_limit)) Inf else 60
  } else
    positive_resource(memory_budget, "--memory-budget-gib")
  if (!is.null(memory_budget) && requested > maximum_memory_gib) {
    stop("--memory-budget-gib must not exceed 100 GB (93.132 GiB)")
  }
  list(threads = as.integer(threads),
       memory_budget_gib = min(requested, min(memory_limit, maximum_memory_gib) * 0.85),
       memory_limit_gib = if (is.finite(memory_limit)) min(memory_limit, maximum_memory_gib) else NULL)
}

# Auto mode discovers supported SmartTile attributes. Explicit requests are
# strict, and even auto mode requires at least one segmentation dimension.
production_instance_dimensions <- function(available, requested = NULL) {
  supported <- c("PredInstance", "PredInstance_SAT", "PredInstance_FM")
  if (is.null(requested)) {
    dimensions <- supported[supported %in% available]
  } else {
    dimensions <- unique(trimws(unlist(strsplit(as.character(requested), ",", fixed = TRUE))))
    if (!length(dimensions) || any(!nzchar(dimensions))) stop("empty instance-dimension request")
    if (any(!dimensions %in% supported)) stop("unsupported instance dimension requested")
    missing <- setdiff(dimensions, available)
    if (length(missing)) {
      stop(paste("missing requested instance dimension(s):", paste(missing, collapse = ", ")))
    }
  }
  if (!length(dimensions)) stop("no supported instance dimensions found in the point cloud")
  dimensions
}
