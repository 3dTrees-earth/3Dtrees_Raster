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
  memory_limit <- min(declared("GALAXY_MEMORY_MB") / 1024,
                      read_limit("memory.max") / 1024^3,
                      read_limit("memory/memory.limit_in_bytes") / 1024^3)
  requested <- if (is.null(memory_budget)) 60 else
    positive_resource(memory_budget, "--memory-budget-gib")
  if (requested > 70) stop("--memory-budget-gib must be at most 70")
  list(threads = as.integer(threads),
       memory_budget_gib = min(requested, memory_limit * 0.85),
       memory_limit_gib = if (is.finite(memory_limit)) memory_limit else NULL)
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
