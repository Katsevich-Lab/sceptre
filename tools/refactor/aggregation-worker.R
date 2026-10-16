#!/usr/bin/env Rscript
# One installed aggregation child. Profiling/correctness are separate from timing.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop("Usage: aggregation-worker.R library input.rds result.json phase.json mode iterations")
.libPaths(c(normalizePath(args[1L]), .libPaths()))
suppressPackageStartupMessages(library(sceptre))
data.table::setDTthreads(1L)
stopifnot(identical(normalizePath(find.package("sceptre")),
                    file.path(normalizePath(args[1L]), "sceptre")))
mode <- match.arg(args[5L], c("timing", "correctness", "profile", "batch"))
iterations <- as.integer(args[6L])
stopifnot(iterations > 0L, mode == "batch" || iterations == 1L)
phase_path <- args[4L]
events <- list()
mark <- function(phase, iteration = 0L) {
    state <- list(pid = Sys.getpid(), phase = phase, iteration = iteration,
                  unix_seconds = as.numeric(Sys.time()))
    events[[length(events) + 1L]] <<- state
    tmp <- paste0(phase_path, ".tmp")
    jsonlite::write_json(state, tmp, auto_unbox = TRUE)
    if (!file.rename(tmp, phase_path)) stop("Cannot publish aggregation phase marker.")
}
await_sample <- function(phase, iteration) {
    marker <- paste0(phase_path, ".", phase, "-", iteration, ".ack")
    deadline <- Sys.time() + 10
    while (!file.exists(marker)) {
        if (Sys.time() > deadline) stop("RSS monitor did not acknowledge phase: ", phase)
        Sys.sleep(.005)
    }
}
mark("input_start")
input <- readRDS(args[2L])
input_hash <- digest::digest(input, algo = "sha256")
mark("input_ready")
fn <- getFromNamespace("collapse_grna_matrix", "sceptre")
warnings <- character()
elapsed <- numeric(iterations)
user <- system <- numeric(iterations)
profile_path <- paste0(args[3L], ".Rprofmem")
profile <- NULL
heap <- NULL
for (iteration in seq_len(iterations)) {
    if (iteration > 1L) input <- readRDS(args[2L])
    gc(reset = TRUE) |> invisible()
    mark("pre_call", iteration)
    await_sample("pre_call", iteration)
    mark("aggregation_start", iteration)
    if (mode == "profile") utils::Rprofmem(profile_path, threshold = 65536)
    start <- proc.time()
    output <- withCallingHandlers(fn(input$m, input$d), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    })
    duration <- proc.time() - start
    if (mode == "profile") utils::Rprofmem(NULL)
    elapsed[iteration] <- unname(duration["elapsed"])
    user[iteration] <- unname(duration["user.self"])
    system[iteration] <- unname(duration["sys.self"])
    mark("aggregation_end", iteration)
    await_sample("aggregation_end", iteration)
    heap <- gc()
    stopifnot(identical(input_hash, digest::digest(input, algo = "sha256")))
    stopifnot(is.matrix(output), typeof(output) == "double",
              nrow(output) == length(unique(input$d$vector_id)), ncol(output) == ncol(input$m))
    structure <- list(class = class(output), typeof = typeof(output), dim = dim(output),
        dimnames_sha256 = digest::digest(dimnames(output), algo = "sha256"),
        value_sum = sum(output), input_nonmutation = TRUE)
    if (mode == "correctness") saveRDS(output, paste0(args[3L], ".output.rds"), compress = FALSE)
    if (mode == "profile") {
        lines <- readLines(profile_path, warn = FALSE)
        lines <- lines[grepl("^[0-9]+ ", lines)]
        bytes <- as.double(sub(" .*", "", lines))
        index <- order(bytes, decreasing = TRUE)[seq_len(min(12L, length(bytes)))]
        profile <- list(threshold_bytes = 65536L, traced_allocation_count = length(bytes),
            traced_allocation_bytes = sum(bytes), largest_R_allocation_bytes = max(c(0, bytes)),
            dense_input_counterfactual_bytes = input$recipe$dense_input_counterfactual_bytes,
            allocations_at_least_dense_input = sum(bytes >= input$recipe$dense_input_counterfactual_bytes),
            largest_stack_examples = lines[index], trace_file = profile_path,
            heap_high_water = unclass(heap), heap_rows = rownames(heap),
            heap_columns = colnames(heap),
            scope = "R allocations >=64KiB; excludes smaller allocations and may omit native malloc. Profile child is excluded from timing/RSS gates.")
    }
    if (iteration < iterations) { rm(output, input); gc() |> invisible() }
}
mark("output_checked")
dependencies <- c("Matrix", "dplyr", "data.table", "Rcpp", "BH", "ondisc", "digest", "jsonlite")
environment <- list(R = R.version.string, platform = R.version$platform,
    external = as.list(extSoftVersion()), RNG = as.list(RNGkind()),
    threads = as.list(Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"))),
    data_table_threads = data.table::getDTthreads(),
    dependencies = stats::setNames(lapply(dependencies, function(package) {
        if (!requireNamespace(package, quietly = TRUE)) return(NULL)
        list(version = as.character(packageVersion(package)), path = find.package(package))
    }), dependencies))
result <- list(schema_version = 1L, mode = mode, iterations = iterations,
    input_sha256 = digest::digest(file = args[2L], algo = "sha256"), input_object_sha256 = input_hash,
    elapsed_seconds = sum(elapsed), elapsed_intervals = elapsed, user_seconds = sum(user),
    system_seconds = sum(system), output_structure = structure, warnings = warnings,
    environment = environment, package_path = find.package("sceptre"), profile = profile,
    events = events, R_MAX_VSIZE = Sys.getenv("R_MAX_VSIZE"))
jsonlite::write_json(result, args[3L], pretty = TRUE, auto_unbox = TRUE, null = "null", digits = NA)
mark("exit")
