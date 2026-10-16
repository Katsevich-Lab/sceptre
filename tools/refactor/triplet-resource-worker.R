#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8L) stop("Supply library fixture orientation mode result phase clock repetitions.")
library <- normalizePath(args[[1L]], mustWork = TRUE)
.libPaths(c(library, .libPaths()))
fixture <- args[[2L]]
row_access <- identical(args[[3L]], "row")
mode <- match.arg(args[[4L]], c("correctness", "performance", "profile"))
result_path <- args[[5L]]
phase_path <- args[[6L]]
clock <- dyn.load(args[[7L]], local = TRUE)
repetitions <- as.integer(args[[8L]])
stopifnot(repetitions >= 1L)
threads <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
stopifnot(all(Sys.getenv(threads) == "1"))
data.table::setDTthreads(1L)
stopifnot(data.table::getDTthreads() == 1L)
namespace <- loadNamespace("sceptre")
stopifnot(identical(normalizePath(getNamespaceInfo(namespace, "path")),
                    normalizePath(file.path(library, "sceptre"))))
convert <- getFromNamespace("set_matrix_accessibility", "sceptre")
monotonic <- function() .Call("sceptre_triplet_monotonic", PACKAGE = clock[["name"]])
triplet_resource_mark <- function(phase) {
    value <- list(phase = phase, pid = Sys.getpid(), unix_seconds = as.numeric(Sys.time()),
                  monotonic_seconds = monotonic())
    line <- jsonlite::toJSON(value, auto_unbox = TRUE, digits = 16)
    cat(line, "\n", file = paste0(result_path, ".events.jsonl"), append = TRUE, sep = "")
    temporary <- paste0(phase_path, ".new")
    writeLines(line, temporary)
    stopifnot(file.rename(temporary, phase_path))
}
stream_hash <- function(object) {
    path <- tempfile(tmpdir = dirname(result_path))
    connection <- file(path, "wb")
    serialize(object, connection, version = 3)
    close(connection)
    on.exit(unlink(path))
    digest::digest(file = path, algo = "sha256")
}
if (mode == "profile") {
    for (function_name in c("data.table", "copy", "setDT", "setorderv")) {
        before <- substitute(get("triplet_resource_mark", .GlobalEnv)(PHASE),
            list(PHASE = paste0(function_name, "_start")))
        after <- substitute(get("triplet_resource_mark", .GlobalEnv)(PHASE),
            list(PHASE = paste0(function_name, "_end")))
        trace(function_name, where = asNamespace("data.table"), tracer = before,
              exit = after, print = FALSE)
    }
}
intervals <- list()
gc_record <- function(value) {
    if (is.null(value)) return(NULL)
    list(rows = rownames(value), columns = colnames(value), values = unname(value))
}
for (iteration in seq_len(repetitions)) {
    triplet_resource_mark("load_start")
    input <- readRDS(fixture)
    before_hash <- if (mode == "correctness") stream_hash(input) else NULL
    triplet_resource_mark("input_ready_pre_gc")
    pre_gc <- gc(reset = TRUE)
    # Untimed marker padding lets the external 20 ms sampler observe the
    # post-GC input-only and post-conversion live-input/output states directly.
    triplet_resource_mark("pre_conversion_ready")
    Sys.sleep(0.06)
    if (mode == "profile") Rprofmem(paste0(result_path, ".rprofmem"))
    triplet_resource_mark("conversion_start")
    start_cpu <- proc.time()
    start <- monotonic()
    output <- convert(input, row_access)
    elapsed <- monotonic() - start
    cpu <- proc.time() - start_cpu
    triplet_resource_mark("conversion_end")
    if (mode == "profile") Rprofmem(NULL)
    Sys.sleep(0.06)
    stopifnot(identical(dim(output), dim(input)), length(output@x) == length(input@x))
    # Input and output remain live through the conversion-end/RSS marker.
    # Expensive hashes/serialization happen only in untimed correctness children.
    correctness <- if (mode == "correctness") {
        after_hash <- stream_hash(input)
        saveRDS(output, paste0(result_path, ".output.rds"), compress = FALSE, version = 3)
        list(input_before_sha256 = before_hash, input_after_sha256 = after_hash,
             caller_unchanged = identical(before_hash, after_hash),
             output_sha256 = digest::digest(file = paste0(result_path, ".output.rds"), algo = "sha256"),
             class = class(output)[[1L]], stored = length(output@x),
             zeros = sum(output@x == 0), fractions = sum(output@x != floor(output@x)))
    } else NULL
    post_gc <- if (mode == "profile") gc() else NULL
    intervals[[iteration]] <- list(elapsed_seconds = elapsed, user_seconds = unname(cpu[["user.self"]]),
        system_seconds = unname(cpu[["sys.self"]]), pre_gc = gc_record(pre_gc), post_gc = gc_record(post_gc),
        correctness = correctness)
    rm(input, output)
    gc()
}
triplet_resource_mark("exit")
result <- list(mode = mode, orientation = args[[3L]], fixture = normalizePath(fixture),
    fixture_sha256 = digest::digest(file = fixture, algo = "sha256"), package = find.package("sceptre"),
    repetitions = repetitions, elapsed_seconds = sum(vapply(intervals, `[[`, numeric(1), "elapsed_seconds")),
    intervals = intervals, threads = Sys.getenv(threads), data_table_threads = data.table::getDTthreads(),
    R = R.version.string, platform = R.version$platform,
    dependencies = lapply(c("Matrix", "data.table", "Rcpp"), function(name) list(name = name,
        version = as.character(utils::packageVersion(name)), path = find.package(name))),
    dependency_dll_sha256 = lapply(getLoadedDLLs(), function(dll) {
        path <- dll[["path"]]
        if (file.exists(path)) list(path = path, sha256 = digest::digest(file = path, algo = "sha256")) else NULL
    }), BLAS = extSoftVersion()[["BLAS"]])
jsonlite::write_json(result, result_path, pretty = TRUE, auto_unbox = TRUE, digits = 16)
