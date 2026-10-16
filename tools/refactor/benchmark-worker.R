#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("Usage: benchmark-worker.R <library> <fixtures> <case> <cold|warm> <output.json>")
library_path <- normalizePath(args[[1L]])
fixtures <- normalizePath(args[[2L]])
.libPaths(c(library_path, .libPaths()))
library(sceptre)
stopifnot(identical(normalizePath(find.package("sceptre")), file.path(library_path, "sceptre")))
source(file.path(dirname(dirname(fixtures)), "helper-refactor.R"))
inputs <- readRDS(file.path(fixtures, "backed-inputs.rds"))
cases <- readRDS(file.path(fixtures, "cases.rds"))
name <- args[[3L]]
case <- cases[[if (name == "backed_perm") "high_perm" else name]]
if (is.null(case)) stop("Unknown case")
mode <- match.arg(args[[4L]], c("cold", "warm"))
directory <- tempfile("refactor-bench-")
object <- refactor_capture(refactor_prepare(inputs, case, backed = name == "backed_perm", directory = directory))$value
if (mode == "warm") object <- refactor_capture(sceptre::run_discovery_analysis(object, print_progress = FALSE, output_amount = 2L))$value
gc()
before_fits <- names(object@response_precomputations)
elapsed <- system.time({
    result <- refactor_capture(sceptre::run_discovery_analysis(object, print_progress = FALSE, output_amount = 2L))
})
report <- list(
    case = name, cache = mode, elapsed_seconds = unname(elapsed[["elapsed"]]),
    user_seconds = unname(elapsed[["user.self"]]), system_seconds = unname(elapsed[["sys.self"]]),
    cells = ncol(inputs$response_matrix), responses = nrow(inputs$response_matrix),
    passing_pairs = result$value@n_ok_discovery_pairs, fits_before = before_fits,
    fits_after = names(result$value@response_precomputations), warnings = result$warnings,
    package_path = find.package("sceptre"), package_version = as.character(packageVersion("sceptre")),
    ondisc_version = if (requireNamespace("ondisc", quietly = TRUE)) as.character(packageVersion("ondisc")) else NULL,
    R = R.version.string, external = extSoftVersion(), threads = Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")),
    parallel = FALSE, memory_scope = "fresh serial R process; excludes process-tree/fork-worker validation"
)
jsonlite::write_json(report, args[[5L]], pretty = TRUE, auto_unbox = TRUE, null = "null")
unlink(directory, recursive = TRUE)
