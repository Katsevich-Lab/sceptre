#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop("Usage: resource-worker.R <library> <case.rds> <cold|warm> <workers> <output.json> <phase.json> <fixtures>")
library_path <- normalizePath(args[[1L]])
.libPaths(c(library_path, .libPaths()))
library(sceptre)
stopifnot(identical(normalizePath(find.package("sceptre")), file.path(library_path, "sceptre")))
source(file.path(dirname(dirname(normalizePath(args[[7L]]))), "helper-refactor.R"))
mode <- match.arg(args[[3L]], c("cold", "warm"))
workers <- as.integer(args[[4L]])
output <- args[[5L]]
phase_path <- args[[6L]]
phase <- function(value) {
    tmp <- paste0(phase_path, ".tmp")
    jsonlite::write_json(list(pid = Sys.getpid(), phase = value), tmp, auto_unbox = TRUE)
    file.rename(tmp, phase_path)
    invisible(NULL)
}
phase("preparation")
set.seed(271828)
case <- readRDS(args[[2L]])
data <- case$data
parameters <- case$parameters
recipe <- case$recipe
input_sha256 <- digest::digest(file = args[[2L]], algo = "sha256")
rm(case)
directory <- tempfile("sceptre-resource-")
dir.create(directory)
preparation <- refactor_capture({
    object <- sceptre::import_data(
        response_matrix = data$response_matrix, grna_matrix = data$grna_matrix,
        grna_target_data_frame = data$grna_target_data_frame,
        extra_covariates = data$extra_covariates, moi = parameters$moi,
        use_ondisc = parameters$backed,
        directory_to_write = if (parameters$backed) directory else NULL
    )
    object <- sceptre::set_analysis_parameters(object,
        discovery_pairs = parameters$discovery_pairs, formula_object = parameters$formula,
        control_group = parameters$control, resampling_mechanism = parameters$mechanism,
        resampling_approximation = "skew_normal", side = "both"
    )
    object <- sceptre::assign_grnas(object, method = "thresholding",
        threshold = parameters$assignment_threshold, print_progress = FALSE)
    sceptre::run_qc(object, response_n_umis_range = c(0, 1),
        response_n_nonzero_range = c(0, 1), n_nonzero_trt_thresh = 7L,
        n_nonzero_cntrl_thresh = 7L)
})
object <- preparation$value
preparation_warnings <- preparation$warnings
rm(preparation, data)
run <- function(x) refactor_capture(sceptre::run_discovery_analysis(
    x, print_progress = FALSE, parallel = workers > 1L,
    n_processors = max(2L, workers), output_amount = 2L, log_dir = directory
))
warm_warnings <- character()
if (mode == "warm") {
    warm <- run(object)
    object <- warm$value
    warm_warnings <- warm$warnings
    rm(warm)
}
before_fits <- names(object@response_precomputations)
gc()
phase("analysis")
elapsed <- system.time({ result <- run(object) })
phase("finalization")
plain <- refactor_plain(result$value@discovery_result)
report <- list(
    case = sub("\\.rds$", "", basename(args[[2L]])), cache = mode, workers = workers,
    elapsed_seconds = unname(elapsed[["elapsed"]]),
    user_self_seconds = unname(elapsed[["user.self"]]),
    system_self_seconds = unname(elapsed[["sys.self"]]),
    user_children_seconds = unname(elapsed[["user.child"]]),
    system_children_seconds = unname(elapsed[["sys.child"]]),
    cells_in_use = length(object@cells_in_use), requested_pairs = nrow(parameters$discovery_pairs),
    passing_pairs = result$value@n_ok_discovery_pairs,
    result_rows = nrow(result$value@discovery_result),
    B = c(B1 = object@B1, B2 = object@B2, B3 = object@B3),
    fits_before = before_fits, fits_after = names(result$value@response_precomputations),
    preparation_warnings = preparation_warnings, warm_warnings = warm_warnings,
    analysis_warnings = result$warnings,
    scientific_sha256 = digest::digest(plain, algo = "sha256"),
    input_sha256 = input_sha256, recipe = recipe,
    package_path = find.package("sceptre"), package_version = as.character(packageVersion("sceptre")),
    ondisc_version = if (requireNamespace("ondisc", quietly = TRUE)) as.character(packageVersion("ondisc")) else NULL,
    R = R.version.string, external = extSoftVersion(), platform = R.version$platform,
    threads = Sys.getenv(c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"))
)
jsonlite::write_json(report, output, pretty = TRUE, auto_unbox = TRUE, null = "null")
saveRDS(plain, sub("\\.json$", "-scientific.rds", output))
unlink(directory, recursive = TRUE)
phase("complete")
