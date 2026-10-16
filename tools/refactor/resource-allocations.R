#!/usr/bin/env Rscript
# Separate diagnostic run: profiling/tracing never contaminates elapsed trials.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop("Usage: resource-allocations.R <library> <case.rds> <fixtures> <cold|warm> <output.json> <original|candidate>")
source_role <- match.arg(args[[6L]], c("original", "candidate"))
.libPaths(c(normalizePath(args[[1L]]), .libPaths()))
library(sceptre)
stopifnot(identical(normalizePath(find.package("sceptre")), file.path(normalizePath(args[[1L]]), "sceptre")))
source(file.path(dirname(dirname(normalizePath(args[[3L]]))), "helper-refactor.R"))
set.seed(271828)
case <- readRDS(args[[2L]])
parameters <- case$parameters
if (!capabilities("profmem")) stop("This R build does not support allocation profiling.")
data <- case$data
backing_directory <- if (parameters$backed) tempfile("sceptre-allocation-backing-") else NULL
if (!is.null(backing_directory)) dir.create(backing_directory)
object <- refactor_capture({
    object <- import_data(data$response_matrix, data$grna_matrix,
        data$grna_target_data_frame, moi = parameters$moi,
        extra_covariates = data$extra_covariates,
        use_ondisc = parameters$backed, directory_to_write = backing_directory)
    object <- set_analysis_parameters(object, discovery_pairs = parameters$discovery_pairs,
        formula_object = parameters$formula, control_group = parameters$control,
        resampling_mechanism = parameters$mechanism, resampling_approximation = "skew_normal")
    object <- assign_grnas(object, method = "thresholding", threshold = parameters$assignment_threshold, print_progress = FALSE)
    run_qc(object, response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0, 1),
        n_nonzero_trt_thresh = 7L, n_nonzero_cntrl_thresh = 7L)
})$value
run <- function(x) refactor_capture(run_discovery_analysis(x, print_progress = FALSE, output_amount = 2L))$value
mode <- match.arg(args[[4L]], c("cold", "warm"))
if (mode == "warm") object <- run(object)
calls <- new.env(parent = emptyenv())
trace_one <- function(name, namespace) {
    calls[[name]] <- 0L
    tracer <- substitute({
        .sceptre_trace_counter <- get("calls", envir = globalenv())
        .sceptre_trace_counter[[KEY]] <- .sceptre_trace_counter[[KEY]] + 1L
    }, list(KEY = name))
    trace(name, tracer = tracer, where = asNamespace(namespace), print = FALSE)
}
trace_one("glm.fit", "stats")
# Enumerate the three actual production samplers explicitly for each source.
# The original-only slow CRT implementation is archival, not an active bank
# producer. A missing mapped function is a hard trace error, never filtered out.
sampler_functions <- switch(source_role,
    original = c(fisher_yates = "fisher_yates_sampler",
        hybrid = "hybrid_fisher_iwor_sampler", crt_fast_production = "crt_index_sampler_fast"),
    candidate = c(fisher_yates = "fisher_yates_sampler",
        hybrid = "hybrid_fisher_iwor_sampler", crt_fast_production = "sample_crt_indices"))
samplers <- unname(sampler_functions)
for (sampler in samplers) trace_one(sampler, "sceptre")
profile <- tempfile("sceptre-rprofmem-")
gc()
utils::Rprofmem(profile)
result <- run(object)
utils::Rprofmem(NULL)
untrace("glm.fit", where = asNamespace("stats"))
for (sampler in samplers) untrace(sampler, where = asNamespace("sceptre"))
lines <- readLines(profile)
allocations <- as.double(sub(" .*", "", lines[grepl("^[0-9]+ ", lines)]))
jsonlite::write_json(list(
    case = sub("\\.rds$", "", basename(args[[2L]])), cache = mode,
    source_role = source_role, sampler_functions = as.list(sampler_functions),
    passing_pairs = result@n_ok_discovery_pairs,
    current_response_fit_expectation = if (mode == "cold") length(unique(object@discovery_pairs_with_info$response_id[object@discovery_pairs_with_info$pass_qc])) else 0L,
    assignment_fit_and_draw_bank_expectation = if (parameters$mechanism == "crt") length(unique(object@discovery_pairs_with_info$grna_group[object@discovery_pairs_with_info$pass_qc])) else 0L,
    permutation_draw_bank_expectation = if (parameters$mechanism == "permutations") 1L else 0L,
    scope = "Separate serial diagnostic; R allocations exclude native C++/HDF5 malloc and fork memory; tracing not present in timing trials",
    input_sha256 = digest::digest(file = args[[2L]], algo = "sha256"),
    R_allocation_count = length(allocations), R_allocation_bytes = sum(allocations),
    largest_R_allocation_bytes = max(allocations), calls = as.list(calls),
    scientific_sha256 = digest::digest(refactor_plain(result@discovery_result), algo = "sha256"),
    B = c(object@B1, object@B2, object@B3), package_path = find.package("sceptre")
), args[[5L]], pretty = TRUE, auto_unbox = TRUE)
unlink(profile)
if (!is.null(backing_directory)) unlink(backing_directory, recursive = TRUE)
