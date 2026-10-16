#!/usr/bin/env Rscript
# Separate diagnostic run: profiling/tracing never contaminates elapsed trials.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("Usage: resource-allocations.R <library> <compact-case.rds> <fixtures> <cold|warm> <output.json>")
.libPaths(c(normalizePath(args[[1L]]), .libPaths()))
library(sceptre)
stopifnot(identical(normalizePath(find.package("sceptre")), file.path(normalizePath(args[[1L]]), "sceptre")))
source(file.path(dirname(dirname(normalizePath(args[[3L]]))), "helper-refactor.R"))
set.seed(271828)
case <- readRDS(args[[2L]])
parameters <- case$parameters
data <- case$data
object <- refactor_capture({
    object <- import_data(data$response_matrix, data$grna_matrix,
        data$grna_target_data_frame, moi = parameters$moi,
        extra_covariates = data$extra_covariates)
    object <- set_analysis_parameters(object, discovery_pairs = parameters$discovery_pairs,
        formula_object = parameters$formula, control_group = parameters$control,
        resampling_mechanism = parameters$mechanism, resampling_approximation = "skew_normal")
    object <- assign_grnas(object, method = "thresholding", threshold = 40, print_progress = FALSE)
    run_qc(object, response_n_umis_range = c(0, 1), response_n_nonzero_range = c(0, 1),
        n_nonzero_trt_thresh = 7L, n_nonzero_cntrl_thresh = 7L)
})$value
run <- function(x) refactor_capture(run_discovery_analysis(x, print_progress = FALSE, output_amount = 2L))$value
mode <- match.arg(args[[4L]], c("cold", "warm"))
if (mode == "warm") object <- run(object)
calls <- new.env(parent = emptyenv())
trace_one <- function(name, namespace) {
    calls[[name]] <- 0L
    count <- local({ key <- name; function() calls[[key]] <- calls[[key]] + 1L })
    trace(name, tracer = count, where = asNamespace(namespace), print = FALSE)
}
trace_one("glm.fit", "stats")
samplers <- c("fisher_yates_sampler", "hybrid_fisher_iwor_sampler", "crt_index_sampler", "crt_index_sampler_fast")
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
    cache = mode, scope = "Separate serial diagnostic; R allocations exclude native C++/HDF5 malloc and fork memory; tracing not present in timing trials",
    input_sha256 = digest::digest(file = args[[2L]], algo = "sha256"),
    R_allocation_count = length(allocations), R_allocation_bytes = sum(allocations),
    largest_R_allocation_bytes = max(allocations), calls = as.list(calls),
    scientific_sha256 = digest::digest(refactor_plain(result@discovery_result), algo = "sha256"),
    B = c(object@B1, object@B2, object@B3), package_path = find.package("sceptre")
), args[[5L]], pretty = TRUE, auto_unbox = TRUE)
unlink(profile)
