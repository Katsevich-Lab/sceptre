#!/usr/bin/env Rscript
# Explicit regeneration only:
# Rscript tools/refactor/generate-native.R <reference-source-dir> <output-dir>
# Tests never invoke this script. Production files must match the pinned commit.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply reference source directory and fixture output directory.")
reference <- normalizePath(args[[1L]], mustWork = TRUE)
output <- normalizePath(args[[2L]], mustWork = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
repository <- normalizePath(file.path(dirname(script), "../.."), mustWork = TRUE)
pinned <- "21f9ea098b69c444a884a49bc254de711968e3b5"
production <- c("DESCRIPTION", "NAMESPACE", "R", "src")
git <- function(arguments, directory = repository) {
    result <- system2("git", c("-C", shQuote(directory), arguments), stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(result, "status"))) stop(paste(result, collapse = "\n"))
    result
}
invisible(git(c("cat-file", "-e", paste0(pinned, "^{commit}"))))
tree <- git(c("ls-tree", "-r", pinned, "--", production))
paths <- sub("^[^\t]+\t", "", tree)
hashes <- sub("^[0-9]+ blob ([0-9a-f]+)\t.*$", "\\1", tree)
if (!all(file.exists(file.path(reference, paths)))) stop("Reference production files are missing.")
actual_hashes <- git(c("hash-object", "--no-filters", "--", shQuote(file.path(reference, paths))))
if (!identical(actual_hashes, hashes)) stop("Reference production files differ from the pinned source.")
# Also reject added source files, including ignored source files. Compiled build
# products in src/ are deliberately outside the production source inventory.
actual_r <- paste0("R/", list.files(file.path(reference, "R"), recursive = TRUE, all.files = TRUE, no.. = TRUE))
actual_src <- paste0("src/", list.files(file.path(reference, "src"), recursive = TRUE, all.files = TRUE, no.. = TRUE))
actual_src <- actual_src[grepl("\\.(c|cc|cpp|cxx|h|hpp|hxx|f|f90|f95|m|mm|inc)$|/Make(vars|file)", actual_src)]
added <- setdiff(c(actual_r, actual_src), paths)
if (length(added)) stop("Additional production source invalidates the reference: ", paste(added, collapse = ", "))
if (!identical(unname(read.dcf(file.path(reference, "DESCRIPTION"))[1L, "Version"]), "0.99.0")) {
    stop("Unexpected reference package version.")
}
source(file.path(repository, "tests/testthat/helper-refactor-native.R"))
pkgload::load_all(reference, export_all = FALSE, helpers = FALSE, quiet = TRUE)
bridge <- refactor_native_bridge(file.path(repository, "tests/testthat/fixtures/refactor/native-bridge.cpp"))

# One treated cell, unit weights, and zero nuisance projection make a selected
# a-value its exact statistic. y/mu retain nontrivial effect/SE calculations.
scores <- c(8, stats::qnorm((seq_len(1000L) - 0.5) / 1000), -8, 0)
normal_bank <- function(n) as.list(as.integer(1L + ((seq_len(n) - 1L) %% 1000L)))
balanced_bank <- function(n) as.list(as.integer(1L + (2L * (seq_len(n) - 1L) %% 1000L)))
case <- function(trt = 1L, side = 1L, fit = FALSE, B2 = 0L, B3 = 7L,
                 first = balanced_bank(499L), second = normal_bank(B2),
                 third = as.list(c(17L, 201L, 498L, 500L, 701L, 901L, 999L))[seq_len(B3)],
                 expected_stage, expected_sn = FALSE) {
    list(
        y = as.numeric(2L + seq_along(scores) %% 3L),
        mu = 1.2 + (seq_along(scores) %% 5L) / 10,
        a = scores, w = rep(1, length(scores)), D = matrix(0, 1L, length(scores)),
        trt_idxs = as.integer(trt), n_trt = 1L, use_all_cells = FALSE,
        B1 = 499L, B2 = as.integer(B2), B3 = as.integer(B3),
        fit_parametric_curve = fit, return_resampling_dist = TRUE,
        side_code = as.integer(side), draws = c(first, second, third),
        expected_stage = as.integer(expected_stage), expected_sn = expected_sn
    )
}
threshold_bank <- function(exceedances, trt = 1L) {
    c(rep(list(as.integer(trt - 1L)), exceedances), rep(list(1002L), 499L - exceedances))
}
cases <- list(
    stage1_right = case(trt = 502L, expected_stage = 1L),
    stage1_left = case(trt = 502L, side = -1L, expected_stage = 1L),
    stage1_both = case(trt = 502L, side = 0L, expected_stage = 1L),
    threshold_right = case(first = threshold_bank(9L), expected_stage = 3L),
    above_threshold_right = case(first = threshold_bank(10L), expected_stage = 1L),
    threshold_left = case(trt = 1002L, side = -1L,
                          first = threshold_bank(9L, 1002L), expected_stage = 3L),
    threshold_both = case(side = 0L, first = threshold_bank(4L), expected_stage = 3L),
    above_threshold_both = case(side = 0L, first = threshold_bank(5L), expected_stage = 1L),
    accepted_sn_right = case(fit = TRUE, B2 = 4999L, B3 = 0L,
                             expected_stage = 2L, expected_sn = TRUE),
    accepted_sn_left = case(trt = 1002L, side = -1L, fit = TRUE, B2 = 4999L, B3 = 0L,
                            expected_stage = 2L, expected_sn = TRUE),
    accepted_sn_both = case(side = 0L, fit = TRUE, B2 = 4999L, B3 = 0L,
                            expected_stage = 2L, expected_sn = TRUE),
    failed_sn_stage3 = case(fit = TRUE, B2 = 4999L,
                            second = rep(list(1002L), 4999L), expected_stage = 3L),
    failed_sn_B3_zero = case(fit = TRUE, B2 = 4999L, B3 = 0L,
                             second = rep(list(1002L), 4999L), expected_stage = 3L),
    failed_sn_outlier_B3_zero = case(fit = TRUE, B2 = 4999L, B3 = 0L,
                                     second = c(normal_bank(4998L), list(0L)), expected_stage = 3L),
    empirical_B3_zero = case(B3 = 0L, expected_stage = 3L)
)
# Include a nonzero nuisance projection and two observed cells so both the
# studentizer and the effect standard error have nontrivial reference values.
projected <- case(trt = c(502L, 700L), expected_stage = 1L)
projected$n_trt <- 2L
projected$D <- rbind(0.01 + (seq_along(scores) %% 7L) / 10000,
                    0.015 * sin(seq_along(scores)))
projected$draws <- lapply(projected$draws, function(draw) c(draw, 1002L))
cases$projected_multi_cell <- projected
variable <- case(trt = 502L, expected_stage = 1L)
variable$use_all_cells <- TRUE
variable$draws <- lapply(seq_along(variable$draws), function(i) {
    if (i %% 2L) variable$draws[[i]] else c(variable$draws[[i]], 1002L)
})
cases$variable_crt_cardinality <- variable
samplers <- list(
    fisher_yates = list(function_name = "fisher_yates_sampler", args = list(n_tot = 13L, M = 4L, B = 17L)),
    hybrid = list(function_name = "hybrid_fisher_iwor_sampler", args = list(N = 9L, m = 2L, M = 5L, B = 17L)),
    crt_fast_production = list(function_name = "crt_index_sampler_fast", args = list(
        fitted_probabilities = c(0, 0.15, 0.4, 0.85, 1), B = 499L + 4999L + 0L
    )),
    crt_legacy = list(function_name = "crt_index_sampler", args = list(
        fitted_probabilities = c(0, 0.15, 0.4, 0.85, 1), B = 17L
    ))
)
inputs <- list(index_base = list(observed_treatment = 1L, synthetic_draws = 0L),
               numeric_tolerance = 1e-8, cases = cases, samplers = samplers)
outputs <- lapply(cases, refactor_native_run_case, bridge = bridge)
for (name in names(cases)) {
    stopifnot(identical(outputs[[name]]$stage, cases[[name]]$expected_stage))
    stopifnot(identical(all(is.finite(outputs[[name]]$sn_params)), cases[[name]]$expected_sn))
}
draws <- lapply(samplers, refactor_native_sample, bridge = bridge)
expected <- list(
    provenance = list(production_commit = pinned, reference_checkout = tryCatch(
                          git(c("rev-parse", "HEAD"), reference),
                          error = function(e) "archive verified against pinned Git blobs"),
                      platform = refactor_native_platform(), tolerance = 1e-8),
    cases = outputs, sampler_draws = draws
)
saveRDS(inputs, file.path(output, "native-inputs.rds"), version = 3)
saveRDS(expected, file.path(output, "native-expected.rds"), version = 3)
diagnostic <- do.call(rbind, lapply(names(draws), function(name) {
    bank <- draws[[name]]
    do.call(rbind, lapply(seq_len(min(8L, length(bank))), function(draw) {
        data.frame(sampler = name, requested_B = samplers[[name]]$args$B,
                   draw = draw - 1L, position = seq_along(bank[[draw]]) - 1L,
                   native_cell = bank[[draw]])
    }))
}))
utils::write.table(diagnostic, file.path(output, "native-draws.tsv"), sep = "\t", quote = FALSE,
                   row.names = FALSE, na = "NA")
cat("Generated", length(cases), "native cases and", length(draws), "sampler banks from", pinned, "\n")
