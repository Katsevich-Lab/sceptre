#!/usr/bin/env Rscript
# Recreate mathematical inputs only. Never run an old package or save its output.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply a fresh output directory outside this checkout.")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
repo <- normalizePath(file.path(dirname(script), "../.."), mustWork = TRUE)
output <- path.expand(args[[1L]])
if (dir.exists(output) && length(list.files(output, all.files = TRUE, no.. = TRUE)))
    stop("Use a fresh empty output directory.")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output, mustWork = TRUE)
if (startsWith(paste0(output, "/"), paste0(repo, "/"))) stop("Write outside the checkout.")
env <- new.env(parent = globalenv())
sys.source(file.path(repo, "tests/testthat/helper-refactor.R"), env)
workflow <- env$refactor_make_inputs()
saveRDS(workflow, file.path(output, "inputs.rds"), version = 3)
saveRDS(env$refactor_case_definitions(), file.path(output, "cases.rds"), version = 3)
saveRDS(env$refactor_expand_inputs(workflow), file.path(output, "backed-inputs.rds"), version = 3)
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
saveRDS(inputs, file.path(output, "native-inputs.rds"), version = 3)
cat("Wrote input-only fixtures to ", output, "\n", sep = "")
