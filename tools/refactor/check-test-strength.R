#!/usr/bin/env Rscript
# Temporary adversarial check of the new permanent tests, not an old/new test.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: check-test-strength.R /installed-library /fresh/output")
.libPaths(c(normalizePath(args[[1L]], mustWork = TRUE), .libPaths()))
suppressPackageStartupMessages(library(sceptre))
suppressPackageStartupMessages(library(testthat))
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
repo <- normalizePath(file.path(dirname(script), "../.."), mustWork = TRUE)
out <- args[[2L]]
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Use a fresh empty output directory.")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source(file.path(repo, "tests/testthat/helper-refactor.R"))
source(file.path(repo, "tests/testthat/helper-refactor-native.R"))
options(sceptre.refactor.fixture_dir = file.path(repo, "tests/testthat/fixtures/refactor"))
records <- list()
reject <- function(name, code, pattern) {
    condition <- tryCatch({ force(code); NULL }, expectation_failure = function(e) e)
    if (!inherits(condition, "expectation_failure") || !grepl(pattern, conditionMessage(condition)))
        stop("Negative control was not rejected by its intended assertion: ", name, " / ", paste(class(condition), collapse = ","), " / ", if (is.null(condition)) "NO FAILURE" else conditionMessage(condition))
    records[[name]] <<- list(rejected = TRUE, assertion = conditionMessage(condition))
}
inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
case <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))$low_complement_crt
actual <- refactor_run_case(inputs, case)
expect_workflow_results(actual, inputs, case) # The genuine output must pass first.
for (kind in c("p_value", "pair_label", "warning", "empty_table", "qc_flag")) {
    bad <- actual
    if (kind == "p_value") bad$value$final$discovery_result$columns$p_value[1L] <- 2
    if (kind == "pair_label") bad$value$final$discovery_result$columns$response_id[1L] <- "r8"
    if (kind == "warning") bad$warnings <- "Unexpected workflow warning"
    if (kind == "empty_table") bad$value$final$discovery_result$columns <- lapply(bad$value$final$discovery_result$columns, function(x) x[FALSE])
    if (kind == "qc_flag") {
        z <- bad$value$final$discovery_result$columns
        i <- which(z$response_id == "r4")
        z$pass_qc[i] <- TRUE; z$p_value[i] <- .1
        bad$value$final$discovery_result$columns <- z
    }
    pattern <- switch(kind, p_value = "eligible p-values", warning = "unexpected workflow warnings", "pair identities and independent QC")
    reject(kind, expect_workflow_results(bad, inputs, case), pattern)
}
bridge <- refactor_native_bridge(file.path(repo, "tests/testthat/fixtures/refactor/native-bridge.cpp"))
p <- c(0, .15, .4, .85, 1); B <- 10000L
joint <- function(pairs) p[pairs[1L, ]] * p[pairs[2L, ]]
real <- refactor_native_sample(list(args = list(fitted_probabilities = p, B = B)), "sample_crt_indices", bridge)
expect_sampler_probabilities(real, p, joint)
reject("constant_CRT", expect_sampler_probabilities(rep(list(4L), B), p, joint), "cell inclusion")
u <- (seq_len(B) - .5) / B
correlated <- lapply(u, function(x) as.integer(which(p > x) - 1L))
reject("correlated_CRT", expect_sampler_probabilities(correlated, p, joint), "joint inclusion")
real <- refactor_native_sample(list(args = list(n_tot = 5L, M = 2L, B = B)), "fisher_yates_sampler", bridge)
wor_joint <- function(pairs) rep(.1, ncol(pairs))
expect_sampler_probabilities(real, rep(.4, 5L), wor_joint)
reject("constant_permutation", expect_sampler_probabilities(rep(list(0:1), B), rep(.4, 5L), wor_joint), "cell inclusion")
jsonlite::write_json(list(positive_controls_passed = TRUE, rejected = records),
    file.path(out, "negative-controls.json"), pretty = TRUE, auto_unbox = TRUE)
cat("Genuine outputs pass; all ", length(records), " deliberately wrong outputs are rejected.\n", sep = "")
