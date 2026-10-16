#!/usr/bin/env Rscript
# Compile actual original/candidate estimator bodies with test-only boundary
# interposition. No package source or installed namespace is modified.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L)
    stop("Supply candidate source and optionally an evidence directory.")
root <- normalizePath(args[[1L]], mustWork = TRUE)
output <- if (length(args) == 2L) args[[2L]] else
    file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "sceptre-dispersion-exceptions")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output, mustWork = TRUE)
original <- file.path(root, "tools/refactor/results/dispersion-characterization-21f9ea0/original.cpp")
candidate <- file.path(root, "src/compute_nb_size_param_functions.cpp")
probe <- file.path(root, "tools/refactor/dispersion-exceptions.cpp")
runner <- file.path(root, "tools/refactor/test-dispersion-exceptions.R")
inputs <- file.path(root, "tools/refactor/results/dispersion-characterization-21f9ea0/report.rds")
sha <- function(path) digest::digest(file = path, algo = "sha256")
files <- c(original, candidate, probe, runner, inputs)
before <- vapply(files, sha, character(1))
stopifnot(sha(original) == "4add3e36820f9af6b794b715862c0981f44cb51f11750cf5107d39c2cd4ee776")

# This bounded correction changes only catches, includes and source contracts.
# Compare all remaining tokens before adding any test-only interposition.
tokens <- function(path) {
    text <- paste(readLines(path), collapse = "\n")
    text <- gsub("//[^\n]*|/\\*[\\s\\S]*?\\*/", "", text, perl = TRUE)
    text <- gsub("(?m)^\\s*#include[^\n]*", "", text, perl = TRUE)
    gsub("[[:space:]]+", "", text)
}
candidate_tokens <- tokens(candidate)
typed <- "catch(conststd::domain_error&){}catch(conststd::overflow_error&){}"
stopifnot(grepl(typed, candidate_tokens, fixed = TRUE),
          identical(sub(typed, "catch(...){}", candidate_tokens, fixed = TRUE), tokens(original)))

instrument <- function(path, destination) {
    source <- readLines(path)
    source <- source[!grepl("^\\s*(#include|// \\[\\[Rcpp::)", source)]
    text <- paste(source, collapse = "\n")
    mle <- "std::vector<double> out = nb_theta_mle(t0, y, mu, limit, eps);"
    mm <- "out = nb_theta_mm(t0, y, mu, dfr, limit, eps);"
    for (call in c(mle, mm))
        stopifnot(length(gregexpr(call, text, fixed = TRUE)[[1L]]) == 1L,
                  grepl(call, text, fixed = TRUE))
    text <- sub(mle, paste0("dispersion_probe::before(1);\n    ", mle), text, fixed = TRUE)
    text <- sub(mm, paste0("dispersion_probe::before(2);\n      ", mm), text, fixed = TRUE)
    writeLines(text, destination)
}
original_instrumented <- file.path(output, "original-instrumented.inc")
candidate_instrumented <- file.path(output, "candidate-instrumented.inc")
instrument(original, original_instrumented)
instrument(candidate, candidate_instrumented)
bridge <- file.path(output, "dispersion-instrumented.cpp")
writeLines(c(
    sprintf('#define SCEPTRE_ORIGINAL_DISPERSION %s', encodeString(original_instrumented, quote = '"')),
    sprintf('#define SCEPTRE_CANDIDATE_DISPERSION %s', encodeString(candidate_instrumented, quote = '"')),
    readLines(probe)
), bridge)
environment <- new.env(parent = globalenv())
Rcpp::sourceCpp(bridge, env = environment, rebuild = TRUE, showOutput = TRUE)
facts <- environment$refactor_dispersion_policy_facts()
stopifnot(facts$double_radix == 2L, facts$double_digits == 53L,
          facts$effective_precision == 53L,
          facts$domain == facts$throw_on_error, facts$pole == facts$throw_on_error,
          facts$overflow == facts$throw_on_error,
          facts$underflow == facts$ignore_error, facts$denorm == facts$ignore_error)

call <- function(candidate, y, mu, dfr, limit = 50L, stage = 0L, type = 0L) {
    environment$refactor_dispersion_case(candidate, y, mu, dfr, limit,
                                         .Machine$double.eps^(1/4), stage, type)
}
numerical <- lapply(readRDS(inputs)$cases, function(case) {
    dfr <- if (is.null(case$input$dfr)) 2 else case$input$dfr
    old <- call(FALSE, case$input$y, case$input$mu, dfr)
    new <- call(TRUE, case$input$y, case$input$mu, dfr)
    stopifnot(identical(old, new), identical(new$escaped, ""))
    list(input = case$input, original = old, candidate = new)
})

# limit=2 follows the actual warning-driven MLE -> MM branch. It is not a
# fabricated fit: the original MLE computation and warning flag both execute.
y <- c(0, 1, 5)
mu <- c(1, 1, 1)
ordinary <- call(TRUE, y, mu, 2)
warning_branch <- call(TRUE, y, mu, 2, limit = 2L)
stopifnot(identical(ordinary$reached, 1L),
          identical(warning_branch$reached, c(1L, 2L)))
failure_types <- c("domain", "overflow", "allocation", "system", "length", "Rcpp", "unknown")
controlled <- list()
for (stage in 1:2) for (type in seq_along(failure_types)) {
    old <- call(FALSE, y, mu, 2, limit = 2L, stage = stage, type = type)
    new <- call(TRUE, y, mu, 2, limit = 2L, stage = stage, type = type)
    expected_stages <- seq_len(stage)
    stopifnot(identical(old$escaped, ""), identical(old$reached, expected_stages),
              identical(new$reached, expected_stages))
    stage_reference <- call(FALSE, y, mu, 2, limit = 2L, stage = stage, type = 1L)
    stopifnot(identical(old$value, stage_reference$value))
    if (type <= 2L) {
        stopifnot(identical(old, new),
                  identical(new$value[[2L]], if (stage == 1L) 3L else 1L))
        if (stage == 1L) stopifnot(identical(new$value[[1L]], 3/sum((y/mu - 1)^2)))
    } else {
        stopifnot(identical(new$escaped, failure_types[[type]]))
    }
    controlled[[paste(stage, failure_types[[type]], sep = "-")]] <-
        list(original = old, candidate = new)
}
stopifnot(identical(before, vapply(files, sha, character(1))))
headers <- file.path(system.file("include", package = "BH"), c(
    "boost/math/policies/policy.hpp", "boost/math/policies/error_handling.hpp",
    "boost/math/special_functions/digamma.hpp", "boost/math/special_functions/trigamma.hpp"))
config <- system2(file.path(R.home("bin"), "R"), c("CMD", "config", "--all"),
                  stdout = TRUE, stderr = TRUE)
stopifnot(is.null(attr(config, "status")))
config <- config[grepl("^(CC|CXX|CPP|SHLIB|CFLAGS|LDFLAGS|LTO)[A-Z0-9_]*[[:space:]]*=", config)]
report <- list(source_files_sha256 = before, original_commit = "21f9ea098b69c444a884a49bc254de711968e3b5",
    generated_sha256 = vapply(c(original_instrumented, candidate_instrumented, bridge), sha, character(1)),
    header_sha256 = vapply(headers, sha, character(1)), policy_facts = facts,
    compiler_configuration = config,
    compiler_environment = Sys.getenv(c("CPATH", "LIBRARY_PATH", "SDKROOT", "PKG_CPPFLAGS", "PKG_CXXFLAGS", "CXXFLAGS")),
    dependencies = list(BH = as.character(packageVersion("BH")), Rcpp = as.character(packageVersion("Rcpp"))),
    session = sessionInfo(), numerical = numerical, controlled = controlled,
    ordinary = ordinary, warning_branch = warning_branch,
    scope = "Fresh test-only compilation of actual source, not installed-package or universal exception qualification.")
saveRDS(report, file.path(output, "dispersion-exception-results.rds"))
print(facts)
cat("Passed", length(numerical), "exact original/candidate numerical cases and",
    length(controlled), "paired controlled MLE/MM failure cases.\n")
