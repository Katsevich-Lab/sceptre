#!/usr/bin/env Rscript
# Check the current estimator's declared exception policy, without old source.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply current immutable source and fresh evidence directory.")
root <- normalizePath(args[[1L]], mustWork = TRUE)
dir.create(args[[2L]], recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(args[[2L]], mustWork = TRUE)
if (length(list.files(output, all.files = TRUE, no.. = TRUE)))
    stop("Use a fresh evidence directory; previous attempts are never overwritten.")
candidate <- file.path(root, "src/compute_nb_size_param_functions.cpp")
probe <- file.path(root, "tools/checks/dispersion-exceptions.cpp")
runner <- file.path(root, "tools/checks/test-dispersion-exceptions.R")
files <- c(candidate, probe, runner)
before <- tools::md5sum(files)
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

instrumented <- file.path(output, "estimator-instrumented.inc")
instrument(candidate, instrumented)
bridge <- file.path(output, "dispersion-instrumented.cpp")
writeLines(c(sprintf('#define SCEPTRE_CANDIDATE_DISPERSION %s', encodeString(instrumented, quote = '"')),
             readLines(probe)), bridge)
environment <- new.env(parent = globalenv())
Rcpp::sourceCpp(bridge, env = environment, rebuild = TRUE, showOutput = TRUE)
facts <- environment$refactor_dispersion_policy_facts()
stopifnot(facts$double_radix == 2L, facts$double_digits == 53L,
          facts$effective_precision == 53L,
          facts$domain == facts$throw_on_error, facts$pole == facts$throw_on_error,
          facts$overflow == facts$throw_on_error,
          facts$underflow == facts$ignore_error, facts$denorm == facts$ignore_error)


call <- function(y, mu, dfr, limit = 50L, stage = 0L, type = 0L) {
    environment$refactor_dispersion_case(y, mu, dfr, limit,
        .Machine$double.eps^(1/4), stage, type)
}
y <- c(0, 1, 5)
mu <- rep(1, 3)
ordinary <- call(y, mu, 2)
warning_branch <- call(y, mu, 2, limit = 2L)
stopifnot(identical(ordinary$reached, 1L),
          identical(warning_branch$reached, c(1L, 2L)))
types <- c("domain", "overflow", "allocation", "system", "length", "Rcpp", "unknown")
results <- list()
for (stage in 1:2) for (type in seq_along(types)) {
    actual <- call(y, mu, 2, limit = 2L, stage = stage, type = type)
    stopifnot(identical(actual$reached, seq_len(stage)))
    if (type <= 2L) {
        stopifnot(identical(actual$escaped, ""),
                  identical(actual$value[[2L]], if (stage == 1L) 3L else 1L))
        if (stage == 1L) stopifnot(identical(actual$value[[1L]], 3/sum((y/mu - 1)^2)))
        domain <- call(y, mu, 2, limit = 2L, stage = stage, type = 1L)
        stopifnot(identical(actual$value, domain$value))
    } else stopifnot(identical(actual$escaped, types[[type]]))
    results[[paste(stage, types[[type]])]] <- actual
}
stopifnot(identical(before, tools::md5sum(files)))
saveRDS(list(scope="Current estimator exception policy; no historical equality check.",
             source_files_md5=before, policy=facts, results=results),
        file.path(output, "dispersion-policy.rds"))
cat("Passed", length(results), "controlled current-estimator policy cases.\n")
