#!/usr/bin/env Rscript
# Run separately from package tests: controlled R errors intentionally bypass
# C++ unwinding. No OOM, production hook, or installed sampler replacement.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || length(args) > 2L)
    stop("Supply immutable candidate source and optionally an evidence output directory.")
root <- normalizePath(args[[1L]], mustWork = TRUE)
output <- if (length(args) == 2L) args[[2L]] else
    file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "sceptre-sampling-ownership")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
source <- file.path(root, "src/generate_samples_functions.cpp")
probe <- file.path(root, "tools/refactor/sampling-ownership.cpp")
header <- file.path(root, "src/sampling_safety.h")
runner <- file.path(root, "tools/refactor/test-sampling-ownership.R")
files <- c(source, probe, header, runner)
before <- tools::md5sum(files)
bridge <- tempfile("sceptre-ownership-", fileext = ".cpp")
writeLines(c(sprintf('#define SCEPTRE_SAMPLING_SOURCE %s', encodeString(source, quote = '"')),
    # Load the guard header globally so its R/Rcpp types never enter the private
    # namespace used for the real candidate sampler functions.
    sprintf('#define SCEPTRE_SAMPLING_HEADER %s', encodeString(header, quote = '"')),
    readLines(probe)), bridge)
environment <- new.env(parent = globalenv())
Rcpp::sourceCpp(bridge, env = environment, rebuild = TRUE, showOutput = TRUE)
results <- list()
for (sampler in 0:3) for (phase in 1:3) {
    error <- tryCatch({
        environment$refactor_sampling_pointer_failure(sampler, phase)
        NA_character_
    }, error = conditionMessage)
    state <- environment$refactor_sampling_failure_state()
    stopifnot(identical(error, "Controlled R pointer-boundary failure."),
              identical(state$phase, phase), isTRUE(state$address_was_null))
    results[[length(results) + 1L]] <- list(sampler = sampler, phase = phase,
        error = error, state = state)
}
stopifnot(identical(before, tools::md5sum(files)))
saveRDS(list(source_files_md5 = before, generated_bridge_md5 = tools::md5sum(bridge),
    dependencies = list(Rcpp = as.character(utils::packageVersion("Rcpp")),
                        BH = as.character(utils::packageVersion("BH"))),
    session = utils::sessionInfo(), results = results),
    file.path(output, "sampling-ownership-results.rds"))
stopifnot(all(file.copy(files, output, overwrite = TRUE)))
stopifnot(file.copy(bridge, file.path(output, "sampling-ownership-instrumented.cpp"), overwrite = TRUE))
unlink(bridge)
cat("Passed", length(results), "controlled R-boundary ownership cases.\n")
