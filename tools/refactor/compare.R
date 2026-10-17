#!/usr/bin/env Rscript
# Explicit developer command. Never called by package tests or CI.
# Usage: Rscript tools/refactor/compare.R /fresh/output [candidate-git-ref]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L)
    stop("Usage: compare.R /fresh/output [candidate-git-ref; default HEAD]")
for (package in c("pkgbuild", "pkgload", "testthat", "digest", "jsonlite")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Developer dependency missing: ", package)
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
repo <- normalizePath(file.path(dirname(script), "../.."), winslash = "/", mustWork = TRUE)
out <- path.expand(args[[1L]])
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE)))
    stop("Use a fresh empty output directory; previous runs are never overwritten.")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, winslash = "/", mustWork = TRUE)
if (startsWith(paste0(out, "/"), paste0(repo, "/"))) stop("Write comparison output outside the checkout.")
for (name in c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"))
    do.call(Sys.setenv, stats::setNames(list("1"), name))
Sys.setenv(SCEPTRE_REFACTOR_ORACLE_DEBUG = "false")
reference_pin <- "21f9ea098b69c444a884a49bc254de711968e3b5"
git_read <- function(arguments) {
    result <- system2("git", c("-C", shQuote(repo), arguments), stdout = TRUE, stderr = TRUE)
    if (!is.null(attr(result, "status"))) stop(paste(result, collapse = "\n"))
    result
}
candidate_ref <- if (length(args) == 2L) args[[2L]] else "HEAD"
candidate_pin <- git_read(c("rev-parse", "--verify", shQuote(paste0(candidate_ref, "^{commit}"))))
commands <- list()
report <- list(schema_version = 1L, scope = "Manual refactor consistency only; not permanent CI, performance or whole-consumer qualification.",
    reference_pin = reference_pin, candidate_pin = candidate_pin, requested_ref = candidate_ref,
    working_tree_included = FALSE, output = out, status = "running", commands = commands)
write_report <- function() jsonlite::write_json(report, file.path(out, "comparison.json"),
    pretty = TRUE, auto_unbox = TRUE, null = "null", digits = NA)
run <- function(stage, command, arguments) {
    log <- file.path(out, paste0(stage, ".log"))
    started <- Sys.time()
    status <- system2(command, arguments, stdout = log, stderr = log)
    report$commands[[length(report$commands) + 1L]] <<- list(stage = stage,
        command = command, arguments = arguments, status = status, log = basename(log),
        seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
    write_report()
    if (status != 0L) {
        report$status <<- "failed"
        report$failed_stage <<- stage
        write_report()
        stop(stage, " failed; inspect ", log, call. = FALSE)
    }
}
# Fetch only if this clone does not yet contain the exact historical object.
has_reference <- system2("git", c("-C", shQuote(repo), "cat-file", "-e", reference_pin),
    stdout = FALSE, stderr = FALSE) == 0L
if (!has_reference) run("fetch-original", "git", c("-C", shQuote(repo), "fetch", "--no-tags", "origin", reference_pin))
for (role in c("original", "candidate")) {
    pin <- if (role == "original") reference_pin else candidate_pin
    archive <- file.path(out, paste0(role, ".tar"))
    run(paste0("archive-", role), "git", c("-C", shQuote(repo), "-c", "core.autocrlf=false",
        "-c", "core.eol=lf", "archive", "--format=tar", paste0("--output=", shQuote(archive)), pin))
    directory <- file.path(out, role)
    dir.create(directory)
    utils::untar(archive, exdir = directory)
}
candidate <- file.path(out, "candidate")
suite <- file.path(candidate, "tools/refactor/consistency/testthat")
if (!dir.exists(suite)) stop("Candidate ref predates the developer consistency suite; choose a newer committed ref.")
library_dir <- file.path(out, "library")
dir.create(library_dir)
run("install-candidate", file.path(R.home("bin"), "R"), c("CMD", "INSTALL", "--preclean",
    paste0("--library=", shQuote(library_dir)), shQuote(candidate)))
oracle <- file.path(out, "oracle")
run("run-original", file.path(R.home("bin"), "Rscript"), c("--vanilla",
    shQuote(file.path(candidate, "tools/refactor/generate-runtime-oracle.R")),
    shQuote(file.path(out, "original")), shQuote(oracle)))
Sys.setenv(SCEPTRE_REFACTOR_ORACLE = oracle)
# All generated R literals use dput(), preserving spaces/backslashes safely.
literal <- function(x) paste(capture.output(dput(x)), collapse = "\n")
worker <- file.path(out, "compare-installed.R")
writeLines(c(
    paste0(".libPaths(c(", literal(library_dir), ", .libPaths()))"),
    "library(sceptre)",
    paste0("stopifnot(identical(normalizePath(find.package('sceptre')), normalizePath(", literal(file.path(library_dir, "sceptre")), "))))"),
    "data.table::setDTthreads(1L)",
    paste0("results <- testthat::test_dir(", literal(suite), ", package='sceptre', load_package='installed', reporter='summary', stop_on_failure=FALSE, stop_on_warning=FALSE)"),
    paste0("saveRDS(results, ", literal(file.path(out, "results.rds")), ")"),
    "counts <- as.data.frame(results)",
    paste0("write.csv(counts[,setdiff(names(counts),'result')], ", literal(file.path(out, "checks.csv")), ", row.names=FALSE)"),
    "summary <- list(blocks=nrow(counts), checks=sum(counts$nb), passed=sum(counts$passed), failed=sum(counts$failed), errors=sum(counts$error), warnings=sum(counts$warning), skipped=sum(counts$skipped))",
    paste0("jsonlite::write_json(summary, ", literal(file.path(out, "counts.json")), ", pretty=TRUE, auto_unbox=TRUE)"),
    "if (summary$failed || summary$errors || summary$warnings) quit(status=1L)"
), worker)
run("compare", file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(worker)))
report$counts <- jsonlite::read_json(file.path(out, "counts.json"), simplifyVector = TRUE)
report$status <- "passed"
report$artifacts <- list(log = "compare.log", detailed_results = "results.rds", checks = "checks.csv")
report$limits <- c("A recorded platform-specific stream test can skip; its full same-environment source comparison still runs.",
    "Only supported historical domains are compared; reviewed bug-fix regressions remain in permanent package tests.",
    "No resource, actual Nextflow or whole-book claim follows from this command.")
write_report()
writeLines(c("# Refactor consistency comparison", "",
    paste("Original:", reference_pin), paste("Candidate:", candidate_pin), "",
    paste("Result:", report$status), paste("Passing expectations:", report$counts$passed),
    paste("Platform-qualified skips:", report$counts$skipped), "",
    "See checks.csv for per-contract results and compare.log/results.rds for test differences.",
    "No expected output was regenerated from candidate source. This command is not CI."),
    file.path(out, "README.md"))
cat("Comparison passed. Read ", file.path(out, "README.md"), " and comparison.json\n", sep = "")
