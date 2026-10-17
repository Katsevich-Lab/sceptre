#!/usr/bin/env Rscript
# Historical bootstrap, no longer invoked by CI. Keep both reference source and results outside the candidate
# checkout and outside R's session tempdir, which disappears when R exits.
args <- commandArgs(trailingOnly = TRUE)
persistent <- if (length(args)) args[[1L]] else Sys.getenv("RUNNER_TEMP")
if (!nzchar(persistent)) stop("Supply a persistent temporary root or RUNNER_TEMP.")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
repo <- normalizePath(file.path(dirname(script), "../.."), winslash = "/", mustWork = TRUE)
dir.create(persistent, recursive = TRUE, showWarnings = FALSE)
root <- tempfile("sceptre-runtime-", tmpdir = persistent)
dir.create(root)
reference <- file.path(root, "original")
output <- file.path(root, "oracle")
dir.create(reference)
baseline <- "21f9ea098b69c444a884a49bc254de711968e3b5"
git <- function(arguments) {
    status <- system2("git", c("-C", shQuote(repo), arguments))
    if (!identical(status, 0L)) stop("Cannot fetch/archive the pinned original package.")
}
git(c("fetch", "--no-tags", "origin", baseline))
archive <- file.path(root, "original.tar")
git(c("-c", "core.autocrlf=false", "-c", "core.eol=lf", "archive", "--format=tar",
      paste0("--output=", shQuote(archive)), baseline))
utils::untar(archive, exdir = reference)
status <- system2(file.path(R.home("bin"), "Rscript"), c(
    "--vanilla", shQuote(file.path(repo, "tools/refactor/generate-runtime-oracle.R")),
    shQuote(reference), shQuote(output)
))
if (!identical(status, 0L)) stop("The separate-process original-package oracle failed.")
output <- normalizePath(output, winslash = "/", mustWork = TRUE)
reference <- normalizePath(reference, winslash = "/", mustWork = TRUE)
env_file <- Sys.getenv("GITHUB_ENV")
if (nzchar(env_file)) {
    cat("SCEPTRE_REFACTOR_ORACLE=", output, "\n", sep = "", file = env_file, append = TRUE)
    cat("SCEPTRE_REFACTOR_ORIGINAL_SOURCE=", reference, "\n", sep = "", file = env_file, append = TRUE)
}
cat("SCEPTRE_REFACTOR_ORACLE=", output, "\n", sep = "")
cat("SCEPTRE_REFACTOR_ORIGINAL_SOURCE=", reference, "\n", sep = "")
