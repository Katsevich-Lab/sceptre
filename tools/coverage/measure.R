#!/usr/bin/env Rscript
# Run against a fresh isolated installed package, never the user's library.
args <- commandArgs(TRUE)
if (length(args) != 2L) stop("Usage: Rscript tools/coverage/measure.R /source /fresh-output")
root <- normalizePath(args[[1L]], mustWork = TRUE)
out <- normalizePath(args[[2L]], mustWork = FALSE)
if (dir.exists(out) && length(list.files(out, all.files = TRUE, no.. = TRUE))) stop("Output must be fresh")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
library(covr)
# Apple gcov is a wrapper for llvm-cov. Use the matching explicit interface.
if (Sys.info()[["sysname"]] == "Darwin") {
    llvm <- system2("xcrun", c("--find", "llvm-cov"), stdout = TRUE)
    if (length(llvm) == 1L && file.exists(llvm)) options(covr.gcov = llvm, covr.gcov_args = "gcov")
}
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1,
    VECLIB_MAXIMUM_THREADS = 1, MAKEFLAGS = "-j2", NOT_CRAN = "true")
input_files <- function() {
    # Bind runtime/test data, installation controls, native make/configure inputs,
    # and the measurement code itself, not merely handwritten scientific code.
    paths <- unlist(lapply(c("R", "tests", "data", "inst"), function(folder)
        list.files(file.path(root, folder), full.names = TRUE, recursive = TRUE)), use.names = FALSE)
    native <- list.files(file.path(root, "src"), full.names = TRUE, recursive = TRUE)
    native <- native[!grepl("[.](o|so|dll|a|gcda|gcno|gcov)$", native)]
    controls <- file.path(root, c("DESCRIPTION", "NAMESPACE", ".Rbuildignore", ".Rprofile",
        "configure", "configure.win", "cleanup", "cleanup.win",
        "tools/coverage/measure.R", "tools/coverage/render.R"))
    paths <- sort(unique(c(paths, native, controls[file.exists(controls)])))
    stats::setNames(as.list(unname(tools::md5sum(paths))), substring(paths, nchar(root) + 2L))
}
input_hashes <- input_files()
saveRDS(input_hashes, file.path(out, "input-hashes.rds"))
saveRDS(list(commit = system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE),
    dirty = system2("git", c("-C", shQuote(root), "status", "--porcelain"), stdout = TRUE)),
    file.path(out, "source-state.rds"))
lib <- file.path(out, "library")
coverage <- covr::package_coverage(root, type = "tests", quiet = FALSE,
    clean = FALSE, install_path = lib,
    code = paste0("stopifnot(identical(normalizePath(find.package('sceptre')), normalizePath(",
        deparse(file.path(lib, "sceptre")), ")))"))
saveRDS(coverage, file.path(out, "coverage.rds"))
source(file.path(root, "tools/coverage/render.R"))
if (!identical(input_hashes, input_files())) stop("Source/tests changed during measurement; rerun on a frozen checkout")
render_coverage(coverage, root, out)
if (identical(Sys.getenv("SCEPTRE_REQUIRE_NATIVE_COVERAGE"), "true")) {
    measured <- unique(covr::tally_coverage(coverage)$filename)
    native_sources <- paste0("src/", list.files(file.path(root, "src"), pattern = "\\.cpp$"))
    native_sources <- setdiff(native_sources, "src/RcppExports.cpp")
    missing <- setdiff(native_sources, measured)
    if (length(missing) || !any(grepl("^src/", measured)))
        stop("Native instrumentation incomplete: ", paste(missing, collapse = ", "))
}
