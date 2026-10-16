#!/usr/bin/env Rscript
# Temporary same-environment expectations from the original package. This never
# writes committed fixtures and never makes inputs or runs candidate source.
# Usage: Rscript --vanilla tools/refactor/generate-runtime-oracle.R <reference> <output>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Supply an archived pinned reference source and a fresh output directory.")
if ("sceptre" %in% loadedNamespaces()) stop("Run the reference generator in a fresh R process, before loading sceptre.")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
repo <- normalizePath(file.path(dirname(script), "../.."), mustWork = TRUE)
reference <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
output <- args[[2L]]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output, winslash = "/", mustWork = TRUE)
fixture <- normalizePath(file.path(repo, "tests/testthat/fixtures/refactor"), winslash = "/")
if (identical(output, fixture) || length(list.files(output, all.files = TRUE, no.. = TRUE))) {
    stop("Runtime oracle output must be a fresh directory separate from committed fixtures.")
}
helper <- file.path(repo, "tests/testthat/helper-refactor.R")
source(helper)
options(sceptre.refactor.fixture_dir = fixture)
manifest <- readRDS(file.path(fixture, "manifest.rds"))
baseline <- "21f9ea098b69c444a884a49bc254de711968e3b5"
if (!identical(manifest$source_commit, baseline) || !isTRUE(manifest$reference_source_verified)) {
    stop("The immutable fixture manifest is not the supported pinned baseline.")
}
input_files <- refactor_oracle_input_files()
sha256 <- function(files) stats::setNames(vapply(files, function(file) {
    digest::digest(file = file, algo = "sha256")
}, character(1)), names(files))
input_paths <- stats::setNames(file.path(fixture, input_files), input_files)
input_hashes <- sha256(input_paths)
if (!identical(input_hashes, unlist(manifest$artifacts_sha256[input_files], use.names = TRUE))) {
    stop("Immutable runtime reference inputs failed SHA256 verification.")
}
source_paths <- stats::setNames(file.path(reference, names(manifest$source_files)), names(manifest$source_files))
source_md5 <- refactor_oracle_hashes(source_paths)
verify_source <- function() {
    refactor_oracle_verify_source(reference, manifest$source_files, source_md5)
    actual <- sha256(source_paths)
    expected <- unlist(manifest$source_files, use.names = TRUE)
    if (!identical(actual, expected)) {
        mismatched <- names(expected)[is.na(actual) | actual != expected]
        detail <- vapply(mismatched, function(file) {
            paste0(file, ": expected ", expected[[file]], ", actual ", actual[[file]])
        }, character(1))
        stop("Runtime reference differs from immutable pinned source SHA256 hashes:\n",
             paste(detail, collapse = "\n"))
    }
}
verify_source()
debug_setting <- Sys.getenv("SCEPTRE_REFACTOR_ORACLE_DEBUG", "true")
if (!debug_setting %in% c("true", "false")) stop("SCEPTRE_REFACTOR_ORACLE_DEBUG must be true or false.")
debug <- identical(debug_setting, "true")
# The pinned archive already contains its Rcpp registration/wrappers. Automatic
# compileAttributes() rewrites those source files (including platform-specific
# line endings), so compile their frozen bytes without regenerating them.
pkgbuild::clean_dll(reference)
pkgbuild::compile_dll(reference, force = TRUE, compile_attributes = FALSE,
                      quiet = TRUE, debug = debug)
verify_source()
pkgload::load_all(reference, export_all = FALSE, helpers = FALSE, quiet = TRUE,
                  compile = FALSE, debug = debug)
if (!identical(normalizePath(getNamespaceInfo("sceptre", "path"), winslash = "/"), reference)) {
    stop("Loaded package namespace is not the verified original source archive.")
}
dll <- getNamespaceInfo("sceptre", "DLLs")[["sceptre"]]
if (is.null(dll) || !startsWith(normalizePath(dll[["path"]], winslash = "/", mustWork = TRUE), paste0(reference, "/"))) {
    stop("Loaded native library is not from the cleanly compiled original source archive.")
}
verify_source()
inputs <- readRDS(file.path(fixture, "inputs.rds"))
cases <- readRDS(file.path(fixture, "cases.rds"))
backed_inputs <- readRDS(file.path(fixture, "backed-inputs.rds"))
save <- function(value, name) saveRDS(value, file.path(output, name), version = 3)
save(lapply(cases, function(case) refactor_run_case(inputs, case)), "workflows.rds")
save(refactor_cache_observations(inputs, cases$low_complement_perm), "cache.rds")
save(refactor_edge_observations(inputs, cases$low_complement_perm), "edges.rds")
saved <- readRDS(file.path(fixture, "saved-object.rds"))
continued <- refactor_capture(sceptre::run_discovery_analysis(saved, print_progress = FALSE, output_amount = 2L))
save(list(value = refactor_observe(continued$value), warnings = continued$warnings), "saved-object-result.rds")
save(refactor_run_case(backed_inputs, cases$high_perm), "backed-workflow.rds")
save(refactor_capture({
    object <- refactor_prepare(backed_inputs, cases$high_perm)
    refactor_observe(sceptre::run_discovery_analysis(object, print_progress = FALSE, output_amount = 3L))
}), "backed-discovery.rds")
has_backed <- requireNamespace("ondisc", quietly = TRUE)
if (has_backed) {
    directory <- tempfile("refactor-runtime-backed-")
    backed <- refactor_run_case(backed_inputs, cases$high_perm, backed = TRUE, directory = directory)
    save(backed, "backed-observations.rds")
    save(list(ondisc_version = as.character(utils::packageVersion("ondisc")),
              warnings = backed$warnings, rng_after = backed$value$rng_after), "backed-metadata.rds")
    unlink(directory, recursive = TRUE)
}
# Verify again after execution, not merely before compilation. The artifact
# manifest is written last, so an interrupted/partial run is never accepted.
verify_source()
if (!identical(input_hashes, sha256(input_paths))) {
    stop("Reference execution changed immutable inputs.")
}
outputs <- refactor_oracle_output_files(has_backed)
runtime <- list(
    schema_version = 1L, type = "sceptre-same-environment-reference", source_commit = baseline,
    source_verified_before_and_after = TRUE, clean_native_rebuild = TRUE, fresh_process = TRUE,
    reference_source = reference, source_files = manifest$source_files, source_md5 = source_md5,
    process_id = Sys.getpid(), debug_compile = debug,
    inputs_sha256 = input_hashes, inputs_md5 = refactor_oracle_hashes(input_paths),
    harness_sha256 = unname(sha256(helper)), harness_md5 = unname(refactor_oracle_hashes(helper)),
    environment = refactor_oracle_environment(),
    outputs_sha256 = sha256(stats::setNames(file.path(output, outputs), outputs)),
    outputs_md5 = refactor_oracle_hashes(stats::setNames(file.path(output, outputs), outputs)),
    policy = "Same-environment original source reference; frozen inputs/API/maps/native oracle unchanged; tolerance 1e-8 and exact identities/decisions unchanged."
)
save(runtime, "runtime-manifest.rds")
cat("Built verified temporary workflow oracle from", baseline, "with debug =", debug, "in", output, "\n")
