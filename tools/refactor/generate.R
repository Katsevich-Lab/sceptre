#!/usr/bin/env Rscript
# Deliberate baseline generation. Never called by package tests or CI.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript tools/refactor/generate.R <archived-reference-source>")
repo <- normalizePath(".")
reference <- normalizePath(args[[1L]])
baseline <- "21f9ea098b69c444a884a49bc254de711968e3b5"
production <- c("R", "src", "data", "DESCRIPTION", "NAMESPACE")
files <- system2("git", c("ls-tree", "-r", "--name-only", baseline, "--", production), stdout = TRUE)
if (!length(files)) stop("Cannot enumerate pinned baseline files")
for (file in files) {
    expected <- system2("git", c("rev-parse", paste0(baseline, ":", file)), stdout = TRUE)
    actual <- system2("git", c("hash-object", shQuote(file.path(reference, file))), stdout = TRUE)
    if (!identical(actual, expected)) stop("Reference differs from pinned source: ", file)
}
for (directory in c("R", "src", "data")) {
    actual <- list.files(file.path(reference, directory), recursive = TRUE, full.names = FALSE)
    actual <- actual[!grepl("\\.(o|so|dll|dylib|d)$", actual)]
    expected <- substring(files[startsWith(files, paste0(directory, "/"))], nchar(directory) + 2L)
    if (!setequal(actual, expected)) stop("Unexpected reference files under ", directory)
}
pkgload::load_all(reference, export_all = FALSE, helpers = FALSE, quiet = TRUE,
                  compile = TRUE, recompile = TRUE)
for (file in files) {
    expected <- system2("git", c("rev-parse", paste0(baseline, ":", file)), stdout = TRUE)
    actual <- system2("git", c("hash-object", shQuote(file.path(reference, file))), stdout = TRUE)
    if (!identical(actual, expected)) stop("Compilation changed pinned source: ", file)
}
source(file.path(repo, "tests/testthat/helper-refactor.R"))
output <- file.path(repo, "tests/testthat/fixtures/refactor")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
options(sceptre.refactor.fixture_dir = output)
inputs <- refactor_make_inputs()
cases <- refactor_case_definitions()
saveRDS(inputs, file.path(output, "inputs.rds"), version = 3)
saveRDS(cases, file.path(output, "cases.rds"), version = 3)
saveRDS(refactor_api(), file.path(output, "api.rds"), version = 3)
observations <- lapply(cases, function(case) refactor_run_case(inputs, case))
saveRDS(observations, file.path(output, "workflows.rds"), version = 3)
saveRDS(refactor_index_observations(), file.path(output, "indexing.rds"), version = 3)
saveRDS(refactor_cache_observations(inputs, cases$low_complement_perm), file.path(output, "cache.rds"), version = 3)
saveRDS(refactor_edge_observations(inputs, cases$low_complement_perm), file.path(output, "edges.rds"), version = 3)
saved <- refactor_capture(refactor_prepare(inputs, cases$low_complement_perm))$value
saveRDS(saved, file.path(output, "saved-object.rds"), version = 3)
result <- refactor_capture(sceptre::run_discovery_analysis(saved, print_progress = FALSE, output_amount = 2L))
saveRDS(list(value = refactor_observe(result$value), warnings = result$warnings), file.path(output, "saved-object-result.rds"), version = 3)
backed_inputs <- refactor_expand_inputs(inputs)
saveRDS(backed_inputs, file.path(output, "backed-inputs.rds"), version = 3)
saveRDS(refactor_run_case(backed_inputs, cases$high_perm), file.path(output, "backed-workflow.rds"), version = 3)
backed_discovery <- refactor_capture({
    object <- refactor_prepare(backed_inputs, cases$high_perm)
    refactor_observe(sceptre::run_discovery_analysis(object, print_progress = FALSE, output_amount = 3L))
})
saveRDS(backed_discovery, file.path(output, "backed-discovery.rds"), version = 3)
if (requireNamespace("ondisc", quietly = TRUE)) {
    directory <- tempfile("refactor-baseline-ondisc-")
    backed <- refactor_run_case(backed_inputs, cases$high_perm, backed = TRUE, directory = directory)
    saveRDS(list(ondisc_version = as.character(utils::packageVersion("ondisc")),
                 warnings = backed$warnings, rng_after = backed$value$rng_after),
            file.path(output, "backed-metadata.rds"), version = 3)
    unlink(directory, recursive = TRUE)
}

write.table(inputs$response_matrix, file.path(output, "response-counts.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(inputs$low_grna_matrix, file.path(output, "low-guide-counts.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(inputs$high_grna_matrix, file.path(output, "high-guide-counts.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(inputs$extra_covariates, file.path(output, "covariates.tsv"), sep = "\t", quote = FALSE, col.names = NA)
jsonlite::write_json(refactor_api(), file.path(output, "api.json"), pretty = TRUE, auto_unbox = TRUE)
deps <- sort(unique(c("sceptre", "ondisc", "Matrix", "Rcpp", "BH", "dplyr", "data.table", "testthat", "pkgload", "withr")))
versions <- stats::setNames(lapply(deps, function(p) {
    if (requireNamespace(p, quietly = TRUE)) list(version = as.character(utils::packageVersion(p)), path = find.package(p)) else NULL
}), deps)
manifest <- list(
    schema_version = 1L, source_commit = baseline, generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    reference_source_verified = TRUE, clean_native_rebuild = TRUE,
    source_files = stats::setNames(lapply(files, function(f) digest::digest(file = file.path(reference, f), algo = "sha256")), files),
    R = R.version.string, platform = R.version$platform, os = Sys.info(), external = extSoftVersion(),
    dependencies = versions, RNGkind = RNGkind(), input_seed = 20261016L, workflow_seed = 271828L,
    native_rng = "Boost MT19937(4); production sampler families/draw order unchanged",
    indices = list(public_and_R = "1-based", native_draw_arrays = "0-based"),
    comparison = list(floating_tolerance = 1e-8, exact = c("API/formals/methods/slots", "memberships/order", "budgets", "QC/significance/stage decisions")),
    resource_validation = "Small characterization fixtures only; realistic performance and process-tree RSS gates outstanding"
)
artifact_files <- sort(list.files(output, full.names = TRUE))
artifact_files <- artifact_files[!grepl("manifest\\.(json|rds)$", artifact_files)]
manifest$artifacts_sha256 <- stats::setNames(lapply(artifact_files, function(f) digest::digest(file = f, algo = "sha256")), basename(artifact_files))
manifest$artifacts_md5 <- stats::setNames(unname(tools::md5sum(artifact_files)), basename(artifact_files))
saveRDS(manifest, file.path(output, "manifest.rds"), version = 3)
jsonlite::write_json(manifest, file.path(output, "manifest.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")
cat("Generated explicit reference fixtures from", baseline, "\n")
