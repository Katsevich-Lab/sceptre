# Run at the package root. Optionally select a separately installed candidate.
root <- normalizePath('.')
output <- Sys.getenv('SCEPTRE_TEST_OUTPUT', file.path(root, 'tools/hardening/validation'))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
library_path <- Sys.getenv('SCEPTRE_TEST_LIBRARY')
if (nzchar(library_path)) {
    .libPaths(c(normalizePath(library_path), .libPaths()))
    library(sceptre)
    mode <- 'installed'
} else {
    pkgload::load_all(root, recompile = FALSE)
    mode <- 'none'
}
result <- testthat::test_dir(file.path(root, 'tests/testthat'), package = 'sceptre',
    load_package = mode, reporter = 'summary', stop_on_failure = FALSE)
saveRDS(result, file.path(output, 'full-tests.rds'))
rows <- as.data.frame(result)
summary <- list(cases = nrow(rows), passed = sum(rows$passed), failed = sum(rows$failed),
    errors = sum(rows$error), warnings = sum(rows$warning), skips = sum(rows$skipped),
    R = R.version.string, platform = R.version$platform,
    testthat = as.character(packageVersion('testthat')), package_path = find.package('sceptre'),
    mode = mode)
jsonlite::write_json(summary, file.path(output, 'full-tests-summary.json'),
    pretty = TRUE, auto_unbox = TRUE)
print(summary)
stopifnot(summary$failed == 0L, summary$errors == 0L, summary$warnings == 0L)
