# Fresh process, run from the copied package root. Source must be installed in
# the isolated floor library first. This script never installs a dependency.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- normalizePath(args[[1L]], mustWork = TRUE)
lib <- normalizePath(Sys.getenv("SCEPTRE_FLOOR_LIBRARY"), mustWork = TRUE)
minor <- Sys.getenv("SCEPTRE_EXPECTED_R_MINOR", unset = "4.4")
stopifnot(startsWith(as.character(getRversion()), paste0(minor, ".")),
    identical(as.character(packageVersion("testthat")), "3.1.8"),
    startsWith(normalizePath(find.package("testthat")), paste0(lib, "/")),
    startsWith(normalizePath(find.package("sceptre")), paste0(lib, "/")))
library(sceptre)
results <- testthat::test_dir("tests/testthat", package = "sceptre",
    load_package = "installed", stop_on_failure = FALSE, reporter = "summary")
saveRDS(results, file.path(out, "floor-suite-results.rds"))
d <- as.data.frame(results)
summary <- list(cases = nrow(d), expectations = sum(d$nb), passed = sum(d$passed),
    failed = sum(d$failed), errors = sum(d$error), warnings = sum(d$warning),
    skipped = sum(d$skipped), elapsed_seconds = sum(d$real),
    failures = d[d$failed > 0L | d$error | d$warning > 0L,
        c("file", "test", "failed", "error", "warning"), drop = FALSE],
    skips = d[d$skipped, c("file", "test"), drop = FALSE])
jsonlite::write_json(summary, file.path(out, "floor-suite-summary.json"),
    pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out, "floor-suite-session.txt"))
write.table(installed.packages()[, c("Package", "Version", "LibPath", "Built")],
    file.path(out, "floor-installed-packages.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summary[setdiff(names(summary), c("failures", "skips"))])
stopifnot(summary$failed == 0L, summary$errors == 0L, summary$warnings == 0L)
