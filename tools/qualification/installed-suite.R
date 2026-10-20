# Test the archive installation, retaining diagnostics before failing.
args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
lib <- normalizePath(args[[1L]])
out <- normalizePath(args[[2L]])
backend <- args[[3L]]
.libPaths(c(lib, .libPaths()))
stopifnot(backend %in% c("present", "absent", "either"))
has_backend <- requireNamespace("ondisc", quietly = TRUE)
if (backend != "either") stopifnot(has_backend == (backend == "present"))
stopifnot(identical(normalizePath(find.package("sceptre")),
    normalizePath(file.path(lib, "sceptre"))))
library(sceptre)
writeLines(capture.output(sessionInfo()), file.path(out, "session.txt"))
write.table(installed.packages()[, c("Package", "Version", "LibPath", "Built")],
    file.path(out, "installed.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
jsonlite::write_json(list(R = R.version.string, platform = R.version$platform,
    library = find.package("sceptre"), ondisc = has_backend,
    ondisc_policy = backend, testthat = as.character(packageVersion("testthat"))),
    file.path(out, "environment.json"), pretty = TRUE, auto_unbox = TRUE)
result <- testthat::test_dir(system.file("tests", "testthat", package = "sceptre"),
    package = "sceptre", load_package = "installed", reporter = "summary",
    stop_on_failure = FALSE)
saveRDS(result, file.path(out, "tests.rds"))
d <- as.data.frame(result)
summary <- list(cases = nrow(d), expectations = sum(d$nb), passed = sum(d$passed),
    failed = sum(d$failed), errors = sum(d$error), warnings = sum(d$warning),
    skipped = sum(d$skipped), elapsed_seconds = sum(d$real),
    failures = d[d$failed > 0L | d$error | d$warning > 0L,
        c("file", "test", "failed", "error", "warning"), drop = FALSE],
    skips = d[d$skipped, c("file", "test"), drop = FALSE])
jsonlite::write_json(summary, file.path(out, "tests.json"),
    pretty = TRUE, auto_unbox = TRUE)
stopifnot(summary$failed == 0L, summary$errors == 0L, summary$warnings == 0L)
