library(testthat)
validation_dir <- file.path("tools", "components", "fits", "validation")
results <- readRDS(file.path(validation_dir, "source-suite-results.rds"))
d <- as.data.frame(results)
summary <- list(cases = nrow(d), expectations = sum(d$nb),
    failed = sum(d$failed), errors = sum(d$error), warnings = sum(d$warning),
    skipped = sum(d$skipped), passed = sum(d$passed),
    elapsed_seconds = sum(d$real),
    failures = d[d$failed > 0L | d$error, c("file", "test", "failed", "error"), drop = FALSE])
jsonlite::write_json(summary, file.path(validation_dir, "source-suite-summary.json"),
    pretty = TRUE, auto_unbox = TRUE)
print(summary[setdiff(names(summary), "failures")])
