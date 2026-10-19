# Run in a fresh R process from the repository root.
pkgload::load_all(".", recompile = FALSE)
validation_dir <- file.path("tools", "components", "runplan", "validation")
results <- testthat::test_dir("tests/testthat", package = "sceptre",
    load_package = "none", stop_on_failure = FALSE, reporter = "summary")
saveRDS(results, file.path(validation_dir, "source-suite-results.rds"))
source(file.path(validation_dir, "summarize-source-suite.R"))
