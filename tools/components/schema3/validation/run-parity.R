# Run in a fresh R process from the repository root.
pkgload::load_all(".", recompile = FALSE)
source("tests/testthat/helper-refactor.R")
validation_dir <- file.path("tools", "components", "schema3", "validation")
inputs <- readRDS("tests/testthat/fixtures/refactor/inputs.rds")
cases <- readRDS("tests/testthat/fixtures/refactor/cases.rds")
outputs <- lapply(cases, function(case) refactor_run_case(inputs, case))
saveRDS(outputs, file.path(validation_dir, "parity-current.rds"))
source(file.path(validation_dir, "compare-parity.R"))
