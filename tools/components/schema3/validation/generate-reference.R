# Run in a fresh R process from the repository root. The reference library must
# contain the archived flat-slot refactor, not the current component namespace.
reference_library <- Sys.getenv("SCEPTRE_REFERENCE_LIBRARY")
if (!nzchar(reference_library)) stop("Set SCEPTRE_REFERENCE_LIBRARY to the archived reference package library.")
.libPaths(c(reference_library, .libPaths()))
library(sceptre)
validation_dir <- file.path("tools", "components", "schema3", "validation")
source(file.path(validation_dir, "reference-flat-helper.R"))
inputs <- readRDS("tests/testthat/fixtures/refactor/inputs.rds")
cases <- readRDS("tests/testthat/fixtures/refactor/cases.rds")
outputs <- lapply(cases, function(case) withr::with_seed(271828, refactor_capture({
    object <- refactor_prepare(inputs, case)
    prepared <- refactor_observe(object)
    object <- sceptre::run_calibration_check(object, n_calibration_pairs = 4L,
        calibration_group_size = 2L, print_progress = FALSE,
        n_processors = 2L, output_amount = case$output)
    calibrated <- refactor_observe(object)
    # Match the established discovery/power QC budget rather than the legacy
    # accidental calibration overwrite. Numerical reference code is untouched.
    object@B3 <- prepared$B3
    object <- sceptre::run_power_check(object, print_progress = FALSE,
        n_processors = 2L, output_amount = case$output)
    object <- sceptre::run_discovery_analysis(object, print_progress = FALSE,
        n_processors = 2L, output_amount = case$output)
    list(prepared = prepared, calibrated = calibrated, final = refactor_observe(object),
        rng_after = get(".Random.seed", envir = .GlobalEnv))
})))
saveRDS(outputs, file.path(validation_dir, "reference-matched-budgets.rds"))
