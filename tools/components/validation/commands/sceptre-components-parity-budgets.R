.libPaths(c('/private/tmp/sceptre-architecture-audit-xqkki89z/library', .libPaths()))
library(sceptre)
source('/private/tmp/sceptre-base-helper.R')
inputs <- readRDS('tests/testthat/fixtures/refactor/inputs.rds')
cases <- readRDS('tests/testthat/fixtures/refactor/cases.rds')
outputs <- lapply(cases, function(case) withr::with_seed(271828, refactor_capture({
 object <- refactor_prepare(inputs, case)
 prepared <- refactor_observe(object)
 object <- sceptre::run_calibration_check(object, n_calibration_pairs=4L,
  calibration_group_size=2L, print_progress=FALSE, n_processors=2L, output_amount=case$output)
 calibrated <- refactor_observe(object)
 # The implementation restores the existing QC policy after calibration. Apply
 # that same budget to the untouched numerical reference, recording raw legacy
 # outputs separately in parity-base.rds rather than overwriting them.
 object@B3 <- prepared$B3
 object <- sceptre::run_power_check(object, print_progress=FALSE,
  n_processors=2L, output_amount=case$output)
 object <- sceptre::run_discovery_analysis(object, print_progress=FALSE,
  n_processors=2L, output_amount=case$output)
 list(prepared=prepared, calibrated=calibrated, final=refactor_observe(object),
  rng_after=get('.Random.seed',envir=.GlobalEnv))
})))
saveRDS(outputs,'/private/tmp/sceptre-components-parity-budgets.rds')
