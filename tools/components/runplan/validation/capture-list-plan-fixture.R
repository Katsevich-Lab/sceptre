# Run from the repository root. Point SCEPTRE_LIST_PLAN_LIBRARY at a library
# containing the package built before RunPlan (the schema-three 2a0938e tree).
# The assertions below prevent overwriting the fixture with typed plans.
old_library <- Sys.getenv("SCEPTRE_LIST_PLAN_LIBRARY", unset =
    "/private/tmp/sceptre-schema3-final-check/sceptre.Rcheck")
.libPaths(c(old_library, .libPaths()))
library(sceptre)
source('tests/testthat/helper-refactor.R')
options(sceptre.refactor.fixture_dir = normalizePath('tests/testthat/fixtures/refactor'))
inputs <- readRDS('tests/testthat/fixtures/refactor/inputs.rds')
case <- readRDS('tests/testthat/fixtures/refactor/cases.rds')$low_complement_perm
object <- suppressMessages(refactor_prepare(inputs, case))
set.seed(20261023)
object <- suppressMessages(run_discovery_analysis(object, print_progress = FALSE, output_amount = 2L))
stopifnot(sceptre:::result_record_is_current(object, 'run_discovery_analysis'))
first <- object@results@records$run_discovery_analysis
stopifnot(is.list(first$provenance$plan), !isS4(first$provenance$plan))
object <- suppressMessages(run_discovery_analysis(object, print_progress = FALSE, output_amount = 2L))
stopifnot(sceptre:::result_record_is_current(object, 'run_discovery_analysis'))
record <- object@results@records$run_discovery_analysis
stopifnot(length(record$history) == 1L, identical(record$history[[1L]]$table, first$table),
          identical(record$history[[1L]]$identity, first$identity))
path <- 'tests/testthat/fixtures/saved-objects/sceptre-schema3-list-plans-2a0938e.rds'
saveRDS(object, path)
roundtrip <- readRDS(path)
stopifnot(sceptre:::result_record_is_current(roundtrip, 'run_discovery_analysis'))
cat('Library:',find.package('sceptre'),'\nFixture:',path,'\nCurrent: TRUE\nHistory records: 1\nPlan representation: list\n')
