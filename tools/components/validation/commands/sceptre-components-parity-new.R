.libPaths(c('/private/tmp/sceptre-components-lib', .libPaths()))
library(sceptre)
source('tests/testthat/helper-refactor.R')
inputs <- readRDS('tests/testthat/fixtures/refactor/inputs.rds')
cases <- readRDS('tests/testthat/fixtures/refactor/cases.rds')
outputs <- lapply(cases, function(case) refactor_run_case(inputs, case))
saveRDS(outputs, '/private/tmp/sceptre-components-parity-new.rds')
base <- readRDS('/private/tmp/sceptre-components-parity-base.rds')
for(name in names(outputs)) {
 cat(name,':',all.equal(base[[name]], outputs[[name]], tolerance=0),'\n')
}
