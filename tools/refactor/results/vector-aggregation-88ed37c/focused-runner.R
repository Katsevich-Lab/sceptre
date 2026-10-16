.libPaths(c('/private/tmp/sceptre-b6-library-88ed37c', .libPaths()))
library(sceptre)
library(testthat)
data.table::setDTthreads(1L)
root <- '/Users/ekatsevi/.codex/worktrees/refactor-execution-cleanup/sceptre'
source(file.path(root, 'tests/testthat/helper-refactor.R'))
results <- list()
for (name in c('test-refactor-vector-aggregation.R', 'test-refactor-import-coordinates.R', 'test-import_functs.R')) {
  cat('\nFILE', name, '\n')
  results[[name]] <- testthat::test_file(file.path(root, 'tests/testthat', name), reporter = 'summary', stop_on_failure = FALSE)
}
saveRDS(results, '/private/tmp/sceptre-b6-focused-88ed37c.rds')
frames <- lapply(results, as.data.frame)
counts <- data.frame(file = names(frames), tests = vapply(frames, nrow, integer(1)),
  expectations = vapply(frames, function(x) sum(x$nb), numeric(1)),
  failed = vapply(frames, function(x) sum(x$failed), numeric(1)),
  errors = vapply(frames, function(x) sum(x$error), numeric(1)),
  warnings = vapply(frames, function(x) sum(x$warning), numeric(1)),
  skipped = vapply(frames, function(x) sum(x$skipped), numeric(1)))
print(counts)
jsonlite::write_json(counts, '/private/tmp/sceptre-b6-focused-88ed37c.json', pretty=TRUE, auto_unbox=TRUE)
stopifnot(all(counts$failed == 0), all(counts$errors == 0), all(counts$warnings == 0), all(counts$skipped == 0))
a <- readRDS('/private/tmp/sceptre-b6-characterization-reference-88ed37c.rds')
b <- readRDS('/private/tmp/sceptre-b6-characterization-candidate-88ed37c.rds')
stopifnot(identical(a$outputs,b$outputs),identical(a$generated_input_sha256,b$generated_input_sha256),a$input_nonmutation,b$input_nonmutation)
cat('\nEXACT: all 228 installed original/candidate outputs identical; both input nonmutation checks passed.\n')
for(file in c('aggregation-worker.R','aggregation-data.R','aggregation-characterization.R')) parse(file.path(root,'tools/refactor',file))
input <- readRDS(file.path(root,'tools/refactor/results/import-aggregation-characterization-6a43fa5/fixtures.rds'))[[1L]]
input$recipe <- list(dense_input_counterfactual_bytes = prod(as.double(dim(input$m))) * 8)
saveRDS(input, '/private/tmp/sceptre-b6-smoke-input-88ed37c.rds', compress=FALSE)
