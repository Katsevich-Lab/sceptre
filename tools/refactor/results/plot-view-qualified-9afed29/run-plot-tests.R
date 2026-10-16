out <- '/private/tmp/sceptre-plot-view-9afed29'
.libPaths(c('/private/tmp/sceptre-cl2-dense-final-909aa3a/candidate/library',
            '/private/tmp/sceptre-rf0-ondisc-library', .libPaths()))
suppressPackageStartupMessages(library(testthat))
namespace <- loadNamespace('sceptre')
source_env <- new.env(parent = namespace)
sys.source(file.path(out, 'source/R/qq_plot_helpers.R'), envir = source_env)
sys.source(file.path(out, 'source/R/plotting_functions.R'), envir = source_env)
test_env <- new.env(parent = globalenv())
test_env$refactor_private <- function(name) get(name, envir = source_env, inherits = TRUE)
results <- testthat::test_file(file.path(out, 'source/tests/testthat/test-refactor-plot-diagnostic-na.R'),
    env = test_env, reporter = 'summary')
saveRDS(results, file.path(out, 'source-plot-tests.rds'))
counts <- as.data.frame(results)
count_columns <- c('file', 'context', 'test', 'nb', 'failed', 'skipped', 'error', 'warning', 'passed')
counts_scalar <- counts[, intersect(count_columns, names(counts)), drop = FALSE]
utils::write.csv(counts_scalar, file.path(out, 'source-plot-counts.csv'), row.names = FALSE)
print(counts[, c('test', 'nb', 'failed', 'skipped', 'error', 'warning', 'passed')])
cat('Namespace support:', find.package('sceptre'), '\n')
print(sessionInfo())
if (any(counts$failed > 0 | counts$error | counts$warning > 0)) quit(status = 1L)
