commandArgs <- function(trailingOnly = FALSE) "/private/tmp/sceptre-cl1-final-pt2wn9_o/source"
reporter <- testthat::ListReporter$new()
testthat::with_reporter(reporter, {
    reporter$start_file("test-native-components.R")
    source("/private/tmp/sceptre-cl1-final-pt2wn9_o/source/tools/refactor/test-native-components.R")
    reporter$end_file()
})
results <- reporter$get_results()
saveRDS(results, "/private/tmp/sceptre-cl1-final-pt2wn9_o/native-component-test-results.rds")
counts <- as.data.frame(results)
write.csv(counts[, setdiff(names(counts), "result")], "/private/tmp/sceptre-cl1-final-pt2wn9_o/native-component-test-counts.csv", row.names = FALSE)
summary <- list(blocks=nrow(counts), checks=sum(counts$nb), passed=sum(counts$passed), failed=sum(counts$failed), errors=sum(counts$error), warnings=sum(counts$warning), skipped=sum(counts$skipped))
jsonlite::write_json(summary, "/private/tmp/sceptre-cl1-final-pt2wn9_o/native-component-test-summary.json", pretty=TRUE, auto_unbox=TRUE)
stopifnot(summary$checks >= 112L, summary$failed == 0L, summary$errors == 0L, summary$warnings == 0L)
