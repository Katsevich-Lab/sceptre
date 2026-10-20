# Args and exact checker/options follow the frozen original-wrapper runner.
# Use R_LIBS=/isolated/BiocCheck1.49.30/library and pass archive output 1.49.30 3.24.
# This is the matched submission lane; run unsuppressed diagnostics separately.
source("tools/hardening/submission-cleanup/qualification/run-bioccheck.R")
jsonlite::write_json(list(
    counts = list(errors = length(result$error), warnings = length(result$warning),
        notes = length(result$note)),
    metadata = result$metadata, error = result$error,
    warning = result$warning, note = result$note),
    "diagnostics-summary.json", pretty = TRUE, auto_unbox = TRUE)
stopifnot(length(result$error) == 0L)
