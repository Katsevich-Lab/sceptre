#!/usr/bin/env Rscript
# Run in a separate R process: controlled R errors can skip C++ destructors.
# Test cleanup without exhausting memory or changing installed package functions.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || length(args) > 2L)
    stop("Supply immutable candidate source and optionally an evidence output directory.")
root <- normalizePath(args[[1L]], mustWork = TRUE)
output <- if (length(args) == 2L) args[[2L]] else
    file.path(Sys.getenv("RUNNER_TEMP", tempdir()), "sceptre-count-ownership")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(output, mustWork = TRUE)
paths <- file.path(root, c("src/sparse_counts.cpp", "src/sparse_matrix_functs.cpp",
    "src/RcppExports.cpp", "tools/checks/count-ownership.cpp", "tools/checks/test-count-ownership.R"))
before <- tools::md5sum(paths)
count <- readLines(paths[[1L]])
covariate <- readLines(paths[[2L]])
# Substitute only the test wrappers, retaining the actual arithmetic and buffer
# lifetimes. Use a smaller integer limit to reach the overflow error without
# allocating billions of members; ordinary cases still use INT_MAX.
needle <- "if (nt_total > std::numeric_limits<int>::max())"
stopifnot(sum(grepl(needle, count, fixed = TRUE)) == 1L)
count <- sub(needle, "if (nt_total > refactor_count_overflow_limit)", count, fixed = TRUE)
transform <- function(lines) gsub("std::vector<", "refactor_count_vector<", lines, fixed = TRUE)
count_path <- file.path(output, "sparse_counts-instrumented.inc")
covariate_path <- file.path(output, "sparse_matrix-instrumented.inc")
writeLines(transform(count), count_path)
writeLines(transform(covariate), covariate_path)
bridge <- file.path(output, "count-ownership-instrumented.cpp")
writeLines(c(sprintf('#define SCEPTRE_COUNT_SOURCE %s', encodeString(count_path, quote = '"')),
    sprintf('#define SCEPTRE_COVARIATE_SOURCE %s', encodeString(covariate_path, quote = '"')),
    readLines(paths[[4L]])), bridge)
environment <- new.env(parent = globalenv())
Rcpp::sourceCpp(bridge, env = environment, rebuild = TRUE, showOutput = TRUE)
inputs <- function() list(c(0L, 2L, 1L, 3L), c(0L, 2L, 4L),
    list(1:4, c(4, 2), integer()), list(1:2, 2:4), 1:2, 2:1, 1:4,
    0:3, c(0L, 4L), as.numeric(1:4), 1:4)
expected <- list(list(n_nonzero_mat = matrix(c(1L, 1L, 1L, 2L), 2L),
    n_nonzero_tot = c(2L, 3L), n_nonzero_trt = c(0L, 2L), n_nonzero_cntrl = c(2L, 3L)),
    matrix(c(2L, 0L, 0L, 2L, 2L, 0L), 3L),
    list(n_nonzero = 4L, n_umi = 10, p_mito = 1, max_feature = 3L, frac_umis_max_feature = 0.4))
results <- list()
record <- function(kernel, phase, position, error, state) {
    stopifnot(identical(state$live_blocks, 0L), identical(state$unsafe_boundaries, 0L),
              all(state$boundary_live_blocks == 0L))
    results[[length(results) + 1L]] <<- list(kernel = kernel, phase = phase,
        position = position, error = error, state = state)
}
for (kernel in 0:2) {
    # Run once so Rcpp caches are initialized. Later numbered allocation sites
    # can then be tested without counting one-time setup calls.
    environment$refactor_count_ownership_call(kernel, inputs(), 0L, 0L, 0L)
    environment$refactor_count_ownership_state()
    original <- inputs()
    snapshot <- serialize(original, NULL)
    actual <- environment$refactor_count_ownership_call(kernel, original, 0L, 0L, 0L)
    baseline <- environment$refactor_count_ownership_state()
    stopifnot(identical(actual, expected[[kernel + 1L]]),
              identical(serialize(original, NULL), snapshot), baseline$native_count > 0L,
              baseline$max_live_blocks > 0L)
    record(kernel, "success", 0L, NA_character_, baseline)
    for (boundary in seq_len(baseline$boundary_count)) {
        error <- tryCatch({
            environment$refactor_count_ownership_call(kernel, inputs(), boundary, 0L, 0L)
            NA_character_
        }, error = conditionMessage)
        state <- environment$refactor_count_ownership_state()
        stopifnot(identical(error, "Controlled count R-boundary failure."),
                  identical(state$failed_boundary, boundary))
        record(kernel, "raw-R-error", boundary, error, state)
    }
    for (kind in 0:1) for (position in seq_len(baseline$native_count)) {
        error <- tryCatch({
            environment$refactor_count_ownership_call(kernel, inputs(), 0L, position, kind)
            NA_character_
        }, error = conditionMessage)
        state <- environment$refactor_count_ownership_state()
        stopifnot(!is.na(error), identical(state$failed_native, position))
        record(kernel, c("bad-alloc", "native-overflow")[[kind + 1L]], position, error, state)
    }
}
error <- tryCatch({
    environment$refactor_count_ownership_call(0L, inputs(), 0L, 0L, 0L, TRUE)
    NA_character_
}, error = conditionMessage)
state <- environment$refactor_count_ownership_state()
stopifnot(identical(error, "NT group sum exceeds integer count capacity."), state$max_live_blocks > 0L)
record(0L, "actual-NT-overflow-branch-bounded-threshold", 0L, error, state)
observed <- unique(unlist(lapply(results, function(x) x$state$boundary_kinds)))
saveRDS(list(results = results, observed_boundary_kinds = sort(observed)),
    file.path(output, "count-ownership-partial-results.rds"))
cat("Observed R boundary kinds:", paste(sort(observed), collapse = ", "), "\n")
# Depending on Rcpp version, a matrix is made with a vector plus dimensions or
# Rf_allocMatrix(). Wrap both routes and require errors at the route actually used.
stopifnot(all(c(1L, 3L, 4L, 5L, 6L, 7L, 8L, 9L) %in% observed),
          any(c(2L, 11L) %in% observed),
          any(vapply(results, function(x) x$phase == "bad-alloc" && x$state$max_live_blocks > 0L, logical(1))))
negative_error <- tryCatch({
    environment$refactor_count_ownership_negative_control()
    NA_character_
}, error = conditionMessage)
negative_state <- environment$refactor_count_ownership_state()
stopifnot(identical(negative_error, "R boundary observed with live native count ownership."),
          identical(negative_state$unsafe_boundaries, 1L),
          identical(negative_state$live_blocks, 0L), negative_state$max_live_blocks > 0L)
stopifnot(identical(before, tools::md5sum(paths)))
saveRDS(list(source_files_md5 = before, instrumented_files_md5 = tools::md5sum(c(bridge, count_path, covariate_path)),
    transformations = list(vector_allocator = "std::vector< -> refactor_count_vector<",
        overflow_threshold = "sole NT-total threshold -> test-only INT_MAX or2"),
    boundary_kind_names = c("allocVector", "allocMatrix", "R_alloc", "coerceVector",
        "INTEGER-materialization", "REAL-materialization", "mkChar", "setAttrib",
        "precious-preserve", "get-callable", "matrix-dim-attribute"),
    scope = "new count/mt native-copy scratch lifetimes; not universal Rcpp preservation/OOM leak freedom",
    dependencies = list(Rcpp = as.character(utils::packageVersion("Rcpp"))),
    session = utils::sessionInfo(), results = results,
    negative_control = list(error = negative_error, state = negative_state)),
    file.path(output, "count-ownership-results.rds"))
stopifnot(all(file.copy(paths, output, overwrite = TRUE)))
cat("Passed", length(results), "controlled count ownership cases.\n")
