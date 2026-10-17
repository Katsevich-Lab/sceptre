#!/usr/bin/env Rscript
# Separate bounded actual-source probe; no actual OOM or production hooks.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L)
    stop("Supply immutable candidate source, verified compiled original source, and fresh evidence directory.")
root <- normalizePath(args[[1L]], mustWork = TRUE)
reference <- normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
dir.create(args[[3L]], recursive = TRUE, showWarnings = FALSE)
output <- normalizePath(args[[3L]], mustWork = TRUE)
if (length(list.files(output, all.files = TRUE, no.. = TRUE)))
    stop("Ownership evidence directory must be fresh; previous attempts are never overwritten.")
dir.create(file.path(output, "transforms"))
dir.create(file.path(output, "originals"))
relative <- c("src/inference.h", "src/inference_scores.cpp", "src/inference_tail.cpp",
    "src/inference_stages.cpp", "src/negative_control_functions.cpp", "src/low_level_full_test.cpp",
    "src/result_lifetime.h", "src/sampling_safety.h", "src/pair_eligibility.h",
    "src/RcppExports.cpp", "tools/refactor/result-ownership.cpp",
    "tools/refactor/test-result-ownership.R",
    "tools/refactor/consistency/testthat/fixtures/refactor/native-inputs.rds",
    "tools/refactor/consistency/testthat/fixtures/refactor/native-expected.rds",
    "tools/refactor/consistency/testthat/helper-refactor.R", "tools/refactor/consistency/testthat/fixtures/refactor/manifest.rds")
paths <- file.path(root, relative)
before <- tools::md5sum(paths)
stopifnot(all(file.exists(paths)),
    identical(Sys.getenv("OPENBLAS_NUM_THREADS"), "1"),
    identical(Sys.getenv("OMP_NUM_THREADS"), "1"))
source_files <- paths[seq_len(6L)]
generated <- file.path(output, "transforms", basename(source_files))
for (i in seq_along(source_files)) {
    text <- readLines(source_files[[i]])
    text <- gsub("std::vector<", "result_vector<", text, fixed = TRUE)
    text <- gsub("std::unordered_set<", "result_set<", text, fixed = TRUE)
    # Test-only lifetime markers have copy/move-aware accounting, no production
    # mutation and no changed native/R allocation or owner scope.
    for (name in c("TailResult", "CalibrationPairResult")) {
        needle <- paste0("struct ", name, " {")
        if (any(grepl(needle, text, fixed = TRUE))) {
            stopifnot(sum(grepl(needle, text, fixed = TRUE)) == 1L)
            text <- sub(needle, paste0(needle, "\n  ResultPayloadLifetime probe_lifetime;"), text, fixed = TRUE)
        }
    }
    for (header in c("inference.h", "result_lifetime.h", "sampling_safety.h", "pair_eligibility.h")) {
        resolved <- if (header == "inference.h") generated[[1L]] else file.path(root, "src", header)
        text <- sub(paste0('#include "', header, '"'),
                    paste("#include", encodeString(resolved, quote = '"')), text, fixed = TRUE)
    }
    writeLines(text, generated[[i]])
}
macros <- c("INFERENCE", "SCORES", "TAIL", "STAGES", "CALIBRATION", "ENTRY")
bridge <- file.path(output, "result-ownership-instrumented.cpp")
writeLines(c(sprintf("#define SCEPTRE_RESULT_%s %s", macros,
    vapply(generated, encodeString, character(1), quote = '"')),
    readLines(file.path(root, "tools/refactor/result-ownership.cpp"))), bridge)
environment <- new.env(parent = globalenv())
Rcpp::sourceCpp(bridge, env = environment, rebuild = TRUE, showOutput = TRUE)
options(warn = 2) # warnings are failures during controlled qualification
# The CI oracle bootstrap already compiled this verified original archive.
# Load those bytes without building or looking for an installed/candidate fallback.
if ("sceptre" %in% loadedNamespaces()) stop("Use a fresh process before loading sceptre.")
source(file.path(root, "tools/refactor/consistency/testthat/helper-refactor.R"))
manifest <- readRDS(file.path(root, "tools/refactor/consistency/testthat/fixtures/refactor/manifest.rds"))
stopifnot(identical(manifest$source_commit, "21f9ea098b69c444a884a49bc254de711968e3b5"),
          isTRUE(manifest$reference_source_verified))
original_paths <- stats::setNames(file.path(reference, names(manifest$source_files)),
                                  names(manifest$source_files))
original_md5 <- refactor_oracle_hashes(original_paths)
verify_original <- function() {
    refactor_oracle_verify_source(reference, manifest$source_files, original_md5)
    actual <- vapply(original_paths, function(p) digest::digest(file = p, algo = "sha256"), "")
    stopifnot(identical(actual, unlist(manifest$source_files, use.names = TRUE)))
}
verify_original()
compiled_dll <- file.path(reference, "src", paste0("sceptre", .Platform$dynlib.ext))
if (!file.exists(compiled_dll)) stop("Verified original source must already have its compiled DLL.")
original_dll_md5 <- unname(tools::md5sum(compiled_dll))
pkgload::load_all(reference, compile = FALSE, helpers = FALSE, export_all = FALSE, quiet = TRUE)
namespace <- asNamespace("sceptre")
loaded_dll <- getNamespaceInfo("sceptre", "DLLs")[["sceptre"]]
stopifnot(identical(normalizePath(getNamespaceInfo("sceptre", "path"), winslash = "/"), reference),
          !is.null(loaded_dll),
          identical(unname(tools::md5sum(loaded_dll[["path"]])), original_dll_md5))
verify_original()
group_expected <- get("sample_combinations_v2", namespace)(2L, 2, 6, 4L, 4, 3L, .5)
pair_original <- get("sample_undercover_pairs_v2", namespace)
pair <- list(matrix(c(2L, 1L, 3L, 0L, 1L, 2L), 2L), c(4L, 5L, 4L),
             matrix(0:1, 2L), 3L, 4L, 1L, 1L)
rejected <- list(matrix(0L, 1L, 500L), rep(2L, 500L), matrix(0L, 100L, 1L),
                 500L, 50000L, 1L, 1L)
fixture <- readRDS(paths[[13L]])
expected <- readRDS(paths[[14L]])$cases
cases <- list(groups = list(kernel = 0L, input = list(), expected = group_expected),
    pairs_accepted = list(kernel = 1L, input = pair, expected = do.call(pair_original, pair)),
    pairs_all_rejected = list(kernel = 1L, input = rejected, expected = do.call(pair_original, rejected)))
tail_cases <- c("stage1_right", "accepted_sn_right", "accepted_sn_left",
                "failed_sn_stage3", "failed_sn_B3_zero", "empirical_B3_zero", "variable_crt_cardinality")
for (name in tail_cases)
    cases[[name]] <- list(kernel = 2L, input = fixture$cases[[name]], expected = expected[[name]])
unused <- fixture$cases$stage1_right
unused$B2 <- 0L
unused$B3 <- .Machine$integer.max - unused$B1
unused$draws <- unused$draws[seq_len(unused$B1)]
cases$stage1_huge_unused_B3 <- list(kernel = 2L, input = unused, expected = expected$stage1_right)
compact <- fixture$cases$stage1_right
compact$return_resampling_dist <- FALSE
cases$stage1_compact <- list(kernel = 2L, input = compact,
    expected = expected$stage1_right[setdiff(names(expected$stage1_right), "resampling_dist")])
results <- list()
record <- function(name, phase, position, error, state) {
    stopifnot(identical(state$live_blocks, 0L), identical(state$live_payloads, 0L),
        identical(state$protected_depth, 0L), identical(state$unsafe_boundaries, 0L),
        all((state$boundary_blocks == 0L & state$boundary_payloads == 0L) |
             state$boundary_protected == 1L))
    results[[length(results) + 1L]] <<- list(case = name, phase = phase,
        position = position, error = error, state = state)
}
baselines <- list()
set.seed(37911L)
rng <- .Random.seed
for (name in names(cases)) {
    case <- cases[[name]]
    bank <- if (case$kernel == 2L) environment$result_ownership_bank(case$input$draws) else NULL
    call <- function(boundary = 0L, native = 0L)
        environment$result_ownership_call(case$kernel, case$input, bank, boundary, native)
    # Warm Rcpp callable caches before enumerating stable boundary ordinals.
    call(); environment$result_ownership_state()
    snapshot <- serialize(case$input, NULL)
    actual <- call()
    state <- environment$result_ownership_state()
    stopifnot(isTRUE(all.equal(actual, case$expected, tolerance = 1e-8, check.attributes = TRUE)),
        identical(names(actual), names(case$expected)),
        identical(serialize(case$input, NULL), snapshot),
        state$native_count > 0L, state$max_live_blocks > 0L,
        identical(state$cleanup_calls, if (case$kernel == 0L) 0L else 1L))
    if (case$kernel < 2L) stopifnot(identical(actual, case$expected))
    if (case$kernel == 2L)
        stopifnot(identical(vapply(actual, typeof, ""), vapply(case$expected, typeof, "")))
    baselines[[name]] <- state
    record(name, "success", 0L, NA_character_, state)
    for (boundary in seq_len(state$boundary_count)) {
        error <- tryCatch({ call(boundary); NA_character_ }, error = conditionMessage)
        after <- environment$result_ownership_state()
        stopifnot(identical(error, "Controlled result R-boundary failure."),
                  identical(after$failed_boundary, boundary), after$cleanup_calls <= 1L)
        if (after$boundary_payloads[[boundary]] > 0L)
            stopifnot(identical(after$cleanup_calls, 1L))
        record(name, "raw-R-error", boundary, error, after)
    }
    for (native in seq_len(state$native_count)) {
        error <- tryCatch({ call(native = native); NA_character_ }, error = conditionMessage)
        after <- environment$result_ownership_state()
        stopifnot(!is.na(error), grepl("bad_alloc", error, fixed = TRUE),
                  identical(after$failed_native, native),
                  identical(after$cleanup_calls, 0L))
        record(name, "bad-alloc", native, error, after)
    }
    stopifnot(identical(serialize(case$input, NULL), snapshot))
}
stopifnot(identical(.Random.seed, rng),
    identical(baselines$stage1_right$peak_native_vector_bytes,
              baselines$stage1_huge_unused_B3$peak_native_vector_bytes),
    identical(baselines$pairs_all_rejected$peak_native_vector_bytes, 50000 * 4),
    baselines$pairs_all_rejected$max_live_payloads == 1L)
negative_error <- tryCatch({
    environment$result_ownership_negative_control(); NA_character_
}, error = conditionMessage)
negative <- environment$result_ownership_state()
stopifnot(identical(negative_error, "Unprotected R boundary with native result ownership."),
    identical(negative$unsafe_boundaries, 1L), identical(negative$live_blocks, 0L),
    identical(negative$live_payloads, 0L), negative$max_live_payloads > 0L)
stopifnot(identical(before, tools::md5sum(paths)))
verify_original()
stopifnot(identical(unname(tools::md5sum(compiled_dll)), original_dll_md5))
evidence <- list(source_files_md5 = before,
    source_files_md5_after = tools::md5sum(paths),
    original_source = list(path = reference, commit = manifest$source_commit,
        source_sha256 = manifest$source_files, source_md5 = original_md5,
        compiled_dll = compiled_dll, compiled_dll_md5 = original_dll_md5,
        loaded_dll_md5 = unname(tools::md5sum(loaded_dll[["path"]]))),
    generated_files_md5 = tools::md5sum(c(bridge, generated)),
    transformations = c("std::vector allocator", "unordered_set allocator",
        "copy/move-aware lifetime marker in TailResult and CalibrationPairResult",
        "absolute source header paths"),
    boundary_kinds = c("allocVector", "allocMatrix", "R_alloc", "coerceVector",
        "INTEGER", "REAL", "mkChar", "setAttrib", "precious-preserve", "get-callable",
        "make-unwind-continuation", "preserve-object", "protect"),
    scope = paste("native calibration scratch/results and selected tail payload only;",
        "not universal Rcpp wrapper/OOM safety; controlled errors, no actual OOM"),
    dependencies = list(Rcpp = as.character(utils::packageVersion("Rcpp")),
                         BH = as.character(utils::packageVersion("BH"))),
    threads = Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")),
    session = utils::sessionInfo(), baselines = baselines, results = results,
    negative_control = list(error = negative_error, state = negative))
saveRDS(evidence, file.path(output, "result-ownership-results.rds"))
for (i in seq_along(paths)) {
    destination <- file.path(output, "originals", relative[[i]])
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    stopifnot(file.copy(paths[[i]], destination, overwrite = FALSE))
}
stopifnot(identical(evidence$generated_files_md5, tools::md5sum(c(bridge, generated))))
cat("Passed", length(results), "controlled result ownership cases across", length(cases), "workloads.\n")
