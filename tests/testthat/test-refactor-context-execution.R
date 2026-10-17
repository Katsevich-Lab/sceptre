test_that("NT intercept-only calibration and discovery preserve explicit matrix views", {
    counts <- matrix(c(1, 2, 5, 0, 7, 1, 2, 4, 1, 2, 0, 1), 1L,
                     dimnames = list("r", NULL))
    object <- methods::new("sceptre_object")
    store <- methods::as(Matrix::Matrix(counts, sparse = TRUE), "RsparseMatrix")
    object <- refactor_private("set_response_matrix")(object, store)
    object@cells_in_use <- c(2L, 4L, 5L, 8L, 11L)
    object@covariate_matrix <- matrix(1, 12L, 1L,
                                      dimnames = list(NULL, "(Intercept)"))
    object@low_moi <- TRUE
    object@control_group_complement <- FALSE
    object@grna_assignments <- list(grna_group_idxs = list(A = c(2L, 3L)),
        indiv_nt_grna_idxs = list(ntA = 1L, ntB = c(2L, 3L)),
        all_nt_idxs = c(4L, 1L, 5L))
    object@B1 <- 9L
    object@B2 <- object@B3 <- 0L
    object@side_code <- 1L
    object@calibration_group_size <- 1L
    object@resampling_approximation <- "no_approximation"

    # Write the original matrix columns explicitly, without using the selection
    # helpers to calculate expected values. Calibration uses NT cells (8,2,11);
    # discovery places treatment cells before the NT cells.
    for (calibration in c(TRUE, FALSE)) {
        source <- if (calibration) c(8L, 2L, 11L) else c(4L, 5L, 8L, 2L, 11L)
        treatment <- if (calibration) c(2L, 3L) else 1:2
        role <- if (calibration) "calibration_check" else "discovery_analysis"
        pairs <- data.frame(response_id = "r", grna_group = if (calibration) "ntB" else "A",
                            pass_qc = TRUE)
        y <- as.numeric(counts[1L, source])
        mu <- mean(y)
        # The response dispersion estimator is unchanged; the independent
        # intercept score below does not call the new fitter/score/view helpers.
        theta <- max(min(refactor_private("estimate_theta")(
            y, rep(mu, length(y)), length(y) - 1L, limit = 50,
            eps = .Machine$double.eps^(1 / 4))[[1L]], 1000), 0.01)
        denom <- 1 + mu / theta
        k <- length(treatment)
        expected_z <- sum((y[treatment] - mu) / denom) /
            sqrt(mu / denom * (k - k^2 / length(y)))
        for (permutations in c(TRUE, FALSE)) {
            object@run_permutations <- permutations
            before <- serialize(object, NULL)
            context <- refactor_private("prepare_analysis_context")(object, role, pairs)
            expect_identical(dim(context$design), c(length(source), 1L))
            expect_identical(context$common_view$source_cells,
                if (calibration) source else object@cells_in_use)
            actual <- suppressMessages(refactor_private("analyze_association_pairs")(
                object, pairs, role, output_amount = 3L, print_progress = FALSE,
                parallel = FALSE, n_processors = 1L, log_dir = tempdir()))
            expect_equal(actual$result$z_orig, expected_z, tolerance = 1e-8)
            expect_identical(actual$result$response_id, "r")
            expect_identical(actual$result$grna_group, pairs$grna_group)
            expect_identical(dim(object@covariate_matrix), c(12L, 1L))
            expect_identical(serialize(object, NULL), before)
            expect_length(actual$response_precomputations, if (calibration) 1L else 0L)
        }
    }
})

test_that("prepared integer partitions retain label and caller RNG order", {
    partition <- refactor_private("partition_response_ids")
    labels <- factor(c("z", "a", "b", "x", "q", "c"),
                     levels = c("unused", "c", "z", "a", "b", "x", "q"))
    rows <- c(8L, 3L, 7L, 1L, 9L, 4L)
    for (parallel in c(FALSE, TRUE)) {
        set.seed(182)
        before <- .Random.seed
        by_label <- partition(labels, parallel, 2L)
        by_row <- partition(rows, parallel, 2L)
        positions <- match(unlist(by_row, use.names = FALSE), rows)
        expect_identical(as.character(labels[positions]),
                         as.character(unlist(by_label, use.names = FALSE)))
        expect_identical(.Random.seed, before)
    }
})

test_that("progress bytes and periodic maintenance retain independent cadences", {
    report <- refactor_private("report_worker_progress")
    positions <- c(1L, 4L, 5L, 199L, 200L, 201L)
    for (feature in c("response", "gRNA group", "gRNA")) {
        prefix <- if (feature == "gRNA") "Performing gRNA-to-cell assignments for " else
            "Analyzing pairs containing "
        messages <- character()
        withCallingHandlers(for (i in positions) report(i, "id", 201L, TRUE, FALSE,
            NULL, feature, prefix), message = function(cnd) {
                messages <<- c(messages, conditionMessage(cnd))
                invokeRestart("muffleMessage")
            })
        expected <- paste0(prefix, feature, " id (", c(1L, 5L, 200L), " of 201)\n")
        expect_identical(messages, expected)
        file <- tempfile()
        on.exit(unlink(file), add = TRUE)
        file.create(file)
        for (i in positions) report(i, "id", 201L, TRUE, TRUE, file, feature, prefix)
        bytes <- readBin(file, "raw", n = file.info(file)$size)
        expected_text <- paste0(paste0(expected, "\n"), collapse = "")
        # Existing write(..., file = path) uses native text-file endings. Keep
        # comparing raw bytes, including its extra newline, on both platforms.
        if (.Platform$OS.type == "windows") {
            expected_text <- gsub("\n", "\r\n", expected_text, fixed = TRUE)
        }
        expect_identical(bytes, charToRaw(expected_text))
        absent <- tempfile()
        for (i in positions) report(i, "id", 201L, FALSE, TRUE, absent, feature, prefix)
        expect_false(file.exists(absent))
    }
    # Track gc() calls in this local copy of the helper, leaving base and package
    # functions unchanged. Garbage collection must still run with progress hidden.
    collect <- refactor_private("collect_worker_garbage")
    calls <- 0L
    scope <- new.env(parent = environment(collect))
    scope$gc <- function() { calls <<- calls + 1L; NULL }
    environment(collect) <- scope
    for (i in c(1L, 199L, 200L, 201L, 400L)) collect(i)
    expect_identical(calls, 2L)
})

test_that("public workflows retain fit, bank and row-read sharing scopes", {
    inputs <- readRDS(file.path(refactor_fixture_dir(), "inputs.rds"))
    base_case <- readRDS(file.path(refactor_fixture_dir(), "cases.rds"))$low_nt_perm
    originals <- lapply(c("fit_response", "fit_assignment_probabilities",
        "sample_crt_indices", "get_synthetic_permutation_idxs",
        "storage_read_response"), refactor_private)
    names(originals) <- c("fit", "assignment", "crt", "permutation", "read")
    events <- new.env(parent = emptyenv())
    reset <- function() {
        events$fit <- events$assignment <- events$read <- 0L
        events$crt <- events$permutation <- numeric()
    }
    testthat::local_mocked_bindings(
        fit_response = function(...) {
            events$fit <- events$fit + 1L; originals$fit(...)
        },
        fit_assignment_probabilities = function(...) {
            events$assignment <- events$assignment + 1L; originals$assignment(...)
        },
        sample_crt_indices = function(fitted_probabilities, B) {
            events$crt <- c(events$crt, B); originals$crt(fitted_probabilities, B)
        },
        get_synthetic_permutation_idxs = function(..., B) {
            events$permutation <- c(events$permutation, B); originals$permutation(..., B = B)
        },
        storage_read_response = function(...) {
            events$read <- events$read + 1L; originals$read(...)
        }, .package = "sceptre")
    for (control in c("nt_cells", "complement")) {
        for (mechanism in c("permutations", "crt")) {
            case <- base_case
            case$control <- control
            case$mechanism <- mechanism
            object <- suppressMessages(refactor_prepare(inputs, case))
            pairs <- object@discovery_pairs_with_info
            pairs <- pairs[pairs$pass_qc, , drop = FALSE]
            nr <- length(unique(pairs$response_id))
            ng <- length(unique(pairs$grna_group))
            np <- nrow(pairs)
            for (warm in c(FALSE, TRUE)) {
                reset()
                before_seed <- withr::with_seed(817L, {
                    seed <- .Random.seed
                    result <- suppressMessages(sceptre::run_discovery_analysis(object,
                        output_amount = 2L, print_progress = FALSE))
                    expect_identical(.Random.seed, seed)
                    result
                })
                expect_identical(events$fit,
                    if (control == "nt_cells") np else if (warm) 0L else nr)
                if (mechanism == "permutations") {
                    expect_identical(events$read, nr)
                    expect_identical(events$assignment, 0L)
                    expect_equal(events$permutation, sum(as.double(c(object@B1, object@B2, object@B3))))
                    expect_length(events$crt, 0L)
                } else {
                    expect_identical(events$read, np + if (control == "complement") nr else 0L)
                    expect_identical(events$assignment, ng)
                    expect_length(events$permutation, 0L)
                    expect_equal(events$crt, rep(sum(as.double(c(object@B1, object@B2, object@B3))), ng))
                }
                if (warm) expect_equal(before_seed@discovery_result,
                                       object@discovery_result, tolerance = 1e-8)
                object <- before_seed
            }
        }
    }
})

test_that("mixture keeps actual EM success, raw-count fallback and caller RNG", {
    g <- c(rep(c(0, 1, 0, 2), 25L), rep(c(25, 35, 30, 40), 2L))
    design <- matrix(1, length(g), 1L)
    guesses <- refactor_private("get_random_starting_guesses")()
    fit_em <- refactor_private("run_reduced_em_algo_cpp")
    pois <- suppressWarnings(stats::glm.fit(design, g, family = stats::poisson()))
    fit <- fit_em(guesses$pi_guesses, guesses$g_pert_guesses, g,
                  pois$fitted.values, lgamma(g + 1))
    expect_true(fit$outer_converged && fit$outer_log_lik != -Inf)
    assign <- refactor_private("obtain_em_assignments")
    args <- list(pi_guesses = guesses$pi_guesses, g_pert_guesses = guesses$g_pert_guesses,
        g = g, covariate_matrix = design, n_nonzero_cells_cutoff = 1L,
        backup_threshold = 20, probability_threshold = 0.8)
    before <- serialize(args, NULL)
    set.seed(271L)
    seed <- .Random.seed
    expect_identical(do.call(assign, args), which(fit$outer_Ti1s >= 0.8))
    expect_identical(.Random.seed, seed)
    expect_identical(serialize(args, NULL), before)
    calls <- 0L
    testthat::local_mocked_bindings(run_reduced_em_algo_cpp = function(...) {
        calls <<- calls + 1L
        fit_em(...)
    }, .package = "sceptre")
    # Fit this case directly in C++ to show that no EM solution is accepted.
    # Assignment must therefore use the existing raw-count threshold fallback.
    failed_g <- c(rep(c(0, 1, 0, 2), 10L), rep(c(25, 35, 30, 40), 10L))
    failed_design <- matrix(1, length(failed_g), 1L)
    failed_pois <- suppressWarnings(stats::glm.fit(failed_design, failed_g, family = stats::poisson()))
    failed_fit <- fit_em(guesses$pi_guesses, guesses$g_pert_guesses, failed_g,
                        failed_pois$fitted.values, lgamma(failed_g + 1))
    expect_false(failed_fit$outer_converged)
    args$g <- failed_g
    args$covariate_matrix <- failed_design
    expect_identical(do.call(assign, args), which(failed_g >= 20))
    expect_identical(calls, 1L)
    args$n_nonzero_cells_cutoff <- length(failed_g) + 1L
    expect_identical(do.call(assign, args), which(failed_g >= 20))
    expect_identical(calls, 1L)
    expect_identical(.Random.seed, seed)
})


test_that("mixture block traversal preserves progress files, assignments and RNG", {
    ids <- paste0("g", seq_len(402L))
    counts <- matrix(rep(c(0, 25), each = length(ids)), nrow = length(ids),
                     dimnames = list(ids, NULL))
    covariates <- data.frame(z = c(1, 2))
    hyperparameters <- list(n_em_rep = 2L, pi_guess_range = c(1e-5, 0.1),
        g_pert_guess_range = c(10, 50), formula_object = ~ 1,
        n_nonzero_cells_cutoff = 3L, backup_threshold = 20,
        probability_threshold = 0.8)
    run <- refactor_private("assign_grnas_to_cells_mixture")
    log_root <- tempfile()
    on.exit(unlink(log_root, recursive = TRUE), add = TRUE)
    set.seed(1829)
    seed <- .Random.seed
    serial <- suppressMessages(run(counts, covariates, hyperparameters,
        FALSE, FALSE, 1L, log_root, ids))
    expect_identical(names(serial), ids)
    expect_true(all(vapply(serial, identical, logical(1), 2L)))
    expect_identical(.Random.seed, seed)
    expect_false(dir.exists(log_root))
    skip_on_os("windows")
    parallel <- suppressMessages(run(counts, covariates, hyperparameters,
        TRUE, TRUE, 2L, log_root, ids))
    expect_identical(parallel[ids], serial)
    expected_ids <- withr::with_seed(4, sample(ids))
    expect_identical(names(parallel), expected_ids)
    expect_identical(.Random.seed, seed)
    for (block in 1:2) {
        block_ids <- expected_ids[(seq_len(201L)) + (block - 1L) * 201L]
        positions <- c(1L, seq.int(5L, 200L, by = 5L))
        expected <- paste0("Performing gRNA-to-cell assignments for gRNA ",
            block_ids[positions], " (", positions, " of 201)\n\n", collapse = "")
        path <- file.path(log_root, "sceptre_logs", paste0("assign_grnas_", block, ".out"))
        expect_identical(readBin(path, "raw", n = file.info(path)$size), charToRaw(expected))
    }
})
