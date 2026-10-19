perturbation_test_design <- function(n = 120L) {
    cbind(intercept = 1, z = sin(seq_len(n)), w = cos(seq_len(n) / 3))
}

perturbation_test_bytes <- function(x) writeBin(as.double(x), raw())

test_that("compact perturbation models reproduce original glm probabilities exactly", {
    withr::local_seed(20261019)
    design <- perturbation_test_design(1000L)
    rownames(design) <- paste0("cell-", seq_len(nrow(design)))
    treated <- c(seq(3L, nrow(design), by = 5L), 3L)
    indicator <- integer(nrow(design))
    indicator[treated] <- 1L
    before <- .Random.seed
    reference <- stats::glm.fit(x = design, y = indicator, family = stats::binomial())
    model <- fit_perturbation_model(treated, design)
    expect_identical(model$probabilities, reference$fitted.values)
    expect_named(model$fit, c("fitted_coefs", "diagnostics"))
    expect_identical(model$fit$fitted_coefs, reference$coefficients)
    expect_named(model$fit$diagnostics,
        c("rank", "converged", "boundary", "iterations", "warnings"))
    expect_identical(perturbation_test_bytes(reconstruct_perturbation_probabilities(model$fit, design)),
        perturbation_test_bytes(reference$fitted.values))
    expect_identical(fit_assignment_probabilities(treated, design), reference$fitted.values)
    expect_identical(.Random.seed, before)
    expect_false(any(c("probabilities", "fitted.values", "qr", "linear.predictors", "cells") %in% names(model$fit)))
    expect_lt(as.numeric(object.size(model$fit)), as.numeric(object.size(model$probabilities)))
})

test_that("inverse-link reconstruction preserves glm clipping rather than plogis", {
    design <- matrix(1, 2, 1, dimnames = list(NULL, "intercept"))
    fit <- list(fitted_coefs = c(intercept = -31), diagnostics = list())
    actual <- reconstruct_perturbation_probabilities(fit, design)
    expect_identical(actual, stats::binomial()$linkinv(rep(-31, 2)))
    expect_false(identical(actual, stats::plogis(rep(-31, 2))))
    fit$fitted_coefs <- c(intercept = 31)
    expect_identical(reconstruct_perturbation_probabilities(fit, design),
        stats::binomial()$linkinv(rep(31, 2)))
})

test_that("aliased or inexact proposals decline caching without refitting", {
    design <- perturbation_test_design()
    design <- cbind(design, duplicate_z = design[, "z"])
    treated <- seq(3L, nrow(design), by = 5L)
    ordinary_fit <- fit_perturbation_glm
    calls <- 0L
    testthat::local_mocked_bindings(fit_perturbation_glm = function(...) {
        calls <<- calls + 1L
        ordinary_fit(...)
    })
    reference <- ordinary_fit(treated, design)
    expect_true(anyNA(reference$coefficients))
    model <- fit_perturbation_model(treated, design)
    expect_identical(calls, 1L)
    expect_null(model$fit)
    expect_identical(model$probabilities, reference$fitted.values)

    testthat::local_mocked_bindings(reconstruct_perturbation_probabilities = function(fit, design) {
        rep(0.37, nrow(design))
    })
    full_rank <- perturbation_test_design()
    reference <- ordinary_fit(treated, full_rank)
    model <- fit_perturbation_model(treated, full_rank)
    expect_identical(calls, 2L)
    expect_null(model$fit)
    expect_identical(model$probabilities, reference$fitted.values)
})

test_that("nonfinite coefficients decline compact proposals with original probabilities", {
    calls <- 0L
    testthat::local_mocked_bindings(fit_perturbation_glm = function(...) {
        calls <<- calls + 1L
        list(coefficients = c(intercept = Inf), fitted.values = c(0.2, 0.8))
    })
    model <- fit_perturbation_model(1L,
        matrix(1, 2, 1, dimnames = list(NULL, "intercept")))
    expect_identical(calls, 1L)
    expect_null(model$fit)
    expect_identical(model$probabilities, c(0.2, 0.8))
})

test_that("warnings are emitted once on fitting and stored without replay on reconstruction", {
    design <- perturbation_test_design()
    seen <- character()
    model <- withCallingHandlers(fit_perturbation_model(integer(), design), warning = function(w) {
        seen <<- c(seen, conditionMessage(w))
        invokeRestart("muffleWarning")
    })
    expect_length(seen, 1L)
    expect_match(seen, "did not converge")
    expect_false(model$fit$diagnostics$converged)
    expect_identical(model$fit$diagnostics$warnings, seen)
    expect_silent(actual <- reconstruct_perturbation_probabilities(model$fit, design))
    expect_identical(actual, model$probabilities)
    expect_identical(model$fit$diagnostics$iterations, 25L)
})

test_that("reconstruction rejects incompatible coefficient columns", {
    design <- perturbation_test_design()
    model <- fit_perturbation_model(seq(3L, nrow(design), by = 5L), design)
    expect_error(reconstruct_perturbation_probabilities(model$fit,
        design[, c("z", "intercept", "w")]), "ordered design columns")
    expect_error(reconstruct_perturbation_probabilities(model$fit,
        design[, c("intercept", "z")]), "aligned to its design")
})

test_that("CRT workers preserve exact native outputs with supplied compact probabilities", {
    withr::local_seed(20261019)
    store <- methods::as(Matrix::Matrix(rbind(r1 = rep(c(2, 3, 5), 4),
        r2 = c(1, 3, 5, 4, 6, 2, 7, 1, 4, 2, 5, 3)), sparse = TRUE), "RsparseMatrix")
    source <- c(2L, 4L, 5L, 8L, 11L)
    design <- matrix(1, length(source), 1, dimnames = list(NULL, "intercept"))
    treatment <- list(trt_idxs = 1:2, n_trt = 2L)
    model <- fit_perturbation_model(treatment$trt_idxs, design)
    probabilities <- reconstruct_perturbation_probabilities(model$fit, design)
    options <- list(B1 = 9L, B2 = 0L, B3 = 0L, total_resamples = 9L,
        side_code = 1L, fit_parametric_curve = FALSE)
    fits <- lapply(1:2, function(i) fit_response(as.numeric(store[i, source]), design))
    names(fits) <- rownames(store)
    args <- list(response_store = store, response_rows = 1:2,
        response_ids = rownames(store), design = design, source_cells = source,
        treatment = treatment, response_fits = fits, test_options = options,
        output_amount = 3L)
    before <- .Random.seed
    cold <- do.call(test_crt_shared_fit, args)
    supplied <- do.call(test_crt_shared_fit, c(args, list(probabilities = probabilities)))
    expect_identical(supplied, cold)

    view <- methods::new("CellView", fit_to_analysis = c(2L, 3L, 4L, 1L, 5L),
        source_cells = source[c(2L, 3L, 4L, 1L, 5L)], trt_fit_pos = 1:2, n_trt = 2L)
    args <- list(response_store = store, response_rows = 1:2, design = design,
        fit_view = view, test_options = options, output_amount = 3L)
    cold <- do.call(test_crt_group_fit, args)
    supplied <- do.call(test_crt_group_fit, c(args, list(probabilities = probabilities)))
    expect_identical(supplied, cold)
    expect_identical(.Random.seed, before)

    testthat::local_mocked_bindings(fit_assignment_probabilities = function(...) stop("unexpected refit"))
    expect_identical(do.call(test_crt_group_fit, c(args, list(probabilities = probabilities))), supplied)
})
