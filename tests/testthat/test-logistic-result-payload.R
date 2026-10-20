# This fixture has an exact finite solution, so packaging is checked against
# likelihood mathematics rather than a serialization of the implementation.
logistic_expected_payload <- function() {
    list(
        coefficients = 0,
        fitted.values = rep(0.5, 8L),
        score = 0,
        converged = TRUE,
        iterations = 0L,
        reason = "ok",
        loglik = -8 * log(2),
        deviance = 16 * log(2),
        gradient_max = 0,
        score_tolerance = 4e-8,
        conditioned_gradient_max = 0,
        conditioned_score_tolerance = 4e-8,
        relative_deviance_change = NaN,
        max_abs_eta = 0,
        rank = 1L,
        information_condition_ratio = NaN,
        full_cell_scans = 0L,
        step_halvings = 0L,
        conditioned = TRUE,
        separation_detected = FALSE,
        used_analytic_prevalence_state = TRUE
    )
}

expect_logistic_payload <- function(actual, expected) {
    # In particular, the appended 21st field must remain a named logical scalar.
    expect_identical(names(actual), names(expected))
    expect_identical(
        vapply(actual, typeof, character(1)),
        vapply(expected, typeof, character(1))
    )
    expect_identical(lengths(actual), lengths(expected))
    expect_true(all(vapply(
        actual,
        function(value) is.null(attributes(value)),
        logical(1)
    )))
    expect_equal(actual, expected, tolerance = 1e-14)
}

test_that("native logistic success preserves all 21 result fields", {
    design <- matrix(1, 8L, 1L, dimnames = list(NULL, "(Intercept)"))
    prepared <- prepare_grna_logistic_design(design)
    actual <- fit_grna_logistic_prepared_cpp(prepared, c(1L, 3L, 5L, 7L))
    expect_logistic_payload(actual, logistic_expected_payload())
})

test_that("native logistic numerical failure preserves all 21 result fields", {
    design <- matrix(1, 8L, 1L, dimnames = list(NULL, "(Intercept)"))
    prepared <- prepare_grna_logistic_design(design)
    actual <- fit_grna_logistic_prepared_cpp(prepared, integer())
    expected <- logistic_expected_payload()
    expected$fitted.values <- numeric()
    expected$score <- NaN
    expected$converged <- FALSE
    expected$reason <- "degenerate_binary_response"
    for (field in c(
        "loglik",
        "deviance",
        "gradient_max",
        "score_tolerance",
        "conditioned_gradient_max",
        "conditioned_score_tolerance",
        "max_abs_eta"
    )) {
        expected[[field]] <- NaN
    }
    expected$separation_detected <- TRUE
    expect_logistic_payload(actual, expected)
})
