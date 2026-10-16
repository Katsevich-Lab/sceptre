# Current response fitting is Poisson glm.fit plus a clipped dispersion estimate.
# It does not execute a second NB regression. Keep arithmetic/warnings here;
# this module has no object, storage, scheduling, or cache responsibilities.
#' Fit the existing response nuisance model on one aligned cell view
#'
#' @param expressions Numeric response vector of length n.
#' @param covariate_matrix Numeric n-by-p design, ordered with `expressions`;
#'   callers supply the validated design and select the fit-cell view.
#' @return Compact list `fitted_coefs` (p coefficients in design-column order) and
#'   scalar `theta`, clipped to [0.01, 1000] using the existing arithmetic.
#' @details Fits Poisson `glm.fit`, then estimates dispersion at its fitted means;
#'   no second negative-binomial regression is performed. Existing fit warnings
#'   propagate and missing coefficients produce the existing covariate error.
#'   Does not mutate inputs, cache results, access storage, or use RNG/files.
#' @noRd
fit_response_current <- function(expressions, covariate_matrix) {
    pois_fit <- stats::glm.fit(
        y = expressions,
        x = covariate_matrix,
        family = stats::poisson()
    )
    response_theta_list <- estimate_theta(
        y = expressions,
        mu = pois_fit$fitted.values,
        dfr = pois_fit$df.residual,
        limit = 50,
        eps = (.Machine$double.eps)^(1 / 4)
    )
    # check that NAs are absent from the fitted coefficient vector
    if (any(is.na(pois_fit$coefficients))) {
        problem_covariates <- paste0(names(which(is.na(pois_fit$coefficients))))
        stop(
            "The coefficients corresponding to the following covariates cannot be estimated in the regression model: ",
            paste0(problem_covariates, collapse = ", "),
            ". Consider removing these covariates from the model (by updating `formula_object` in `set_analysis_parameters()`)."
        )
    }
    theta <- max(min(response_theta_list[[1]], 1000), 0.01)
    result <- list(fitted_coefs = pois_fit$coefficients, theta = theta)
    return(result)
}

#' Fit the existing binary assignment model on one aligned cell view
#'
#' @param trt_idxs Valid R 1-based design-row positions of treated cells, not source
#'   columns. Repetitions set the same binary indicator entry to one.
#' @param covariate_matrix Numeric n-by-p design in the selected fit-cell order.
#' @param return_fitted_values Logical; return probabilities rather than coefficients.
#' @return Numeric vector of n fitted probabilities in fit-cell order, or p
#'   coefficients in design-column order. Uses binomial `glm.fit` with existing
#'   warnings and coefficient behavior; no additional convergence checks occur.
#' @details Does not mutate inputs, persist a cache, access storage, or use RNG/files.
#' @noRd
fit_assignment_current <- function(
    trt_idxs,
    covariate_matrix,
    return_fitted_values
) {
    indicator <- integer(length = nrow(covariate_matrix))
    indicator[trt_idxs] <- 1L
    logistic_fit <- stats::glm.fit(
        y = indicator,
        x = covariate_matrix,
        family = stats::binomial()
    )
    if (return_fitted_values) {
        out <- logistic_fit$fitted.values
    } else {
        out <- logistic_fit$coefficients
    }
    return(out)
}
