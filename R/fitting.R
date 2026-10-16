# Current response fitting is Poisson glm.fit plus a clipped dispersion estimate.
# It does not execute a second NB regression. Keep arithmetic/warnings here;
# this module has no object, storage, scheduling, or cache responsibilities.
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
