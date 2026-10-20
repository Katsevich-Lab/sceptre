# Fit the response using Poisson glm.fit(), then estimate dispersion at those
# fitted means. There is no second negative-binomial regression. These functions
# take vectors and design matrices directly; their callers select cells, read
# counts and decide whether a previous fit can be reused.
#' Fit the response model for the selected cells
#' @param expressions Numeric response counts, one per selected cell.
#' @param covariate_matrix Numeric design matrix with one row per selected cell,
#'   in exactly the same order as `expressions`.
#' @return Coefficients (`fitted_coefs`) in design-column order and dispersion
#'   (`theta`), restricted to `[0.01, 1000]`. These are small enough to save for
#'   reuse.
#' @details Fit Poisson `glm.fit()`, then estimate dispersion at its fitted
#'   means.
#'   Warnings from the fit are passed through. Missing coefficients trigger the
#'   covariate error below. The caller handles saving or reusing this fit.
#' @noRd
fit_response <- function(expressions, covariate_matrix) {
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

#' Fit each selected cell's probability of belonging to the treatment group
#' @param trt_idxs Treated cells' row numbers in `covariate_matrix`, starting at
#'   1.
#'   These are positions in this fit, not columns of the original count matrix.
#'   Repeated positions mark the same cell as treated once.
#' @param covariate_matrix Numeric design matrix for the cells in this fit.
#' @return One fitted probability per design row, from binomial `glm.fit()`.
#'   Existing glm warnings and convergence behavior are passed through; this
#'   helper adds no separate convergence check.
#' @noRd
fit_assignment_probabilities <- function(trt_idxs, covariate_matrix) {
    fit_perturbation_glm(trt_idxs, covariate_matrix)$fitted.values
}

# Fit exactly once with the established indicator encoding and glm defaults.
# The caller may retain compact coefficients; the full glm object is transient.
#' @noRd
fit_perturbation_glm <- function(trt_idxs, design) {
    indicator <- integer(length = nrow(design))
    indicator[trt_idxs] <- 1L
    stats::glm.fit(y = indicator, x = design, family = stats::binomial())
}
