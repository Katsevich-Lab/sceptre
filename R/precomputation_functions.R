#' Run response precomputation
#'
#' Runs the precomputation for a given response vector, which consists of
#' regressing the response expressions onto the covariate matrix.
#'
#' We (i) perform an initial Poisson regression, (ii) estimate the size
#' parameter using the Poisson residuals, and then (iii) perform an NB
#' regression, treating the size parameter as estimated in step (ii) as fixed
#' and known.
#'
#' @param expressions the numeric vector of response expressions
#' @param covariate_matrix the covariate matrix on which to regress the
#' expressions (NOTE: the matrix should contain an interecept term)
#'
#' @return a list containing the following elements: (i) "fitted_coefs": a
#' vector of fitted coefficients; (ii) "theta": the fitted theta.
#' @noRd
perform_response_precomputation <- function(expressions, covariate_matrix) {
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


perform_grna_precomputation <- function(
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


#' Obtain the gRNA fitted probabilities, reusing a cached fit when available
#'
#' The logistic regression underlying the gRNA precomputation depends only on
#' which cells contain the gRNA target and on the covariate matrix; it does not
#' depend on the response counts. Analyses that share a set of cells and a set
#' of gRNA assignments can therefore share the fit, in the same way that they
#' already share the response precomputations.
#'
#' We cache the coefficients rather than the fitted values. The latter carry one
#' entry per cell, which for a large screen would make the cache far bigger than
#' the data it summarizes, whereas the coefficients are a handful of numbers per
#' gRNA target. Nothing is lost: [reconstruct_grna_fitted_probabilities()]
#' recovers the fitted values from the coefficients exactly.
#'
#' @param grna_precomputations the named list of cached precomputations
#' @param curr_grna_group the gRNA target being analyzed
#' @param trt_idxs the indices of the treated cells, into the rows of
#' `covariate_matrix`
#' @param covariate_matrix the covariate matrix to regress on
#' @param fit_idxs the index vector that identifies this fit; see
#' [grna_precomputation_is_current()]
#'
#' @return a list containing the following elements: (i) "fitted_probabilities":
#' the fitted probabilities, one per row of `covariate_matrix`; (ii) "to_cache":
#' the precomputation to add to the cache, or `NULL` if the cache already holds
#' it or the fit cannot be cached.
#' @noRd
get_grna_fitted_probabilities <- function(
    grna_precomputations,
    curr_grna_group,
    trt_idxs,
    covariate_matrix,
    fit_idxs
) {
    # `NULL[["x"]]` is NULL, so an absent cache falls through to the fit
    grna_precomp <- grna_precomputations[[curr_grna_group]]
    if (grna_precomputation_is_current(grna_precomp, covariate_matrix, fit_idxs)) {
        return(list(
            fitted_probabilities = reconstruct_grna_fitted_probabilities(
                fitted_coefs = grna_precomp$fitted_coefs,
                covariate_matrix = covariate_matrix
            ),
            to_cache = NULL
        ))
    }
    fitted_coefs <- perform_grna_precomputation(
        trt_idxs = trt_idxs,
        covariate_matrix = covariate_matrix,
        return_fitted_values = FALSE
    )
    # a rank-deficient fit yields NA coefficients for the aliased columns while
    # still producing valid fitted values, so it cannot be replayed from its
    # coefficients; fall back to the fitted values and decline to cache
    if (anyNA(fitted_coefs)) {
        return(list(
            fitted_probabilities = perform_grna_precomputation(
                trt_idxs = trt_idxs,
                covariate_matrix = covariate_matrix,
                return_fitted_values = TRUE
            ),
            to_cache = NULL
        ))
    }
    list(
        fitted_probabilities = reconstruct_grna_fitted_probabilities(
            fitted_coefs = fitted_coefs,
            covariate_matrix = covariate_matrix
        ),
        to_cache = list(
            fitted_coefs = fitted_coefs,
            fit_idxs = as.integer(fit_idxs),
            covariate_names = colnames(covariate_matrix)
        )
    )
}


#' Recover the gRNA fitted probabilities from the fitted coefficients
#'
#' `stats::glm.fit()` computes the fitted values it returns as
#' `linkinv(drop(x %*% coefficients))` with a zero offset, and then strips the
#' names of the result. Reproducing both steps here makes the reconstruction
#' bit-identical to a fresh fit; without the `unname()` the values still agree
#' but `identical()` does not, because `drop(x %*% coefs)` inherits the row names
#' of the covariate matrix. Exactness matters: a one-ulp difference here would
#' flip a Bernoulli draw in `crt_index_sampler_fast()` and change a p-value.
#'
#' @param fitted_coefs the fitted coefficients of the logistic regression
#' @param covariate_matrix the covariate matrix the fit was performed on
#'
#' @return a numeric vector of fitted probabilities, one per cell
#' @noRd
reconstruct_grna_fitted_probabilities <- function(
    fitted_coefs,
    covariate_matrix
) {
    linear_predictor <- drop(covariate_matrix %*% fitted_coefs)
    unname(stats::binomial()$linkinv(linear_predictor))
}


#' Determine whether a cached gRNA precomputation applies to the current fit
#'
#' A cached fit is valid only for the design it was produced on, and the two CRT
#' workhorse functions use different row spaces: `crt_glm_factored_out()`
#' regresses on every cell in use, whereas `discovery_ntcells_crt()` regresses on
#' the treated cells followed by the non-targeting cells. An entry cached by one
#' is silently wrong in the other, and the two analyses share a cache, so a
#' length check alone would not separate them.
#'
#' We therefore fingerprint each entry with `fit_idxs`, the index vector that
#' determines the fit: the treated cells for `crt_glm_factored_out()`, and the
#' treated-plus-non-targeting cells for `discovery_ntcells_crt()`. Comparing it
#' identifies the code path, the gRNA assignment, and the cell set in a single
#' check. A stale entry is refit and replaced rather than raising an error, since
#' the cache is maintained by the package rather than supplied by the user.
#'
#' @param grna_precomp one element of a `grna_precomputations` list, or `NULL`
#' @param covariate_matrix the covariate matrix for the current fit
#' @param fit_idxs the index vector identifying the current fit
#'
#' @return `TRUE` if the cached entry can be replayed, and `FALSE` otherwise
#' @noRd
grna_precomputation_is_current <- function(
    grna_precomp,
    covariate_matrix,
    fit_idxs
) {
    !is.null(grna_precomp) &&
        identical(
            grna_precomp$covariate_names,
            colnames(covariate_matrix)
        ) &&
        identical(grna_precomp$fit_idxs, as.integer(fit_idxs))
}


#' Read the gRNA precomputations from a `sceptre_object`
#'
#' The `grna_precomputations` slot was added to the class after the first
#' release, so a `sceptre_object` deserialized from an earlier version does not
#' carry it. Reading through this accessor lets such an object be analyzed as
#' usual; the slot is created the first time it is assigned.
#'
#' @param sceptre_object a `sceptre_object`
#'
#' @return the named list of cached gRNA precomputations
#' @noRd
get_grna_precomputations <- function(sceptre_object) {
    if (methods::.hasSlot(sceptre_object, "grna_precomputations")) {
        sceptre_object@grna_precomputations
    } else {
        list()
    }
}


compute_D_matrix <- function(Zt_wZ, wZ) {
    P_decomp <- eigen(Zt_wZ, symmetric = TRUE)
    U <- P_decomp$vectors
    Lambda_minus_half <- 1 / sqrt(P_decomp$values)
    D <- (Lambda_minus_half * t(U)) %*% t(wZ)
    return(D)
}


compute_precomputation_pieces <- function(
    expression_vector,
    covariate_matrix,
    fitted_coefs,
    theta,
    full_test_stat
) {
    mu <- exp(as.numeric(covariate_matrix %*% fitted_coefs))
    if (full_test_stat) {
        denom <- 1 + mu / theta
        w <- mu / denom
        a <- (expression_vector - mu) / denom
        wZ <- w * covariate_matrix
        Zt_wZ <- t(covariate_matrix) %*% wZ
        D <- compute_D_matrix(Zt_wZ, wZ)
        out <- list(mu = mu, w = w, a = a, D = D)
    } else {
        a <- expression_vector -
            (expression_vector * mu + theta * mu) / (theta + mu)
        b <- (theta * mu) / (theta + mu)
        out <- list(mu = mu, a = a, b = b)
    }
    return(out)
}
