# Score arrays share one ordered fit-cell view: y/mu/a/w are length n;
# D is p by n: rows are eigenbasis coordinates in eigen() order, with p equal
# to the design coefficient count; columns retain the ordered fit-cell view.
# Keep the legacy arithmetic/eigen decomposition order for numerical parity.
#' Form the existing score denominator factor without changing eigensystem order
#'
#' @param Zt_wZ Numeric symmetric p-by-p design crossproduct.
#' @param wZ Numeric n-by-p design multiplied rowwise by weights w, not sqrt(w).
#' @return Numeric p-by-n matrix `(Lambda^(-1/2) * t(U)) %*% t(wZ)` from the symmetric
#'   eigendecomposition. Rows are eigenbasis coordinates in `eigen()` order,
#'   with row count equal to the design coefficient count; columns retain the
#'   fit-cell order in `wZ`.
#' @details No ridge, eigenvalue clipping, or pseudoinverse is introduced; existing
#'   nonpositive-eigenvalue behavior propagates. Inputs, RNG and files are untouched.
#' @noRd
compute_D_matrix <- function(Zt_wZ, wZ) {
    P_decomp <- eigen(Zt_wZ, symmetric = TRUE)
    U <- P_decomp$vectors
    Lambda_minus_half <- 1 / sqrt(P_decomp$values)
    D <- (Lambda_minus_half * t(U)) %*% t(wZ)
    return(D)
}


#' Prepare the current score arrays on the nuisance fit's ordered cell view
#' @param expression_vector Numeric response vector of length n.
#' @param covariate_matrix Numeric n-by-p matrix aligned with response and fit.
#' @param fit Compact record with p fitted_coefs in design-column order and
#'   scalar theta, computed for this same fit-cell view.
#' @return List mu, w, a (length n) and D (p by n), preserving the existing full
#'   statistic formulas and numerical warnings. Does not refit, mutate input,
#'   persist fitted vectors, resample or access storage/RNG/files.
#' @noRd
prepare_score_components <- function(expression_vector, covariate_matrix, fit) {
    fitted_coefs <- fit$fitted_coefs
    theta <- fit$theta
    mu <- exp(as.numeric(covariate_matrix %*% fitted_coefs))
    denom <- 1 + mu / theta
    w <- mu / denom
    a <- (expression_vector - mu) / denom
    wZ <- w * covariate_matrix
    Zt_wZ <- t(covariate_matrix) %*% wZ
    D <- compute_D_matrix(Zt_wZ, wZ)
    list(mu = mu, w = w, a = a, D = D)
}
