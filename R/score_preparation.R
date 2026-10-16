# Score arrays share one ordered fit-cell view: y/mu/a/w are length n;
# D is p by n, with design coefficients in rows and fit cells in columns.
# Keep the legacy arithmetic/eigen decomposition order for numerical parity.
#' Form the existing score denominator factor without changing eigensystem order
#'
#' @param Zt_wZ Numeric symmetric p-by-p design crossproduct.
#' @param wZ Numeric n-by-p design multiplied rowwise by weights w, not sqrt(w).
#' @return Numeric p-by-n matrix `(Lambda^(-1/2) * t(U)) %*% t(wZ)` from the symmetric
#'   eigendecomposition. Coefficient rows and fit-cell columns retain caller order.
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


#' Prepare ephemeral score arrays on the nuisance fit's cell view
#'
#' @param expression_vector Numeric response vector of length n.
#' @param covariate_matrix Numeric n-by-p design aligned with the response and fit.
#' @param fit Compact response record with p `fitted_coefs` in design-column order
#'   and scalar `theta`, computed for this same fit-cell view.
#' @param full_test_stat Logical; select the full or legacy reduced score arrays.
#' @return Full mode: list `mu`, `w`, `a` (length n), `D` (p by n). Reduced mode:
#'   list `mu`, `a`, `b` (length n), preserving the compatibility wrapper's formulas.
#' @details Keeps the existing arithmetic and numerical warnings. Does not refit,
#'   mutate inputs, persist fitted vectors, resample, or access storage/RNG/files.
#' @noRd
prepare_score_current <- function(
    expression_vector,
    covariate_matrix,
    fit,
    full_test_stat = TRUE
) {
    fitted_coefs <- fit$fitted_coefs
    theta <- fit$theta
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
