# Every vector below uses the same selected cells in the same order.
# y, mu, a and w have one value per cell; D has one column per cell.
# D's rows describe eigenvectors, in eigen() order, rather than individual
# covariates. There are as many rows as design coefficients.
# Keep the calculation and eigenvector order to preserve numerical results.
#' Build the matrix used to adjust the score variance for fitted covariates
#' @param Zt_wZ Symmetric p-by-p matrix `t(Z) %*% (w * Z)`, where Z is the design.
#' @param wZ n-by-p design matrix with each cell's row multiplied by its weight w
#'   (not by sqrt(w)).
#' @return A p-by-n matrix, `(Lambda^(-1/2) * t(U)) %*% t(wZ)`, where U and Lambda
#'   come from `eigen(Zt_wZ)`. Columns follow the cells in `wZ`; rows follow the
#'   eigenvectors in eigen() order, rather than individual covariates.
#' @details Use the eigenvalues directly, without adding a ridge term, clipping
#'   them or using a pseudoinverse. Nonpositive eigenvalues keep their existing
#'   numerical behavior.
#' @noRd
compute_D_matrix <- function(Zt_wZ, wZ) {
    P_decomp <- eigen(Zt_wZ, symmetric = TRUE)
    U <- P_decomp$vectors
    Lambda_minus_half <- 1 / sqrt(P_decomp$values)
    D <- (Lambda_minus_half * t(U)) %*% t(wZ)
    return(D)
}


#' Turn a fitted response model into the arrays used by the score test
#' @param expression_vector Response counts, one per selected cell.
#' @param covariate_matrix Design matrix whose rows follow those same cells.
#' @param fit Coefficients (`fitted_coefs`) and dispersion (`theta`) fitted on
#'   exactly these cells and this design.
#' @return Fitted means `mu`, weights `w`, residual contributions `a` (one value
#'   per cell) and the covariate-adjustment matrix `D` (one column per cell).
#'   The formulas below keep the current test's calculations and warnings.
#'   This helper prepares arrays from an existing fit; it does not fit again.
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
