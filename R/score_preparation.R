# Score arrays share one ordered fit-cell view: y/mu/a/w are length n;
# D is p by n, with design coefficients in rows and fit cells in columns.
# Keep the legacy arithmetic/eigen decomposition order for numerical parity.
compute_D_matrix <- function(Zt_wZ, wZ) {
    P_decomp <- eigen(Zt_wZ, symmetric = TRUE)
    U <- P_decomp$vectors
    Lambda_minus_half <- 1 / sqrt(P_decomp$values)
    D <- (Lambda_minus_half * t(U)) %*% t(wZ)
    return(D)
}


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
