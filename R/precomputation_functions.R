#' Run response precomputation
#'
#' Runs the precomputation for a given response vector, which consists of regressing the response expressions onto the covariate matrix.
#'
#' We (i) perform an initial Poisson regression, (ii) estimate the size parameter using the Poisson residuals, and then (iii) perform an NB regression, treating the size parameter as estimated in step (ii) as fixed and known.
#'
#' @param expressions the numeric vector of response expressions
#' @param covariate_matrix the covariate matrix on which to regress the expressions (NOTE: the matrix should contain an interecept term)
#'
#' @return a list containing the following elements: (i) "precomp_str": a string summarizing the method used to fit the GLM and calculate the NB size parameter. The part of the string before the colon is either "nb" or "pois," indicating whether NB regression or Poisson regression was used. ("pois_{warn}" indicates that Poisson regression was used and that a warning was generated). The part of the string after the colon indicates the method used to estimate the size parameter ("mass" for the "mass" package, "resid_mle" for MLE on the residuals, and "resid_mm" for method of moments on the residuals.) (ii) "fitted_coefs": a vector of fitted coefficients; the final entry of this vector is the fitted theta.
#' @noRd
perform_response_precomputation <- function(expressions, covariate_matrix, regression_method) {
  # backup: return fitted coefficients from Poisson regression
  backup_2 <- function(pois_fit, pois_warn) {
    list(fitted_coef_str = paste0("pois", if (pois_warn) "_(warn)" else NULL),
         fitted_coefs = pois_fit$coefficients)
  }

  # fit Poisson model, tracking warning
  pois_warn <- FALSE
  wHandler <- function(w) {pois_warn <<- TRUE; invokeRestart("muffleWarning")}
  withCallingHandlers(expr = {
    pois_fit <- stats::glm.fit(y = expressions, x = covariate_matrix, family = stats::poisson())
  }, warning = wHandler)

  # get theta; save theta itself (in response_theta) and theta_fit_string (i.e., procedure used to fit theta)
  response_theta_list <- estimate_theta(y = expressions, mu = pois_fit$fitted.values, dfr = pois_fit$df.residual,
                                        limit = 50, eps = (.Machine$double.eps)^(1/4))
  theta_fit_str <- c("mle", "mm", "pilot")[response_theta_list[[2]]]
  theta <- max(min(response_theta_list[[1]], 1000), 0.01)

  # obtain the fitted coefficients
  if (regression_method == "nb_glm") {
    fitted_coefs_list <- tryCatch({
      fit_nb <- stats::glm.fit(y = expressions,
                               x = covariate_matrix,
                               family = MASS::negative.binomial(theta),
                               mustart = pois_fit$fitted.values)
      list(fitted_coef_str = "nb",
           fitted_coefs = fit_nb$coefficients)
    }, error = function(e) backup_2(pois_fit, pois_warn), warning = function(w) backup_2(pois_fit, pois_warn))
  } else {
    fitted_coefs_list <- backup_2(pois_fit, pois_warn)
  }

  # set the precomp str
  precomp_str <- paste0(fitted_coefs_list$fitted_coef_str, ":", theta_fit_str)

  # output the result
  result <- list(precomp_str = precomp_str,
                 fitted_coefs = fitted_coefs_list$fitted_coefs,
                 theta = theta)
  return(result)
}


perform_grna_precomputation <- function(trt_idxs, covariate_matrix, return_fitted_values) {
  indicator <- integer(length = nrow(covariate_matrix))
  indicator[trt_idxs] <- 1L
  logistic_fit <- stats::glm.fit(y = indicator, x = covariate_matrix, family = stats::binomial())
  if (return_fitted_values) {
    out <- logistic_fit$fitted.values
  } else {
    out <- logistic_fit$coefficients
  }
  return(out)
}


compute_D_matrix <- function(Zt_wZ, wZ) {
  P_decomp <- eigen(Zt_wZ, symmetric = TRUE)
  U <- P_decomp$vectors
  Lambda_minus_half <- 1/sqrt(P_decomp$values)
  D <- (Lambda_minus_half * t(U)) %*% t(wZ)
  return(D)
}


compute_precomputation_pieces <- function(expression_vector, covariate_matrix, fitted_coefs, theta, full_test_stat) {
  mu <- exp(as.numeric(covariate_matrix %*% fitted_coefs))
  if (full_test_stat) {
    denom <- 1 + mu/theta
    w <- mu/denom
    a <- (expression_vector - mu)/denom
    wZ <- w * covariate_matrix
    Zt_wZ <- t(covariate_matrix) %*% wZ
    D <- compute_D_matrix(Zt_wZ, wZ)
    out <- list(mu = mu, w = w, a = a, D = D)
  } else {
    a <- expression_vector - (expression_vector * mu + theta * mu)/(theta + mu)
    b <- (theta * mu)/(theta + mu)
    out <- list(mu = mu, a = a, b = b)
  }
  return(out)
}
