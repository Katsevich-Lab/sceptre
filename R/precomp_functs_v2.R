#' Run gene precomputation (v2)
#'
#' Runs precomputation on a gene
#'
#' @param expressions the numeric vector of gene expressions
#' @param covariate_matrix the covariate matrix on which to regress (NOTE: should contain an interecept term)
#' @param fam a string indicating the family to use in the regression (either "nb" or "poisson")
#'
#' @return a named vector of fitted coefficients, alongside the fitted size parameter (named "gene_theta") if fam == "nb"
run_gene_precomputation_v2 <- function(expressions, covariate_matrix, fam) {
  # second backup: method of moments
  backup_2 <- function(pois_fit) {
    MASS::theta.mm(y = expressions, mu = pois_fit$fitted.values, dfr = pois_fit$df.residual)
  }

  # first backup: MLE on poisson reg
  backup <- function() {
    pois_fit <- stats::glm(expressions ~ . + 0, data = covariate_matrix, family = stats::poisson())
    gene_theta <- tryCatch({
      MASS::theta.ml(expressions, pois_fit$fitted.values, limit = 50)[1]
    }, error = function(e) backup_2(pois_fit), warning = function(w) backup_2(pois_fit))
    gene_theta <- max(gene_theta, 0.1)
    fit_nb <- VGAM::vglm(formula = expressions ~ . + 0, family = VGAM::negbinomial.size(gene_theta), data = covariate_matrix)
    fitted_coefs <- stats::coef(fit_nb)
    return(c(fitted_coefs = fitted_coefs, gene_theta = gene_theta))
  }

  # try to fit a negative binomial GLM with unknown dispersion
  if (fam == "nb") {
    result <- tryCatch({
      fit_nb_init <- MASS::glm.nb(formula = expressions ~ . + 0, data = covariate_matrix)
      gene_theta <- max(fit_nb_init$theta, 0.1)
      fit_nb <- VGAM::vglm(formula = expressions ~ . + 0, family = VGAM::negbinomial.size(gene_theta), data = covariate_matrix)
      fitted_coefs <- stats::coef(fit_nb)
      return(c(fitted_coefs, gene_theta = gene_theta))
    }, error = function(e) backup(), warning = function(w) backup())
  }
  if (fam == "poisson") {
    pois_fit <- stats::glm(expressions ~ . + 0, data = covariate_matrix, family = stats::poisson())
    result <- stats::coef(pois_fit)
  }
  return(result)
}


#' Run grna precomputation (v2)
#'
#' Runs precomputation on a vector of grna indicators
#'
#' @param indicators the binary vector of grna indicators
#' @param covariate_matrix the covariate matrix on which to regress (NOTE: should contain an interecept term)
#'
#' @return a named vector of fitted coefficients
run_grna_precomputation_v2 <- function(indicators, covariate_matrix) {
  fit <- stats::glm(formula = indicators ~ . + 0, family = "binomial", data = covariate_matrix)
  ret <- stats::coef(fit)
  return(ret)
}


#' Generate synthetic data
#'
#' Generate a (sparse) matrix of synthetic grna indicators from a vector of fitted probabilitites.
#'
#' @param fitted_probs a vector of fitted grna presence/absence indicators
#' @param B the number of resampled vectors to create
#'
#' @return a sparse matrix (of dimension n x B) containing the resampled indicators
generate_synthetic_grna_data <- function(fitted_probs, B) {
  draws_idxs <- lapply(X = fitted_probs, FUN = function(prob) {
    draws <- stats::rbinom(n = B, size = 1, prob = prob)
    which(draws == 1)
  })
  n_cells <- length(fitted_probs)
  row_idxs <- rep(x = seq(1L, n_cells), times = sapply(draws_idxs, length))
  col_idxs <- unlist(draws_idxs)
  synth_data <- Matrix::sparseMatrix(i = row_idxs, j = col_idxs, dims = c(n_cells, B))
  return(synth_data)
}
