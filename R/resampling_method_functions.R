#' Fit skew-t
#'
#' Fits a skew-t distribution on a set of resampled test statistics.
#'
#' @param t_nulls sampled statistics that form a null distribution
#' @param t_star the observed test statistic
#' @param side the side of the test (left, right, of both)
#'
#' @return
#' A list containing (i) skew_t_fit_success (boolean), (ii) out_p (the p-value), and (iii) skew-t mle (a vector containing the fitted MLEs, NA if fit failed).
#' @noRd
fit_skew_t <- function(t_nulls, t_star, side) {
  p_value_skew_t <- NA
  skew_t_fit <- tryCatch(
    R.utils::withTimeout(expr = sn::selm(t_nulls ~ 1, family = "ST"), timeout = 3),
    error = function(e) return(NA),
    warning = function(w) return(NA))
  if (is(skew_t_fit, "selm")) { # If the fit worked,
    dp <- skew_t_fit@param$dp # then extract the parameters.
    if (!any(is.na(dp))) { # If all the fitted parameters are numbers,
      p_value_skew_t <- switch(side,
                               'left' = pmax(.Machine$double.eps, sn::pst(x = t_star, dp = dp, method = 2, rel.tol = .Machine$double.eps)), # then compute the skew t-based p-value. pst(x = t_star, dp = dp)
                               'right' = pmax(.Machine$double.eps, 1 - sn::pst(x = t_star, dp = dp, method = 2, rel.tol = .Machine$double.eps)),
                               'both' = pmax(.Machine$double.eps, sn::pst(x = -abs(t_star), dp = dp, method = 2, rel.tol = .Machine$double.eps) +
                                               (1 - sn::pst(x = abs(t_star), dp = dp, method = 2, rel.tol = .Machine$double.eps)))
      )
    }
  }
  # check if the skew-t fit worked
  skew_t_fit_success <- !is.na(p_value_skew_t)
  if (skew_t_fit_success) {
    out_p <- p_value_skew_t
    skew_t_mle <- dp
  } else {
    out_p <- switch(side,
                    'left' = mean(c(-Inf, t_nulls) <= t_star),
                    'right' = mean(c(Inf, t_nulls) >= t_star),
                    'both' = mean(c(Inf, abs(t_nulls)) >= abs(t_star)))
    skew_t_mle <- c(xi = NA, omega = NA, alpha = NA, nu = NA)
  }
  return(list(skew_t_fit_success = skew_t_fit_success, out_p = out_p, skew_t_mle = skew_t_mle))
}

#' Run sceptre using precomputations for gRNAs and genes.
#'
#' This function is the workhorse function of the sceptre package. It runs a distilled CRT using a negative binomial test statistic based on an expression vector, a gRNA indicator vector, an offset vector (from the distillation step), gRNA conditional probabilities, an estimate of the negative binomial size parameter, and the number of resampling replicates.
#'
#' @param expressions a vector of gene expressions (in UMI counts)
#' @param gRNA_indicators a vector of gRNA indicators
#' @param gRNA_precomp a n_cells x B matrix of synthetic indicators
#' @param gene_precomp_size the pre-computed size
#' @param gene_precomp_offsets the pre-computed distillation offsets
#' @param side an argument to set test for left-sided, right-sided or two-tailed. Default as 'left' and can take 'left', 'right' or 'both'.#'
#' @noRd
#' @return if reduced_output, a dataframe with the p-value, test statistic, and other helpful values; if !reduced_output, a list containing all of the above, plus the resampled statistics.
run_sceptre_using_precomp_fast <- function(expressions, gRNA_indicators, gRNA_precomp, side, gene_precomp_size, gene_precomp_offsets, full_output = FALSE) {
  exp_gene_offsets <- exp(gene_precomp_offsets)
  # compute test statistic on real data
  y <- expressions[gRNA_indicators == 1]
  exp_o <- exp_gene_offsets[gRNA_indicators == 1]
  log_fold_change <- log(mean(y)) - log(mean(exp_o))
  z_star <- compute_nb_test_stat_fast_score(y, exp_o, gene_precomp_size)
  # compute test statistics on the resampled statistics
  z_null <- apply(X = gRNA_precomp, MARGIN = 2, FUN = function(col) {
    compute_nb_test_stat_fast_score(expressions[col], exp_gene_offsets[col], gene_precomp_size)
  })
  # fit skew-t
  skew_t_fit <- fit_skew_t(z_null, z_star, side)
  # create output
  if (full_output) {
    B <- length(z_null)
    z_df <- as.data.frame(matrix(z_null, nrow = 1, ncol = B))
    colnames(z_df) <- paste0("z_null_", seq(1, B))
    out <- data.frame(p_value = skew_t_fit$out_p, skew_t_fit_success = skew_t_fit$skew_t_fit_success,
                      xi = skew_t_fit$skew_t_mle[["xi"]], omega = skew_t_fit$skew_t_mle[["omega"]],
                      alpha = skew_t_fit$skew_t_mle[["alpha"]], nu = skew_t_fit$skew_t_mle[["nu"]],
                      z_value = z_star, log_fold_change = log_fold_change) %>% dplyr::mutate(z_df)
  } else {
    out <- data.frame(p_value = skew_t_fit$out_p, z_value = z_star, log_fold_change = log_fold_change)
  }
  return(out)
}


compute_nb_test_stat_fast_score <- function(y, exp_o, gene_precomp_size) {
  r_exp_o <- gene_precomp_size * exp_o
  y_exp_o <- y * exp_o
  r_plus_exp_o <- gene_precomp_size + exp_o
  sum_y <- sum(y)
  top <- (y_exp_o + r_exp_o)/r_plus_exp_o
  bottom <- r_exp_o/r_plus_exp_o
  z <- (sum_y - sum(top))/sqrt(sum(bottom))
  return(z)
}


#' Run gRNA precomputation
#'
#' This function runs the precomputation for a given gRNA.
#'
#' @param gRNA_indicators a vector of gRNA indicators
#' @param covariate_matrix the cell-specific covariate matrix
#'
#' @noRd
#' @return the fitted probabilities
run_gRNA_precomputation <- function(gRNA_indicators, covariate_matrix, B, seed) {
  set.seed(seed)
  # first, fit a logistic regression model to estimate perturbation probabilities
  fit_model_grna <- stats::glm(gRNA_indicators ~ ., family = stats::binomial(), data = covariate_matrix)
  out <- as.numeric(stats::fitted(fit_model_grna))
  # generate synthetic data
  draws_idxs <- lapply(X = out, FUN = function(prob) {
    draws <- stats::rbinom(n = B, size = 1, prob = prob)
    which(draws == 1)
  })
  n_cells <- length(gRNA_indicators)
  row_idxs <- rep(x = seq(1L, n_cells), times = sapply(draws_idxs, length))
  col_idxs <- unlist(draws_idxs)
  synth_data <- Matrix::sparseMatrix(i = row_idxs, j = col_idxs, dims = c(n_cells, B))
  return(synth_data)
}


#' Run gene precomputation
#'
#' This function runs the precomputation for a given gene. It fits an NB regression of expression on covariates. The estimate of theta (i.e., the NB size parameter) is obtained from the glm.nb function. Offsets are obtained by log-transforming the fitted values.
#'
#' @param expressions the vector of gene expressions
#' @param covariate_matrix the cell-specific covariate matrix
#' @param gene_precomp_size the pre-computed size parameter
#'
#' @return a named list containing two items: offsets and size.
#' @noRd
run_gene_precomputation <- function(expressions, covariate_matrix, gene_precomp_size) {
  # cases on gene_precomp_size
  if (is.null(gene_precomp_size)) {
    # no size supplied; use glm.nb to estimate size and fit model
    # second backup: method of moments on poisson reg
    backup_2 <- function(pois_fit) {
      MASS::theta.mm(y = expressions, mu = pois_fit$fitted.values, dfr = pois_fit$df.residual)
    }

    # first backup: MLE on poisson reg
    backup <- function() {
      pois_fit <- stats::glm(expressions ~ ., data = covariate_matrix, family = stats::poisson())
      gene_precomp_size_out <- tryCatch({
        MASS::theta.ml(expressions, pois_fit$fitted.values, limit = 50)[1]
      }, error = function(e) backup_2(pois_fit), warning = function(w) backup_2(pois_fit))
      fit_nb <- VGAM::vglm(formula = expressions ~ ., family = VGAM::negbinomial.size(gene_precomp_size_out), data = covariate_matrix)
      fitted_vals <- as.numeric(VGAM::fittedvlm(fit_nb))
      list(fitted_vals = fitted_vals, gene_precomp_size_out = gene_precomp_size_out)
    }

    # try to fit a negative binomial GLM with unknown dispersion
    result <- tryCatch({
      fit_nb <- MASS::glm.nb(formula = expressions ~ ., data = covariate_matrix)
      fitted_vals <- as.numeric(fit_nb$fitted.values)
      gene_precomp_size_out <- fit_nb$theta
      list(fitted_vals = fitted_vals, gene_precomp_size_out = gene_precomp_size_out)
    }, error = function(e) backup(), warning = function(w) backup())

    fitted_vals <- result$fitted_vals; gene_precomp_size_out <- result$gene_precomp_size_out

  } else { # size supplied; use vglm to fit model
    gene_precomp_size_out <- gene_precomp_size
    fit_nb <- VGAM::vglm(formula = expressions ~ ., family = VGAM::negbinomial.size(gene_precomp_size_out), data = covariate_matrix)
    fitted_vals <- as.numeric(VGAM::fittedvlm(fit_nb))
  }

  gene_precomp_offsets_out <- log(fitted_vals)
  out <- list(gene_precomp_offsets = gene_precomp_offsets_out, gene_precomp_size = gene_precomp_size_out)
  return(out)
}


#' Run `sceptre` on a gRNA-gene pair
#'
#' Runs `sceptre` on a given gRNA-gene pair.
#'
#' This function is for demonstration purposes only. The function `run_sceptre_in_memory` should be used when testing for association between multiple genes and gRNAs.
#'
#' @param gene_expressions (numeric vector) a vector of gene expressions
#' @param gRNA_expressions (numeric vector) a vector of gRNA expressions (or, optionally, a user-thresholded vector of binary gRNA indicators)
#' @param covariate_matrix (data frame) the cell-specific matrix of covariates
#' @param side sidedness of the test, one of "left," "right," and "both"
#' @param B number of resamples for the conditional randomization test
#' @param full_output return the full output (TRUE) or a streamlined, reduced output (FALSE)?
#' @param seed seed to the random number generator
#'
#' @examples
#' \dontrun{
#' # load the example data
#' data(gene_matrix); data(gRNA_matrix); data(covariate_matrix)
#' gene_expressions <- gene_matrix[1,]
#' gRNA_expressions <- gRNA_matrix[1,]
#' # run method
#' result <- run_sceptre_gRNA_gene_pair(gene_expressions, gRNA_expressions, covariate_matrix, "left")
#' }
run_sceptre_gRNA_gene_pair <- function(gene_expressions, gRNA_expressions, covariate_matrix, side = "both", B = 1500, full_output = FALSE, seed = 4) {
  THRESHOLD <- 3

  cat("Running perturbation precomputation.")
  # threshold gRNA_expressions, if necessary
  if (max(gRNA_expressions) >= 2) gRNA_expressions <- as.integer(gRNA_expressions >= THRESHOLD)
  gene_precomp <- run_gene_precomputation(gene_expressions, covariate_matrix, NULL)
  cat(crayon::green(' \u2713\n'))

  cat("Running gRNA precomputation.")
  gRNA_precomp <- run_gRNA_precomputation(gRNA_expressions, covariate_matrix, B, seed)
  cat(crayon::green(' \u2713\n'))

  cat("Running perturbation-to-gene association analysis.")
  out <- run_sceptre_using_precomp_fast(gene_expressions,
                                        gRNA_expressions,
                                        gRNA_precomp,
                                        side,
                                        gene_precomp$gene_precomp_size,
                                        gene_precomp$gene_precomp_offsets,
                                        full_output)
  cat(crayon::green(' \u2713\n'))
  return(out)
}
