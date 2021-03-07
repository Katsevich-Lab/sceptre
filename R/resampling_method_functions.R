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
  skew_t_fit <- tryCatch(sn::selm(t_nulls ~ 1, family = "ST"), error = function(e) return(NA))
  if (class(skew_t_fit) == "selm") { # If the fit worked,
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
#' @param gRNA_precomp a vector of conditional probabilities for gRNA assignments
#' @param gene_precomp_size the pre-computed size
#' @param gene_precomp_offsets the pre-computed distillation offsets
#' @param B the number of resamples to make (default 500)
#' @param seed an arguement to set.seed; if null, no seed is set
#' @param side an argument to set test for left-sided, right-sided or two-tailed. Default as 'left' and can take 'left', 'right' or 'both'.
#' @param reduced_output logical indicating whether the function should the function return the full output or the reduced output.
#' @param verbose print status updates to the console?
#'
#' @noRd
#' @return if reduced_output, a dataframe with the p-value, test statistic, and other helpful values; if !reduced_output, a list containing all of the above, plus the resampled statistics.
run_sceptre_using_precomp <- function(expressions, gRNA_indicators, gRNA_precomp, side, gene_precomp_size, gene_precomp_offsets, B, seed, reduced_output, verbose) {
  if (!is.null(seed)) set.seed(seed)

  # compute the test statistic on the real data
  fit_star <- VGAM::vglm(formula = expressions[gRNA_indicators == 1] ~ 1, family = VGAM::negbinomial.size(gene_precomp_size), offset = gene_precomp_offsets[gRNA_indicators == 1])
  t_star <- VGAM::summaryvglm(fit_star)@coef3["(Intercept)", "z value"]

  # Define a closure to resample B times (omitting the NAs)
  resample_B_times <- function(my_B) {
    t_nulls <- sapply(1:my_B, function(i) {
      if (verbose && (i %% 100 == 0)) cat(paste0("Running resample ", i ,"/", my_B, ".\n"))
      gRNA_indicators_null <- stats::rbinom(n = length(gRNA_precomp), size = 1, prob = gRNA_precomp)
      tryCatch({
        fit_null <- VGAM::vglm(formula = expressions[gRNA_indicators_null == 1] ~ 1, family = VGAM::negbinomial.size(gene_precomp_size), offset = gene_precomp_offsets[gRNA_indicators_null == 1])
        VGAM::summaryvglm(fit_null)@coef3["(Intercept)", "z value"]},
        error = function(e) return(NA),
        warning = function(w) return(NA)
      )
    })
    t_nulls[!is.na(t_nulls)]
  }

  # resample B times
  t_nulls <- resample_B_times(B)

  # obtain a p-value
  skew_t_fit <- fit_skew_t(t_nulls, t_star, side)

  # determine if the fit was successful
  if (!skew_t_fit$skew_t_fit_success || length(t_nulls) <= floor(0.95 * B)) { # If the skew-t fit failed, then try again.
    t_nulls_second_set <- resample_B_times(9 * B)
    t_nulls <- c(t_nulls, t_nulls_second_set)
    skew_t_fit <- fit_skew_t(t_nulls, t_star, side)
  }

  # Prepare the output
  if (reduced_output) {
    out <- data.frame(p_value = skew_t_fit$out_p, skew_t_fit_success = skew_t_fit$skew_t_fit_success,
                      xi = skew_t_fit$skew_t_mle[["xi"]], omega = skew_t_fit$skew_t_mle[["omega"]],
                      alpha = skew_t_fit$skew_t_mle[["alpha"]], nu = skew_t_fit$skew_t_mle[["nu"]],
                      z_value = t_star, n_successful_resamples = length(t_nulls))
  } else {
    out <- list(p_value = skew_t_fit$out_p,
                skew_t_fit_success = skew_t_fit$skew_t_fit_success,
                skew_t_mle = skew_t_fit$skew_t_mle,
                z_value = t_star,
                resampled_z_values = t_nulls)
  }
  return(out)
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
run_gRNA_precomputation <- function(gRNA_indicators, covariate_matrix) {
  fit_model_grna <- stats::glm(gRNA_indicators ~ ., family = stats::binomial(), data = covariate_matrix)
  out <- as.numeric(stats::fitted(fit_model_grna))
  return(out)
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


#' Run SCEPTRE on a gRNA-gene pair
#'
#' Runs SCEPTRE on a given gRNA-gene pair. The function requires as arguments the gene expression vector, the gRNA indicator vector, and the covariate matrix.
#'
#' @param expressions a vector a gene expressions
#' @param gRNA_indicators a vector of gRNA inicators
#' @param covariate_matrix the data frame of cell-specific covariates
#' @param side (optional; default "left") side of the test; one of "left," "right," and "both."
#' @param gene_precomp_size (optional) the pre-computed size of the gene NB distribution; if not supplied, will be estimated
#' @param B (optional; default 500) number of conditional randomization test resamples
#' @param reduced_output (optional; default TRUE) return the reduced output?
#' @param verbose (optional; default TRUE) print status updates?
#' @param seed (optional; default 4) seed to the random number generator
#' @return If reduced_output is TRUE, a data frame containing the following columns: p-value, skew-t fit success, skew-t fit MLEs (xi, omega, alpha, nu), the original signed test statistic, and the number of successful resamples; if reduced_output is FALSE, return a list containing all of the above, plus the vector of resampled test statistics.
#' @export
#' @examples
#' data(expressions)
#' data(gRNA_indicators)
#' data(covariate_matrix)
#' result <- run_sceptre_gRNA_gene_pair(expressions,
#' gRNA_indicators,
#' covariate_matrix)
run_sceptre_gRNA_gene_pair <- function(expressions, gRNA_indicators, covariate_matrix, side = "left", gene_precomp_size = NULL, B = 500, seed = NULL, reduced_output = TRUE, verbose = TRUE) {
  cat(paste0("Running gRNA precomputation.\n"))
  gRNA_precomp <- run_gRNA_precomputation(gRNA_indicators, covariate_matrix)

  cat(paste0("Running gene precomputation.\n"))
  gene_precomp <- run_gene_precomputation(expressions, covariate_matrix, gene_precomp_size)

  out <- run_sceptre_using_precomp(expressions,
                                   gRNA_indicators,
                                   gRNA_precomp,
                                   side,
                                   gene_precomp$gene_precomp_size,
                                   gene_precomp$gene_precomp_offsets,
                                   B,
                                   if(is.null(seed)) 4 else seed,
                                   reduced_output,
                                   verbose)
  return(out)
}
