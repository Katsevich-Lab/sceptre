assign_grnas_to_cells_mixture <- function(grna_matrix, grna_group_data_frame, cell_covariate_data_frame) {
  # a. get random starting guesses for pi and g_pert
  set.seed(4)
  n_em_rep <- 10L
  pi_guess_range <- c(1e-5, 0.1)
  exp_g_pert_guess_ranges = list(c(10, 500), c(500, 10000))
  n_guesses_per_range <- ceiling(n_em_rep/2)
  n_total_em_rep <- n_guesses_per_range * 2
  pi_guesses <- stats::runif(n = n_guesses_per_range * 2, min = pi_guess_range[1], max = pi_guess_range[2])
  g_pert_guesses <- lapply(X = exp_g_pert_guess_ranges,
                                   FUN = function(l) stats::runif(n = n_guesses_per_range, min = l[1], max = l[2])) |>
    unlist() |> log()

  # a. obtain the covariate matrix
  formula_object <- auto_construct_formula_object(cell_covariates = cell_covariate_data_frame,
                                                  include_grna_covariates = TRUE)
  covariate_matrix <- convert_covariate_df_to_design_matrix(covariate_data_frame = cell_covariate_data_frame,
                                                            formula_object = formula_object)
  # b. make the grna expression matrix row-accessible
  grna_matrix <- set_matrix_accessibility(grna_matrix, make_row_accessible = TRUE)
  grna_ids <- rownames(grna_matrix)

  # c. iterate over grna ids, fitting a mixture model to each and obtaining the assignments
  assignment_list <- sapply(grna_ids, function(grna_id) {
    print(grna_id)
    # 1. load the grna expressions
    g <- load_csr_row(j = grna_matrix@j,
                      p = grna_matrix@p,
                      x = grna_matrix@x,
                      row_idx = which(grna_id == grna_ids),
                      n_cells = ncol(grna_matrix))
    head(g, 400)
    fit <- fit_latent_variable_model(g, covariate_matrix, pi_guesses, g_pert_guesses)
    which(fit$Ti1s >= 0.5)
  }, simplify = FALSE)
}


fit_latent_variable_model <- function(g, covariate_matrix, pi_guesses, g_pert_guesses) {
  pois_fit <- suppressWarnings(stats::glm.fit(y = g, x = covariate_matrix, family = stats::poisson()))
  g_fitted <- pois_fit$fitted.values
  fit <- run_reduced_em_algo_v2(pi_guesses, g_pert_guesses, g, g_fitted)
  return(best_fit)
}


run_reduced_em_algo_v2 <- function(pi_guesses, g_pert_guesses, g, g_fitted) {
  # define options as well as "outer" Ti1, g_pert, and pi
  ep_tol <- 0.5 * 1e-4
  max_it <- 50
  min_it <- 3
  outer_Ti1s <- numeric()
  outer_g_pert <- NA_integer_
  outer_pi <- NA_integer_
  outer_log_lik <- -Inf
  outer_i <- NA_integer_

  # loop over initial guesses
  for (i in seq_along(pi_guesses)) {
    converged <- FALSE
    prev_log_lik <- -Inf
    n <- length(g)
    curr_g_pert <- g_pert_guesses[i]
    curr_pi <- pi_guesses[i]
    iteration <- 1L

    # Iterate through E and M steps until convergence.
    while (TRUE) {
      ##############
      # E step start
      ##############
      g_mus_pert0 <- g_fitted
      g_mus_pert1 <- g_fitted + curr_g_pert

      # compute log likelihood
      p0 <- exp(log(1 - curr_pi) + stats::dpois(g, g_mus_pert0, log = TRUE))
      p1 <- exp(log(curr_pi) + stats::dpois(g, g_mus_pert1, log = TRUE))
      s <- p0 + p1
      if (0 %in% s) {
        s_wo_0 <- s[s != 0]
        s[s == 0] <- min(s_wo_0)
      }
      curr_log_lik <- sum(log(s))

      # compute membership probabilities
      quotient <- log(1 - curr_pi) - log(curr_pi) + g * (log(g_mus_pert0) - log(g_mus_pert1)) + g_mus_pert1 - g_mus_pert0
      Ti1s <- 1/(exp(quotient) + 1)

      # verify that Ti1s are not NA or 0
      if (any(is.na(Ti1s)) || all(Ti1s == 0)) {
        curr_log_lik <- -Inf
        break()
      }

      # estimate new_pi and ensure new_pi is less than 0.5
      curr_pi <- sum(Ti1s)/n
      if (curr_pi > 0.5) {
        Ti1s <- 1 - Ti1s
        curr_pi <- 1 - curr_pi
      }

      # Check for convergence
      tol <- compute_tolerance(curr_log_lik, prev_log_lik)
      if (tol < ep_tol && iteration >= min_it) {
        converged <- TRUE
        break()
      } else {
        prev_log_lik <- curr_log_lik
        iteration <- iteration + 1L
        if (iteration >= max_it) break()
      }
      ############
      # E step end
      ############

      # M step
      e1 <- sum(Ti1s * g)
      e2 <- sum(Ti1s * g_fitted)
      curr_g_pert <- log(e1) - log(e2)
    }

    if (curr_log_lik > outer_log_lik && converged) {
      outer_g_pert <- curr_g_pert
      outer_pi <- curr_pi
      outer_Ti1s <- Ti1s
      outer_i <- i
    }
  }
  return(list(g_pert = outer_g_pert, pi = outer_pi, Ti1s = outer_Ti1s))
}


compute_tolerance <- function(curr_log_lik, prev_log_lik) {
  if (curr_log_lik == -Inf) {
    tol <- Inf
  } else {
    tol <- abs(curr_log_lik - prev_log_lik)/min(abs(curr_log_lik), abs(prev_log_lik))
  }
  return(tol)
}
