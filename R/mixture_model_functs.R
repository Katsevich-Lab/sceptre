assign_grnas_to_cells_mixture <- function(grna_matrix, grna_group_data_frame, cell_covariate_data_frame, n_em_rep = 5L) {
  # 0. get random starting guesses for pi and g_pert
  starting_guesses <- get_random_starting_guesses(n_em_rep)
  pi_guesses <- starting_guesses$pi_guesses
  g_pert_guesses <- starting_guesses$g_pert_guesses

  # 1. obtain the covariate matrix
  formula_object <- auto_construct_formula_object(cell_covariates = cell_covariate_data_frame,
                                                  include_grna_covariates = TRUE)
  covariate_matrix <- convert_covariate_df_to_design_matrix(covariate_data_frame = cell_covariate_data_frame,
                                                            formula_object = formula_object)
  # 2. make the grna expression matrix row-accessible
  grna_matrix <- set_matrix_accessibility(grna_matrix, make_row_accessible = TRUE)
  grna_ids <- rownames(grna_matrix)

  # 3. iterate over grna ids, obtaining the assignments for each
  assignment_list <- sapply(grna_ids, function(grna_id) {
    print(grna_id)
    g <- load_csr_row(j = grna_matrix@j, p = grna_matrix@p, x = grna_matrix@x,
                      row_idx = which(grna_id == grna_ids), n_cells = ncol(grna_matrix))
    assignments <- obtain_em_assignments(pi_guesses, g_pert_guesses, g, covariate_matrix,
                                         use_glm = TRUE, n_nonzero_cells_cutoff = 10L, backup_threshold = 5)
    return(assignments)
  }, simplify = FALSE)
}


obtain_em_assignments <- function(pi_guesses, g_pert_guesses, g, covariate_matrix, use_glm,
                                  n_nonzero_cells_cutoff, backup_threshold) {
  threshold_em_probabilities <- TRUE
  n_nonzero <- sum(g >= 1)
  if (n_nonzero >= n_nonzero_cells_cutoff) {
    # 1. estimate the conditional means
    if (use_glm) {
      pois_fit <- suppressWarnings(stats::glm.fit(y = g, x = covariate_matrix, family = stats::poisson()))
      g_mus_pert0 <- pois_fit$fitted.values
    } else {
      stop("not yet implemented.")
    }
    # 2. compute log(g!)
    log_g_factorial <- lgamma(g+1)
    # 3. run the em algorithm
    fit_cpp <- run_reduced_em_algo_cpp(pi_guesses, g_pert_guesses, g, g_mus_pert0, log_g_factorial)
    # obtain the assignments
    if (fit_cpp$outer_converged && fit_cpp$outer_log_lik != -Inf) {
      assignments <- which(fit_cpp$outer_Ti1s >= 0.5)
    } else {
      threshold_em_probabilities <- FALSE
    }
  } else {
    threshold_em_probabilities <- FALSE
  }

  if (!threshold_em_probabilities) {
    assignments <- which(g >= backup_threshold)
  }

  return(assignments)
}


get_random_starting_guesses <- function(n_em_rep = 5, pi_guess_range = c(1e-5, 0.1), exp_g_pert_guess_range = c(10, 5000)) {
  set.seed(4)
  pi_guesses <- stats::runif(n = n_em_rep, min = pi_guess_range[1], max = pi_guess_range[2])
  g_pert_guesses <- stats::runif(n = n_em_rep, min = log(exp_g_pert_guess_range[1]), max = log(exp_g_pert_guess_range[2]))
  list(pi_guesses = pi_guesses, g_pert_guesses = g_pert_guesses)
}
