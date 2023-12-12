assign_grnas_to_cells_mixture <- function(grna_matrix, cell_covariate_data_frame, grna_assignment_hyperparameters, print_progress, parallel, n_processors, log_dir) {
  if (!parallel) cat(crayon::red("Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`.\n\n"))
  # 0. get random starting guesses for pi and g_pert
  starting_guesses <- get_random_starting_guesses(n_em_rep = grna_assignment_hyperparameters$n_em_rep,
                                                  pi_guess_range = grna_assignment_hyperparameters$pi_guess_range,
                                                  g_pert_guess_range = grna_assignment_hyperparameters$g_pert_guess_range)

  # 1. obtain the covariate matrix
  covariate_matrix <- convert_covariate_df_to_design_matrix(covariate_data_frame = cell_covariate_data_frame,
                                                            formula_object = grna_assignment_hyperparameters$formula_object)

  # 2. make the grna expression matrix row-accessible
  grna_matrix <- set_matrix_accessibility(grna_matrix, make_row_accessible = TRUE)
  grna_ids <- rownames(grna_matrix)

  # 3. define the function to perform assignments for a set of grnas
  analyze_given_grna_ids <- function(curr_grna_ids, proc_id = NULL) {
    if (parallel && print_progress) {
      f_name <- paste0(get_log_dir(log_dir), "assign_grnas_", proc_id, ".out")
      file.create(f_name) |> invisible()
    } else {
      f_name <- NULL
    }
    initial_assignment_list <- sapply(seq_along(curr_grna_ids), function(i) {
      grna_id <- get_id_from_idx(response_idx = i, print_progress = print_progress, response_ids = curr_grna_ids,
                                 feature = "gRNA", str = "Performing gRNA-to-cell assignments for ", parallel = parallel, f_name = f_name)
      g <- load_csr_row(j = grna_matrix@j, p = grna_matrix@p, x = grna_matrix@x,
                        row_idx = which(grna_id == grna_ids), n_cells = ncol(grna_matrix))
      assignments <- obtain_em_assignments(pi_guesses = starting_guesses$pi_guesses,
                                           g_pert_guesses = starting_guesses$g_pert_guesses,
                                           g = g, covariate_matrix = covariate_matrix, use_glm = TRUE,
                                           n_nonzero_cells_cutoff = grna_assignment_hyperparameters$n_nonzero_cells_cutoff,
                                           backup_threshold = grna_assignment_hyperparameters$backup_threshold,
                                           probability_threshold = grna_assignment_hyperparameters$probability_threshold)
      return(assignments)
    }, simplify = FALSE) |> stats::setNames(curr_grna_ids)
  }

  # 4. run the analysis
  if (!parallel) {
    initial_assignment_list <- analyze_given_grna_ids(grna_ids)
  } else {
    cat("Running gRNA assignments in parallel. ")
    if (print_progress) {
      cat(paste0("Change directories to ", crayon::blue(get_log_dir(log_dir)), " and view the files ",
                 crayon::blue("assign_grnas_*.out"), " for progress updates.\n"))
    }
    grna_ids_partitioned <- partition_response_ids(grna_ids, parallel, n_processors)
    res <- parallel::mclapply(seq_along(grna_ids_partitioned),
                              function(proc_id) analyze_given_grna_ids(grna_ids_partitioned[[proc_id]], proc_id),
                              mc.cores = length(grna_ids_partitioned))
    initial_assignment_list <- res |> purrr::flatten()
  }

  # return
  return (initial_assignment_list)
}

obtain_em_assignments <- function(pi_guesses, g_pert_guesses, g, covariate_matrix, use_glm,
                                  n_nonzero_cells_cutoff, backup_threshold, probability_threshold) {
  threshold_em_probabilities <- TRUE
  n_nonzero <- sum(g >= 1)
  if (n_nonzero >= n_nonzero_cells_cutoff) {
    # 1. estimate the conditional means
    if (use_glm) {
      pois_fit <- suppressWarnings(stats::glm.fit(y = g, x = covariate_matrix, family = stats::poisson()))
      g_mus_pert0 <- pois_fit$fitted.values
    }
    # 2. compute log(g!)
    log_g_factorial <- lgamma(g+1)
    # 3. run the em algorithm
    fit_cpp <- run_reduced_em_algo_cpp(pi_guesses, g_pert_guesses, g, g_mus_pert0, log_g_factorial)
    # obtain the assignments
    if (fit_cpp$outer_converged && fit_cpp$outer_log_lik != -Inf) {
      assignments <- which(fit_cpp$outer_Ti1s >= probability_threshold)
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


get_random_starting_guesses <- function(n_em_rep = 5, pi_guess_range = c(1e-5, 0.1), g_pert_guess_range = c(10, 5000)) {
  set.seed(4)
  pi_guesses <- stats::runif(n = n_em_rep, min = pi_guess_range[1], max = pi_guess_range[2])
  g_pert_guesses <- stats::runif(n = n_em_rep, min = g_pert_guess_range[1], max = g_pert_guess_range[2])
  list(pi_guesses = pi_guesses, g_pert_guesses = g_pert_guesses)
}
