#' Assign guides using the existing covariate-adjusted count mixture.
#'
#' @param grna_matrix Supported row-readable guide-by-source-cell counts.
#' @param cell_covariate_data_frame Covariates in the same source-cell order.
#' @param grna_assignment_hyperparameters Existing EM, formula and threshold options.
#' @param print_progress,parallel,n_processors,log_dir Worker/reporting options.
#' @param grna_ids Ordered guide IDs to process, retained as output names.
#' @return Named lists of one-based source-cell assignments, in serial or existing
#'   seeded block order. Reads guide rows, may create progress files and runs GC
#'   every 200 block items. Fixed starting guesses preserve caller RNG; active
#'   count/nonconvergence fallback and formula-error context remain unchanged.
#' @noRd
assign_grnas_to_cells_mixture <- function(
    grna_matrix,
    cell_covariate_data_frame,
    grna_assignment_hyperparameters,
    print_progress,
    parallel,
    n_processors,
    log_dir,
    grna_ids
) {
    if (!parallel) {
        message(crayon::red(
            "Note: If you are on a Mac laptop or desktop, consider setting `parallel = TRUE` to improve speed. Otherwise, keep `parallel = FALSE`."
        ))
    }
    # 0. get random starting guesses for pi and g_pert
    starting_guesses <- get_random_starting_guesses(
        n_em_rep = grna_assignment_hyperparameters$n_em_rep,
        pi_guess_range = grna_assignment_hyperparameters$pi_guess_range,
        g_pert_guess_range = grna_assignment_hyperparameters$g_pert_guess_range
    )

    # 1. obtain the covariate matrix
    covariate_matrix <- tryCatch(
        convert_covariate_df_to_design_matrix(
            covariate_data_frame = cell_covariate_data_frame,
            formula_object = grna_assignment_hyperparameters$formula_object
        ),
        error = function(cnd) {
            underlying_error <- conditionMessage(cnd)
            cnd$message <- paste0(
                "The `formula_object` used for gRNA assignment could not be ",
                "converted to a valid design matrix. Check the formula and ",
                "covariate data. Possible causes include missing variables or ",
                "values, non-finite values, categorical covariates with fewer ",
                "than two observed levels, and perfectly collinear or ",
                "redundant covariates. Supply a corrected `formula_object` to ",
                "`assign_grnas()`; the formula used in ",
                "`set_analysis_parameters()` may be appropriate.",
                "\nUnderlying error: ",
                underlying_error
            )
            stop(cnd)
        }
    )

    # 2. make the grna expression matrix row-accessible
    grna_matrix <- set_matrix_accessibility(
        grna_matrix,
        make_row_accessible = TRUE
    )

    #' Assign one ordered block of guide IDs.
    #' @param curr_grna_ids Character guide IDs; proc_id selects the worker log.
    #' @return Named assignment list. Reports then collects garbage before each
    #'   row read; the shared design and fixed guesses are read-only captures.
    #' @noRd
    assign_guide_block <- function(curr_grna_ids, proc_id = NULL) {
        f_name <- analysis_log_file("assign_grnas", proc_id, parallel,
                                    print_progress, log_dir)
        initial_assignment_list <- lapply(
            seq_along(curr_grna_ids),
            function(i) {
                grna_id <- as.character(curr_grna_ids[i])
                report_worker_progress(i, grna_id, length(curr_grna_ids),
                    print_progress, parallel, f_name, feature = "gRNA",
                    prefix = "Performing gRNA-to-cell assignments for ")
                collect_worker_garbage(i)
                g <- load_row(grna_matrix, grna_id)
                assignments <- obtain_em_assignments(
                    pi_guesses = starting_guesses$pi_guesses,
                    g_pert_guesses = starting_guesses$g_pert_guesses,
                    g = g,
                    covariate_matrix = covariate_matrix,
                    n_nonzero_cells_cutoff = grna_assignment_hyperparameters$n_nonzero_cells_cutoff,
                    backup_threshold = grna_assignment_hyperparameters$backup_threshold,
                    probability_threshold = grna_assignment_hyperparameters$probability_threshold
                )
                return(assignments)
            }
        ) |>
            stats::setNames(curr_grna_ids)
    }

    # 4. run the analysis
    if (!parallel) {
        initial_assignment_list <- assign_guide_block(grna_ids)
    } else {
        message("Running gRNA assignments in parallel.")
        if (print_progress) {
            message(
                "Change directories to ",
                crayon::blue(get_log_dir(log_dir)),
                " and view the files ",
                crayon::blue("assign_grnas_*.out"),
                " for progress updates."
            )
        }
        grna_ids_partitioned <- partition_response_ids(
            grna_ids,
            parallel,
            n_processors
        )
        res <- parallel::mclapply(
            seq_along(grna_ids_partitioned),
            function(proc_id) {
                assign_guide_block(grna_ids_partitioned[[proc_id]], proc_id)
            },
            mc.cores = length(grna_ids_partitioned)
        )
        initial_assignment_list <- res |> purrr::flatten()
    }

    # return
    return(initial_assignment_list)
}


#' Assign one guide using EM or the existing raw-count fallback.
#'
#' @param pi_guesses,g_pert_guesses Numeric starting values in matched order.
#' @param g Numeric counts in source-cell order.
#' @param covariate_matrix Matrix aligned to g, including one-column designs.
#' @param n_nonzero_cells_cutoff Minimum count of cells with g >= 1 for EM.
#' @param backup_threshold Raw-count fallback threshold.
#' @param probability_threshold Posterior threshold on an accepted EM fit.
#' @return Increasing one-based cell positions. Uses the existing Poisson
#'   initialization and native EM; expected initialization warnings are muffled.
#'   Low counts or unsuccessful EM retain raw-count thresholding. No input mutation.
#' @noRd
obtain_em_assignments <- function(
    pi_guesses,
    g_pert_guesses,
    g,
    covariate_matrix,
    n_nonzero_cells_cutoff,
    backup_threshold,
    probability_threshold
) {
    threshold_em_probabilities <- TRUE
    n_nonzero <- sum(g >= 1)
    if (n_nonzero >= n_nonzero_cells_cutoff) {
        # 1. estimate the conditional means
        # Initial Poisson warnings are expected; native EM refines this fit.
        pois_fit <- withCallingHandlers(
            stats::glm.fit(y = g, x = covariate_matrix, family = stats::poisson()),
            warning = function(w) invokeRestart("muffleWarning")
        )
        g_mus_pert0 <- pois_fit$fitted.values
        # 2. compute log(g!)
        log_g_factorial <- lgamma(g + 1)
        # 3. run the em algorithm
        fit_cpp <- run_reduced_em_algo_cpp(
            pi_guesses,
            g_pert_guesses,
            g,
            g_mus_pert0,
            log_g_factorial
        )
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


#' Generate the existing deterministic mixture starting guesses.
#'
#' @param n_em_rep Number of paired guesses.
#' @param pi_guess_range,g_pert_guess_range Numeric uniform bounds.
#' @return Lists pi_guesses and g_pert_guesses, drawn in that order under seed 4.
#'   with_seed restores the caller's RNG state, including its absence.
#' @noRd
get_random_starting_guesses <- function(
    n_em_rep = 5,
    pi_guess_range = c(1e-5, 0.1),
    g_pert_guess_range = c(10, 5000)
) {
    withr::with_seed(4, {
        pi_guesses <- stats::runif(
            n = n_em_rep,
            min = pi_guess_range[1],
            max = pi_guess_range[2]
        )
        g_pert_guesses <- stats::runif(
            n = n_em_rep,
            min = g_pert_guess_range[1],
            max = g_pert_guess_range[2]
        )
        list(pi_guesses = pi_guesses, g_pert_guesses = g_pert_guesses)
    })
}
