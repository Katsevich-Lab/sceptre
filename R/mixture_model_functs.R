#' Assign guides with a count mixture adjusted for cell covariates
#' @param grna_matrix Guide-by-cell counts prepared for row access.
#' @param cell_covariate_data_frame Covariates in matching original-cell order.
#' @param grna_assignment_hyperparameters EM, formula and assignment-threshold settings.
#' @param print_progress,parallel,n_processors,log_dir Worker and reporting settings.
#' @param grna_ids Guide names to process, in order.
#' @return Named lists of original cell-column numbers, starting at 1. Serial work
#'   follows input order; parallel work follows the seeded block order. Read one
#'   guide at a time, report progress and collect garbage every 200 items. Starting
#'   guesses use a fixed seed and restore R's random state. Low counts or failed
#'   mixture fits use the existing raw-count fallback.
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

    #' Assign cells for one worker's block of guides
    #' @param curr_grna_ids Guide names in processing order; proc_id names the worker log.
    #' @return Named assignment list. Report progress and collect garbage before
    #'   reading each guide, using the shared design and fixed starting guesses.
    #' @noRd
    assign_guide_block <- function(curr_grna_ids, proc_id = NULL) {
        f_name <- analysis_log_file(
            "assign_grnas",
            proc_id,
            parallel,
            print_progress,
            log_dir
        )
        initial_assignment_list <- lapply(
            seq_along(curr_grna_ids),
            function(i) {
                grna_id <- as.character(curr_grna_ids[i])
                report_worker_progress(
                    i,
                    grna_id,
                    length(curr_grna_ids),
                    print_progress,
                    parallel,
                    f_name,
                    feature = "gRNA",
                    prefix = "Performing gRNA-to-cell assignments for "
                )
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


#' Assign one guide's cells using a fitted mixture or a raw-count fallback
#' @param pi_guesses,g_pert_guesses Paired starting values for EM fits.
#' @param g Guide counts in original-cell order.
#' @param covariate_matrix Design matrix with rows in that same order.
#' @param n_nonzero_cells_cutoff Minimum number of cells with g >= 1 to attempt EM.
#' @param backup_threshold Raw-count threshold used when EM is unsuitable.
#' @param probability_threshold Posterior-probability cutoff after a successful EM fit.
#' @return Original cell positions, starting at 1, in increasing order. Initialize
#'   with a Poisson fit, silencing its expected warnings, then run C++ EM. If too
#'   few cells have counts or EM fails, use raw-count thresholding instead.
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
            stats::glm.fit(
                y = g,
                x = covariate_matrix,
                family = stats::poisson()
            ),
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


#' Generate repeatable starting values for the mixture fits
#' @param n_em_rep Number of paired starting guesses.
#' @param pi_guess_range,g_pert_guess_range Uniform sampling bounds.
#' @return Lists pi_guesses and g_pert_guesses, sampled in that order with seed 4.
#'   Restore the caller's previous random state afterward, including when no
#'   random seed existed before the call.
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
