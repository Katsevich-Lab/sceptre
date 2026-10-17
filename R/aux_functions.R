#' Choose the default model formula from the available cell covariates
#'
#' Keep covariates in column order, excluding response mitochondrial fraction and,
#' when requested, guide count summaries. Use log(x) for names matching n_umis or
#' n_nonzero, adding one if any value is zero. Other covariates are included only
#' if they have fewer than 15 distinct values. Return a formula for model.matrix();
#' validation of the supplied covariates happens elsewhere.
#' @noRd
auto_construct_formula_object <- function(
    cell_covariates,
    include_grna_covariates
) {
    MAX_N_LEVELS_ALLOWED <- 15L
    cell_covariate_names <- colnames(cell_covariates)
    cell_covariate_names <- cell_covariate_names[
        cell_covariate_names != "response_p_mito"
    ]
    if (!include_grna_covariates) {
        # by default, do not use grna count-based covariates in low moi
        cell_covariate_names <- cell_covariate_names[
            !(cell_covariate_names %in% c("grna_n_umis", "grna_n_nonzero"))
        ]
    }
    form_str <- vapply(
        cell_covariate_names,
        function(curr_name) {
            count_based_covariate <- grepl(
                pattern = "n_umis|n_nonzero",
                x = curr_name
            )
            if (count_based_covariate) {
                if (any(cell_covariates[[curr_name]] == 0)) {
                    out <- paste0("log(", curr_name, "+1)")
                } else {
                    out <- paste0("log(", curr_name, ")")
                }
            } else {
                if (
                    length(unique(cell_covariates[[curr_name]])) >=
                        MAX_N_LEVELS_ALLOWED
                ) {
                    out <- NA_character_
                } else {
                    out <- curr_name
                }
            }
            return(out)
        },
        FUN.VALUE = character(1)
    ) |>
        stats::na.omit() |>
        paste0(collapse = " + ")
    form <- paste0("~ ", form_str) |> stats::as.formula()
    return(form)
}


#' Combine response summaries, guide summaries and extra cell covariates
#'
#' Response and guide matrices must have cells in the same order. Calculate their
#' summaries, prefix the column names with response_ or grna_, and append any extra
#' covariates. Use response-matrix cell IDs as row names when supplied. Matrix
#' storage may be converted for column access; cell order is kept.
#' @noRd
auto_compute_cell_covariates <- function(
    response_matrix,
    grna_matrix,
    extra_covariates,
    response_names
) {
    # compute the response covariates
    covariate_df <- compute_cell_covariates(
        matrix_in = response_matrix,
        feature_names = response_names,
        compute_p_mito = TRUE,
        compute_max_feature = FALSE
    )
    colnames(covariate_df) <- paste0("response_", colnames(covariate_df))

    # compute the grna covariates
    grna_covariate_df <- compute_cell_covariates(
        matrix_in = grna_matrix,
        feature_names = character(0),
        compute_p_mito = TRUE,
        compute_max_feature = TRUE
    )
    colnames(grna_covariate_df) <- paste0("grna_", colnames(grna_covariate_df))
    covariate_df <- cbind(covariate_df, grna_covariate_df)

    # if extra covariates have been provided, add those as well
    if (nrow(extra_covariates) >= 1L) {
        covariate_df <- cbind(covariate_df, extra_covariates)
    }
    cell_ids <- colnames(response_matrix)
    if (!is.null(cell_ids)) {
        rownames(covariate_df) <- cell_ids
    }
    return(covariate_df)
}


#' Divide gene or guide IDs into blocks for workers
#'
#' When parallel work has at least twice as many items as workers, shuffle using
#' seed 4 and split into blocks. Otherwise return one block in input order.
#' For n_processors = "auto", use the preferred fraction of available physical
#' cores once. Restore the caller's random state and temporary options afterward.
#' This prepares the blocks; it does not start workers.
#' @noRd
partition_response_ids <- function(response_ids, parallel, n_processors) {
    groups_set <- FALSE
    if (parallel) {
        if (identical(n_processors, "auto")) {
            # Apply the user's core fraction once, defaulting to half the cores.
            # availableCores() can apply that same option itself, so temporarily
            # set it to 1 while counting cores, then multiply by frac below.
            # Use an option rather than fraction = 1: that argument was added in
            # parallelly 1.47.0, while supported installations can be older.
            frac <- getOption("parallelly.availableCores.fraction", 0.5)
            n_processors <- max(
                1L,
                floor(
                    withr::with_options(
                        list(parallelly.availableCores.fraction = 1),
                        parallelly::availableCores(logical = FALSE)
                    ) *
                        frac
                )
            )
        }
        if (length(response_ids) >= 2 * n_processors) {
            s <- withr::with_seed(4, sample(response_ids))
            out <- split(
                s,
                cut(
                    seq_along(s),
                    n_processors,
                    labels = paste0("group_", seq(1, n_processors))
                )
            )
            groups_set <- TRUE
        }
    }
    if (!groups_set) {
        out <- list(group_1 = response_ids)
    }
    return(out)
}


#' Find or create the sceptre_logs subdirectory
#'
#' Return log_dir/sceptre_logs/ with a trailing slash. Create parent directories
#' as needed; individual worker log files are created by their callers.
#' @noRd
get_log_dir <- function(log_dir) {
    log_dir <- paste0(log_dir, "/sceptre_logs/")
    if (!dir.exists(log_dir)) {
        dir.create(log_dir, recursive = TRUE)
    }
    return(log_dir)
}
